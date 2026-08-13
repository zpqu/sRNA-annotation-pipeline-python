"""Step-02 annotation engine (port of ``02_annotation_smallRNA.R``).

Reads are annotated as non-redundant unique reads. A sense pass walks the
genomic-feature priority list (matmiRNA > piRNA > snoRNA > tRNA > RM >
NM.exon > NM.intron > lincRNA.exon), assigning each overlapping read a category,
feature id, the number of overlapping features, and (for the four
small-RNA categories) a gene context region. An antisense pass then annotates
the remaining reads with any-overlap/ignore-strand features. Everything left
over is ``other``.

Work is split per chromosome so the ~83M piRNA loci never need to be resident
in memory at once.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from loguru import logger

from srna.annotation.pairs import overlap_pairs
from srna.config import Strategy
from srna.features.base import (
    FEATURE_META,
    GENE_CONTEXT_CATEGORIES,
    GENE_FEATURES,
    GENE_REGION,
    GENOMIC_FEATURES,
    element_id_of,
)
from srna.features.store import FeatureStore

_SENSE_MODE: dict[Strategy, str] = {
    Strategy.FULLY_CONTAINED: "within",
    Strategy.UNION: "union",
    Strategy.ANY: "any",
}

#: feature types whose hits receive a gene-context region in step 02.
_PIRNA = "piRNA"


@dataclass
class AnnotationResult:
    """Step-02 outputs for one sample."""

    sample: str
    per_read: pd.DataFrame
    counts_unique: list[tuple[str, str, int]] = field(default_factory=list)
    counts_reads: list[tuple[str, str, int]] = field(default_factory=list)


class _FeatureCache:
    """Per-chromosome feature cache with lazy whole-table loading."""

    def __init__(self, store: FeatureStore):
        self._store = store
        self._whole: dict[str, pd.DataFrame] = {}
        self._chrom: dict[tuple[str, str], pd.DataFrame] = {}
        self._prepped: dict[tuple[str, str], pd.DataFrame] = {}

    def chrom(self, name: str, chrom: str) -> pd.DataFrame:
        key = (name, chrom)
        df = self._chrom.get(key)
        if df is None:
            if name == _PIRNA:
                df = self._store.chrom_table(name, chrom)
            else:
                whole = self._whole.get(name)
                if whole is None:
                    whole = self._store.table(name)
                    self._whole[name] = whole
                df = whole[whole["chrom"] == chrom]
            self._chrom[key] = df
        return df

    def prepped(self, name: str, chrom: str) -> pd.DataFrame:
        """Return the bookkeeping-enriched feature slice (memoized)."""
        key = (name, chrom)
        df = self._prepped.get(key)
        if df is None:
            df = _prep_features(self.chrom(name, chrom), name)
            self._prepped[key] = df
        return df


def _prep_features(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Add overlap bookkeeping columns to a feature table slice."""
    meta = FEATURE_META[name]
    df = df.copy()
    df["feat_idx"] = np.arange(len(df))
    df["elem_id"] = element_id_of(df, meta)
    if meta.kind == "grl":
        df["elem_rank"] = df.groupby("feature_id")["feat_idx"].transform("min")
    else:
        df["elem_rank"] = df["feat_idx"]
    return df


def _pairs(
    reads: pd.DataFrame,
    features: pd.DataFrame,
    mode: str,
    ignore_strand: bool = False,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Return ``(pairs, primary)`` where *primary* flags within-direction pairs."""
    if mode == "union" and not ignore_strand:
        p1 = overlap_pairs(reads, features, "within")
        p2 = overlap_pairs(reads, features, "reverse_within")
        if len(p1) == 0:
            return p2, np.zeros(len(p2), dtype=bool)
        if len(p2) == 0:
            return p1, np.ones(len(p1), dtype=bool)
        return (
            np.vstack([p1, p2]),
            np.concatenate([np.ones(len(p1), dtype=bool), np.zeros(len(p2), dtype=bool)]),
        )
    pairs = overlap_pairs(reads, features, mode, ignore_strand=ignore_strand)
    return pairs, None


def _select(
    pairs: np.ndarray,
    features: pd.DataFrame,
    is_pirna: bool,
    primary: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, pd.Series]:
    """Pick one feature id per hit read and count overlapping features.

    Returns ``(read_idx, feature_id, n_features)`` aligned to ascending read
    index. piRNA picks the highest-scoring overlapping locus (ties broken by
    feature order); other features pick the first feature in feature order.
    """
    if len(pairs) == 0:
        return (
            np.array([], dtype=np.int64),
            np.array([], dtype=object),
            pd.Series(dtype=np.int64),
        )
    pr = pd.DataFrame(
        {"read_idx": pairs[:, 0].astype(np.int64), "feat_idx": pairs[:, 1].astype(np.int64)}
    )
    pr = pr.merge(
        features[["feat_idx", "elem_rank", "elem_id", "score"]], on="feat_idx", how="left"
    )
    if is_pirna:
        pr["_score"] = pr["score"].fillna(-np.inf)
        pr = pr.sort_values(
            ["read_idx", "_score", "feat_idx"], ascending=[True, False, True], kind="stable"
        )
    else:
        pr["_primary"] = True if primary is None else primary
        pr = pr.sort_values(
            ["read_idx", "_primary", "feat_idx"], ascending=[True, True, True], kind="stable"
        )
    per_element = pr.drop_duplicates(subset=["read_idx", "elem_rank"], keep="first")
    pick = per_element.drop_duplicates(subset=["read_idx"], keep="first").sort_values(
        "read_idx", kind="stable"
    )
    n_features = per_element.groupby("read_idx", sort=False).size()
    n_features = n_features.reindex(pick["read_idx"]).astype(np.int64)
    return (
        pick["read_idx"].to_numpy(),
        pick["elem_id"].to_numpy(),
        n_features.reset_index(drop=True),
    )


def _assign_regions(
    hit_df: pd.DataFrame,
    cache: _FeatureCache,
    mode: str,
    chroms: list[str],
) -> pd.Series:
    """Port of ``mygeneFeature.R``: assign a gene-context region per hit read.

    Returns a Series keyed by the global read index. Sense regions are assigned
    in priority order (CDS > 5UTR > 3UTR > intron > up1k > down1k > RM), then
    antisense regions (``AS.*``) to what remains, then ``intergenic``.
    """
    region = pd.Series("", index=hit_df["read_idx"], dtype=object)
    remaining = hit_df

    def assign_round(features: list[str], prefix: str, ignore_strand: bool) -> None:
        nonlocal remaining
        for feat in features:
            reg = prefix + GENE_REGION[feat]
            for chrom in chroms:
                sub = remaining[remaining["chrom"] == chrom]
                if sub.empty:
                    continue
                f = cache.prepped(feat, chrom)
                if f.empty:
                    continue
                if ignore_strand:
                    pairs, _ = _pairs(sub, f, "any", ignore_strand=True)
                else:
                    pairs, primary = _pairs(sub, f, mode)
                if len(pairs) == 0:
                    continue
                ridx, _, _ = _select(pairs, f, False, None if ignore_strand else primary)
                if len(ridx):
                    region.loc[ridx] = reg
                    remaining = remaining[~remaining["read_idx"].isin(set(ridx.tolist()))]

    assign_round(GENE_FEATURES, "", False)
    assign_round(GENE_FEATURES, "AS.", True)
    return region.where(region != "", "intergenic")


def annotate_sample(
    reads: pd.DataFrame,
    store: FeatureStore,
    strategy: Strategy,
    sample: str,
) -> AnnotationResult:
    """Annotate one sample's non-redundant reads and build step-02 outputs."""
    mode = _SENSE_MODE[strategy]
    df = reads.copy()
    df["read_idx"] = np.arange(len(df))
    df["poskey"] = (
        df["chrom"].astype(str)
        + ":"
        + (df["start"] + 1).astype(str)
        + ":"
        + df["end"].astype(str)
        + ":"
        + df["strand"].astype(str)
    )
    chroms = [c for c in pd.unique(df["chrom"]) if c is not None]
    cache = _FeatureCache(store)
    active = np.ones(len(df), dtype=bool)
    records: list[pd.DataFrame] = []

    for name in GENOMIC_FEATURES:
        hit_ridx: list[np.ndarray] = []
        eids: list[np.ndarray] = []
        nfs: list[pd.Series] = []
        for chrom in chroms:
            sel = np.flatnonzero((df["chrom"].to_numpy() == chrom) & active)
            if len(sel) == 0:
                continue
            r = df.iloc[sel].copy()
            f = cache.prepped(name, chrom)
            if f.empty:
                continue
            pairs, primary = _pairs(r, f, mode)
            if len(pairs) == 0:
                continue
            ridx, eid, nf = _select(pairs, f, name == _PIRNA, primary)
            hit_ridx.append(ridx)
            eids.append(eid)
            nfs.append(nf)
        if len(hit_ridx) == 0:
            continue
        ridx_all = np.concatenate(hit_ridx)
        eid_all = np.concatenate(eids)
        nf_all = pd.concat(nfs, ignore_index=True)
        hits = df.iloc[ridx_all].copy()
        hits["category"] = name
        hits["feature_id"] = eid_all
        hits["n_features"] = nf_all.to_numpy()
        if name in GENE_CONTEXT_CATEGORIES:
            regions = _assign_regions(hits, cache, mode, chroms)
            hits["region"] = regions.loc[ridx_all].to_numpy()
        else:
            hits["region"] = "NA"
        records.append(hits)
        active[ridx_all] = False
        logger.info(
            "sample '{sample}': {name} annotated {n} unique reads",
            sample=sample,
            name=name,
            n=len(ridx_all),
        )

    for name in GENOMIC_FEATURES:
        category = "AS." + name
        hit_ridx: list[np.ndarray] = []
        eids: list[np.ndarray] = []
        nfs: list[pd.Series] = []
        for chrom in chroms:
            sel = np.flatnonzero((df["chrom"].to_numpy() == chrom) & active)
            if len(sel) == 0:
                continue
            r = df.iloc[sel].copy()
            f = cache.prepped(name, chrom)
            if f.empty:
                continue
            pairs, _ = _pairs(r, f, "any", ignore_strand=True)
            if len(pairs) == 0:
                continue
            ridx, eid, nf = _select(pairs, f, name == _PIRNA, None)
            hit_ridx.append(ridx)
            eids.append(eid)
            nfs.append(nf)
        if len(hit_ridx) == 0:
            continue
        ridx_all = np.concatenate(hit_ridx)
        eid_all = np.concatenate(eids)
        nf_all = pd.concat(nfs, ignore_index=True)
        hits = df.iloc[ridx_all].copy()
        hits["category"] = category
        hits["feature_id"] = eid_all
        hits["n_features"] = nf_all.to_numpy()
        hits["region"] = "NA"
        records.append(hits)
        active[ridx_all] = False
        logger.info(
            "sample '{sample}': {category} annotated {n} unique reads",
            sample=sample,
            category=category,
            n=len(ridx_all),
        )

    remaining = df.iloc[np.flatnonzero(active)].copy()
    remaining["category"] = "other"
    remaining["feature_id"] = np.nan
    remaining["n_features"] = pd.NA
    remaining["region"] = "NA"
    records.append(remaining)

    annotated = pd.concat(records, ignore_index=True)
    total = int(annotated["count"].sum())
    per_read = pd.DataFrame(
        {
            "sample": sample,
            "read_id": np.arange(len(annotated)) + 1,
            "chr": annotated["chrom"],
            "start": annotated["start"] + 1,
            "end": annotated["end"],
            "strand": annotated["strand"],
            "size": annotated["end"] - annotated["start"],
            "count": annotated["count"],
            "cpm": (annotated["count"] / total * 1e6).round(3),
            "category": annotated["category"],
            "feature_id": annotated["feature_id"].astype("string"),
            "gene_context": annotated["region"].astype("string"),
            "n_features": annotated["n_features"],
        }
    )
    per_read = per_read.sort_values(
        ["count", "chr", "start", "end"], ascending=[False, True, True, True], kind="stable"
    ).reset_index(drop=True)

    counts_unique, counts_reads = _consolidated_counts(annotated)
    return AnnotationResult(
        sample=sample,
        per_read=per_read,
        counts_unique=counts_unique,
        counts_reads=counts_reads,
    )


def _consolidated_counts(
    annotated: pd.DataFrame,
) -> tuple[list[tuple[str, str, int]], list[tuple[str, str, int]]]:
    """Build the Table_02a/02b rows, mirroring ``add.counts`` call order."""
    unique_rows: list[tuple[str, str, int]] = []
    read_rows: list[tuple[str, str, int]] = []
    counts = annotated["count"].to_numpy()

    def add(category: str, items) -> None:
        items = np.asarray(items)
        seen: dict[str, int] = {}
        rows_unique: list[tuple[str, int, int]] = []
        for item, c in zip(items, counts, strict=False):
            key = str(item)
            if key not in seen:
                seen[key] = len(rows_unique)
                rows_unique.append((key, 0, 0))
            row = rows_unique[seen[key]]
            rows_unique[seen[key]] = (row[0], row[1] + 1, row[2] + int(c))
        unique_rows.extend((category, k, u) for k, u, _ in rows_unique)
        read_rows.extend((category, k, r) for k, _, r in rows_unique)

    by_type = {t: annotated["category"] == t for t in GENE_CONTEXT_CATEGORIES}
    for t in GENE_CONTEXT_CATEGORIES:
        mask = by_type[t].to_numpy()
        add(f"{t}.annotation", annotated.loc[mask, "region"])
        add(f"{t}.size", (annotated.loc[mask, "end"] - annotated.loc[mask, "start"]).to_numpy())
    add("read.annotation", annotated["category"])
    add("read.size", (annotated["end"] - annotated["start"]).to_numpy())
    return unique_rows, read_rows
