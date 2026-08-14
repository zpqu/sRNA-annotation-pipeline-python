"""Build the refGene-derived feature tables (``refGene.gtf``).

Port of the ``refGene GTF`` section of ``00_build_annotation_DB.R``:

* ``NM_`` transcripts -> exon / CDS / 5UTR / 3UTR (per transcript), introns,
  transcript spans (``NM.mRNA``) and +/- 1 kb flanks (``NM.up1k``/``NM.down1k``).
* ``NR_`` transcripts -> snoRNAs (gene name matching ``Snor*``), ``NR.exon`` and
  ``lincRNA.exon``.

All returned tables use 0-based half-open intervals with the unified feature
schema (see :mod:`srna.features.schema`).
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from .schema import FEATURE_COLUMNS

_SNORD_RE = re.compile(r"^Snor", re.IGNORECASE)


def parse_refgene_gtf(path: Any, genome: str) -> dict[str, pd.DataFrame]:
    """Parse a UCSC refGene GTF and return all derived feature tables.

    Args:
        path: path to the ``refGene.gtf`` file.
        genome: genome assembly id (used for provenance only).

    Returns:
        A mapping of canonical feature type -> feature DataFrame.

    """
    from .readers import extract_gtf_fields, read_gtf

    gtf = read_gtf(path)
    fields = extract_gtf_fields(gtf, ["gene_id", "transcript_id", "gene_name"])
    gtf = pd.concat([gtf, fields], axis=1)

    transcript_id = gtf["transcript_id"].astype(str)
    is_nm = transcript_id.str.startswith("NM_")
    is_nr = transcript_id.str.startswith("NR_")

    tables: dict[str, pd.DataFrame] = {}

    # --- NM_ (protein-coding) -------------------------------------------------
    nm = gtf[is_nm]
    nm_exon = nm[nm["feature"] == "exon"]
    nm_cds = nm[nm["feature"] == "CDS"]
    nm_5utr = nm[nm["feature"] == "5UTR"]
    nm_3utr = nm[nm["feature"] == "3UTR"]

    tables["NM.exon"] = _emit(nm_exon, "NM.exon", "refGene", genome)
    tables["NM.CDS"] = _emit(nm_cds, "NM.CDS", "refGene", genome)
    tables["NM.5UTR"] = _emit(nm_5utr, "NM.5UTR", "refGene", genome)
    tables["NM.3UTR"] = _emit(nm_3utr, "NM.3UTR", "refGene", genome)

    # transcript spans / introns / +/- 1kb flanks per NM_ transcript
    spans = nm_exon.groupby("transcript_id", sort=False).agg(
        chrom=("chrom", "first"),
        strand=("strand", "first"),
        start=("start", "min"),
        end=("end", "max"),
    )
    spans = spans.reset_index()
    spans["feature_type"] = "NM.mRNA"
    spans["feature_id"] = spans["transcript_id"]
    spans["gene_id"] = spans["transcript_id"].map(_first(nm_exon, "gene_id"))
    spans["gene_name"] = spans["transcript_id"].map(_first(nm_exon, "gene_name"))
    spans["source"] = "refGene"
    spans["score"] = np.nan
    tables["NM.mRNA"] = spans[FEATURE_COLUMNS].copy()

    introns = _intron_rows(nm_exon)
    tables["NM.intron"] = introns[FEATURE_COLUMNS].copy()

    for name, sign in (("NM.up1k", 1), ("NM.down1k", -1)):
        flanks = _flank_rows(spans, upstream=sign > 0, name=name)
        tables[name] = flanks[FEATURE_COLUMNS].copy()

    # --- NR_ (noncoding): snoRNAs + other ncRNA -------------------------------
    nr = gtf[is_nr]
    nr_exon = nr[nr["feature"] == "exon"]
    nr_snord = nr_exon["gene_name"].str.match(_SNORD_RE.pattern)
    tables["NR.exon"] = _emit(nr_exon, "NR.exon", "refGene", genome)
    tables["snoRNA"] = _emit(nr_exon[nr_snord], "snoRNA", "refGene", genome)
    tables["lincRNA.exon"] = _emit(nr_exon[~nr_snord], "lincRNA.exon", "refGene", genome)

    return tables


def _emit(df: pd.DataFrame, feature_type: str, source: str, genome: str) -> pd.DataFrame:
    """Convert a parsed GTF subset into the unified feature schema."""
    if df.empty:
        return _empty(feature_type, source, genome)
    return pd.DataFrame(
        {
            "chrom": df["chrom"].astype(str),
            "start": df["start"].astype(np.int64),
            "end": df["end"].astype(np.int64),
            "strand": df["strand"].astype(str),
            "feature_type": feature_type,
            "feature_id": df["transcript_id"].astype(str),
            "gene_id": df["gene_id"].astype(str),
            "gene_name": df["gene_name"].astype(str),
            "source": source,
            "score": np.nan,
        }
    )[FEATURE_COLUMNS].reset_index(drop=True)


def _empty(feature_type: str, source: str, genome: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": pd.Series(dtype="string"),
            "start": pd.Series(dtype=np.int64),
            "end": pd.Series(dtype=np.int64),
            "strand": pd.Series(dtype="string"),
            "feature_type": pd.Series(dtype="string"),
            "feature_id": pd.Series(dtype="string"),
            "gene_id": pd.Series(dtype="string"),
            "gene_name": pd.Series(dtype="string"),
            "source": pd.Series(dtype="string"),
            "score": pd.Series(dtype="float64"),
        }
    )


def _first(df: pd.DataFrame, column: str) -> dict[str, str]:
    """Map each transcript_id to its first value of *column*."""
    return df.groupby("transcript_id", sort=False)[column].first().to_dict()


def _intron_rows(exon_df: pd.DataFrame) -> pd.DataFrame:
    """Compute intron intervals (gaps between exons within a transcript)."""
    if exon_df.empty:
        return _empty("NM.intron", "refGene", "")
    rows: list[dict[str, Any]] = []
    for _tx, exons in exon_df.sort_values(["transcript_id", "start"]).groupby(
        "transcript_id", sort=False
    ):
        starts = exons["start"].to_numpy(dtype=np.int64)
        ends = exons["end"].to_numpy(dtype=np.int64)
        for i in range(len(starts) - 1):
            gap_start = ends[i]
            gap_end = starts[i + 1]
            if gap_end > gap_start:
                rows.append(
                    {
                        "chrom": exons["chrom"].iloc[0],
                        "start": gap_start,
                        "end": gap_end,
                        "strand": exons["strand"].iloc[0],
                        "feature_type": "NM.intron",
                        "feature_id": _tx,
                        "gene_id": exons["gene_id"].iloc[0],
                        "gene_name": exons["gene_name"].iloc[0],
                        "source": "refGene",
                        "score": np.nan,
                    }
                )
    if not rows:
        return _empty("NM.intron", "refGene", "")
    return pd.DataFrame(rows)[FEATURE_COLUMNS]


def _flank_rows(spans: pd.DataFrame, upstream: bool, name: str) -> pd.DataFrame:
    """Build +/- 1 kb flank intervals from transcript spans.

    ``upstream=True`` gives ``NM.up1k`` (TSS - 1000 .. TSS - 1), otherwise
    ``NM.down1k`` (TTS + 1 .. TTS + 1000).
    """
    if spans.empty:
        return _empty(name, "refGene", "")
    start = spans["start"].to_numpy(dtype=np.int64)
    end = spans["end"].to_numpy(dtype=np.int64)
    plus = spans["strand"].to_numpy() == "+"
    minus = ~plus
    flank_start = np.zeros(len(spans), dtype=np.int64)
    flank_end = np.zeros(len(spans), dtype=np.int64)
    if upstream:
        flank_start[plus] = start[plus] - 1000
        flank_end[plus] = start[plus]
        flank_end[minus] = end[minus] + 1000
        flank_start[minus] = end[minus]
    else:
        flank_start[plus] = end[plus]
        flank_end[plus] = end[plus] + 1000
        flank_end[minus] = start[minus]
        flank_start[minus] = start[minus] - 1000
    out = spans.copy()
    out["start"] = flank_start
    out["end"] = flank_end
    out["feature_type"] = name
    return out[FEATURE_COLUMNS]
