"""Build the piRNA feature table from piRBase (BED) + piRNAdb (GTF).

The piRBase BED is the reference set and can contain ~80 M loci, so it is
streamed and kept with memory-efficient (pyarrow-backed) dtypes. piRNAdb
entries whose coordinates (chrom, start, end, strand) exactly match a piRBase
entry are dropped; loci duplicated within piRNAdb itself are dropped as well.
The remaining piRNAdb loci are appended to the piRBase set, mirroring the R
DB-build script.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv

_PIRBASE_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]
_PIRBASE_TYPES = {
    "chrom": pa.string(),
    "start": pa.int32(),
    "end": pa.int32(),
    "name": pa.string(),
    "score": pa.int32(),
    "strand": pa.string(),
}

_PRIMARY_DROP_PATTERN = r"_alt|_random|_fix|_hap|chrUn|^Un|^HSCHR"


def parse_pirna(
    pirbase_path: Any,
    pirnadb_path: Any,
    genome: str,
    chr_style: str = "chr",
) -> pa.Table:
    """Parse piRBase BED + piRNAdb GTF into a merged piRNA feature table.

    Returns a pyarrow Table in the unified feature schema (order preserved:
    piRBase rows first, then the appended non-redundant piRNAdb loci).
    """
    pirbase = _read_pirbase(pirbase_path, chr_style)
    pirnadb = _read_pirnadb(pirnadb_path, chr_style)

    # drop piRNAdb loci whose coordinates exactly match a piRBase locus
    pirnadb = _drop_matching(pirnadb, pirbase)
    # drop loci duplicated within piRNAdb itself (same coords, different codes)
    dups = pirnadb["start"].duplicated() & pirnadb["end"].duplicated()
    pirnadb = pirnadb[~dups]

    return _merge(pirbase, pirnadb, genome)


def _read_pirbase(path: Any, chr_style: str) -> pd.DataFrame:
    """Stream the piRBase BED into a normalized, primary-only DataFrame."""
    reader = pacsv.open_csv(
        str(path),
        read_options=pacsv.ReadOptions(column_names=_PIRBASE_COLUMNS),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=_PIRBASE_TYPES),
    )
    batches: list[pd.DataFrame] = []
    for batch in reader:
        df = batch.to_pandas()
        df["chrom"] = _normalize(df["chrom"], chr_style)
        df = df[_is_primary(df["chrom"])]
        df["score"] = df["score"].fillna(0).astype(np.int32)
        batches.append(df)
    out = pd.concat(batches, ignore_index=True)
    out["chrom"] = out["chrom"].astype("string[pyarrow]")
    out["name"] = out["name"].astype("string[pyarrow]")
    out["strand"] = out["strand"].astype("string[pyarrow]")
    return out


def _read_pirnadb(path: Any, chr_style: str) -> pd.DataFrame:
    """Read and normalize the piRNAdb GTF (non-comment lines only)."""
    from .readers import extract_gtf_fields, read_gtf

    gtf = read_gtf(path)
    gtf = gtf[gtf["feature"] == "piRNA"]
    code = extract_gtf_fields(gtf, ["piRNA_code"])["piRNA_code"]
    out = pd.DataFrame(
        {
            "chrom": gtf["chrom"],
            "start": gtf["start"].astype(np.int32),
            "end": gtf["end"].astype(np.int32),
            "strand": gtf["strand"],
            "name": code.astype("string[pyarrow]"),
        }
    )
    out["chrom"] = _normalize(out["chrom"], chr_style)
    out["chrom"] = out["chrom"].astype("string[pyarrow]")
    out["strand"] = out["strand"].astype("string[pyarrow]")
    return out[_is_primary(out["chrom"])].reset_index(drop=True)


def _drop_matching(pirnadb: pd.DataFrame, pirbase: pd.DataFrame) -> pd.DataFrame:
    """Return *pirnadb* rows whose coords are absent from *pirbase* (per chrom)."""
    keep: list[pd.DataFrame] = []
    for chrom, group in pirnadb.groupby("chrom", sort=False, dropna=False):
        base = pirbase[pirbase["chrom"] == chrom]
        if base.empty:
            keep.append(group)
            continue
        match = np.isin(
            _coords_keys(
                group["start"].to_numpy(), group["end"].to_numpy(), group["strand"].to_numpy()
            ),
            _coords_keys(
                base["start"].to_numpy(), base["end"].to_numpy(), base["strand"].to_numpy()
            ),
        )
        keep.append(group[~match])
    if not keep:
        return pirnadb.iloc[0:0]
    return pd.concat(keep, ignore_index=True)


def _coords_keys(start: np.ndarray, end: np.ndarray, strand: np.ndarray) -> np.ndarray:
    """Build a structured-array key over (start, end, strand) for isin lookups."""
    dtype = np.dtype([("start", "<i4"), ("end", "<i4"), ("strand", "S1")])
    out = np.empty(len(start), dtype=dtype)
    out["start"] = start.astype(np.int32)
    out["end"] = end.astype(np.int32)
    out["strand"] = np.char.encode(np.asarray(strand, dtype="U1"))
    return out


def _merge(pirbase: pd.DataFrame, pirnadb: pd.DataFrame, genome: str) -> pa.Table:
    """Combine piRBase + appended piRNAdb rows into the unified schema."""
    n_base = len(pirbase)
    n_add = len(pirnadb)
    chrom = pa.concat_arrays(
        [pa.array(pirbase["chrom"].to_numpy()), pa.array(pirnadb["chrom"].to_numpy())]
    )
    start = pa.concat_arrays(
        [
            pa.array(pirbase["start"].to_numpy().astype(np.int64)),
            pa.array(pirnadb["start"].to_numpy().astype(np.int64)),
        ]
    )
    end = pa.concat_arrays(
        [
            pa.array(pirbase["end"].to_numpy().astype(np.int64)),
            pa.array(pirnadb["end"].to_numpy().astype(np.int64)),
        ]
    )
    strand = pa.concat_arrays(
        [pa.array(pirbase["strand"].to_numpy()), pa.array(pirnadb["strand"].to_numpy())]
    )
    names = pa.concat_arrays(
        [pa.array(pirbase["name"].to_numpy()), pa.array(pirnadb["name"].to_numpy())]
    )
    score = pa.concat_arrays(
        [
            pa.array(pirbase["score"].to_numpy(dtype=np.float64)),
            pa.array(np.full(n_add, np.nan, dtype=np.float64)),
        ]
    )
    source = pa.concat_arrays(
        [
            pa.array(np.full(n_base, "piRBase", dtype=object)),
            pa.array(np.full(n_add, "piRNAdb", dtype=object)),
        ]
    )
    n_total = n_base + n_add
    return pa.Table.from_arrays(
        [
            chrom,
            start,
            end,
            strand,
            pa.array(np.full(n_total, "piRNA", dtype=object)),
            names,
            pa.array(np.full(n_total, "", dtype=object)),
            pa.array(np.full(n_total, "", dtype=object)),
            source,
            score,
        ],
        names=[
            "chrom",
            "start",
            "end",
            "strand",
            "feature_type",
            "feature_id",
            "gene_id",
            "gene_name",
            "source",
            "score",
        ],
    )


def _normalize(values: pd.Series, chr_style: str) -> pd.Series:
    """Normalize a pandas string column of chromosome names."""
    from ..chrom import normalize_seqnames

    return pd.Series(
        normalize_seqnames(values.astype("string").fillna("").tolist(), chr_style=chr_style),
        index=values.index,
    )


def _is_primary(chroms: pd.Series) -> pd.Series:
    """Return a boolean mask keeping primary chromosomes (regex search)."""
    matches = pc.match_substring_regex(pa.array(chroms.fillna("")), _PRIMARY_DROP_PATTERN)
    return ~pd.Series(matches.to_numpy(zero_copy_only=False), index=chroms.index)
