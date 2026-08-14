"""Reading raw annotation files (BED / GTF / GFF3) into DataFrames.

All parsers return 0-based half-open intervals (``start`` inclusive,
``end`` exclusive), matching pysam and pyranges. Large files are read with the
pyarrow CSV engine for speed and memory efficiency.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

_BED6_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


def read_bed6(path: Path) -> pd.DataFrame:
    """Read a 6-column BED file (chrom, start, end, name, score, strand)."""
    return _read_delim(path, _BED6_COLUMNS, {"chrom": "string", "start": "int32", "end": "int32"})


def read_bed12(path: Path) -> pd.DataFrame:
    """Read a 12-column BED file, returning the first six columns."""
    columns = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "item_rgb",
        "block_count",
        "block_sizes",
        "block_starts",
    ]
    df = _read_delim(path, columns, {"chrom": "string", "start": "int32", "end": "int32"})
    return df[_BED6_COLUMNS]


def read_gtf(path: Path) -> pd.DataFrame:
    """Read a GTF file into a DataFrame with 0-based intervals.

    Returns columns ``chrom, feature, start, end, strand, score,
    attributes`` (attributes preserved for key extraction).
    """
    df = _read_delim(
        path,
        ["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        {
            "chrom": "string",
            "feature": "string",
            "start": "int32",
            "end": "int32",
            "strand": "string",
        },
    )
    df["start"] = df["start"] - 1
    return df


def read_gff3(path: Path) -> pd.DataFrame:
    """Read a GFF3 file into a DataFrame with 0-based intervals."""
    return read_gtf(path)


def extract_gtf_fields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """Extract quoted attribute values (``key "value"``) from a GTF attributes column."""
    attrs = df["attributes"].astype("string")
    out: dict[str, pd.Series] = {}
    for field in fields:
        pattern = rf'{field}\s+"([^"]*)"'
        out[field] = attrs.str.extract(pattern, expand=False).astype("string").fillna("")
    return pd.DataFrame(out)


def extract_gff_fields(df: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """Extract ``key=value`` attributes from a GFF3 attributes column."""
    attrs = df["attributes"].astype("string")
    out: dict[str, pd.Series] = {}
    for field in fields:
        pattern = rf"(?:^|;){field}=([^;]*)"
        out[field] = attrs.str.extract(pattern, expand=False).astype("string").fillna("")
    return pd.DataFrame(out)


def _read_delim(path: Path, columns: list[str], dtypes: dict[str, Any]) -> pd.DataFrame:
    """Read a tab-delimited file with pyarrow CSV (skips ``#`` comment lines)."""
    try:
        import pyarrow.csv as pacsv

        table = pacsv.read_csv(
            str(path),
            parse_options=pacsv.ParseOptions(delimiter="\t", comment_char="#"),
            convert_options=pacsv.ConvertOptions(
                column_names=columns,
                include_columns=columns,
                string_columns=[c for c in columns if c in dtypes],
            ),
        )
        df = table.to_pandas()
    except (ImportError, Exception):
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            names=columns,
            usecols=range(len(columns)),
            dtype=dtypes,
        )
    return df
