"""Parse tRNA (tRNAscan-SE BED12) and RepeatMasker (BED6) feature files."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .base import FEATURE_COLUMNS


def parse_trna_bed(path: Any, genome: str) -> pd.DataFrame:
    """Parse a tRNAscan-SE BED12 file into a feature table.

    Interval coordinates from BED are 0-based half-open already.
    """
    from .parsing import read_bed12

    bed = read_bed12(path)
    return _emit(bed, "tRNA", f"{genome}_tRNAs")


def parse_rm_bed(path: Any, genome: str) -> pd.DataFrame:
    """Parse a RepeatMasker BED6 file into a feature table."""
    from .parsing import read_bed6

    bed = read_bed6(path)
    return _emit(bed, "RM", f"{genome}_rmsk")


def _emit(bed: pd.DataFrame, feature_type: str, source: str) -> pd.DataFrame:
    if bed.empty:
        return _empty(feature_type, source)
    return pd.DataFrame(
        {
            "chrom": bed["chrom"].astype(str),
            "start": bed["start"].astype(np.int64),
            "end": bed["end"].astype(np.int64),
            "strand": bed["strand"].astype(str),
            "feature_type": feature_type,
            "feature_id": bed["name"].astype(str),
            "gene_id": bed["name"].astype(str),
            "gene_name": bed["name"].astype(str),
            "source": source,
            "score": pd.to_numeric(bed["score"], errors="coerce"),
        }
    )[FEATURE_COLUMNS].reset_index(drop=True)


def _empty(feature_type: str, source: str) -> pd.DataFrame:
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
