"""Overlap-pair finding between reads and features.

Thin wrapper around pyranges joins. pyranges' ``containment`` mode has
inconsistent strand semantics, so every mode is implemented as a strand-aware
(or strand-ignoring) *any-overlap* join plus an explicit containment filter,
which is exactly what the R pipeline's ``findOverlaps(type="within")`` computes.

All coordinates are 0-based half-open intervals. Both inputs must carry the
``read_idx`` / ``feat_idx`` identifier columns.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyranges as pr

#: identifier columns used by the engine.
_READ_COL = "read_idx"
_FEAT_COL = "feat_idx"


def overlap_pairs(
    reads: pd.DataFrame,
    features: pd.DataFrame,
    mode: str,
    ignore_strand: bool = False,
) -> np.ndarray:
    """Return the ``(read_idx, feat_idx)`` overlap pairs as an ``(n, 2)`` array.

    Args:
        reads: read table with ``chrom, start, end, strand, read_idx`` columns.
        features: feature table with ``chrom, start, end, strand, feat_idx``.
        mode: ``"any"`` (>=1 bp overlap), ``"within"`` (read inside feature) or
            ``"reverse_within"`` (feature inside read).
        ignore_strand: when True the strand condition is dropped.

    Returns:
        Array of pairs sorted by ``(read_idx, feat_idx)``.

    """
    keep_strand = not ignore_strand
    rpr = _to_pyranges(reads, _READ_COL, keep_strand)
    fpr = _to_pyranges(features, _FEAT_COL, keep_strand)
    joined = rpr.join(fpr, how="left").df
    joined = joined[joined[_FEAT_COL] >= 0]
    if joined.empty:
        return np.empty((0, 2), dtype=np.int64)
    if mode == "within":
        contained = (joined["Start_b"] <= joined["Start"]) & (joined["End_b"] >= joined["End"])
        joined = joined[contained]
    elif mode == "reverse_within":
        contained = (joined["Start"] <= joined["Start_b"]) & (joined["End"] >= joined["End_b"])
        joined = joined[contained]
    elif mode != "any":
        raise ValueError(f"unknown overlap mode: {mode!r}")
    pairs = joined[[_READ_COL, _FEAT_COL]].to_numpy(dtype=np.int64)
    return pairs


def _to_pyranges(df: pd.DataFrame, idx_col: str, keep_strand: bool) -> pr.PyRanges:
    """Convert an engine table into a PyRanges preserving *idx_col*."""
    cols = {
        "Chromosome": df["chrom"].to_numpy(),
        "Start": df["start"].to_numpy(),
        "End": df["end"].to_numpy(),
    }
    if keep_strand:
        cols["Strand"] = df["strand"].to_numpy()
    cols[idx_col] = df[idx_col].to_numpy()
    return pr.PyRanges(pd.DataFrame(cols))
