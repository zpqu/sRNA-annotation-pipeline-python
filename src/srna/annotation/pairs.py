"""Overlap-pair finding between reads and features.

Thin wrapper around pyranges joins. pyranges' ``containment`` mode has
inconsistent strand semantics, so every mode is implemented as a strand-aware
(or strand-ignoring) *any-overlap* join plus an explicit containment filter,
which is exactly what the R pipeline's ``findOverlaps(type="within")`` computes.

pyranges 0.1.4 is unstable on multi-chromosome joins with pandas 3 (it chokes
when a chromosome is present in only one of the two inputs), so work is split
per chromosome exactly like the step-02 engine does.

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
    common_chroms = set(reads["chrom"].astype(str)) & set(features["chrom"].astype(str))
    if not common_chroms:
        return np.empty((0, 2), dtype=np.int64)
    parts: list[np.ndarray] = []
    for chrom in common_chroms:
        sub = reads[reads["chrom"].astype(str) == chrom]
        f = features[features["chrom"].astype(str) == chrom]
        if sub.empty or f.empty:
            continue
        parts.append(_join_chrom(sub, f, mode, ignore_strand))
    if not parts:
        return np.empty((0, 2), dtype=np.int64)
    pairs = np.vstack(parts)
    # each per-chromosome join is internally sorted, so re-sort globally to
    # honour the documented (read_idx, feat_idx) ordering.
    return pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]


def _join_chrom(
    reads: pd.DataFrame,
    features: pd.DataFrame,
    mode: str,
    ignore_strand: bool,
) -> np.ndarray:
    """Run the overlap join for a single chromosome."""
    keep_strand = not ignore_strand
    rpr = _to_pyranges(reads, _READ_COL, keep_strand)
    fpr = _to_pyranges(features, _FEAT_COL, keep_strand)
    joined = rpr.join(fpr, how="left").df
    joined = joined[joined[_FEAT_COL] >= 0]
    if joined.empty:
        return np.empty((0, 2), dtype=np.int64)
    if keep_strand:
        # pyranges 0.1.4 join() ignores the Strand column even when both inputs
        # are stranded (and `strandedness="same"` crashes), so enforce the
        # same-strand condition ourselves (matches findOverlaps' default).
        joined = joined[joined["Strand"] == joined["Strand_b"]]
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
