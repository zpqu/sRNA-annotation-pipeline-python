"""Collapse BAM reads into non-redundant reads with per-position counts.

The BAM is streamed contig by contig (never loaded whole). Reads sharing the
same (chrom, start, end, strand) are reported once with a ``count`` column.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from loguru import logger
from tqdm import tqdm

_STRAND_DT = np.dtype([("start", "<i8"), ("end", "<i8"), ("strand", "i1")])


def collapse_bam(bam_path: Path, thread: int = 4) -> tuple[pd.DataFrame, int]:
    """Collapse a BAM into non-redundant reads.

    Args:
        bam_path: indexed BAM file.
        thread: pysam fetch thread count.

    Returns:
        ``(unique_reads, total_reads)`` where ``unique_reads`` is a DataFrame
        with columns ``chrom, start, end, strand, count`` (0-based half-open
        intervals) and ``total_reads`` is the number of alignments read.

    """
    frames: list[pd.DataFrame] = []
    total = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        contigs = [c for c in bam.references]
        for contig in tqdm(contigs, desc=str(bam_path.name), unit="contig", disable=None):
            starts: list[int] = []
            ends: list[int] = []
            strands: list[int] = []
            for read in bam.fetch(contig):
                starts.append(read.reference_start)
                ends.append(read.reference_end)
                strands.append(1 if read.is_reverse else 0)
            if not starts:
                continue
            total += len(starts)
            u_start, u_end, u_strand, counts = _unique_counts(
                np.asarray(starts, dtype=np.int64),
                np.asarray(ends, dtype=np.int64),
                np.asarray(strands, dtype=np.int8),
            )
            frames.append(
                pd.DataFrame(
                    {
                        "chrom": np.full(len(u_start), contig, dtype=object),
                        "start": u_start,
                        "end": u_end,
                        "strand": np.where(u_strand == 1, "-", "+"),
                        "count": counts.astype(np.int64),
                    }
                )
            )
    if not frames:
        logger.warning("no reads found in {bam_path}", bam_path=bam_path)
        return (
            pd.DataFrame(columns=["chrom", "start", "end", "strand", "count"]),
            0,
        )
    unique = pd.concat(frames, ignore_index=True)
    unique["strand"] = unique["strand"].astype("string")
    return unique, total


def _unique_counts(
    start: np.ndarray, end: np.ndarray, strand: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return unique (start, end, strand) rows and their counts (sorted)."""
    keys = np.empty(len(start), dtype=_STRAND_DT)
    keys["start"] = start
    keys["end"] = end
    keys["strand"] = strand
    uniq, counts = np.unique(keys, return_counts=True)
    return uniq["start"], uniq["end"], uniq["strand"], counts


def add_poskey(df: pd.DataFrame) -> pd.DataFrame:
    """Add the ``poskey`` column (chrom:start:end:strand, 1-based)."""
    out = df.copy()
    out["poskey"] = (
        out["chrom"].astype(str)
        + ":"
        + (out["start"] + 1).astype(str)
        + ":"
        + out["end"].astype(str)
        + ":"
        + out["strand"].astype(str)
    )
    return out
