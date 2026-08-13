"""Step-01 summaries: read-size / read-count distribution tables (Table_01a-d)."""

from __future__ import annotations

import numpy as np
import pandas as pd

_COUNT_COL = "count"
_WIDTH_COL = "width"


def build_table_01a(unique: pd.DataFrame, total_reads: int, sample: str) -> dict[str, object]:
    """Build the Table_01a sample-summary row (mirrors the R logic).

    Args:
        unique: non-redundant reads (``count`` column required).
        total_reads: number of alignments in the BAM.
        sample: sample label.

    """
    plot = pd.DataFrame({_WIDTH_COL: unique["end"] - unique["start"], _COUNT_COL: unique["count"]})
    width_counts = plot.groupby(_WIDTH_COL)[_COUNT_COL].count().sort_values(ascending=False)
    width_reads = plot.groupby(_WIDTH_COL)[_COUNT_COL].sum().sort_values(ascending=False)
    top5 = ";".join(str(w) for w in width_counts.index[:5])
    return {
        "sample": sample,
        "total_reads": int(total_reads),
        "unique_reads": len(plot),
        "singleton_reads": int((plot[_COUNT_COL] == 1).sum()),
        "pct_singletons": round(100 * (plot[_COUNT_COL] == 1).mean(), 2),
        "median_count": float(plot[_COUNT_COL].median()),
        "max_count": int(plot[_COUNT_COL].max()),
        "dominant_size_by_unique_nt": (
            int(width_counts.index[0]) if not width_counts.empty else np.nan
        ),
        "dominant_size_by_unique_n": (
            int(width_counts.iloc[0]) if not width_counts.empty else np.nan
        ),
        "dominant_size_by_reads_nt": int(width_reads.index[0]) if not width_reads.empty else np.nan,
        "dominant_size_by_reads_n": int(width_reads.iloc[0]) if not width_reads.empty else np.nan,
        "top5_sizes_nt": top5,
    }


def build_table_01b(unique: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Read-size distribution: ``n_unique`` and ``n_reads`` per read size."""
    plot = pd.DataFrame({_WIDTH_COL: unique["end"] - unique["start"], _COUNT_COL: unique["count"]})
    out = plot.groupby(_WIDTH_COL, as_index=False).agg(
        n_unique=(_COUNT_COL, "count"), n_reads=(_COUNT_COL, "sum")
    )
    out.insert(0, "sample", sample)
    return out


def build_table_01c(unique: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Read-count distribution: ``n_unique`` per rounded ``log2(count)`` bin."""
    plot = pd.DataFrame({_COUNT_COL: unique["count"]})
    plot["lc"] = np.round(np.log2(plot[_COUNT_COL]))
    out = plot.groupby("lc", as_index=False).agg(n_unique=(_COUNT_COL, "count"))
    out.insert(0, "sample", sample)
    return out


def build_table_01d(unique: pd.DataFrame, sample: str) -> pd.DataFrame:
    """2D read-size x read-count grid: ``n_unique`` per (width, log2-count) bin."""
    plot = pd.DataFrame({_WIDTH_COL: unique["end"] - unique["start"], _COUNT_COL: unique["count"]})
    plot["lc"] = np.round(np.log2(plot[_COUNT_COL]))
    out = plot.groupby([_WIDTH_COL, "lc"], as_index=False).agg(n_unique=(_COUNT_COL, "count"))
    out.insert(0, "sample", sample)
    return out


def subsample_counts(
    unique: pd.DataFrame, max_points: int = 1_000_000, seed: int = 123
) -> np.ndarray:
    """Return the ``log2(count)`` values of a seeded subsample of unique reads."""
    rng = np.random.default_rng(seed)
    n = len(unique)
    keep = rng.choice(n, size=min(n, max_points), replace=False)
    return np.log2(unique["count"].to_numpy()[keep])
