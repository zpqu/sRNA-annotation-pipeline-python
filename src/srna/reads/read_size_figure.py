"""Step-01 figure: read size / read count overview (Figure_01)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from srna.plotting import fig_dims, save_figure, small_fonts
from srna.reads.summaries import build_table_01b, subsample_counts

_ROW_H_MM = 69.0


def read_size_overview_figure(samples: dict[str, pd.DataFrame], base: Path) -> Path:
    """Render the read-size x read-count overview for all samples.

    Args:
        samples: sample label -> non-redundant reads (``count`` column).
        base: output path without extension (``Figure_01.read_size_vs_count``).

    Returns:
        *base* as a Path (both ``.pdf`` and ``.png`` are written).

    """  # noqa: D205
    n = len(samples)
    width, height = fig_dims(n, ncol=1, per_h=_ROW_H_MM / 25.4)
    fig, axes = plt.subplots(nrows=n, ncols=3, figsize=(width, height), squeeze=False)
    fig.suptitle("Read size vs count", fontsize=8)

    for row, (sample, unique) in enumerate(samples.items()):
        plot = pd.DataFrame({"width": unique["end"] - unique["start"], "count": unique["count"]})
        lc = subsample_counts(unique)
        width_counts = build_table_01b(unique, sample)

        ax_size, ax_dens, ax_2d = axes[row]

        _panel_size(ax_size, width_counts)
        _panel_density(ax_dens, lc)
        _panel_2d(ax_2d, plot)

        ax_size.set_title(sample, fontsize=8)
        is_bottom = row == n - 1
        for col_idx, ax in enumerate([ax_size, ax_dens, ax_2d]):
            is_left = col_idx == 0
            small_fonts(ax)
            if not is_bottom:
                ax.tick_params(axis="x", labelbottom=False)

    fig.tight_layout()
    save_figure(fig, base)
    return base


def _panel_size(ax, width_counts: pd.DataFrame) -> None:
    """Bar plot of n_unique reads per size."""
    ax.bar(width_counts["width"], width_counts["n_unique"], color="steelblue", width=0.9)
    ax.set_xlabel("size (nt)")
    ax.set_ylabel("n_unique")


def _panel_density(ax, lc: np.ndarray) -> None:
    """Density histogram of the log2 read count."""
    if len(lc) == 0:
        ax.text(0.5, 0.5, "no reads", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlabel("log2 count")
        return
    ax.hist(lc, bins=60, density=True, color="grey", alpha=0.6)
    ax.set_xlabel("log2 count")
    ax.set_ylabel("density")


def _panel_2d(ax, plot: pd.DataFrame) -> None:
    """2D heatmap of n_unique across size x log2 count."""
    if plot.empty:
        ax.text(0.5, 0.5, "no reads", ha="center", va="center", transform=ax.transAxes)
        return
    widths = plot["width"].to_numpy()
    lcs = np.round(np.log2(plot["count"].to_numpy()))
    counts, xedges, yedges = np.histogram2d(widths, lcs, bins=(60, 60))
    counts = np.ma.masked_where(counts == 0, counts)
    mesh = ax.pcolormesh(xedges, yedges, counts.T, cmap="YlGnBu", norm=LogNorm(vmin=1))
    ax.set_xlabel("size (nt)")
    ax.set_ylabel("log2 count")
    return mesh
