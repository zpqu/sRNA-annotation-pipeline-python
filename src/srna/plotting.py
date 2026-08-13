"""Shared plotting helpers: publication figure sizing and PDF/PNG saving.

Port of the R helpers in ``scripts/R/lib/init.R`` (``fig.dims``, ``small.font``)
and a consistent save routine used by every figure-producing step.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

WIDTH_MM = 180.0
HEIGHT_MM = 220.0

#: orientation colour scheme used by step-04 figures.
ORIENT_COLS = {"sense": "#1f78b4", "antisense": "#e31a1c", "other": "grey50"}

#: per-category colour scheme used by step-03 / step-06 figures.
CATEGORY_COLS = {
    "matmiRNA": "#1f78b4",
    "piRNA": "#33a02c",
    "tRNA": "#e31a1c",
    "snoRNA": "#ff7f00",
}

#: overlap-rule colour scheme used by step-s01 figures.
STRATEGY_COLS = {"fully-contained": "#33a02c", "union": "#1f78b4", "any": "#e31a1c"}


def fig_dims(n: int, ncol: int, per_h: float = 3.2, stack: int = 1) -> tuple[float, float]:
    """Return ``(width, height)`` in inches for a faceted figure.

    Mirrors ``fig.dims`` in the R bootstrap: fixed 180 mm width, height derived
    from the number of facet rows and capped at 220 mm.
    """
    width = WIDTH_MM / 25.4
    ncol = max(1, min(ncol, n))
    rows = math.ceil(n / ncol)
    height = min(stack * rows * per_h, HEIGHT_MM / 25.4)
    return width, height


def save_figure(fig, base: Path, dpi: int = 300, close: bool = True) -> None:
    """Save a matplotlib figure as both PDF and PNG at *base*."""
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(base) + ".pdf")
    fig.savefig(str(base) + ".png", dpi=dpi)
    if close:
        plt.close(fig)


def small_fonts(ax) -> None:
    """Apply the small publication fonts used on 180 mm figures."""
    for item in (
        [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()
    ):
        item.set_fontsize(7)
    if ax.legend_ is not None:
        ax.legend_.set_fontsize(7)


def rotate_xticks(ax, angle: int = 45) -> None:
    """Rotate x tick labels to avoid overlap."""
    ax.tick_params(axis="x", rotation=angle)
