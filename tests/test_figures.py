"""Regression tests for figure rendering (samples on rows).

Step-04/05/06 figures place samples on rows and read flavor on columns, so a
grid with *n* samples has ``2 * n`` panels in ``n`` rows. Two past bugs:

* ``fig_dims`` was called with ``n = len(samples)`` instead of the total panel
  count, sizing a 2-row grid for a single row and truncating the bottom row.
* a bogus cleanup loop turned ``axis("off")`` on the second row of axes,
  hiding every tick mark and tick label for the second sample.

Both are guarded here: figures must be two rows tall and every panel must have
fully visible tick marks and tick labels within the canvas.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import srna.analysis.annotation_and_size_barplots as barplots
from srna.analysis import pirna_position_windows
from srna.analysis.annotation_and_size_barplots import (
    FLAVORS,
    READ_LEVELS,
    _figure_04_class,
    _figure_05_class,
)
from srna.analysis.pirna_position_windows import _overlap_summary_figure, _position_figure

SAMPLES = ["Cumulus-cells", "Granulosa-cells"]

_ROW_PX = int(4.1 * 300)  # one 4.1 in row at the 300 dpi save resolution


def _captured_figs(func, *args, **kwargs) -> dict[str, plt.Figure]:
    """Run *func* with ``save_figure`` stubbed to keep the figures open."""
    captured: dict[str, plt.Figure] = {}
    module = barplots if func in (_figure_04_class, _figure_05_class) else pirna_position_windows

    def fake_save(fig, base, dpi=300, close=False):
        base.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(base) + ".png", dpi=dpi)
        captured[base.name] = fig

    saved = module.save_figure
    module.save_figure = fake_save
    try:
        func(*args, **kwargs)
    finally:
        module.save_figure = saved
    assert captured, "no figure was saved"
    return captured


def _assert_two_rows_and_full_axes(captured: dict[str, plt.Figure]) -> None:
    for name, fig in captured.items():
        fig.set_dpi(300)
        fig.canvas.draw()
        width, height = fig.get_size_inches() * fig.get_dpi()
        assert len(fig.axes) == 2 * len(SAMPLES), f"{name}: panel count"
        for ax in fig.axes:
            xt = [t for t in ax.get_xticklabels() if t.get_text() != ""]
            yt = [t for t in ax.get_yticklabels() if t.get_text() != ""]
            assert len(xt) >= 1, f"{name}: {ax.get_title()}: no x tick labels"
            assert len(yt) >= 1, f"{name}: {ax.get_title()}: no y tick labels"
            assert all(t.get_visible() for t in xt), f"{name}: hidden x labels"
            assert all(t.get_visible() for t in yt), f"{name}: hidden y labels"
            for t in xt:
                bb = t.get_window_extent()
                assert bb.y0 >= 0, f"{name}: x label clipped below canvas"
            for t in yt:
                bb = t.get_window_extent()
                assert bb.x0 >= 0, f"{name}: y label clipped left of canvas"
        assert height == 2 * _ROW_PX, f"{name}: expected 2 rows, got {height}px"


def _png_height(path: Path) -> int:
    from PIL import Image

    with Image.open(path) as im:
        return im.size[1]


def _t2_block(items: list[str]) -> pd.DataFrame:
    rows = []
    for s in SAMPLES:
        for flavor in FLAVORS:
            for i, item in enumerate(items):
                rows.append(
                    {
                        "sample": s,
                        "item": item,
                        "flavor": flavor,
                        "Freq": 10 + i + (0 if s == "Cumulus-cells" else 1),
                    }
                )
    return pd.DataFrame(rows)


class TestFigure04:
    def test_two_sample_grid_renders_both_rows(self, tmp_path):
        sub = _t2_block(READ_LEVELS[:6])
        captured = _captured_figs(_figure_04_class, "matmiRNA", sub, READ_LEVELS, SAMPLES, tmp_path)
        assert set(captured) == {
            "Figure_04b.matmiRNA_annotation_count_barplot",
            "Figure_04b.matmiRNA_annotation_percentage_barplot",
        }
        _assert_two_rows_and_full_axes(captured)
        for name, fig in captured.items():
            for ax in fig.axes:
                labels = {t.get_text() for t in ax.get_xticklabels() if t.get_text() != ""}
                assert labels == set(READ_LEVELS), (
                    f"{name}: {ax.get_title()}: x axis must show every category "
                    f"including zero-count ones (missing {set(READ_LEVELS) - labels})"
                )
        assert (
            _png_height(tmp_path / "Figure_04b.matmiRNA_annotation_count_barplot.png")
            == 2 * _ROW_PX
        )


class TestFigure05:
    def test_two_sample_grid_renders_both_rows(self, tmp_path):
        sub = _t2_block([str(sz) for sz in range(16, 26)])
        sub["item"] = sub["item"].astype(int)
        breaks = np.asarray([16, 21, 26])
        captured = _captured_figs(_figure_05_class, "matmiRNA", sub, SAMPLES, tmp_path, breaks)
        assert set(captured) == {
            "Figure_05a.matmiRNA_size_barplot",
            "Figure_05b.matmiRNA_size_barplot.percentage",
        }
        _assert_two_rows_and_full_axes(captured)
        assert _png_height(tmp_path / "Figure_05a.matmiRNA_size_barplot.png") == 2 * _ROW_PX


class TestOverlapSummaryFigure:
    def test_legend_does_not_overlap_title(self, tmp_path):
        cats = list(pirna_position_windows._OVERLAP_LEVELS)
        ov = pd.DataFrame(
            [
                {"sample": s, "category": c, "n_unique": 1, "n_reads": 2}
                for s in SAMPLES
                for c in cats
            ]
        )
        captured = _captured_figs(_overlap_summary_figure, ov, SAMPLES, tmp_path)
        assert set(captured) == {"Figure_06.piRNA_overlap_summary"}
        fig = captured["Figure_06.piRNA_overlap_summary"]
        fig.set_dpi(300)
        fig.canvas.draw()
        width, height = fig.get_size_inches() * fig.get_dpi()
        bb = fig.legends[0].get_window_extent()
        st = fig._suptitle.get_window_extent()
        lt = fig.legends[0].get_title().get_window_extent()
        assert bb.x0 >= 0 and bb.x1 <= width and bb.y0 >= 0 and bb.y1 <= height
        assert min(st.y1, lt.y1) - max(st.y0, lt.y0) <= 0, "legend title overlaps suptitle"


class TestPositionFigure:
    def test_two_sample_grid_renders_both_rows(self, tmp_path):
        profile = pd.DataFrame(
            {
                "position": [1, 2] * 2,
                "sample": [s for s in SAMPLES for _ in range(2)],
                "flavor": ["unique reads"] * 4,
                "sense": [1.0, 2.0, 3.0, 4.0],
                "antisense": [0.5, 1.0, 1.5, 2.0],
            }
        )
        captured = _captured_figs(_position_figure, profile, SAMPLES, "tRNA", tmp_path)
        assert set(captured) == {
            "Figure_06.piRNA_vs_tRNA_pos_barplot",
            "Figure_06.piRNA_vs_tRNA_pos_barplot_AS",
        }
        _assert_two_rows_and_full_axes(captured)
        assert _png_height(tmp_path / "Figure_06.piRNA_vs_tRNA_pos_barplot.png") == 2 * _ROW_PX
