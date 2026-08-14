"""Steps 04-05 figures (ports of ``04_figure_annotation.R`` / ``05_figure_size.R``).

Both steps read the consolidated count tables (``Table_02a`` unique reads,
``Table_02b`` all reads) produced by step 02 and draw faceted barplots with
samples on the rows and read flavor (unique/all) on the columns:

* step 04 -- per-class annotation composition: ``Figure_04{a-e}.<class>_..._barplot``
* step 05 -- per-class read-size distribution: ``Figure_05a/05b.<class>_size_barplot``
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from srna.plotting import ORIENT_COLS, fig_dims, save_figure, small_fonts

CLASS_MAP = {
    "read": "read.annotation",
    "matmiRNA": "matmiRNA.annotation",
    "snoRNA": "snoRNA.annotation",
    "piRNA": "piRNA.annotation",
    "tRNA": "tRNA.annotation",
}

READ_LEVELS = [
    "matmiRNA",
    "piRNA",
    "snoRNA",
    "tRNA",
    "RM",
    "refGene.NM.exon",
    "refGene.NM.intron",
    "lincRNA.exon",
    "AS.matmiRNA",
    "AS.piRNA",
    "AS.snoRNA",
    "AS.tRNA",
    "AS.RM",
    "AS.refGene.NM.exon",
    "AS.refGene.NM.intron",
    "AS.lincRNA.exon",
    "other",
]

REGION_LEVELS = [
    "CDS",
    "5UTR",
    "3UTR",
    "intron",
    "up1k",
    "down1k",
    "RM",
    "AS.CDS",
    "AS.5UTR",
    "AS.3UTR",
    "AS.intron",
    "AS.up1k",
    "AS.down1k",
    "AS.RM",
    "intergenic",
]

FLAVORS = ["unique reads", "all reads"]

_SIZE_CLASSES = ["read", "matmiRNA", "piRNA", "snoRNA", "tRNA"]


def _orientation(item: str) -> str:
    """Map an item label to sense / antisense / other (R ``orient.of``)."""
    if item.startswith("AS."):
        return "antisense"
    if item in ("intergenic", "other"):
        return "other"
    return "sense"


def _load_t2(tables_dir: Path) -> pd.DataFrame:
    """Load and combine Table_02a (unique) + Table_02b (all reads)."""
    unique = pd.read_csv(tables_dir / "Table_02a_annotation_count_unique_reads.csv")
    all_reads = pd.read_csv(tables_dir / "Table_02b_annotation_count_all_reads.csv")
    for df, flavor in ((unique, "unique reads"), (all_reads, "all reads")):
        df["flavor"] = flavor
    return pd.concat([unique, all_reads], ignore_index=True)


def _axis_bar(
    ax,
    items: list[str],
    values: np.ndarray,
    groups: list[str],
    log_scale: bool,
) -> None:
    """Draw a grouped-by-orientation bar chart on *ax*."""
    width = 0.7
    pos = 0.0
    for _item, value, _group in zip(items, values, groups, strict=False):
        ax.bar(
            pos,
            value,
            width=width,
            color=[ORIENT_COLS[g] for g in groups],
            log=log_scale,
        )
        pos += 1.0
    ax.set_xticks(np.arange(len(items)))
    ax.set_xticklabels(items, rotation=45, ha="right")
    small_fonts(ax)


def _figure_04_class(
    class_name: str,
    sub: pd.DataFrame,
    levels: list[str],
    samples: list[str],
    figures_dir: Path,
) -> None:
    """Draw count + percentage barplots for one step-04 class."""
    import matplotlib.pyplot as plt

    samples = list(sub["sample"].unique())
    ncol = 2
    nrow = len(samples)
    w, h = fig_dims(len(samples) * ncol, ncol, per_h=4.1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        row_sub = sub[sub["sample"] == s]
        for j, flavor in enumerate(FLAVORS):
            ax = axes[i][j]
            f = row_sub[row_sub["flavor"] == flavor]
            counts = {str(it): val for it, val in zip(f["item"], f["Freq"], strict=False)}
            items = list(levels)
            values = np.array([counts.get(lv, 0) for lv in items], dtype=float)
            groups = [_orientation(lv) for lv in items]
            if len(items):
                _axis_bar(ax, items, np.maximum(values, 1.0), groups, log_scale=True)
            ax.set_yscale("log")
            ax.set_title(f"{s} | {flavor}", fontsize=8)
            ax.set_ylabel("Count (log10 scale)", fontsize=8)
    fig.suptitle(f"{class_name} annotation - count", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(
        fig, figures_dir / f"Figure_04{_letter(class_name)}.{class_name}_annotation_count_barplot"
    )

    # percentage panel: complete all item x sample x flavor combos (zero-filled)
    full = sub.copy()
    full = full.groupby(["sample", "item", "flavor"], as_index=False)["Freq"].sum()
    full["item"] = full["item"].astype(str)
    for s in samples:
        for flavor in FLAVORS:
            present = set(full[(full["sample"] == s) & (full["flavor"] == flavor)]["item"])
            for lv in levels:
                if lv not in present:
                    full = pd.concat(
                        [
                            full,
                            pd.DataFrame(
                                {"sample": [s], "item": [lv], "flavor": [flavor], "Freq": [0]}
                            ),
                        ],
                        ignore_index=True,
                    )
    full["per"] = full.groupby(["sample", "flavor"], group_keys=False)["Freq"].transform(
        lambda x: x / x.sum()
    )
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        for j, flavor in enumerate(FLAVORS):
            ax = axes[i][j]
            f = full[(full["sample"] == s) & (full["flavor"] == flavor)]
            counts = dict(zip(f["item"], f["per"], strict=False))
            items = list(levels)
            values = np.array([counts.get(lv, 0) for lv in items])
            groups = [_orientation(lv) for lv in items]
            if len(items):
                _axis_bar(ax, items, values, groups, log_scale=False)
            ax.set_ylim(0, 1)
            ax.set_yticklabels([f"{tick:.0%}" for tick in ax.get_yticks()])
            ax.set_title(f"{s} | {flavor}", fontsize=8)
            ax.set_ylabel("Percentage", fontsize=8)
    fig.suptitle(f"{class_name} annotation - percentage", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(
        fig,
        figures_dir / f"Figure_04{_letter(class_name)}.{class_name}_annotation_percentage_barplot",
    )


def _letter(class_name: str) -> str:
    return "abcde"[list(CLASS_MAP).index(class_name)]


def _figure_05_class(
    class_name: str,
    sub: pd.DataFrame,
    samples: list[str],
    figures_dir: Path,
    size_breaks: np.ndarray,
) -> None:
    """Draw count + percentage read-size barplots for one step-05 class."""
    import matplotlib.pyplot as plt

    sub = sub.copy()
    sub["size"] = sub["item"].astype(int)
    nrow = len(samples)
    w, h = fig_dims(len(samples) * 2, 2, per_h=4.1)
    for kind, _yval, title, ylabel in (
        ("counts", "Freq", "counts", "Count"),
        ("percentage", None, "percentage", "Percentage"),
    ):
        fig, axes = plt.subplots(nrow, 2, figsize=(w, h), squeeze=False)
        for i, s in enumerate(samples):
            for j, flavor in enumerate(FLAVORS):
                ax = axes[i][j]
                f = sub[(sub["sample"] == s) & (sub["flavor"] == flavor)]
                if kind == "percentage":
                    denom = f["Freq"].sum()
                    y = f["Freq"].to_numpy() / denom if denom else np.zeros(len(f))
                else:
                    y = f["Freq"].to_numpy(dtype=float)
                ax.bar(f["size"].to_numpy(), y, width=0.8)
                ax.set_xticks(size_breaks)
                ax.set_title(f"{s} | {flavor}", fontsize=8)
                ax.set_ylabel(ylabel, fontsize=8)
        fig.suptitle(f"{class_name} read-size distribution ({title})", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        base = (
            figures_dir / f"Figure_05{'a' if kind == 'counts' else 'b'}.{class_name}_size_barplot"
        )
        if kind == "percentage":
            base = base.with_name(base.name + ".percentage")
        save_figure(fig, base)


def run_annotation_barplots(
    tables_dir: Path,
    figures_dir: Path,
) -> None:
    """Draw all step-04 annotation-composition figures."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    all_df = _load_t2(tables_dir)
    for class_name, category in CLASS_MAP.items():
        sub = all_df[all_df["category"] == category]
        if sub.empty:
            logger.warning("step 04: no rows for category '{}', skipping", category)
            continue
        levels = READ_LEVELS if class_name == "read" else REGION_LEVELS
        _figure_04_class(class_name, sub, levels, list(sub["sample"].unique()), figures_dir)


def run_read_size_barplots(
    tables_dir: Path,
    figures_dir: Path,
) -> None:
    """Draw all step-05 read-size figures."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    all_df = _load_t2(tables_dir)
    size_sub = all_df[all_df["category"].isin([f"{c}.size" for c in _SIZE_CLASSES])]
    if size_sub.empty:
        logger.warning("step 05: no size rows found, skipping")
        return
    sizes = sorted(int(x) for x in set(size_sub["item"].astype(str)))
    min_size = min(sizes) if sizes else 0
    max_size = max(sizes) if sizes else 0
    first = (min_size // 5) * 5 + 5
    breaks = list(range(first, max_size + 1, 5))
    if sizes and breaks and breaks[0] > min_size:
        breaks = [min_size] + breaks
    breaks = np.asarray(breaks)
    for class_name in _SIZE_CLASSES:
        sub = all_df[all_df["category"] == f"{class_name}.size"]
        if sub.empty:
            continue
        samples = list(sub["sample"].unique())
        _figure_05_class(class_name, sub, samples, figures_dir, breaks)
