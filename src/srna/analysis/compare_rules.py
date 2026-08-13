"""Step s01: three-way comparison of overlap-rule strategies.

Port of ``s01_compare_overlap_rules.R``.

Only meaningful in comparison mode. The three strategies (fully-contained,
union, any) differ only in the sense overlap rule; this step compares their
step-02 outputs: annotation composition, per-category read sizes, category
totals, read movement between strategies, mature-miRNA expression, and strand
specificity of the ``any`` strategy.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from srna.features.store import FeatureStore
from srna.plotting import STRATEGY_COLS, fig_dims, save_figure, small_fonts

STRATEGIES = ("fully-contained", "union", "any")
_STRATEGY_DIRS = {"fully-contained": "fully_contained", "union": "union", "any": "any"}

_ITEM_FACTOR = [
    "matmiRNA",
    "piRNA",
    "snoRNA",
    "tRNA",
    "RM",
    "gene.exon",
    "gene.intron",
    "lincRNA.exon",
    "AS_matmiRNA",
    "AS_piRNA",
    "AS_snoRNA",
    "AS_tRNA",
    "AS_RM",
    "AS_gene.exon",
    "AS_gene.intron",
    "AS_lincRNA.exon",
    "other",
]

_SIZE_CATEGORIES = ["matmiRNA", "piRNA", "snoRNA", "tRNA"]


def _strategy_dir(base: Path, strategy: str) -> Path:
    return base / _STRATEGY_DIRS[strategy]


def _composition(
    base: Path,
) -> pd.DataFrame:
    """read.annotation composition for every strategy (d.all)."""
    frames: list[pd.DataFrame] = []
    for st in STRATEGIES:
        tables = _strategy_dir(base, st) / "tables"
        t2a = pd.read_csv(tables / "Table_02a_annotation_count_unique_reads.csv")
        t2b = pd.read_csv(tables / "Table_02b_annotation_count_all_reads.csv")
        u = t2a[t2a["category"] == "read.annotation"][["sample", "item", "Freq"]].rename(
            columns={"Freq": "n_unique"}
        )
        r = t2b[t2b["category"] == "read.annotation"][["sample", "item", "Freq"]].rename(
            columns={"Freq": "n_reads"}
        )
        m = u.merge(r, on=["sample", "item"], how="outer").fillna(0)
        m["pct_unique"] = 100 * m["n_unique"] / m.groupby("sample")["n_unique"].transform("sum")
        m["pct_reads"] = 100 * m["n_reads"] / m.groupby("sample")["n_reads"].transform("sum")
        m["strategy"] = st
        frames.append(m)
    d = pd.concat(frames, ignore_index=True)
    d["item"] = d["item"].astype(str).str.replace(r"^AS\.", "AS_", regex=True)
    d["item"] = d["item"].str.replace(r"^refGene\.NM\.", "gene.", regex=True)
    return d


def table_s01a(d: pd.DataFrame) -> pd.DataFrame:
    """Ordered overlap-rule composition table."""
    return d.sort_values(["strategy", "pct_reads"], ascending=[True, False]).reset_index(drop=True)


def table_s01b(d: pd.DataFrame) -> pd.DataFrame:
    """Per-category read totals (pct) with any/union deltas vs fully-contained."""
    cat = d.pivot_table(
        index=["sample", "item"],
        columns="strategy",
        values="pct_reads",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    cat["dany"] = cat["any"] - cat["fully-contained"]
    cat["dunion"] = cat["union"] - cat["fully-contained"]
    return cat.sort_values(["sample", "fully-contained"], ascending=[True, False]).reset_index(
        drop=True
    )


def _read_movement(
    base: Path,
    per_reads: dict[str, dict[str, pd.DataFrame]],
    to_col: str,
) -> pd.DataFrame:
    """Cross-tabulate read movement fully-contained -> *to_col* strategy.

    The movement is keyed on read position (poskey) with the *to_col* strategy
    table as the base (its counts are used); the fully-contained category is
    carried along as ``typeB``. Reads present in only one strategy get the
    missing category as ``other``.
    """
    rows: list[pd.DataFrame] = []
    for _s, tables in per_reads.items():
        fc = tables["fully-contained"]
        to = tables[to_col]
        fc["poskey"] = (
            fc["chr"].astype(str)
            + ":"
            + fc["start"].astype(str)
            + ":"
            + fc["end"].astype(str)
            + ":"
            + fc["strand"].astype(str)
        )
        to["poskey"] = (
            to["chr"].astype(str)
            + ":"
            + to["start"].astype(str)
            + ":"
            + to["end"].astype(str)
            + ":"
            + to["strand"].astype(str)
        )
        m = to[["poskey", "count", "category"]].merge(
            fc[["poskey", "category"]], on="poskey", how="left", suffixes=("", "_B")
        )
        m["typeB"] = m["category_B"].fillna("other")
        rows.append(m[["typeB", "category", "count"]])
    mov = pd.concat(rows, ignore_index=True)
    agg = mov.groupby(["typeB", "category"], sort=False)["count"].sum().reset_index()
    return (
        agg.sort_values("count", ascending=False)
        .rename(columns={"category": to_col})
        .reset_index(drop=True)
    )


def table_s01e(
    base: Path,
    store: FeatureStore,
    per_reads: dict[str, dict[str, pd.DataFrame]],
) -> pd.DataFrame:
    """Mature-miRNA expression per strategy (any-overlap remapping to loci)."""
    feat = store.table("matmiRNA").drop_duplicates(subset=["chrom", "start", "end", "strand"])
    if feat.empty:
        return pd.DataFrame(columns=["sample", "strategy", "name", "n_unique", "n_reads"])
    f = feat.copy()
    f["feat_idx"] = np.arange(len(f))
    records: list[pd.DataFrame] = []
    for s, tables in per_reads.items():
        for st in STRATEGIES:
            g = tables[st]
            g = g[g["category"] == "matmiRNA"]
            if g.empty:
                continue
            reads = pd.DataFrame(
                {
                    "chrom": g["chr"].astype(str).to_numpy(),
                    "start": (g["start"] - 1).to_numpy(),
                    "end": g["end"].to_numpy(),
                    "strand": g["strand"].astype(str).to_numpy(),
                    "read_idx": np.arange(len(g)),
                    "count": g["count"].to_numpy(dtype=np.int64),
                }
            )
            from srna.annotation.pairs import overlap_pairs

            pairs = overlap_pairs(reads, f, "any")
            if len(pairs) == 0:
                continue
            hit = pd.DataFrame(
                {
                    "sample": s,
                    "strategy": st,
                    "name": feat["gene_name"].to_numpy()[pairs[:, 1]],
                    "read_idx": pairs[:, 0],
                    "count": reads["count"].to_numpy()[pairs[:, 0]],
                }
            )
            agg = (
                hit.groupby(["sample", "strategy", "name"], sort=False)
                .agg(n_unique=("read_idx", "nunique"), n_reads=("count", "sum"))
                .reset_index()
            )
            records.append(agg)
    if not records:
        return pd.DataFrame(columns=["sample", "strategy", "name", "n_unique", "n_reads"])
    expr = pd.concat(records, ignore_index=True)
    return expr.sort_values(
        ["sample", "strategy", "n_reads"], ascending=[True, True, False]
    ).reset_index(drop=True)


def _load_per_reads(base: Path, samples: list[str]) -> dict[str, dict[str, pd.DataFrame]]:
    """Load per-read tables for every strategy (missing -> empty)."""
    out: dict[str, dict[str, pd.DataFrame]] = {}
    for s in samples:
        out[s] = {}
        for st in STRATEGIES:
            path = _strategy_dir(base, st) / "tables" / f"Table_02_{s}_unique_reads_annotation.csv"
            if path.exists():
                out[s][st] = pd.read_csv(path)
            else:
                logger.warning("s01: missing per-read table {}", path)
                out[s][st] = pd.DataFrame(
                    columns=[
                        "sample",
                        "read_id",
                        "chr",
                        "start",
                        "end",
                        "strand",
                        "size",
                        "count",
                        "cpm",
                        "category",
                        "feature_id",
                        "gene_context",
                        "n_features",
                    ]
                )
    return out


def figure_s01a(d: pd.DataFrame, samples: list[str], figures_dir: Path) -> None:
    """Annotation-composition dodge barplot, faceted by sample."""
    import matplotlib.pyplot as plt

    items = [it for it in _ITEM_FACTOR if it in set(d["item"])]
    ncol = 1
    nrow = len(samples)
    w, h = fig_dims(len(samples), 1, per_h=8.2)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        ax = axes[i][0]
        sub = d[d["sample"] == s]
        x = np.arange(len(items))
        width = 0.7 / 3
        for k, st in enumerate(STRATEGIES):
            vals = []
            for it in items:
                row = sub[(sub["strategy"] == st) & (sub["item"] == it)]
                vals.append(row["pct_reads"].sum() if len(row) else 0.0)
            ax.bar(x + (k - 1) * width, vals, width=width, color=STRATEGY_COLS[st], label=st)
        ax.set_xticks(x)
        ax.set_xticklabels(items, rotation=45, ha="right")
        ax.set_title(s, fontsize=8)
        ax.set_ylabel("% of reads", fontsize=8)
        small_fonts(ax)
    axes[0][0].legend(fontsize=7, ncol=3, loc="upper right")
    for j in range(len(samples), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("read.annotation composition - all reads (overlap-rule strategies)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, figures_dir / "Figure_s01a_overlap_rules_composition")


def figure_s01b(base: Path, samples: list[str], figures_dir: Path) -> None:
    """Per-category read-size distributions, faceted by category x sample."""
    import matplotlib.pyplot as plt

    frames: list[pd.DataFrame] = []
    for st in STRATEGIES:
        t2b = pd.read_csv(
            _strategy_dir(base, st) / "tables" / "Table_02b_annotation_count_all_reads.csv"
        )
        sub = t2b[t2b["category"].isin([f"{c}.size" for c in _SIZE_CATEGORIES])]
        if sub.empty:
            continue
        sub = sub.copy()
        sub["strategy"] = st
        frames.append(sub)
    if not frames:
        return
    size_all = pd.concat(frames, ignore_index=True)
    size_all["size"] = size_all["item"].astype(int)
    size_all["pct"] = (
        100
        * size_all["Freq"]
        / size_all.groupby(["sample", "category", "strategy"], sort=False)["Freq"].transform("sum")
    )
    sizes = sorted(set(size_all["size"]))
    cats = [c for c in _SIZE_CATEGORIES if c in set(size_all["category"])]
    ncol = len(samples)
    nrow = len(cats)
    w, h = fig_dims(4 * len(samples), len(samples), per_h=4.1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, cat in enumerate(cats):
        for j, s in enumerate(samples):
            ax = axes[i][j]
            sub = size_all[(size_all["category"] == f"{cat}.size") & (size_all["sample"] == s)]
            x = np.arange(len(sizes))
            width = 0.8 / 3
            for k, st in enumerate(STRATEGIES):
                vals = []
                for sz in sizes:
                    row = sub[(sub["strategy"] == st) & (sub["size"] == sz)]
                    vals.append(row["pct"].sum() if len(row) else 0.0)
                ax.bar(x + (k - 1) * width, vals, width=width, color=STRATEGY_COLS[st])
            ax.set_xticks(x)
            ax.set_xticklabels(sizes, rotation=45, ha="right", fontsize=7)
            ax.set_title(f"{cat} | {s}", fontsize=8)
            ax.set_ylabel("%", fontsize=8)
            small_fonts(ax)
    for j in range(len(cats) * len(samples), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("per-category read-size distribution (all reads)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, figures_dir / "Figure_s01b_overlap_rules_category_size")


def run_s01(
    base: Path,
    store: FeatureStore,
    samples: list[str],
    figures_dir: Path,
    tables_dir: Path,
) -> None:
    """Run the full step s01 for a comparison layout."""
    for st in STRATEGIES:
        if not _strategy_dir(base, st).exists():
            raise FileNotFoundError(
                f"strategy dir missing: {_strategy_dir(base, st)} (run steps 02-06 first)"
            )

    d = _composition(base)
    table_s01a(d).to_csv(tables_dir / "Table_s01a_overlap_rule_composition.csv", index=False)
    table_s01b(d).to_csv(tables_dir / "Table_s01b_overlap_rule_category_totals.csv", index=False)

    per_reads = _load_per_reads(base, samples)
    mov_ba = _read_movement(base, per_reads, "typeA")
    mov_ba.to_csv(tables_dir / "Table_s01c_read_movement_contained_vs_any.csv", index=False)
    mov_bc = _read_movement(base, per_reads, "typeC")
    mov_bc.to_csv(tables_dir / "Table_s01d_read_movement_contained_vs_union.csv", index=False)

    expr = table_s01e(base, store, per_reads)
    expr.to_csv(tables_dir / "Table_s01e_mature_miRNA_expression_strategies.csv", index=False)

    figure_s01a(d, samples, figures_dir)
    figure_s01b(base, samples, figures_dir)
