"""Step-03 abundance analysis (port of ``03_abundance_by_category.R``).

For each genome-feature category the per-locus read abundance is computed and
described:

* ``matmiRNA`` -- per mature miRNA locus (miRBase Name)
* ``snoRNA``   -- per snoRNA gene (refGene gene_name)
* ``tRNA``     -- per tRNA gene (gene_id)
* ``piRNA``    -- per read position (chr:start:end:strand); the piRNA feature
  set is huge so read positions are used directly

Outputs: Table_03a (per-category summary statistics incl. cross-sample tests),
Table_03b (long-format per-locus abundance) and Figure_03a/b/c.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from srna.annotation.overlap_pairs import overlap_pairs
from srna.features.feature_store import FeatureStore
from srna.plotting import CATEGORY_COLS, fig_dims, save_figure

_ABUNDANCE_CATEGORIES = ("matmiRNA", "snoRNA", "tRNA")
_PIRNA = "piRNA"

#: plot order / colours for the four abundance categories (R step-03).
CATEGORY_ORDER = ("matmiRNA", "piRNA", "tRNA", "snoRNA")

#: feature id column used as the locus name per category (mirrors the R switch).
_LOCUS_COL: dict[str, str] = {"matmiRNA": "gene_name", "snoRNA": "gene_id", "tRNA": "gene_id"}


def _locus_reads(
    per_read: pd.DataFrame,
    sample: str,
    store: FeatureStore,
) -> pd.DataFrame:
    """Compute per-locus read abundance for one sample.

    matmiRNA/snoRNA/tRNA reads are re-mapped to every overlapping feature locus
    (any overlap, strand-aware, as in the R script); piRNA reads use their read
    position as the locus.
    """
    rows: list[pd.DataFrame] = []
    for cat in _ABUNDANCE_CATEGORIES:
        g = per_read[per_read["category"] == cat]
        if g.empty:
            continue
        feat = store.table(cat)
        if feat.empty:
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
        feat = feat.copy()
        feat["feat_idx"] = np.arange(len(feat))
        pairs = overlap_pairs(reads, feat, "any")
        if len(pairs) == 0:
            continue
        ridx = pairs[:, 0]
        f_idx = pairs[:, 1]
        locus = feat[_LOCUS_COL[cat]].to_numpy()[f_idx]
        count = reads["count"].to_numpy()[ridx]
        dt = pd.DataFrame({"category": cat, "locus": locus, "count": count})
        dt = dt.groupby(["category", "locus"], sort=False)["count"].sum().reset_index()
        dt.insert(0, "sample", sample)
        rows.append(dt.rename(columns={"count": "n_reads"}))
    g = per_read[per_read["category"] == _PIRNA]
    if not g.empty:
        rows.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "category": _PIRNA,
                    "locus": (
                        g["chr"].astype(str)
                        + " "
                        + g["start"].astype(str)
                        + " "
                        + g["end"].astype(str)
                        + " "
                        + g["strand"].astype(str)
                    ),
                    "n_reads": g["count"].to_numpy(dtype=np.int64),
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["sample", "category", "locus", "n_reads"])
    return pd.concat(rows, ignore_index=True)


def _gini(x: np.ndarray) -> float:
    """Gini coefficient on per-locus read counts (sorted ascending)."""
    x = np.asarray(x, dtype=np.float64)
    if len(x) < 2 or x.sum() == 0:
        return np.nan
    x = np.sort(x)
    n = len(x)
    return float((2 * np.sum(np.arange(1, n + 1) * x)) / (n * x.sum()) - (n + 1) / n)


def _spearman(log_a: pd.Series, log_b: pd.Series) -> float:
    if len(log_a) < 2:
        return np.nan
    return float(log_a.rank().corr(log_b.rank()))


def _wilcoxon_paired(a: pd.Series, b: pd.Series) -> float:
    """Two-sided Wilcoxon signed-rank test, normal approximation (like R)."""
    d = (np.asarray(a) - np.asarray(b)).astype(np.float64)
    d = d[d != 0]
    n = len(d)
    if n == 0:
        return np.nan
    r = pd.Series(np.abs(d)).rank().to_numpy()
    w_plus = float(r[d > 0].sum())
    mean = n * (n + 1) / 4
    sd = np.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    if sd == 0:
        return np.nan
    z = (w_plus - mean - 0.5) / sd
    from math import erf

    p = 2 * (1 - 0.5 * (1 + erf(abs(z) / np.sqrt(2))))
    return float(p)


def table_03a(
    locus_tab: pd.DataFrame,
    samples: list[str],
) -> pd.DataFrame:
    """Per-category abundance summary statistics + cross-sample tests."""
    summary: list[dict[str, Any]] = []
    for cat in ("matmiRNA", "snoRNA", "tRNA", _PIRNA):
        sub = locus_tab[locus_tab["category"] == cat]
        rho: float | None = None
        p_w: float | None = None
        if len(samples) == 2:
            cs = sub.pivot_table(
                index="locus", columns="sample", values="n_reads", aggfunc="sum", fill_value=0
            )
            if set(samples) <= set(cs.columns) and len(cs) > 0:
                v1 = np.log10(cs[samples[0]].to_numpy() + 1)
                v2 = np.log10(cs[samples[1]].to_numpy() + 1)
                rho = _spearman(pd.Series(v1), pd.Series(v2))
                p_w = _wilcoxon_paired(pd.Series(v1), pd.Series(v2))
        for s in samples:
            x = sub[sub["sample"] == s].sort_values("n_reads", ascending=False)
            if x.empty:
                continue
            reads = x["n_reads"].to_numpy(dtype=np.float64)
            tot = reads.sum()
            n = len(reads)
            log10r = np.log10(reads + 1)
            top5 = float(reads[:5].sum())
            top10 = float(reads[:10].sum())
            summary.append(
                {
                    "category": cat,
                    "sample": s,
                    "n_loci": n,
                    "total_reads": int(tot),
                    "mean_log10": round(float(log10r.mean()), 3),
                    "median_log10": round(float(np.median(log10r)), 3),
                    "gini": round(_gini(reads), 3),
                    "top1_pct": round(100 * reads[0] / tot, 2),
                    "top5_pct": round(100 * top5 / tot, 2),
                    "top10_pct": round(100 * top10 / tot, 2),
                    "top_locus": str(x["locus"].iloc[0]),
                    "spearman_cross_sample": round(rho, 3) if rho is not None else np.nan,
                    "wilcoxon_p_2sample": (_format_p(p_w) if p_w is not None else pd.NA),
                }
            )
    return pd.DataFrame(summary)


def _format_p(p: float) -> str:
    if p < 2.2e-16:
        return "<2.2e-16"
    return f"{p:.3g}"


def table_03b(locus_tab: pd.DataFrame) -> pd.DataFrame:
    """Long-format per-locus abundance, ordered by category/sample/reads."""
    return locus_tab.sort_values(
        ["category", "sample", "n_reads"], ascending=[True, True, False], kind="stable"
    ).reset_index(drop=True)


def figure_03a(
    locus_tab: pd.DataFrame,
    samples: list[str],
    figures_dir: Path,
) -> None:
    """Per-locus log10 read distribution (violin + boxplot) per sample."""
    import matplotlib.pyplot as plt

    plot = locus_tab.copy()
    plot["log10_reads"] = np.log10(plot["n_reads"].to_numpy() + 1)
    order = list(CATEGORY_ORDER)
    cats = [c for c in order if c in plot["category"].unique()]
    ncol = 2 if len(samples) <= 8 else 4
    nrow = int(np.ceil(len(samples) / ncol))
    w, h = fig_dims(len(samples), ncol, per_h=4.1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        ax = axes[i // ncol][i % ncol]
        sub = plot[plot["sample"] == s]
        pos = 0
        for c in cats:
            data = sub.loc[sub["category"] == c, "log10_reads"].to_numpy()
            if len(data) == 0:
                continue
            idx = len(data) // 100
            if idx < 2:
                idx = 2
            if idx > 20:
                idx = 20
            vp = ax.violinplot([data], positions=[pos], widths=0.8, showextrema=False)
            vp["bodies"][0].set_facecolor(CATEGORY_COLS[c])
            vp["bodies"][0].set_alpha(0.7)
            q1, med, q3 = np.percentile(data, [25, 50, 75])
            ax.plot([pos - 0.1, pos + 0.1], [med, med], color="black", lw=0.8)
            ax.plot([pos - 0.07, pos + 0.07], [q1, q1], color="black", lw=0.5)
            ax.plot([pos - 0.07, pos + 0.07], [q3, q3], color="black", lw=0.5)
            ax.plot([pos, pos], [q1, q3], color="black", lw=0.5)
            pos += 1
        ax.set_title(s, fontsize=9)
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(r"$\log_{10}(reads+1)$", fontsize=9)
    for j in range(len(samples), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("per-locus read abundance by annotation category", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, figures_dir / "Figure_03a_per_locus_distribution")
    plt.close(fig)


def figure_03b(
    locus_tab: pd.DataFrame,
    samples: list[str],
    figures_dir: Path,
) -> None:
    """Rank-abundance (Whittaker) curves per category and sample."""
    import matplotlib.pyplot as plt

    ra = (
        locus_tab.groupby(["sample", "category", "locus"], sort=False)["n_reads"]
        .sum()
        .reset_index()
    )
    ra["rank"] = (
        ra.groupby(["sample", "category"], sort=False)["n_reads"]
        .rank(ascending=False, method="first")
        .to_numpy()
    )
    rng = np.random.default_rng(0)
    ra["rank"] = ra["rank"].to_numpy() + rng.normal(0, 0.08, len(ra))
    order = list(CATEGORY_ORDER)
    cats = [c for c in order if c in ra["category"].unique()]
    ncol = 2 if len(samples) <= 8 else 4
    nrow = int(np.ceil(len(samples) / ncol))
    w, h = fig_dims(len(samples), ncol, per_h=4.1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        ax = axes[i // ncol][i % ncol]
        sub = ra[ra["sample"] == s]
        for c in cats:
            x = sub[sub["category"] == c]
            if x.empty:
                continue
            x = x.sort_values("rank")
            ax.plot(x["rank"], x["n_reads"], color=CATEGORY_COLS[c], lw=0.4, alpha=0.8)
        ax.set_title(s, fontsize=9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("rank of locus (log10)", fontsize=8)
        ax.set_ylabel("reads (log10)", fontsize=8)
    for j in range(len(samples), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("rank-abundance curves per annotation category", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, figures_dir / "Figure_03b_rank_abundance")
    plt.close(fig)


def figure_03c(
    locus_tab: pd.DataFrame,
    samples: list[str],
    figures_dir: Path,
) -> None:
    """Lorenz curves for mature-miRNA loci."""
    import matplotlib.pyplot as plt

    ncol = 2 if len(samples) <= 8 else 4
    nrow = int(np.ceil(len(samples) / ncol))
    w, h = fig_dims(len(samples), ncol, per_h=4.1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    for i, s in enumerate(samples):
        ax = axes[i // ncol][i % ncol]
        x = (
            locus_tab[(locus_tab["sample"] == s) & (locus_tab["category"] == "matmiRNA")]["n_reads"]
            .sort_values(ascending=False)
            .to_numpy(dtype=np.float64)
        )
        n = len(x)
        ax.plot([0, 1], [0, 1], linestyle="--", color="gray", lw=1)
        if n > 0:
            cum = np.concatenate([[0], np.cumsum(x) / x.sum()])
            i_axis = np.concatenate([[0], np.arange(1, n + 1) / n])
            ax.plot(i_axis, cum, color="steelblue", lw=0.6)
        ax.set_title(s, fontsize=9)
        ax.set_xlabel("cumulative fraction of loci", fontsize=8)
        ax.set_ylabel("cumulative fraction of reads", fontsize=8)
        ax.set_aspect("equal")
    for j in range(len(samples), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Lorenz curves - mature miRNA loci", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, figures_dir / "Figure_03c_lorenz_matmiRNA")
    plt.close(fig)


def run_step03(
    per_reads: dict[str, pd.DataFrame],
    store: FeatureStore,
    samples: list[str],
    out_dir: Path,
) -> None:
    """Run the full step-03 for the given samples."""
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    locus_tab = pd.concat(
        [_locus_reads(per_reads[s], s, store) for s in samples], ignore_index=True
    )
    logger.info("step 03: {} per-locus abundance rows", len(locus_tab))

    t3a = table_03a(locus_tab, samples)
    t3a.to_csv(tables_dir / "Table_03a_category_abundance_summary.csv", index=False)
    table_03b(locus_tab).to_csv(tables_dir / "Table_03b_per_locus_abundance.csv", index=False)

    figure_03a(locus_tab, samples, figures_dir)
    figure_03b(locus_tab, samples, figures_dir)
    figure_03c(locus_tab, samples, figures_dir)
