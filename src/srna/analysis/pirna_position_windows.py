"""Step-06 analysis: piRNA reads on tRNA/snoRNA gene windows.

Port of ``06_tRNA_snoRNA_position.windows.R``.

piRNA-annotated reads (sense to a piRNA locus) may be degradation products of
tRNA/snoRNA genes. For each tRNA/snoRNA gene a 21 bp window is slid in 1 bp
steps along the gene body, the mean read count per window position is plotted
per sample (sense and antisense panels), and a summary reports the percentage
of piRNA reads overlapping snoRNA/tRNA genes (sense/antisense) or neither.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from srna.annotation.overlap_pairs import overlap_pairs
from srna.features.feature_store import FeatureStore
from srna.plotting import fig_dims, save_figure, small_fonts

_OVERLAP_LEVELS = [
    "snoRNA_sense",
    "snoRNA_antisense",
    "tRNA_sense",
    "tRNA_antisense",
    "none",
]

_OVERLAP_COLS = {
    "snoRNA_sense": "#1f78b4",
    "snoRNA_antisense": "#a6cee3",
    "tRNA_sense": "#e31a1c",
    "tRNA_antisense": "#fb9a99",
    "none": "grey",
}


def sliding_windows(features: pd.DataFrame) -> pd.DataFrame:
    """Build 21 bp sliding windows (1 bp step) along each gene body.

    Position is counted 5' -> 3' (position 1 = the 5' end for both strands).
    Genes of length <= 20 bp produce no windows. The output is 0-based.
    """
    records: list[dict] = []
    for row in features.itertuples(index=False):
        length = int(row.end) - int(row.start)
        if length <= 20:
            continue
        if row.strand == "+":
            starts = [row.start + k for k in range(length - 20)]
        else:
            starts = [row.end - 21 - k for k in range(length - 20)]
        for k, start in enumerate(starts):
            records.append(
                {
                    "chrom": row.chrom,
                    "start": start,
                    "end": start + 21,
                    "strand": row.strand,
                    "position": k + 1,
                    "gene_id": row.gene_id,
                }
            )
    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "position", "gene_id"])
    return pd.DataFrame(records)


def _window_counts(
    windows: pd.DataFrame,
    reads: pd.DataFrame,
    ignore_strand: bool,
) -> np.ndarray:
    """Sum read counts overlapping each window (any overlap)."""
    out = np.zeros(len(windows), dtype=np.int64)
    if reads.empty or windows.empty:
        return out
    win = windows.copy()
    win["read_idx"] = np.arange(len(win))
    r = reads.copy()
    r["feat_idx"] = np.arange(len(r))
    pairs = overlap_pairs(win, r, "any", ignore_strand=ignore_strand)
    if len(pairs) == 0:
        return out
    win_idx = pairs[:, 0]
    read_idx = pairs[:, 1]
    counts = r["count"].to_numpy()[read_idx]
    np.add.at(out, win_idx, counts)
    return out


def _pirna_reads(per_read: pd.DataFrame) -> pd.DataFrame:
    """Extract sense-piRNA reads with 0-based coordinates + count column."""
    g = per_read[per_read["category"] == "piRNA"]
    if g.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "count"])
    return pd.DataFrame(
        {
            "chrom": g["chr"].astype(str).to_numpy(),
            "start": (g["start"] - 1).to_numpy(),
            "end": g["end"].to_numpy(),
            "strand": g["strand"].astype(str).to_numpy(),
            "count": g["count"].to_numpy(dtype=np.int64),
        }
    )


def _overlap_table(per_read: pd.DataFrame, store: FeatureStore) -> pd.DataFrame:
    """Categorise every sense-piRNA read by snoRNA/tRNA sense/antisense/neither."""
    g = _pirna_reads(per_read)
    if g.empty:
        return pd.DataFrame(
            columns=["sample", "category", "n_unique", "n_reads", "pct_unique", "pct_reads"]
        )
    reads = g.copy()
    reads["read_idx"] = np.arange(len(reads))
    labels = ["none"] * len(reads)

    def assign(feature_name: str, sense_label: str, anti_label: str) -> None:
        feat = store.table(feature_name)
        if feat.empty:
            return
        f = feat.copy()
        f["feat_idx"] = np.arange(len(f))
        sense = overlap_pairs(reads, f, "any", ignore_strand=False)
        for ridx in sense[:, 0]:
            labels[ridx] = sense_label
        remain = [i for i, lab in enumerate(labels) if lab == "none"]
        if remain:
            sub = reads.iloc[remain].copy()
            sub["read_idx"] = np.arange(len(sub))
            anti = overlap_pairs(sub, f, "any", ignore_strand=True)
            for pos in anti[:, 0]:
                labels[remain[pos]] = anti_label

    assign("snoRNA", "snoRNA_sense", "snoRNA_antisense")
    assign("tRNA", "tRNA_sense", "tRNA_antisense")

    tab = pd.DataFrame(
        {
            "category": labels,
            "count": g["count"].to_numpy(dtype=np.int64),
        }
    )
    tab = (
        tab.groupby("category", sort=False)
        .agg(n_unique=("count", "size"), n_reads=("count", "sum"))
        .reset_index()
    )
    tab["sample"] = ""
    tab["pct_unique"] = (100 * tab["n_unique"] / tab["n_unique"].sum()).round(2)
    tab["pct_reads"] = (100 * tab["n_reads"] / tab["n_reads"].sum()).round(2)
    return tab[["sample", "category", "n_unique", "n_reads", "pct_unique", "pct_reads"]]


def _window_profile(
    windows: pd.DataFrame,
    pirna: pd.DataFrame,
    sample: str,
) -> pd.DataFrame:
    """sense/antisense counts per window for one sample, both read flavors."""
    base = windows[["chrom", "start", "end", "strand", "position", "gene_id"]].copy()
    out: list[pd.DataFrame] = []
    for flavor, counts in (("all reads", pirna["count"]), ("unique reads", np.ones(len(pirna)))):
        df = base.copy()
        pir = pirna.copy()
        pir["count"] = counts
        df["sense"] = _window_counts(windows, pir, ignore_strand=False)
        flipped = windows.copy()
        flipped["strand"] = np.where(flipped["strand"] == "+", "-", "+")
        df["antisense"] = _window_counts(flipped, pir, ignore_strand=False)
        df["sample"] = sample
        df["flavor"] = flavor
        out.append(df)
    return pd.concat(out, ignore_index=True)


def run_pirna_windows(
    per_reads: dict[str, pd.DataFrame],
    store: FeatureStore,
    samples: list[str],
    out_dir: Path,
) -> None:
    """Run the full step-06 analysis."""
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    rdata_dir = out_dir / "rdata"
    for d in (tables_dir, figures_dir, rdata_dir):
        d.mkdir(parents=True, exist_ok=True)

    trna_windows = sliding_windows(store.table("tRNA"))
    snorna_feat = store.table("snoRNA")
    snorna_feat = snorna_feat[snorna_feat["end"] - snorna_feat["start"] > 20]
    snorna_windows = sliding_windows(snorna_feat)
    logger.info(
        "step 06: tRNA windows={}, snoRNA windows={}", len(trna_windows), len(snorna_windows)
    )

    trna_profiles: list[pd.DataFrame] = []
    snorna_profiles: list[pd.DataFrame] = []
    overlap_tabs: list[pd.DataFrame] = []
    for s in samples:
        pirna = _pirna_reads(per_reads[s])
        logger.info("step 06: sample '{}' has {} sense-piRNA reads", s, len(pirna))
        if len(trna_windows):
            trna_profiles.append(_window_profile(trna_windows, pirna, s))
        if len(snorna_windows):
            snorna_profiles.append(_window_profile(snorna_windows, pirna, s))
        ot = _overlap_table(per_reads[s], store)
        ot["sample"] = s
        overlap_tabs.append(ot)

    if trna_profiles:
        pd.concat(trna_profiles, ignore_index=True).to_parquet(
            rdata_dir / "piRNA_on_tRNA.20bp.dis.all.parquet"
        )
    if snorna_profiles:
        pd.concat(snorna_profiles, ignore_index=True).to_parquet(
            rdata_dir / "piRNA_on_snoRNA.20bp.dis.all.parquet"
        )

    ov = pd.concat(overlap_tabs, ignore_index=True)
    ov.to_csv(tables_dir / "Table_06_piRNA_overlap_summary.csv", index=False)

    if trna_profiles:
        _position_figure(pd.concat(trna_profiles, ignore_index=True), samples, "tRNA", figures_dir)
    if snorna_profiles:
        _position_figure(
            pd.concat(snorna_profiles, ignore_index=True), samples, "snoRNA", figures_dir
        )
    _overlap_summary_figure(ov, samples, figures_dir)


def _position_figure(
    profile: pd.DataFrame,
    samples: list[str],
    gene_name: str,
    figures_dir: Path,
) -> None:
    """Barplots of mean sense/antisense counts per window position."""
    import matplotlib.pyplot as plt

    for panel, col in (("", "sense"), ("_AS", "antisense")):
        agg = (
            profile.groupby(["position", "sample", "flavor"], sort=False)[col].mean().reset_index()
        )
        nrow = len(samples)
        w, h = fig_dims(len(samples), 2, per_h=4.1)
        fig, axes = plt.subplots(nrow, 2, figsize=(w, h), squeeze=False)
        for i, s in enumerate(samples):
            for j, flavor in enumerate(["unique reads", "all reads"]):
                ax = axes[i][j]
                sub = agg[(agg["sample"] == s) & (agg["flavor"] == flavor)]
                ax.bar(sub["position"], sub[col], width=1.0)
                ax.set_title(f"{s} | {flavor}", fontsize=8)
                ax.set_ylabel("Mean count", fontsize=8)
                small_fonts(ax)
        for j in range(len(samples), nrow * 2):
            axes[j // 2][j % 2].axis("off")
        orient = "sense" if col == "sense" else "antisense"
        fig.suptitle(
            f"{gene_name} - piRNA reads, {orient} position distribution (20 bp windows)",
            fontsize=10,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        save_figure(fig, figures_dir / f"Figure_06.piRNA_vs_{gene_name}_pos_barplot{panel}")


def _overlap_summary_figure(
    ov: pd.DataFrame,
    samples: list[str],
    figures_dir: Path,
) -> None:
    """Stacked percentage bar of piRNA overlap categories (samples x flavor)."""
    import matplotlib.pyplot as plt

    long = pd.melt(
        ov,
        id_vars=["sample", "category"],
        value_vars=["n_unique", "n_reads"],
        var_name="flavor",
        value_name="n",
    )
    long["flavor"] = long["flavor"].map({"n_unique": "unique reads", "n_reads": "all reads"})
    long["pct"] = 100 * long["n"] / long.groupby(["sample", "flavor"])["n"].transform("sum")
    long["category"] = pd.Categorical(long["category"], categories=_OVERLAP_LEVELS, ordered=True)
    flavors = ["unique reads", "all reads"]
    w, h = fig_dims(len(samples), 2, per_h=2.7)
    fig, axes = plt.subplots(1, 2, figsize=(w, h), squeeze=True)
    for j, flavor in enumerate(flavors):
        ax = axes[j]
        sub = long[long["flavor"] == flavor]
        bottom = np.zeros(len(samples))
        for cat in _OVERLAP_LEVELS:
            vals = np.zeros(len(samples))
            for i, s in enumerate(samples):
                row = sub[(sub["sample"] == s) & (sub["category"] == cat)]
                vals[i] = row["pct"].sum() if len(row) else 0.0
            ax.bar(np.arange(len(samples)), vals, bottom=bottom, color=_OVERLAP_COLS[cat])
            bottom += vals
        ax.set_xticks(np.arange(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("% of piRNA reads", fontsize=8)
        ax.set_title(flavor, fontsize=8)
        small_fonts(ax)
    handles = [plt.Rectangle((0, 0), 1, 1, color=_OVERLAP_COLS[c]) for c in _OVERLAP_LEVELS]
    fig.legend(
        handles, _OVERLAP_LEVELS, title="piRNA overlap", loc="upper center", ncol=5, fontsize=7
    )
    fig.suptitle("Percentage of piRNA reads overlapping snoRNA/tRNA genes", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save_figure(fig, figures_dir / "Figure_06.piRNA_overlap_summary")
