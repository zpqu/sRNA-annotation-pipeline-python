"""Pipeline orchestration (port of ``run_smallRNA_annotation.sh``).

Runs the whole small-RNA annotation pipeline for a :class:`PipelineConfig`:

* step 00 -- build/cache the feature DB
* step 01 -- collapse BAMs to non-redundant reads; write Table_01a-d + Figure_01
* steps 02-06 -- annotate reads (per strategy) and produce annotation/abundance
  tables and figures
* step s01 -- overlap-rule comparison (comparison mode only)
* step 10 -- end-of-pipeline summary report
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from loguru import logger

from srna.config import PipelineConfig, Strategy
from srna.features.feature_store import FeatureStore
from srna.paths import PathResolver
from srna.reads.collapse import collapse_bam
from srna.reads.read_size_figure import read_size_overview_figure
from srna.reads.summaries import (
    build_table_01a,
    build_table_01b,
    build_table_01c,
    build_table_01d,
)

_STEP02_COLUMNS = [
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


@dataclass
class CollapseResult:
    """Step-01 outputs: per-sample unique reads and alignment totals."""

    reads: dict[str, pd.DataFrame] = field(default_factory=dict)
    total_reads: dict[str, int] = field(default_factory=dict)


def run_build_feature_db(store: FeatureStore) -> None:
    """Build the feature DB (no-op when cached and not forced)."""
    store.build()
    logger.info("step 00: feature DB ready at {}", store.cache_dir)


def run_read_collapse(resolver: PathResolver, samples: list[str]) -> CollapseResult:
    """Collapse BAMs and write shared read tables + Figure_01."""
    tables_dir = resolver.output_base / "tables"
    rdata_dir = resolver.output_base / "rdata"
    figures_dir = resolver.output_base / "figures"
    for d in (tables_dir, rdata_dir, figures_dir):
        d.mkdir(parents=True, exist_ok=True)

    result = CollapseResult()
    t1a: list[dict] = []
    t1b: list[pd.DataFrame] = []
    t1c: list[pd.DataFrame] = []
    t1d: list[pd.DataFrame] = []
    for sample in samples:
        bam = resolver.bam_file(sample)
        logger.info("step 01: collapsing {bam}", bam=bam.name)
        unique, total = collapse_bam(bam)
        result.reads[sample] = unique
        result.total_reads[sample] = total
        unique.to_parquet(rdata_dir / f"{sample}.bam.unique.parquet")
        t1a.append(build_table_01a(unique, total, sample))
        t1b.append(build_table_01b(unique, sample))
        t1c.append(build_table_01c(unique, sample))
        t1d.append(build_table_01d(unique, sample))
        logger.info("step 01: sample '{}' total={} unique={}", sample, total, len(unique))

    pd.DataFrame(t1a).to_csv(tables_dir / "Table_01a_sample_summary.csv", index=False)
    pd.concat(t1b, ignore_index=True).to_csv(
        tables_dir / "Table_01b_read_size_distribution.csv", index=False
    )
    pd.concat(t1c, ignore_index=True).to_csv(
        tables_dir / "Table_01c_read_count_distribution.csv", index=False
    )
    pd.concat(t1d, ignore_index=True).to_csv(
        tables_dir / "Table_01d_read_size_vs_count.csv", index=False
    )
    read_size_overview_figure(result.reads, figures_dir / "Figure_01.read_size_vs_count")
    return result


def run_annotation_steps(
    config: PipelineConfig,
    resolver: PathResolver,
    store: FeatureStore,
    step01: CollapseResult,
) -> None:
    """Run steps 02-06 for one concrete strategy into its output directory."""
    from srna.analysis.annotation_and_size_barplots import (
        run_annotation_barplots,
        run_read_size_barplots,
    )
    from srna.analysis.locus_abundance import run_abundance
    from srna.analysis.pirna_position_windows import run_pirna_windows
    from srna.annotation.annotator import annotate_sample

    if resolver.config.strategy.is_comparison:
        sdir = resolver.strategy_dir(config.strategy)
    else:
        sdir = resolver.output_base
    for d in ("tables", "figures", "rdata"):
        (sdir / d).mkdir(parents=True, exist_ok=True)

    samples = list(step01.reads)
    per_reads: dict[str, pd.DataFrame] = {}
    for sample in samples:
        result = annotate_sample(step01.reads[sample], store, config.strategy, sample)
        per_reads[sample] = result.per_read
        out = result.per_read.copy()
        out.to_csv(
            sdir / "tables" / f"Table_02_{sample}_unique_reads_annotation.csv",
            index=False,
            na_rep="NA",
        )
        out.to_parquet(sdir / "rdata" / f"{sample}.bam.annotated.parquet")

        rows_u = [
            {"sample": sample, "category": c, "item": i, "Freq": n}
            for c, i, n in result.counts_unique
        ]
        rows_r = [
            {"sample": sample, "category": c, "item": i, "Freq": n}
            for c, i, n in result.counts_reads
        ]
        if rows_u:
            pd.DataFrame(rows_u).to_csv(
                sdir / "tables" / "Table_02a_annotation_count_unique_reads.csv",
                index=False,
                mode="a" if sample != samples[0] else "w",
                header=(sample == samples[0]),
            )
        if rows_r:
            pd.DataFrame(rows_r).to_csv(
                sdir / "tables" / "Table_02b_annotation_count_all_reads.csv",
                index=False,
                mode="a" if sample != samples[0] else "w",
                header=(sample == samples[0]),
            )
        logger.info(
            "step 02: sample '{}' annotated ({} unique reads)", sample, len(per_reads[sample])
        )

    run_abundance(per_reads, store, samples, sdir)
    run_annotation_barplots(sdir / "tables", sdir / "figures")
    run_read_size_barplots(sdir / "tables", sdir / "figures")
    run_pirna_windows(per_reads, store, samples, sdir)


def run_pipeline(
    config: PipelineConfig,
    root: Path,
) -> None:
    """Run the full pipeline for *config* with project root *root*."""
    resolver = PathResolver(root, config)
    resolver.require_raw_db_files()
    resolver.output_base.mkdir(parents=True, exist_ok=True)

    store = FeatureStore(resolver.db_raw_dir, resolver.db_cache_dir, config)
    run_build_feature_db(store)

    samples = resolver.samples()
    if not samples:
        raise FileNotFoundError(f"no *.bam files found in {resolver.bam_dir}")
    logger.info("pipeline: samples = {}", samples)

    step01 = run_read_collapse(resolver, samples)

    if config.strategy.is_comparison:
        for sub in (Strategy.FULLY_CONTAINED, Strategy.UNION, Strategy.ANY):
            logger.info("pipeline: running steps 02-06 for strategy '{}'", sub.value)
            run_annotation_steps(
                PipelineConfig.for_substrategy(config, sub), resolver, store, step01
            )
        from srna.analysis.overlap_rule_comparison import run_overlap_rule_comparison

        run_overlap_rule_comparison(
            resolver.output_base,
            store,
            samples,
            resolver.output_base / "figures",
            resolver.output_base / "tables",
        )
    else:
        run_annotation_steps(config, resolver, store, step01)

    from srna.summary.pipeline_report import run_pipeline_summary

    run_pipeline_summary(resolver, samples, store)
