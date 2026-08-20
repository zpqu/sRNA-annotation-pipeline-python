# srna-annotation

A production-quality, end-to-end small-RNA sequencing annotation pipeline in
Python. It collapses aligned reads (BAM), annotates each unique read against
small-RNA and genomic feature sets with a configurable sense-overlap rule, and
produces abundance analyses, publication-ready figures and a self-contained
report.

## Overview

The pipeline processes small-RNA sequencing data from aligned reads to an
annotated, quantified, and visually summarised dataset:

- Collapse BAM alignments into non-redundant reads with per-sample size and
  abundance summaries.
- Annotate each read by priority against mature miRNA, piRNA, snoRNA, tRNA,
  repeat-masker and refGene (exon / intron / lincRNA) features, on both strands,
  plus antisense and intergenic categories.
- Quantify per-locus abundance with Gini, rank-abundance and Lorenz analyses.
- Draw annotation-composition and read-size figures, piRNA degradation profiles
  over tRNA/snoRNA gene windows, and — in comparison mode — a three-way
  comparison of sense-overlap rules.
- Finish with an end-of-pipeline markdown report and sanity checks.

## Features

- Non-redundant read collapsing with abundance (count, CPM) weighting.
- Priority-ordered annotation with gene-context region assignment
  (CDS / 5UTR / 3UTR / intron / ±1 kb / repeat-masker) and antisense detection.
- Configurable sense-overlap strategy: `fully-contained`, `union`, `any`, or a
  `comparison` of all three.
- Streaming, index-aware BAM processing and cached feature databases (parquet)
  so reruns skip work.
- PDF + PNG figures and CSV + parquet tables at every step.
- End-of-pipeline markdown report with sanity checks.

## Requirements

- Python ≥ 3.12
- A Linux/macOS environment with the BAM files and raw feature annotations
  described below.

## Installation

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

The only executable is the `srna` CLI entry point.

## Input data layout

The pipeline reads everything from a project root (default: the current
directory):

```
data/
  bam/        aligned reads, one indexed BAM per sample (<sample>.bwa.bam + .bai)
  DB/         raw feature annotations
  cache/      per-genome feature DB caches (built and reused automatically)
```

`data/DB/` is expected to contain the feature files for the reference genome:

| File | Content |
|---|---|
| `miRNA.gff3` | mature miRNA coordinates |
| `piRBase.bed` | piRNA loci |
| `piRNAdb.gtf` | piRNA loci (additional source) |
| `refGene.gtf` | refGene transcripts (exons, CDS, UTRs, introns, ±1 kb) |
| `RM.bed` | RepeatMasker repeats |
| `tRNAs.bed` | tRNA genes |

Each BAM in `data/bam/` is treated as one sample; its basename (without
`.bwa.bam`/`.bam`) becomes the sample id used in outputs.

## Usage

```bash
.venv/bin/srna --root . --genome mm39 --strategy fully-contained
```

### Command-line options

| Option | Default | Description |
|---|---|---|
| `--root PATH` | current dir | project root containing `data/` and `output/` |
| `--genome ID` | `mm39` | assembly id; raw files read from `data/DB/`, caches written to `data/cache/<ID>/` |
| `--strategy RULE` | `fully-contained` | sense-overlap rule (see below) |
| `--chr-style STYLE` | `chr` | chromosome naming convention, `chr` or `none` |
| `--force-rebuild-db` | off | rebuild the feature DB even when caches are valid |
| `--log-dir PATH` | `<root>/logs` | directory for the pipeline log |

The pipeline is also configurable through environment variables
(`SMALLRNA_GENOME`, `SMALLRNA_STRATEGY`, `SMALLRNA_CHR_STYLE`,
`SMALLRNA_FORCE_REBUILD_DB`).

### Sense-overlap strategies

| Strategy | Rule |
|---|---|
| `fully-contained` | a read is annotated to a feature only if the read is fully contained in it |
| `union` | read within feature **or** feature within read (captures long reads spanning a small feature) |
| `any` | any overlap (≥ 1 bp) between read and feature |
| `comparison` | runs the three rules above and then the overlap-rule comparison analysis |

## Pipeline steps

| Step | Purpose |
|---|---|
| 00 | Build the genomic-feature DB from `data/DB/` and cache to `data/cache/<GENOME>/` |
| 01 | Collapse BAMs to non-redundant reads; read-size / count tables and overview figure |
| 02 | Annotate non-redundant reads by priority (sense + antisense + other) |
| 03 | Per-locus abundance, Gini, rank-abundance and Lorenz analyses |
| 04 | Annotation count / percentage barplots |
| 05 | Read-size distribution barplots per class |
| 06 | piRNA degradation profiles over tRNA/snoRNA gene windows + overlap summary |
| s01 | Three-way comparison of the sense-overlap rules (comparison mode) |
| 10 | End-of-pipeline markdown summary + sanity checks |

## Outputs

All outputs are written under `output/` (single strategy) or
`output/comparison/` (comparison mode).

```
output/
  tables/     CSV tables (Table_01a–d, Table_02a/b, Table_03a/b, Table_06, ...)
  figures/    PDF + PNG figures (Figure_01, Figure_03a–c, Figure_04a–e,
                                  Figure_05a/b, Figure_06, Figure_s01a/b)
  parquet/    parquet intermediates (unique reads, annotations, position profiles)
  pipeline_summary.md   end-of-pipeline report
output/comparison/
  tables/     overlap-rule comparison tables (Table_s01a–e) + shared tables
  figures/    comparison figures + Figure_01
  fully_contained/  union/  any/
                per-strategy copies of tables/, figures/ and parquet/
```

Key tables:

| Table | Content |
|---|---|
| `Table_01a_sample_summary.csv` | per-sample total / unique reads, singletons, dominant sizes |
| `Table_01b/c/d` | read-size distribution, read-count distribution, size × count grid |
| `Table_02_<sample>_unique_reads_annotation.csv` | per-read annotation + abundance |
| `Table_02a/b` | annotation composition for unique reads / all reads |
| `Table_03a/b` | per-category and per-locus abundance statistics |
| `Table_06_piRNA_overlap_summary.csv` | piRNA reads overlapping snoRNA/tRNA genes (sense / antisense / none) |
| `Table_s01a–e` | overlap-rule comparison outputs (comparison mode) |
| `Table_10_sanity_checks.csv` | end-of-pipeline sanity checks |

## Development

```bash
.venv/bin/ruff check src tests
.venv/bin/black --check src tests
.venv/bin/pytest
```

## License

MIT.
