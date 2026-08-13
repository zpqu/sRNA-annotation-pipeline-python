# srna-annotation

Production-quality small-RNA annotation pipeline (Python), refactored from the
original R implementation in `scripts/R/`.

## Overview

The pipeline annotates small-RNA sequencing reads (BAM) against a set of genomic
features (mature miRNAs, piRNAs, snoRNAs, tRNAs, RepeatMasker, refGene
exons/introns/lincRNAs) with a configurable sense-overlap rule, then produces
abundance analyses, figures and an end-of-pipeline report.

Pipeline steps (mirroring `scripts/run_smallRNA_annotation.sh`):

| Step | Purpose |
|---|---|
| 00 | Build the genomic-feature DB from raw annotation files (`data/DB/`) and cache to parquet |
| 01 | Preprocess BAM files: read-size / count distribution tables and figure |
| 02 | Annotate non-redundant reads by priority (sense + antisense + other) |
| 03 | Per-locus abundance, Gini, rank-abundance and Lorenz analyses |
| 04 | Annotation count / percentage barplots |
| 05 | Read-size distribution barplots per class |
| 06 | piRNA degradation profiles over tRNA/snoRNA gene windows |
| s01 | Three-way comparison of the sense-overlap rules (comparison mode) |
| 10 | End-of-pipeline markdown summary + sanity checks |

## Installation

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

## Usage

```bash
.venv/bin/srna run [GENOME] [STRATEGY]
```

- `GENOME`: assembly id, default `mm39` (raw files are read from
  `data/DB/`, caches are written to `data/cache/<GENOME>/`).
- `STRATEGY`: sense-overlap rule, one of `fully-contained` (default),
  `union`, `any`, or `comparison`.

Outputs are written to `output/` (single strategy) or `output/comparison/`
(comparison mode).

## Development

```bash
.venv/bin/ruff check . && .venv/bin/black --check .
.venv/bin/pytest
```
