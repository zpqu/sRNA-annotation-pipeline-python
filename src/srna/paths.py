"""Filesystem layout: every path used by the pipeline, derived from one root.

Centralises the directory layout:

* ``data/bam/``   -- input BAM files
* ``data/DB/``    -- raw annotation files for the reference genome
* ``data/cache/<genome>/`` -- rebuilt feature tables (parquet + meta.json)
* ``output/``     -- single-strategy run results
* ``output/comparison/`` -- comparison-mode run results
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, Strategy

_SAMPLES_RE = (".bam",)
_ALIGNER_SUFFIXES = (".bwa", ".bowtie2")

_RAW_DB_FILES = (
    "refGene.gtf",
    "miRNA.gff3",
    "tRNAs.bed",
    "RM.bed",
    "piRBase.bed",
    "piRNAdb.gtf",
)


def derive_sample_label(bam_name: str) -> str:
    """Derive a sample label from a BAM file name.

    Strips ``.bam`` and then a trailing ``.bwa``/``.bowtie2`` aligner suffix,
    e.g. ``Cumulus-cells.bwa.bam -> Cumulus-cells``.
    """
    name = bam_name
    for suffix in _SAMPLES_RE:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in _ALIGNER_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


@dataclass(frozen=True)
class PathResolver:
    """Resolve all input/output paths for a pipeline configuration."""

    root: Path
    config: PipelineConfig

    # -- inputs -------------------------------------------------------------
    @property
    def bam_dir(self) -> Path:
        """Directory containing input BAM files."""
        return self.root / "data" / "bam"

    @property
    def db_raw_dir(self) -> Path:
        """Directory containing the raw annotation files for the genome."""
        return self.root / "data" / "DB"

    @property
    def db_cache_dir(self) -> Path:
        """Directory where rebuilt feature tables are cached (per genome)."""
        return self.root / "data" / "cache" / self.config.genome

    # -- outputs ------------------------------------------------------------
    @property
    def output_base(self) -> Path:
        """Root output directory (``output/`` or ``output/comparison/``)."""
        if self.config.strategy.is_comparison:
            return self.root / "output" / "comparison"
        return self.root / "output"

    @property
    def out_dir(self) -> Path:
        """Directory for the current strategy (comparison sub-folder or base)."""
        if self.config.strategy.is_comparison:
            return self.output_base
        return self.output_base

    def strategy_dir(self, strategy: Strategy) -> Path:
        """Return the output directory for a concrete (non-comparison) strategy."""
        return self.output_base / strategy.dir_name

    # -- logs ---------------------------------------------------------------
    @property
    def log_dir(self) -> Path:
        """Directory where pipeline logs are written."""
        return self.root / "logs"

    # -- helpers ------------------------------------------------------------
    def samples(self) -> list[str]:
        """Auto-detect sample labels from ``bams/*.bam`` (sorted, unique)."""
        labels = sorted({derive_sample_label(f.name) for f in self.bam_dir.glob("*.bam")})
        return labels

    def bam_file(self, sample: str) -> Path:
        """Locate the BAM file for a sample (tolerates ``.bwa``/``.bowtie2`` names)."""
        candidates = sorted(self.bam_dir.glob(f"{sample}*.bam"))
        if not candidates:
            raise FileNotFoundError(f"no BAM file for sample {sample!r} in {self.bam_dir}")
        return candidates[0]

    def require_raw_db_files(self) -> None:
        """Raise if any raw annotation file needed to build the DB is missing."""
        missing = [f for f in _RAW_DB_FILES if not (self.db_raw_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"raw feature files missing in {self.db_raw_dir}: {', '.join(missing)}"
            )
