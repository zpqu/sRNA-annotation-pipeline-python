"""Feature DB store: builds all feature tables and caches them to parquet.

Each feature table is rebuilt from its raw source file(s) in ``data/DB/`` and
cached under ``data/cache/<genome>/``. Caches are invalidated when any source
file's mtime or size changes (recorded in a ``.meta.json`` sidecar).
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger

from ..chrom import filter_primary, normalize_seqnames
from ..config import PipelineConfig
from .mirna import parse_mirna_gff3
from .pirna import parse_pirna
from .refgene import parse_refgene_gtf
from .schema import CANONICAL_TYPES
from .trna_rm import parse_rm_bed, parse_trna_bed

_CACHE_VERSION = 1

#: canonical table name -> raw source file(s) it depends on
_SOURCES: dict[str, tuple[str, ...]] = {
    "matmiRNA": ("miRNA.gff3",),
    "primiRNA": ("miRNA.gff3",),
    "piRNA": ("piRBase.bed", "piRNAdb.gtf"),
    "tRNA": ("tRNAs.bed",),
    "RM": ("RM.bed",),
    "NM.exon": ("refGene.gtf",),
    "NM.CDS": ("refGene.gtf",),
    "NM.5UTR": ("refGene.gtf",),
    "NM.3UTR": ("refGene.gtf",),
    "NM.intron": ("refGene.gtf",),
    "NM.mRNA": ("refGene.gtf",),
    "NM.up1k": ("refGene.gtf",),
    "NM.down1k": ("refGene.gtf",),
    "NR.exon": ("refGene.gtf",),
    "snoRNA": ("refGene.gtf",),
    "lincRNA.exon": ("refGene.gtf",),
}


class FeatureStore:
    """Builds and loads the cached feature tables for one genome."""

    def __init__(self, raw_dir: Path, cache_dir: Path, config: PipelineConfig):
        self.raw_dir = raw_dir
        self.cache_dir = cache_dir
        self.config = config

    # -- public API -----------------------------------------------------------
    def build(self) -> None:
        """Rebuild missing or stale caches for every feature table."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        for name in CANONICAL_TYPES:
            self._ensure_table(name)

    def table(self, name: str) -> pd.DataFrame:
        """Load a whole feature table as a DataFrame."""
        self._ensure_table(name)
        return self._read(name)

    def chrom_table(self, name: str, chrom: str) -> pd.DataFrame:
        """Load the features of *name* restricted to one chromosome."""
        self._ensure_table(name)
        path = self._parquet_path(name)
        if _is_partitioned(name):
            df = pq.read_table(path, filters=[("chrom", "=", chrom)]).to_pandas()
        else:
            df = pd.read_parquet(path)
            df = df[df["chrom"] == chrom]
        return df

    def all_tables(self) -> dict[str, pd.DataFrame]:
        """Load every feature table (used by the annotation engine)."""
        return {name: self.table(name) for name in CANONICAL_TYPES}

    # -- cache management -----------------------------------------------------
    def _ensure_table(self, name: str) -> None:
        if self.config.force_rebuild_db or not self._cache_valid(name):
            logger.info("building feature table '{name}' ...", name=name)
            self._build_one(name)

    def _cache_valid(self, name: str) -> bool:
        path = self._parquet_path(name)
        meta = self._meta_path(name)
        exists = path.is_dir() if _is_partitioned(name) else path.is_file()
        if not exists or not meta.exists():
            return False
        try:
            current = self._source_signatures(_SOURCES[name])
            return json.loads(meta.read_text()) == {
                "version": _CACHE_VERSION,
                "sources": current,
            }
        except (json.JSONDecodeError, OSError):
            return False

    def _build_one(self, name: str) -> None:
        builder = _builders()[name]
        data = builder(self.raw_dir, self.config)
        self._write(name, data)
        self._write_meta(name)

    def _write(self, name: str, data: pa.Table | pd.DataFrame) -> None:
        path = self._parquet_path(name)
        if isinstance(data, pa.Table):
            table = data.sort_by([("chrom", "ascending")])
        else:
            table = pa.Table.from_pandas(data, preserve_index=False)
        n = table.num_rows
        if _is_partitioned(name):
            pq.write_to_dataset(table, root_path=path, partition_cols=["chrom"])
            logger.info(
                "cached '{name}' -> {path}/ (partitioned by chromosome, {n} rows)",
                name=name,
                path=path,
                n=n,
            )
        else:
            pq.write_table(table, path)
            logger.info("cached '{name}' -> {path} ({n} rows)", name=name, path=path, n=n)

    def _write_meta(self, name: str) -> None:
        meta = {
            "version": _CACHE_VERSION,
            "sources": self._source_signatures(_SOURCES[name]),
        }
        self._meta_path(name).write_text(json.dumps(meta, indent=2))

    def _source_signatures(self, sources: tuple[str, ...]) -> dict[str, dict[str, int]]:
        out: dict[str, dict[str, int]] = {}
        for source in sources:
            stat = (self.raw_dir / source).stat()
            out[source] = {"mtime": int(stat.st_mtime_ns), "size": stat.st_size}
        return out

    def _parquet_path(self, name: str) -> Path:
        return self.cache_dir / f"{name}.parquet"

    def _meta_path(self, name: str) -> Path:
        return self.cache_dir / f"{name}.meta.json"

    def _read(self, name: str) -> pd.DataFrame:
        path = self._parquet_path(name)
        if _is_partitioned(name):
            return pq.read_table(path).to_pandas()
        return pd.read_parquet(path)


def _is_partitioned(name: str) -> bool:
    """Return True for tables stored as a per-chromosome parquet dataset."""
    return name == "piRNA"


def _builders() -> dict[str, Callable[[Path, PipelineConfig], pa.Table | pd.DataFrame]]:
    """Return the canonical-name -> build-callable mapping."""

    def refgene(raw: Path, cfg: PipelineConfig) -> dict[str, pd.DataFrame]:
        return parse_refgene_gtf(raw / "refGene.gtf", cfg.genome)

    def mirna(raw: Path, cfg: PipelineConfig) -> dict[str, pd.DataFrame]:
        return parse_mirna_gff3(raw / "miRNA.gff3")

    def trna(raw: Path, cfg: PipelineConfig) -> pd.DataFrame:
        return _normalize_frame(parse_trna_bed(raw / "tRNAs.bed", cfg.genome), cfg)

    def rm(raw: Path, cfg: PipelineConfig) -> pd.DataFrame:
        return _normalize_frame(parse_rm_bed(raw / "RM.bed", cfg.genome), cfg)

    def pirna(raw: Path, cfg: PipelineConfig) -> pa.Table:
        return parse_pirna(
            raw / "piRBase.bed", raw / "piRNAdb.gtf", cfg.genome, chr_style=cfg.chr_style
        )

    builders: dict[str, Callable[[Path, PipelineConfig], pa.Table | pd.DataFrame]] = {}
    for name in _SOURCES:
        if name in ("matmiRNA", "primiRNA"):
            builders[name] = _select(name, mirna)
        elif name == "piRNA":
            builders[name] = pirna
        elif name == "tRNA":
            builders[name] = trna
        elif name == "RM":
            builders[name] = rm
        else:
            builders[name] = _select(name, refgene)
    return builders


def _select(name: str, builder: Callable[[Path, PipelineConfig], dict[str, pd.DataFrame]]):
    def selected(raw: Path, cfg: PipelineConfig) -> pd.DataFrame:
        tables = builder(raw, cfg)
        return _normalize_frame(tables[name], cfg)

    return selected


def _normalize_frame(df: pd.DataFrame, cfg: PipelineConfig) -> pd.DataFrame:
    """Apply chromosome normalization + primary filter to a feature table."""
    if df.empty:
        return df
    df = df.copy()
    df["chrom"] = pd.Series(
        normalize_seqnames(df["chrom"].astype(str).tolist(), chr_style=cfg.chr_style),
        index=df.index,
    )
    df = df[df["chrom"].isin(filter_primary(df["chrom"].unique().tolist()))]
    return df.reset_index(drop=True)
