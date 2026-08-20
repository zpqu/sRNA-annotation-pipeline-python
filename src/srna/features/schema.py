"""Unified genomic-feature table and the feature catalog.

Every feature set in the DB is represented as a :class:`pandas.DataFrame` with
a single, fixed column schema (0-based half-open intervals, matching pysam and
pyranges):

* ``chrom``  -- normalized chromosome name
* ``start``  -- 0-based start coordinate
* ``end``    -- exclusive end coordinate
* ``strand`` -- ``"+"`` or ``"-"``
* ``feature_type`` -- canonical feature type name
* ``feature_id``   -- element/transcript/locus id
* ``gene_id``, ``gene_name`` -- gene annotations
* ``source`` -- provenance (``refGene``, ``miRBase``, ...)
* ``score``  -- optional numeric score (piRNA loci)

The catalog in :data:`FEATURE_META` declares, per canonical type, the output
label, whether the feature is a flat DataFrame or a grouped DataFrame
whose elements must be deduplicated by ``feature_id`` during overlap
counting, and which column provides the feature id emitted into
the annotated-read output.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd

FEATURE_COLUMNS = [
    "chrom",
    "start",
    "end",
    "strand",
    "feature_type",
    "feature_id",
    "gene_id",
    "gene_name",
    "source",
    "score",
]

CANONICAL_TYPES = (
    "matmiRNA",
    "primiRNA",
    "piRNA",
    "snoRNA",
    "tRNA",
    "RM",
    "NM.exon",
    "NM.CDS",
    "NM.5UTR",
    "NM.3UTR",
    "NM.intron",
    "NM.mRNA",
    "NM.up1k",
    "NM.down1k",
    "NR.exon",
    "lincRNA.exon",
)


@dataclass(frozen=True)
class FeatureMeta:
    """Metadata describing one canonical feature type."""

    name: str
    out_label: str
    kind: str  # "gr" (flat) or "grl" (grouped by feature_id)
    id_col: str  # column emitted as the read's feature id


def _meta(name: str, out_label: str, kind: str, id_col: str) -> FeatureMeta:
    return FeatureMeta(name=name, out_label=out_label, kind=kind, id_col=id_col)


FEATURE_META: dict[str, FeatureMeta] = {
    t.name: t
    for t in (
        _meta("matmiRNA", "matmiRNA", "gr", "gene_name"),
        _meta("primiRNA", "primiRNA", "gr", "feature_id"),
        _meta("piRNA", "piRNA", "gr", "feature_id"),
        _meta("snoRNA", "snoRNA", "gr", "gene_id"),
        _meta("tRNA", "tRNA", "gr", "gene_id"),
        _meta("RM", "RM", "gr", "gene_id"),
        _meta("NM.exon", "refGene.NM.exon", "grl", "gene_name"),
        _meta("NM.CDS", "refGene.NM.CDS", "grl", "gene_name"),
        _meta("NM.5UTR", "refGene.NM.5UTR", "grl", "gene_name"),
        _meta("NM.3UTR", "refGene.NM.3UTR", "grl", "gene_name"),
        _meta("NM.intron", "refGene.NM.intron", "grl", "gene_name"),
        _meta("NM.mRNA", "refGene.NM.mRNA", "grl", "gene_name"),
        _meta("NM.up1k", "refGene.NM.up1k", "grl", "gene_name"),
        _meta("NM.down1k", "refGene.NM.down1k", "grl", "gene_name"),
        _meta("NR.exon", "refGene.NR.exon", "grl", "gene_name"),
        _meta("lincRNA.exon", "lincRNA.exon", "grl", "gene_name"),
    )
}

#: Order in which reads are annotated (sense pass) in step 02.
GENOMIC_FEATURES = [
    "matmiRNA",
    "piRNA",
    "snoRNA",
    "tRNA",
    "RM",
    "NM.exon",
    "NM.intron",
    "lincRNA.exon",
]

#: Gene-context features probed inside the matmiRNA/snoRNA/piRNA/tRNA
#: annotation, in region-priority order.
GENE_FEATURES = [
    "NM.CDS",
    "NM.5UTR",
    "NM.3UTR",
    "NM.intron",
    "NM.up1k",
    "NM.down1k",
    "RM",
]

#: region labels emitted for each gene-context feature.
GENE_REGION = {
    "NM.CDS": "CDS",
    "NM.5UTR": "5UTR",
    "NM.3UTR": "3UTR",
    "NM.intron": "intron",
    "NM.up1k": "up1k",
    "NM.down1k": "down1k",
    "RM": "RM",
}

#: Feature types whose hits get a gene-context region in step 02.
GENE_CONTEXT_CATEGORIES = ("matmiRNA", "snoRNA", "piRNA", "tRNA")


def empty_frame() -> pd.DataFrame:
    """Return an empty feature DataFrame with the canonical schema."""
    return pd.DataFrame({c: pd.Series(dtype=object) for c in FEATURE_COLUMNS})


def make_frame(
    chrom: Any,
    start: Any,
    end: Any,
    strand: Any,
    feature_type: str,
    feature_id: Any,
    gene_id: Any = None,
    gene_name: Any = None,
    source: str = "",
    score: Any = None,
) -> pd.DataFrame:
    """Build a feature DataFrame from equal-length array-likes."""
    return pd.DataFrame(
        {
            "chrom": list(chrom),
            "start": list(start),
            "end": list(end),
            "strand": list(strand),
            "feature_type": feature_type,
            "feature_id": list(feature_id),
            "gene_id": "" if gene_id is None else list(gene_id),
            "gene_name": "" if gene_name is None else list(gene_name),
            "source": source,
            "score": pd.Series(score, dtype="float64"),
        }
    )


def element_id_of(frame: pd.DataFrame, meta: FeatureMeta) -> pd.Series:
    """Return the per-element display id used as the read's ``feature_id``.

    For ``gr`` features this is the catalog ``id_col``. For ``grl`` features it
    is the gene name, falling back to the gene id.
    """
    if meta.kind == "gr":
        return frame[meta.id_col].astype(str)
    name = frame["gene_name"].fillna("").astype(str)
    gid = frame["gene_id"].fillna("").astype(str)
    return name.where(name != "", gid).astype(str)
