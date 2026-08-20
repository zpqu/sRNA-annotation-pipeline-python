"""Pipeline configuration: genome, annotation strategy and chromosome style.

The behaviour of every pipeline step is controlled by a :class:`PipelineConfig`
object, which controls the genome, sense-overlap rule and chromosome naming
convention.

Sense-overlap rules
-------------------
``fully-contained``
    a read is annotated to a feature only if the read is fully contained in the
    feature.
``union``
    read within feature **or** feature within read (captures long reads
    spanning a small feature).
``any``
    any overlap (>= 1 bp) between read and feature.
``comparison``
    not a rule itself; runs the three rules above and then the overlap-rule
    comparison analysis (step s01).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import Literal


class Strategy(Enum):
    """Canonical sense-overlap strategies supported by the pipeline."""

    FULLY_CONTAINED = "fully-contained"
    UNION = "union"
    ANY = "any"
    COMPARISON = "comparison"

    @property
    def is_comparison(self) -> bool:
        """Return True for the comparison mode."""
        return self is Strategy.COMPARISON

    @property
    def overlap_rule(self) -> str:
        """Return the overlap rule used for the sense pass."""
        if self is Strategy.UNION:
            return "union"
        if self is Strategy.ANY:
            return "any"
        return "within"

    @property
    def dir_name(self) -> str:
        """Return the output sub-directory name used for this strategy."""
        return self.value.replace("-", "_")


_RULE_ALIASES: dict[str, str] = {
    "fully-contained": "fully-contained",
    "fully_contained": "fully-contained",
    "within": "fully-contained",
    "union": "union",
    "any": "any",
    "comparison": "comparison",
    "compar": "comparison",
}


def canon_strategy(value: str) -> Strategy:
    """Canonicalize a strategy string to a :class:`Strategy`.

    Normalizes: lowercased, whitespace/underscores collapsed to dashes,
    then matched by prefix.
    """
    normalized = re.sub(r"[ _]+", "-", value.strip().lower())
    for rule, canonical in _RULE_ALIASES.items():
        if normalized == rule or normalized.startswith(rule):
            return Strategy(canonical)
    if "union" in normalized:
        return Strategy.UNION
    if "any" in normalized:
        return Strategy.ANY
    if "compar" in normalized:
        return Strategy.COMPARISON
    return Strategy.FULLY_CONTAINED


_CHR_STYLES = ("chr", "", "none")


@dataclass(frozen=True)
class PipelineConfig:
    """Immutable runtime configuration for one pipeline invocation.

    Attributes:
        genome: assembly id used to locate raw feature files and caches
            (e.g. ``"mm39"``, ``"hg38"``).
        strategy: canonical sense-overlap strategy.
        chr_style: chromosome naming convention (``"chr"`` or ``""``/``"none"``).
        force_rebuild_db: force a full feature-DB rebuild even when cached.

    """

    genome: str = "mm39"
    strategy: Strategy = Strategy.FULLY_CONTAINED
    chr_style: Literal["chr", ""] = "chr"
    force_rebuild_db: bool = False

    def __post_init__(self) -> None:
        if self.chr_style not in _CHR_STYLES:
            raise ValueError(f"chr_style must be one of {_CHR_STYLES}, got {self.chr_style!r}")

    @classmethod
    def from_env(cls, env: dict[str, str] | None = None) -> PipelineConfig:
        """Build a config from environment variables.

        Reads ``SMALLRNA_GENOME``, ``SMALLRNA_STRATEGY``,
        ``SMALLRNA_CHR_STYLE`` and ``SMALLRNA_FORCE_REBUILD_DB``.
        """
        env = env if env is not None else os_environ()
        strategy = canon_strategy(env.get("SMALLRNA_STRATEGY", "fully-contained"))
        chr_style = env.get("SMALLRNA_CHR_STYLE", "chr")
        if chr_style == "none":
            chr_style = ""
        return cls(
            genome=env.get("SMALLRNA_GENOME", "mm39"),
            strategy=strategy,
            chr_style=chr_style,  # type: ignore[arg-type]
            force_rebuild_db=env.get("SMALLRNA_FORCE_REBUILD_DB", "0") == "1",
        )

    @classmethod
    def for_substrategy(cls, config: PipelineConfig, strategy: Strategy) -> PipelineConfig:
        """Return a copy of *config* running *strategy* (used inside comparison mode)."""
        return cls(
            genome=config.genome,
            strategy=strategy,
            chr_style=config.chr_style,
            force_rebuild_db=config.force_rebuild_db,
        )


def os_environ() -> dict[str, str]:
    """Return ``os.environ`` as a plain dict (imported lazily to keep module pure)."""
    import os

    return dict(os.environ)
