"""Chromosome-name normalization and primary-chromosome filtering.

Port of the helpers in ``scripts/R/00_build_DB/00_build_annotation_DB.R``:

* ``normalize_seqnames`` -- make every chromosome name agree with the reads
  (``chr1``/``1``/``Chr1`` -> ``chr1``, ``MT`` -> ``chrM``).
* ``is_primary`` / ``primary`` -- keep only primary chromosomes / scaffolds,
  dropping alts, randoms, fixations, haplotypes and unplaced contigs.
"""

from __future__ import annotations

import re

from numpy.typing import NDArray

_PRIMARY_DROP = re.compile(r"_alt|_random|_fix|_hap|chrUn|^Un|^HSCHR")


def normalize_seqnames(names: list[str], chr_style: str = "chr") -> list[str]:
    """Normalize a list of chromosome names to *chr_style*.

    For ``chr_style == "chr"``: strip a leading ``Chr``, prepend ``chr`` when
    missing and map ``chrMT`` -> ``chrM``. For ``chr_style == ""``: strip a
    leading ``chr``. Any other style returns the names unchanged.
    """
    out: list[str] = []
    for raw in names:
        name = re.sub(r"^Chr", "chr", raw)
        if chr_style == "chr":
            if not name.startswith("chr"):
                name = "chr" + name
            if name == "chrMT":
                name = "chrM"
        elif chr_style == "":
            name = re.sub(r"^chr", "", name)
        out.append(name)
    return out


def is_primary(name: str) -> bool:
    """Return True when *name* is a primary chromosome / scaffold."""
    return not bool(_PRIMARY_DROP.search(name))


def filter_primary(names: NDArray | list[str]) -> list[str]:
    """Filter *names* down to primary chromosomes (order preserved)."""
    return [n for n in names if is_primary(n)]
