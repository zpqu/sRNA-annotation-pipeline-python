"""Parse the miRNA GFF3 (miRBase) into mature and primary-transcript tables."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .base import FEATURE_COLUMNS


def parse_mirna_gff3(path: Any) -> dict[str, pd.DataFrame]:
    """Parse a miRBase GFF3 file.

    Returns:
        ``{"matmiRNA": ..., "primiRNA": ...}`` tables keyed by mature miRNA
        (``Name``, e.g. ``mmu-miR-206-5p``) and primary transcript.

    """
    from .parsing import extract_gff_fields, read_gff3

    gff = read_gff3(path)
    fields = extract_gff_fields(gff, ["ID", "Name"])
    gff = pd.concat([gff, fields], axis=1)

    name = gff["Name"].astype(str)
    gid = gff["ID"].astype(str)

    mature = gff[gff["feature"] == "miRNA"]
    primary = gff[gff["feature"] == "miRNA_primary_transcript"]

    tables: dict[str, pd.DataFrame] = {}
    for key, df, label in (
        ("matmiRNA", mature, "matmiRNA"),
        ("primiRNA", primary, "primiRNA"),
    ):
        if df.empty:
            tables[key] = _empty(label)
            continue
        tables[key] = pd.DataFrame(
            {
                "chrom": df["chrom"].astype(str),
                "start": df["start"].astype(np.int64),
                "end": df["end"].astype(np.int64),
                "strand": df["strand"].astype(str),
                "feature_type": label,
                "feature_id": df["Name"].fillna(df["ID"]).astype(str),
                "gene_id": gid.loc[df.index].astype(str),
                "gene_name": name.loc[df.index].astype(str),
                "source": "miRBase",
                "score": np.nan,
            }
        )[FEATURE_COLUMNS].reset_index(drop=True)
    return tables


def _empty(feature_type: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": pd.Series(dtype="string"),
            "start": pd.Series(dtype=np.int64),
            "end": pd.Series(dtype=np.int64),
            "strand": pd.Series(dtype="string"),
            "feature_type": pd.Series(dtype="string"),
            "feature_id": pd.Series(dtype="string"),
            "gene_id": pd.Series(dtype="string"),
            "gene_name": pd.Series(dtype="string"),
            "source": pd.Series(dtype="string"),
            "score": pd.Series(dtype="float64"),
        }
    )
