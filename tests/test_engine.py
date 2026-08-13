"""End-to-end annotation engine tests on synthetic data (no real DB required)."""

from __future__ import annotations

import pandas as pd

from srna.annotation.engine import annotate_sample
from srna.config import Strategy
from srna.features.base import FEATURE_COLUMNS, make_frame

SAMPLE = "test-sample"


class FakeStore:
    """Minimal FeatureStore stand-in exposing ``chrom_table`` / ``table``."""

    def __init__(self, tables: dict[str, pd.DataFrame]):
        self._tables = tables

    def table(self, name: str) -> pd.DataFrame:
        return self._tables[name]

    def chrom_table(self, name: str, chrom: str) -> pd.DataFrame:
        df = self._tables[name]
        return df[df["chrom"] == chrom]


def _feature_table(rows: list[tuple], feature_type: str) -> pd.DataFrame:
    """Build a feature table from rows of feature fields."""
    frame = make_frame(
        chrom=[r[0] for r in rows],
        start=[r[1] for r in rows],
        end=[r[2] for r in rows],
        strand=[r[3] for r in rows],
        feature_type=feature_type,
        feature_id=[r[4] for r in rows],
        gene_id=[r[5] for r in rows],
        gene_name=[r[6] for r in rows],
        source="test",
        score=[r[7] if len(r) > 7 else None for r in rows],
    )
    return frame[FEATURE_COLUMNS]


def _reads() -> pd.DataFrame:
    """Synthetic unique reads: ``(chrom, start, end, strand, count)``."""
    rows = [
        ("chr1", 12, 18, "+", 5),
        ("chr1", 16, 20, "+", 3),
        ("chr1", 31, 39, "+", 2),
        ("chr1", 61, 69, "-", 1),
        ("chr1", 85, 90, "+", 7),
        ("chr1", 210, 290, "+", 4),
        ("chr1", 26, 29, "+", 6),
        ("chr1", 17, 24, "+", 8),
    ]
    return pd.DataFrame(
        {
            "chrom": [r[0] for r in rows],
            "start": [r[1] for r in rows],
            "end": [r[2] for r in rows],
            "strand": [r[3] for r in rows],
            "count": [r[4] for r in rows],
        }
    )


def _store() -> FakeStore:
    tables: dict[str, pd.DataFrame] = {
        "matmiRNA": _feature_table(
            [("chr1", 10, 20, "+", "mir-1", "MIMAT0001", "mir-1", None)], "matmiRNA"
        ),
        "piRNA": _feature_table([("chr1", 30, 40, "+", "piR-1", "", "", 5)], "piRNA"),
        "snoRNA": _feature_table(
            [("chr1", 60, 70, "-", "NR_x", "SnordX", "SnordX", None)], "snoRNA"
        ),
        "tRNA": _feature_table([("chr1", 80, 90, "+", "tRNA-A", "tRNA-A", "tRNA-A", None)], "tRNA"),
        "RM": _feature_table([("chr1", 200, 300, "+", "L1", "L1", "L1", None)], "RM"),
        "NM.exon": _feature_table([("chr1", 15, 25, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.exon"),
        "NM.intron": _feature_table(
            [("chr1", 40, 50, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.intron"
        ),
        "NM.CDS": _feature_table([("chr1", 15, 25, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.CDS"),
        "NM.5UTR": _feature_table([("chr1", 5, 15, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.5UTR"),
        "NM.3UTR": _feature_table([("chr1", 25, 30, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.3UTR"),
        "NM.up1k": _feature_table([("chr1", 0, 5, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.up1k"),
        "NM.down1k": _feature_table(
            [("chr1", 30, 35, "+", "NM_1", "Xkr4", "Xkr4", None)], "NM.down1k"
        ),
        "lincRNA.exon": _feature_table([], "lincRNA.exon"),
    }
    return FakeStore(tables)


def test_fully_contained_annotation():
    result = annotate_sample(_reads(), _store(), Strategy.FULLY_CONTAINED, SAMPLE)
    tab = result.per_read.set_index("read_id")
    assert len(tab) == 8
    assert tab.loc[1, "category"] == "matmiRNA"
    assert tab.loc[1, "feature_id"] == "mir-1"
    assert tab.loc[1, "gene_context"] == "AS.CDS"
    assert tab.loc[1, "n_features"] == 1
    assert tab.loc[2, "category"] == "matmiRNA"
    assert tab.loc[2, "gene_context"] == "CDS"
    assert tab.loc[3, "category"] == "piRNA"
    assert tab.loc[3, "gene_context"] == "AS.down1k"
    assert tab.loc[4, "category"] == "snoRNA"
    assert tab.loc[4, "gene_context"] == "intergenic"
    assert tab.loc[5, "category"] == "tRNA"
    assert tab.loc[5, "gene_context"] == "intergenic"
    assert tab.loc[6, "category"] == "RM"
    assert tab.loc[6, "gene_context"] == "NA"
    assert tab.loc[7, "category"] == "refGene.NM.exon"
    assert tab.loc[7, "gene_context"] == "NA"
    assert tab.loc[7, "feature_id"] == "Xkr4"
    assert tab.loc[8, "category"] == "other"
    assert tab.loc[8, "gene_context"] == "NA"


def test_consolidated_counts_unique_reads():
    result = annotate_sample(_reads(), _store(), Strategy.FULLY_CONTAINED, SAMPLE)
    counts = {(a, b): c for a, b, c in result.counts_unique}
    assert counts[("matmiRNA.annotation", "AS.CDS")] == 1
    assert counts[("matmiRNA.annotation", "CDS")] == 1
    assert counts[("piRNA.annotation", "AS.down1k")] == 1
    assert counts[("snoRNA.annotation", "intergenic")] == 1
    assert counts[("tRNA.annotation", "intergenic")] == 1
    assert counts[("read.annotation", "matmiRNA")] == 2
    assert counts[("read.annotation", "other")] == 1


def test_consolidated_counts_all_reads():
    result = annotate_sample(_reads(), _store(), Strategy.FULLY_CONTAINED, SAMPLE)
    counts = {(a, b): c for a, b, c in result.counts_reads}
    assert counts[("matmiRNA.annotation", "AS.CDS")] == 5
    assert counts[("matmiRNA.annotation", "CDS")] == 3
    assert counts[("read.annotation", "matmiRNA")] == 8
    assert counts[("read.annotation", "refGene.NM.exon")] == 8
    assert counts[("read.annotation", "other")] == 6


def test_cpm_column():
    result = annotate_sample(_reads(), _store(), Strategy.FULLY_CONTAINED, SAMPLE)
    total = _reads()["count"].sum()
    for _, row in result.per_read.iterrows():
        assert row["cpm"] == round(row["count"] / total * 1e6, 3)


def test_strategy_union_captures_long_reads():
    long_read = pd.DataFrame(
        [("chr1", 5, 45, "+", 1)],
        columns=["chrom", "start", "end", "strand", "count"],
    )
    result = annotate_sample(long_read, _store(), Strategy.UNION, SAMPLE)
    tab = result.per_read
    assert len(tab) == 1
    assert tab.loc[0, "category"] == "matmiRNA"


def test_strategy_any_overlap():
    read = pd.DataFrame(
        [("chr1", 8, 15, "+", 1)],
        columns=["chrom", "start", "end", "strand", "count"],
    )
    result = annotate_sample(read, _store(), Strategy.ANY, SAMPLE)
    tab = result.per_read
    assert tab.loc[0, "category"] == "matmiRNA"


def test_pirna_max_score_selection():
    store = _store()
    pirna = _feature_table(
        [
            ("chr1", 30, 40, "+", "piR-low", "", "", 2),
            ("chr1", 30, 40, "+", "piR-high", "", "", 9),
        ],
        "piRNA",
    )
    store._tables["piRNA"] = pirna
    read = pd.DataFrame(
        [("chr1", 33, 37, "+", 1)],
        columns=["chrom", "start", "end", "strand", "count"],
    )
    result = annotate_sample(read, store, Strategy.FULLY_CONTAINED, SAMPLE)
    assert result.per_read.loc[0, "feature_id"] == "piR-high"
    assert result.per_read.loc[0, "n_features"] == 2


def test_empty_table_handling():
    empty = pd.DataFrame(columns=["chrom", "start", "end", "strand", "count"])
    result = annotate_sample(empty, _store(), Strategy.FULLY_CONTAINED, SAMPLE)
    assert len(result.per_read) == 0
