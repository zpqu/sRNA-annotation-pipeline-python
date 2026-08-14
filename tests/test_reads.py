"""Unit tests for step-01 read collapsing and summary tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from srna.reads.collapse import _unique_counts, add_poskey, collapse_bam
from srna.reads.summaries import (
    build_table_01a,
    build_table_01b,
    build_table_01c,
    build_table_01d,
    subsample_counts,
)


def _write_bam(path: Path) -> None:
    """Write a tiny coordinate-sorted BAM with 4 alignments over 2 contigs."""
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
    }
    rows = [
        ("chr1", 10, 18, False),
        ("chr1", 10, 18, False),
        ("chr1", 50, 70, True),
        ("chr2", 5, 9, False),
    ]
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (chrom, start, end, rev) in enumerate(rows):
            a = pysam.AlignedSegment()
            a.query_name = f"read{i}"
            a.reference_id = 0 if chrom == "chr1" else 1
            a.reference_start = start
            a.cigarstring = f"{end - start}M"
            a.query_sequence = "A" * (end - start)
            a.flag = 16 if rev else 0
            out.write(a)
    pysam.index(str(path))


class TestUniqueCounts:
    def test_groups_identical_reads(self):
        start = np.array([10, 10, 20, 10], dtype=np.int64)
        end = np.array([18, 18, 28, 18], dtype=np.int64)
        strand = np.array([0, 0, 1, 0], dtype=np.int8)
        u_start, u_end, u_strand, counts = _unique_counts(start, end, strand)
        assert list(u_start) == [10, 20]
        assert list(u_end) == [18, 28]
        assert list(u_strand) == [0, 1]
        assert list(counts) == [3, 1]

    def test_empty_input(self):
        u_start, u_end, u_strand, counts = _unique_counts(
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int8),
        )
        assert len(u_start) == 0 and len(counts) == 0
        assert u_strand.dtype.kind == "i"


class TestAddPoskey:
    def test_one_based_poskey(self):
        df = pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [18], "strand": ["+"]})
        out = add_poskey(df)
        assert out.loc[0, "poskey"] == "chr1:11:18:+"


class TestCollapseBam:
    def test_collapse_matches_r_semantics(self, tmp_path):
        bam = tmp_path / "reads.bam"
        _write_bam(bam)
        unique, total = collapse_bam(bam)
        assert total == 4
        assert len(unique) == 3
        counts = dict(zip(unique["start"], unique["count"], strict=True))
        assert counts == {5: 1, 10: 2, 50: 1}
        reverse = unique[unique["strand"] == "-"]
        assert reverse["start"].tolist() == [50]
        assert reverse["end"].tolist() == [70]


class TestSummaries:
    def _unique(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "start": [10, 20, 50],
                "end": [18, 28, 70],
                "strand": ["+", "+", "-"],
                "count": [3, 1, 1],
            }
        )

    def test_table_01a(self):
        row = build_table_01a(self._unique(), total_reads=5, sample="S")
        assert row["sample"] == "S"
        assert row["total_reads"] == 5
        assert row["unique_reads"] == 3
        assert row["singleton_reads"] == 2
        assert row["pct_singletons"] == 66.67
        assert row["median_count"] == 1.0
        assert row["max_count"] == 3
        assert row["dominant_size_by_unique_nt"] == 8
        assert row["dominant_size_by_unique_n"] == 2
        assert row["dominant_size_by_reads_nt"] == 8
        assert row["top5_sizes_nt"] == "8;20"

    def test_table_01b(self):
        tab = build_table_01b(self._unique(), "S")
        assert list(tab["sample"]) == ["S", "S"]
        assert list(tab["width"]) == [8, 20]
        assert list(tab["n_unique"]) == [2, 1]
        assert list(tab["n_reads"]) == [4, 1]

    def test_table_01c(self):
        tab = build_table_01c(self._unique(), "S")
        assert list(tab["lc"]) == [0, 2]
        assert list(tab["n_unique"]) == [2, 1]

    def test_table_01d(self):
        tab = build_table_01d(self._unique(), "S")
        pairs = {
            (w, lc): n for w, lc, n in tab[["width", "lc", "n_unique"]].itertuples(index=False)
        }
        assert pairs == {(8, 0): 1, (8, 2): 1, (20, 0): 1}

    def test_subsample_counts(self):
        values = subsample_counts(self._unique(), max_points=100)
        assert len(values) == 3
        assert set(np.round(values, 6)) == {0.0, float(np.round(np.log2(3), 6))}
