"""Unit tests for step-06 (piRNA-on-tRNA/snoRNA window analysis)."""

from __future__ import annotations

import pandas as pd

from srna.analysis.step06 import _overlap_table, _window_counts, sliding_windows


def _windows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 150],
            "end": [121, 171],
            "strand": ["+", "+"],
            "position": [1, 1],
            "gene_id": ["geneA", "geneB"],
        }
    )


def _reads() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 150],
            "end": [105, 151],
            "strand": ["+", "-"],
            "count": [5, 2],
        }
    )


class TestSlidingWindows:
    def test_plus_strand_windows(self):
        feat = pd.DataFrame(
            {"chrom": ["chr1"], "start": [100], "end": [130], "strand": ["+"], "gene_id": ["g"]}
        )
        out = sliding_windows(feat)
        assert len(out) == 10
        assert list(out["position"]) == list(range(1, 11))
        assert list(out["start"]) == list(range(100, 110))
        assert list(out["end"]) == list(range(121, 131))
        assert (out["strand"] == "+").all()
        assert (out["gene_id"] == "g").all()

    def test_minus_strand_windows_reversed(self):
        feat = pd.DataFrame(
            {"chrom": ["chr1"], "start": [100], "end": [130], "strand": ["-"], "gene_id": ["g"]}
        )
        out = sliding_windows(feat)
        # position 1 is the 5' end (right-hand side for minus strand)
        assert out.loc[out["position"] == 1, "start"].iloc[0] == 109
        assert out.loc[out["position"] == 10, "start"].iloc[0] == 100

    def test_short_gene_produces_no_windows(self):
        feat = pd.DataFrame(
            {"chrom": ["chr1"], "start": [100], "end": [120], "strand": ["+"], "gene_id": ["g"]}
        )
        out = sliding_windows(feat)
        assert len(out) == 0

    def test_empty_features(self):
        feat = pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_id"])
        out = sliding_windows(feat)
        assert len(out) == 0


class TestWindowCounts:
    def test_sense_counts_only_same_strand(self):
        counts = _window_counts(_windows(), _reads(), ignore_strand=False)
        assert counts.tolist() == [5, 0]

    def test_ignore_strand_counts_cross_strand(self):
        counts = _window_counts(_windows(), _reads(), ignore_strand=True)
        assert counts.tolist() == [5, 2]

    def test_empty_reads(self):
        empty = pd.DataFrame(columns=["chrom", "start", "end", "strand", "count"])
        assert _window_counts(_windows(), empty, ignore_strand=False).tolist() == [0, 0]


class TestOverlapTable:
    def _per_read(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chr": ["chr1", "chr1", "chr1"],
                "start": [5, 6, 50],
                "end": [15, 16, 60],
                "strand": ["+", "-", "+"],
                "category": ["piRNA", "piRNA", "piRNA"],
                "count": [3, 4, 1],
            }
        )

    def _store(self):
        class Store:
            def table(self, name: str) -> pd.DataFrame:
                if name == "snoRNA":
                    return pd.DataFrame(
                        {
                            "chrom": ["chr1"],
                            "start": [0],
                            "end": [20],
                            "strand": ["+"],
                            "gene_id": ["SnordX"],
                        }
                    )
                return pd.DataFrame(columns=["chrom", "start", "end", "strand", "gene_id"])

        return Store()

    def test_classification(self):
        tab = _overlap_table(self._per_read(), self._store())
        rows = {(r.category): r for r in tab.itertuples(index=False)}
        assert rows["snoRNA_sense"].n_unique == 1
        assert rows["snoRNA_sense"].n_reads == 3
        assert rows["snoRNA_antisense"].n_unique == 1
        assert rows["snoRNA_antisense"].n_reads == 4
        assert rows["none"].n_reads == 1
        assert rows["snoRNA_sense"].pct_reads == 37.5
        assert round(rows["none"].pct_reads, 1) == 12.5

    def test_empty_per_read(self):
        empty = pd.DataFrame(columns=["chr", "start", "end", "strand", "category", "count"])
        tab = _overlap_table(empty, self._store())
        assert len(tab) == 0
