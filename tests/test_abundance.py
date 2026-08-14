"""Unit tests for step-03 abundance analysis functions."""

from __future__ import annotations

import numpy as np
import pandas as pd

from srna.analysis.locus_abundance import (
    _gini,
    _locus_reads,
    _spearman,
    _wilcoxon_paired,
    abundance_summary_table,
    per_locus_abundance_table,
)


class FakeStore:
    def __init__(self, tables: dict[str, pd.DataFrame]):
        self._tables = tables

    def table(self, name: str) -> pd.DataFrame:
        return self._tables[name]


class TestGini:
    def test_perfect_equality(self):
        assert _gini(np.array([1.0, 1.0])) == 0.0

    def test_known_value(self):
        assert round(_gini(np.array([1.0, 2.0, 3.0])), 4) == 0.2222

    def test_degenerate_inputs(self):
        assert np.isnan(_gini(np.array([])))
        assert np.isnan(_gini(np.array([5.0])))
        assert np.isnan(_gini(np.array([0.0, 0.0])))


class TestSpearman:
    def test_perfect_positive(self):
        assert _spearman(pd.Series([1.0, 2.0, 3.0]), pd.Series([1.0, 2.0, 3.0])) == 1.0

    def test_perfect_negative(self):
        assert _spearman(pd.Series([1.0, 2.0, 3.0]), pd.Series([3.0, 2.0, 1.0])) == -1.0

    def test_too_few_points(self):
        assert np.isnan(_spearman(pd.Series([1.0]), pd.Series([2.0])))


class TestWilcoxonPaired:
    def test_identical_data_is_nan(self):
        assert np.isnan(_wilcoxon_paired(pd.Series([3.0, 3.0, 3.0]), pd.Series([3.0, 3.0, 3.0])))

    def test_mixed_signs_weaker_than_consistent(self):
        all_same = _wilcoxon_paired(
            pd.Series([1.0, 2.0, 3.0, 4.0, 5.0]), pd.Series([0.0, 0.0, 0.0, 0.0, 0.0])
        )
        mixed = _wilcoxon_paired(
            pd.Series([1.0, 2.0, 3.0, 0.0, 0.0]), pd.Series([0.0, 0.0, 0.0, 4.0, 5.0])
        )
        assert all_same < mixed

    def test_matches_normal_approximation(self):
        p = _wilcoxon_paired(
            pd.Series([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]), pd.Series([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        )
        # manual normal-approximation value for n=6 with all-positive diffs
        assert round(p, 3) == 0.036


class TestLocusReads:
    def _per_read(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [12, 100],
                "end": [18, 110],
                "strand": ["+", "+"],
                "category": ["matmiRNA", "piRNA"],
                "count": [5, 3],
            }
        )

    def _store(self) -> FakeStore:
        mat = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [10, 15],
                "end": [20, 25],
                "strand": ["+", "+"],
                "gene_name": ["mir-A", "mir-B"],
            }
        )
        return FakeStore({"matmiRNA": mat})

    def test_any_overlap_locus_remap_and_pirna_positions(self):
        tab = _locus_reads(self._per_read(), "S1", self._store())
        rows = {(r.sample, r.category, r.locus): r.n_reads for r in tab.itertuples(index=False)}
        assert rows == {
            ("S1", "matmiRNA", "mir-A"): 5,
            ("S1", "matmiRNA", "mir-B"): 5,
            ("S1", "piRNA", "chr1 100 110 +"): 3,
        }


class TestAbundanceSummaryTable:
    def _locus_tab(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "sample": ["A", "A", "B", "B"],
                "category": ["matmiRNA", "matmiRNA", "matmiRNA", "matmiRNA"],
                "locus": ["L1", "L2", "L1", "L2"],
                "n_reads": [100, 50, 80, 70],
            }
        )

    def test_summary_values(self):
        tab = abundance_summary_table(self._locus_tab(), ["A", "B"])
        assert len(tab) == 2
        row = tab[tab["sample"] == "A"].iloc[0]
        assert row["n_loci"] == 2
        assert row["total_reads"] == 150
        assert row["top_locus"] == "L1"
        assert row["top1_pct"] == 66.67
        assert row["gini"] == 0.167
        assert row["spearman_cross_sample"] == 1.0
        assert not pd.isna(row["wilcoxon_p_2sample"])

    def test_single_sample_no_cross_test(self):
        tab = abundance_summary_table(self._locus_tab().query("sample == 'A'"), ["A"])
        assert np.isnan(tab.iloc[0]["spearman_cross_sample"])


class TestPerLocusAbundanceTable:
    def test_sorted_by_category_sample_reads(self):
        tab = per_locus_abundance_table(
            pd.DataFrame(
                {
                    "sample": ["B", "A", "A"],
                    "category": ["matmiRNA", "matmiRNA", "piRNA"],
                    "locus": ["x", "y", "z"],
                    "n_reads": [1, 10, 5],
                }
            )
        )
        assert list(tab["sample"]) == ["A", "B", "A"]
        assert list(tab["category"]) == ["matmiRNA", "matmiRNA", "piRNA"]
        assert list(tab["n_reads"]) == [10, 1, 5]
