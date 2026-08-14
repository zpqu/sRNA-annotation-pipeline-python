"""Unit tests for :mod:`srna.annotation.overlap_pairs` (overlap-pair finding)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from srna.annotation.overlap_pairs import overlap_pairs


def _reads():
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [10, 30, 60, 10],
            "end": [20, 40, 70, 20],
            "strand": ["+", "-", "+", "+"],
            "read_idx": [0, 1, 2, 3],
        }
    )


def _features():
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2"],
            "start": [5, 25, 65, 12],
            "end": [25, 50, 68, 18],
            "strand": ["+", "-", "+", "+"],
            "feat_idx": [0, 1, 2, 3],
        }
    )


def _pairs_by_read(pairs: np.ndarray) -> dict[int, list[int]]:
    out: dict[int, list[int]] = {}
    for r, f in pairs:
        out.setdefault(int(r), []).append(int(f))
    return out


class TestOverlapPairs:
    def test_within_mode(self):
        pairs = _pairs_by_read(overlap_pairs(_reads(), _features(), "within"))
        # read 0 (10-20) within feat 0 (5-25); read 1 (30-40) within feat 1 (25-50)
        assert pairs[0] == [0]
        assert pairs[1] == [1]
        assert 2 not in pairs  # 60-70 not within 65-68

    def test_any_mode(self):
        pairs = _pairs_by_read(overlap_pairs(_reads(), _features(), "any"))
        assert pairs[0] == [0]
        assert pairs[1] == [1]
        assert pairs[2] == [2]  # 60-70 any-overlaps 65-68
        assert pairs[3] == [3]

    def test_reverse_within(self):
        # read 3 (10-20, chr2) contains feature 3 (12-18, chr2) fully.
        pairs = _pairs_by_read(overlap_pairs(_reads(), _features(), "reverse_within"))
        assert pairs[3] == [3]
        assert 0 not in pairs  # read 10-20 does not contain feat 5-25

    def test_strand_respected_by_default(self):
        # add a minus-strand read overlapping a plus-strand feature.
        reads = _reads().copy()
        reads = pd.concat(
            [
                reads,
                pd.DataFrame(
                    [{"chrom": "chr1", "start": 6, "end": 14, "strand": "-", "read_idx": 9}]
                ),
            ],
            ignore_index=True,
        )
        pairs = overlap_pairs(reads, _features(), "any")
        read_ids = {int(r) for r, _ in pairs}
        assert 9 not in read_ids
        pairs_ignore = overlap_pairs(reads, _features(), "any", ignore_strand=True)
        assert 9 in {int(r) for r, _ in pairs_ignore}

    def test_no_common_chromosomes(self):
        reads = _reads()[_reads()["chrom"] == "chr1"]
        features = _features()[_features()["chrom"] == "chr2"]
        assert len(overlap_pairs(reads, features, "any")) == 0

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError, match="unknown overlap mode"):
            overlap_pairs(_reads(), _features(), "bogus")

    def test_sorted_by_read_then_feature(self):
        pairs = overlap_pairs(_reads(), _features(), "any")
        assert (pairs[:, 0] == np.sort(pairs[:, 0])).all()

    def test_dtype_int64(self):
        pairs = overlap_pairs(_reads(), _features(), "within")
        assert pairs.dtype == np.int64
