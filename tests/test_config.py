"""Unit tests for :mod:`srna.config` (strategies and pipeline config)."""

from __future__ import annotations

import pytest

from srna.config import PipelineConfig, Strategy, canon_strategy


class TestStrategy:
    def test_overlap_rule(self):
        assert Strategy.FULLY_CONTAINED.overlap_rule == "within"
        assert Strategy.UNION.overlap_rule == "union"
        assert Strategy.ANY.overlap_rule == "any"
        assert Strategy.COMPARISON.overlap_rule == "within"

    def test_dir_name_uses_underscores(self):
        assert Strategy.FULLY_CONTAINED.dir_name == "fully_contained"
        assert Strategy.COMPARISON.dir_name == "comparison"

    def test_is_comparison(self):
        assert Strategy.COMPARISON.is_comparison
        assert not Strategy.FULLY_CONTAINED.is_comparison


class TestCanonStrategy:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("fully-contained", Strategy.FULLY_CONTAINED),
            ("fully_contained", Strategy.FULLY_CONTAINED),
            ("within", Strategy.FULLY_CONTAINED),
            ("fully", Strategy.FULLY_CONTAINED),
            ("UNION", Strategy.UNION),
            ("any", Strategy.ANY),
            ("comparison", Strategy.COMPARISON),
            ("compar", Strategy.COMPARISON),
            ("garbage", Strategy.FULLY_CONTAINED),
        ],
    )
    def test_canonicalization(self, raw, expected):
        assert canon_strategy(raw) is expected


class TestPipelineConfig:
    def test_invalid_chr_style_rejected(self):
        with pytest.raises(ValueError):
            PipelineConfig(chr_style="nope")  # type: ignore[arg-type]

    def test_from_env_defaults(self):
        cfg = PipelineConfig.from_env({})
        assert cfg.genome == "mm39"
        assert cfg.strategy is Strategy.FULLY_CONTAINED
        assert cfg.chr_style == "chr"
        assert cfg.force_rebuild_db is False

    def test_from_env_parses_values(self):
        cfg = PipelineConfig.from_env(
            {
                "SMALLRNA_GENOME": "hg38",
                "SMALLRNA_STRATEGY": "any",
                "SMALLRNA_CHR_STYLE": "none",
                "SMALLRNA_FORCE_REBUILD_DB": "1",
            }
        )
        assert cfg.genome == "hg38"
        assert cfg.strategy is Strategy.ANY
        assert cfg.chr_style == ""
        assert cfg.force_rebuild_db is True

    def test_for_substrategy_preserves_config(self):
        base = PipelineConfig(genome="mm39", strategy=Strategy.COMPARISON)
        sub = PipelineConfig.for_substrategy(base, Strategy.UNION)
        assert sub.strategy is Strategy.UNION
        assert sub.genome == base.genome
        assert sub.chr_style == base.chr_style
