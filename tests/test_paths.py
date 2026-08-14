"""Unit tests for :mod:`srna.paths` (sample labels and directory layout)."""

from __future__ import annotations

from pathlib import Path

import pytest

from srna.config import PipelineConfig, Strategy
from srna.paths import PathResolver, derive_sample_label


class TestDeriveSampleLabel:
    @pytest.mark.parametrize(
        ("name", "expected"),
        [
            ("Cumulus-cells.bwa.bam", "Cumulus-cells"),
            ("Granulosa-cells.bowtie2.bam", "Granulosa-cells"),
            ("Sample.bam", "Sample"),
            ("plain.txt", "plain.txt"),
        ],
    )
    def test_labels(self, name, expected):
        assert derive_sample_label(name) == expected


class TestPathResolver:
    def _resolver(self, strategy: Strategy) -> PathResolver:
        return PathResolver(
            root=Path("/tmp/root"),
            config=PipelineConfig(genome="mm39", strategy=strategy),
        )

    def test_input_dirs(self):
        r = self._resolver(Strategy.FULLY_CONTAINED)
        assert r.bam_dir == Path("/tmp/root/data/bam")
        assert r.db_raw_dir == Path("/tmp/root/data/DB")
        assert r.db_cache_dir == Path("/tmp/root/data/cache/mm39")

    def test_output_base_comparison(self):
        assert self._resolver(Strategy.FULLY_CONTAINED).output_base == Path("/tmp/root/output")
        assert self._resolver(Strategy.COMPARISON).output_base == Path(
            "/tmp/root/output/comparison"
        )

    def test_strategy_dir(self, tmp_path):
        r = PathResolver(root=tmp_path, config=PipelineConfig())
        assert r.strategy_dir(Strategy.FULLY_CONTAINED).name == "fully_contained"
        assert r.strategy_dir(Strategy.ANY).name == "any"

    def test_bam_file_raises_when_missing(self, tmp_path):
        r = PathResolver(root=tmp_path, config=PipelineConfig())
        with pytest.raises(FileNotFoundError):
            r.bam_file("Sample")

    def test_bam_file_picks_aligner_suffix(self, tmp_path):
        (tmp_path / "data" / "bam").mkdir(parents=True)
        bam = tmp_path / "data" / "bam" / "Sample.bwa.bam"
        bam.write_bytes(b"")
        r = PathResolver(root=tmp_path, config=PipelineConfig())
        assert r.bam_file("Sample") == bam

    def test_samples_detection(self, tmp_path):
        bam_dir = tmp_path / "data" / "bam"
        bam_dir.mkdir(parents=True)
        (bam_dir / "A.bwa.bam").write_bytes(b"")
        (bam_dir / "B.bam").write_bytes(b"")
        r = PathResolver(root=tmp_path, config=PipelineConfig())
        assert r.samples() == ["A", "B"]

    def test_require_raw_db_files(self, tmp_path):
        r = PathResolver(root=tmp_path, config=PipelineConfig())
        with pytest.raises(FileNotFoundError):
            r.require_raw_db_files()
