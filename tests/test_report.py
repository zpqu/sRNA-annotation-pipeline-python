"""Unit tests for step-10 report helpers (:mod:`srna.summary.pipeline_report`)."""

from __future__ import annotations

from srna.summary.pipeline_report import describe


class TestDescribe:
    def test_known_files(self):
        assert "sample summary" in describe("Table_01a_sample_summary.csv").lower()
        assert "read size distribution" in describe("Table_01b_read_size_distribution.csv").lower()
        assert (
            "annotation composition"
            in describe("Table_02a_annotation_count_unique_reads.csv").lower()
        )
        assert "abundance statistics" in describe("Table_03a_category_abundance_summary.csv")
        assert "sanity-check" in describe("Table_10_sanity_checks.csv")
        assert (
            "read size / count distribution" in describe("Figure_01.read_size_vs_count.png").lower()
        )
        assert "lorenz curves" in describe("Figure_03c_lorenz_matmiRNA.pdf").lower()

    def test_unknown_file_defaults(self):
        assert describe("mystery_file.txt") == "See producing step for details"

    def test_pirna_overlap_summary(self):
        assert "piRNA reads" in describe("Table_06_piRNA_overlap_summary.csv")

    def test_s01_tables(self):
        assert "overlap-rule" in describe("Table_s01a_overlap_rule_composition.csv").lower()
        assert "read movement" in describe("Table_s01c_read_movement_contained_vs_any.csv").lower()

    def test_s01_figures(self):
        assert "overlap-rule" in describe("Figure_s01a_overlap_rules_composition.png").lower()
