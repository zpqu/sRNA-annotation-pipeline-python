"""Unit tests for the feature-DB parsers (:mod:`srna.features`)."""

from __future__ import annotations

import numpy as np
import pyarrow as pa

from srna.features.mirna import parse_mirna_gff3
from srna.features.pirna import parse_pirna
from srna.features.refgene import parse_refgene_gtf
from srna.features.trna_rm import parse_rm_bed, parse_trna_bed


class TestParseMirna:
    def test_mature_and_primary_tables(self, tmp_path):
        p = tmp_path / "miRNA.gff3"
        p.write_text(
            "##gff-version 3\n"
            "chr1\tsrc\tmiRNA\t1\t22\t.\t+\t.\tID=mir-1;Name=mmu-miR-1-5p;\n"
            "chr1\tsrc\tmiRNA_primary_transcript\t1\t50\t.\t-\t.\tID=pri-1;Name=pri-mir-1;\n"
        )
        tables = parse_mirna_gff3(p)
        mat = tables["matmiRNA"]
        assert len(mat) == 1
        assert mat.loc[0, "feature_id"] == "mmu-miR-1-5p"
        assert mat.loc[0, "gene_id"] == "mir-1"
        assert mat.loc[0, "source"] == "miRBase"
        assert mat.loc[0, "start"] == 0
        pri = tables["primiRNA"]
        assert len(pri) == 1
        assert pri.loc[0, "gene_id"] == "pri-1"
        assert pri.loc[0, "strand"] == "-"

    def test_empty_file(self, tmp_path):
        p = tmp_path / "miRNA.gff3"
        p.write_text("##gff-version 3\n")
        tables = parse_mirna_gff3(p)
        assert len(tables["matmiRNA"]) == 0
        assert len(tables["primiRNA"]) == 0


class TestParseTrnaRm:
    def test_trna_bed12(self, tmp_path):
        p = tmp_path / "tRNAs.bed"
        p.write_text("chr1\t100\t150\ttRNA-A\t0\t+\t100\t150\t0\t1\t50\t0,\n")
        tab = parse_trna_bed(p, "mm39")
        assert tab.loc[0, "feature_type"] == "tRNA"
        assert tab.loc[0, "feature_id"] == "tRNA-A"
        assert tab.loc[0, "source"] == "mm39_tRNAs"
        assert tab.loc[0, "start"] == 100
        assert tab.loc[0, "end"] == 150

    def test_rm_bed6(self, tmp_path):
        p = tmp_path / "RM.bed"
        p.write_text("chr1\t200\t300\tL1\tabc\t-\n")
        tab = parse_rm_bed(p, "mm39")
        assert tab.loc[0, "feature_type"] == "RM"
        assert tab.loc[0, "feature_id"] == "L1"
        assert tab.loc[0, "source"] == "mm39_rmsk"
        assert np.isnan(tab.loc[0, "score"])

    def test_empty_bed(self, tmp_path):
        p = tmp_path / "RM.bed"
        p.write_text("# nothing\n")
        tab = parse_rm_bed(p, "mm39")
        assert len(tab) == 0


class TestParsePirna:
    def _write(self, tmp_path):
        bed = tmp_path / "piRBase.bed"
        bed.write_text(
            "chr1\t100\t120\tpiR-1\t10\t+\n"
            "chr2\t200\t230\tpiR-2\t5\t-\n"
            "chr1_alt\t1\t10\tpiR-alt\t0\t+\n"
        )
        gtf = tmp_path / "piRNAdb.gtf"
        gtf.write_text(
            'chr1\tsrc\tpiRNA\t101\t120\t.\t+\t.\tpiRNA_code "code-dup";\n'
            'chr1\tsrc\tpiRNA\t130\t150\t.\t+\t.\tpiRNA_code "code-new";\n'
            'chr1\tsrc\tpiRNA\t101\t120\t.\t+\t.\tpiRNA_code "code-dup2";\n'
        )
        return bed, gtf

    def test_merge_order_and_dedup(self, tmp_path):
        bed, gtf = self._write(tmp_path)
        tab = parse_pirna(bed, gtf, "mm39")
        assert isinstance(tab, pa.Table)
        assert tab.column_names == [
            "chrom",
            "start",
            "end",
            "strand",
            "feature_type",
            "feature_id",
            "gene_id",
            "gene_name",
            "source",
            "score",
        ]
        assert tab.num_rows == 3  # 2 piRBase + 1 piRNAdb
        sources = tab["source"].to_pylist()
        assert sources[:2] == ["piRBase", "piRBase"]
        assert sources[2] == "piRNAdb"
        assert tab["feature_id"].to_pylist() == ["piR-1", "piR-2", "code-new"]
        # alt scaffold dropped, and the piRNAdb dup/matching rows dropped
        assert "piR-alt" not in tab["feature_id"].to_pylist()
        assert "code-dup" not in tab["feature_id"].to_pylist()

    def test_score_kept_for_pirbase(self, tmp_path):
        bed, gtf = self._write(tmp_path)
        tab = parse_pirna(bed, gtf, "mm39")
        assert tab["score"].to_pylist()[:2] == [10.0, 5.0]


class TestParseRefgene:
    def _write(self, tmp_path):
        p = tmp_path / "refGene.gtf"

        def attr(tx: str) -> str:
            return f'gene_id "{tx}"; transcript_id "NM_1"; gene_name "{tx}";'

        def attr_nr(tx: str, nr: str) -> str:
            return f'gene_id "{tx}"; transcript_id "{nr}"; gene_name "{tx}";'

        gtf = (
            f"chr1\tsrc\texon\t100\t150\t.\t+\t.\t{attr('Xkr4')}\n"
            f"chr1\tsrc\texon\t200\t250\t.\t+\t.\t{attr('Xkr4')}\n"
            f"chr1\tsrc\tCDS\t140\t250\t.\t+\t.\t{attr('Xkr4')}\n"
            f"chr1\tsrc\t5UTR\t100\t139\t.\t+\t.\t{attr('Xkr4')}\n"
            f"chr1\tsrc\t3UTR\t249\t260\t.\t+\t.\t{attr('Xkr4')}\n"
            f"chr1\tsrc\texon\t500\t560\t.\t-\t.\t{attr_nr('SnordX', 'NR_1')}\n"
            f"chr1\tsrc\texon\t700\t750\t.\t+\t.\t{attr_nr('Gm1', 'NR_2')}\n"
        )
        p.write_text(gtf)
        return p

    def test_nm_tables(self, tmp_path):
        tables = parse_refgene_gtf(self._write(tmp_path), "mm39")
        assert len(tables["NM.exon"]) == 2
        assert tables["NM.CDS"].loc[0, "start"] == 139  # 0-based
        assert tables["NM.5UTR"].loc[0, "start"] == 99
        assert tables["NM.3UTR"].loc[0, "start"] == 248
        assert tables["NM.mRNA"].loc[0, "start"] == 99
        assert tables["NM.mRNA"].loc[0, "end"] == 250
        # intron = gap between exon 100-150 and exon 200-250
        intron = tables["NM.intron"]
        assert len(intron) == 1
        assert intron.loc[0, "start"] == 150
        assert intron.loc[0, "end"] == 199

    def test_flanks_respect_strand(self, tmp_path):
        tables = parse_refgene_gtf(self._write(tmp_path), "mm39")
        up = tables["NM.up1k"]
        down = tables["NM.down1k"]
        plus = up[up["gene_name"] == "Xkr4"]
        assert plus.loc[0, "start"] == 99 - 1000
        assert plus.loc[0, "end"] == 99
        down_plus = down[down["gene_name"] == "Xkr4"]
        assert down_plus.loc[0, "start"] == 250
        assert down_plus.loc[0, "end"] == 1250

    def test_nr_snorna_vs_lincrna(self, tmp_path):
        tables = parse_refgene_gtf(self._write(tmp_path), "mm39")
        assert tables["snoRNA"]["gene_name"].tolist() == ["SnordX"]
        assert tables["lincRNA.exon"]["gene_name"].tolist() == ["Gm1"]
        assert tables["NR.exon"]["gene_name"].tolist() == ["SnordX", "Gm1"]
