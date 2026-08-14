"""Unit tests for raw annotation-file readers (:mod:`srna.features.readers`)."""

from __future__ import annotations

import pandas as pd

from srna.features.readers import (
    extract_gff_fields,
    extract_gtf_fields,
    read_bed6,
    read_bed12,
    read_gff3,
    read_gtf,
)


def test_read_bed6(tmp_path):
    p = tmp_path / "features.bed"
    p.write_text("chr1\t10\t20\tfeat1\t0\t+\nchr2\t30\t40\tfeat2\t1\t-\n")
    df = read_bed6(p)
    assert list(df.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert df.loc[0, "chrom"] == "chr1"
    assert df.loc[0, "start"] == 10
    assert df.loc[1, "strand"] == "-"


def test_read_bed6_skips_comments(tmp_path):
    p = tmp_path / "features.bed"
    p.write_text("# header\nchr1\t10\t20\tfeat1\t0\t+\n")
    df = read_bed6(p)
    assert len(df) == 1


def test_read_bed12_keeps_first_six_columns(tmp_path):
    p = tmp_path / "features.bed12"
    p.write_text("chr1\t0\t100\tg\t0\t+\t0\t100\t0\t1\t100\t0,\n")
    df = read_bed12(p)
    assert list(df.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert df.loc[0, "end"] == 100


def test_read_gtf_zero_based(tmp_path):
    p = tmp_path / "genes.gtf"
    p.write_text('chr1\tsrc\texon\t11\t20\t.\t+\t.\tgene_id "Xkr4"; transcript_id "NM_1";\n')
    df = read_gtf(p)
    assert list(df.columns) == [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    assert df.loc[0, "start"] == 10  # 1-based 11 -> 0-based 10
    assert df.loc[0, "end"] == 20
    assert df.loc[0, "feature"] == "exon"


def test_read_gff3(tmp_path):
    p = tmp_path / "mirna.gff3"
    p.write_text("##gff-version 3\nchr1\tsrc\tmiRNA\t1\t10\t.\t+\t.\tID=mir-1;\n")
    df = read_gff3(p)
    assert len(df) == 1
    assert df.loc[0, "start"] == 0


def test_extract_gtf_fields(tmp_path):
    p = tmp_path / "genes.gtf"
    p.write_text('chr1\ts\tgene\t1\t10\t.\t+\t.\tgene_id "A"; transcript_id "NM_1";\n')
    out = extract_gtf_fields(read_gtf(p), ["gene_id", "transcript_id"])
    assert out.loc[0, "gene_id"] == "A"
    assert out.loc[0, "transcript_id"] == "NM_1"


def test_extract_gff_fields(tmp_path):
    p = tmp_path / "mirna.gff3"
    p.write_text("chr1\ts\tmiRNA\t1\t10\t.\t+\t.\tID=mir-1;Name=mmu-miR-1-5p;\n")
    out = extract_gff_fields(read_gff3(p), ["ID", "Name"])
    assert out.loc[0, "ID"] == "mir-1"
    assert out.loc[0, "Name"] == "mmu-miR-1-5p"


def test_gtf_missing_field_is_empty(tmp_path):
    p = tmp_path / "genes.gtf"
    p.write_text('chr1\ts\tgene\t1\t10\t.\t+\t.\tgene_id "A";\n')
    out = extract_gtf_fields(read_gtf(p), ["transcript_id"])
    assert out.loc[0, "transcript_id"] == ""
    assert pd.isna(out.loc[0, "transcript_id"]) or out.loc[0, "transcript_id"] == ""
