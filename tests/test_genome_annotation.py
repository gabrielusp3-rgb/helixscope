"""GFF/GTF context, wrong-assembly reject, untrusted input."""

from __future__ import annotations

import os

import pytest

from modules import crispr_assemblies, genome_annotation, genome_store

FIXTURE_GFF = os.path.join(
    os.path.dirname(__file__), "fixtures", "helixscope_test_reference.gff3"
)


def test_parse_gff_and_classify_gene_exon_intron_intergenic():
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    parsed = genome_annotation.parse_gff_text(text, declared_assembly="HELIXSCOPE_TEST_REF")
    index = genome_annotation.build_interval_index(parsed["features"])
    cds = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_1", start_0based=40, end_0based=50
    )
    assert cds["primary_label"] == "CDS"
    assert cds["coding"] is True
    assert cds["effect"] is None
    intron = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_1", start_0based=55, end_0based=60
    )
    assert intron["primary_label"] == "intron"
    inter = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_1", start_0based=0, end_0based=5
    )
    assert inter["primary_label"] == "intergenic"
    utr = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_1", start_0based=80, end_0based=90
    )
    assert utr["primary_label"] == "UTR"


def test_negative_strand_and_multiple_overlaps():
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    parsed = genome_annotation.parse_gff_text(text, declared_assembly="HELIXSCOPE_TEST_REF")
    index = genome_annotation.build_interval_index(parsed["features"])
    minus = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_2", start_0based=0, end_0based=10
    )
    assert minus["genes"]
    assert any(g.get("strand") == "-" for g in minus["genes"])
    overlap = genome_annotation.classify_hit_context(
        index, contig="HS_TEST_2", start_0based=20, end_0based=30
    )
    ids = {g.get("gene_id") for g in overlap["genes"]}
    assert "HSG2" in ids
    assert "HSG3" in ids
    assert len(overlap["transcripts"]) >= 2
    assert overlap["transcript_selection"] == "all_overlapping"


def test_wrong_assembly_annotation_rejected():
    with pytest.raises(genome_annotation.AnnotationError) as exc:
        genome_annotation.bind_annotation_to_assembly("GRCh38.p14", "GRCh37")
    assert exc.value.category == "CORRUPT"
    with pytest.raises(genome_annotation.AnnotationError):
        genome_annotation.bind_annotation_to_assembly("GRCh38.p14", "hg19")
    genome_annotation.bind_annotation_to_assembly("GRCh38.p14", "GCF_000001405.40")
    genome_annotation.bind_annotation_to_assembly(
        crispr_assemblies.TEST_REFERENCE_ID, "HELIXSCOPE_TEST_REF"
    )


def test_malicious_seqid_and_oversize_line_rejected():
    with pytest.raises(genome_annotation.AnnotationError):
        genome_annotation.parse_gff_text(
            "chr1/../etc\ttestdb\tgene\t1\t10\t.\t+\t.\tID=x\n",
            declared_assembly="GRCh38.p14",
        )
    with pytest.raises(genome_annotation.AnnotationError):
        genome_annotation.parse_gff_text(
            "chr1\ttestdb\tgene\t1\t10\t.\t+\t.\t" + ("A" * 40000) + "\n",
            declared_assembly="GRCh38.p14",
        )


def test_install_test_annotation_and_annotate_hits(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIXSCOPE_REFERENCE_DIR", str(tmp_path / "refs"))
    monkeypatch.setenv("HELIXSCOPE_JOBS_DIR", str(tmp_path / "jobs"))
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    genome_store.install_test_reference()
    manifest = genome_annotation.install_annotation_from_text(
        crispr_assemblies.TEST_REFERENCE_ID,
        text,
        declared_assembly="HELIXSCOPE_TEST_REF",
    )
    assert manifest["status"] == "READY"
    loaded = genome_annotation.load_ready_annotation(crispr_assemblies.TEST_REFERENCE_ID)
    assert loaded is not None
    hits = [
        {
            "chromosome_or_contig": "HS_TEST_1",
            "start_0based": 40,
            "end_0based": 60,
            "site_role": "on_target",
        }
    ]
    annotated = genome_annotation.annotate_hits(
        hits,
        loaded["index"],
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        annotation_sha256=loaded["file_sha256"],
    )
    assert annotated[0]["feature_annotation"] == "CDS"
    assert annotated[0]["genomic_context"]["effect"] is None
    raw_identity = "search-key"
    ctx_a = genome_annotation.context_cache_identity(
        search_identity=raw_identity, annotation_sha256="aaa"
    )
    ctx_b = genome_annotation.context_cache_identity(
        search_identity=raw_identity, annotation_sha256="bbb"
    )
    assert ctx_a != ctx_b


def test_catalog_has_official_gff_md5():
    human = crispr_assemblies.catalog_by_id("GRCh38.p14")
    assert human["official_gff_md5"] == "24b731562b9d4cae9e37b23404b5be16"
    assert human["official_gtf_md5"] == "a561524a1ac438a2a95aed54d00e490e"
    assert "GCF_000001405.40" in human["official_gff_url"]
    urls = genome_annotation.annotation_catalog_urls("GRCh38.p14")
    assert urls["gff_md5"] == human["official_gff_md5"]
