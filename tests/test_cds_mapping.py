"""CDS to codon/residue mapping. Ambiguous transcripts preserved."""

from __future__ import annotations

import os

from modules import cds_mapping, genome_annotation

FIXTURE_GFF = os.path.join(
    os.path.dirname(__file__), "fixtures", "helixscope_test_reference.gff3"
)


def test_cds_mapping_plus_strand_first_codon():
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    parsed = genome_annotation.parse_gff_text(text, declared_assembly="HELIXSCOPE_TEST_REF")
    result = cds_mapping.map_genomic_base_to_codon(
        parsed["features"], contig="HS_TEST_1", position_0based=40
    )
    assert result["unmapped"] is False
    assert result["effect"] is None
    assert result["n_mappings"] >= 1
    hit = result["mappings"][0]
    assert hit["codon_index_0based"] == 0
    assert hit["position_in_codon_0based"] == 0
    assert hit["protein_residue_number_1based"] == 1
    assert hit["transcript_id"] == "tx1"


def test_cds_mapping_unmapped_intergenic():
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    parsed = genome_annotation.parse_gff_text(text, declared_assembly="HELIXSCOPE_TEST_REF")
    result = cds_mapping.map_genomic_base_to_codon(
        parsed["features"], contig="HS_TEST_1", position_0based=2
    )
    assert result["unmapped"] is True
    assert result["mappings"] == []


def test_cds_mapping_preserves_multiple_transcripts():
    text = open(FIXTURE_GFF, encoding="utf-8").read()
    parsed = genome_annotation.parse_gff_text(text, declared_assembly="HELIXSCOPE_TEST_REF")
    result = cds_mapping.map_genomic_base_to_codon(
        parsed["features"], contig="HS_TEST_2", position_0based=20
    )
    tx_ids = {row["transcript_id"] for row in result["mappings"]}
    assert "tx2" in tx_ids
    assert "tx3" in tx_ids
