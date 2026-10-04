"""Testes da referencia CRISPR fornecida (FASTA, hashes, duplicatas, gzip)."""

import gzip

import pytest

from modules import crispr_reference


def test_load_reference_fasta_computes_hash_and_rejects_filename_as_assembly():
    fasta = ">locus1\nGAGTCCGAGCAGAAGAAGAAGGG\n"
    loaded = crispr_reference.load_reference_from_text(
        fasta,
        source="user FASTA",
        upload_filename="hg38.fa",
    )
    assert loaded["genome_wide"] is False
    assert loaded["filename_trusted_as_assembly"] is False
    assert loaded["assembly_declared"] == ""
    assert loaded["n_contigs"] == 1
    assert loaded["identity_hash"]
    assert loaded["contigs"][0]["sequence"] == "GAGTCCGAGCAGAAGAAGAAGGG"
    assert "hg38" not in loaded["assembly_declared"]
    assert any("filename" in warning.lower() for warning in loaded["warnings"])


def test_duplicate_contig_identifiers_are_not_merged():
    fasta = ">dup\nATGCATGCATGC\n>dup\nGGCCGGCCGGCC\n"
    with pytest.raises(crispr_reference.ReferenceError, match="Duplicate") as exc:
        crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    assert exc.value.category == "PARSING_ERROR"


def test_empty_header_in_multifasta_is_parsing_error():
    fasta = ">ok\nATGCATGCATGC\n>\nGGCCGGCCGGCC\n"
    with pytest.raises(crispr_reference.ReferenceError, match="empty identifier") as exc:
        crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    assert exc.value.category == "PARSING_ERROR"


def test_rna_symbols_are_rejected_as_dna_reference():
    fasta = ">rna\nACGUACGUACGUACGU\n"
    with pytest.raises(crispr_reference.ReferenceError, match="non-DNA") as exc:
        crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    assert exc.value.category == "INVALID_INPUT"


def test_corrupt_fasta_header_without_records():
    with pytest.raises(crispr_reference.ReferenceError) as exc:
        crispr_reference.load_reference_from_text(">", source="user FASTA")
    assert exc.value.category in {"PARSING_ERROR", "INVALID_INPUT"}


def test_same_content_same_hash_different_declared_assembly_changes_cache_key():
    fasta = ">c1\nATGCATGCATGCATGCATGC\n"
    a = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", assembly="asmA"
    )
    b = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", assembly="asmB"
    )
    assert a["identity_hash"] == b["identity_hash"]
    assert crispr_reference.reference_cache_key(a) != crispr_reference.reference_cache_key(
        b
    )


def test_oversized_reference_is_resource_limit():
    huge = ">c\n" + ("A" * (crispr_reference.MAX_REFERENCE_NT + 1)) + "\n"
    with pytest.raises(crispr_reference.ReferenceError) as exc:
        crispr_reference.load_reference_from_text(huge, source="user FASTA")
    assert exc.value.category == "RESOURCE_LIMIT"


def test_gzip_roundtrip_and_decompression_bomb():
    fasta = b">c1\nATGCATGCATGCATGCATGC\n"
    compressed = gzip.compress(fasta)
    loaded = crispr_reference.load_reference_from_bytes(
        compressed, source="user FASTA", upload_filename="locus.fa.gz"
    )
    assert loaded["n_contigs"] == 1
    bomb = gzip.compress(b"A" * (crispr_reference.MAX_UNCOMPRESSED_BYTES + 1))
    with pytest.raises(crispr_reference.ReferenceError) as exc:
        crispr_reference.decode_reference_payload(bomb)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_zip_archive_is_rejected():
    with pytest.raises(crispr_reference.ReferenceError, match="ZIP") as exc:
        crispr_reference.decode_reference_payload(b"PK\x03\x04" + b"\x00" * 20)
    assert exc.value.category == "INVALID_INPUT"


def test_ncbi_protein_record_is_not_a_dna_reference():
    record = {
        "sequence_available": True,
        "sequence": "MKTFF",
        "database": "protein",
        "molecule": "protein",
        "accession": "P12345",
    }
    with pytest.raises(crispr_reference.ReferenceError, match="Protein") as exc:
        crispr_reference.reference_from_ncbi_record(record)
    assert exc.value.category == "INVALID_INPUT"


def test_ncbi_record_without_sequence_is_invalid():
    record = {
        "sequence_available": False,
        "sequence": "",
        "database": "nucleotide",
        "accession": "NC_000001",
    }
    with pytest.raises(crispr_reference.ReferenceError) as exc:
        crispr_reference.reference_from_ncbi_record(record)
    assert exc.value.category == "INVALID_INPUT"


def test_workspace_dna_becomes_single_contig_not_genome_wide():
    loaded = crispr_reference.reference_from_workspace_dna(
        "GAGTCCGAGCAGAAGAAGAAGGG", identifier="workspace"
    )
    assert loaded["scope"] == crispr_reference.SEARCH_SCOPE_WORKSPACE
    assert loaded["genome_wide"] is False
    assert loaded["n_contigs"] == 1
