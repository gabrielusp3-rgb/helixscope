"""Pasted-sequence contract: FASTA headers never reach molecule alphabets."""

from __future__ import annotations

import pytest

from helixscope_core.alignment import pairwise_align
from helixscope_core.crispr import design_guides
from helixscope_core.dna import analyze_dna
from helixscope_core.motif import search_motifs
from helixscope_core.protein import analyze_protein
from helixscope_core.rna import analyze_rna
from helixscope_core.sequence_input import resolve_pasted_sequence


INSULIN_FASTA = (
    ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n"
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n"
)

DNA_FASTA = ">dna_fixture\nATGCATGCATGC\n"


def test_resolve_strips_single_fasta_header() -> None:
    parsed = resolve_pasted_sequence(DNA_FASTA)
    assert parsed["format"] == "fasta"
    assert parsed["sequence"] == "ATGCATGCATGC"
    assert ">" not in parsed["sequence"]
    assert parsed["identifier"]


def test_resolve_rejects_multifasta() -> None:
    with pytest.raises(ValueError, match="records"):
        resolve_pasted_sequence(">one\nATGC\n>two\nGGCC\n")


def test_analyze_dna_fasta_computed() -> None:
    result = analyze_dna(DNA_FASTA)
    assert result["status"] == "COMPUTED"
    assert result["length"] == 12
    assert result["gc_content"] == pytest.approx(50.0)
    assert result["input_format"] == "fasta"
    assert ">" not in str(result["sequence"])


def test_analyze_dna_invalid_alphabet_stays_error() -> None:
    result = analyze_dna("NNNN")
    assert result["status"] == "ERROR"
    assert result["validation"]["is_valid"] is False


def test_analyze_protein_fasta_does_not_see_punctuation() -> None:
    result = analyze_protein(INSULIN_FASTA)
    assert result["status"] == "COMPUTED"
    assert result["length"] == 110
    assert result["input_format"] == "fasta"
    assert ">" not in str(result["sequence"])
    assert "|" not in str(result["sequence"])


def test_analyze_protein_reports_real_invalid_residues_only() -> None:
    with pytest.raises(ValueError, match="O, U") as caught:
        analyze_protein("MKTAYIOUAK")
    message = str(caught.value)
    assert "Found unsupported residues: O, U." in message
    assert ">" not in message
    assert "|" not in message


def test_analyze_protein_header_punctuation_is_not_invalid_after_parse() -> None:
    with pytest.raises(ValueError) as caught:
        analyze_protein("MKWVTFISLLLLFSSAYSROU")
    message = str(caught.value)
    assert "O" in message and "U" in message
    assert ">" not in message
    assert "|" not in message


def test_analyze_rna_fasta() -> None:
    result = analyze_rna(">rna\nAUGCAUGC\n")
    assert result["status"] == "COMPUTED"
    assert result["length"] == 8
    assert result["gc_content"] == pytest.approx(50.0)
    assert result["au_content"] == pytest.approx(50.0)


def test_alignment_and_motif_accept_fasta() -> None:
    aligned = pairwise_align(">a\nACGTACGT\n", ">b\nACGTACGG\n")
    assert aligned["status"] in {"COMPUTED", "HEURISTIC"} or "aligned_seq1" in aligned
    hits = search_motifs(">m\nACGTACGTACGT\n", "ACGT")
    assert hits


def test_crispr_fasta_target_yields_guides() -> None:
    guides = design_guides(">target\nACGTACGTACGTACGTACGTAGG\n", "SpCas9")
    assert guides
    assert guides[0]["pam_sequence"]
