"""Testes unitarios do modulo modules.motif_search.

Cobrem busca IUPAC, overlaps, fita reversa e resultado vazio. Nenhum teste faz
chamadas de rede.
"""

import pytest

from modules import motif_search


def test_find_motif_overlapping_forward():
    hits = motif_search.find_motif("AAAA", "AA")
    starts = [item["start"] for item in hits if item["strand"] == "+"]
    assert starts == [0, 1, 2]


def test_find_motif_reverse_strand():
    hits = motif_search.find_motif("GAATTC", "GAATTC")
    strands = {item["strand"] for item in hits}
    assert "+" in strands
    assert "-" in strands


def test_dickerson_ecoRI_is_one_locus_on_two_strands():
    hits = motif_search.find_motif("CGCGAATTCGCG", "GAATTC")
    summary = motif_search.summarize_motif_hits(hits)
    assert summary["literal_hits"] == 2
    assert summary["forward_hits"] == 1
    assert summary["reverse_complement_hits"] == 1
    assert summary["unique_loci"] == 1
    assert {(item["start"], item["end"]) for item in hits} == {(3, 9)}
    assert "CGCGAATTCGCG"[3:9] == "GAATTC"


def test_find_motif_iupac():
    hits = motif_search.find_motif("GAATTC", "GAATTY")
    assert any(item["strand"] == "+" and item["start"] == 0 for item in hits)


def test_find_motif_empty_result_is_empty_list():
    assert motif_search.find_motif("ACGT", "GGGG") == []


def test_find_motif_rejects_empty_inputs():
    with pytest.raises(ValueError):
        motif_search.find_motif("", "AA")
    with pytest.raises(ValueError):
        motif_search.find_motif("ACGT", "   ")


def test_find_motif_rna_uses_uracil():
    hits = motif_search.find_motif_rna("AUGCAUGC", "AUG")
    forward = [item for item in hits if item["strand"] == "+"]
    assert [item["start"] for item in forward] == [0, 4]
    assert all(item["match"] == "AUG" for item in forward)


def test_find_motif_protein_exact_and_x():
    hits = motif_search.find_motif_protein("MKVLTESTPEP", "TEST")
    assert hits == [
        {"start": 4, "end": 8, "strand": "+", "match": "TEST"},
    ]
    wildcard = motif_search.find_motif_protein("MKVL", "MXVL")
    assert wildcard[0]["match"] == "MKVL"
    with pytest.raises(ValueError):
        motif_search.find_motif_protein("MKVL", "M*VL")


def test_find_motif_hit_cap_raises_instead_of_partial_list():
    sequence = "A" * (motif_search.MAX_MOTIF_HITS + 1)
    with pytest.raises(ValueError, match="exceeded"):
        motif_search.find_motif(sequence, "A")


def test_parse_and_find_multiple_motifs():
    patterns = motif_search.parse_motif_patterns("AA; GC\nAT")
    assert patterns == ["AA", "GC", "AT"]
    hits = motif_search.find_motifs("AAGC", patterns, "DNA")
    assert all("pattern" in item for item in hits)
    plus = [item for item in hits if item["strand"] == "+"]
    assert any(item["pattern"] == "AA" and item["start"] == 0 for item in plus)
    filtered = motif_search.filter_motif_hits(hits, strand="+")
    assert all(item["strand"] == "+" for item in filtered)
    protein = motif_search.find_motifs("MKVLTESTPEP", ["TEST"], "PROTEIN")
    assert protein[0]["match"] == "TEST"
    assert protein[0]["pattern"] == "TEST"
    with pytest.raises(ValueError):
        motif_search.parse_motif_patterns("   ")
    too_many = ",".join(["AA"] * (motif_search.MAX_MOTIF_PATTERNS + 1))
    with pytest.raises(ValueError, match="At most"):
        motif_search.parse_motif_patterns(too_many)
    with pytest.raises(ValueError):
        motif_search.find_motifs("ACGT", ["AA"], "VIRUS")
