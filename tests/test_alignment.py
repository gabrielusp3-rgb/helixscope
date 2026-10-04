"""Testes unitarios do modulo modules.alignment.

Cobrem alinhamento global e local com os parametros reais do modulo e a matriz
de dotplot. Nenhum teste faz chamadas de rede.
"""

import math

import pytest

from modules import alignment


def test_pairwise_global_known_score():
    result = alignment.pairwise_global("ACGT", "AGT")
    assert result["score"] == pytest.approx(1.0)
    assert result["gaps"] >= 1
    assert result["identity_pct"] == pytest.approx(75.0)
    assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])


def test_pairwise_global_identical():
    result = alignment.pairwise_global("ACGTACGT", "ACGTACGT")
    assert result["score"] == pytest.approx(16.0)
    assert result["identity_pct"] == pytest.approx(100.0)
    assert result["gaps"] == 0


def test_pairwise_local_finds_substring():
    result = alignment.pairwise_local("ACGTACGT", "CGT")
    assert result["score"] > 0
    assert "CGT" in result["aligned_seq1"] or "CGT" in result["aligned_seq2"]


def test_pairwise_rejects_empty():
    with pytest.raises(ValueError):
        alignment.pairwise_global("", "ACGT")
    with pytest.raises(ValueError):
        alignment.pairwise_local("ACGT", "   ")


def test_pairwise_local_zero_when_no_positive_alignment():
    result = alignment.pairwise_local("AAAA", "TTTT")
    assert result["score"] == pytest.approx(0.0)
    assert result["aligned_seq1"] == ""
    assert result["aligned_seq2"] == ""


def test_pairwise_rejects_mixed_molecule_types():
    with pytest.raises(ValueError, match="Cannot align"):
        alignment.pairwise_global("ACGTACGT", "ACGUACGU")
    with pytest.raises(ValueError, match="Cannot align"):
        alignment.pairwise_local("ACGTACGT", "MKLPQERT")


def _affine_global_score(
    seq_a: str,
    seq_b: str,
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -5.0,
    gap_extend: float = -0.5,
) -> float:
    """Needleman-Wunsch-Gotoh independente, so para conferir o score."""
    negative = -1e12
    n = len(seq_a)
    m = len(seq_b)
    match_m = [[negative] * (m + 1) for _ in range(n + 1)]
    gap_x = [[negative] * (m + 1) for _ in range(n + 1)]
    gap_y = [[negative] * (m + 1) for _ in range(n + 1)]
    match_m[0][0] = 0.0
    for i in range(1, n + 1):
        gap_x[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, m + 1):
        gap_y[0][j] = gap_open + (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            match_m[i][j] = pair + max(
                match_m[i - 1][j - 1], gap_x[i - 1][j - 1], gap_y[i - 1][j - 1]
            )
            gap_x[i][j] = max(match_m[i - 1][j] + gap_open, gap_x[i - 1][j] + gap_extend)
            gap_y[i][j] = max(match_m[i][j - 1] + gap_open, gap_y[i][j - 1] + gap_extend)
    return max(match_m[n][m], gap_x[n][m], gap_y[n][m])


@pytest.mark.parametrize(
    "seq_a,seq_b",
    [
        ("ACGT", "AGT"),
        ("ACGTACGT", "ACGTACGT"),
        ("ATGAAATTTGGG", "ATGAAATTCGGG"),
        ("GATTACA", "GCATGCT"),
        ("AAACCCGGGTTTAAACCCGGGTTT", "AAACCCGGTTTTAACCCGGGTTT"),
    ],
)
def test_pairwise_global_matches_independent_gotoh(seq_a, seq_b):
    expected = _affine_global_score(seq_a, seq_b)
    result = alignment.pairwise_global(seq_a, seq_b)
    assert result["score"] == pytest.approx(expected)
    assert len(result["aligned_seq1"]) == len(result["aligned_seq2"])


def test_pairwise_local_real_coding_and_protein_fragments():
    dna = alignment.pairwise_local(
        "ATGAAATTTGGGAAACCCGGGTTT",
        "AAATTTGGGAAA",
    )
    assert dna["score"] == pytest.approx(24.0)
    assert "AAATTTGGGAAA" in dna["aligned_seq1"] or "AAATTTGGGAAA" in dna["aligned_seq2"]
    protein = alignment.pairwise_local("MKVLTESTPEPTIDE", "TESTPEP")
    assert protein["score"] > 0
    assert protein["scoring"] == "BLOSUM62"
    assert protein["molecule"] == "PROTEIN"
    assert "TEST" in protein["aligned_seq1"] or "TEST" in protein["aligned_seq2"]
    matrix = alignment.dotplot_matrix("ACGT", "ACGT")
    assert matrix.shape == (4, 4)
    assert int(matrix[0, 0]) == 1
    assert int(matrix[0, 1]) == 0


def test_pairwise_protein_identical_matches_blosum62_diagonal():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    from Bio.Align import substitution_matrices

    blosum = substitution_matrices.load("BLOSUM62")
    expected = sum(float(blosum[residue, residue]) for residue in seq)
    result = alignment.pairwise_global(seq, seq)
    assert result["scoring"] == "BLOSUM62"
    assert result["identity_pct"] == pytest.approx(100.0)
    assert result["score"] == pytest.approx(expected)
    assert result["status"] == "COMPUTED"


def test_format_pairwise_export_and_msa_unavailable():
    result = alignment.pairwise_global("ACGTACGT", "ACGTACGT")
    text = alignment.format_pairwise_export(result, "alpha", "beta", width=8)
    assert "BLOSUM62" not in text
    assert "identity match=2 mismatch=-1" in text
    assert "ACGTACGT" in text
    assert "ungapped_identity_pct=" in text
    assert "coverage_seq1_pct=" in text
    availability = alignment.multiple_alignment_availability()
    assert availability["available"] is True
    assert availability["blast_is_not_msa"] is True
    assert "Needleman" in availability["reason"] or "not MSA" in availability["reason"]
    blast = alignment.blast_search_availability()
    assert blast["available"] is False
    assert "BLAST" in blast["reason"]
    with pytest.raises(ValueError):
        alignment.format_pairwise_export({"aligned_seq1": "", "aligned_seq2": ""})


def test_alignment_coverage_and_ungapped_identity_separate_from_identity():
    identical = alignment.pairwise_global("ACGTACGT", "ACGTACGT")
    assert identical["identity_pct"] == pytest.approx(100.0)
    assert identical["ungapped_identity_pct"] == pytest.approx(100.0)
    assert identical["coverage_seq1_pct"] == pytest.approx(100.0)
    assert identical["coverage_seq2_pct"] == pytest.approx(100.0)
    assert identical["method"] == "Needleman-Wunsch"
    gapped = alignment.pairwise_global("ACGT", "AGT")
    assert gapped["identity_pct"] == pytest.approx(75.0)
    assert gapped["gaps"] >= 1
    ungapped_denom = gapped["matches"] + gapped["mismatches"]
    expected_ungapped = round((gapped["matches"] / ungapped_denom) * 100.0, 2)
    assert gapped["ungapped_identity_pct"] == pytest.approx(expected_ungapped)
    assert gapped["coverage_seq1_pct"] == pytest.approx(100.0)
    assert gapped["coverage_seq2_pct"] == pytest.approx(100.0)
    empty_local = alignment.pairwise_local("AAAA", "TTTT")
    assert empty_local["identity_pct"] == pytest.approx(0.0)
    assert math.isnan(empty_local["ungapped_identity_pct"])
    assert empty_local["coverage_seq1_pct"] == pytest.approx(0.0)
    classes = alignment.classify_alignment_columns("AC-GT", "ACAGT")
    assert classes == ["match", "match", "gap", "match", "match"]
    with pytest.raises(ValueError):
        alignment.classify_alignment_columns("AC", "AGT")
