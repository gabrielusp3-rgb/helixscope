"""Property-based checks for sequence and alignment invariants."""

from __future__ import annotations

import math

import pytest

from modules import alignment, dna_analysis, explain, rna_analysis

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st


DNA_ALPHA = st.text(alphabet="ACGT", min_size=1, max_size=80)


@settings(max_examples=40, deadline=None)
@given(DNA_ALPHA)
def test_reverse_complement_involution(seq: str) -> None:
    twice = dna_analysis.reverse_complement(dna_analysis.reverse_complement(seq))
    assert twice == seq


@settings(max_examples=40, deadline=None)
@given(DNA_ALPHA)
def test_gc_at_percentages_use_canonical_denominator(seq: str) -> None:
    gc = dna_analysis.gc_content(seq)
    at = dna_analysis.at_content(seq) if hasattr(dna_analysis, "at_content") else None
    if math.isnan(gc):
        return
    assert 0.0 <= gc <= 100.0
    if at is not None and not math.isnan(at):
        assert abs((gc + at) - 100.0) < 0.15


@settings(max_examples=25, deadline=None)
@given(DNA_ALPHA, DNA_ALPHA)
def test_pairwise_alignment_rows_equal_length(seq_a: str, seq_b: str) -> None:
    result = alignment.pairwise_global(seq_a, seq_b)
    row_a = str(result.get("aligned_seq1") or "")
    row_b = str(result.get("aligned_seq2") or "")
    assert len(row_a) == len(row_b)
    identity = float(result.get("identity_pct") or 0)
    assert 0.0 <= identity <= 100.0


@settings(max_examples=20, deadline=None)
@given(st.sampled_from(["PREDICTED", "ILLUSTRATIVE", "UNAVAILABLE"]))
def test_status_invariants_in_explanations(status: str) -> None:
    if status == "PREDICTED":
        rec = explain.explain(
            "rna_fold",
            {"status": status, "mfe_kcal_mol": -1.0, "backend": "viennarna_python"},
        )
        assert rec.status == "PREDICTED"
        blob = (rec.plain_meaning + rec.limitations).lower()
        assert "not an experimental" in blob or "not an experimentally" in blob
    elif status == "ILLUSTRATIVE":
        rec = explain.explain("structure", {"kind": "illustrative", "structure_id": "x"})
        assert rec.status == "ILLUSTRATIVE"
        assert rec.status != "EXPERIMENTAL"
    else:
        rec = explain.explain("sasa", {"status": "UNAVAILABLE", "sasa_angstrom2": None})
        assert "successfully computed" not in rec.plain_meaning.lower()


def test_transcription_length_preserved() -> None:
    dna = "ACGTAC"
    rna = rna_analysis.transcribe(dna) if hasattr(rna_analysis, "transcribe") else dna.replace("T", "U")
    assert len(rna) == len(dna)


@settings(max_examples=25, deadline=None)
@given(DNA_ALPHA)
def test_transcription_alphabet_and_length(seq: str) -> None:
    rna = rna_analysis.transcribe(seq)
    assert len(rna) == len(seq)
    assert "T" not in rna
    assert set(rna) <= set("ACGU")
