"""Contrato de 100k, 125k e 150k. Os esperados sao contados neste arquivo."""

from __future__ import annotations

import hashlib
import math

from modules import dna_analysis, scale_contract, scale_profile


def _independent_counts(sequence: str) -> dict[str, int]:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for base in sequence:
        counts[base] = counts[base] + 1
    return counts


def _independent_gc_percent(counts: dict[str, int], length: int) -> float:
    return ((counts["G"] + counts["C"]) / length) * 100.0


def test_official_fixtures_keep_exact_length_and_composition() -> None:
    for length in scale_contract.FULL_SCALE_LENGTHS:
        sequence = scale_contract.deterministic_nucleotide_fixture(length)
        assert len(sequence) == length
        assert sequence[:4] == "ACGT"
        assert "N" not in sequence
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        assert scale_contract.fixture_digest(sequence) == digest
        counts = _independent_counts(sequence)
        assert sum(counts.values()) == length
        assert counts["A"] == counts["C"] == counts["G"] == counts["T"] == length // 4
        composition = dna_analysis.nucleotide_composition(sequence)
        for base, count in counts.items():
            assert composition[base]["count"] == count
        gc = dna_analysis.gc_content(sequence)
        assert gc == _independent_gc_percent(counts, length)
        assert dna_analysis.at_content(sequence) == 100.0 - gc
        entropy = dna_analysis.shannon_entropy(sequence)
        assert entropy == 2.0
        assert dna_analysis.gc_skew_global(sequence) == 0.0
        assert dna_analysis.at_skew(sequence) == 0.0
        cumulative = dna_analysis.cumulative_skew_summary(sequence)
        assert cumulative["input_length"] == length
        assert cumulative["analyzed_length"] == length
        assert cumulative["final_gc"] == 0
        assert cumulative["final_at"] == 0
        assert cumulative["visualized_length"] <= scale_profile.MAX_PLOT_POINTS + 1
        assert cumulative["visualized_length"] < length
        assert cumulative["series"][-1]["position"] == length
        kmers = dna_analysis.kmer_counts(sequence, k=3)
        assert sum(item["count"] for item in kmers.values()) == length - 2
        info = dna_analysis.validate_for_molecule(sequence, "DNA")
        assert info["is_valid"] is True
        assert info["length"] == length


def test_cumulative_plot_refuses_silent_full_draw_past_the_cap() -> None:
    sequence = "ACGT" * 800
    with __import__("pytest").raises(ValueError, match="Technical plot limit"):
        dna_analysis.cumulative_skew_summary(sequence, plot_step=1)
    summary = dna_analysis.cumulative_skew_summary("ACGTACGT")
    assert summary["analyzed_length"] == 8
    assert summary["visualized_length"] == 8
    assert summary["visualization"] == "FULL"
    assert math.isfinite(summary["global_gc_skew_ratio"])


def test_length_contract_rejects_a_short_analysis() -> None:
    with __import__("pytest").raises(ValueError, match="analyzed_length"):
        scale_contract.length_contract(
            input_length=150_000,
            analyzed_length=149_999,
            visualized_length=100,
            visualization="SAMPLED",
            full_scale=True,
        )
