"""Testes da fundacao de proveniencia, evidencia e invariantes cientificas.

Nenhum teste faz chamadas de rede. Nenhum teste inventa um resultado analitico:
eles conferem classificacao, redacao de segredos e consistencia de coordenadas.
"""

import pytest

from modules import crispr, provenance, scientific_checks


def test_build_provenance_computed_and_unknown_status():
    record = provenance.build_provenance(
        status="computed",
        algorithm="gc_content",
        parameters={"window": 100},
        input_identifier="user-paste",
    )
    assert record["status"] == "COMPUTED"
    assert record["software"] == "HelixScope"
    assert record["software_version"] == provenance.HELIXSCOPE_VERSION
    assert record["algorithm"] == "gc_content"
    assert record["parameters"]["window"] == 100
    assert record["computed_timestamp"]
    assert record["retrieval_timestamp"] == ""
    unknown = provenance.build_provenance(status="made-up", algorithm="none")
    assert unknown["status"] == "ERROR"


def test_build_provenance_retrieved_sets_retrieval_timestamp():
    record = provenance.build_provenance(
        status="RETRIEVED",
        algorithm="Entrez efetch",
        source="NCBI Entrez",
        database="nucleotide",
        accession="NM_000518.5",
    )
    assert record["status"] == "RETRIEVED"
    assert record["retrieval_timestamp"]
    assert record["computed_timestamp"] == ""
    assert record["accession"] == "NM_000518.5"


def test_structure_disclaimer_and_invalid_kind():
    text = provenance.structure_disclaimer("experimental")
    assert "Experimental structure" in text
    predicted = provenance.structure_disclaimer("predicted")
    assert "not an experimental measurement" in predicted
    unavailable = provenance.structure_disclaimer("unavailable")
    assert "unavailable" in unavailable.lower()
    with pytest.raises(ValueError):
        provenance.structure_disclaimer("crystal-like")


def test_redact_secrets_omits_env_and_key_names(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "secret-ncbi-key-value")
    payload = {
        "ncbi_api_key": "should-hide",
        "message": "secret-ncbi-key-value",
        "module": "ncbi_fetch",
    }
    redacted = provenance.redact_secrets(payload)
    assert redacted["ncbi_api_key"] == "[redacted]"
    assert redacted["message"] == "[redacted]"
    assert redacted["module"] == "ncbi_fetch"
    event = provenance.diagnostic_event(
        module="ncbi_fetch",
        operation="fetch",
        status="error",
        exception="NCBIQueryError",
    )
    assert event["module"] == "ncbi_fetch"
    assert "NCBI_API_KEY" not in str(event)


def test_percent_and_composition_invariants():
    assert scientific_checks.percent_is_valid(50.0) is True
    assert scientific_checks.percent_is_valid(float("nan")) is True
    assert scientific_checks.percent_is_valid(-1.0) is False
    assert scientific_checks.percent_is_valid(101.0) is False
    assert scientific_checks.composition_counts_match_length({"A": 2, "T": 2}, 4) is True
    assert scientific_checks.composition_counts_match_length({"A": 2}, 4) is False
    assert scientific_checks.composition_counts_match_length({"A": 1}, -1) is False


def test_span_alignment_chart_and_motif_invariants():
    assert scientific_checks.span_is_valid(0, 4, 4, base=0) is True
    assert scientific_checks.span_is_valid(1, 4, 4, base=1) is True
    assert scientific_checks.span_is_valid(0, 5, 4, base=0) is False
    with pytest.raises(ValueError):
        scientific_checks.span_is_valid(0, 1, 4, base=2)
    assert scientific_checks.alignment_rows_are_consistent("AC-GT", "ACAGT") is True
    assert scientific_checks.alignment_rows_are_consistent("ACGT", "AC") is False
    assert scientific_checks.alignment_rows_are_consistent("", "AC") is False
    assert scientific_checks.chart_xy_lengths_match([1, 2], [3, 4]) is True
    assert scientific_checks.chart_xy_lengths_match([1], [3, 4]) is False
    assert scientific_checks.motif_hit_matches_sequence("GAATTC", 0, 6, "GAATTC") is True
    assert scientific_checks.motif_hit_matches_sequence("GAATTC", 0, 6, "TTTTTT") is False


def test_crispr_guide_in_sequence_plus_and_minus():
    sequence = crispr.EXAMPLE_EMX1_SPCAS9
    plus = {
        "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
        "position": 1,
        "strand": "+",
    }
    assert scientific_checks.crispr_guide_in_sequence(sequence, plus) is True
    found = crispr.find_guides(sequence, "SpCas9")
    assert found
    assert all(scientific_checks.crispr_guide_in_sequence(sequence, guide) for guide in found)
    outside = {
        "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
        "position": 50,
        "strand": "+",
    }
    assert scientific_checks.crispr_guide_in_sequence(sequence, outside) is False
    assert scientific_checks.crispr_guide_in_sequence(sequence, {}) is False


def test_sequence_digest_envelope_and_csv_cells():
    digest = provenance.sequence_digest("ACGT")
    assert digest == provenance.sequence_digest("ACGT")
    assert digest != provenance.sequence_digest("ACGU")
    envelope = provenance.analysis_envelope(
        module="DNA",
        payload={"gc_percent": float("nan"), "count": 0},
        status="COMPUTED",
        algorithm="gc_content",
        sequence="NNNN",
    )
    assert envelope["gc_percent"] is None
    assert envelope["count"] == 0
    assert envelope["input_hash"] == provenance.sequence_digest("NNNN")
    assert envelope["status"] == "COMPUTED"
    assert provenance.csv_cell(0) == 0
    assert provenance.csv_cell(0.0) == 0.0
    assert provenance.csv_cell(float("nan")) == "N/A"
    assert provenance.csv_cell(None) == "N/A"
    assert provenance.csv_cell(float("inf")) == "ERROR"
    snap = provenance.workspace_snapshot(
        sequence="ACGT",
        molecule="DNA",
        source="user sequence",
        identifier="seq1",
    )
    assert snap["length"] == 4
    assert snap["input_hash"] == provenance.sequence_digest("ACGT")
    assert "sequence" not in snap


def test_make_span_and_protein_composition_check():
    span = scientific_checks.make_span(
        start=1,
        end=4,
        sequence="ACGTAA",
        strand="+",
        label="Predicted ORF",
        kind="orf",
        status="PREDICTED",
    )
    assert span["sequence"] == "CGT"
    assert span["status"] == "PREDICTED"
    assert span["type"] == "orf"
    assert span["kind"] == "orf"
    with pytest.raises(ValueError):
        scientific_checks.make_span(start=0, end=9, sequence="ACGT")
    assert scientific_checks.protein_composition_matches_length({"A": 2, "K": 2}, 4)
    assert scientific_checks.protein_composition_matches_length({"A": 1}, 4) is False
    valid_hsp = {
        "identities": 20,
        "alignment_length": 20,
        "identity_pct": 100.0,
        "query_coverage_pct": 100.0,
        "query_from": 1,
        "query_to": 20,
        "hit_from": 10,
        "hit_to": 29,
        "qseq": "A" * 20,
        "hseq": "A" * 20,
        "midline": "|" * 20,
        "positives": 20,
        "gaps": 0,
    }
    assert scientific_checks.blast_hsp_is_valid(valid_hsp, query_len=20, hit_len=100)
    invalid = dict(valid_hsp)
    invalid["identities"] = 21
    assert scientific_checks.blast_hsp_is_valid(invalid, query_len=20, hit_len=100) is False
    assert scientific_checks.rna_structure_length_matches("GGGAAACCC", "(((...)))")
    assert scientific_checks.rna_structure_length_matches("GGG", "..") is False
    assert scientific_checks.rna_parentheses_are_balanced("(((...)))")
    assert scientific_checks.rna_parentheses_are_balanced("((...)") is False
    assert scientific_checks.mfe_kcal_mol_is_valid(0.0)
    assert scientific_checks.mfe_kcal_mol_is_valid(-3.3)
    assert scientific_checks.mfe_kcal_mol_is_valid(float("nan")) is False
    assert scientific_checks.mfe_kcal_mol_is_valid(float("inf")) is False
    assert scientific_checks.atom_coordinate_is_finite(1.0, 2.0, 3.0)
    assert scientific_checks.atom_coordinate_is_finite(0.0, 0.0, 0.0)
    assert scientific_checks.atom_coordinate_is_finite(float("nan"), 0.0, 0.0) is False
    assert scientific_checks.protein_mapping_covers_query(46, 46)
    assert scientific_checks.protein_mapping_covers_query(10, 46) is False
    assert scientific_checks.mapping_status_is_declared("EXACT")
    assert scientific_checks.mapping_status_is_declared("BEST_EFFORT")
    assert scientific_checks.mapping_status_is_declared("exact-enough") is False
    assert provenance.hashes_match("abc", "abc")
    assert provenance.hashes_match("", "") is False
    assert provenance.hashes_match("abc", "def") is False
    assert scientific_checks.to_1based_inclusive(0, 1) == (1, 1)
    assert scientific_checks.to_1based_inclusive(0, 4) == (1, 4)
    last = scientific_checks.to_1based_inclusive(3, 4)
    assert last == (4, 4)
    assert scientific_checks.to_0based_half_open(1, 1) == (0, 1)
    assert scientific_checks.to_0based_half_open(1, 4) == (0, 4)
    assert scientific_checks.to_0based_half_open(*scientific_checks.to_1based_inclusive(2, 3)) == (2, 3)
    with pytest.raises(ValueError):
        scientific_checks.to_1based_inclusive(0, 0)
    with pytest.raises(ValueError):
        scientific_checks.to_0based_half_open(0, 1)
    assert scientific_checks.bootstrap_support_is_valid(None, replicates=0)
    assert scientific_checks.bootstrap_support_is_valid(0.0, replicates=0) is False
    assert scientific_checks.bootstrap_support_is_valid(0.0, replicates=10)
    assert scientific_checks.branch_length_is_valid(None)
    assert scientific_checks.branch_length_is_valid(float("nan")) is False
    assert scientific_checks.phylogenetic_leaf_set_matches_msa(["a", "b"], ["b", "a"])
    assert scientific_checks.phylogenetic_leaf_set_matches_msa(["a", "a"], ["a", "b"]) is False

