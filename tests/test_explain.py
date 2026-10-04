"""Deterministic scientific explanations: status honesty and missing values."""

from __future__ import annotations

import math

import pytest

from modules.explain import (
    DIAGNOSIS_MARKERS,
    display_number,
    explain,
    explanation_as_dict,
)


def test_display_number_never_turns_missing_into_zero() -> None:
    assert display_number(None) == "Unavailable"
    assert display_number(float("nan")) == "Unavailable"
    assert display_number(float("inf")) == "Unavailable"
    assert display_number("") == "Unavailable"
    assert display_number(0) == "0.00"
    assert display_number(48.43, suffix="%") == "48.43%"


def test_gc_explanation_uses_actual_value() -> None:
    rec = explain("dna", {"status": "COMPUTED", "length": 100, "gc_percent": 48.43, "at_percent": 51.57})
    assert "48.43%" in rec.title or "48.43%" in rec.plain_meaning
    assert "100" in rec.plain_meaning
    assert "does not identify the organism" in rec.limitations.lower()


def test_predicted_fold_cannot_claim_experimental_structure() -> None:
    rec = explain(
        "rna_fold",
        {
            "status": "PREDICTED",
            "mfe_kcal_mol": -1.2,
            "backend": "viennarna_python",
            "tool": "ViennaRNA RNAlib",
            "method": "RNA.fold MFE",
        },
    )
    blob = (rec.plain_meaning + rec.limitations).lower()
    assert rec.status == "PREDICTED"
    assert "predicted" in blob
    assert "not an experimental" in blob or "not converted into rna 3d" in blob
    assert "experimental coordinates" not in blob


def test_illustrative_structure_is_not_experimental() -> None:
    rec = explain(
        "structure",
        {"kind": "illustrative", "structure_id": "helix", "source": "HelixScope"},
    )
    assert rec.status == "ILLUSTRATIVE"
    assert "illustrative" in rec.plain_meaning.lower()
    assert "not coordinates measured in an experiment" in rec.plain_meaning.lower()


def test_alphafold_is_predicted() -> None:
    rec = explain(
        "structure",
        {"kind": "predicted", "structure_id": "AF-P01308-F1", "source": "AlphaFold DB"},
    )
    assert rec.status == "PREDICTED"
    assert "computational prediction" in rec.plain_meaning.lower()


def test_unavailable_does_not_imply_success() -> None:
    rec = explain("rna_fold", {"status": "UNAVAILABLE", "reason": "ViennaRNA not installed"})
    blob = (rec.plain_meaning + rec.interpretation).lower()
    assert "successfully computed" not in blob
    assert "unavailable" in rec.plain_meaning.lower() or "not available" in rec.plain_meaning.lower()


def test_clinvar_is_not_diagnosis() -> None:
    rec = explain(
        "variant",
        {
            "hgvs": "17:43093557 C>G",
            "assembly": "GRCh38",
            "most_severe_consequence": "missense_variant",
            "clinvar_significance": "Pathogenic",
        },
    )
    blob = (rec.plain_meaning + rec.interpretation + rec.limitations).lower()
    for banned in DIAGNOSIS_MARKERS:
        assert banned not in blob
    assert "clinvar reports" in rec.interpretation.lower()
    assert "not a helixscope diagnosis" in rec.interpretation.lower()


def test_conservation_is_not_function() -> None:
    rec = explain("evolution", {"conservation_label": "HIGH CONSERVATION", "status": "COMPUTED"})
    assert "does not prove" in rec.interpretation.lower()
    assert "essential" in rec.limitations.lower() or "pathogenic" in rec.limitations.lower()


def test_test_reference_is_not_whole_genome() -> None:
    rec = explain(
        "crispr",
        {
            "n_guides": 3,
            "pam": "NGG",
            "reference_scope": "TEST REFERENCE",
            "test_reference": True,
        },
    )
    assert "not a genome-wide result" in rec.limitations.lower()
    assert rec.status == "PREDICTED"


def test_block_and_global_rmsd_explanations_are_distinct() -> None:
    values = {
        "engine": "RCSB PDB Structure Alignment API",
        "engine_location": "remote",
        "rmsd_global_angstrom": 0.99,
        "rmsd_block0_angstrom": 0.8,
        "n_aligned_residue_pairs": 20,
        "tm_scores": [{"type": "TM-score", "value": 0.23}],
        "method_display": "jFATCAT-rigid",
        "status": "RETRIEVED",
    }
    global_rec = explain("compare", values)
    block_rec = explain("compare_rmsd_block", values)
    assert "0.99" in global_rec.plain_meaning
    assert "0.8" in block_rec.plain_meaning
    assert "block 0" in block_rec.plain_meaning.lower() or "block rmsd" in block_rec.title.lower()
    assert "do not substitute" in block_rec.interpretation.lower()


def test_local_blast_not_merged_with_remote() -> None:
    rec = explain(
        "blast",
        {
            "backend": "NCBI BLAST+ local",
            "program": "blastn",
            "n_hits": 1,
            "not_remote_ncbi": True,
            "status": "READY",
        },
    )
    assert "local" in rec.plain_meaning.lower()
    assert "not ncbi nt/nr" in rec.plain_meaning.lower() or "not ncbi nt" in rec.limitations.lower()


def test_motif_occurrence_is_not_function() -> None:
    rec = explain("motif", {"pattern": "GAATTC", "n_hits": 2, "status": "COMPUTED"})
    assert "does not prove biological function" in rec.limitations.lower()
    assert "2" in rec.plain_meaning


def test_orf_explanation_is_predicted_not_gene() -> None:
    rec = explain("dna_orf", {"n_orfs": 3, "n_codons": 120, "status": "PREDICTED"})
    blob = (rec.plain_meaning + rec.limitations).lower()
    assert rec.status == "PREDICTED"
    assert "3" in rec.plain_meaning
    assert "not experimental gene" in blob or "not experimental" in blob


def test_contacts_do_not_prove_interaction() -> None:
    rec = explain("contacts", {"n_contacts": 12, "status": "COMPUTED"})
    assert "does not prove interaction" in rec.limitations.lower()
    assert "12" in rec.plain_meaning


def test_hostile_external_text_stays_in_explanation_values() -> None:
    rec = explain(
        "ncbi",
        {
            "accession": "NM_000000.1",
            "organism": "<script>alert(1)</script>",
            "status": "RETRIEVED",
        },
    )
    assert "<script>" in rec.plain_meaning
    from ui.components import html_escape

    escaped = html_escape(rec.plain_meaning)
    assert "<script>" not in escaped
    assert "&lt;script&gt;" in escaped


def test_unknown_kind_raises() -> None:
    with pytest.raises(ValueError):
        explain("not_a_real_kind", {})


def test_explanation_dict_is_json_safe() -> None:
    rec = explain("protein", {"length": 10, "molecular_weight_kda": 1.2, "isoelectric_point": 6.5, "gravy": -0.3})
    payload = explanation_as_dict(rec)
    assert payload["kind"] == "protein"
    assert math.isfinite(float(1.2))
    assert "software_version" in payload
