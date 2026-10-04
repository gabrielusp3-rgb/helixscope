"""Restoration payloads: Core must still emit the Streamlit-era scientific arrays."""

from __future__ import annotations

import math

import pytest

from helixscope_core.alignment import pairwise_align
from helixscope_core.dna import analyze_dna
from helixscope_core.protein import analyze_protein
from helixscope_core.rna import analyze_rna
from helixscope_core.structure import load_catalog_complex
from helixscope_core.errors import StructureError


def test_analyze_dna_includes_dinucleotides_and_restriction_rows() -> None:
    result = analyze_dna(
        "ATGC" * 40,
        include_orfs=True,
        include_restriction=True,
    )
    assert result["status"] == "COMPUTED"
    assert "CG" in result["dinucleotides"]
    oe = result["dinucleotides"]["CG"]["observed_expected"]
    assert isinstance(oe, float)
    assert result["composition_consistent"] is True
    assert isinstance(result["restriction_enzyme_rows"], list)
    assert isinstance(result["restriction_hits"], list)
    assert isinstance(result["sequence_features"], list)


def test_analyze_dna_dinucleotide_nan_is_not_zero() -> None:
    result = analyze_dna("A" * 12)
    oe = result["dinucleotides"]["CG"]["observed_expected"]
    assert math.isnan(oe)


def test_analyze_protein_categories_and_charge_profile() -> None:
    result = analyze_protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
    assert result["amino_acid_categories_status"] == "COMPUTED"
    assert "charged_total" in result["amino_acid_categories"]
    assert result["charge_profile_status"] == "COMPUTED"
    assert result["charge_profile"][0]["method"]


def test_analyze_rna_includes_dinucleotides() -> None:
    result = analyze_rna("AUGGCUGCUGCUGCUUAA")
    assert result["status"] == "COMPUTED"
    assert "dinucleotides" in result
    assert "AT" in result["dinucleotides"]


def test_pairwise_align_exposes_column_classes() -> None:
    result = pairwise_align("ACGTACGTACGT", "ACGTACGGACGT", mode="global")
    classes = result["column_classes"]
    assert len(classes) == len(result["aligned_seq1"])
    assert set(classes) <= {"match", "mismatch", "gap"}


def test_catalog_complex_rejects_unknown_pdb() -> None:
    with pytest.raises(StructureError):
        load_catalog_complex("1CRN")


def test_dna_3d_viewer_reports_lod_and_inspect_mapping() -> None:
    from helixscope_core.dna import dna_3d_viewer

    payload = dna_3d_viewer(
        "ATGC" * 20,
        prefer="illustrative",
        lod="low",
        region_start=0,
        region_end=40,
        inspect_position=2,
    )
    assert payload["visual_lod"] == "low"
    assert payload["inspect_base"] == "G"
    assert payload["helix_limit_nt"] == 80
    assert "position_mapping" in payload


def test_analyze_rna_attaches_circular_layout_when_fold_has_sequence() -> None:
    from modules.rna_folding import detect_python_bindings

    if not detect_python_bindings().get("available"):
        pytest.skip("ViennaRNA Python module RNA is not installed")
    result = analyze_rna("GGGGAAAACCCC", fold=True)
    fold = result.get("fold")
    if not isinstance(fold, dict):
        return
    if fold.get("status") in {"PREDICTED", "COMPUTED"} and fold.get("sequence"):
        layout = fold.get("circular_layout")
        assert isinstance(layout, dict)
        assert len(layout.get("x") or []) == len(str(fold.get("sequence") or ""))


def test_design_guides_keeps_local_off_target_hits() -> None:
    from helixscope_core.crispr import design_guides

    guides = design_guides(
        "ACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGT",
        "SpCas9",
        run_off_target=True,
    )
    assert guides
    assert isinstance(guides[0].get("off_targets"), list)
    assert "n_off_targets" in guides[0]
