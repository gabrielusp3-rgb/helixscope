"""Testes do scene model 3D: coordenadas reais, mapping, stale, imutabilidade.

1CRN.cif e rcsb_entry_1CRN.json sao copias publicas do RCSB PDB (crambin,
X-ray, 1.5 A). parser_* nao sao deposicoes. Nenhum teste inventa pLDDT,
cobertura ou um PDB/AlphaFold vivo.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from modules import protein_structure, provenance, scientific_checks, structure_scene
from ui import structure_viewer

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _fail_network(*_args, **_kwargs):
    raise AssertionError("network")


def _1crn_result() -> dict:
    meta = protein_structure.parse_rcsb_entry_metadata(json.loads(_load("rcsb_entry_1CRN.json")))
    return protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
    )


def _parser_result(query: str, filename: str, structure_id: str) -> dict:
    return protein_structure.load_protein_structure(
        query,
        source="parser fixture",
        structure_id=structure_id,
        structure_text=_load(filename),
        metadata={
            "kind": "experimental",
            "method": "NMR" if "nmr" in filename else "X-ray",
            "structure_id": structure_id,
        },
        urlopen_fn=_fail_network,
    )


def test_1crn_scene_uses_deposited_ca_and_is_deterministic():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    snapshot = structure_scene.scientific_snapshot(envelope)
    started = time.perf_counter()
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        color_mode="chain",
    )
    elapsed = time.perf_counter() - started
    assert elapsed < 2.0
    assert scene["status"] == "READY"
    assert scene["engine"] == "plotly"
    assert scene["kind"] == "experimental"
    assert scene["kind_label"] == "Experimental Structure"
    assert scene["mapping_status"] == "EXACT"
    assert scene["source"] == "RCSB PDB"
    traces = scene["objects"]["backbone_traces"]
    assert traces
    first_ca = result["residue_mapping"][0]["ca"]
    assert traces[0]["x"][0] == pytest.approx(first_ca["x"])
    assert traces[0]["y"][0] == pytest.approx(first_ca["y"])
    assert traces[0]["z"][0] == pytest.approx(first_ca["z"])
    assert len(traces[0]["x"]) == 46
    again = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        color_mode="chain",
    )
    assert structure_scene.scene_fingerprint(scene) == structure_scene.scene_fingerprint(again)
    assert structure_scene.scientific_snapshot(envelope) == snapshot
    traces[0]["x"][0] = 999.0
    assert envelope["atoms"][0]["x"] != 999.0
    assert result["atoms"][0]["x"] != 999.0
    assert result["sequence_hash"] == envelope["sequence_hash"]
    assert result["content_hash"] == envelope["structure_hash"]
    dumped = str(scene).lower()
    assert "http://" not in dumped
    assert "https://" not in dumped
    report = structure_scene.structure_3d_report(scene, envelope)
    assert report["structure_id"] == "1CRN"
    assert "quality" not in str(report).lower()
    assert "realism" not in str(report).lower()


def test_stale_sequence_hash_does_not_render():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    other = provenance.sequence_digest("ACDEFGHIKLMNPQRSTVWY")
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=other,
        representation="backbone",
    )
    assert scene["status"] == "STALE"
    assert scene["n_atoms_rendered"] == 0
    assert scene["objects"]["backbone_traces"] == []
    assert "stale" in str(scene.get("message") or "").lower()


def test_auth_offset_mapping_and_reverse_lookup():
    result = _parser_result("AK", "parser_auth_offset.pdb", "parser_auth_offset")
    assert result["source"] == "parser fixture"
    assert result["mapping_status"] == "PARTIAL"
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="ca",
    )
    assert scene["status"] == "READY"
    assert scene["mapping_status"] == "PARTIAL"
    row = structure_scene.query_to_structure_residue(envelope, 0)
    assert row["auth_seq_id"] == 10
    assert row["query_index_0based"] == 0
    reverse = structure_scene.structure_residue_to_query(
        envelope, chain_id="A", auth_seq_id=10, label_seq_id=10
    )
    assert reverse["query_index_0based"] == 0
    assert structure_scene.structure_residue_to_query(envelope, chain_id="A", auth_seq_id=1) is None


def test_nmr_models_missing_residue_and_chain_groups():
    result = _parser_result("AKA", "parser_nmr_altloc.cif", "parser_nmr_altloc")
    assert result["models"] == [1, 2]
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        model_filter=1,
        selected_query_index=1,
    )
    assert scene["status"] == "READY"
    assert 1 in scene["residues_without_coordinates"]
    marker = scene["objects"]["selection_marker"]
    assert marker is None
    mapped = structure_scene.query_to_structure_residue(envelope, 1)
    assert mapped["has_coordinates"] is False
    model_two = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="ca",
        model_filter=2,
    )
    assert model_two["model_number"] == 2
    assert model_two["partial_view"] is True
    assert model_two["view_label"] == "Partial structure view"


def test_nan_coordinates_are_parsing_error_not_drawn():
    envelope = protein_structure.structure_3d_input(_1crn_result())
    envelope["atoms"][0]["x"] = float("nan")
    scene = structure_scene.build_scene(envelope, representation="atoms")
    assert scene["status"] == "PARSING_ERROR"
    assert scene["n_atoms_rendered"] == 0


def test_resource_limit_does_not_render_partial_as_complete():
    stub_atom = {
        "group": "ATOM",
        "atom_name": "CA",
        "auth_asym_id": "A",
        "label_asym_id": "A",
        "x": 1.0,
        "y": 2.0,
        "z": 3.0,
        "model": 1,
        "comp_id": "ALA",
    }
    envelope = {
        "source": "test stub",
        "kind": "experimental",
        "structure_id": "limit-stub",
        "sequence_hash": "stub",
        "atoms": [dict(stub_atom) for _ in range(protein_structure.MAX_ATOMS + 1)],
        "residues": [],
        "mappings": [],
        "chain_id": "A",
        "model_number": 1,
        "renderer": "data-only",
    }
    scene = structure_scene.build_scene(envelope, expected_sequence_hash="stub", representation="atoms")
    assert scene["status"] == "RESOURCE_LIMIT"
    assert scene["n_atoms_rendered"] == 0
    assert scene["objects"]["backbone_traces"] == []


def test_empty_atoms_are_unavailable_not_fake_geometry():
    envelope = {
        "source": "test stub",
        "kind": "experimental",
        "structure_id": "empty-stub",
        "sequence_hash": "stub",
        "atoms": [],
        "residues": [],
        "mappings": [],
        "chain_id": "A",
        "model_number": 1,
        "renderer": "data-only",
    }
    scene = structure_scene.build_scene(envelope, expected_sequence_hash="stub")
    assert scene["status"] == "UNAVAILABLE"
    assert scene["objects"]["backbone_traces"] == []


def test_plddt_color_refused_for_experimental_and_url_rejected():
    envelope = protein_structure.structure_3d_input(_1crn_result())
    scene = structure_scene.build_scene(envelope, color_mode="plddt")
    assert scene["status"] == "READY"
    assert scene["color_mode"] == "chain"
    assert "pLDDT" in scene["color_note"]
    envelope["download_url"] = "https://files.rcsb.org/download/1CRN.cif"
    blocked = structure_scene.build_scene(envelope)
    assert blocked["status"] == "ERROR"


def test_predicted_kind_uses_deposited_b_factor_bins_on_parser_fixture():
    result = _parser_result("AKA", "parser_nmr_altloc.cif", "parser_nmr_altloc")
    envelope = protein_structure.structure_3d_input(result)
    envelope["kind"] = "predicted"
    envelope["source"] = "parser fixture"
    envelope["kind_label"] = "Predicted Structure"
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        color_mode="plddt",
        representation="ca",
    )
    assert scene["status"] == "READY"
    assert scene["color_mode"] == "plddt"
    assert "not a probability" in scene["color_note"].lower()
    assert "accuracy percentage" in scene["color_note"].lower()
    assert scene["source"] == "parser fixture"


def test_selection_and_span_highlight_ids():
    envelope = protein_structure.structure_3d_input(_1crn_result())
    indices = structure_scene.highlight_query_indices_for_span(0, 3)
    assert indices == [0, 1, 2]
    scene = structure_scene.build_scene(
        envelope,
        representation="ca",
        color_mode="selection",
        selected_query_index=0,
        highlight_query_indices=indices,
    )
    assert scene["selected_query_index"] == 0
    marker = scene["objects"]["selection_marker"]
    assert marker is not None
    parsed = structure_scene.parse_selection_payload(marker["customdata"][0])
    assert parsed["query_index_0based"] == 0
    assert parsed["chain_id"] == "A"
    event = structure_scene.make_structure_event(
        "residueSelected",
        sequence_hash=envelope["sequence_hash"],
        structure_id="1CRN",
        chain_id="A",
        model=1,
        query_index=0,
        label_seq_id=1,
        auth_seq_id=1,
    )
    assert event["kind"] == "residueSelected"
    assert event["query_index_0based"] == 0


def test_figure_from_scene_matches_backbone_coordinates():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(envelope, representation="backbone")
    figure = structure_viewer.figure_from_scene(scene)
    assert list(figure.data[0].x) == scene["objects"]["backbone_traces"][0]["x"]
    assert "<script" not in str(figure.data[0].text)
    empty = structure_scene.build_scene(
        {"source": "test stub", "atoms": [], "sequence_hash": "x", "kind": "unavailable"}
    )
    with pytest.raises(ValueError, match="Insufficient data"):
        structure_viewer.figure_from_scene(empty)


def test_webgl_note_does_not_claim_gpu_detection():
    note = structure_scene.webgl_capability_note()
    assert "user-agent" in note.lower()
    assert "unavailable" in note.lower()
    log = structure_scene.renderer_log("renderer_init", status="READY", n_atoms=327)
    assert "token" not in str(log).lower()
    assert log["engine"] == "plotly"


def test_atom_and_ca_scene_counts_for_1crn():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    atoms = structure_scene.build_scene(envelope, representation="atoms")
    assert atoms["n_atoms_source"] == 327
    assert atoms["n_atoms_rendered"] == 327
    ca = structure_scene.build_scene(envelope, representation="ca")
    assert ca["n_atoms_rendered"] == 46
    assert ca["view_scope"] == "selected_chain"
    assert ca["view_label"] != "Full structure"


def test_hydropathy_overlay_uses_provided_values_not_invented_scores():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    overlay = {0: -4.5, 1: 4.5}
    scene = structure_scene.build_scene(
        envelope,
        color_mode="hydropathy",
        overlay=overlay,
        representation="ca",
    )
    assert scene["color_mode"] == "hydropathy"
    assert "Kyte-Doolittle" in scene["color_note"]
    missing = structure_scene.build_scene(envelope, color_mode="conservation")
    assert missing["color_mode"] == "chain"


def test_charge_overlay_uses_protein_engine_values():
    from modules import protein_analysis

    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    charges = protein_analysis.per_residue_formal_charge(CRAMBIN)
    overlay = {index: float(value) for index, value in enumerate(charges)}
    scene = structure_scene.build_scene(
        envelope,
        color_mode="charge",
        overlay=overlay,
        representation="ca",
    )
    assert scene["color_mode"] == "charge"
    assert "net_charge" in scene["color_note"]
    assert "Henderson-Hasselbalch" in scene["color_note"]
    assert scene["timings_ms"]["scene_construction"] >= 0
    missing = structure_scene.build_scene(envelope, color_mode="charge")
    assert missing["color_mode"] == "chain"


def test_two_chains_visible_subset_is_partial_not_full_structure():
    result = _parser_result("AK", "parser_two_chains.cif", "parser_two_chains")
    envelope = protein_structure.structure_3d_input(result)
    snapshot = structure_scene.scientific_snapshot(envelope)
    hidden = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="atoms",
        chain_filter="A",
        visible_chains=["A"],
    )
    assert hidden["status"] == "READY"
    assert hidden["partial_view"] is True
    assert hidden["view_label"] == "Partial structure view"
    assert "Full structure" not in hidden["view_label"]
    both = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="atoms",
        chain_filter="A",
        visible_chains=["A", "B"],
    )
    assert both["n_atoms_rendered"] >= hidden["n_atoms_rendered"]
    assert both["visible_chains"] == ["A", "B"]
    mutated = dict(hidden["objects"]["atom_points"])
    mutated["x"] = [999.0] * len(mutated["x"])
    assert structure_scene.scientific_snapshot(envelope) == snapshot
    report = structure_scene.structure_3d_report(hidden, envelope)
    assert report["partial_view"] is True
    assert report["kind"] == "experimental"


def test_mutating_coordinates_or_kind_is_not_silently_accepted():
    original = _1crn_result()
    envelope = protein_structure.structure_3d_input(original)
    envelope["atoms"][0]["x"] = float("inf")
    scene = structure_scene.build_scene(envelope, representation="atoms")
    assert scene["status"] == "PARSING_ERROR"
    assert original["kind"] == "experimental"
    experimental = protein_structure.structure_3d_input(original)
    plddt = structure_scene.build_scene(experimental, color_mode="plddt")
    assert plddt["color_mode"] == "chain"
    assert original["atoms"][0]["x"] != float("inf")


def test_rna_and_crispr_modules_do_not_emit_fake_3d():
    import inspect

    from modules import crispr, rna_folding

    rna_src = inspect.getsource(rna_folding).lower()
    assert "structure_scene" not in rna_src
    assert "scatter3d" not in rna_src
    assert "fake 3d" not in rna_src
    crispr_src = inspect.getsource(crispr).lower()
    assert "structure_scene" not in crispr_src
    assert "scatter3d" not in crispr_src
    assert "cas9 complex" not in crispr_src or "unavailable" in crispr_src

