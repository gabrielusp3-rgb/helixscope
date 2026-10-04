"""Fase 10B: validacao de figura 3D, metadados mmCIF e isolamento do renderer."""

from __future__ import annotations

import ast
import copy
from pathlib import Path

import pytest

from modules import dna_3d, protein_structure, provenance, structure_scene
from ui import structure_viewer

FIXTURES = Path(__file__).parent / "fixtures"
DICKERSON = "CGCGAATTCGCG"
RNA_DUPLEX = "UUAUAUAUAUAUAA"


def test_1rna_mmcif_metadata_comes_from_file() -> None:
    parsed = protein_structure.parse_mmcif((FIXTURES / "1RNA.cif").read_text(encoding="utf-8"))
    assert parsed["entry_id"] == "1RNA"
    assert parsed["method"] == "X-RAY DIFFRACTION"
    assert parsed["resolution_angstrom"] == pytest.approx(2.25)
    chains = {item["chain_id"]: item for item in parsed["chains"]}
    assert "A" in chains
    assert RNA_DUPLEX in str(chains["A"].get("sequence") or "")


def test_1crn_mmcif_metadata_comes_from_file() -> None:
    parsed = protein_structure.parse_mmcif((FIXTURES / "1CRN.cif").read_text(encoding="utf-8"))
    assert parsed["entry_id"] == "1CRN"
    assert parsed["method"] == "X-RAY DIFFRACTION"
    assert parsed["resolution_angstrom"] == pytest.approx(1.5)


def test_figure_rejects_ready_scene_with_zero_traces() -> None:
    scene = {
        "status": "READY",
        "objects": {"backbone_traces": [], "pair_traces": [], "atom_points": None},
        "materials": {},
        "camera": {},
        "kind": "experimental",
        "structure_id": "none",
    }
    with pytest.raises(ValueError, match="Insufficient data"):
        structure_viewer.figure_from_scene(scene)


def test_figure_rejects_empty_trace_arrays() -> None:
    scene = {
        "status": "READY",
        "objects": {
            "backbone_traces": [{"x": [], "y": [], "z": [], "name": "empty"}],
            "pair_traces": [],
            "atom_points": None,
        },
        "materials": {"line_width": 4, "marker_size": 5, "opacity": 0.9},
        "camera": {},
        "kind": "experimental",
        "structure_id": "none",
    }
    with pytest.raises(ValueError, match="Insufficient data"):
        structure_viewer.figure_from_scene(scene)


def test_structure_viewer_has_no_network_imports() -> None:
    source = Path(structure_viewer.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module.split(".")[0])
    banned = {"urllib", "requests", "httpx", "aiohttp", "http", "socket"}
    assert not banned.intersection(imported)


def test_structure_scene_has_no_cas_assembler() -> None:
    source = Path(structure_scene.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    assert "place_cas_guide_dna" not in names
    assert "synthetic_cas_complex" not in names


def test_backbone_trace_name_passes_through_mapping_role() -> None:
    auto = dna_3d.representation_for_sequence(DICKERSON)
    envelope = copy.deepcopy(auto["envelope"])
    for item in list(envelope.get("mappings") or []):
        if str(item.get("chain_id") or "") == "A":
            item["role"] = "target_dna"
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=str(auto["sequence_hash"]),
        representation="backbone",
    )
    assert scene["status"] == "READY"
    names = " ".join(str(trace.get("name") or "") for trace in scene["objects"]["backbone_traces"])
    assert "target dna" in names.lower()
    hover = " ".join(
        " ".join(str(text) for text in (trace.get("hovertext") or []))
        for trace in scene["objects"]["backbone_traces"]
    )
    assert "role target_dna" in hover


def test_lod_changes_materials_not_coordinates() -> None:
    auto = dna_3d.representation_for_sequence(DICKERSON, prefer="illustrative")
    envelope = auto["envelope"]
    high = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=str(auto["sequence_hash"]),
        lod="high",
    )
    low = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=str(auto["sequence_hash"]),
        lod="low",
    )
    hx = list(high["objects"]["backbone_traces"][0]["x"])
    lx = list(low["objects"]["backbone_traces"][0]["x"])
    assert lx[0] == pytest.approx(hx[0])
    assert lx[-1] == pytest.approx(hx[-1])
    assert len(lx) < len(hx)
    assert int(high["materials"]["marker_size"]) > int(low["materials"]["marker_size"])
    assert int(high["materials"]["line_width"]) > int(low["materials"]["line_width"])


def test_illustrative_dna_roundtrip_matches_sequence() -> None:
    seq = "ATGCGAATTCGG"
    rep = dna_3d.representation_for_sequence(seq, prefer="illustrative")
    residues = list((rep.get("model") or {}).get("residues") or [])
    drawn = "".join(str(item.get("base") or "") for item in residues)
    start = int(rep.get("view_start") or 0)
    end = int(rep.get("view_end") or len(seq))
    assert drawn == seq[start:end]


def test_version_is_10b() -> None:
    assert provenance.HELIXSCOPE_VERSION.startswith("0.")
