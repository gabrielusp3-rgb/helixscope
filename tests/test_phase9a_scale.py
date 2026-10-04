"""Fase 9A: escala 50k, LOD visual, janela de viewer, orcamento de grelha.

Nenhum teste inventa coordenadas experimentais, arvores ou scores. A sequencia
de 50 kb e sintetica apenas como input de comprimento; as metricas sao calculadas
sobre ela. parser fixtures nao sao deposicoes PDB.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from modules import (
    dna_analysis,
    nucleic_geometry,
    protein_structure,
    provenance,
    region_nav,
    scale_profile,
    structure_scene,
)
from ui import charts, components, structure_viewer

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _fail_network(*_args, **_kwargs):
    raise AssertionError("network")


def _seq_50k() -> str:
    return ("ACGT" * 12_500)


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


def test_50k_dna_analysis_is_valid_and_complete():
    seq = _seq_50k()
    assert len(seq) == scale_profile.TARGET_DNA_NT
    info = dna_analysis.validate_for_molecule(seq, "DNA")
    assert info["is_valid"] is True
    assert info["length"] == 50_000
    started = time.perf_counter()
    gc = dna_analysis.gc_content(seq)
    composition = dna_analysis.nucleotide_composition(seq)
    skew = dna_analysis.gc_skew(seq, window=1000)
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    assert gc == pytest.approx(50.0)
    total = sum(int(item["count"]) for item in composition.values())
    assert total == 50_000
    assert len(skew) == 50
    assert elapsed_ms < 15_000
    assert "instant" not in str(elapsed_ms)


def test_50k_sliding_window_refuses_silent_step_one():
    seq = _seq_50k()
    with pytest.raises(ValueError, match="Technical display/memory limit"):
        dna_analysis.gc_sliding_window(seq, window=100, step=1)
    with pytest.raises(ValueError, match="Technical display/memory limit"):
        dna_analysis.windowed_profiles(seq, window=100, step=1)
    plan = scale_profile.analysis_plan(50_000)
    profile = dna_analysis.gc_sliding_window(
        seq,
        window=int(plan["suggested_window"]),
        step=int(plan["suggested_step"]),
    )
    assert 0 < len(profile) <= scale_profile.MAX_WINDOW_PROFILE_ROWS


def test_contiguous_view_on_50k_caps_helix_without_skipping_bases():
    view = region_nav.contiguous_view(50_000, start=0, end=50_000, lod="high")
    assert view["n"] == region_nav.MAX_HELIX_NT_HIGH
    assert view["partial_view"] is True
    assert view["end"] - view["start"] == view["n"]
    seq = _seq_50k()
    model = nucleic_geometry.illustrative_bdna_model(seq)
    assert model["kind"] == nucleic_geometry.KIND_ILLUSTRATIVE
    assert model["partial_view"] is True
    assert model["n_residues"] == region_nav.MAX_HELIX_NT_HIGH
    indices = [int(item["index_0based"]) for item in model["residues"]]
    assert indices == list(range(indices[0], indices[0] + len(indices)))
    assert model["sequence_hash"] == provenance.sequence_digest(seq)
    assert "experimental" not in str(model["kind"]).lower()


def test_illustrative_region_keeps_global_indices():
    seq = "ATGC" * 50
    model = nucleic_geometry.illustrative_bdna_model(
        seq, region_start=40, region_end=80, lod="high"
    )
    assert model["view_start"] == 40
    assert model["residues"][0]["index_0based"] == 40
    assert model["residues"][0]["base"] == seq[40]


def test_lod_stride_omits_points_without_changing_remaining_xyz():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    snapshot = structure_scene.scientific_snapshot(envelope)
    high = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        lod="high",
    )
    low = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        lod="low",
    )
    assert high["status"] == "READY"
    assert low["status"] == "READY"
    assert high["partial_view"] is False
    assert low["partial_view"] is True
    assert high["visual_stride"] == 1
    assert low["visual_stride"] == 4
    hx = list(high["objects"]["backbone_traces"][0]["x"])
    hy = list(high["objects"]["backbone_traces"][0]["y"])
    hz = list(high["objects"]["backbone_traces"][0]["z"])
    lx = list(low["objects"]["backbone_traces"][0]["x"])
    ly = list(low["objects"]["backbone_traces"][0]["y"])
    lz = list(low["objects"]["backbone_traces"][0]["z"])
    assert lx[0] == pytest.approx(hx[0])
    assert lx[-1] == pytest.approx(hx[-1])
    assert len(lx) < len(hx)
    high_set = {
        (round(float(x), 6), round(float(y), 6), round(float(z), 6))
        for x, y, z in zip(hx, hy, hz)
    }
    for x, y, z in zip(lx, ly, lz):
        assert (round(float(x), 6), round(float(y), 6), round(float(z), 6)) in high_set
    assert structure_scene.scientific_snapshot(envelope) == snapshot
    assert high["source_snapshot"] == snapshot
    assert low["source_snapshot"] == snapshot


def test_1crn_high_lod_still_draws_all_mapped_ca():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        lod="high",
    )
    assert scene["status"] == "READY"
    assert len(scene["objects"]["backbone_traces"][0]["x"]) == 46
    figure = structure_viewer.figure_from_scene(scene)
    assert figure.layout.legend.y < 0.5


def test_sequence_display_region_window():
    seq = "ACGT" * 40
    html = components.sequence_display(seq, "DNA", view_start=8, view_end=16)
    assert "Showing residues 8" in html
    assert "16 of 160" in html


def test_downsample_keeps_first_and_last():
    kept = region_nav.downsample_indices(10, 4)
    assert kept[0] == 0
    assert kept[-1] == 9
    assert kept == sorted(set(kept))


def test_phylogenetic_figure_hides_labels_when_many_leaves():
    nodes = [{"id": "root", "x": 0.0, "y": 10.0, "is_leaf": False}]
    edges = []
    for index in range(scale_profile.MAX_TREE_LEAF_LABELS + 4):
        leaf_id = f"leaf{index}"
        nodes.append(
            {
                "id": leaf_id,
                "x": 1.0,
                "y": float(index),
                "is_leaf": True,
                "label": leaf_id,
            }
        )
        edges.append(
            {
                "parent": "root",
                "child": leaf_id,
                "x0": 0.0,
                "y0": 10.0,
                "x1": 1.0,
                "y1": float(index),
                "length": 1.0,
            }
        )
    figure = charts.phylogenetic_tree_figure(
        {"nodes": nodes, "edges": edges, "scale_bar": {"present": False}}
    )
    leaf_traces = [trace for trace in figure.data if trace.name == "Leaves"]
    assert leaf_traces
    assert leaf_traces[0].mode == "markers"
    assert "hover only" in str(figure.layout.title.text)


def test_inspector_and_meta_wrap_long_text():
    html = components.inspector_panel(
        {
            "IDENTITY": [("Hash", "a" * 80)],
            "SOURCE": [("Note", "experimental")],
        }
    )
    assert "helix-inspector" in html
    assert "a" * 80 in html
    grid = components.meta_grid([("Structure ID", "1CRN experimental mapping status EXACT")])
    assert "helix-meta-grid" in grid
    row = components.badge_row(["EXPERIMENTAL", "PARTIAL_VIEW", "ILLUSTRATIVE"])
    assert "helix-badge-row" in row
    assert "EXPERIMENTAL" in row
