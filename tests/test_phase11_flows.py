"""Fase 11: selecao 3D, erros classificados, 1CRN bundled, stale, LOD, MSA."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from streamlit.testing.v1 import AppTest

from modules import dna_3d, msa, molecule_selection, protein_structure, provenance, region_nav, structure_scene
from ui import components, structure_viewer
from tests.helix_apptest import open_module, switch_module

ROOT = Path(__file__).resolve().parents[1]
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
DICKERSON = "CGCGAATTCGCG"
SEQ_50K = "ATGC" * 12_500


def _ui_text(app: AppTest) -> str:
    return " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption", "info", "error", "warning")
        for item in list(getattr(app, bucket, []) or [])
    )


def test_plotly_selection_payload_roundtrip() -> None:
    event = SimpleNamespace(
        selection=SimpleNamespace(
            points=[{"customdata": [4, "A", 1, 5, 5, "", "P"]}]
        )
    )
    parsed = structure_viewer.parse_plotly_selection(event)
    assert parsed is not None
    assert parsed["query_index_0based"] == 4
    assert parsed["chain_id"] == "A"
    assert parsed["atom_name"] == "P"
    empty = structure_viewer.parse_plotly_selection(SimpleNamespace(selection=None))
    assert empty is None


def test_plotly_multi_point_indices_are_sorted_and_unique() -> None:
    event = SimpleNamespace(
        selection=SimpleNamespace(
            points=[
                {"customdata": [7, "A", 1, 8, 8, "", "P"]},
                {"customdata": [2, "A", 1, 3, 3, "", "P"]},
                {"customdata": [7, "A", 1, 8, 8, "", "P"]},
                {"customdata": [None, "A", 1, 1, 1, "", "P"]},
            ]
        )
    )
    indices = structure_viewer.parse_plotly_selection_indices(event)
    assert indices == [2, 7]


def test_stale_hash_rejects_structure_scene() -> None:
    rep = dna_3d.representation_for_sequence(DICKERSON)
    scene = structure_scene.build_scene(
        rep["envelope"],
        expected_sequence_hash=provenance.sequence_digest("ATGCATGCATGC"),
    )
    assert scene["status"] == "STALE"


def test_lod_low_points_are_subset_of_high_coordinates() -> None:
    rep = dna_3d.representation_for_sequence(DICKERSON)
    envelope = rep["envelope"]
    digest = str(rep["sequence_hash"])
    high = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=digest,
        lod=region_nav.LOD_HIGH,
    )
    low = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=digest,
        lod=region_nav.LOD_LOW,
    )
    assert high["status"] == "READY"
    assert low["status"] == "READY"
    assert low["visual_lod"] == region_nav.LOD_LOW
    assert high["visual_stride"] == 1
    assert low["visual_stride"] > 1
    assert int(low["n_backbone_points"]) <= int(high["n_backbone_points"])
    high_xyz = {
        (round(x, 5), round(y, 5), round(z, 5))
        for trace in list((high.get("objects") or {}).get("backbone_traces") or [])
        for x, y, z in zip(trace.get("x") or [], trace.get("y") or [], trace.get("z") or [])
        if x is not None and y is not None and z is not None
    }
    low_xyz = {
        (round(x, 5), round(y, 5), round(z, 5))
        for trace in list((low.get("objects") or {}).get("backbone_traces") or [])
        for x, y, z in zip(trace.get("x") or [], trace.get("y") or [], trace.get("z") or [])
        if x is not None and y is not None and z is not None
    }
    assert low_xyz
    assert low_xyz.issubset(high_xyz)
    assert int(high["n_atoms_source"]) == len(list(envelope.get("atoms") or []))
    assert int(low["n_atoms_source"]) == int(high["n_atoms_source"])


def test_error_badges_are_distinct() -> None:
    html_not_found = components.status_badge("NOT_FOUND")
    html_network = components.status_badge("NETWORK_ERROR")
    html_timeout = components.status_badge("TIMEOUT")
    html_parse = components.status_badge("PARSING_ERROR")
    html_limit = components.status_badge("RESOURCE_LIMIT")
    html_none = components.status_badge("NO_STRUCTURE")
    html_empty = components.status_badge("INSUFFICIENT_DATA")
    assert "NOT FOUND" in html_not_found
    assert "NETWORK ERROR" in html_network
    assert "TIMEOUT" in html_timeout
    assert "PARSING ERROR" in html_parse
    assert "RESOURCE LIMIT" in html_limit
    assert "NO STRUCTURE" in html_none
    assert "INSUFFICIENT DATA" in html_empty
    assert html_not_found != html_network
    assert html_not_found != html_none


def test_msa_gap_column_is_unmapped() -> None:
    members = [
        msa.build_collection_member(
            sequence="ACGT",
            identifier="a",
            source="test fixture",
            molecule="DNA",
        ),
        msa.build_collection_member(
            sequence="ACCGT",
            identifier="b",
            source="test fixture",
            molecule="DNA",
        ),
    ]
    result = msa.build_msa_result(
        raw_alignment=">a\nAC-GT\n>b\nACCGT\n",
        members=members,
        tool="test_fixture",
        tool_version="phase11",
        method="fixture aligned FASTA (not Clustal Omega)",
        parameters={"source": "tests/test_phase11_flows.py"},
        source="test fixture",
    )
    gap = molecule_selection.selection_from_msa_column(
        result, 2, sequence="ACGT", molecule="DNA"
    )
    assert gap["status"] == "UNMAPPED"
    assert gap["query_indices"] == []
    mapped = molecule_selection.selection_from_msa_column(
        result, 0, sequence="ACGT", molecule="DNA"
    )
    assert mapped["status"] in {"READY", "UNMAPPED"}
    assert mapped["query_indices"] == [0]


def test_apptest_protein_bundled_1crn_search_and_view() -> None:
    app = open_module("protein", timeout=120)
    assert not app.exception
    app.text_area(key="prot_text").set_value(CRAMBIN).run()
    app.button(key="prot_analyze").click().run()
    assert not app.exception
    app.toggle(key="prot_structure").set_value(True).run()
    assert not app.exception
    app.text_input(key="prot_struct_pdb").set_value("1CRN").run()
    app.button(key="prot_struct_search_btn").click().run()
    assert not app.exception
    assert "prot_struct_hits" in app.session_state
    hits = app.session_state["prot_struct_hits"]
    assert isinstance(hits, dict)
    bundled = [
        item
        for item in list(hits.get("hits") or [])
        if item.get("bundled_fixture") == "1CRN.cif"
    ]
    assert bundled, hits
    app.button(key="prot_struct_load").click().run()
    assert not app.exception
    assert "prot_struct_result" in app.session_state
    loaded = app.session_state["prot_struct_result"]
    assert isinstance(loaded, dict)
    assert loaded.get("structure_id") == "1CRN"
    assert loaded.get("kind") == "experimental"
    app.toggle(key="prot_view_3d").set_value(True).run()
    assert not app.exception
    combined = _ui_text(app)
    assert "1CRN" in combined
    assert "EXPERIMENTAL" in combined or "experimental" in combined.lower()
    assert "Insufficient data" not in combined


def test_apptest_dna_3d_selection_widget_still_renders() -> None:
    app = open_module("dna", timeout=90)
    app.text_area(key="dna_text").set_value(DICKERSON).run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_3d").set_value(True).run()
    app.toggle(key="dna_illust_view_3d").set_value(True).run()
    assert not app.exception
    combined = _ui_text(app)
    assert "1BNA" in combined
    assert "3D viewport (Plotly Scatter3d)" in combined
    assert "Click a mapped backbone point" in combined
    assert provenance.HELIXSCOPE_VERSION.startswith("0.")
    assert "3D scene summary" in combined or "helixscope_dna_illust" in combined or "sequence_hash" in combined


def test_apptest_dna_stale_replaces_1bna_with_illustrative() -> None:
    app = open_module("dna", timeout=90)
    app.text_area(key="dna_text").set_value(DICKERSON).run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_3d").set_value(True).run()
    app.toggle(key="dna_illust_view_3d").set_value(True).run()
    assert "1BNA" in _ui_text(app)
    app.text_area(key="dna_text").set_value("ATGCATGCATGC").run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_3d").set_value(True).run()
    app.toggle(key="dna_illust_view_3d").set_value(True).run()
    assert not app.exception
    combined = _ui_text(app)
    assert "1BNA" not in combined
    assert "ILLUSTRATIVE" in combined


def test_apptest_50k_dna_analyze_keeps_full_length() -> None:
    app = open_module("dna", timeout=180)
    app.text_area(key="dna_text").set_value(SEQ_50K).run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    stored = app.session_state["dna_input"] if "dna_input" in app.session_state else ""
    assert len(str(stored)) == 50_000
    combined = _ui_text(app)
    assert "50,000" in combined or "50000" in combined
    assert "silently truncated" not in combined.lower()


def _protein_search_error_text(monkeypatch: pytest.MonkeyPatch, category: str, message: str) -> str:
    def boom(*_args, **_kwargs):
        raise protein_structure.StructureError(message, category)

    monkeypatch.setattr(protein_structure, "fetch_rcsb_entry", boom)
    app = open_module("protein", timeout=90)
    app.text_area(key="prot_text").set_value(CRAMBIN).run()
    app.button(key="prot_analyze").click().run()
    app.toggle(key="prot_structure").set_value(True).run()
    app.text_input(key="prot_struct_pdb").set_value("9ZZZ").run()
    app.button(key="prot_struct_search_btn").click().run()
    assert not app.exception
    return _ui_text(app)


def test_apptest_not_found_is_not_no_structure(monkeypatch: pytest.MonkeyPatch) -> None:
    text = _protein_search_error_text(
        monkeypatch, "NOT_FOUND", "Structure identifier was not found."
    )
    assert "NOT_FOUND" in text or "Not found" in text
    assert "Structure retrieval failed" in text
    collapsed = "No validated structure found" in text and "NOT_FOUND" not in text and "Not found" not in text
    assert not collapsed


def test_apptest_network_error(monkeypatch: pytest.MonkeyPatch) -> None:
    text = _protein_search_error_text(monkeypatch, "NETWORK_ERROR", "connection refused")
    assert "NETWORK_ERROR" in text or "Network error" in text
    assert "NOT_FOUND" not in text


def test_apptest_timeout_error(monkeypatch: pytest.MonkeyPatch) -> None:
    text = _protein_search_error_text(monkeypatch, "TIMEOUT", "timed out")
    assert "TIMEOUT" in text or "Timeout" in text


def test_apptest_parsing_error(monkeypatch: pytest.MonkeyPatch) -> None:
    text = _protein_search_error_text(monkeypatch, "PARSING_ERROR", "Remote JSON is malformed.")
    assert "PARSING_ERROR" in text or "Parsing error" in text


def test_apptest_resource_limit_error(monkeypatch: pytest.MonkeyPatch) -> None:
    text = _protein_search_error_text(
        monkeypatch, "RESOURCE_LIMIT", "Structure exceeds the byte limit."
    )
    assert "RESOURCE_LIMIT" in text or "Resource limit" in text


def test_apptest_crispr_stale_complex_hash() -> None:
    app = open_module("crispr", timeout=90)
    target = "GAGTCCGAGCAGAAGAAGAAGGG"
    app.text_area(key="crispr_text").set_value(target).run()
    app.session_state["crispr_complex"] = {
        "status": "RETRIEVED",
        "bound_guide_hash": "",
        "bound_target_hash": provenance.sequence_digest("AAAAAAAAAAAAAAAAAAAAAAA"),
        "disclaimer": "test fixture stale complex",
        "chains": [],
    }
    app.run()
    assert not app.exception
    combined = _ui_text(app)
    assert "STALE" in combined
    assert "stale/incompatible structure" in combined.lower()


def test_apptest_msa_column_highlight_and_gap() -> None:
    members = [
        msa.build_collection_member(
            sequence=CRAMBIN,
            identifier="crambin",
            source="test fixture",
            molecule="PROTEIN",
        ),
        msa.build_collection_member(
            sequence=CRAMBIN[1:],
            identifier="short",
            source="test fixture",
            molecule="PROTEIN",
        ),
    ]
    aligned = f">crambin\n{CRAMBIN}\n>short\n-{CRAMBIN[1:]}\n"
    msa_result = msa.build_msa_result(
        raw_alignment=aligned,
        members=members,
        tool="test_fixture",
        tool_version="phase11",
        method="fixture aligned FASTA (not Clustal Omega)",
        parameters={"source": "tests/test_phase11_flows.py"},
        source="test fixture",
    )
    app = open_module("protein", timeout=120)
    app.text_area(key="prot_text").set_value(CRAMBIN).run()
    app.button(key="prot_analyze").click().run()
    app.session_state["msa_collection"] = members
    app.session_state["msa_result"] = msa_result
    app.session_state["msa_column"] = 0
    switch_module(app, "msa")
    assert not app.exception
    keys = [item.key for item in app.button]
    assert "msa_3d_highlight" in keys
    app.button(key="msa_3d_highlight").click().run()
    assert not app.exception
    assert "molecule_selection" in app.session_state
    selection = app.session_state["molecule_selection"]
    assert selection.get("source") == "msa"
    assert selection.get("query_indices") == [0]
    app.session_state["msa_column"] = 0
    app.run()
    gap_sel = molecule_selection.selection_from_msa_column(
        msa_result, 0, sequence=CRAMBIN[1:], molecule="PROTEIN"
    )
    assert gap_sel["status"] == "UNMAPPED"


def test_apptest_phylogeny_leaf_uses_dna_3d_hash() -> None:
    members = [
        msa.build_collection_member(
            sequence=DICKERSON,
            identifier="dickerson",
            source="test fixture",
            molecule="DNA",
        ),
        msa.build_collection_member(
            sequence="CGCGAATTCGCA",
            identifier="mut_a",
            source="test fixture",
            molecule="DNA",
        ),
        msa.build_collection_member(
            sequence="CGCGAATTCGCT",
            identifier="mut_t",
            source="test fixture",
            molecule="DNA",
        ),
    ]
    aligned = (
        f">dickerson\n{DICKERSON}\n>mut_a\nCGCGAATTCGCA\n>mut_t\nCGCGAATTCGCT\n"
    )
    msa_result = msa.build_msa_result(
        raw_alignment=aligned,
        members=members,
        tool="test_fixture",
        tool_version="phase11",
        method="fixture aligned FASTA (not Clustal Omega)",
        parameters={"source": "tests/test_phase11_flows.py"},
        source="test fixture",
    )
    app = open_module("dna", timeout=120)
    app.text_area(key="dna_text").set_value(DICKERSON).run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_3d").set_value(True).run()
    app.session_state["msa_collection"] = members
    app.session_state["msa_result"] = msa_result
    switch_module(app, "phylogeny")
    assert not app.exception
    app.button(key="phylo_build").click().run()
    assert not app.exception
    assert "phylo_result" in app.session_state
    tree = app.session_state["phylo_result"]
    assert isinstance(tree, dict)
    digest = provenance.sequence_digest(DICKERSON)
    leaves = list(tree.get("leaves") or [])
    hashes = [str(item.get("sequence_hash") or "") for item in leaves]
    assert digest in hashes
    target_id = next(
        str(item.get("tree_id") or "")
        for item in leaves
        if str(item.get("sequence_hash") or "") == digest
    )
    app.selectbox(key="phylo_leaf_select").set_value(target_id).run()
    assert not app.exception
    combined = _ui_text(app)
    assert "Validated structure in this session" in combined or "1BNA" in combined
    assert "taxonomic tree" in combined.lower() or "Not a taxonomic" in combined


def test_ncbi_nucleotide_feature_is_not_protein_residue() -> None:
    spans = [{"start": 0, "end": 4, "status": "AVAILABLE"}]
    protein_sel = molecule_selection.selection_from_ncbi_spans(
        spans, sequence=CRAMBIN, molecule="PROTEIN", accession="NM_000000"
    )
    assert protein_sel["molecule"] == "PROTEIN"
    dna_sel = molecule_selection.selection_from_ncbi_spans(
        spans, sequence=DICKERSON, molecule="DNA", accession="NM_000000"
    )
    assert dna_sel["molecule"] == "DNA"
    assert dna_sel["query_indices"] == [0, 1, 2, 3]
    assert not molecule_selection.molecules_are_compatible("DNA", "PROTEIN")
