"""Fase 13: MSA prealigned, UI da arvore, BLAST+ local sem DB, SASA cache."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from modules import blast_search, msa, phylogeny, protein_structure, provenance, runtime, tool_paths
from tests.helix_apptest import open_module

ROOT = Path(__file__).resolve().parents[1]


def _ui_text(app: AppTest) -> str:
    return " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption", "info", "error", "warning")
        for item in list(getattr(app, bucket, []) or [])
    )


def test_import_prealigned_fasta_is_completed_not_clustal():
    text = ">seq_a\nAAAAAAAAAA\n>seq_b\nAAAAAAAATT\n>seq_c\nTTTTTTTTTT\n"
    result = msa.import_prealigned_fasta(text)
    assert result["status"] == "COMPLETED"
    assert "prealigned" in str(result.get("method") or "").lower()
    assert "user-supplied" in str(result.get("tool") or "").lower()
    assert result["n_sequences"] == 3
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        rooting=phylogeny.ROOTING_UNROOTED,
    )
    assert tree["status"] == "COMPUTED"
    assert tree["n_leaves"] == 3


def test_prealigned_rejects_mixed_molecules():
    text = ">dna\nACGTACGTACGTACGT\n>prot\nMKTFFVLVVL\n"
    with pytest.raises(msa.MsaError) as exc:
        msa.import_prealigned_fasta(text)
    assert exc.value.category == "INVALID_INPUT"


def test_apptest_tree_button_hidden_without_completed_msa():
    app = open_module("msa", timeout=90)
    keys = [item.key for item in app.button]
    assert "phylo_build" not in keys
    combined = _ui_text(app).lower()
    assert "no completed msa" in combined
    assert "phylogeny engines are not missing" in combined or "not invented" in combined


def test_apptest_tree_button_visible_with_completed_msa():
    text = ">seq_a\nAAAAAAAAAA\n>seq_b\nAAAAAAAATT\n>seq_c\nTTTTTTTTTT\n"
    result = msa.import_prealigned_fasta(text)
    app = open_module("msa", timeout=90)
    app.session_state["msa_collection"] = list(result.get("input_members") or [])
    app.session_state["msa_result"] = result
    app.run()
    assert not app.exception
    keys = [item.key for item in app.button]
    assert "phylo_build" in keys
    assert "msa_load_prealigned" in keys


def test_local_blast_without_db_raises_unavailable(monkeypatch):
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    monkeypatch.setattr(blast_search, "declared_blast_prefix", lambda program="blastn": "")
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.run_local_blast(program="blastn", query="ACGTACGTACGTACGT", molecule="DNA")
    assert exc.value.category in {"UNAVAILABLE", "INVALID_DATABASE"}


def test_sasa_cache_hits_on_structure_hash():
    parsed = {
        "content_hash": "phase13-sasa-cache",
        "atoms": [
            {
                "group": "ATOM",
                "atom_name": "C",
                "element": "C",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "auth_asym_id": "A",
                "auth_seq_id": 1,
            }
        ],
    }
    first = protein_structure.shrake_rupley_sasa(parsed)
    second = protein_structure.shrake_rupley_sasa(parsed)
    assert first["status"] == "COMPUTED"
    assert second["cache_status"] == "cached"
    assert first["sasa_angstrom2"] == second["sasa_angstrom2"]


def test_sasa_neighbor_list_matches_isolated_carbon():
    parsed = {
        "atoms": [
            {
                "group": "ATOM",
                "atom_name": "C",
                "element": "C",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "auth_asym_id": "A",
                "auth_seq_id": 1,
            }
        ]
    }
    result = protein_structure.shrake_rupley_sasa(parsed)
    radius = 1.70 + 1.4
    expected = 4.0 * math.pi * radius * radius
    assert result["sasa_angstrom2"] == pytest.approx(round(expected, 3), abs=0.001)
    assert result["neighbor_list"] is True


def test_tool_paths_scan_is_bounded():
    assert tool_paths.MAX_SCAN_DEPTH <= 3
    source = Path(tool_paths.__file__).read_text(encoding="utf-8")
    assert "C:\\\\" not in source or "no" in source.lower()
    assert "os.walk" in source
    assert "repository_tools_dir" in source


def test_runtime_wired_in_app_source():
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "runtime.ensure_runtime" in text
    assert "importlib.reload(structure_viewer)" not in text
    assert "dna_analysis = importlib.reload(dna_analysis)" not in text


def test_dssp_cache_hits_on_structure_hash_and_version(monkeypatch):
    from types import SimpleNamespace

    fixture = (ROOT / "tests" / "fixtures" / "dssp_legacy_min.txt").read_text(encoding="utf-8")
    calls = {"n": 0}

    def fake_run(argv, **kwargs):
        calls["n"] += 1
        return SimpleNamespace(returncode=0, stdout=fixture, stderr="mkdssp 4.6.1")

    monkeypatch.setattr(
        protein_structure,
        "dssp_availability",
        lambda: {"available": True, "version": "4.6.1", "reason": "mock detector"},
    )
    protein_structure._DSSP_CACHE.clear()
    text = "data_mock\n_atom_site.id 1\n"
    first = protein_structure.assign_secondary_structure_dssp(text, run_fn=fake_run)
    second = protein_structure.assign_secondary_structure_dssp(text, run_fn=fake_run)
    assert first["status"] == "EXPERIMENTAL"
    assert first["cache_status"] == "live"
    assert second["cache_status"] == "cached"
    assert calls["n"] == 1
    assert first["version"] == "4.6.1"


def test_blast_cache_key_changes_with_backend_db_and_tool_version():
    remote = blast_search.cache_key(
        query_hash="abc",
        program="blastn",
        database="core_nt",
        expect=10.0,
        hitlist_size=20,
    )
    local_a = blast_search.cache_key(
        query_hash="abc",
        program="blastn",
        database="core_nt",
        expect=10.0,
        hitlist_size=20,
        backend="blastplus_local",
        database_identity="tiny-test-db",
        tool_version="2.17.0",
    )
    local_b = blast_search.cache_key(
        query_hash="abc",
        program="blastn",
        database="core_nt",
        expect=10.0,
        hitlist_size=20,
        backend="blastplus_local",
        database_identity="tiny-test-db",
        tool_version="2.16.0",
    )
    assert remote != local_a
    assert local_a != local_b
