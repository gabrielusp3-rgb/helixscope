"""Phase 19 Part 1: navigation registry, shell, session state, smoke."""

from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest

from tests.helix_apptest import open_module, switch_module
from ui.navigation import (
    ACTIVE_MODULE_KEY,
    GROUP_ORDER,
    MODULES,
    MODULE_BY_ID,
    grouped_modules,
    module_ids,
)

ROOT = Path(__file__).resolve().parents[1]
DICKERSON = "CGCGAATTCGCG"


def _ui_text(app) -> str:
    return " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption", "info", "error", "warning")
        for item in list(getattr(app, bucket, []) or [])
    )


def test_cold_start_back_returns_to_overview() -> None:
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=90)
    app.run()
    assert not app.exception
    app.button(key="nav_mod_dna").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "dna"
    assert app.session_state["helix_module_history"] == ["overview"]
    app.button(key="helix_back").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "overview"
    assert not app.exception


def test_back_button_returns_to_the_previous_module() -> None:
    app = open_module("overview", timeout=90)
    assert not app.exception
    assert app.button(key="helix_back") is not None
    assert app.button(key="helix_forward") is not None
    app.button(key="nav_mod_dna").click().run()
    app.button(key="nav_mod_protein").click().run()
    app.button(key="nav_mod_rna").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "rna"
    app.button(key="helix_back").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "protein"
    app.button(key="helix_back").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "dna"
    app.button(key="helix_forward").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "protein"
    assert not app.exception


def test_module_switch_keeps_the_previous_analysis() -> None:
    app = open_module("dna", timeout=90)
    app.session_state["dna_input"] = "ATGCATGCATGC"
    app.session_state["protein_input"] = "ACDEFGHIKL"
    app.button(key="nav_mod_protein").click().run()
    assert app.session_state["dna_input"] == "ATGCATGCATGC"
    app.button(key="nav_mod_rna").click().run()
    app.button(key="helix_back").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "protein"
    assert app.session_state["protein_input"] == "ACDEFGHIKL"
    app.button(key="helix_back").click().run()
    assert app.session_state[ACTIVE_MODULE_KEY] == "dna"
    assert app.session_state["dna_input"] == "ATGCATGCATGC"
    assert not app.exception


def test_registry_ids_are_unique_and_grouped() -> None:
    ids = module_ids()
    assert len(ids) == len(set(ids))
    assert ids[0] == "overview"
    assert "dna" in ids
    assert "variant" in ids
    assert "compare" in ids
    assert "phylogeny" in ids
    groups = [group for group, _ in grouped_modules()]
    assert groups == [item for item in GROUP_ORDER if item in groups]
    for spec in MODULES:
        assert spec.renderer
        assert spec.title
        assert spec.method
        assert spec.scope
        assert "AI Homology" not in spec.title
        assert spec.id in MODULE_BY_ID


def test_ui_layer_does_not_recompute_science() -> None:
    forbidden = (
        "gc_content(",
        "melting_temperature(",
        "infer_phylogeny(",
        "run_local_blast(",
        "submit_alignment(",
    )
    for path in (ROOT / "ui").glob("*.py"):
        text = path.read_text(encoding="utf-8").lower()
        for token in forbidden:
            assert token not in text, f"{path.name} contains {token}"


def test_app_source_retired_global_tab_strip() -> None:
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "st.tabs(" not in text
    assert "_render_hero(" not in text
    assert "helix-hero-modules" not in text
    assert "render_sidebar_navigation" in text
    assert "persist_orphaned_widget_keys" in text


def test_apptest_overview_and_contextual_header() -> None:
    app = open_module("overview", timeout=90)
    assert not app.exception
    combined = _ui_text(app)
    assert "Workspace overview" in combined
    assert "Active scientific objects" in combined
    assert "Scientific method" in combined or "Scope and limitations" in combined
    keys = [item.key for item in app.button]
    assert "nav_mod_dna" in keys
    assert "nav_mod_compare" in keys


def test_apptest_module_switch_preserves_dna_input() -> None:
    app = open_module("dna", timeout=90)
    assert not app.exception
    app.text_area(key="dna_text").set_value(DICKERSON).run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    stored = str(app.session_state["dna_input"])
    switch_module(app, "protein")
    assert not app.exception
    assert str(app.session_state["dna_input"]) == stored
    switch_module(app, "dna")
    assert not app.exception
    assert str(app.session_state["dna_input"]) == stored
    combined = _ui_text(app)
    assert "DNA analysis" in combined or "Independent of the RNA" in combined


def test_apptest_variant_compare_variant_keeps_controls() -> None:
    app = open_module("variant", timeout=90)
    assert not app.exception
    toggles = [item.key for item in app.toggle]
    assert "variant_use_vep" in toggles
    assert "variant_use_clinvar" in toggles
    switch_module(app, "compare")
    assert not app.exception
    radios = [item.key for item in app.radio]
    assert "compare_mode" in radios
    switch_module(app, "variant")
    assert not app.exception
    toggles = [item.key for item in app.toggle]
    assert "variant_use_vep" in toggles


def test_apptest_msa_phylogeny_msa_roundtrip() -> None:
    app = open_module("msa", timeout=90)
    assert not app.exception
    assert "msa_fasta_text" in [item.key for item in app.text_area]
    switch_module(app, "phylogeny")
    assert not app.exception
    combined = _ui_text(app).lower()
    assert "no completed msa" in combined
    assert "phylogeny engines are not missing" in combined or "not invented" in combined
    switch_module(app, "msa")
    assert not app.exception
    assert "msa_add_fasta" in [item.key for item in app.button]


def test_apptest_dna_analyze_smoke() -> None:
    app = open_module("dna", timeout=90)
    app.text_area(key="dna_text").set_value("ATGGCATTACGTACGTACGT").run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    assert "dna_input" in app.session_state
    combined = _ui_text(app)
    assert "COMPUTED" in combined


def test_apptest_alignment_smoke() -> None:
    app = open_module("alignment", timeout=90)
    assert not app.exception
    app.text_area(key="align_seq1").set_value("ATGCATGCATGC").run()
    app.text_area(key="align_seq2").set_value("ATGCATGCATGC").run()
    app.button(key="align_button").click().run()
    assert not app.exception
    assert "align_result" in app.session_state
    result = app.session_state["align_result"]
    assert isinstance(result, dict)


def test_apptest_motif_smoke() -> None:
    app = open_module("motif", timeout=90)
    assert not app.exception
    widget = next(item.key for item in app.text_area if str(item.key).startswith("motif_text_"))
    app.text_area(key=widget).set_value("AAAGAATTCAAA").run()
    app.text_input(key="motif_pattern").set_value("GAATTC").run()
    app.button(key="motif_search_button").click().run()
    assert not app.exception


def test_apptest_compare_keeps_alignment_method() -> None:
    app = open_module("compare", timeout=90)
    assert not app.exception
    combined = _ui_text(app)
    assert "fatcat" in combined.lower() or "Alignment method" in combined
    assert "compare_mode" in [item.key for item in app.radio]


def test_apptest_all_registry_modules_open() -> None:
    app = open_module("overview", timeout=90)
    assert not app.exception
    for module_id in module_ids():
        switch_module(app, module_id)
        assert not app.exception, module_id
        assert app.session_state[ACTIVE_MODULE_KEY] == module_id
        keys = [item.key for item in app.button]
        assert f"nav_mod_{module_id}" in keys
