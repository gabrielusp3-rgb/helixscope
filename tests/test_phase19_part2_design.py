"""Phase 19 Part 2: design tokens, command palette, workspace HTML, no frontend science."""

from __future__ import annotations

from pathlib import Path

from ui.command_palette import registered_command_targets
from ui.tokens import STATUS_COLORS, canonical_status, css_variables, status_style
from ui.workspace import (
    chain_role_legend_html,
    compare_hero_html,
    job_state_html,
    unavailable_text,
    variant_journey_html,
)
from tests.helix_apptest import open_module

ROOT = Path(__file__).resolve().parents[1]


def test_status_tokens_distinguish_experimental_predicted_illustrative() -> None:
    exp = status_style("EXPERIMENTAL")
    pred = status_style("PREDICTED")
    ill = status_style("ILLUSTRATIVE")
    unav = status_style("UNAVAILABLE")
    err = status_style("ERROR")
    assert exp["fg"] != pred["fg"]
    assert pred["fg"] != ill["fg"]
    assert unav["fg"] != err["fg"]
    assert canonical_status("EXPERIMENTAL") == "EXPERIMENTAL"
    assert "LIVE_VALIDATED" in STATUS_COLORS
    assert "REMOTE_VALIDATED" in STATUS_COLORS
    assert "TEST_ONLY" in STATUS_COLORS
    assert canonical_status("TEST_ONLY") == "TEST_ONLY"
    assert canonical_status("NOT_INSTALLED") == "NOT_INSTALLED"
    css = css_variables()
    assert "--hs-cyan:" in css
    assert "--hs-status-experimental-fg:" in css
    assert "--hs-status-test-only-fg:" in css


def test_unavailable_text_does_not_show_zero() -> None:
    assert unavailable_text(None) == "Unavailable"
    assert unavailable_text("") == "Unavailable"
    assert unavailable_text(0) == "0"


def test_compare_hero_keeps_distinct_rmsd_fields() -> None:
    html = compare_hero_html(
        {
            "reference": {"entry_id": "8HSK", "kind_label": "EXPERIMENTAL"},
            "target": {"entry_id": "8HSF", "kind_label": "EXPERIMENTAL"},
            "method_display": "jFATCAT-rigid",
            "engine": "RCSB",
            "status": "REMOTE_VALIDATED",
            "rmsd_global_angstrom": 0.99,
            "rmsd_block0_angstrom": 0.8,
            "n_aligned_residue_pairs": 20,
            "aln_coverage_percent": [95, 100],
            "tm_scores": [{"type": "TM-score", "value": 0.23}],
        }
    )
    assert "0.99" in html
    assert "0.8" in html
    assert "0.9 A" not in html
    assert "20" in html
    assert "Aligned residue pairs" in html
    assert "Block RMSD" in html
    assert "0.23" in html
    assert "95%" in html


def test_variant_journey_breaks_on_unmapped() -> None:
    html = variant_journey_html(
        (
            ("Genome", "GRCh38 17:1 A>G", "COMPUTED"),
            ("Protein", "", "UNMAPPED"),
        )
    )
    assert "UNMAPPED" in html
    assert "hs-journey-break" in html
    assert "hs-journey-gap" in html


def test_chain_role_legend_uses_deposited_roles_only() -> None:
    html = chain_role_legend_html(
        (
            {"chain_id": "B", "role": "cas_protein"},
            {"chain_id": "A", "role": "guide_rna"},
            {"chain_id": "C", "role": "target_dna"},
        )
    )
    assert "Cas protein" in html
    assert "Guide RNA" in html
    assert "Target DNA" in html
    assert "B" in html
    assert chain_role_legend_html([]) == ""


def test_job_state_html_keeps_status_text() -> None:
    html = job_state_html("RUNNING", "RID ABC123")
    assert "RUNNING" in html
    assert "ABC123" in html
    waiting = job_state_html("WAITING", "")
    assert "WAITING" in waiting


def test_command_palette_targets_match_registry() -> None:
    ids = registered_command_targets()
    assert "dna" in ids
    assert "variant" in ids
    assert "compare" in ids
    assert "phylogeny" in ids


def test_pointer_script_migrated_off_components_v1() -> None:
    styles = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    pointer = (ROOT / "ui" / "pointer.py").read_text(encoding="utf-8")
    assert "components.v1.html" not in styles
    assert "st.components.v2.component" in pointer
    assert "Trusted application script" in pointer


def test_ui_tokens_are_not_scientific_engines() -> None:
    text = (ROOT / "ui" / "tokens.py").read_text(encoding="utf-8").lower()
    assert "gc_content(" not in text
    assert "infer_phylogeny(" not in text


def test_apptest_command_button_and_dna_header() -> None:
    app = open_module("overview", timeout=90)
    assert not app.exception
    keys = [item.key for item in app.button]
    assert "helix_cmd_open" in keys
    labels = [str(getattr(item, "label", "") or getattr(item, "value", "") or "") for item in app.button]
    combined_labels = " ".join(labels)
    assert "Pesquisar" in combined_labels or any("Pesquisar" in str(item) for item in app.button)
    app = open_module("dna", timeout=90)
    assert not app.exception
    app.text_area(key="dna_text").set_value("ATGGCATTACGTACGTACGT").run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption")
        for item in list(getattr(app, bucket, []) or [])
    )
    assert "COMPUTED" in combined
    assert "bp" in combined.lower() or "DNA analysis" in combined
