"""Monochrome liquid-glass identity: tokens, chrome, frozen science."""

from __future__ import annotations

from pathlib import Path

from ui.tokens import (
    HS_CYAN,
    HS_ELECTRIC_BLUE,
    HS_ICE_BLUE,
    HS_PRIMARY,
    HS_VOID,
    canonical_status,
    css_variables,
    status_style,
)
from ui.workspace import metric_grid_html

ROOT = Path(__file__).resolve().parents[1]


def test_monochrome_palette_tokens_are_centralized() -> None:
    css = css_variables()
    assert HS_VOID == "#000000"
    assert HS_PRIMARY == "#FFFFFF"
    assert HS_ELECTRIC_BLUE == "#2355A0"
    assert HS_ICE_BLUE == "#91B4E4"
    assert HS_CYAN == HS_ICE_BLUE
    assert HS_CYAN.lower() != "#22d3ee"
    for token in (
        "--hs-void:",
        "--hs-bg-deep:",
        "--hs-bg-raised:",
        "--hs-blue-deep:",
        "--hs-blue:",
        "--hs-blue-mid:",
        "--hs-ice:",
        "--hs-frost:",
        "--hs-white-blue:",
        "--hs-glass-thin:",
        "--hs-glass-medium:",
        "--hs-glass-thick:",
        "--hs-glass-fallback:",
        "--hs-glass-border:",
        "--hs-glass-highlight:",
        "--hs-glass-shadow:",
        "--hs-highlight:",
        "--hs-border-hover:",
        "--hs-cyan:",
        "--hs-text:",
        "--hs-primary:",
    ):
        assert token in css
    assert "--hs-void: #000000;" in css
    assert "--hs-primary: #FFFFFF;" in css
    assert "--hs-blue: #2355A0;" in css


def test_status_semantics_remain_distinguishable() -> None:
    assert status_style("EXPERIMENTAL")["fg"] != status_style("PREDICTED")["fg"]
    assert status_style("PREDICTED")["fg"] != status_style("ILLUSTRATIVE")["fg"]
    assert status_style("UNAVAILABLE")["fg"] != status_style("ERROR")["fg"]
    assert status_style("LIVE_VALIDATED")["fg"] != status_style("REMOTE_VALIDATED")["fg"]
    assert canonical_status("RESOURCE_LIMIT") == "RESOURCE_LIMIT"


def test_metric_grid_geometry_is_preserved() -> None:
    html = metric_grid_html(
        [
            ("Length", "12", "bp"),
            ("GC", "50.00", "%"),
            ("AT", "50.00", "%"),
            ("Tm", "Unavailable", ""),
            ("MW", "1.00", "kDa"),
            ("Entropy", "2.00", "bits"),
            ("AT skew", "0.00", ""),
            ("GC skew", "0.00", ""),
            ("SantaLucia Tm", "46.79", "C"),
        ],
        columns=4,
    )
    assert 'data-cols="4"' in html
    assert "--hs-metric-cols:4" in html
    assert "auto-fit" not in html
    assert html.count('class="hs-metric"') == 9
    assert "Unavailable" in html
    styles = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    assert "repeat(var(--hs-metric-cols, 4), minmax(0, 1fr))" in styles


def test_chrome_keeps_command_key_and_content_solid() -> None:
    styles = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    pointer = (ROOT / "ui" / "pointer.py").read_text(encoding="utf-8")
    assert ".st-key-helix_cmd_open" in styles
    assert "justify-content: flex-end" in styles
    assert "backdrop-filter: none" in styles
    assert "@supports not" in styles
    assert "prefers-reduced-motion" in styles
    assert "prefers-reduced-transparency" in styles
    assert "particle-js" not in styles.lower()
    assert "animated dna wallpaper" not in styles.lower()
    assert "st.components.v2.component" in pointer
    assert "isolate_styles=False" in pointer
    assert "components.v1.html" not in styles
    glass_card_block = styles.split(".glass-card, .helix-panel")[1].split(".hs-glass")[0]
    assert "backdrop-filter: none" in glass_card_block


def test_streamlit_production_pin_unchanged() -> None:
    req = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    cfg = (ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8")
    assert "streamlit==1.58.0" in req
    assert "biopython==1.87" in req
    assert 'toolbarMode = "viewer"' in cfg
    assert 'backgroundColor = "#000000"' in cfg
    assert 'primaryColor = "#FFFFFF"' in cfg
    assert "#00f2ff" not in cfg.lower()


def test_tokens_are_not_scientific_engines() -> None:
    text = (ROOT / "ui" / "tokens.py").read_text(encoding="utf-8").lower()
    assert "gc_content(" not in text
    assert "infer_phylogeny(" not in text
    assert "doench" not in text
    assert "santalucia" not in text


def test_hand_control_defaults_remain_off_and_camera_only() -> None:
    hand = (ROOT / "ui" / "hand_control.py").read_text(encoding="utf-8")
    js = (ROOT / "ui" / "hand_control.js").read_text(encoding="utf-8")
    assert 'value=False' in hand
    assert "Gestures move the 3D camera, not atomic coordinates." in hand
    assert "trace.x" not in js
    assert ".st-key-" in js
    assert "plotHostKey" in js
    assert "Allow open-palm camera reset" in hand
    assert "value=False" in hand


def test_crispr_cold_start_copy_not_reintroduced() -> None:
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    command = (ROOT / "ui" / "command_palette.py").read_text(encoding="utf-8")
    assert "placeholder=crispr.EXAMPLE_EMX1_SPCAS9" not in text
    assert "Paste target DNA" in text
    assert 'key="helix_cmd_open"' in command


def test_3d_highlight_is_ice_not_turquoise() -> None:
    from modules.structure_scene import HIGHLIGHT_COLOR, CHAIN_PALETTE

    assert HIGHLIGHT_COLOR.lower() != "#00f2ff"
    assert "#22d3ee" not in {item.lower() for item in CHAIN_PALETTE}
    assert HIGHLIGHT_COLOR == "#AFCBF1"
