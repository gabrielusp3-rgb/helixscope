"""CRISPR cold start: no automatic EMX1 target or name."""

from __future__ import annotations

from pathlib import Path

from modules import crispr
from tests.helix_apptest import open_module, switch_module

ROOT = Path(__file__).resolve().parents[1]


def test_emx1_fixture_remains_available_for_science() -> None:
    assert crispr.EXAMPLE_EMX1_SPCAS9 == "GAGTCCGAGCAGAAGAAGAAGGG"
    guides = crispr.find_guides(crispr.EXAMPLE_EMX1_SPCAS9, "SpCas9")
    assert guides


def test_crispr_placeholders_are_not_emx1() -> None:
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "placeholder=crispr.EXAMPLE_EMX1_SPCAS9" not in text
    assert 'placeholder="EMX1"' not in text
    assert "Paste target DNA" in text
    assert "Optional name" in text


def test_crispr_fresh_session_target_and_name_are_empty() -> None:
    app = open_module("crispr", timeout=90)
    assert not app.exception
    target = app.text_area(key="crispr_text").value
    name = app.text_input(key="crispr_name").value
    assert not str(target or "").strip()
    assert not str(name or "").strip()
    assert crispr.EXAMPLE_EMX1_SPCAS9 not in str(target or "")
    assert str(name or "").upper() != "EMX1"


def test_crispr_user_input_survives_module_switch() -> None:
    app = open_module("crispr", timeout=90)
    assert not app.exception
    pasted = "ATGCTAGTCGGATCCTGAATGCGTACGACTAG"
    app.text_area(key="crispr_text").set_value(pasted).run()
    app.text_input(key="crispr_name").set_value("user-target").run()
    assert not app.exception
    switch_module(app, "dna")
    assert not app.exception
    switch_module(app, "crispr")
    assert not app.exception
    assert app.text_area(key="crispr_text").value == pasted
    assert app.text_input(key="crispr_name").value == "user-target"
