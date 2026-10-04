"""Propriedades baratas de DNA e o escopo de produto dos engines ausentes."""

from __future__ import annotations

from pathlib import Path

from modules import dna_analysis

ROOT = Path(__file__).resolve().parents[1]
_ALPHABET = "ACGT"


def _dna_from_index(index: int, length: int) -> str:
    chars = []
    value = index
    for _ in range(length):
        chars.append(_ALPHABET[value & 3])
        value >>= 2
    return "".join(chars)


def test_reverse_complement_involution_on_unique_short_dna() -> None:
    """20_000 sequencias unicas de 8 nt. A enumeracao de 4^10 fica fora da suíte."""
    for index in range(20_000):
        sequence = _dna_from_index(index, 8)
        complemented = dna_analysis.reverse_complement(sequence)
        assert len(complemented) == len(sequence)
        assert dna_analysis.reverse_complement(complemented) == sequence
        counts = {base: sequence.count(base) for base in _ALPHABET}
        assert sum(counts.values()) == len(sequence)
        gc = dna_analysis.gc_content(sequence)
        assert 0.0 <= gc <= 100.0


def test_cumulative_skew_final_matches_an_independent_walk() -> None:
    sequence = "GGCCAATTACGTN"  # N is rejected by canonical validation of some paths
    canonical = "GGCCAATTACGT"
    summary = dna_analysis.cumulative_skew_summary(canonical, plot_step=1)
    gc_total = 0
    at_total = 0
    for base in canonical:
        if base == "G":
            gc_total += 1
        elif base == "C":
            gc_total -= 1
        elif base == "A":
            at_total += 1
        elif base == "T":
            at_total -= 1
    assert summary["final_gc"] == gc_total
    assert summary["final_at"] == at_total
    assert summary["analyzed_length"] == len(canonical)
    assert summary["visualized_length"] == len(canonical)
    assert summary["input_length"] == summary["analyzed_length"]


def test_sidebar_css_scrolls_instead_of_clipping() -> None:
    styles = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    assert 'section[data-testid="stSidebar"] [data-testid="stSidebarContent"]' in styles
    content_rule = styles.split(
        'section[data-testid="stSidebar"] [data-testid="stSidebarContent"] {'
    )[1].split("}")[0]
    assert "overflow-y: auto" in content_rule
    assert "min-height: 0" in content_rule


def test_overview_does_not_print_unvalidated_engine_names() -> None:
    from tests.helix_apptest import open_module

    app = open_module("overview", timeout=90)
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption")
        for item in list(getattr(app, bucket, []) or [])
    )
    assert "STRIDE" not in combined
    assert "MSMS" not in combined
    assert "RNAfold" not in combined
    assert "NOT INSTALLED" not in combined.upper()
    assert "NOT_INSTALLED" not in combined.upper()
    assert ">DSSP<" not in combined
    from modules.phylogeny import detect_fasttree
    from modules.rna_folding import detect_python_bindings

    if detect_fasttree().get("available"):
        assert "FastTree" in combined
    if detect_python_bindings().get("available"):
        assert "ViennaRNA Python" in combined


def test_unvalidated_engines_stay_off_the_engines_screen() -> None:
    from ui.shell import engine_row_is_product_visible

    assert engine_row_is_product_visible(
        {"status": "not_installed", "available": False, "name": "DSSP"}
    ) is False
    assert engine_row_is_product_visible(
        {"status": "not_installed", "available": False, "name": "STRIDE"}
    ) is False
    assert engine_row_is_product_visible(
        {"status": "not_installed", "available": False, "name": "RNAfold"}
    ) is False
    assert engine_row_is_product_visible(
        {
            "status": "live_validated",
            "available": True,
            "live": {"ok": True, "status": "live_validated"},
        }
    ) is True
    assert engine_row_is_product_visible(
        {
            "status": "remote_validated",
            "available": True,
            "engine_location": "remote",
            "live": {"ok": True, "status": "remote_validated"},
        }
    ) is True


def test_missing_local_structure_engines_are_not_offered_as_buttons() -> None:
    source = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "Assign DSSP (mkdssp if installed)" not in source
    assert "Assign STRIDE (if installed)" not in source
    assert 'dssp.get("available") and st.button("Assign DSSP (local mkdssp)"' in source
    assert "Cas-OFFinder (if installed)" not in source
    assert "engine_choices.append(\"Cas-OFFinder\")" in source
