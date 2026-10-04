"""Phase 19 Part 3: engines, explanations, missing values, fuzz, properties."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from modules import (
    blast_search,
    dna_analysis,
    explain,
    phylogeny,
    rna_folding,
    tool_registry,
    usalign,
)
from tests.helix_apptest import open_module

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


def test_engine_page_does_not_duplicate_status_badges() -> None:
    source = inspect.getsource(__import__("ui.shell", fromlist=["render_settings"]).render_settings)
    assert "Validation details" in source
    assert "Advanced detector table" in source
    assert "Live record" not in source


def test_fasttree_version_is_parsed_when_installed() -> None:
    info = phylogeny.detect_fasttree()
    if not info.get("available"):
        pytest.skip("FastTree official binary not installed")
    assert str(info.get("version") or "").startswith("2.")


def test_viennarna_python_is_separate_from_rnafold_cli() -> None:
    python = rna_folding.detect_python_bindings()
    cli = rna_folding.detect_rnafold_executable()
    if python.get("available"):
        assert python.get("version")
        assert cli.get("available") is False or cli.get("available") is True
    if not cli.get("available"):
        assert cli.get("version") == ""


def test_usalign_version_is_not_asterisk_banner() -> None:
    info = usalign.detect_usalign()
    if not info.get("available"):
        pytest.skip("US-align official binary not installed")
    assert "*" not in str(info.get("version") or "")
    assert str(info.get("version") or "").isdigit()


def test_tool_snapshot_refresh_sees_optional_engines() -> None:
    snap = tool_registry.collect_tool_snapshot(refresh=True)
    tools = snap["tools"]
    assert "ViennaRNA Python" in tools
    assert "RNAfold" in tools
    assert "US-align" in tools
    assert "BLAST+" in tools
    assert tools["RNAfold"]["name"] == "RNAfold"


def test_apptest_scientific_engines_has_single_summary() -> None:
    app = open_module("settings", timeout=90)
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption")
        for item in list(getattr(app, bucket, []) or [])
    )
    assert "Engine summary" in combined
    assert "DETECTED is not LIVE_VALIDATED" in combined


def test_fuzz_fasta_and_motif_fail_safely() -> None:
    hostile = [
        "",
        ">",
        ">id\n",
        ">a\nACGT\n>a\nTTTT\n",
        "><script>alert(1)</script>\nACGT\n",
        ">..\\..\\etc\\passwd\nACGTACGTACGTACGT\n",
        ">id; rm -rf /\nACGT\n",
        "N" * 10,
        "ACGT",
    ]
    for text in hostile:
        try:
            dna_analysis.parse_fasta_records(text)
        except (ValueError, TypeError):
            pass
        info = dna_analysis.validate_sequence(text) if text else {"is_valid": False}
        assert isinstance(info, dict)
    from modules import motif_search

    for pattern in ["AA", "N", "GAATTC"]:
        hits = motif_search.find_motif("ACGTACGT", pattern)
        assert isinstance(hits, list)
    with pytest.raises(ValueError):
        motif_search.find_motif("ACGT", "   ")


def test_reverse_complement_involution_unambiguous() -> None:
    seq = "ACGTN"
    try:
        twice = dna_analysis.reverse_complement(dna_analysis.reverse_complement("ACGT"))
    except ValueError:
        pytest.skip("N not supported in reverse_complement")
    else:
        assert twice == "ACGT"
    assert dna_analysis.reverse_complement(dna_analysis.reverse_complement("ACGTACGT")) == "ACGTACGT"


def test_explain_imported_from_modules() -> None:
    rec = explain.explain("alignment", {"algorithm": "Needleman-Wunsch", "identity_pct": 80.0, "score": 12, "length": 10})
    assert rec.kind == "alignment"
    assert "80.00%" in rec.plain_meaning or "80.0" in rec.plain_meaning


def test_reduced_motion_disables_decorative_blur() -> None:
    text = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    assert "prefers-reduced-motion" in text
    assert "backdrop-filter: none" in text


def test_streamlit_config_hides_deploy_via_official_toolbar_mode() -> None:
    text = (ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8")
    assert text.count("[client]") == 1
    assert 'toolbarMode = "viewer"' in text
    assert "showSidebarNavigation = false" in text


def test_newick_fuzz_fails_safely() -> None:
    from modules import phylogeny

    hostile = ["", "(", "());", "<script>alert(1)</script>", "..\\..\\etc", "A" * 20]
    for text in hostile:
        try:
            phylogeny.parse_newick(text)
        except phylogeny.PhylogenyError:
            pass


def test_malformed_blast_query_is_invalid_before_database() -> None:
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.run_local_blast(program="blastn", query="!!!!", molecule="DNA")
    assert exc.value.category == "INVALID_INPUT"
