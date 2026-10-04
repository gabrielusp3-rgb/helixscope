"""Testes de folding de RNA: parser, pares, proveniencia e deteccao de ferramenta.

Nenhum teste inventa uma estrutura MFE e a atribui ao ViennaRNA. Fixtures de
dot-bracket testam o parser. Folding ao vivo so corre se RNA/RNAfold existirem.
A sequencia hsa-let-7a-5p e o registro publicado MIMAT0000062 (miRBase).
"""

from __future__ import annotations

import inspect
import subprocess
from pathlib import Path

import pytest

from modules import provenance, rna_analysis, rna_folding, scientific_checks
from tests.helix_apptest import open_module

FIXTURES = Path(__file__).parent / "fixtures"
LET7A = "UGAGGUAGUAGGUUGUAUAGUU"
HAIRPIN_SEQ = "GGGAAACCC"
HAIRPIN_DB = "(((...)))"
VIENNA_DOC_SEQ = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"


def _hairpin_engine(sequence: str) -> tuple[str, float]:
    if sequence != HAIRPIN_SEQ:
        raise rna_folding.FoldingError("unexpected sequence", "TOOL_FAILED")
    return HAIRPIN_DB, -3.3


def test_availability_does_not_fake_viennarna():
    info = rna_folding.tool_availability()
    delegated = rna_analysis.secondary_structure_availability()
    assert delegated["available"] is info["available"]
    assert "ViennaRNA" in info["reason"] or "RNAfold" in info["reason"]
    if not info["available"]:
        assert info["python_rna"]["available"] is False
        assert info["rnafold"]["available"] is False
        assert info["python_rna"]["version"] == ""
        assert info["rnafold"]["version"] == ""
        assert "unavailable" in info["reason"].lower() or "not detected" in info["reason"].lower()


def test_python_bindings_absent_are_not_faked(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "RNA":
            raise ImportError("blocked")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    detected = rna_folding.detect_python_bindings()
    assert detected["available"] is False
    assert detected["version"] == ""
    assert detected["version_status"] == "unavailable"


def test_rnafold_absent_is_not_faked(monkeypatch):
    monkeypatch.setattr(rna_folding, "_resolve_allowlisted_executable", lambda _names: None)
    detected = rna_folding.detect_rnafold_executable()
    assert detected["available"] is False
    assert detected["version"] == ""
    assert "not installed" in detected["reason"]


def test_dna_and_ambiguous_rna_are_rejected():
    with pytest.raises(rna_folding.FoldingError) as dna:
        rna_folding.validate_rna_for_folding("ACGTACGT")
    assert dna.value.category == "INVALID_INPUT"
    assert "T" in str(dna.value)
    with pytest.raises(rna_folding.FoldingError) as ambig:
        rna_folding.validate_rna_for_folding("ACGN")
    assert ambig.value.category == "INVALID_INPUT"
    with pytest.raises(rna_folding.FoldingError):
        rna_folding.validate_rna_for_folding("")
    one = rna_folding.validate_rna_for_folding("A")
    assert one["sequence"] == "A"
    info = rna_folding.validate_rna_for_folding("  acguacgu  ")
    assert info["sequence"] == "ACGUACGU"
    assert info["normalization"]["whitespace_removed"] is True
    assert info["normalization"]["uppercased"] is True
    assert info["normalization"]["nucleotides_changed"] is False
    assert info["hash"] == provenance.sequence_digest("ACGUACGU")


def test_resource_limit_is_not_empty_structure(monkeypatch):
    monkeypatch.setattr(rna_folding, "MAX_FOLD_NT", 4)
    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.validate_rna_for_folding("ACGUAC")
    assert exc.value.category == "RESOURCE_LIMIT"
    assert "empty" not in str(exc.value).lower() or "RESOURCE_LIMIT" in str(exc.value)


def test_parse_rnafold_output_fixture_is_format_not_live_mfe():
    body = (FIXTURES / "rnafold_output_hairpin.txt").read_text(encoding="utf-8")
    parsed = rna_folding.parse_rnafold_output(body, HAIRPIN_SEQ)
    assert parsed["structure"] == HAIRPIN_DB
    assert parsed["mfe_kcal_mol"] == pytest.approx(-3.30)
    with pytest.raises(rna_folding.FoldingError) as mismatch:
        rna_folding.parse_rnafold_output(body, "AAAAAAAAA")
    assert mismatch.value.category == "PARSING_ERROR"
    with pytest.raises(rna_folding.FoldingError) as empty:
        rna_folding.parse_rnafold_output("   ", HAIRPIN_SEQ)
    assert empty.value.category == "PARSING_ERROR"
    huge = "A" * (rna_folding.MAX_OUTPUT_BYTES + 1)
    with pytest.raises(rna_folding.FoldingError) as limited:
        rna_folding.parse_rnafold_output(huge, HAIRPIN_SEQ)
    assert limited.value.category == "RESOURCE_LIMIT"


def test_dot_bracket_parser_round_trip_and_balance():
    pairs = rna_folding.parse_dot_bracket(HAIRPIN_DB)
    assert pairs == [(0, 8), (1, 7), (2, 6)]
    rebuilt = rna_folding.reconstruct_dot_bracket(9, pairs)
    assert rebuilt == HAIRPIN_DB
    assert scientific_checks.rna_parentheses_are_balanced(HAIRPIN_DB)
    assert scientific_checks.rna_parentheses_are_balanced("((...)") is False
    with pytest.raises(rna_folding.FoldingError) as unbalanced:
        rna_folding.parse_dot_bracket("((...)")
    assert unbalanced.value.category == "PARSING_ERROR"
    with pytest.raises(rna_folding.FoldingError) as knots:
        rna_folding.parse_dot_bracket("([)]")
    assert knots.value.category == "PARSING_ERROR"


def test_impossible_pairs_are_rejected():
    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.validate_pairs_against_sequence("AAAAAA", [(0, 5)])
    assert exc.value.category == "PARSING_ERROR"
    rna_folding.validate_pairs_against_sequence(HAIRPIN_SEQ, [(0, 8), (1, 7), (2, 6)])


def test_fold_pipeline_with_injected_engine_validates_output():
    result = rna_folding.fold_rna(HAIRPIN_SEQ, source="unit test", fold_fn=_hairpin_engine)
    assert result["status"] == "PREDICTED"
    assert result["kind"] == "predicted"
    assert result["sequence"] == HAIRPIN_SEQ
    assert result["dot_bracket"] == HAIRPIN_DB
    assert result["mfe_kcal_mol"] == pytest.approx(-3.3)
    assert result["mfe_unit"] == "kcal/mol"
    assert result["sequence_hash"] == provenance.sequence_digest(HAIRPIN_SEQ)
    assert result["n_pairs"] == 3
    assert result["base_pairs"][0]["position_0based"] == 0
    assert result["base_pairs"][0]["base"] == "G"
    assert result["base_pairs"][0]["paired_position_0based"] == 8
    assert result["base_pairs"][0]["paired_base"] == "C"
    assert result["coordinate_system"] == "0-based"
    types = {item["type"] for item in result["elements"]}
    assert "stem" in types
    assert "hairpin_loop" in types
    assert "experimental" not in result["disclaimer"].lower() or "not an experimental" in result["disclaimer"].lower()
    assert result["cache_status"] == "live"
    rebuilt = rna_folding.reconstruct_dot_bracket(9, [(p["position_0based"], p["paired_position_0based"]) for p in result["base_pairs"]])
    assert rebuilt == result["dot_bracket"]


def test_injected_malformed_structure_is_parsing_error():
    def bad(_sequence: str) -> tuple[str, float]:
        return "((...", -1.0

    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=bad)
    assert exc.value.category == "PARSING_ERROR"


def test_injected_length_mismatch_is_parsing_error():
    def short(_sequence: str) -> tuple[str, float]:
        return "...", 0.0

    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=short)
    assert exc.value.category == "PARSING_ERROR"


def test_mfe_missing_or_non_finite_is_not_zero():
    def nan_energy(_sequence: str) -> tuple[str, float]:
        return ".........", float("nan")

    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.fold_rna("AAAAAAAAA", fold_fn=nan_energy)
    assert exc.value.category == "PARSING_ERROR"
    assert scientific_checks.mfe_kcal_mol_is_valid(0.0)
    assert scientific_checks.mfe_kcal_mol_is_valid(float("nan")) is False


def test_all_unpaired_zero_mfe_is_kept():
    def unpaired(sequence: str) -> tuple[str, float]:
        return "." * len(sequence), 0.0

    result = rna_folding.fold_rna("AAAAAAAAAA", fold_fn=unpaired)
    assert result["mfe_kcal_mol"] == 0.0
    assert result["n_pairs"] == 0
    csv_text = rna_folding.export_pairs_csv(result)
    assert "0.0" in csv_text or "0" in csv_text


def test_region_folding_records_subsequence_disclaimer():
    parent = "AAAA" + HAIRPIN_SEQ + "UUUU"
    result = rna_folding.fold_rna(
        parent,
        source="unit test",
        region_start=4,
        region_end=13,
        full_sequence=parent,
        fold_fn=_hairpin_engine,
    )
    assert result["region_mode"] == "selected subsequence"
    assert result["region_start"] == 4
    assert result["region_end"] == 13
    assert result["sequence"] == HAIRPIN_SEQ
    assert result["sequence_hash"] == provenance.sequence_digest(HAIRPIN_SEQ)
    assert result["full_sequence_hash"] == provenance.sequence_digest(parent)
    assert "subsequence" in result["disclaimer"].lower()
    assert "native local structure" in result["disclaimer"].lower()


def test_ncbi_sequence_source_stays_separate_from_predicted_structure():
    result = rna_folding.fold_rna(
        HAIRPIN_SEQ,
        source="NCBI Entrez",
        accession="NM_000000",
        version="NM_000000.1",
        organism="Homo sapiens",
        retrieved_at="2026-01-01T00:00:00Z",
        fold_fn=_hairpin_engine,
    )
    assert result["sequence_source"] == "NCBI Entrez"
    assert result["structure_source"]
    assert result["status"] == "PREDICTED"
    assert result["kind"] == "predicted"
    assert "NCBI experimental" not in result["disclaimer"]
    assert result["accession"] == "NM_000000"


def test_cache_key_changes_with_sequence_parameter_and_version():
    params = rna_folding.default_fold_parameters("full sequence")
    assert params["temperature_c"] == rna_folding.VIENNA_DEFAULT_TEMPERATURE_C
    assert params["temperature_passed_by_helixscope"] is True
    assert params["salt_passed_by_helixscope"] is False
    a = rna_folding.cache_key(sequence_hash="aaa", backend="viennarna_python", tool_version="", parameters=params)
    b = rna_folding.cache_key(sequence_hash="bbb", backend="viennarna_python", tool_version="", parameters=params)
    c = rna_folding.cache_key(
        sequence_hash="aaa",
        backend="viennarna_python",
        tool_version="",
        parameters=rna_folding.default_fold_parameters("selected subsequence"),
    )
    d = rna_folding.cache_key(sequence_hash="aaa", backend="viennarna_python", tool_version="2.7.2", parameters=params)
    assert a != b
    assert a != c
    assert a != d
    live = rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=_hairpin_engine)
    cached = rna_folding.mark_cached_result(live)
    assert cached["cache_status"] == "cached"
    assert cached["computed_at_utc"] == live["computed_at_utc"]
    assert cached["served_from_cache_at_utc"]


def test_tool_not_installed_is_not_a_structure(monkeypatch):
    monkeypatch.setattr(rna_folding, "detect_python_bindings", lambda: {"available": False, "version": ""})
    monkeypatch.setattr(rna_folding, "detect_rnafold_executable", lambda: {"available": False, "version": "", "path": ""})
    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.fold_rna(HAIRPIN_SEQ)
    assert exc.value.category == "TOOL_NOT_INSTALLED"


def test_timeout_is_timeout_not_empty_structure(monkeypatch, tmp_path):
    fake = tmp_path / "RNAfold.exe"
    fake.write_text("", encoding="utf-8")
    monkeypatch.setattr(
        rna_folding,
        "detect_rnafold_executable",
        lambda: {"available": True, "path": str(fake), "version": "RNAfold 0"},
    )

    def boom(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd="RNAfold", timeout=1)

    monkeypatch.setattr(rna_folding.subprocess, "run", boom)
    monkeypatch.setattr(rna_folding, "detect_python_bindings", lambda: {"available": False, "version": ""})
    with pytest.raises(rna_folding.FoldingError) as exc:
        rna_folding.fold_rna(HAIRPIN_SEQ, backend="rnafold_exe")
    assert exc.value.category == "TIMEOUT"


def test_subprocess_uses_nops_noconv_and_no_shell():
    source = inspect.getsource(rna_folding._fold_executable)
    assert "shell=False" in source
    assert "shell=True" not in source
    assert "--noPS" in source
    assert "--noconv" in source
    assert "--temp=" in source


def test_pair_selection_and_msa_column_map():
    result = rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=_hairpin_engine)
    pair = rna_folding.pair_at_position(result, 0)
    assert pair["paired_position_0based"] == 8
    assert rna_folding.pair_at_position(result, 4) is None
    span = rna_folding.spans_for_position(result, 0)
    assert span["start"] == 0
    assert span["end"] == 1
    assert span["status"] == "PREDICTED"
    assert rna_folding.map_msa_column_to_fold_position([None, 0, 1], 0) is None
    assert rna_folding.map_msa_column_to_fold_position([None, 0, 1], 2) == 1


def test_exports_include_provenance():
    result = rna_folding.fold_rna(HAIRPIN_SEQ, source="unit test", fold_fn=_hairpin_engine)
    dbn = rna_folding.export_dot_bracket(result)
    assert result["sequence_hash"] in dbn
    assert "PREDICTED" in dbn
    assert HAIRPIN_DB in dbn
    bundle = rna_folding.export_result_bundle(result)
    assert bundle["mfe_kcal_mol"] == pytest.approx(-3.3)
    assert bundle["status"] == "PREDICTED"


def test_let7a_fixture_is_published_sequence_not_an_expected_structure():
    text = (FIXTURES / "rna_hsa_let7a.fa").read_text(encoding="utf-8")
    assert "MIMAT0000062" in text
    assert LET7A in text.replace("\n", "")
    validated = rna_folding.validate_rna_for_folding(LET7A)
    assert validated["sequence"] == LET7A
    assert validated["hash"] == provenance.sequence_digest(LET7A)


def test_reproducible_injected_fold():
    first = rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=_hairpin_engine)
    second = rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=_hairpin_engine)
    assert first["dot_bracket"] == second["dot_bracket"]
    assert first["mfe_kcal_mol"] == second["mfe_kcal_mol"]
    assert first["sequence_hash"] == second["sequence_hash"]
    assert first["method"] == second["method"]


@pytest.mark.skipif(
    not rna_folding.tool_availability()["available"],
    reason="ViennaRNA RNAlib/RNAfold not installed in this environment",
)
def test_live_viennarna_on_published_let7a_and_docs_example():
    let7 = rna_folding.fold_rna(LET7A, source="miRBase MIMAT0000062")
    assert let7["status"] == "PREDICTED"
    assert len(let7["dot_bracket"]) == len(LET7A)
    assert scientific_checks.mfe_kcal_mol_is_valid(let7["mfe_kcal_mol"])
    again = rna_folding.fold_rna(LET7A, source="miRBase MIMAT0000062")
    assert again["dot_bracket"] == let7["dot_bracket"]
    assert again["mfe_kcal_mol"] == let7["mfe_kcal_mol"]
    docs = rna_folding.fold_rna(VIENNA_DOC_SEQ, source="ViennaRNA Python examples")
    assert len(docs["dot_bracket"]) == len(VIENNA_DOC_SEQ)
    assert docs["kind"] == "predicted"


@pytest.mark.skipif(
    not rna_folding.detect_rnafold_executable().get("available"),
    reason="official RNAfold CLI not installed",
)
def test_live_rnafold_cli_on_published_let7a():
    result = rna_folding.fold_rna(
        LET7A,
        source="miRBase MIMAT0000062",
        backend="rnafold_exe",
    )
    assert result["status"] == "PREDICTED"
    assert result["backend"] == "rnafold_exe" or "rnafold" in str(result.get("backend") or "").lower()
    assert len(result["dot_bracket"]) == len(LET7A)
    assert scientific_checks.mfe_kcal_mol_is_valid(result["mfe_kcal_mol"])


def test_arc_and_circular_layouts_use_only_real_pairs():
    result = rna_folding.fold_rna(HAIRPIN_SEQ, fold_fn=_hairpin_engine)
    paths = rna_folding.arc_paths(9, result["base_pairs"])
    assert len(paths) == 3
    assert paths[0]["i"] == 0
    assert paths[0]["j"] == 8
    layout = rna_folding.circular_layout(HAIRPIN_SEQ, result["base_pairs"])
    assert layout["n"] == 9
    assert len(layout["chords"]) == 3


def test_all_u_and_invalid_character_are_classified():
    info = rna_folding.validate_rna_for_folding("U" * 20)
    assert info["sequence"] == "U" * 20
    with pytest.raises(rna_folding.FoldingError) as invalid:
        rna_folding.validate_rna_for_folding("ACGX")
    assert invalid.value.category == "INVALID_INPUT"


def test_json_export_keeps_zero_and_does_not_turn_nan_into_zero():
    result = rna_folding.fold_rna("AAAAAAAAAA", fold_fn=lambda sequence: ("." * len(sequence), 0.0))
    bundle = rna_folding.export_result_bundle(result)
    assert bundle["mfe_kcal_mol"] == 0.0
    assert bundle["n_pairs"] == 0
    sanitized = provenance.json_safe({"mfe_kcal_mol": float("nan"), "n_pairs": 0})
    assert sanitized["mfe_kcal_mol"] is None
    assert sanitized["n_pairs"] == 0
    csv_text = rna_folding.export_pairs_csv(result)
    assert "N/A" not in csv_text.splitlines()[0]
    assert "0.0" in csv_text or ",0," in csv_text or csv_text.endswith("0\n") or ",0\n" in csv_text


def test_allowlist_rejects_wrong_executable_basename(monkeypatch, tmp_path):
    fake = tmp_path / "malware.exe"
    fake.write_text("", encoding="utf-8")
    monkeypatch.setattr(rna_folding.shutil, "which", lambda _name: str(fake))
    resolved = rna_folding._resolve_allowlisted_executable(rna_folding.RNAFOLD_NAMES)
    assert resolved is None or "malware" not in resolved.lower()
    if resolved:
        assert rna_folding._normalized_tool_basename(resolved) == "rnafold"


def test_parser_and_layout_complete_within_one_second_at_cap():
    import time

    sequence = "A" * rna_folding.MAX_FOLD_NT
    structure = "." * rna_folding.MAX_FOLD_NT
    started = time.perf_counter()
    pairs = rna_folding.parse_dot_bracket(structure)
    rna_folding.arc_paths(len(sequence), [])
    rna_folding.circular_layout(sequence, [])
    elapsed = time.perf_counter() - started
    assert pairs == []
    assert elapsed < 1.0


def test_rna_tab_reaches_structure_before_cds_abort():
    source = Path(__file__).resolve().parents[1] / "app.py"
    text = source.read_text(encoding="utf-8")
    start = text.index("def render_rna_analysis")
    end = text.index("def render_protein_analysis")
    body = text[start:end]
    assert body.index("_render_rna_structure_section") < body.index("if resolution is None:")
    assert "fold_rna(" not in body


def test_short_rna_analyze_reaches_structure_toggle_without_fake_fold():
    app = open_module("rna", timeout=60)
    assert not app.exception
    app.text_area(key="rna_text").set_value(LET7A).run()
    assert not app.exception
    app.button(key="rna_analyze").click().run()
    assert not app.exception
    app.toggle(key="rna_structure").set_value(True).run()
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for item in list(app.info) + list(app.error) + list(app.markdown)
    )
    fold_keys = [item.key for item in app.button]
    assert "Codon-usage sections require a resolved coding region" in combined or any(
        "coding" in str(getattr(item, "value", "")).lower() for item in app.error
    )
    if not rna_folding.tool_availability()["available"]:
        assert "RNA folding unavailable in this environment." in combined
        assert "rna_fold_run" not in fold_keys
        assert "This is a computational prediction" in combined or "not an experimental" in combined.lower() or "will not draw an invented structure" in combined
    else:
        assert "rna_fold_run" in fold_keys
