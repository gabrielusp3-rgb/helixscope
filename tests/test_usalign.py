"""US-align: deteccao, parser, subprocesso sem shell, live se instalado."""

from __future__ import annotations

import inspect
import subprocess
from pathlib import Path

import pytest

from modules import usalign

FIXTURES = Path(__file__).parent / "fixtures"


def test_usalign_not_installed_keeps_nucleic_unavailable(monkeypatch):
    monkeypatch.setattr(
        usalign.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: None,
    )
    detected = usalign.detect_usalign()
    assert detected["available"] is False
    assert detected["status"] == "not_installed" or "NOT_INSTALLED" in str(detected["status"]).upper() or detected["status"] in {"not_installed"}
    nucleic = usalign.nucleic_structural_comparison_availability()
    assert nucleic["available"] is False
    assert nucleic["status"] == "UNAVAILABLE"
    assert "C-alpha" in nucleic["rcsb_alignment_api"]


@pytest.mark.parametrize(
    "chain_id",
    [
        " ",
        " A",
        "A ",
        "A B",
        "\"A\"",
        "A;rm",
        "A&B",
        "A|B",
        "$A",
        "`A`",
        "A\n",
        "A\r",
        "-A",
        "--chain",
        "../A",
        "..\\A",
        "%2e%2e",
        "Ａ",
        "A" * 5,
        "A/B",
        "A\x00",
    ],
)
def test_chain_id_rejects_argument_injection(chain_id: str) -> None:
    with pytest.raises(usalign.USAlignError) as exc:
        usalign.normalize_chain_id(chain_id)
    assert exc.value.category == "INVALID_INPUT"


@pytest.mark.parametrize("chain_id", ["A", "a", "AA", "A1", "1", "ABCD"])
def test_chain_id_keeps_structure_token(chain_id: str) -> None:
    assert usalign.normalize_chain_id(chain_id) == chain_id


def test_structure_filename_rejects_path_syntax(tmp_path) -> None:
    with pytest.raises(usalign.USAlignError) as exc:
        usalign.write_structure_text("data_x\n", str(tmp_path), "foo..bar.cif")
    assert exc.value.category == "INVALID_INPUT"
    written = usalign.write_structure_text("data_x\n", str(tmp_path), "../secret.cif")
    assert Path(written).name == "secret.cif"
    assert Path(written).resolve().parent == tmp_path.resolve()
    kept = usalign.write_structure_text("data_x\n", str(tmp_path), "1CRN.cif")
    assert Path(kept).is_file()


def test_single_chain_write_keeps_only_that_chain(tmp_path) -> None:
    dest = tmp_path / "onlyA.cif"
    usalign._write_single_chain(str(FIXTURES / "1CRN.cif"), str(dest), "A")
    assert dest.is_file()
    assert dest.stat().st_size > 0
    with pytest.raises(usalign.USAlignError) as exc:
        usalign._write_single_chain(str(FIXTURES / "1CRN.cif"), str(tmp_path / "none.cif"), "Z")
    assert exc.value.category == "INVALID_INPUT"


def test_empty_chain_id_means_no_filter() -> None:
    assert usalign.normalize_chain_id("") == ""


def test_bad_chain_never_reaches_subprocess(monkeypatch) -> None:
    def fail_run(*_args, **_kwargs):
        raise AssertionError("subprocess must not start")

    monkeypatch.setattr(usalign.subprocess, "run", fail_run)
    with pytest.raises(usalign.USAlignError) as exc:
        usalign.align_structure_files("missing.cif", "missing.cif", chain1="-mol")
    assert exc.value.category == "INVALID_INPUT"


def test_usalign_subprocess_is_shell_false():
    probe_src = inspect.getsource(usalign._probe_version)
    align_src = inspect.getsource(usalign.align_structure_files)
    assert "shell=False" in probe_src
    assert "shell=True" not in probe_src
    assert "shell=False" in align_src
    assert "shell=True" not in align_src


def test_parse_usalign_1crn_self_fixture():
    text = (FIXTURES / "usalign_1crn_self.txt").read_text(encoding="utf-8")
    parsed = usalign.parse_usalign_stdout(text)
    assert parsed["n_aligned"] == 46
    assert parsed["rmsd_angstrom"] == pytest.approx(0.0)
    assert parsed["sequence_identity"] == pytest.approx(1.0)
    assert len(parsed["tm_scores"]) == 2
    assert parsed["tm_scores"][0]["value"] == pytest.approx(1.0)
    assert parsed["aligned_seq1"].startswith("TTCCPS")
    record = usalign.build_alignment_record(
        parsed,
        path1="1CRN.cif",
        path2="1CRN.cif",
        mol="prot",
        chain1="A",
        chain2="A",
        matrix=None,
        executable_version=str(parsed.get("version") or ""),
    )
    assert record["engine"] == "US-align"
    assert record["engine_location"] == "local"
    assert record["n_aligned_residue_pairs"] == 46
    assert record["rmsd_global_angstrom"] == record["rmsd_block0_angstrom"]
    assert record["superposition_visual"] == "UNAVAILABLE"
    assert "RCSB" not in record["engine"]


def test_parse_usalign_matrix_identity():
    text = (
        "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n"
        "m               t[m]        u[m][0]        u[m][1]        u[m][2]\n"
        "0       0.0000000000   1.0000000000  -0.0000000000   0.0000000000\n"
        "1       0.0000000000   0.0000000000   1.0000000000  -0.0000000000\n"
        "2      -0.0000000000  -0.0000000000   0.0000000000   1.0000000000\n"
    )
    matrix = usalign.parse_usalign_matrix(text)
    assert matrix is not None
    assert matrix[0] == pytest.approx(1.0)
    assert matrix[5] == pytest.approx(1.0)
    assert matrix[10] == pytest.approx(1.0)
    assert matrix[15] == pytest.approx(1.0)


def test_usalign_missing_output_is_parsing_error():
    with pytest.raises(usalign.USAlignError) as exc:
        usalign.parse_usalign_stdout("not an alignment")
    assert exc.value.category == "PARSING_ERROR"


def test_live_usalign_protein_self_and_nucleic_pair():
    info = usalign.detect_usalign()
    if not info.get("available"):
        pytest.skip("US-align official binary not installed")
    assert str(info.get("version") or "").isdigit()
    crn = FIXTURES / "1CRN.cif"
    bna = FIXTURES / "1BNA.cif"
    rna = FIXTURES / "1RNA.cif"
    self_aln = usalign.align_structure_files(
        str(crn), str(crn), mol="prot", record_validation=True
    )
    assert self_aln["n_aligned_residue_pairs"] == 46
    assert self_aln["rmsd_global_angstrom"] == pytest.approx(0.0, abs=1e-6)
    tms = [row["value"] for row in self_aln["tm_scores"]]
    assert tms
    assert all(value == pytest.approx(1.0, abs=1e-4) for value in tms)
    nucleic = usalign.align_structure_files(
        str(bna), str(rna), mol="RNA", record_validation=True
    )
    assert int(nucleic["n_aligned_residue_pairs"] or 0) >= 1
    assert nucleic["rmsd_global_angstrom"] is not None
    avail = usalign.nucleic_structural_comparison_availability()
    assert avail["status"] == "LIVE_VALIDATED"


def test_live_usalign_timeout_is_timeout(monkeypatch):
    info = usalign.detect_usalign()
    if not info.get("available"):
        pytest.skip("US-align official binary not installed")

    def boom(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd="USalign", timeout=1)

    monkeypatch.setattr(usalign.subprocess, "run", boom)
    crn = str(FIXTURES / "1CRN.cif")
    with pytest.raises(usalign.USAlignError) as exc:
        usalign.align_structure_files(crn, crn, timeout_s=1)
    assert exc.value.category == "TIMEOUT"
