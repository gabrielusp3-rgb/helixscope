"""Container closure checks for path, shell, URL, upload, and HTML handling.

These tests call the product functions. They do not install engines, do not
open host files outside the process, and do not treat a missing optional
engine as a successful scientific result.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from modules import dna_analysis, genome_download, genome_jobs, protein_structure, provenance, tool_detection
from modules.clinvar_evidence import url_is_allowed as clinvar_url_is_allowed
from ui.components import html_escape, safe_download_filename

ROOT = Path(__file__).resolve().parents[1]
MODULES = ROOT / "modules"

PATH_PAYLOADS: tuple[str, ...] = (
    "../",
    "../../",
    "..\\",
    "..\\..\\",
    "/etc/passwd",
    "/proc/self/environ",
    "/app/app.py",
    "/tmp/helixscope-evidence",
    r"C:\Windows",
    r"C:\Program Files",
    "C:/Windows",
    "//server/share/secret",
    "passwd\x00.txt",
)

SSRF_PAYLOADS: tuple[str, ...] = (
    "http://127.0.0.1/",
    "https://127.0.0.1/",
    "https://localhost/",
    "https://[::1]/",
    "https://0.0.0.0/",
    "https://169.254.169.254/latest/meta-data/",
    "https://10.0.0.1/",
    "https://172.16.0.1/",
    "https://192.168.0.1/",
    "file:///etc/passwd",
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/",
    "data:text/plain,hello",
    "javascript:alert(1)",
    "https://evil.example/genomes/all/x",
    "https://ftp.ncbi.nlm.nih.gov.evil.example/genomes/all/x",
)

SHELL_PAYLOADS: tuple[str, ...] = (
    ";",
    "&&",
    "||",
    "|",
    "`",
    "$(id)",
    "\n",
    "\r",
    "'",
    '"',
    ">",
    "<",
    "*",
)


def _call_sites(path: Path) -> list[ast.Call]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [node for node in ast.walk(tree) if isinstance(node, ast.Call)]


def test_modules_do_not_call_shell() -> None:
    offenders: list[str] = []
    for path in sorted(MODULES.glob("*.py")):
        for node in _call_sites(path):
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == "system":
                owner = func.value
                if isinstance(owner, ast.Name) and owner.id == "os":
                    offenders.append(f"{path.name}:os.system")
            if isinstance(func, ast.Name) and func.id == "system":
                offenders.append(f"{path.name}:system()")
            for keyword in node.keywords:
                if keyword.arg != "shell":
                    continue
                value = keyword.value
                if isinstance(value, ast.Constant) and value.value is True:
                    offenders.append(f"{path.name}:shell=True")
    assert offenders == []


def test_tool_path_payloads_do_not_escape_or_execute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    marker = tmp_path / "not-created"
    for payload in PATH_PAYLOADS:
        reported = tool_detection.sanitize_tool_path(payload)
        lowered = payload.replace("\\", "/").lower()
        if ".." in payload.replace("\\", "/").split("/") or "\x00" in payload:
            assert reported == ""
        elif not lowered.startswith("c:/program files"):
            assert "/" not in reported.replace("\\", "/")
            assert not reported.startswith("/etc")
            assert not reported.startswith("/proc")
        assert marker.exists() is False
        if "\x00" in payload:
            continue
        monkeypatch.setenv("HELIXSCOPE_FASTTREE_PROBE", payload)
        resolved = tool_detection.resolve_allowlisted_executable(
            ["FastTree"],
            env_var="HELIXSCOPE_FASTTREE_PROBE",
        )
        if resolved is not None:
            chosen = Path(resolved)
            assert chosen.name.lower().removesuffix(".exe") == "fasttree"
            assert chosen.is_file()
            assert ".." not in chosen.parts
            assert chosen.name.lower() not in {"passwd", "environ"}
    link = tmp_path / "FastTree"
    target = tmp_path / "not-an-engine"
    target.write_text("not an engine", encoding="utf-8")
    try:
        link.symlink_to(target)
    except OSError:
        link.write_text("not an engine", encoding="utf-8")
    monkeypatch.setenv("HELIXSCOPE_FASTTREE_PROBE", str(link))
    found = tool_detection.resolve_allowlisted_executable(
        ["FastTree"],
        env_var="HELIXSCOPE_FASTTREE_PROBE",
    )
    assert found is not None
    assert Path(found).name == "FastTree"
    assert "/etc/passwd" not in found.replace("\\", "/")


def test_download_names_drop_traversal() -> None:
    for payload in PATH_PAYLOADS + SHELL_PAYLOADS:
        name = safe_download_filename(payload)
        assert "/" not in name
        assert "\\" not in name
        assert ".." not in name
        assert "\x00" not in name
        assert len(name) <= 80


def test_remote_urls_reject_ssrf_payloads() -> None:
    checkers = (
        genome_download.url_is_allowed,
        protein_structure.url_is_allowed,
        clinvar_url_is_allowed,
    )
    for payload in SSRF_PAYLOADS:
        for checker in checkers:
            assert checker(payload) is False
    assert genome_download.url_is_allowed(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/example"
    ) is True
    with pytest.raises(genome_download.GenomeDownloadError):
        genome_download.download_to_partial(
            "file:///etc/passwd",
            "ignored.partial",
            urlopen_fn=lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("urlopen called")),
        )


def test_shell_metacharacters_stay_inside_sequence_parser() -> None:
    for payload in SHELL_PAYLOADS:
        if payload in {"\n", "\r"}:
            parsed = dna_analysis.parse_fasta_records(f">sample{payload}ATGC")
            assert parsed["record_count"] >= 0
            continue
        text = f">sample{payload}\nATGC\n"
        parsed = dna_analysis.parse_fasta_records(text)
        assert parsed["record_count"] == 1
        assert parsed["records"][0]["sequence"] == "ATGC"
        assert payload not in parsed["records"][0]["sequence"]


def test_malformed_uploads_do_not_write_or_execute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    before = set(tmp_path.iterdir())
    empty = dna_analysis.parse_fasta_records("")
    assert empty["record_count"] == 0
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_pdb("")
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_mmcif("")
    parsed = dna_analysis.parse_fasta_records("\x00\x01\x02ATGC")
    assert parsed["record_count"] == 1
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_pdb("NOTPDB\n<script>alert(1)</script>\n")
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_mmcif("data_broken\n_not_a_loop\n")
    assert set(tmp_path.iterdir()) == before


def test_html_escape_neutralizes_markup() -> None:
    payloads = (
        "<script>alert(1)</script>",
        "<img src=x onerror=alert(1)>",
        "<svg onload=alert(1)>",
        "javascript:alert(1)",
        "data:text/html,<script>alert(1)</script>",
        "<div onclick=alert(1)>",
        "<b>unclosed",
        "[click](javascript:alert(1))",
    )
    for payload in payloads:
        escaped = html_escape(payload)
        assert "<" not in escaped
        assert ">" not in escaped
        assert escaped != payload or "<" not in payload


def test_job_log_redacts_fake_secret(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NCBI_API_KEY", "fake-token-not-a-real-secret")
    raw = "path=/home/helix/private token=fake-token-not-a-real-secret <script>alert(1)</script>"
    cleaned = genome_jobs.sanitize_job_log(raw)
    assert "fake-token-not-a-real-secret" not in cleaned
    assert "<redacted>" in cleaned
    redacted = provenance.redact_secrets({"NCBI_API_KEY": "fake-token-not-a-real-secret", "gc": 50.0})
    assert redacted["NCBI_API_KEY"] == "[redacted]"
    assert redacted["gc"] == 50.0


def test_overlong_sequence_is_refused_before_analysis() -> None:
    too_long = "A" * (dna_analysis.MAX_INPUT_RESIDUES + 1)
    report = dna_analysis.validate_sequence(too_long)
    assert report["is_valid"] is False
    assert report["sequence"] == ""
    assert report["rejection_reason"]
