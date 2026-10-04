"""Deteccao de executaveis: allowlist, sem varrer o disco, diagnostico de plataforma."""

from __future__ import annotations

import os

from modules import tool_detection


def test_resolve_rejects_wrong_basename(monkeypatch, tmp_path):
    fake = tmp_path / "malware.exe"
    fake.write_text("", encoding="utf-8")
    monkeypatch.setattr(tool_detection.shutil, "which", lambda _name: str(fake))
    monkeypatch.setattr(tool_detection, "VIENNARNA_WINDOWS_CANDIDATES", ())
    assert tool_detection.resolve_allowlisted_executable(("RNAfold", "RNAfold.exe")) is None


def test_resolve_accepts_env_file_with_allowlisted_basename(monkeypatch, tmp_path):
    exe = tmp_path / "RNAfold.exe"
    exe.write_text("", encoding="utf-8")
    monkeypatch.setenv("HELIXSCOPE_RNAFOLD", str(exe))
    monkeypatch.setattr(tool_detection.shutil, "which", lambda _name: None)
    found = tool_detection.resolve_allowlisted_executable(
        ("RNAfold", "RNAfold.exe"),
        env_var="HELIXSCOPE_RNAFOLD",
        extra_file_candidates=(),
    )
    assert found == os.path.abspath(str(exe))


def test_env_wrong_basename_is_not_used(monkeypatch, tmp_path):
    exe = tmp_path / "python.exe"
    exe.write_text("", encoding="utf-8")
    monkeypatch.setenv("HELIXSCOPE_RNAFOLD", str(exe))
    monkeypatch.setattr(tool_detection.shutil, "which", lambda _name: None)
    assert (
        tool_detection.resolve_allowlisted_executable(
            ("RNAfold",),
            env_var="HELIXSCOPE_RNAFOLD",
            extra_file_candidates=(),
        )
        is None
    )


def test_environment_report_does_not_embed_home_path():
    report = tool_detection.environment_report(
        tools={
            "RNAfold": {
                "available": False,
                "version": "",
                "path": os.path.join(os.path.expanduser("~"), "secret", "RNAfold.exe"),
            }
        }
    )
    assert report["python_version"]
    assert report["platform"]
    dumped = str(report)
    assert os.path.expanduser("~") not in dumped
    assert report["tools"]["RNAfold"]["status"] == tool_detection.TOOL_STATUS_NOT_INSTALLED
    assert report["tools"]["RNAfold"]["path"] == "RNAfold.exe"


def test_classify_installed_without_version():
    status = tool_detection.classify_tool_record(
        available=True, path="C:\\RNAfold.exe", version=""
    )
    assert status == tool_detection.TOOL_STATUS_VERSION_UNAVAILABLE
