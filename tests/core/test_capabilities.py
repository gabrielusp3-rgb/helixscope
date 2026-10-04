"""Capability registry and NCBI client construction (no live HTTP)."""

from __future__ import annotations

import json

import pytest

from helixscope_core.capabilities import get_capabilities
from helixscope_core.references import list_references
from helixscope_core.remote import NCBIClient


def test_capabilities_redact_paths_and_omit_streamlit() -> None:
    snap = get_capabilities()
    dumped = json.dumps(snap, default=str)
    assert "C:\\" not in dumped
    assert "/Users/" not in dumped
    assert "streamlit" not in dumped.lower()
    assert "tools" in snap
    assert "opencl" in snap
    assert "ready" in snap["opencl"]
    for row in snap["tools"].values():
        binary = str(row.get("binary") or "")
        assert "\\" not in binary
        assert "/" not in binary


def test_list_references_has_no_filesystem_path_keys() -> None:
    rows = list_references()
    for row in rows:
        assert "fasta_path" not in row
        assert "path" not in row
        if row.get("not_a_public_assembly"):
            assert "TEST" in str(row.get("assembly_id") or "").upper() or "test" in str(
                row.get("label") or ""
            ).lower()


def test_ncbi_client_requires_email() -> None:
    with pytest.raises(ValueError):
        NCBIClient(email="")
    client = NCBIClient(email="helixscope-tests@example.invalid")
    assert client.email == "helixscope-tests@example.invalid"
    assert client.api_key is None
