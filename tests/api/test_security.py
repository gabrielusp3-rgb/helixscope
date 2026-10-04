"""Security bounds: no Streamlit/modules bypass, SSRF, path, oversized, XSS-as-data."""

from __future__ import annotations

import ast
from pathlib import Path

from fastapi.testclient import TestClient

API_ROOT = Path(__file__).resolve().parents[2] / "services" / "api" / "helixscope_api"


def test_api_source_has_no_forbidden_imports() -> None:
    offenders: list[str] = []
    for path in sorted(API_ROOT.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            for name in names:
                if name == "streamlit" or name.startswith("streamlit."):
                    offenders.append(f"{path.name}: streamlit")
                if name == "app" or name.startswith("app."):
                    offenders.append(f"{path.name}: app.py")
                if name == "modules" or name.startswith("modules."):
                    offenders.append(f"{path.name}: modules")
    assert offenders == []


def test_core_does_not_import_api() -> None:
    core = Path(__file__).resolve().parents[2] / "helixscope_core"
    offenders: list[str] = []
    for path in sorted(core.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        if "helixscope_api" in text:
            offenders.append(str(path))
    assert offenders == []


def test_malicious_inputs_are_data_or_classified(client: TestClient) -> None:
    script = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": "<script>alert(1)</script>", "include_explanation": False},
    )
    assert script.status_code in {200, 400, 422}
    body = script.json()
    if script.status_code == 200:
        assert body["result"]["status"] == "ERROR"
    html = client.post(
        "/api/v1/motif/search",
        json={"sequence": "ATGC", "pattern": "<b>not-a-motif", "molecule": "DNA"},
    )
    assert html.status_code in {200, 400, 422}
    traversal = client.post(
        "/api/v1/structures/fixture",
        json={"structure_id": "../etc/passwd"},
    )
    assert traversal.status_code == 422
    ssrf = client.post(
        "/api/v1/ncbi/fetch",
        json={"accession": "https://evil.example/seq", "email": "dev@example.com"},
    )
    assert ssrf.status_code == 400
    assert ssrf.json()["error"]["code"] == "INVALID_INPUT"
    exe = client.post(
        "/api/v1/msa/jobs",
        json={"fasta": ">a\nATGC\n>b\nATGG\n", "email": "dev@example.com", "backend": "C:\\\\evil.exe"},
    )
    assert exe.status_code == 400
    oversized = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": "A" * 250_001},
    )
    assert oversized.status_code == 422
    variant = client.post(
        "/api/v1/variants/identify",
        json={"text": "not a variant ;;; rm -rf /", "assembly": "GRCh38.p14"},
    )
    assert variant.status_code in {200, 400, 422}
    if variant.status_code >= 400:
        assert "traceback" not in variant.text.lower()
        assert "C:\\\\Users" not in variant.text
    unicode_ok = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": "ATGCΑ", "include_explanation": False},
    )
    assert unicode_ok.status_code in {200, 400, 422}


def test_casoffinder_does_not_queue_without_engine(client: TestClient) -> None:
    response = client.post(
        "/api/v1/crispr/casoffinder/jobs",
        json={"guide_sequence": "ACGTACGTACGTACGTACGT", "assembly_id": "GRCh38", "pam": "NGG"},
    )
    assert response.status_code in {409, 503}
    code = response.json()["error"]["code"]
    assert code in {"ENGINE_NOT_INSTALLED", "RESOURCE_LIMIT"}
