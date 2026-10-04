"""Health, capabilities, references, CORS."""

from __future__ import annotations

import json

from fastapi.testclient import TestClient


def test_live_and_ready(client: TestClient) -> None:
    live = client.get("/health/live")
    assert live.status_code == 200
    assert live.json()["status"] == "alive"
    ready = client.get("/health/ready")
    assert ready.status_code == 200
    body = ready.json()
    assert body["status"] == "ready"
    assert body["product_version"] == "0.24.3-19"
    assert body["api_version"] == "v1"
    assert body["job_runner"] is True
    json.dumps(body, allow_nan=False)


def test_api_root(client: TestClient) -> None:
    response = client.get("/api/v1")
    assert response.status_code == 200
    body = response.json()
    assert body["name"] == "HelixScope API"
    assert body["product_version"] == "0.24.3-19"


def test_capabilities_have_no_home_paths(client: TestClient) -> None:
    response = client.get("/api/v1/system/capabilities")
    assert response.status_code == 200
    payload = response.json()
    json.dumps(payload, allow_nan=False)
    dumped = json.dumps(payload)
    assert "C:\\\\Users" not in dumped
    assert "C:/Users" not in dumped
    tools = payload["result"]["tools"]
    assert isinstance(tools, dict)
    opencl = payload["result"]["opencl"]
    assert "ready" in opencl


def test_references_preserve_status_vocabulary(client: TestClient) -> None:
    response = client.get("/api/v1/references")
    assert response.status_code == 200
    rows = response.json()["result"]["references"]
    assert isinstance(rows, list)
    for row in rows:
        assert "assembly_id" in row
        assert "status" in row
        assert "ready" in row
        if row.get("not_a_public_assembly"):
            assert "GRCh38" not in str(row.get("assembly_id") or "")


def test_cors_allows_localhost_3000_not_star(client: TestClient) -> None:
    allowed = client.get(
        "/api/v1/system/capabilities",
        headers={"Origin": "http://localhost:3000"},
    )
    assert allowed.headers.get("access-control-allow-origin") == "http://localhost:3000"
    denied = client.get(
        "/api/v1/system/capabilities",
        headers={"Origin": "https://evil.example"},
    )
    assert denied.headers.get("access-control-allow-origin") != "*"
    assert denied.headers.get("access-control-allow-origin") != "https://evil.example"
    assert "access-control-allow-credentials" not in {k.lower() for k in denied.headers.keys()} or (
        denied.headers.get("access-control-allow-credentials") != "true"
    )
