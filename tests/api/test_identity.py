"""Runtime identity and OpenAPI hash must be comparable to the committed spec."""

from __future__ import annotations

import hashlib
import json
from fastapi.testclient import TestClient

from helixscope_api.contract import dump_openapi, openapi_sha256, runtime_identity
from helixscope_api.main import create_app


def test_identity_endpoint_is_sanitized(client: TestClient) -> None:
    response = client.get("/api/v1/system/identity")
    assert response.status_code == 200
    body = response.json()
    assert body["api_version"] == "v1"
    assert body["product_version"]
    assert body["core_version"]
    assert len(body["openapi_sha256"]) == 64
    assert "C:\\" not in json.dumps(body)
    assert "/Users/" not in json.dumps(body)
    dumped = json.dumps(body)
    assert "password" not in dumped.lower()
    assert "secret" not in dumped.lower()


def test_identity_hash_matches_live_openapi(client: TestClient) -> None:
    live = openapi_sha256(client.app)
    identity = client.get("/api/v1/system/identity").json()
    assert identity["openapi_sha256"] == live
    canonical = dump_openapi(client.app)
    assert hashlib.sha256(canonical.encode("utf-8")).hexdigest() == live


def test_runtime_identity_has_no_absolute_paths() -> None:
    app = create_app()
    payload = runtime_identity(app)
    text = json.dumps(payload)
    assert "C:\\\\Users" not in text
    assert "Desktop" not in text
