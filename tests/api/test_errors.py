"""Classified errors, not generic 500 for known conditions."""

from __future__ import annotations

from fastapi.testclient import TestClient


def test_invalid_alignment_mode(client: TestClient) -> None:
    response = client.post(
        "/api/v1/alignment/pairwise",
        json={"seq1": "ATGC", "seq2": "ATGC", "mode": "blast"},
    )
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "INVALID_INPUT"
    assert "traceback" not in response.text.lower()


def test_invalid_protein_is_classified(client: TestClient) -> None:
    response = client.post(
        "/api/v1/protein/analyze",
        json={"sequence": "not-a-peptide-000", "include_explanation": False},
    )
    assert response.status_code in {400, 422}
    body = response.json()
    assert "error" in body
    assert body["error"]["code"] in {"INVALID_INPUT", "PARSE_ERROR"}


def test_request_id_header(client: TestClient) -> None:
    response = client.get("/health/live")
    assert response.headers.get("x-request-id")
    custom = client.get("/health/live", headers={"X-Request-ID": "trace-test-1"})
    assert custom.headers.get("x-request-id") == "trace-test-1"
