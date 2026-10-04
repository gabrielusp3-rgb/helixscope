"""Prompt 2 transport: allowlisted bundled mmCIF for Mol*. No science rewrite."""

from __future__ import annotations

from fastapi.testclient import TestClient


def test_coordinate_file_serves_bundled_1crn(client: TestClient) -> None:
    response = client.post("/api/v1/structures/coordinate-file", json={"structure_id": "1CRN"})
    assert response.status_code == 200, response.text
    result = response.json()["result"]
    assert result["structure_id"] == "1CRN"
    assert result["status"] == "EXPERIMENTAL"
    assert result["format"] == "mmcif"
    assert isinstance(result.get("mmcif"), str) and result["mmcif"].startswith("data_")
    assert "http://" not in result["mmcif"][:80]


def test_coordinate_file_rejects_non_bundled_ids(client: TestClient) -> None:
    rejected = client.post("/api/v1/structures/coordinate-file", json={"structure_id": "4UN3"})
    assert rejected.status_code in {400, 422}
    url = client.post("/api/v1/structures/coordinate-file", json={"structure_id": "https://example.com/x.cif"})
    assert url.status_code in {400, 422}
