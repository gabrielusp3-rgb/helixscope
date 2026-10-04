"""API transport for restored Streamlit-era scientific arrays."""

from __future__ import annotations

from fastapi.testclient import TestClient


def test_api_dna_exposes_dinucleotides_and_feature_tracks(client: TestClient) -> None:
    response = client.post(
        "/api/v1/dna/analyze",
        json={
            "sequence": "ATGC" * 40,
            "include_orfs": True,
            "include_restriction": True,
            "include_cpg": True,
            "include_kmers": True,
            "include_profiles": True,
            "include_santalucia": True,
            "include_explanation": True,
            "profile_window": 100,
            "profile_step": 50,
            "kmer_k": 3,
        },
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["status"] == "COMPUTED"
    assert "CG" in body["dinucleotides"]
    assert isinstance(body["restriction_enzyme_rows"], list)
    assert "kmer_summary" in body
    assert body.get("santalucia_tm")
    assert "sequence_feature_tracks" in body or body.get("sequence_features") is not None


def test_api_experimental_complex_rejects_non_catalog_id(client: TestClient) -> None:
    response = client.post(
        "/api/v1/structures/experimental-complex",
        json={"structure_id": "1CRN"},
    )
    assert response.status_code == 422


def test_api_dna_helix_3d_accepts_lod_window(client: TestClient) -> None:
    response = client.post(
        "/api/v1/dna/helix-3d",
        json={
            "sequence": "ATGC" * 20,
            "prefer": "illustrative",
            "lod": "medium",
            "region_start": 0,
            "region_end": 40,
            "inspect_position": 2,
            "center_on_selected": False,
        },
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["visual_lod"] == "medium"
    assert body["inspect_position"] == 2
    assert "position_mapping" in body
