"""Property-based API vs core DNA parity."""

from __future__ import annotations

import pytest

from helixscope_core import analyze_dna
from helixscope_api.main import create_app
from helixscope_api.serialization import to_transport
from fastapi.testclient import TestClient

hypothesis = pytest.importorskip("hypothesis")
st = pytest.importorskip("hypothesis.strategies")


@hypothesis.given(seq=st.text(alphabet="ACGT", min_size=12, max_size=48))
@hypothesis.settings(max_examples=12, deadline=None)
def test_random_dna_api_matches_core(seq: str) -> None:
    with TestClient(create_app()) as client:
        core = to_transport(analyze_dna(seq))
        response = client.post(
            "/api/v1/dna/analyze",
            json={"sequence": seq, "include_explanation": False},
        )
        assert response.status_code == 200
        api = response.json()["result"]
        assert api["status"] == core["status"]
        assert api["gc_content"] == core["gc_content"]
        assert api["length"] == core["length"]
        assert api["sequence_hash"] == core["sequence_hash"]
