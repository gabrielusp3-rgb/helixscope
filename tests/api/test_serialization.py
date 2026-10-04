"""Strict JSON, DataFrame, and array transport."""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
from fastapi.testclient import TestClient

from helixscope_core.serialize import arrays_equivalent
from helixscope_api.serialization import assert_strict_json, encode_array, to_transport


def test_to_transport_nan_inf_and_zero() -> None:
    assert to_transport(0.0) == 0.0
    assert to_transport(float("nan")) == {"value": None, "value_state": "NAN"}
    assert to_transport(float("inf")) == {"value": None, "value_state": "INF"}
    assert to_transport(float("-inf")) == {"value": None, "value_state": "-INF"}
    assert to_transport(None) is None
    payload = to_transport({"gc": float("nan"), "zero": 0.0})
    assert_strict_json(payload)
    assert payload["zero"] == 0.0
    assert payload["gc"]["value_state"] == "NAN"


def test_dataframe_and_array_transport_do_not_mutate() -> None:
    frame = pd.DataFrame({"value": [0.0, float("nan")], "label": ["a", "b"]})
    original = frame.copy(deep=True)
    encoded = to_transport(frame)
    pd.testing.assert_frame_equal(frame, original)
    assert encoded["columns"] == ["value", "label"]
    assert encoded["records"][0]["value"] == 0.0
    assert encoded["records"][1]["value"]["value_state"] == "NAN"
    array = np.array([[1.0, float("nan")], [0.0, 2.0]])
    snapshot = array.copy()
    packed = encode_array(array)
    assert arrays_equivalent(array, snapshot)
    assert packed["shape"] == [2, 2]
    assert packed["data"][0][1]["value_state"] == "NAN"
    assert packed["data"][1][0] == 0.0
    assert_strict_json(packed)


def test_api_responses_are_strict_json(client: TestClient) -> None:
    responses = [
        client.get("/health/live"),
        client.get("/api/v1"),
        client.get("/api/v1/system/capabilities"),
        client.get("/api/v1/references"),
        client.post("/api/v1/dna/analyze", json={"sequence": "NNNN", "include_explanation": False}),
        client.post("/api/v1/protein/analyze", json={"sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", "include_explanation": False}),
    ]
    for response in responses:
        assert response.status_code == 200, response.text
        body = response.json()
        json.dumps(body, allow_nan=False)
        raw = response.content.decode("utf-8")
        assert "Infinity" not in raw
