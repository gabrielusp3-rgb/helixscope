"""Compare remote jobs restore Core overlay atoms without inventing XYZ."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from helixscope_api.job_handlers import run_compare_remote_job
from helixscope_api.schemas import CompareRemoteJobRequest


def test_compare_remote_schema_accepts_alphafold_id() -> None:
    body = CompareRemoteJobRequest(
        reference_entry="AF-P01308-F1",
        target_entry="8HSF",
        fetch_coordinates=True,
    )
    dumped = body.model_dump()
    assert dumped["fetch_coordinates"] is True
    assert dumped["reference_entry"] == "AF-P01308-F1"


def test_compare_remote_schema_rejects_url_shaped_id() -> None:
    with pytest.raises(Exception):
        CompareRemoteJobRequest(
            reference_entry="https://evil.example/1crn",
            target_entry="8HSF",
        )


def test_compare_job_passes_loaded_coordinates(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}

    def fake_load(entry_id: str, **kwargs):
        return {
            "status": "AVAILABLE",
            "parsed": {"atoms": [{"x": 1.0, "y": 2.0, "z": 3.0, "name": "CA"}]},
            "kind_label": "EXPERIMENTAL",
            "entry_id": entry_id,
        }

    def fake_compare(**kwargs):
        captured.update(kwargs)
        return {
            "mode": "structures",
            "alignment": {"rmsd_global_angstrom": 1.25, "n_aligned_residue_pairs": 10},
            "superposition": {
                "reference_atoms": [{"x": 1.0, "y": 2.0, "z": 3.0}],
                "transformed_target_atoms": [{"x": 4.0, "y": 5.0, "z": 6.0}],
            },
        }

    monkeypatch.setattr(
        "helixscope_api.job_handlers.load_parsed_structure_for_entry", fake_load
    )
    monkeypatch.setattr(
        "helixscope_api.job_handlers.compare_structures_remote", fake_compare
    )
    record = SimpleNamespace(job_id="test")
    out = run_compare_remote_job(
        {
            "reference_entry": "8HSK",
            "target_entry": "8HSF",
            "reference_chain": "A",
            "target_chain": "A",
            "method": "tm-align",
            "fetch_coordinates": True,
        },
        record,  # type: ignore[arg-type]
    )
    assert captured["reference_parsed"]["atoms"][0]["x"] == 1.0
    assert captured["target_parsed"]["atoms"][0]["z"] == 3.0
    assert out["coordinate_fetch"]["overlay_status"] == "AVAILABLE"
    assert out["superposition"]["transformed_target_atoms"][0]["x"] == 4.0


def test_compare_job_skips_coordinate_load_when_disabled(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_load(*args, **kwargs):
        raise AssertionError("coordinate fetch was not requested")

    captured: dict = {}

    def fake_compare(**kwargs):
        captured.update(kwargs)
        return {"mode": "structures", "alignment": {}, "superposition": None}

    monkeypatch.setattr(
        "helixscope_api.job_handlers.load_parsed_structure_for_entry", fail_load
    )
    monkeypatch.setattr(
        "helixscope_api.job_handlers.compare_structures_remote", fake_compare
    )
    out = run_compare_remote_job(
        {"reference_entry": "8HSK", "target_entry": "8HSF", "fetch_coordinates": False},
        SimpleNamespace(job_id="t"),  # type: ignore[arg-type]
    )
    assert captured["reference_parsed"] is None
    assert captured["target_parsed"] is None
    assert out["coordinate_fetch"]["overlay_status"] == "UNAVAILABLE"
    assert out["superposition"] is None


def test_compare_parse_envelope_has_no_overlay_xyz(client: TestClient) -> None:
    from pathlib import Path
    import json

    fixture = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "rcsb_alignment_complete.json"
    payload = json.loads(fixture.read_text(encoding="utf-8"))
    response = client.post("/api/v1/compare/parse", json={"payload": payload})
    assert response.status_code == 200, response.text
    result = response.json()["result"]
    assert result.get("superposition") in (None, {}, [])
    assert "transformed_target_atoms" not in result
    assert result.get("n_aligned_residue_pairs")
