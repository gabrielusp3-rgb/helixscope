"""API overhead samples for the Prompt 2 performance artifact."""

from __future__ import annotations

import json
import time
from pathlib import Path

from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "MIGRATION_API_PERFORMANCE.json"


def _timed(fn) -> tuple[object, float]:
    started = time.perf_counter()
    result = fn()
    return result, (time.perf_counter() - started) * 1000.0


def test_write_performance_artifact(client: TestClient) -> None:
    samples: dict[str, float] = {}
    _, samples["health_live_ms"] = _timed(lambda: client.get("/health/live"))
    _, samples["capabilities_ms"] = _timed(lambda: client.get("/api/v1/system/capabilities"))
    short = "ATGCATGCATGC"
    dna_resp, samples["dna_short_ms"] = _timed(
        lambda: client.post("/api/v1/dna/analyze", json={"sequence": short, "include_explanation": False})
    )
    assert dna_resp.status_code == 200
    prot = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    _, samples["protein_ms"] = _timed(
        lambda: client.post("/api/v1/protein/analyze", json={"sequence": prot, "include_explanation": False})
    )
    _, samples["alignment_ms"] = _timed(
        lambda: client.post(
            "/api/v1/alignment/pairwise",
            json={"seq1": "ACGTACGTACGT", "seq2": "ACGTACGGACGT", "mode": "global", "include_explanation": False},
        )
    )
    seq_50k = "ATGC" * 12_500
    resp_50k, samples["dna_50k_ms"] = _timed(
        lambda: client.post("/api/v1/dna/analyze", json={"sequence": seq_50k, "include_explanation": False})
    )
    assert resp_50k.status_code == 200
    samples["dna_50k_response_bytes"] = float(len(resp_50k.content))
    fixture, samples["structure_1crn_ms"] = _timed(
        lambda: client.post("/api/v1/structures/fixture", json={"structure_id": "1CRN"})
    )
    if fixture.status_code == 200:
        samples["structure_1crn_response_bytes"] = float(len(fixture.content))
    job, samples["job_submit_ms"] = _timed(
        lambda: client.post(
            "/api/v1/phylogeny/jobs",
            json={
                "fasta": ">seq_a\nACGTACGTACGTACGTACGT\n>seq_b\nACGTACGTACGTACGTTTTT\n",
                "method": "neighbor_joining",
            },
        )
    )
    if job.status_code == 202:
        job_id = job.json()["job_id"]
        _, samples["job_status_ms"] = _timed(lambda: client.get(f"/api/v1/jobs/{job_id}"))
    payload = {
        "product_version": "0.24.3-19",
        "api_contract": "v1",
        "note": "Wall times include core science plus FastAPI transport on this machine. Not a SLA.",
        "samples_ms": samples,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assert OUT.is_file()
