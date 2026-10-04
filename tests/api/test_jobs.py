"""Local non-durable job lifecycle."""

from __future__ import annotations

import time

from fastapi.testclient import TestClient

from helixscope_api.config import ApiSettings
from helixscope_api.jobs import FAILED, TIMEOUT, LocalMemoryJobRunner
from helixscope_api.main import create_app
from helixscope_api.serialization import assert_strict_json

PREALIGNED = """>seq_a
ACGTACGTACGTACGTACGT
>seq_b
ACGTACGTACGTACGTTTTT
>seq_c
TTTTACGTACGTACGTACGT
>seq_d
ACGTACGTAAAAACGTACGT
"""


def test_phylogeny_job_completes(client: TestClient) -> None:
    submitted = client.post(
        "/api/v1/phylogeny/jobs",
        json={"fasta": PREALIGNED, "method": "neighbor_joining", "distance_model": "p_distance"},
    )
    assert submitted.status_code == 202, submitted.text
    body = submitted.json()
    assert body["execution_status"] in {"QUEUED", "RUNNING", "COMPLETED"}
    job_id = body["job_id"]
    status = ""
    snapshot = {}
    for _ in range(80):
        got = client.get(f"/api/v1/jobs/{job_id}")
        assert got.status_code == 200, got.text
        snapshot = got.json()
        assert_strict_json(snapshot)
        status = snapshot["execution_status"]
        assert "progress" not in snapshot["result"] or snapshot["result"].get("progress") in (None, "RUNNING")
        if status in {"COMPLETED", "FAILED", "TIMEOUT", "RESOURCE_LIMIT"}:
            break
        time.sleep(0.05)
    assert status == "COMPLETED"
    assert snapshot["scientific_status"] == "COMPUTED"
    assert snapshot["execution_status"] == "COMPLETED"
    assert snapshot["result"]["result"]["n_leaves"] == 4
    assert snapshot["result"]["durable"] is False


def test_job_not_found_and_iqtree_rejected_from_sync(client: TestClient) -> None:
    missing = client.get("/api/v1/jobs/does-not-exist")
    assert missing.status_code == 404
    assert missing.json()["error"]["code"] == "JOB_NOT_FOUND"
    sync = client.post(
        "/api/v1/phylogeny/infer",
        json={"fasta": PREALIGNED, "method": "iqtree_ml"},
    )
    assert sync.status_code == 400
    assert sync.json()["error"]["code"] == "INVALID_INPUT"


def test_job_timeout_and_failure_classification() -> None:
    runner = LocalMemoryJobRunner(max_concurrent=1, queue_capacity=2)

    def timeout_fn(params, record):
        exc = RuntimeError("engine timed out")
        exc.category = "TIMEOUT"
        raise exc

    def fail_fn(params, record):
        raise ValueError("bad tree")

    timed = runner.submit("timeout_fixture", {}, timeout_fn)
    failed = runner.submit("fail_fixture", {}, fail_fn)
    for _ in range(40):
        if timed.execution_status != "QUEUED" and timed.execution_status != "RUNNING":
            if failed.execution_status not in {"QUEUED", "RUNNING"}:
                break
        time.sleep(0.05)
    assert timed.execution_status == TIMEOUT
    assert failed.execution_status == FAILED
    assert failed.error and failed.error.get("message")
    runner.shutdown()


def test_queue_capacity_resource_limit() -> None:
    settings = ApiSettings(job_max_concurrent=1, job_queue_capacity=1)
    with TestClient(create_app(settings)) as client:
        def hang(params, record):
            time.sleep(1.5)
            return {"status": "COMPUTED"}

        client.app.state.jobs.submit("hang", {}, hang)
        blocked = client.post(
            "/api/v1/phylogeny/jobs",
            json={"fasta": PREALIGNED, "method": "neighbor_joining"},
        )
        assert blocked.status_code == 409
        assert blocked.json()["error"]["code"] == "RESOURCE_LIMIT"


def test_cancel_queued_only() -> None:
    runner = LocalMemoryJobRunner(max_concurrent=1, queue_capacity=4)

    def hang(params, record):
        time.sleep(2)
        return {"status": "COMPUTED"}

    first = runner.submit("hang", {}, hang)
    for _ in range(50):
        if first.execution_status == "RUNNING":
            break
        time.sleep(0.02)
    second = runner.submit("queued", {}, hang)
    outcome = runner.cancel(second.job_id)
    assert outcome["code"] == "CANCELLED"
    running = runner.cancel(first.job_id)
    assert running["code"] in {"CANCEL_NOT_SUPPORTED", "CANCELLED"}
    runner.shutdown()
