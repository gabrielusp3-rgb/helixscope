"""LOCAL MIGRATION JOB BACKEND. NON-DURABLE. NOT A FINAL PRODUCTION QUEUE.

Process restart loses in-memory jobs. Later workers can implement JobBackend.
"""

from __future__ import annotations

import logging
import secrets
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable, Optional, Protocol

from fastapi import Request

from helixscope_api.errors import ApiError
from helixscope_api.schemas import JobAccepted

LOGGER = logging.getLogger("helixscope_api.jobs")

QUEUED = "QUEUED"
RUNNING = "RUNNING"
COMPLETED = "COMPLETED"
FAILED = "FAILED"
CANCELLED = "CANCELLED"
TIMEOUT = "TIMEOUT"
RESOURCE_LIMIT = "RESOURCE_LIMIT"

JobFn = Callable[[dict[str, Any], "JobRecord"], dict[str, Any]]


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class JobRecord:
    job_id: str
    kind: str
    execution_status: str
    created_at: str
    params: dict[str, Any]
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    scientific_status: Optional[str] = None
    result: Optional[dict[str, Any]] = None
    error: Optional[dict[str, Any]] = None
    cancel_requested: bool = False
    future: Any = field(default=None, repr=False)


class JobBackend(Protocol):
    def submit(self, kind: str, params: dict[str, Any], handler: JobFn) -> JobRecord: ...
    def get(self, job_id: str) -> Optional[JobRecord]: ...
    def cancel(self, job_id: str) -> dict[str, Any]: ...


class LocalMemoryJobRunner:
    """In-process thread pool. Bounded. Not durable."""

    def __init__(self, *, max_concurrent: int = 2, queue_capacity: int = 8) -> None:
        self.max_concurrent = int(max_concurrent)
        self.queue_capacity = int(queue_capacity)
        self._jobs: dict[str, JobRecord] = {}
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(max_workers=self.max_concurrent, thread_name_prefix="hs-job")

    def submit(self, kind: str, params: dict[str, Any], handler: JobFn) -> JobRecord:
        with self._lock:
            active = sum(
                1
                for job in self._jobs.values()
                if job.execution_status in {QUEUED, RUNNING}
            )
            if active >= self.queue_capacity:
                raise RuntimeError("JOB_QUEUE_FULL")
            job_id = secrets.token_hex(16)
            record = JobRecord(
                job_id=job_id,
                kind=str(kind),
                execution_status=QUEUED,
                created_at=_utc(),
                params=dict(params),
            )
            self._jobs[job_id] = record
            record.future = self._pool.submit(self._run, record, handler)
            return record

    def get(self, job_id: str) -> Optional[JobRecord]:
        with self._lock:
            found = self._jobs.get(job_id)
            return found

    def cancel(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                return {"ok": False, "code": "JOB_NOT_FOUND"}
            if record.execution_status in {COMPLETED, FAILED, CANCELLED, TIMEOUT}:
                return {"ok": False, "code": "CANCEL_NOT_SUPPORTED", "execution_status": record.execution_status}
            if record.execution_status == QUEUED:
                record.cancel_requested = True
                record.execution_status = CANCELLED
                record.finished_at = _utc()
                return {"ok": True, "code": "CANCELLED", "execution_status": CANCELLED}
            record.cancel_requested = True
            return {"ok": False, "code": "CANCEL_NOT_SUPPORTED", "execution_status": record.execution_status}

    def _run(self, record: JobRecord, handler: JobFn) -> None:
        if record.cancel_requested:
            record.execution_status = CANCELLED
            record.finished_at = _utc()
            return
        record.execution_status = RUNNING
        record.started_at = _utc()
        started = time.perf_counter()
        try:
            payload = handler(record.params, record)
            if record.cancel_requested:
                record.execution_status = CANCELLED
                record.finished_at = _utc()
                return
            record.result = payload
            record.scientific_status = str(payload.get("status") or payload.get("scientific_status") or "COMPUTED")
            record.execution_status = COMPLETED
        except Exception as exc:
            LOGGER.warning("job_id=%s kind=%s failed type=%s", record.job_id, record.kind, type(exc).__name__)
            category = str(getattr(exc, "category", "") or "")
            if category.upper() == "TIMEOUT" or "timeout" in str(exc).lower():
                record.execution_status = TIMEOUT
            elif category.upper() == "RESOURCE_LIMIT":
                record.execution_status = RESOURCE_LIMIT
            else:
                record.execution_status = FAILED
            record.error = {
                "code": category.upper() or "ENGINE_FAILED",
                "message": str(exc),
                "type": type(exc).__name__,
            }
        record.finished_at = _utc()
        LOGGER.info(
            "job_id=%s kind=%s execution_status=%s elapsed_ms=%s",
            record.job_id,
            record.kind,
            record.execution_status,
            int((time.perf_counter() - started) * 1000),
        )

    def shutdown(self) -> None:
        self._pool.shutdown(wait=False, cancel_futures=True)

    def snapshot(self, record: JobRecord) -> dict[str, Any]:
        return {
            "job_id": record.job_id,
            "kind": record.kind,
            "execution_status": record.execution_status,
            "scientific_status": record.scientific_status,
            "created_at": record.created_at,
            "started_at": record.started_at,
            "finished_at": record.finished_at,
            "result": record.result,
            "error": record.error,
            "backend": "LOCAL_MIGRATION_JOB_BACKEND",
            "durable": False,
            "note": "NON-DURABLE. Process restart loses local job state. Not a production queue.",
        }


def queue_job(request: Request, kind: str, params: dict[str, Any], handler: JobFn) -> JobAccepted:
    """Submit to the local runner. Queue overflow is RESOURCE_LIMIT, not a hang."""
    runner: LocalMemoryJobRunner = request.app.state.jobs
    try:
        record = runner.submit(kind, params, handler)
    except RuntimeError as exc:
        if str(exc) == "JOB_QUEUE_FULL":
            raise ApiError(
                409,
                "RESOURCE_LIMIT",
                "Local job queue is full.",
            ) from exc
        raise
    return JobAccepted(
        request_id=str(getattr(request.state, "request_id", "") or ""),
        job_id=record.job_id,
        kind=record.kind,
        execution_status=record.execution_status,
        poll_url=f"/api/v1/jobs/{record.job_id}",
    )
