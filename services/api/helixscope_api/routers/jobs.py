"""Common job status. Execution status is not scientific status."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.jobs import LocalMemoryJobRunner
from helixscope_api.schemas import ApiEnvelope

router = APIRouter(prefix="/api/v1/jobs", tags=["jobs"])


@router.get(
    "/{job_id}",
    operation_id="jobGet",
    summary="Get local job execution status and result",
    response_model=ApiEnvelope,
)
def job_get(job_id: str, request: Request) -> ApiEnvelope:
    """NON-DURABLE local runner. Progress percentages are not invented."""
    runner: LocalMemoryJobRunner = request.app.state.jobs
    record = runner.get(job_id)
    if record is None:
        raise ApiError(404, "JOB_NOT_FOUND", "Job was not found in this process.")
    snap = runner.snapshot(record)
    return envelope(
        request,
        snap,
        scientific_status=snap.get("scientific_status"),
        execution_status=str(snap.get("execution_status") or ""),
        drop_heavy=True,
    )


@router.delete(
    "/{job_id}",
    operation_id="jobCancel",
    summary="Request job cancellation if actually supported",
    response_model=ApiEnvelope,
)
def job_cancel(job_id: str, request: Request) -> ApiEnvelope:
    """CANCELLED only if the job was still queued. Running engines are not fake-cancelled."""
    runner: LocalMemoryJobRunner = request.app.state.jobs
    outcome = runner.cancel(job_id)
    code = str(outcome.get("code") or "")
    if code == "JOB_NOT_FOUND":
        raise ApiError(404, "JOB_NOT_FOUND", "Job was not found in this process.")
    if code == "CANCEL_NOT_SUPPORTED":
        raise ApiError(
            409,
            "CANCEL_NOT_SUPPORTED",
            "This job cannot be cancelled safely while running or after completion.",
            {"execution_status": outcome.get("execution_status")},
        )
    return envelope(
        request,
        outcome,
        execution_status=str(outcome.get("execution_status") or "CANCELLED"),
        scientific_status=None,
    )
