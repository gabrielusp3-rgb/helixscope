"""BLAST transport. Local BLAST+ is not nt/nr. Remote NCBI is a separate backend."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from helixscope_core.blast import local_blast_availability, ncbi_blast_availability

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.job_handlers import run_blast_local_job, run_blast_remote_submit
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import ApiEnvelope, BlastLocalJobRequest, BlastRemoteJobRequest, JobAccepted
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/blast", tags=["blast"])


@router.get(
    "/availability",
    operation_id="blastAvailability",
    summary="Local BLAST+ and NCBI remote BLAST capability",
    response_model=ApiEnvelope,
)
def blast_availability(request: Request) -> ApiEnvelope:
    """Local database is the declared tiny test prefix, not nt/nr."""
    local = local_blast_availability()
    remote = ncbi_blast_availability()
    payload = {
        "local": local,
        "remote": remote,
        "local_is_not_nt_nr": True,
    }
    return envelope(request, payload, scientific_status="COMPUTED")


@router.post(
    "/jobs/local",
    operation_id="blastLocalSubmitJob",
    summary="Queue local BLAST+ against the declared test database",
    response_model=JobAccepted,
    status_code=202,
)
def blast_local_job(request: Request, body: BlastLocalJobRequest) -> JSONResponse:
    """Refuses to queue when BLAST+ or the declared database is missing."""
    avail = local_blast_availability()
    if not avail.get("available"):
        raise ApiError(
            503,
            "ENGINE_NOT_INSTALLED",
            str(avail.get("reason") or "Local BLAST+ is not available."),
        )
    accepted = queue_job(
        request,
        "blast_local",
        body.model_dump(),
        run_blast_local_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())


@router.post(
    "/jobs/remote",
    operation_id="blastRemoteSubmitJob",
    summary="Queue NCBI remote BLAST Put (RID)",
    response_model=JobAccepted,
    status_code=202,
)
def blast_remote_job(request: Request, body: BlastRemoteJobRequest) -> JSONResponse:
    """Remote BLAST is explicit. Timeout does not fall back to local BLAST+."""
    reject_url_or_path(body.query, field="query")
    accepted = queue_job(
        request,
        "blast_remote",
        body.model_dump(),
        run_blast_remote_submit,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())
