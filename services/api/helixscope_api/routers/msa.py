"""MSA transport. Prealigned FASTA is not an engine MSA."""

from __future__ import annotations

from fastapi import APIRouter, File, Request, UploadFile
from fastapi.responses import JSONResponse

from helixscope_core.msa import import_alignment, tool_availability

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.job_handlers import run_msa_engine_job
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import ApiEnvelope, JobAccepted, MsaEngineJobRequest, MsaPrealignedRequest

router = APIRouter(prefix="/api/v1/msa", tags=["msa"])


@router.get(
    "/availability",
    operation_id="msaAvailability",
    summary="Declared MSA backends",
    response_model=ApiEnvelope,
)
def msa_availability(request: Request) -> ApiEnvelope:
    """EBI Clustal Omega and local allowlisted binaries. No arbitrary executable."""
    return envelope(request, tool_availability(), scientific_status="COMPUTED")


@router.post(
    "/prealigned",
    operation_id="msaImportPrealigned",
    summary="Import user-declared prealigned FASTA",
    response_model=ApiEnvelope,
)
def msa_prealigned(request: Request, body: MsaPrealignedRequest) -> ApiEnvelope:
    """Does not run Clustal/MAFFT/MUSCLE."""
    payload = import_alignment(body.fasta, fmt=body.format)
    return envelope(request, payload)


@router.post(
    "/prealigned/upload",
    operation_id="msaImportPrealignedUpload",
    summary="Import prealigned FASTA from an uploaded file",
    response_model=ApiEnvelope,
)
def msa_prealigned_upload(request: Request, file: UploadFile = File(...)) -> ApiEnvelope:
    """Filename is never used as a server path. Bytes are size-limited."""
    settings = request.app.state.settings
    raw = file.file.read(int(settings.max_upload_bytes) + 1)
    if len(raw) > int(settings.max_upload_bytes):
        raise ApiError(413, "RESOURCE_LIMIT", "Uploaded FASTA exceeds the configured size limit.")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ApiError(422, "PARSE_ERROR", "Uploaded FASTA is not valid UTF-8.") from exc
    payload = import_alignment(text, fmt="auto")
    return envelope(request, payload)


@router.post(
    "/jobs",
    operation_id="msaEngineSubmitJob",
    summary="Queue an MSA engine job",
    response_model=JobAccepted,
    status_code=202,
)
def msa_engine_job(request: Request, body: MsaEngineJobRequest) -> JSONResponse:
    """Computed MSA, distinct from user-supplied prealigned FASTA."""
    backend = str(body.backend or "")
    if "/" in backend or "\\" in backend or ".." in backend or backend.endswith(".exe"):
        raise ApiError(400, "INVALID_INPUT", "MSA backend must be a registered identifier.")
    accepted = queue_job(
        request,
        "msa_engine",
        body.model_dump(),
        run_msa_engine_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())
