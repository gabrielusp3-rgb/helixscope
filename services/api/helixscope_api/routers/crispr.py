"""CRISPR guide design and Cas-OFFinder jobs. TEST REFERENCE is not GRCh38."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from helixscope_core import get_capabilities
from helixscope_core.crispr import analyze_crispr, list_cas_systems, model_availability
from helixscope_core.genome import detect_cas_offinder, full_genome_search_preflight
from helixscope_core.serialize import redact_filesystem_path

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.job_handlers import run_casoffinder_job
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import ApiEnvelope, CasOffinderJobRequest, CrisprGuidesRequest, JobAccepted
from helixscope_api.science_attach import attach_explanations
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/crispr", tags=["crispr"])


@router.get(
    "/availability",
    operation_id="crisprAvailability",
    summary="Advanced CRISPR model availability",
    response_model=ApiEnvelope,
)
def crispr_availability(request: Request) -> ApiEnvelope:
    """Azimuth/DeepHF/Rule Set 2 remain unavailable offline."""
    payload = model_availability()
    return envelope(request, payload, scientific_status="UNAVAILABLE")


@router.get(
    "/systems",
    operation_id="crisprListSystems",
    summary="Canonical Cas-system registry",
    response_model=ApiEnvelope,
)
def crispr_systems(request: Request) -> ApiEnvelope:
    """UI dropdowns must use canonical_key. Informal tokens are invalid."""
    payload = list_cas_systems()
    return envelope(request, payload, scientific_status="COMPUTED")


@router.post(
    "/guides",
    operation_id="crisprGuides",
    summary="Design guides on a provided sequence",
    response_model=ApiEnvelope,
)
def crispr_guides(request: Request, body: CrisprGuidesRequest) -> ApiEnvelope:
    """Provided-sequence off-target is not genome-wide. None scores stay null."""
    payload = analyze_crispr(
        body.sequence,
        body.cas_system,
        max_mismatches=body.max_mismatches,
        run_off_target=body.run_off_target,
    )
    if body.include_explanation:
        attach_explanations(
            payload,
            include=True,
            kind="crispr",
            domain="crispr",
            values={
                "status": "HEURISTIC",
                "cas_system": body.cas_system,
                "n_guides": payload.get("n_guides"),
            },
        )
    return envelope(request, payload)


@router.get(
    "/casoffinder/availability",
    operation_id="casOffinderAvailability",
    summary="Cas-OFFinder and OpenCL capability",
    response_model=ApiEnvelope,
)
def casoffinder_availability(request: Request) -> ApiEnvelope:
    """Truthful capability. Missing OpenCL is not simulated hits."""
    detected = detect_cas_offinder()
    caps = get_capabilities(refresh=False)
    opencl = caps.get("opencl") if isinstance(caps.get("opencl"), dict) else {}
    sanitized = {
        "available": bool(detected.get("available")),
        "version": detected.get("version") or "",
        "status": detected.get("status"),
        "reason": detected.get("reason"),
        "binary": redact_filesystem_path(detected.get("report_path") or detected.get("path")),
        "opencl_ready": bool(opencl.get("ready")),
        "opencl_blocker": opencl.get("blocker") or "",
        "not_grch38": True,
        "test_reference_is_not_genome_wide": True,
    }
    status = "LIVE_VALIDATED" if sanitized["available"] and sanitized["opencl_ready"] else "UNAVAILABLE"
    return envelope(request, sanitized, scientific_status=status)


@router.post(
    "/casoffinder/jobs",
    operation_id="casOffinderSubmitJob",
    summary="Submit a Cas-OFFinder job after capability checks",
    response_model=JobAccepted,
    status_code=202,
)
def casoffinder_submit(request: Request, body: CasOffinderJobRequest) -> JSONResponse:
    """Do not queue a job that cannot execute. TEST REFERENCE is not FULL GENOME."""
    reject_url_or_path(body.assembly_id, field="assembly_id")
    reject_url_or_path(body.guide_sequence, field="guide_sequence")
    detected = detect_cas_offinder()
    if not detected.get("available"):
        raise ApiError(503, "ENGINE_NOT_INSTALLED", str(detected.get("reason") or "Cas-OFFinder is not installed."))
    caps = get_capabilities(refresh=False)
    opencl = caps.get("opencl") if isinstance(caps.get("opencl"), dict) else {}
    if not opencl.get("ready"):
        raise ApiError(
            503,
            "ENGINE_NOT_INSTALLED",
            str(opencl.get("blocker") or "OpenCL is not ready for Cas-OFFinder."),
        )
    pre = full_genome_search_preflight(body.assembly_id)
    test_ok = bool(pre.get("is_test_reference") and pre.get("test_search_allowed"))
    if not pre.get("allowed") and not test_ok:
        reason = str(pre.get("reason") or "Cas-OFFinder search is not admitted.")
        code = "RESOURCE_LIMIT" if "RESOURCE_LIMIT" in reason.upper() else "ENGINE_NOT_INSTALLED"
        status = 409 if code == "RESOURCE_LIMIT" else 503
        raise ApiError(status, code, reason, {"is_test_reference": bool(pre.get("is_test_reference"))})
    accepted = queue_job(
        request,
        "cas_offinder",
        {
            "guide_sequence": body.guide_sequence,
            "assembly_id": body.assembly_id,
            "pam": body.pam,
            "cas_system": body.cas_system,
        },
        run_casoffinder_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())
