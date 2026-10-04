"""Variant identity and optional remote annotation. No HelixScope diagnosis."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from helixscope_core import identify_variant
from helixscope_core.variants import explore_variant

from helixscope_api.dependencies import envelope
from helixscope_api.job_handlers import run_variant_explore_job
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import ApiEnvelope, JobAccepted, VariantExploreRequest, VariantIdentifyRequest
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/variants", tags=["variants"])


@router.post(
    "/identify",
    operation_id="variantIdentify",
    summary="Parse variant identity without effect prediction",
    response_model=ApiEnvelope,
)
def variant_identify(request: Request, body: VariantIdentifyRequest) -> ApiEnvelope:
    """Assembly is part of identity. effect_status remains NOT_COMPUTED."""
    reject_url_or_path(body.text, field="text")
    reject_url_or_path(body.assembly, field="assembly")
    payload = identify_variant(body.text, assembly=body.assembly)
    return envelope(request, payload)


@router.post(
    "/explore",
    operation_id="variantExplore",
    summary="Identity plus optional remote layers (defaults off)",
    response_model=ApiEnvelope,
)
def variant_explore(request: Request, body: VariantExploreRequest) -> ApiEnvelope:
    """Does not set helixscope_pathogenic. Remote layers are explicit flags."""
    reject_url_or_path(body.text, field="text")
    reject_url_or_path(body.assembly, field="assembly")
    variant = identify_variant(body.text, assembly=body.assembly)
    payload = explore_variant(
        variant=variant,
        email=body.email,
        enable_vep=body.enable_vep,
        enable_clinvar=body.enable_clinvar,
        enable_domains=body.enable_domains,
    )
    return envelope(request, payload)


@router.post(
    "/jobs",
    operation_id="variantExploreSubmitJob",
    summary="Queue remote variant annotation",
    response_model=JobAccepted,
    status_code=202,
)
def variant_explore_job(request: Request, body: VariantExploreRequest) -> JSONResponse:
    """Use when VEP/ClinVar/InterPro are requested. Local runner is non-durable."""
    reject_url_or_path(body.text, field="text")
    reject_url_or_path(body.assembly, field="assembly")
    accepted = queue_job(
        request,
        "variant_explore",
        body.model_dump(),
        run_variant_explore_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())
