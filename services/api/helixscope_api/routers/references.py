"""Local assembly/reference catalog and safe catalogued install."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core.references import admit_catalog_download, list_references, start_catalog_download

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.schemas import ApiEnvelope, ReferenceAdmissionRequest, ReferenceInstallRequest
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1", tags=["references"])


@router.get(
    "/references",
    operation_id="listReferences",
    summary="Local reference catalog with READY/NOT_INSTALLED/RESOURCE_LIMIT",
    response_model=ApiEnvelope,
)
def references_list(request: Request) -> ApiEnvelope:
    """Browser must not decide readiness. Status comes from core."""
    rows = list_references()
    return envelope(request, {"references": rows}, scientific_status="COMPUTED")


@router.post(
    "/references/admission",
    operation_id="referencesAdmitDownload",
    summary="Resource admission for a catalogued public assembly",
    response_model=ApiEnvelope,
)
def references_admission(request: Request, body: ReferenceAdmissionRequest) -> ApiEnvelope:
    """Does not download. TEST REFERENCE is not a public genome."""
    reject_url_or_path(body.assembly_id, field="assembly_id")
    payload = admit_catalog_download(body.assembly_id)
    return envelope(request, payload)


@router.post(
    "/references/install",
    operation_id="referencesInstallCatalog",
    summary="Start a catalogued reference download after admission PROCEED",
    response_model=ApiEnvelope,
)
def references_install(request: Request, body: ReferenceInstallRequest) -> ApiEnvelope:
    """Explicit confirm required. RESOURCE_LIMIT is not bypassed. No arbitrary URL."""
    reject_url_or_path(body.assembly_id, field="assembly_id")
    if not body.confirm:
        raise ApiError(400, "INVALID_INPUT", "Install requires explicit confirm=true from the user.")
    payload = start_catalog_download(body.assembly_id)
    return envelope(request, payload)
