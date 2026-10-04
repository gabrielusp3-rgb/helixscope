"""System metadata and capabilities."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import CORE_PACKAGE_VERSION, get_capabilities

from helixscope_api import API_CONTRACT_VERSION, PRODUCT_VERSION
from helixscope_api.contract import runtime_identity
from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope

router = APIRouter(prefix="/api/v1", tags=["system"])


@router.get(
    "",
    operation_id="apiRoot",
    summary="API contract identity",
)
def api_root() -> dict[str, str]:
    """Identify HelixScope API v1. Not a marketing page."""
    return {
        "name": "HelixScope API",
        "api_version": API_CONTRACT_VERSION,
        "product_version": PRODUCT_VERSION,
        "core_package_version": CORE_PACKAGE_VERSION,
        "docs": "/docs",
        "openapi": "/openapi.json",
    }


@router.get(
    "/system/identity",
    operation_id="systemIdentity",
    summary="Sanitized runtime and OpenAPI identity",
)
def system_identity(request: Request) -> dict[str, str]:
    """Frontend compares this OpenAPI hash to the types it was generated against."""
    payload = runtime_identity(request.app)
    return {
        "product_version": str(payload["product_version"]),
        "core_version": str(payload["core_version"]),
        "api_version": str(payload["api_version"]),
        "openapi_sha256": str(payload["openapi_sha256"]),
        "runtime_build_identity": str(payload["runtime_build_identity"]),
    }


@router.get(
    "/system/capabilities",
    operation_id="systemCapabilities",
    summary="Sanitized engine and reference capabilities",
    response_model=ApiEnvelope,
)
def system_capabilities(request: Request) -> ApiEnvelope:
    """Uses helixscope_core.get_capabilities. Paths remain basenames."""
    payload = get_capabilities(refresh=False)
    return envelope(request, payload, scientific_status="COMPUTED")
