"""Health probes. Live is process-alive. Ready does not require remote science."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_api import API_CONTRACT_VERSION, PRODUCT_VERSION

router = APIRouter(tags=["health"])


@router.get(
    "/health/live",
    operation_id="healthLive",
    summary="Liveness probe",
)
def health_live() -> dict[str, str]:
    """HTTP process is up. No engine or remote calls."""
    return {"status": "alive"}


@router.get(
    "/health/ready",
    operation_id="healthReady",
    summary="Readiness probe",
)
def health_ready(request: Request) -> dict[str, object]:
    """Core import and local job runner. NCBI/ClinVar/RCSB/GRCh38 are not required."""
    from helixscope_core import CORE_PACKAGE_VERSION

    jobs = getattr(request.app.state, "jobs", None)
    settings = getattr(request.app.state, "settings", None)
    return {
        "status": "ready",
        "api_version": API_CONTRACT_VERSION,
        "product_version": PRODUCT_VERSION,
        "core_package_version": CORE_PACKAGE_VERSION,
        "job_backend": "LOCAL_MIGRATION_JOB_BACKEND",
        "job_runner": jobs is not None,
        "environment": getattr(settings, "environment", "") if settings is not None else "",
    }
