"""Deterministic explanation. Does not recompute science."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import explain_result

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope, ExplainRequest

router = APIRouter(prefix="/api/v1", tags=["explain"])


@router.post(
    "/explain",
    operation_id="explainResult",
    summary="Deterministic explanation of an already computed result",
    response_model=ApiEnvelope,
)
def explain_endpoint(request: Request, body: ExplainRequest) -> ApiEnvelope:
    """Core explain_result. FastAPI does not generate text."""
    payload = explain_result(body.kind, body.values)
    return envelope(request, payload, scientific_status=str(payload.get("status") or ""))
