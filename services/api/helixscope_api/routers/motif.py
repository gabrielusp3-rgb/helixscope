"""Motif search transport. Hits are pattern occurrences, not function."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core.motif import analyze_motif

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope, MotifSearchRequest
from helixscope_api.science_attach import attach_explanations

router = APIRouter(prefix="/api/v1/motif", tags=["motif"])


@router.post(
    "/search",
    operation_id="motifSearch",
    summary="IUPAC/exact motif search",
    response_model=ApiEnvelope,
)
def motif_search(request: Request, body: MotifSearchRequest) -> ApiEnvelope:
    """Coordinates are 0-based half-open. Hits do not prove biological function."""
    payload = analyze_motif(body.sequence, body.pattern, molecule=body.molecule)
    if body.include_explanation:
        attach_explanations(
            payload,
            include=True,
            kind="motif",
            domain="motif",
            values={
                "pattern": body.pattern,
                "n_hits": payload.get("n_hits"),
                "status": payload.get("status"),
            },
        )
    return envelope(request, payload)
