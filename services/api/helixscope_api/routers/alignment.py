"""Pairwise alignment transport. Method names stay Needleman-Wunsch / Smith-Waterman."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import pairwise_align
from helixscope_core.alignment import MAX_DOTPLOT_LENGTH, dotplot_matrix
from helixscope_core.sequence_input import resolve_pasted_sequence

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope, PairwiseAlignRequest
from helixscope_api.science_attach import attach_explanations

router = APIRouter(prefix="/api/v1/alignment", tags=["alignment"])


@router.post(
    "/pairwise",
    operation_id="alignmentPairwise",
    summary="Needleman-Wunsch or Smith-Waterman pairwise alignment",
    response_model=ApiEnvelope,
)
def alignment_pairwise(request: Request, body: PairwiseAlignRequest) -> ApiEnvelope:
    """global = Needleman-Wunsch. local = Smith-Waterman. Not BLAST."""
    payload = pairwise_align(body.seq1, body.seq2, mode=body.mode, translate_nucleic=body.translate_nucleic)
    if body.include_dotplot:
        try:
            first = str(resolve_pasted_sequence(body.seq1).get("sequence") or "")
            second = str(resolve_pasted_sequence(body.seq2).get("sequence") or "")
            if max(len(first), len(second)) > MAX_DOTPLOT_LENGTH:
                payload["dotplot_status"] = "UNAVAILABLE"
                payload["dotplot_reason"] = (
                    f"Dot plot is limited to {MAX_DOTPLOT_LENGTH} residues per sequence."
                )
            else:
                payload["dotplot"] = dotplot_matrix(first, second)
                payload["dotplot_status"] = "COMPUTED"
        except (ValueError, TypeError) as exc:
            payload["dotplot_status"] = "UNAVAILABLE"
            payload["dotplot_reason"] = str(exc)
    if body.include_explanation:
        attach_explanations(
            payload,
            include=True,
            kind="alignment",
            domain="alignment",
            values={
                "status": payload.get("status"),
                "method": payload.get("method"),
                "identity_pct": payload.get("identity_pct"),
                "mode": body.mode,
            },
        )
    return envelope(request, payload)
