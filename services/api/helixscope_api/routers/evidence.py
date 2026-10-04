"""Evidence pack transport. No majority-vote confidence."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import evidence_conflict_pack
from helixscope_core.evidence import evidence_item

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope, EvidencePackRequest

router = APIRouter(prefix="/api/v1/evidence", tags=["evidence"])


@router.post(
    "/pack",
    operation_id="evidencePack",
    summary="Build an evidence pack and list source conflicts",
    response_model=ApiEnvelope,
)
def evidence_pack(request: Request, body: EvidencePackRequest) -> ApiEnvelope:
    """confidence_score remains null. Conflicts are listed, not voted."""
    items = [
        evidence_item(
            field=item.field,
            value=item.value,
            source=item.source,
            evidence_status=item.evidence_status,
            retrieved_at_utc=item.retrieved_at_utc,
            identifier=item.identifier,
            mapping_status=item.mapping_status or "AVAILABLE",
            note=item.note,
        )
        for item in body.items
    ]
    pack = evidence_conflict_pack(items, kind=body.kind, inputs=body.inputs)
    return envelope(request, pack)
