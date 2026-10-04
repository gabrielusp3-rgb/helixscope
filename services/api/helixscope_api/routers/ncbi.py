"""NCBI Entrez transport. Uses the locked NCBIClient. Does not log email."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core.remote import NCBIClient, api_key_from_environment

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.schemas import ApiEnvelope, NcbiFetchRequest, NcbiSearchRequest
from helixscope_api.science_attach import attach_explanations
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/ncbi", tags=["ncbi"])

ALLOWED_DB = frozenset({"nucleotide", "protein", "gene", "pubmed"})


def _client(body_email: str, body_key: str | None) -> NCBIClient:
    key = (body_key or "").strip() or api_key_from_environment()
    return NCBIClient(body_email, api_key=key)


@router.post(
    "/fetch",
    operation_id="ncbiFetch",
    summary="Fetch one NCBI record by accession, GeneID, or PMID",
    response_model=ApiEnvelope,
)
def ncbi_fetch(request: Request, body: NcbiFetchRequest) -> ApiEnvelope:
    """Entrez fetch under the process lock. Email and API key are not returned."""
    reject_url_or_path(body.accession, field="accession")
    db = str(body.db or "nucleotide").strip().lower()
    if db not in ALLOWED_DB:
        raise ApiError(400, "INVALID_INPUT", "NCBI database is not offered.")
    record = _client(body.email, body.api_key).fetch_by_accession(body.accession, db=db)
    attach_explanations(
        record,
        include=True,
        kind="ncbi",
        domain="ncbi",
        values={
            "status": "RETRIEVED",
            "accession": record.get("accession"),
            "organism": record.get("organism"),
            "title": record.get("title") or record.get("description"),
            "database": record.get("database") or db,
        },
    )
    return envelope(request, record, scientific_status="RETRIEVED")


@router.post(
    "/search",
    operation_id="ncbiSearch",
    summary="Entrez esearch plus esummary",
    response_model=ApiEnvelope,
)
def ncbi_search(request: Request, body: NcbiSearchRequest) -> ApiEnvelope:
    """Search nucleotide, protein, gene, or pubmed. Not a generic database API."""
    reject_url_or_path(body.term, field="term")
    db = str(body.db or "nucleotide").strip().lower()
    if db not in ALLOWED_DB:
        raise ApiError(400, "INVALID_INPUT", "NCBI database is not offered.")
    payload = _client(body.email, body.api_key).search_records(
        body.term,
        db=db,
        retmax=body.retmax,
    )
    attach_explanations(
        payload,
        include=True,
        kind="ncbi",
        domain="ncbi",
        values={
            "status": payload.get("status") or "RETRIEVED",
            "query": payload.get("query"),
            "database": payload.get("database") or db,
            "count": payload.get("count"),
        },
    )
    return envelope(request, payload, scientific_status="RETRIEVED")
