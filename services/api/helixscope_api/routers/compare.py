"""Comparison transport. Form values reach Core. No silent 1CRN/chain A/tm-align substitution."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from helixscope_core import explain_result, get_capabilities
from helixscope_core.compare import (
    compare_evolution,
    compare_evidence_from_variants,
    compare_guides,
    compare_modes,
    compare_proteins,
    compare_variants_from_text,
    parse_alignment_payload,
    structure_alignment_methods,
)

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.job_handlers import run_compare_remote_job, run_usalign_job
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import (
    ApiEnvelope,
    CompareEvolutionRequest,
    CompareEvidenceRequest,
    CompareGuidesRequest,
    CompareParseRequest,
    CompareProteinsRequest,
    CompareRemoteJobRequest,
    CompareUSalignJobRequest,
    CompareVariantsRequest,
    JobAccepted,
)
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/compare", tags=["compare"])


@router.get(
    "/modes",
    operation_id="compareListModes",
    summary="Declared comparison modes and dedicated endpoints",
    response_model=ApiEnvelope,
)
def compare_list_modes(request: Request) -> ApiEnvelope:
    """Six Core modes. Not one generic JSON endpoint."""
    return envelope(request, compare_modes(), scientific_status="COMPUTED")


@router.get(
    "/methods",
    operation_id="compareListMethods",
    summary="RCSB Alignment API methods HelixScope may submit",
    response_model=ApiEnvelope,
)
def compare_list_methods(request: Request) -> ApiEnvelope:
    """Official RCSB method catalog. No invented aligners."""
    return envelope(
        request,
        {"status": "COMPUTED", "methods": structure_alignment_methods()},
        scientific_status="COMPUTED",
    )


@router.post(
    "/parse",
    operation_id="compareParseAlignment",
    summary="Parse an RCSB COMPLETE alignment payload",
    response_model=ApiEnvelope,
)
def compare_parse(request: Request, body: CompareParseRequest) -> ApiEnvelope:
    """Keeps block RMSD and global RMSD as distinct fields."""
    record = parse_alignment_payload(body.payload)
    payload = dict(record)
    payload["rmsd_note"] = (
        "rmsd_block0_angstrom and rmsd_global_angstrom are distinct measurements. "
        "They must not be merged."
    )
    payload["explanation"] = explain_result(
        "compare",
        {
            "status": payload.get("status"),
            "method": payload.get("method"),
            "rmsd_global_angstrom": payload.get("rmsd_global_angstrom"),
            "rmsd_block0_angstrom": payload.get("rmsd_block0_angstrom"),
        },
    )
    return envelope(request, payload)


@router.post(
    "/jobs/remote",
    operation_id="compareRemoteSubmitJob",
    summary="Queue RCSB Alignment API comparison",
    response_model=JobAccepted,
    status_code=202,
)
def compare_remote_job(request: Request, body: CompareRemoteJobRequest) -> JSONResponse:
    """Explicit RCSB backend. Does not fall back to US-align or rewrite IDs/chains/method."""
    reject_url_or_path(body.reference_entry, field="reference_entry")
    reject_url_or_path(body.target_entry, field="target_entry")
    accepted = queue_job(
        request,
        "compare_rcsb",
        body.model_dump(),
        run_compare_remote_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())


@router.post(
    "/jobs/usalign",
    operation_id="compareUSalignSubmitJob",
    summary="Queue local US-align on bundled structures",
    response_model=JobAccepted,
    status_code=202,
)
def compare_usalign_job(request: Request, body: CompareUSalignJobRequest) -> JSONResponse:
    """Explicit local US-align backend. Does not fall back to RCSB or substitute 1CRN/1BNA."""
    caps = get_capabilities(refresh=False)
    tools = caps.get("tools") if isinstance(caps.get("tools"), dict) else {}
    usa = tools.get("usalign") if isinstance(tools.get("usalign"), dict) else {}
    if usa and not usa.get("available"):
        raise ApiError(503, "ENGINE_NOT_INSTALLED", str(usa.get("reason") or "US-align is not installed."))
    accepted = queue_job(
        request,
        "compare_usalign",
        body.model_dump(),
        run_usalign_job,
    )
    return JSONResponse(status_code=202, content=accepted.model_dump())


@router.post(
    "/variants",
    operation_id="compareVariants",
    summary="Compare two identified variants without pathogenicity ranking",
    response_model=ApiEnvelope,
)
def compare_variants_endpoint(request: Request, body: CompareVariantsRequest) -> ApiEnvelope:
    """Identity side-by-side. ClinVar/VEP layers are not invented here."""
    payload = compare_variants_from_text(body.text_a, body.assembly_a, body.text_b, body.assembly_b)
    return envelope(request, payload)


@router.post(
    "/proteins",
    operation_id="compareProteins",
    summary="Compare two protein sequences",
    response_model=ApiEnvelope,
)
def compare_proteins_endpoint(request: Request, body: CompareProteinsRequest) -> ApiEnvelope:
    """Sequence identity is not structure similarity."""
    payload = compare_proteins(
        sequence_a=body.sequence_a,
        sequence_b=body.sequence_b,
        identifier_a=body.identifier_a,
        identifier_b=body.identifier_b,
    )
    return envelope(request, payload)


@router.post(
    "/guides",
    operation_id="compareGuides",
    summary="Compare two CRISPR guides",
    response_model=ApiEnvelope,
)
def compare_guides_endpoint(request: Request, body: CompareGuidesRequest) -> ApiEnvelope:
    """No invented universal score. Genome-wide MIT/CFD stay unavailable unless supplied."""
    payload = compare_guides(body.guide_a, body.guide_b)
    return envelope(request, payload)


@router.post(
    "/evolution",
    operation_id="compareEvolution",
    summary="Inspect MSA columns and optional group-associated positions",
    response_model=ApiEnvelope,
)
def compare_evolution_endpoint(request: Request, body: CompareEvolutionRequest) -> ApiEnvelope:
    """Conservation among submitted sequences. Not dN/dS."""
    payload = compare_evolution(
        body.fasta,
        group_ids=body.group_ids,
        column=body.column,
        member_id=body.member_id,
    )
    return envelope(request, payload)


@router.post(
    "/evidence",
    operation_id="compareEvidence",
    summary="Evidence pack for two identified variants",
    response_model=ApiEnvelope,
)
def compare_evidence_endpoint(request: Request, body: CompareEvidenceRequest) -> ApiEnvelope:
    """Does not vote, rank, or diagnose."""
    payload = compare_evidence_from_variants(body.text_a, body.assembly_a, body.text_b, body.assembly_b)
    return envelope(request, payload)
