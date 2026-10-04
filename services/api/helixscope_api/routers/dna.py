"""DNA analysis transport. Science stays in helixscope_core."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import analyze_dna
from helixscope_core.dna import dna_3d_viewer

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import ApiEnvelope, DnaAnalyzeRequest, DnaHelix3dRequest
from helixscope_api.science_attach import attach_explanations

router = APIRouter(prefix="/api/v1/dna", tags=["dna"])


@router.post(
    "/analyze",
    operation_id="dnaAnalyze",
    summary="Analyze a DNA sequence",
    response_model=ApiEnvelope,
)
def dna_analyze(request: Request, body: DnaAnalyzeRequest) -> ApiEnvelope:
    """Call analyze_dna. Optional windowed profiles and deterministic explanation."""
    payload = analyze_dna(
        body.sequence,
        min_orf_length=body.min_orf_length,
        include_orfs=body.include_orfs,
        include_restriction=body.include_restriction,
        drop_overlapping_orfs=body.drop_overlapping_orfs,
        include_profiles=body.include_profiles,
        include_cpg=body.include_cpg,
        include_kmers=body.include_kmers,
        include_santalucia=body.include_santalucia,
        profile_window=body.profile_window,
        profile_step=body.profile_step,
        kmer_k=body.kmer_k,
        kmer_top_n=body.kmer_top_n,
        santalucia_na_m=body.santalucia_na_m,
        santalucia_oligo_nm=body.santalucia_oligo_nm,
    )
    computed = payload.get("status") == "COMPUTED"
    if body.include_explanation and computed:
        attach_explanations(
            payload,
            include=True,
            kind="dna",
            domain="dna",
            values={
                "status": payload.get("status"),
                "length": payload.get("length"),
                "gc_percent": payload.get("gc_content"),
                "at_percent": payload.get("at_content"),
                "shannon_entropy": payload.get("shannon_entropy"),
            },
        )
    return envelope(request, payload)


@router.post(
    "/helix-3d",
    operation_id="dnaHelix3d",
    summary="Illustrative B-DNA or deposited DNA 3D coordinates",
    response_model=ApiEnvelope,
)
def dna_helix_3d(request: Request, body: DnaHelix3dRequest) -> ApiEnvelope:
    """ILLUSTRATIVE helix is not experimental. 1BNA remains a deposited duplex."""
    payload = dna_3d_viewer(
        body.sequence,
        prefer=body.prefer,
        region_start=body.region_start,
        region_end=body.region_end,
        lod=body.lod,
        inspect_position=body.inspect_position,
        center_on_selected=body.center_on_selected,
    )
    status = str(payload.get("status") or "UNAVAILABLE")
    return envelope(request, payload, scientific_status=status, drop_heavy=True)
