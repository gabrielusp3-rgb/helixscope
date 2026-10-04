"""RNA analysis transport. Folding remains ViennaRNA in core."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import analyze_rna, transcribe_dna
from helixscope_core.rna import rna_3d_viewer, translate_coding_composition

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import (
    ApiEnvelope,
    RnaAnalyzeRequest,
    RnaHelix3dRequest,
    RnaTranslateRequest,
    TranscribeRequest,
)
from helixscope_api.science_attach import attach_explanations

router = APIRouter(prefix="/api/v1/rna", tags=["rna"])


@router.post(
    "/analyze",
    operation_id="rnaAnalyze",
    summary="Analyze an RNA sequence",
    response_model=ApiEnvelope,
)
def rna_analyze(request: Request, body: RnaAnalyzeRequest) -> ApiEnvelope:
    """Composition plus optional codon metrics and ViennaRNA fold."""
    payload = analyze_rna(
        body.sequence,
        fold=body.fold,
        backend=body.backend,
        include_codon_metrics=body.include_codon_metrics,
        include_profiles=body.include_profiles,
        include_kmers=body.include_kmers,
        profile_window=body.profile_window,
        profile_step=body.profile_step,
        kmer_k=body.kmer_k,
    )
    computed = payload.get("status") == "COMPUTED"
    if body.include_explanation and computed:
        fold = payload.get("fold") if isinstance(payload.get("fold"), dict) else {}
        attach_explanations(
            payload,
            include=True,
            kind="rna_fold" if body.fold else "rna",
            domain="rna",
            values={
                "status": fold.get("status") if body.fold else payload.get("status"),
                "length": payload.get("length"),
                "gc_percent": payload.get("gc_content"),
                "mfe_kcal_mol": fold.get("mfe_kcal_mol"),
                "backend": fold.get("backend"),
                "tool": fold.get("tool"),
                "method": fold.get("method"),
                "shannon_entropy": payload.get("shannon_entropy"),
            },
        )
    return envelope(request, payload)


@router.post(
    "/transcribe",
    operation_id="rnaTranscribe",
    summary="Transcribe DNA to RNA (T to U)",
    response_model=ApiEnvelope,
)
def rna_transcribe(request: Request, body: TranscribeRequest) -> ApiEnvelope:
    """Exact transcription helper. Not splicing."""
    rna = transcribe_dna(body.sequence)
    return envelope(
        request,
        {"rna": rna, "status": "COMPUTED"},
        scientific_status="COMPUTED",
    )


@router.post(
    "/translate-coding",
    operation_id="rnaTranslateCoding",
    summary="Translate resolved coding RNA/CDS and report amino-acid composition",
    response_model=ApiEnvelope,
)
def rna_translate_coding(request: Request, body: RnaTranslateRequest) -> ApiEnvelope:
    """Does not treat arbitrary RNA as a gene. ORF/CDS context is required."""
    payload = translate_coding_composition(body.sequence)
    return envelope(request, payload)


@router.post(
    "/helix-3d",
    operation_id="rnaHelix3d",
    summary="Illustrative A-RNA or deposited RNA 3D coordinates",
    response_model=ApiEnvelope,
)
def rna_helix_3d(request: Request, body: RnaHelix3dRequest) -> ApiEnvelope:
    """ILLUSTRATIVE A-RNA is not experimental. MFE is not 3D. 1RNA is deposited."""
    payload = rna_3d_viewer(
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
