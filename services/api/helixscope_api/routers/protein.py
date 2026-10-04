"""Protein physicochemical and structure-source transport. No invented SS from sequence."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core import analyze_protein, retrieve_uniprot_entry
from helixscope_core.structure import (
    BUNDLED_STRUCTURE_IDS,
    analyze_structure_coordinates,
    fetch_alphafold_prediction,
    load_bundled_deposited,
    load_deposited_macromolecule,
    read_bundled_mmcif,
    search_structure_sources,
    structure_engine_capabilities,
)

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import (
    ApiEnvelope,
    ProteinAlphafoldRequest,
    ProteinAnalyzeRequest,
    ProteinStructureAnalyzeRequest,
    ProteinStructureLoadRequest,
    ProteinStructureSearchRequest,
    ProteinUniprotRequest,
)
from helixscope_api.science_attach import attach_explanations
from helixscope_api.security_inputs import reject_url_or_path

router = APIRouter(prefix="/api/v1/protein", tags=["protein"])


@router.post(
    "/analyze",
    operation_id="proteinAnalyze",
    summary="Physicochemical protein metrics",
    response_model=ApiEnvelope,
)
def protein_analyze(request: Request, body: ProteinAnalyzeRequest) -> ApiEnvelope:
    """ProtParam-style metrics. Sequence SS predictors remain UNAVAILABLE."""
    payload = analyze_protein(body.sequence, ph=body.ph)
    if body.include_explanation:
        attach_explanations(
            payload,
            include=True,
            kind="protein",
            domain="protein",
            values={
                "status": payload.get("status"),
                "length": payload.get("length"),
                "molecular_weight_kda": payload.get("molecular_weight_kda"),
                "isoelectric_point": payload.get("isoelectric_point"),
                "gravy": payload.get("gravy"),
                "instability_index": payload.get("instability_index"),
            },
        )
    return envelope(request, payload)


@router.post(
    "/uniprot",
    operation_id="proteinUniprot",
    summary="Retrieve UniProtKB annotation by accession",
    response_model=ApiEnvelope,
)
def protein_uniprot(request: Request, body: ProteinUniprotRequest) -> ApiEnvelope:
    """Official rest.uniprot.org entry. Does not guess identity from a pasted sequence."""
    reject_url_or_path(body.accession, field="accession")
    payload = retrieve_uniprot_entry(body.accession, include_domains=body.include_domains)
    if body.include_explanation:
        attach_explanations(
            payload,
            include=True,
            kind="ncbi",
            domain="ncbi",
            values={
                "status": "RETRIEVED",
                "accession": payload.get("accession"),
                "organism": payload.get("organism"),
                "title": payload.get("protein_name"),
            },
        )
    return envelope(request, payload, scientific_status="RETRIEVED")


@router.get(
    "/structure/capabilities",
    operation_id="proteinStructureCapabilities",
    summary="Protein-structure engine capability contract",
    response_model=ApiEnvelope,
)
def protein_structure_capabilities(request: Request) -> ApiEnvelope:
    """DSSP, PDB-REDO DSSP, STRIDE, EDTSurf, SASA, AlphaFold retrieval, RCSB search."""
    return envelope(request, structure_engine_capabilities(), scientific_status="COMPUTED")


@router.post(
    "/structure/search",
    operation_id="proteinStructureSearch",
    summary="Search RCSB / UniProt / AlphaFold structure sources",
    response_model=ApiEnvelope,
)
def protein_structure_search(request: Request, body: ProteinStructureSearchRequest) -> ApiEnvelope:
    """Official RCSB Search API when sequence search is requested. No HTML scrape."""
    if body.pdb_id:
        reject_url_or_path(body.pdb_id, field="pdb_id")
    if body.uniprot:
        reject_url_or_path(body.uniprot, field="uniprot")
    if body.alphafold_id:
        reject_url_or_path(body.alphafold_id, field="alphafold_id")
    payload = search_structure_sources(
        body.sequence,
        pdb_id=body.pdb_id,
        uniprot=body.uniprot,
        alphafold_id=body.alphafold_id,
        search_pdb_by_sequence=body.search_pdb_by_sequence,
        lookup_uniprot=body.lookup_uniprot,
        lookup_alphafold=body.lookup_alphafold,
        identity_cutoff=body.identity_cutoff,
        molecule=body.molecule,
    )
    return envelope(request, payload)


@router.post(
    "/structure/load",
    operation_id="proteinStructureLoad",
    summary="Load a bundled fixture or an official RCSB PDB identifier",
    response_model=ApiEnvelope,
)
def protein_structure_load(request: Request, body: ProteinStructureLoadRequest) -> ApiEnvelope:
    """Bundled 1CRN/1BNA/1RNA never fetch. Other IDs use the official RCSB files API."""
    reject_url_or_path(body.structure_id, field="structure_id")
    ident = str(body.structure_id or "").strip().upper()
    if body.bundled or ident in BUNDLED_STRUCTURE_IDS:
        payload = load_bundled_deposited(ident)
    else:
        payload = load_deposited_macromolecule(source="RCSB PDB", structure_id=ident)
    return envelope(request, payload)


@router.post(
    "/structure/analyze",
    operation_id="proteinStructureAnalyze",
    summary="Coordinate analyses on a bundled deposited structure",
    response_model=ApiEnvelope,
)
def protein_structure_analyze(request: Request, body: ProteinStructureAnalyzeRequest) -> ApiEnvelope:
    """SASA/DSSP/STRIDE/EDTSurf on deposited coordinates. No sequence SS invention."""
    parsed = load_bundled_deposited(body.structure_id)
    text = read_bundled_mmcif(f"{body.structure_id}.cif")
    payload = analyze_structure_coordinates(
        parsed,
        run_sasa=body.run_sasa,
        run_dssp=body.run_dssp,
        run_dssp_remote=body.run_dssp_remote,
        run_stride=body.run_stride,
        run_edtsurf=body.run_edtsurf,
        structure_text=text,
    )
    return envelope(request, payload)


@router.post(
    "/structure/alphafold",
    operation_id="proteinAlphafoldRetrieve",
    summary="Retrieve an AlphaFold DB prediction by UniProt accession",
    response_model=ApiEnvelope,
)
def protein_alphafold(request: Request, body: ProteinAlphafoldRequest) -> ApiEnvelope:
    """Deposited AFDB coordinates. Does not run AlphaFold locally."""
    reject_url_or_path(body.uniprot, field="uniprot")
    payload = fetch_alphafold_prediction(body.uniprot)
    return envelope(request, payload)
