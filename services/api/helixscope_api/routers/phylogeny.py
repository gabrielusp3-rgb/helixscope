"""Phylogeny transport. NJ/UPGMA may run synchronously; IQ-TREE/FastTree are jobs."""

from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from helixscope_core.msa import import_prealigned_fasta
from helixscope_core.phylogeny import METHOD_FASTTREE, METHOD_IQTREE, METHOD_NJ, METHOD_UPGMA, annotate_taxonomy, infer_tree

from helixscope_api.dependencies import envelope
from helixscope_api.errors import ApiError
from helixscope_api.job_handlers import run_phylogeny_job
from helixscope_api.jobs import queue_job
from helixscope_api.schemas import ApiEnvelope, JobAccepted, PhylogenyInferRequest, PhylogenyTaxonomyRequest

router = APIRouter(prefix="/api/v1/phylogeny", tags=["phylogeny"])

SYNC_METHODS = frozenset({METHOD_NJ, METHOD_UPGMA})
JOB_METHODS = frozenset({METHOD_IQTREE, METHOD_FASTTREE})
ALL_METHODS = SYNC_METHODS | JOB_METHODS


@router.post(
    "/infer",
    operation_id="phylogenyInfer",
    summary="Infer NJ or UPGMA from prealigned FASTA",
    response_model=ApiEnvelope,
)
def phylogeny_infer(request: Request, body: PhylogenyInferRequest) -> ApiEnvelope:
    """IQ-TREE and FastTree must use the job endpoint."""
    method = str(body.method or "").strip()
    if method in JOB_METHODS:
        raise ApiError(
            400,
            "INVALID_INPUT",
            "IQ-TREE and FastTree must be submitted as jobs at POST /api/v1/phylogeny/jobs.",
        )
    if method not in SYNC_METHODS:
        raise ApiError(400, "INVALID_INPUT", "Phylogeny method must be neighbor_joining or upgma.")
    msa = import_prealigned_fasta(body.fasta)
    payload = infer_tree(msa, method=method, distance_model=body.distance_model)
    payload["source_msa_hash"] = body.alignment_hash or msa.get("alignment_hash")
    payload["source_msa_n_sequences"] = msa.get("n_sequences")
    payload["source_msa_engine"] = msa.get("tool") or msa.get("method")
    return envelope(request, payload)


@router.post(
    "/taxonomy",
    operation_id="phylogenyAttachTaxonomy",
    summary="Attach NCBI taxonomy to an existing tree without changing topology",
    response_model=ApiEnvelope,
)
def phylogeny_taxonomy(request: Request, body: PhylogenyTaxonomyRequest) -> ApiEnvelope:
    """Taxonomy is classification, not the phylogenetic tree."""
    payload = annotate_taxonomy(
        body.tree,
        email=body.email,
        declared_organisms=body.declared_organisms,
    )
    return envelope(request, payload)


@router.post(
    "/jobs",
    operation_id="phylogenySubmitJob",
    summary="Queue phylogeny (NJ, UPGMA, IQ-TREE, FastTree)",
    response_model=JobAccepted,
    status_code=202,
)
def phylogeny_job(request: Request, body: PhylogenyInferRequest) -> JSONResponse:
    """Local non-durable runner. IQ-TREE/FastTree execute off the request path."""
    method = str(body.method or "").strip()
    if method not in ALL_METHODS:
        raise ApiError(
            400,
            "INVALID_INPUT",
            "Phylogeny method must be neighbor_joining, upgma, iqtree_ml, or fasttree_ml.",
        )
    params = {"fasta": body.fasta, "method": method, "alignment_hash": body.alignment_hash}
    if method in SYNC_METHODS:
        params["distance_model"] = body.distance_model
    accepted = queue_job(request, "phylogeny", params, run_phylogeny_job)
    return JSONResponse(status_code=202, content=accepted.model_dump())
