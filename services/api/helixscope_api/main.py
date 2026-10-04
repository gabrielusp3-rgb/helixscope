"""HelixScope FastAPI application factory. No science in this module."""

from __future__ import annotations

import logging

from fastapi import FastAPI
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from starlette.exceptions import HTTPException as StarletteHTTPException

from helixscope_api import API_CONTRACT_VERSION, PRODUCT_VERSION
from helixscope_api.config import ApiSettings
from helixscope_api.errors import (
    ApiError,
    DOMAIN_EXCEPTIONS,
    api_error_handler,
    domain_exception_handler,
    http_exception_handler,
    unhandled_exception_handler,
    validation_exception_handler,
)
from helixscope_api.lifespan import create_lifespan
from helixscope_api.middleware import BodyLimitMiddleware, RequestIdMiddleware, SecurityHeadersMiddleware
from helixscope_api.routers import (
    alignment,
    blast,
    compare,
    crispr,
    dna,
    evidence,
    explain,
    health,
    jobs,
    motif,
    msa,
    ncbi,
    phylogeny,
    protein,
    references,
    rna,
    structures,
    system,
    variants,
)


def create_app(settings: ApiSettings | None = None) -> FastAPI:
    """Build the ASGI app. Scientific work is delegated to helixscope_core."""
    cfg = settings or ApiSettings()
    logging.basicConfig(level=getattr(logging, str(cfg.log_level).upper(), logging.INFO))
    app = FastAPI(
        title="HelixScope API",
        version=API_CONTRACT_VERSION,
        summary="Typed scientific HTTP API for HelixScope core.",
        description=(
            f"Product version {PRODUCT_VERSION}. API contract {API_CONTRACT_VERSION}. "
            "This service validates transport input, calls helixscope_core, and serializes "
            "results. It does not independently calculate biology. "
            "Job backend: LOCAL MIGRATION JOB BACKEND, NON-DURABLE, NOT FINAL PRODUCTION QUEUE."
        ),
        lifespan=create_lifespan(cfg),
        docs_url="/docs",
        redoc_url="/redoc",
        openapi_url="/openapi.json",
    )
    app.state.settings = cfg

    for exc_type in DOMAIN_EXCEPTIONS:
        app.add_exception_handler(exc_type, domain_exception_handler)
    app.add_exception_handler(ApiError, api_error_handler)
    app.add_exception_handler(RequestValidationError, validation_exception_handler)
    app.add_exception_handler(StarletteHTTPException, http_exception_handler)
    app.add_exception_handler(Exception, unhandled_exception_handler)

    app.add_middleware(BodyLimitMiddleware, max_bytes=cfg.max_upload_bytes)
    app.add_middleware(SecurityHeadersMiddleware)
    app.add_middleware(RequestIdMiddleware)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=cfg.cors_origin_list(),
        allow_credentials=False,
        allow_methods=["GET", "POST", "DELETE", "OPTIONS"],
        allow_headers=["Authorization", "Content-Type", "X-Request-ID"],
    )

    app.include_router(health.router)
    app.include_router(system.router)
    app.include_router(references.router)
    app.include_router(dna.router)
    app.include_router(rna.router)
    app.include_router(protein.router)
    app.include_router(alignment.router)
    app.include_router(motif.router)
    app.include_router(explain.router)
    app.include_router(ncbi.router)
    app.include_router(crispr.router)
    app.include_router(variants.router)
    app.include_router(evidence.router)
    app.include_router(structures.router)
    app.include_router(compare.router)
    app.include_router(msa.router)
    app.include_router(phylogeny.router)
    app.include_router(blast.router)
    app.include_router(jobs.router)
    return app


app = create_app()
