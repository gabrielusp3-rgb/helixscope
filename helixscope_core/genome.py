"""Genome-store and Cas-OFFinder job facade. Not a public FastAPI import of modules."""

from __future__ import annotations

from modules.crispr_casoffinder import detect_cas_offinder
from modules.genome_jobs import (
    full_genome_search_preflight,
    load_job,
    load_result,
    load_status,
    request_cancel,
    submit_casoffinder_job,
)

__all__ = (
    "detect_cas_offinder",
    "full_genome_search_preflight",
    "load_job",
    "load_result",
    "load_status",
    "request_cancel",
    "submit_casoffinder_job",
)
