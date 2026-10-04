"""BLAST facade. Local BLAST+ is not NCBI nt/nr. Remote BLAST is a separate backend."""

from __future__ import annotations

from modules.blast_search import (
    local_blast_availability,
    ncbi_blast_availability,
    poll_search,
    retrieve_search,
    run_local_blast,
    submit_search,
)

__all__ = (
    "local_blast_availability",
    "ncbi_blast_availability",
    "poll_search",
    "retrieve_search",
    "run_local_blast",
    "submit_search",
)
