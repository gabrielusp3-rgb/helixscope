"""Remote scientific clients. Email and keys are constructor arguments, not widgets.

Entrez.email is a Biopython process global. All Entrez calls that go through
this module take ``_ENTREZ_LOCK`` so concurrent FastAPI requests cannot mix
contact emails.
"""

from __future__ import annotations

import threading
from contextlib import contextmanager
from typing import Any, Iterator, Optional

from modules import clinvar_evidence, ensembl_vep, ncbi_fetch, protein_domains, rcsb_alignment
from modules.ncbi_fetch import NCBIQueryError, api_key_from_environment, fetch_by_accession, search_records

__all__ = (
    "NCBIClient",
    "NCBIQueryError",
    "api_key_from_environment",
    "entrez_lock",
    "fetch_by_accession_locked",
    "search_records_locked",
)

_ENTREZ_LOCK = threading.Lock()


@contextmanager
def entrez_lock() -> Iterator[None]:
    """Serialize Biopython Entrez configuration and the following network call."""
    _ENTREZ_LOCK.acquire()
    try:
        yield
    finally:
        _ENTREZ_LOCK.release()


def fetch_by_accession_locked(
    accession: str,
    email: str,
    db: str = "nucleotide",
    api_key: Optional[str] = None,
) -> dict[str, Any]:
    """``fetch_by_accession`` under the Entrez lock."""
    with entrez_lock():
        return fetch_by_accession(accession, email, db=db, api_key=api_key)


def search_records_locked(
    term: str,
    email: str,
    db: str = "nucleotide",
    api_key: Optional[str] = None,
    retmax: int = 15,
    retstart: int = 0,
) -> dict[str, Any]:
    """``search_records`` under the Entrez lock."""
    with entrez_lock():
        return search_records(
            term,
            email,
            db=db,
            api_key=api_key,
            retmax=retmax,
            retstart=retstart,
        )


class NCBIClient:
    """Entrez client with explicit contact email.

    Each call holds the process lock for configure + request so two emails
    cannot overlap on ``Bio.Entrez.email``.
    """

    def __init__(self, email: str, api_key: Optional[str] = None) -> None:
        if not email or not str(email).strip():
            raise ValueError("NCBI Entrez requires a contact email.")
        self.email = str(email).strip()
        self.api_key = (str(api_key).strip() if api_key else None) or None

    def fetch_by_accession(self, accession: str, db: str = "nucleotide") -> dict[str, Any]:
        """Fetch one NCBI record. Same science as ``ncbi_fetch.fetch_by_accession``."""
        return fetch_by_accession_locked(
            accession, self.email, db=db, api_key=self.api_key
        )

    def search_records(
        self,
        term: str,
        db: str = "nucleotide",
        retmax: int = 15,
        retstart: int = 0,
    ) -> dict[str, Any]:
        """Entrez esearch + esummary under the Entrez lock."""
        return search_records_locked(
            term,
            self.email,
            db=db,
            api_key=self.api_key,
            retmax=retmax,
            retstart=retstart,
        )


def vep_client_module() -> Any:
    """Ensembl VEP module (fixture-tested). Not a web session object."""
    return ensembl_vep


def clinvar_client_module() -> Any:
    """ClinVar E-utilities module."""
    return clinvar_evidence


def interpro_uniprot_module() -> Any:
    """InterPro / UniProt domain client module."""
    return protein_domains


def rcsb_alignment_module() -> Any:
    """RCSB structure-alignment API module."""
    return rcsb_alignment
