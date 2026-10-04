"""Controlled live compatibility smokes for the scientific rebuild.

These tests are skipped unless HELIXSCOPE_LIVE=1. They must not run in the
default deterministic pytest gate. They contact official providers at most
once per provider and never download genomes.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pytest

LIVE = os.environ.get("HELIXSCOPE_LIVE", "").strip() == "1"
ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "SCIENTIFIC_REBUILD_LIVE.json"

pytestmark = pytest.mark.skipif(not LIVE, reason="Set HELIXSCOPE_LIVE=1 for controlled live smokes")


def _stamp(provider: str, operation: str, **fields: object) -> dict:
    row = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "provider": provider,
        "operation": operation,
        **fields,
    }
    existing: list = []
    if REPORT.exists():
        existing = json.loads(REPORT.read_text(encoding="utf-8"))
        if not isinstance(existing, list):
            existing = []
    existing.append(row)
    REPORT.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")
    return row


def test_live_uniprot_rest_entry() -> None:
    from modules.protein_domains import fetch_uniprot_protein

    entry = fetch_uniprot_protein(accession="P01308")
    assert entry.get("accession")
    _stamp(
        "UniProt REST",
        "uniprotkb/P01308",
        accession=entry.get("accession"),
        organism=entry.get("organism") or None,
        status=entry.get("evidence_status") or entry.get("status"),
        parser_ok=True,
    )


def test_live_ensembl_rest_info() -> None:
    from modules.ensembl_vep import VEP_HOST_GRCH38, _http_json, service_availability

    ping = service_availability()
    rest = _http_json(f"https://{VEP_HOST_GRCH38}/info/rest")
    version = None
    if isinstance(rest, dict):
        version = rest.get("release") or rest.get("version") or rest.get("server")
    _stamp(
        "Ensembl REST",
        "/info/rest",
        provider_version=version,
        ping_available=ping.get("available"),
        data_release=ping.get("release"),
        http_ok=True,
        parser_ok=isinstance(rest, dict),
    )
    assert ping.get("available") is True
    assert isinstance(rest, dict)


def test_live_ncbi_einfo_nucleotide() -> None:
    from Bio import Entrez

    Entrez.email = "helixscope@example.invalid"
    Entrez.tool = "helixscope"
    handle = Entrez.einfo(db="nucleotide")
    try:
        payload = Entrez.read(handle)
    finally:
        handle.close()
    db = payload.get("DbInfo") if isinstance(payload, dict) else None
    _stamp(
        "NCBI E-utilities",
        "einfo nucleotide",
        http_ok=True,
        parser_ok=db is not None,
        db_name=str((db or {}).get("DbName") if isinstance(db, dict) else ""),
    )
    assert db is not None
