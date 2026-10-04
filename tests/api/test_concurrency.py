"""Concurrency: DNA, capabilities, Entrez lock."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor

from fastapi.testclient import TestClient

from helixscope_core.remote import NCBIClient


def test_concurrent_dna_and_capabilities(client: TestClient) -> None:
    def dna() -> None:
        response = client.post(
            "/api/v1/dna/analyze",
            json={"sequence": "ATGCATGCATGC", "include_explanation": False},
        )
        assert response.status_code == 200
        assert response.json()["result"]["gc_content"] == 50.0

    def caps() -> None:
        response = client.get("/api/v1/system/capabilities")
        assert response.status_code == 200
        assert "tools" in response.json()["result"]

    with ThreadPoolExecutor(max_workers=8) as pool:
        futs = [pool.submit(dna) for _ in range(6)] + [pool.submit(caps) for _ in range(6)]
        for fut in futs:
            fut.result()


def test_entrez_lock_serializes_two_emails(monkeypatch) -> None:
    seen: list[str] = []

    def fake_fetch(accession, email, db="nucleotide", api_key=None):
        from Bio import Entrez
        from modules.ncbi_fetch import _configure_entrez

        _configure_entrez(email, api_key)
        time.sleep(0.08)
        seen.append(str(getattr(Entrez, "email", "")))
        return {"accession": accession, "query": accession}

    monkeypatch.setattr("helixscope_core.remote.fetch_by_accession", fake_fetch)

    def worker(email: str, accession: str) -> str:
        return NCBIClient(email).fetch_by_accession(accession)["accession"]

    with ThreadPoolExecutor(max_workers=2) as pool:
        a = pool.submit(worker, "alpha@example.com", "NM_000001")
        b = pool.submit(worker, "beta@example.com", "NM_000002")
        assert {a.result(), b.result()} == {"NM_000001", "NM_000002"}
    assert set(seen) == {"alpha@example.com", "beta@example.com"}


def test_entrez_two_clients_no_email_leak_in_api(client: TestClient, monkeypatch) -> None:
    def fake_fetch(accession, email, db="nucleotide", api_key=None):
        time.sleep(0.02)
        return {"accession": accession, "description": "fixture", "email": email, "sequence": "ATGC"}

    monkeypatch.setattr("helixscope_core.remote.fetch_by_accession", fake_fetch)

    def one(email: str, acc: str) -> dict:
        return client.post(
            "/api/v1/ncbi/fetch",
            json={"accession": acc, "email": email, "db": "nucleotide"},
        )

    with ThreadPoolExecutor(max_workers=2) as pool:
        fa = pool.submit(one, "alpha@example.com", "NM_000001")
        fb = pool.submit(one, "beta@example.com", "NM_000002")
        ra, rb = fa.result(), fb.result()
    assert ra.status_code == 200, ra.text
    assert rb.status_code == 200, rb.text
    dumped = ra.text + rb.text
    assert "alpha@example.com" not in dumped
    assert "beta@example.com" not in dumped
    assert ra.json()["result"]["accession"] in {"NM_000001"}
    assert rb.json()["result"]["accession"] in {"NM_000002"}
