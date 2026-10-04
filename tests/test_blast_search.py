"""Testes offline do cliente NCBI BLAST Common URL API.

Fixtures XML sao copias de esquema BlastOutput rotuladas como fixtures, nao
respostas ao vivo. Nenhum teste inventa E-value, bit score ou accession.
"""

from __future__ import annotations

import io
from email.message import EmailMessage
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request

import pytest

from modules import blast_search, scientific_checks

FIXTURES = Path(__file__).parent / "fixtures"
QUERY = "ACGTACGTACGTACGTACGT"
EMAIL = "user@example.com"


def _load_fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


class _FakeHandle:
    def __init__(self, body: str) -> None:
        self._body = body.encode("utf-8")

    def read(self) -> bytes:
        return self._body

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _dispatcher(mapping):
    def fake(request: Request, timeout=None):
        payload = request.data.decode("utf-8") if request.data else ""
        assert str(request.full_url) == blast_search.BLAST_ENDPOINT
        if "CMD=Put" in payload:
            return _FakeHandle(mapping["put"])
        if "FORMAT_OBJECT=SearchInfo" in payload or "FORMAT_OBJECT=SearchInfo" in payload.replace(
            "%3D", "="
        ):
            return _FakeHandle(mapping["status"])
        if "FORMAT_TYPE=XML" in payload:
            return _FakeHandle(mapping["xml"])
        raise AssertionError(f"Unexpected BLAST payload: {payload}")

    return fake


def test_ncbi_remote_blast_is_independent_of_local_plus():
    remote = blast_search.ncbi_blast_availability()
    assert remote["available"] is True
    assert remote["endpoint"] == blast_search.BLAST_ENDPOINT
    assert remote["backend"] == "NCBI BLAST Common URL API"
    local = blast_search.local_blast_availability()
    assert isinstance(local, dict)
    assert "executables_detected" in local
    assert isinstance(local.get("database"), dict)
    assert remote["available"] is True
    ready = blast_search.deployment_readiness()
    assert ready["authentication"] is False
    assert ready["local_blast_plus"] == bool(local.get("available"))
    assert blast_search.BLAST_ENDPOINT in ready["endpoint_allowlist"]
    if local.get("available"):
        assert local.get("remote_unchanged") is True
        assert local.get("backend") == "NCBI BLAST+"


def test_molecule_program_mismatch_is_invalid_input():
    with pytest.raises(blast_search.BlastError) as protein_to_blastn:
        blast_search.prepare_query("MKTFFV", "blastn")
    assert protein_to_blastn.value.category == "INVALID_INPUT"
    with pytest.raises(blast_search.BlastError) as dna_to_blastp:
        blast_search.prepare_query(QUERY, "blastp")
    assert dna_to_blastp.value.category == "INVALID_INPUT"
    rna = blast_search.prepare_query("ACGUACGUACGUACGUACGU", "blastn")
    assert rna["rna_to_dna"] is True
    assert "U" not in rna["sequence"]
    assert rna["molecule"] == "RNA"


def test_query_length_resource_limit():
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.prepare_query("A" * (blast_search.MAX_QUERY_LENGTH + 1), "blastn")
    assert exc.value.category == "RESOURCE_LIMIT"


def test_invalid_database_and_program():
    with pytest.raises(blast_search.BlastError) as bad_prog:
        blast_search.databases_for_program("blastz")
    assert bad_prog.value.category == "INVALID_INPUT"
    with pytest.raises(blast_search.BlastError) as bad_db:
        blast_search.submit_search(
            QUERY,
            EMAIL,
            program="blastn",
            database="made_up_db",
            urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
        )
    assert bad_db.value.category == "INVALID_DATABASE"


def test_parse_put_and_status_pages():
    rid, rtoe = blast_search.parse_put_response("RID = ABCDEFGH\nRTOE = 22\n")
    assert rid == "ABCDEFGH"
    assert rtoe == 22
    with pytest.raises(blast_search.BlastError) as missing:
        blast_search.parse_put_response("no identifier here")
    assert missing.value.category == "JOB_FAILED"
    assert blast_search.parse_status_response("Status=WAITING\n") == "WAITING"
    assert blast_search.parse_status_response("Status=READY\n") == "READY"
    assert blast_search.parse_status_response(_load_fixture("blast_output_hit.xml")) == "READY"


def test_parse_valid_hit_fixture():
    parsed = blast_search.parse_blast_xml(_load_fixture("blast_output_hit.xml"))
    assert parsed["program"] == "blastn"
    assert parsed["blast_version"] == "BLASTN 2.16.0+"
    assert len(parsed["hits"]) == 1
    hit = parsed["hits"][0]
    assert hit["accession"] == "NM_000000"
    assert hit["evalue"] == pytest.approx(0.00012)
    assert hit["bit_score"] == pytest.approx(40.1)
    assert hit["identities"] == 20
    assert hit["alignment_length"] == 20
    assert hit["identity_pct"] == 100.0
    assert hit["organism"] == "Homo sapiens"
    assert scientific_checks.blast_hsp_is_valid(hit, query_len=20, hit_len=100)


def test_parse_no_hits_fixture_is_not_an_error():
    parsed = blast_search.parse_blast_xml(_load_fixture("blast_output_no_hits.xml"))
    assert parsed["hits"] == []
    assert parsed["skipped_hsps"] == 0


def test_malformed_xml_is_parsing_error():
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.parse_blast_xml("<BlastOutput><not-closed>")
    assert exc.value.category == "PARSING_ERROR"
    with pytest.raises(blast_search.BlastError) as html:
        blast_search.parse_blast_xml("this is not BlastOutput XML")
    assert html.value.category == "PARSING_ERROR"


def test_bad_identity_hsp_is_omitted_not_fabricated():
    parsed = blast_search.parse_blast_xml(_load_fixture("blast_output_bad_identity.xml"))
    assert parsed["hits"] == []
    assert parsed["skipped_hsps"] >= 1


def test_retrieve_all_invalid_hsps_is_parsing_error_not_no_hits():
    job = {
        "rid": "ABCDEFGH",
        "status": "READY",
        "program": "blastn",
        "database": "core_nt",
        "expect": 10.0,
        "hitlist_size": 20,
        "query_hash": "abc",
        "query_length": 20,
    }
    fake = _dispatcher(
        {
            "put": "RID = ABCDEFGH\nRTOE = 10\n",
            "status": "Status=READY\n",
            "xml": _load_fixture("blast_output_bad_identity.xml"),
        }
    )
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.retrieve_search(job, EMAIL, urlopen_fn=fake)
    assert exc.value.category == "PARSING_ERROR"
    assert exc.value.category != "NO_HITS"


def test_submit_poll_retrieve_valid_hit():
    fake = _dispatcher(
        {
            "put": "RID = ABCDEFGH\nRTOE = 10\n",
            "status": "Status=READY\n",
            "xml": _load_fixture("blast_output_hit.xml"),
        }
    )
    job = blast_search.submit_search(
        QUERY, EMAIL, program="blastn", database="core_nt", urlopen_fn=fake
    )
    assert job["rid"] == "ABCDEFGH"
    assert job["status"] == "WAITING"
    job["last_poll_monotonic"] = None
    job["submitted_monotonic"] = 0
    polled = blast_search.poll_search(
        job, EMAIL, urlopen_fn=fake, now=blast_search.MIN_CONTACT_INTERVAL_S + 1
    )
    assert polled["status"] == "READY"
    result = blast_search.retrieve_search(polled, EMAIL, urlopen_fn=fake)
    assert result["status"] == "READY"
    assert result["n_hits"] == 1
    assert result["hits"][0]["evalue"] == pytest.approx(0.00012)
    assert result["cache_status"] == "live"
    assert result["blast_version"] == "BLASTN 2.16.0+"
    csv_text = blast_search.export_hits_table(result)
    assert "NM_000000" in csv_text
    assert "0.00012" in csv_text
    cached = blast_search.mark_cached_result(result)
    assert cached["cache_status"] == "cached"
    assert cached["retrieved_at_utc"] == result["retrieved_at_utc"]
    assert cached["served_from_cache_at_utc"]


def test_retrieve_no_hits_status():
    fake = _dispatcher(
        {
            "put": "RID = ABCDEFGH\nRTOE = 10\n",
            "status": "Status=READY\n",
            "xml": _load_fixture("blast_output_no_hits.xml"),
        }
    )
    job = {
        "rid": "ABCDEFGH",
        "status": "READY",
        "program": "blastn",
        "database": "core_nt",
        "expect": 10.0,
        "hitlist_size": 20,
        "query_hash": "abc",
        "query_length": 20,
        "submitted_at_utc": "2026-01-01T00:00:00Z",
    }
    result = blast_search.retrieve_search(job, EMAIL, urlopen_fn=fake)
    assert result["status"] == "NO_HITS"
    assert result["hits"] == []
    assert result["n_hits"] == 0


def test_http_429_is_rate_limited_not_no_hits():
    def boom(request, timeout=None):
        raise HTTPError(
            blast_search.BLAST_ENDPOINT,
            429,
            "Too Many Requests",
            hdrs=EmailMessage(),
            fp=io.BytesIO(b""),
        )

    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.submit_search(QUERY, EMAIL, urlopen_fn=boom)
    assert exc.value.category == "RATE_LIMITED"
    assert exc.value.category != "NO_HITS"


def test_timeout_is_not_no_hits():
    def boom(request, timeout=None):
        raise TimeoutError("timed out")

    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.submit_search(QUERY, EMAIL, urlopen_fn=boom)
    assert exc.value.category == "TIMEOUT"


def test_service_unavailable_5xx():
    def boom(request, timeout=None):
        raise HTTPError(
            blast_search.BLAST_ENDPOINT,
            503,
            "Service Unavailable",
            hdrs=EmailMessage(),
            fp=io.BytesIO(b""),
        )

    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.submit_search(QUERY, EMAIL, urlopen_fn=boom)
    assert exc.value.category == "SERVICE_UNAVAILABLE"


def test_poll_too_soon_is_rate_limited():
    job = {
        "rid": "ABCDEFGH",
        "status": "WAITING",
        "submitted_monotonic": 1000.0,
        "last_poll_monotonic": None,
        "rtoe_s": 10.0,
    }
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.poll_search(
            job, EMAIL, urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(
                AssertionError("should not contact NCBI")
            ), now=1001.0
        )
    assert exc.value.category == "RATE_LIMITED"


def test_poll_timeout_after_max_wait():
    job = {
        "rid": "ABCDEFGH",
        "status": "WAITING",
        "submitted_monotonic": 0.0,
        "last_poll_monotonic": 0.0,
        "rtoe_s": 10.0,
    }
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.poll_search(
            job,
            EMAIL,
            urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("no net")),
            now=blast_search.MAX_JOB_WAIT_S + 1,
        )
    assert exc.value.category == "TIMEOUT"


def test_invalid_rid_rejected():
    job = {"rid": "https://evil.example/blast", "status": "WAITING", "submitted_monotonic": 0}
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.poll_search(job, EMAIL, now=10_000)
    assert exc.value.category == "INVALID_INPUT"


def test_xml_size_limit(monkeypatch):
    monkeypatch.setattr(blast_search, "MAX_XML_BYTES", 64)
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.parse_blast_xml("<BlastOutput>" + "A" * 200 + "</BlastOutput>")
    assert exc.value.category == "RESOURCE_LIMIT"


def test_paginate_large_hit_list():
    hits = [{"accession": f"ACC{i}", "evalue": float(i)} for i in range(25)]
    page = blast_search.paginate_hits(hits, 2, page_size=10)
    assert page["n_hits"] == 25
    assert page["page"] == 2
    assert page["n_pages"] == 3
    assert len(page["items"]) == 10
    assert page["items"][0]["accession"] == "ACC10"


def test_cache_key_changes_with_parameters():
    key_a = blast_search.cache_key(
        query_hash="abc", program="blastn", database="core_nt", expect=10.0, hitlist_size=20
    )
    key_b = blast_search.cache_key(
        query_hash="abc", program="blastn", database="nt", expect=10.0, hitlist_size=20
    )
    assert key_a != key_b


def test_job_failed_status_is_not_no_hits():
    fake = _dispatcher(
        {
            "put": "RID = ABCDEFGH\nRTOE = 10\n",
            "status": "Status=FAILED\n",
            "xml": "",
        }
    )
    job = {
        "rid": "ABCDEFGH",
        "status": "WAITING",
        "submitted_monotonic": 0.0,
        "last_poll_monotonic": None,
        "rtoe_s": 10.0,
    }
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.poll_search(job, EMAIL, urlopen_fn=fake, now=30.0)
    assert exc.value.category == "JOB_FAILED"


def test_exit_code_zero_is_success_not_failure() -> None:
    class _Ok:
        returncode = 0

    class _Fail:
        returncode = 1

    class _Missing:
        pass

    assert blast_search._subprocess_exited_nonzero(_Ok()) is False
    assert blast_search._subprocess_exited_nonzero(_Fail()) is True
    assert blast_search._subprocess_exited_nonzero(_Missing()) is True

