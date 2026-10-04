"""Cliente RCSB Alignment API: validacao, SSRF, 429 e parse de ticket."""

from __future__ import annotations

import io
import json
from email.message import EmailMessage
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request

import pytest

from modules import rcsb_alignment

FIXTURES = Path(__file__).parent / "fixtures"


class _FakeHandle:
    def __init__(self, body: str | bytes, url: str) -> None:
        self._body = body.encode("utf-8") if isinstance(body, str) else body
        self._url = url

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _http_error(url: str, code: int, retry_after: str | None = None) -> HTTPError:
    headers = EmailMessage()
    if retry_after is not None:
        headers["Retry-After"] = retry_after
    return HTTPError(url, code, "error", headers, io.BytesIO(b""))


def test_supported_methods_are_documented_api_names_only():
    names = {row["api_name"] for row in rcsb_alignment.supported_methods()}
    assert names == {
        "fatcat-rigid",
        "fatcat-flexible",
        "tm-align",
        "ce",
        "ce-cp",
    }
    assert "us-align" not in names
    assert rcsb_alignment.method_kind("tm-align") == "rigid"
    assert rcsb_alignment.method_kind("fatcat-flexible") == "flexible"
    with pytest.raises(rcsb_alignment.AlignmentApiError) as exc:
        rcsb_alignment.method_kind("made-up-aligner")
    assert exc.value.category == "INVALID_INPUT"


def test_entry_ids_accept_pdb_af_hyphen_and_csm():
    assert rcsb_alignment.validate_entry_id("8hsk") == "8HSK"
    assert rcsb_alignment.validate_entry_id("AF-P01308-F1") == "AF_AFP01308F1"
    assert rcsb_alignment.rcsb_csm_id_from_uniprot("B3EWR1") == "AF_AFB3EWR1F1"
    assert rcsb_alignment.parse_rcsb_csm_uniprot("AF_AFP01308F1") == ("P01308", 1)
    with pytest.raises(rcsb_alignment.AlignmentApiError):
        rcsb_alignment.validate_entry_id("https://evil.example/1CRN.cif")
    with pytest.raises(rcsb_alignment.AlignmentApiError):
        rcsb_alignment.validate_entry_id("BRCA1")
    with pytest.raises(rcsb_alignment.AlignmentApiError):
        rcsb_alignment.validate_entry_id("file:///etc/passwd")


def test_url_allowlist_blocks_ssrf_schemes_and_hosts():
    assert rcsb_alignment.url_is_allowed(
        "https://alignment.rcsb.org/api/v1/structures/submit"
    )
    assert rcsb_alignment.url_is_allowed(
        "https://alignment.rcsb.org/api/v1/structures/results?uuid=x"
    )
    for url in (
        "http://alignment.rcsb.org/api/v1/structures/submit",
        "https://evil.example/api/v1/structures/submit",
        "https://alignment.rcsb.org.evil.tld/api/v1/structures/submit",
        "https://alignment.rcsb.org/admin",
        "file:///c:/windows/win.ini",
        "gopher://alignment.rcsb.org/api/v1/structures/submit",
        "ftp://alignment.rcsb.org/api/v1/structures/submit",
        "https://127.0.0.1/api/v1/structures/submit",
        "https://localhost/api/v1/structures/submit",
    ):
        assert rcsb_alignment.url_is_allowed(url) is False, url
    with pytest.raises(rcsb_alignment.AlignmentApiError) as exc:
        rcsb_alignment._request_bytes("https://evil.example/x", method="GET")
    assert exc.value.category == "INVALID_INPUT"


def test_ticket_must_be_uuid():
    with pytest.raises(rcsb_alignment.AlignmentApiError) as exc:
        rcsb_alignment.fetch_results("../etc/passwd", poll=False)
    assert exc.value.category == "INVALID_INPUT"
    with pytest.raises(rcsb_alignment.AlignmentApiError):
        rcsb_alignment.fetch_results("not-a-uuid", poll=False)


def test_submit_and_fetch_with_injected_http():
    complete = (FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8")
    ticket = "095be615-a8ad-4c33-8e9c-c7612fbf6c9f"

    def opener(request: Request, timeout: float = 0) -> _FakeHandle:
        url = request.full_url
        assert rcsb_alignment.url_is_allowed(url)
        if request.get_method() == "POST":
            body = request.data or b""
            assert b"query=" in body
            return _FakeHandle(json.dumps({"ticket": ticket}), url)
        return _FakeHandle(complete, url)

    wrapped = rcsb_alignment.align_pairwise(
        reference_entry="8HSK",
        target_entry="8HSF",
        method="fatcat-rigid",
        urlopen_fn=opener,
        sleep_fn=lambda _s: None,
    )
    assert wrapped["ticket"] == ticket
    assert wrapped["engine_location"] == "remote"
    assert wrapped["raw"]["info"]["status"] == "COMPLETE"


def test_http_429_retries_then_succeeds(monkeypatch):
    calls = {"n": 0}
    ticket = "aaaaaaaa-bbbb-4ccc-8ddd-eeeeeeeeeeee"

    def opener(request: Request, timeout: float = 0) -> _FakeHandle:
        calls["n"] += 1
        if calls["n"] == 1:
            raise _http_error(request.full_url, 429, retry_after="0")
        return _FakeHandle(json.dumps({"ticket": ticket}), request.full_url)

    monkeypatch.setattr(rcsb_alignment.time, "sleep", lambda _s: None)
    submitted = rcsb_alignment.submit_pairwise(
        reference_entry="8HSK",
        target_entry="8HSF",
        urlopen_fn=opener,
    )
    assert submitted["ticket"] == ticket
    assert calls["n"] == 2


def test_fetch_retries_404_until_complete():
    complete = (FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8")
    ticket = "095be615-a8ad-4c33-8e9c-c7612fbf6c9f"
    calls = {"n": 0}

    def opener(request: Request, timeout: float = 0) -> _FakeHandle:
        calls["n"] += 1
        if calls["n"] == 1:
            raise _http_error(request.full_url, 404)
        return _FakeHandle(complete, request.full_url)

    payload = rcsb_alignment.fetch_results(
        ticket, urlopen_fn=opener, sleep_fn=lambda _s: None, poll=True
    )
    assert payload["info"]["status"] == "COMPLETE"
    assert calls["n"] == 2
    with pytest.raises(rcsb_alignment.AlignmentApiError) as exc:
        rcsb_alignment.build_pairwise_query(
            reference_entry="1CRN",
            target_entry="1CRN",
            reference_chain="A",
            target_chain="A",
            method="tm-align",
        )
    assert exc.value.category == "INVALID_INPUT"


def test_network_disclosure_does_not_claim_local_engine():
    note = rcsb_alignment.network_disclosure()
    assert "remote" in note["engine_location"]
    assert "arbitrary user urls" in " ".join(note["not_sent"]).lower()
