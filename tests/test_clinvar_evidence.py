"""Testes da camada de evidencia clinica ClinVar. Offline, com rede injetada.

Os payloads sao recortes reduzidos de respostas reais de
`esummary.fcgi?db=clinvar&retmode=json` para BRCA1, preservando os nomes de campo
originais do NCBI.
"""

from __future__ import annotations

import io
import json
from email.message import EmailMessage
from urllib.error import HTTPError, URLError

import pytest

from modules import clinvar_evidence

EMAIL = "helixscope@example.org"


class _FakeHandle:
    def __init__(self, body, url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/x") -> None:
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


def _summary(
    uid: str,
    *,
    accession: str,
    title: str,
    description: str,
    review_status: str,
    start: str,
    stop: str = "",
    scv: tuple[str, ...] = ("SCV000001.1",),
) -> dict:
    return {
        "uid": uid,
        "accession": accession,
        "accession_version": f"{accession}.1",
        "title": title,
        "obj_type": "single nucleotide variant",
        "record_status": "current",
        "protein_change": "M658I",
        "molecular_consequence_list": ["missense variant"],
        "genes": [{"symbol": "BRCA1", "GeneID": "672"}],
        "germline_classification": {
            "description": description,
            "last_evaluated": "2026/01/31 00:00",
            "review_status": review_status,
            "fda_recognized_database": "",
            "trait_set": [
                {
                    "trait_name": "Hereditary cancer-predisposing syndrome",
                    "trait_xrefs": [{"db_source": "MedGen", "db_id": "C0027672"}],
                }
            ],
        },
        "clinical_impact_classification": {
            "description": "",
            "last_evaluated": "1/01/01 00:00",
            "review_status": "",
            "trait_set": [],
        },
        "oncogenicity_classification": {
            "description": "",
            "last_evaluated": "1/01/01 00:00",
            "review_status": "",
            "trait_set": [],
        },
        "supporting_submissions": {"scv": list(scv), "rcv": ["RCV000001"]},
        "variation_set": [
            {
                "measure_id": "69092",
                "variation_name": title,
                "cdna_change": "c.1974G>C",
                "variation_xrefs": [{"db_source": "dbSNP", "db_id": "55678461"}],
                "variation_loc": [
                    {
                        "status": "current",
                        "assembly_name": "GRCh38",
                        "chr": "17",
                        "band": "17q21.31",
                        "start": start,
                        "stop": stop or start,
                        "assembly_acc_ver": "GCF_000001405.38",
                        "ref": "",
                        "alt": "",
                    },
                    {
                        "status": "previous",
                        "assembly_name": "GRCh37",
                        "chr": "17",
                        "start": "41245574",
                        "stop": "41245574",
                        "assembly_acc_ver": "GCF_000001405.25",
                        "ref": "",
                        "alt": "",
                    },
                ],
            }
        ],
    }


CONFLICTING = _summary(
    "54425",
    accession="VCV000054425",
    title="NM_007294.4(BRCA1):c.1974G>C (p.Met658Ile)",
    description="Conflicting classifications of pathogenicity",
    review_status="criteria provided, conflicting classifications",
    start="43093557",
    scv=tuple(f"SCV{index:09d}.1" for index in range(19)),
)
BENIGN_NEIGHBOUR = _summary(
    "1168433",
    accession="VCV001168433",
    title="NM_007294.4(BRCA1):c.1974G>A (p.Met658Ile)",
    description="Benign",
    review_status="criteria provided, single submitter",
    start="43093557",
)
FAR_AWAY = _summary(
    "868986",
    accession="VCV000868986",
    title="NM_007294.4(BRCA1):c.5499G>C (p.Val1833=)",
    description="Uncertain significance",
    review_status="criteria provided, single submitter",
    start="43045750",
)


def _pipeline(uids: list[str], records: list[dict]):
    """Devolve um opener que responde esearch e depois esummary."""

    def opener(request, timeout=None):
        url = str(request.full_url)
        if "esearch" in url:
            body = {"esearchresult": {"count": str(len(uids)), "idlist": uids,
                                      "querytranslation": "translated"}}
        else:
            result = {"uids": uids}
            for record in records:
                result[record["uid"]] = record
            body = {"result": result}
        return _FakeHandle(json.dumps(body))

    return opener


def test_url_allowlist_blocks_other_hosts() -> None:
    assert clinvar_evidence.url_is_allowed(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar"
    )
    assert not clinvar_evidence.url_is_allowed(
        "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    )
    assert not clinvar_evidence.url_is_allowed("https://evil.example.org/entrez/eutils/x")
    assert not clinvar_evidence.url_is_allowed(
        "https://eutils.ncbi.nlm.nih.gov/robots.txt"
    )
    assert not clinvar_evidence.url_is_allowed("")


def test_network_disclosure_includes_clinical_disclaimer() -> None:
    disclosure = clinvar_evidence.network_disclosure()
    assert disclosure["remote"] is True
    assert "not a diagnosis" in disclosure["disclaimer"]
    assert "not a clinical recommendation" in disclosure["disclaimer"]


def test_coordinate_search_term_requires_assembly() -> None:
    with pytest.raises(clinvar_evidence.ClinVarError) as exc:
        clinvar_evidence.build_search_term(contig="17", position_1based=43093557)
    assert exc.value.category == "INVALID_INPUT"
    term = clinvar_evidence.build_search_term(
        assembly="GRCh38.p14", contig="17", position_1based=43093557
    )
    assert term == "17[chr] AND 43093557[chrpos]"


def test_search_term_prefers_the_most_specific_identifier() -> None:
    assert clinvar_evidence.build_search_term(accession="VCV000054425") == "VCV000054425"
    assert clinvar_evidence.build_search_term(rsid="rs80357382") == "rs80357382"
    assert clinvar_evidence.build_search_term(hgvs="NM_007294.4:c.1974G>C").startswith("NM_")
    with pytest.raises(clinvar_evidence.ClinVarError):
        clinvar_evidence.build_search_term()


def test_no_record_is_not_found_and_never_benign() -> None:
    result = clinvar_evidence.fetch_clinvar_records(
        term="17[chr] AND 1[chrpos]",
        email=EMAIL,
        urlopen_fn=_pipeline([], []),
    )
    assert result["status"] == "NOT_FOUND"
    assert result["records"] == []
    assert "not evidence that the variant is benign" in result["message"]
    assert "benign" not in result["message"].lower().replace(
        "not evidence that the variant is benign", ""
    )


def test_conflicting_classification_is_preserved_not_resolved() -> None:
    result = clinvar_evidence.fetch_clinvar_records(
        term="17[chr] AND 43093557[chrpos]",
        email=EMAIL,
        urlopen_fn=_pipeline(["54425"], [CONFLICTING]),
    )
    record = result["records"][0]
    assert record["has_conflict"] is True
    assert record["germline_descriptions"] == ["Conflicting classifications of pathogenicity"]
    assert "does not resolve it" in record["conflict_note"]
    assert record["submission_count"] == 19
    assert record["accession"] == "VCV000054425"
    assert record["evidence_status"] == "RETRIEVED"
    assert "Not a HelixScope" in record["role"]


def test_classification_is_attributed_to_the_submitters() -> None:
    parsed = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["54425"], "54425": CONFLICTING}}, uid="54425"
    )
    block = parsed["classifications"][0]
    assert block["classification_type"] == "germline"
    assert "NCBI ClinVar aggregates" in block["aggregated_by"]
    assert "HelixScope does not recompute" in block["aggregated_by"]
    assert block["conditions"][0]["trait_name"] == "Hereditary cancer-predisposing syndrome"
    assert parsed["disclaimer"] == clinvar_evidence.CLINICAL_DISCLAIMER


def test_empty_classification_blocks_are_dropped_not_shown_as_blank() -> None:
    parsed = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["54425"], "54425": CONFLICTING}}, uid="54425"
    )
    types = {block["classification_type"] for block in parsed["classifications"]}
    assert types == {"germline"}


def test_fuzzy_text_hit_at_another_position_is_not_evidence() -> None:
    result = clinvar_evidence.fetch_clinvar_records(
        term="NM_007294.4(BRCA1):c.1974G>C",
        email=EMAIL,
        urlopen_fn=_pipeline(["868986"], [FAR_AWAY]),
        assembly="GRCh38.p14",
        contig="17",
        position_1based=43093557,
    )
    record = result["records"][0]
    assert record["identity_match"] is False
    assert result["identity_matched_count"] == 0
    assert "related record" in record["identity_match_reason"]


def test_multiple_variants_at_one_position_are_not_merged() -> None:
    result = clinvar_evidence.fetch_clinvar_records(
        term="17[chr] AND 43093557[chrpos]",
        email=EMAIL,
        urlopen_fn=_pipeline(["54425", "1168433"], [CONFLICTING, BENIGN_NEIGHBOUR]),
        assembly="GRCh38.p14",
        contig="17",
        position_1based=43093557,
        ref="C",
        alt="G",
    )
    assert result["identity_matched_count"] == 2
    assert "must not be merged" in result["identity_note"]
    descriptions = sorted(
        item for record in result["records"] for item in record["germline_descriptions"]
    )
    assert descriptions == ["Benign", "Conflicting classifications of pathogenicity"]
    # Alleles are not published in the summary, so allele identity stays unverified.
    assert result["allele_unverified_count"] == 2
    assert result["allele_verified_count"] == 0


def test_allele_match_is_computed_when_clinvar_publishes_alleles() -> None:
    record = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    record["locations"][0]["ref"] = "C"
    record["locations"][0]["alt"] = "G"
    outcome = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38.p14", contig="17", position_1based=43093557, ref="C", alt="G"
    )
    assert outcome["identity_match"] is True
    assert outcome["allele_match"] is True
    mismatch = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38.p14", contig="17", position_1based=43093557, ref="C", alt="T"
    )
    assert mismatch["allele_match"] is False


def test_grch37_and_grch38_locations_are_not_confused() -> None:
    record = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    on_38 = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38.p14", contig="17", position_1based=43093557
    )
    assert on_38["identity_match"] is True
    on_37_position = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38.p14", contig="17", position_1based=41245574
    )
    assert on_37_position["identity_match"] is False
    on_37 = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh37", contig="17", position_1based=41245574
    )
    assert on_37["identity_match"] is True


def test_patch_level_difference_is_reported_not_hidden() -> None:
    record = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    outcome = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38.p14", contig="17", position_1based=43093557
    )
    assert "GCF_000001405.38" in outcome["reason"]
    assert "patch" in outcome["reason"].lower()


def test_unsupported_assembly_is_not_silently_compared() -> None:
    record = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    outcome = clinvar_evidence.record_matches_variant(
        record, assembly="T2T-CHM13v2.0", contig="17", position_1based=43093557
    )
    assert outcome["identity_match"] is False
    assert "not comparable" in outcome["reason"]


def test_chr_prefix_is_accepted_for_matching() -> None:
    record = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    outcome = clinvar_evidence.record_matches_variant(
        record, assembly="GRCh38", contig="chr17", position_1based=43093557
    )
    assert outcome["identity_match"] is True


def test_invalid_email_is_refused_before_any_request() -> None:
    def explode(request, timeout=None):
        raise AssertionError("network must not be touched")

    for email in ("", "not-an-email", "a@b", "a b@c.org"):
        with pytest.raises(clinvar_evidence.ClinVarError) as exc:
            clinvar_evidence.search_clinvar(term="x", email=email, urlopen_fn=explode)
        assert exc.value.category == "INVALID_INPUT"


def test_oversized_term_is_resource_limit() -> None:
    with pytest.raises(clinvar_evidence.ClinVarError) as exc:
        clinvar_evidence.search_clinvar(
            term="a" * (clinvar_evidence.MAX_TERM_LENGTH + 1), email=EMAIL
        )
    assert exc.value.category == "RESOURCE_LIMIT"


def test_malformed_payloads_are_parsing_errors() -> None:
    with pytest.raises(clinvar_evidence.ClinVarError) as exc:
        clinvar_evidence.parse_clinvar_summary({"nope": 1})
    assert exc.value.category == "PARSING_ERROR"
    with pytest.raises(clinvar_evidence.ClinVarError):
        clinvar_evidence.parse_clinvar_summary("a string")
    with pytest.raises(clinvar_evidence.ClinVarError) as missing:
        clinvar_evidence.parse_clinvar_summary({"result": {"uids": []}})
    assert missing.value.category == "NOT_FOUND"


def test_network_failures_keep_their_category(monkeypatch) -> None:
    monkeypatch.setattr(clinvar_evidence.time, "sleep", lambda *_a: None)

    def timeout_opener(request, timeout=None):
        raise TimeoutError("timed out")

    with pytest.raises(clinvar_evidence.ClinVarError) as timed:
        clinvar_evidence.search_clinvar(term="x", email=EMAIL, urlopen_fn=timeout_opener)
    assert timed.value.category == "TIMEOUT"
    assert "not an absence of clinical submissions" in str(timed.value)

    def net_opener(request, timeout=None):
        raise URLError("refused")

    with pytest.raises(clinvar_evidence.ClinVarError) as net:
        clinvar_evidence.search_clinvar(term="x", email=EMAIL, urlopen_fn=net_opener)
    assert net.value.category == "NETWORK_ERROR"

    def rate_opener(request, timeout=None):
        raise HTTPError(
            str(request.full_url), 429, "slow down", hdrs=EmailMessage(), fp=io.BytesIO(b"")
        )

    with pytest.raises(clinvar_evidence.ClinVarError) as limited:
        clinvar_evidence.search_clinvar(term="x", email=EMAIL, urlopen_fn=rate_opener)
    assert limited.value.category == "RATE_LIMITED"


def test_response_from_unexpected_host_is_refused() -> None:
    def opener(request, timeout=None):
        return _FakeHandle(
            json.dumps({"esearchresult": {"idlist": [], "count": "0"}}),
            "https://evil.example.org/entrez/eutils/x",
        )

    with pytest.raises(clinvar_evidence.ClinVarError) as exc:
        clinvar_evidence.search_clinvar(term="x", email=EMAIL, urlopen_fn=opener)
    assert exc.value.category == "INVALID_INPUT"


def test_hostile_metadata_is_carried_as_text_not_executed() -> None:
    hostile = json.loads(json.dumps(CONFLICTING))
    hostile["title"] = "<script>alert('xss')</script>"
    hostile["genes"] = [{"symbol": "<img src=x onerror=alert(1)>", "GeneID": "672"}]
    parsed = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["54425"], "54425": hostile}}, uid="54425"
    )
    assert parsed["title"] == "<script>alert('xss')</script>"
    assert isinstance(parsed["genes"][0]["symbol"], str)
    assert parsed["accession"] == "VCV000054425"


def test_review_status_ordering_is_presentation_only() -> None:
    assert "practice guideline" in clinvar_evidence.REVIEW_STATUS_ORDER
    parsed = clinvar_evidence.parse_clinvar_summary(
        {"result": {"uids": ["1"], "1": CONFLICTING}}, uid="1"
    )
    ranks = [block["review_status_rank"] for block in parsed["classifications"]]
    assert ranks == sorted(ranks)
