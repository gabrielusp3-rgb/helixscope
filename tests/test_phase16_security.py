"""Red team da Fase 16: SSRF, injecao, traversal, XSS e limites de recurso.

Nenhum teste abre socket. Todos os clientes de rede recebem `urlopen_fn`
injetado e a maioria dos casos afirma que a requisicao e RECUSADA antes de
qualquer chamada, contando as invocacoes do fake.

Nenhum teste ataca um servico de terceiros: os alvos hostis sao hosts
ficticios e cargas locais.
"""

from __future__ import annotations

import io
import json
import re
from email.message import EmailMessage
from pathlib import Path
from urllib.error import HTTPError

import pytest

from modules import (
    clinvar_evidence,
    ensembl_vep,
    protein_domains,
    variant_core,
    variant_explorer,
)
from ui import components

REPO_ROOT = Path(__file__).resolve().parents[1]


class _CountingHandle:
    """Resposta HTTP falsa que registra cada abertura."""

    def __init__(self, body: str, calls: list, url: str) -> None:
        self._body = body.encode("utf-8")
        self._calls = calls
        self._url = url

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_CountingHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _counting_urlopen(calls: list, body: str = "[]", url: str = "https://rest.ensembl.org/vep/x"):
    def _open(request, timeout=None):  # noqa: ANN001 - assinatura de urlopen
        target = getattr(request, "full_url", str(request))
        calls.append(target)
        return _CountingHandle(body, calls, url)

    return _open


def _http_error(code: int, url: str, retry_after: str = "") -> HTTPError:
    headers = EmailMessage()
    if retry_after:
        headers["Retry-After"] = retry_after
    return HTTPError(url, code, "forced", headers, io.BytesIO(b"{}"))


@pytest.fixture(autouse=True)
def _no_backoff_sleep(monkeypatch):
    """Remove a espera exponencial para que o teste meca tentativas, nao tempo."""
    for module in (ensembl_vep, clinvar_evidence, protein_domains):
        monkeypatch.setattr(module.time, "sleep", lambda _seconds: None)


HOSTILE_URLS = (
    "http://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T",
    "https://rest.ensembl.org.attacker.example/vep/homo_sapiens/region/17:1-1:1/T",
    "https://attacker.example/vep/homo_sapiens/region/17:1-1:1/T",
    "https://rest.ensembl.org@attacker.example/vep/x",
    "https://127.0.0.1/vep/x",
    "https://localhost:8501/vep/x",
    "https://169.254.169.254/latest/meta-data/",
    "https://[::1]/vep/x",
    "file:///c:/windows/win.ini",
    "gopher://rest.ensembl.org/vep/x",
    "",
)


def test_vep_allowlist_refuses_hostile_urls():
    for url in HOSTILE_URLS:
        assert ensembl_vep.url_is_allowed(url) is False, url
    assert ensembl_vep.url_is_allowed(
        "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T?content-type=application/json"
    )
    assert ensembl_vep.url_is_allowed("https://grch37.rest.ensembl.org/info/data/")
    assert ensembl_vep.url_is_allowed("https://rest.ensembl.org/lookup/id/ENSG00000012048") is False


def test_clinvar_and_domain_allowlists_refuse_hostile_urls():
    for url in HOSTILE_URLS:
        assert clinvar_evidence.url_is_allowed(url) is False, url
        assert protein_domains.url_is_allowed(url) is False, url
    assert clinvar_evidence.url_is_allowed(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar"
    )
    assert clinvar_evidence.url_is_allowed("https://eutils.ncbi.nlm.nih.gov/robots.txt") is False
    assert protein_domains.url_is_allowed(
        "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398"
    )
    assert protein_domains.url_is_allowed("https://rest.uniprot.org/uniprotkb/P38398.json")
    assert protein_domains.url_is_allowed("https://www.ebi.ac.uk/proteins/api/features") is False


def test_non_allowlisted_url_is_refused_before_any_socket_is_opened():
    calls: list = []
    opener = _counting_urlopen(calls)
    with pytest.raises(ensembl_vep.VepError) as vep_error:
        ensembl_vep._http_json("https://attacker.example/vep/x", urlopen_fn=opener)
    assert vep_error.value.category == "INVALID_INPUT"
    with pytest.raises(clinvar_evidence.ClinVarError) as clinvar_error:
        clinvar_evidence._http_json("https://attacker.example/entrez/eutils/x", urlopen_fn=opener)
    assert clinvar_error.value.category == "INVALID_INPUT"
    with pytest.raises(protein_domains.DomainError) as domain_error:
        protein_domains._http_json("https://attacker.example/interpro/api/x", urlopen_fn=opener)
    assert domain_error.value.category == "INVALID_INPUT"
    assert calls == []


def test_redirect_handlers_refuse_off_allowlist_targets():
    with pytest.raises(ensembl_vep.VepError):
        ensembl_vep._AllowlistRedirect().redirect_request(
            None, None, 302, "found", {}, "https://attacker.example/vep/x"
        )
    with pytest.raises(clinvar_evidence.ClinVarError):
        clinvar_evidence._AllowlistRedirect().redirect_request(
            None, None, 302, "found", {}, "https://attacker.example/entrez/eutils/x"
        )
    with pytest.raises(protein_domains.DomainError):
        protein_domains._AllowlistRedirect().redirect_request(
            None, None, 302, "found", {}, "https://attacker.example/interpro/api/x"
        )


UNSAFE_CONTIGS = (
    "../../etc/passwd",
    "17/../18",
    "17?db=evil",
    "17#frag",
    "17&x=1",
    "17 43093557",
    "17%2f..",
    "17:1",
)


def test_vep_region_query_refuses_unsafe_contig_tokens():
    calls: list = []
    opener = _counting_urlopen(calls)
    for contig in UNSAFE_CONTIGS:
        with pytest.raises(ensembl_vep.VepError) as error:
            ensembl_vep.annotate_region(
                contig=contig, position_1based=43093557, ref="C", alt="G", urlopen_fn=opener
            )
        assert error.value.category == "INVALID_INPUT", contig
    assert calls == []


def test_vep_region_query_refuses_unchecked_alleles_and_oversized_events():
    calls: list = []
    opener = _counting_urlopen(calls)
    with pytest.raises(ensembl_vep.VepError) as ref_error:
        ensembl_vep.annotate_region(
            contig="17", position_1based=43093557, ref="CNX", alt="G", urlopen_fn=opener
        )
    assert ref_error.value.category == "INVALID_INPUT"
    with pytest.raises(ensembl_vep.VepError) as alt_error:
        ensembl_vep.annotate_region(
            contig="17", position_1based=43093557, ref="C", alt="G<script>", urlopen_fn=opener
        )
    assert alt_error.value.category == "INVALID_INPUT"
    with pytest.raises(ensembl_vep.VepError) as size_error:
        ensembl_vep.annotate_region(
            contig="17",
            position_1based=43093557,
            ref="A" * (ensembl_vep.MAX_ALLELE_LENGTH + 1),
            alt="G",
            urlopen_fn=opener,
        )
    assert size_error.value.category == "RESOURCE_LIMIT"
    assert calls == []


def test_vep_hgvs_and_identifier_endpoints_refuse_path_changing_input():
    calls: list = []
    opener = _counting_urlopen(calls)
    for notation in (
        "NC_000017.11:g.43093557C>G/../../info",
        "NC_000017.11:g.43093557C>G?x=1",
        "NC_000017.11:g.43093557C>G#f",
        "NC_000017.11:g.43093557C>G AND pathogenic",
        "A" * 300,
    ):
        with pytest.raises(ensembl_vep.VepError):
            ensembl_vep.annotate_hgvs(notation=notation, urlopen_fn=opener)
    for identifier in ("rs80357382; rm -rf /", "../../etc/passwd", "rs80357382/x", "rs" + "9" * 80):
        with pytest.raises(ensembl_vep.VepError):
            ensembl_vep.annotate_identifier(identifier=identifier, urlopen_fn=opener)
    assert calls == []


def test_vep_refuses_unknown_assembly_instead_of_guessing_a_host():
    for assembly in ("", "GRCm39", "T2T-CHM13v2.0", "hg17"):
        with pytest.raises(ensembl_vep.VepError) as error:
            ensembl_vep.host_for_assembly(assembly)
        assert error.value.category == "INVALID_INPUT"
    assert ensembl_vep.host_for_assembly("GRCh38.p14") == ensembl_vep.VEP_HOST_GRCH38
    assert ensembl_vep.host_for_assembly("GRCh37") == ensembl_vep.VEP_HOST_GRCH37


ENTREZ_INJECTIONS = (
    "17[chr] OR 13",
    "17) OR (1",
    "17 AND pathogenic[clinsig]",
    "17;13",
    "../../etc/passwd",
    "17'",
    '17"',
)


def test_clinvar_term_builder_refuses_entrez_operator_injection():
    for contig in ENTREZ_INJECTIONS:
        with pytest.raises(clinvar_evidence.ClinVarError) as error:
            clinvar_evidence.build_search_term(
                assembly="GRCh38", contig=contig, position_1based=43093557
            )
        assert error.value.category == "INVALID_INPUT", contig
    for accession in ("VCV000017677 OR pathogenic", "VCV; rm -rf /", "NOTANACC123", "VCV"):
        with pytest.raises(clinvar_evidence.ClinVarError):
            clinvar_evidence.build_search_term(accession=accession)
    for rsid in ("rs80357382 OR rs1", "rsX", "rs80357382[chr]"):
        with pytest.raises(clinvar_evidence.ClinVarError):
            clinvar_evidence.build_search_term(rsid=rsid)
    for hgvs in ("NC_000017.11:g.1C>G OR pathogenic[clinsig]", "NC_000017.11:g.1C>G (BRCA1)"):
        with pytest.raises(clinvar_evidence.ClinVarError):
            clinvar_evidence.build_search_term(hgvs=hgvs)
    assert (
        clinvar_evidence.build_search_term(
            assembly="GRCh38", contig="17", position_1based=43093557
        )
        == "17[chr] AND 43093557[chrpos]"
    )
    assert clinvar_evidence.build_search_term(accession="vcv000017677.2") == "VCV000017677.2"
    assert clinvar_evidence.build_search_term(rsid="RS80357382") == "rs80357382"


def test_clinvar_requires_a_valid_contact_email_and_never_leaks_it_into_the_path():
    calls: list = []
    opener = _counting_urlopen(calls)
    for email in ("", "not-an-email", "a@b", "a@b.c\nX-Injected: 1"):
        with pytest.raises(clinvar_evidence.ClinVarError) as error:
            clinvar_evidence.search_clinvar(term="17[chr]", email=email, urlopen_fn=opener)
        assert error.value.category == "INVALID_INPUT", email
    assert calls == []
    ok_calls: list = []
    clinvar_evidence.search_clinvar(
        term="17[chr] AND 43093557[chrpos]",
        email="scientist@example.org",
        urlopen_fn=_counting_urlopen(
            ok_calls,
            body=json.dumps({"esearchresult": {"idlist": [], "count": "0"}}),
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        ),
    )
    assert len(ok_calls) == 1
    assert ok_calls[0].startswith("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?")
    assert "scientist%40example.org" in ok_calls[0]
    assert "\n" not in ok_calls[0]


UNSAFE_ACCESSIONS = (
    "../../etc/passwd",
    "P38398/../../uniprotkb",
    "P38398?fields=all",
    "P38398 OR Q00000",
    "P38398#f",
    "'; DROP TABLE",
    "<script>alert(1)</script>",
    "P" * 200,
)


def test_uniprot_accession_validation_blocks_path_and_query_manipulation():
    calls: list = []
    opener = _counting_urlopen(calls)
    for accession in UNSAFE_ACCESSIONS:
        with pytest.raises(protein_domains.DomainError) as error:
            protein_domains.validate_uniprot_accession(accession)
        assert error.value.category == "INVALID_INPUT", accession
        with pytest.raises(protein_domains.DomainError):
            protein_domains.fetch_protein_domains(accession=accession, urlopen_fn=opener)
        with pytest.raises(protein_domains.DomainError):
            protein_domains.fetch_uniprot_protein(accession=accession, urlopen_fn=opener)
    assert calls == []
    assert protein_domains.validate_uniprot_accession(" p38398 ") == "P38398"


def test_interpro_pagination_pointing_outside_the_allowlist_is_refused():
    pages = [
        json.dumps(
            {
                "count": 2,
                "next": "https://attacker.example/interpro/api/entry/interpro/protein/uniprot/P38398",
                "results": [],
            }
        )
    ]

    def _open(request, timeout=None):  # noqa: ANN001
        return _CountingHandle(
            pages.pop(0),
            [],
            "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398",
        )

    with pytest.raises(protein_domains.DomainError) as error:
        protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=_open)
    assert error.value.category == "INVALID_INPUT"
    assert "allowlist" in str(error.value)


VARIANT_INJECTIONS = (
    "17 43093557 C G; rm -rf /",
    "17 43093557 C $(whoami)",
    "17 43093557 C `id`",
    "../../etc/passwd",
    "17 43093557 C G && curl attacker.example",
    "17 43093557 C G\nX-Injected: 1",
    "17 43093557 C <script>alert(1)</script>",
    "chr17:43093557 A>G",
    "17 0 C G",
    "17 -5 C G",
    "17 43093557 C C",
    "",
    "   ",
)


def test_variant_input_refuses_injection_traversal_and_degenerate_alleles():
    for text in VARIANT_INJECTIONS:
        with pytest.raises(variant_core.VariantInputError):
            variant_core.build_variant(text=text, assembly="GRCh38.p14")


def test_variant_input_refuses_oversized_payloads():
    with pytest.raises(variant_core.VariantInputError) as error:
        variant_core.build_variant(text="17 43093557 " + "A" * 20000 + " G", assembly="GRCh38.p14")
    assert error.value.category in {"INVALID_INPUT", "RESOURCE_LIMIT"}


def test_variant_identity_refuses_missing_assembly_and_unsafe_contig_override():
    with pytest.raises(variant_core.VariantInputError) as missing:
        variant_core.build_variant(text="17 43093557 C G", assembly="")
    assert missing.value.category == "INVALID_INPUT"
    for contig in ("../../etc/passwd", "17/18", "17 18", "17;18", "a" * 200):
        with pytest.raises(variant_core.VariantInputError):
            variant_core.build_variant(
                text="17 43093557 C G", assembly="GRCh38.p14", contig_override=contig
            )


def test_hostile_external_metadata_is_preserved_verbatim_but_escaped_on_render():
    payload = [
        {
            "id": "17_1_C/G",
            "assembly_name": "GRCh38",
            "seq_region_name": "17",
            "start": 1,
            "end": 1,
            "allele_string": "C/G",
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [
                {
                    "transcript_id": "ENST00000357654",
                    "gene_symbol": '<img src=x onerror="alert(1)">',
                    "biotype": "protein_coding",
                    "consequence_terms": ["missense_variant"],
                    "protein_start": 10,
                    "amino_acids": "S/R",
                    "codons": "agT/agA",
                    "hgvsp": "<script>alert(1)</script>",
                }
            ],
        }
    ]
    parsed = ensembl_vep.parse_vep_response(payload, host=ensembl_vep.VEP_HOST_GRCH38)
    row = parsed["transcript_consequences"][0]
    assert row["gene_symbol"] == '<img src=x onerror="alert(1)">'
    rendered = components.meta_grid(
        [("gene", row["gene_symbol"]), ("HGVSp", row["hgvsp"])]
    )
    assert "<img" not in rendered
    assert "<script>" not in rendered
    assert "&lt;img" in rendered
    inspector = components.inspector_panel({"VEP": [("gene", row["gene_symbol"])]})
    assert "<img" not in inspector
    badge = components.status_badge(row["gene_symbol"])
    assert "<img" not in badge


def test_hostile_clinvar_and_interpro_metadata_is_escaped_on_render():
    clinvar_payload = {
        "result": {
            "uids": ["12345"],
            "12345": {
                "uid": "12345",
                "accession": "VCV000012345",
                "accession_version": "VCV000012345.1",
                "title": '<script>alert("clinvar")</script>',
                "germline_classification": {
                    "description": "<b>Pathogenic</b>",
                    "review_status": "criteria provided, single submitter",
                    "trait_set": [{"trait_name": "<img src=x onerror=alert(1)>"}],
                },
            },
        }
    }
    record = clinvar_evidence.parse_clinvar_summary(clinvar_payload, uid="12345")
    assert record["title"] == '<script>alert("clinvar")</script>'
    rendered = components.meta_grid([("title", record["title"])])
    assert "<script>" not in rendered
    interpro_payload = {
        "count": 1,
        "results": [
            {
                "metadata": {
                    "accession": "IPR001841",
                    "name": "<svg onload=alert(1)>",
                    "type": "domain",
                    "source_database": "interpro",
                },
                "proteins": [
                    {
                        "protein_length": 1863,
                        "entry_protein_locations": [
                            {"fragments": [{"start": 24, "end": 65}]}
                        ],
                    }
                ],
            }
        ],
    }
    entries = protein_domains.parse_interpro_entries(interpro_payload, accession="P38398")
    assert entries[0]["name"] == "<svg onload=alert(1)>"
    domain_render = components.meta_grid([("domain", entries[0]["name"])])
    assert "<svg" not in domain_render
    assert "&lt;svg" in domain_render


def test_oversized_responses_are_truncated_and_refused():
    oversized = "[" + "0," * (ensembl_vep.MAX_RESPONSE_BYTES // 2) + "0]"
    with pytest.raises(ensembl_vep.VepError) as error:
        ensembl_vep._http_json(
            "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T",
            urlopen_fn=_counting_urlopen([], body=oversized),
        )
    assert error.value.category in {"RESOURCE_LIMIT", "PARSING_ERROR"}


def test_malformed_json_from_every_service_is_reported_as_parsing_error():
    with pytest.raises(ensembl_vep.VepError) as vep_error:
        ensembl_vep._http_json(
            "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T",
            urlopen_fn=_counting_urlopen([], body="{not json"),
        )
    assert vep_error.value.category == "PARSING_ERROR"
    clinvar_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    with pytest.raises(clinvar_evidence.ClinVarError) as clinvar_error:
        clinvar_evidence._http_json(
            clinvar_url,
            urlopen_fn=_counting_urlopen([], body="<html>not json</html>", url=clinvar_url),
        )
    assert clinvar_error.value.category == "PARSING_ERROR"
    interpro_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398"
    with pytest.raises(protein_domains.DomainError) as domain_error:
        protein_domains._http_json(
            interpro_url,
            urlopen_fn=_counting_urlopen([], body="\x00\x01\x02", url=interpro_url),
        )
    assert domain_error.value.category == "PARSING_ERROR"


def test_a_response_arriving_from_a_non_allowlisted_url_is_refused():
    """Defesa contra redirect fora do allowlist que escape ao handler."""
    for module, error_type, url in (
        (
            ensembl_vep,
            ensembl_vep.VepError,
            "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T",
        ),
        (
            clinvar_evidence,
            clinvar_evidence.ClinVarError,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        ),
        (
            protein_domains,
            protein_domains.DomainError,
            "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398",
        ),
    ):
        with pytest.raises(error_type) as error:
            module._http_json(
                url,
                urlopen_fn=_counting_urlopen(
                    [], body="{}", url="https://attacker.example/collected"
                ),
            )
        assert error.value.category == "INVALID_INPUT"


def test_retries_are_bounded_so_a_failing_service_cannot_be_hammered():
    for module, error_type, url in (
        (
            ensembl_vep,
            ensembl_vep.VepError,
            "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T",
        ),
        (
            clinvar_evidence,
            clinvar_evidence.ClinVarError,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        ),
        (
            protein_domains,
            protein_domains.DomainError,
            "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398",
        ),
    ):
        calls: list = []

        def _open(request, timeout=None, _url=url, _calls=calls):  # noqa: ANN001
            _calls.append(getattr(request, "full_url", str(request)))
            raise _http_error(503, _url)

        with pytest.raises(error_type) as error:
            module._http_json(url, urlopen_fn=_open)
        assert error.value.category == "SERVICE_UNAVAILABLE"
        assert len(calls) == module.MAX_RETRIES


def test_rate_limit_response_is_labelled_and_not_retried_forever():
    calls: list = []

    def _open(request, timeout=None):  # noqa: ANN001
        calls.append(1)
        raise _http_error(429, "https://rest.ensembl.org/vep/x", retry_after="1")

    with pytest.raises(ensembl_vep.VepError) as error:
        ensembl_vep._http_json(
            "https://rest.ensembl.org/vep/homo_sapiens/region/17:1-1:1/T", urlopen_fn=_open
        )
    assert error.value.category == "RATE_LIMITED"
    assert error.value.retry_after_s == 1.0
    assert len(calls) == ensembl_vep.MAX_RETRIES


def test_variant_cache_key_prevents_cross_identity_and_cross_version_poisoning():
    base = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    other_assembly = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh37")
    other_allele = variant_core.build_variant(text="17 43093557 C T", assembly="GRCh38.p14")
    keys = {
        variant_explorer.variant_cache_key(
            variant=base, layer="consequences", tool="Ensembl VEP", tool_version="116"
        ),
        variant_explorer.variant_cache_key(
            variant=other_assembly, layer="consequences", tool="Ensembl VEP", tool_version="116"
        ),
        variant_explorer.variant_cache_key(
            variant=other_allele, layer="consequences", tool="Ensembl VEP", tool_version="116"
        ),
        variant_explorer.variant_cache_key(
            variant=base, layer="clinical_evidence", tool="Ensembl VEP", tool_version="116"
        ),
        variant_explorer.variant_cache_key(
            variant=base, layer="consequences", tool="Ensembl VEP", tool_version="115"
        ),
        variant_explorer.variant_cache_key(
            variant=base,
            layer="consequences",
            tool="Ensembl VEP",
            tool_version="116",
            parameters={"domains": "1"},
        ),
    }
    assert len(keys) == 6
    assert base["identity_hash"] != other_assembly["identity_hash"]


def test_network_disclosure_names_every_remote_host_before_any_query():
    vep = ensembl_vep.network_disclosure("GRCh38.p14")
    assert vep["remote"] is True
    assert vep["host"] == ensembl_vep.VEP_HOST_GRCH38
    assert ensembl_vep.network_disclosure("GRCh37")["host"] == ensembl_vep.VEP_HOST_GRCH37
    assert vep["data_sent"]
    clinvar = clinvar_evidence.network_disclosure()
    assert clinvar["host"] == clinvar_evidence.EUTILS_HOST
    assert clinvar["disclaimer"] == clinvar_evidence.CLINICAL_DISCLAIMER
    domains = protein_domains.network_disclosure()
    hosts = {str(service["host"]) for service in domains["services"]}
    assert hosts == {protein_domains.INTERPRO_HOST, protein_domains.UNIPROT_HOST}
    for disclosure in (vep, clinvar, domains):
        assert disclosure["transport"] == "HTTPS"
        assert disclosure["rate_limit_policy"]


def test_variant_modules_do_not_use_pickle_shell_or_eval():
    forbidden = (
        re.compile(r"\bimport\s+pickle\b"),
        re.compile(r"\bpickle\.loads?\b"),
        re.compile(r"[,(]\s*shell\s*=\s*True"),
        re.compile(r"\beval\s*\("),
        re.compile(r"\bexec\s*\("),
        re.compile(r"\bos\.system\s*\("),
    )
    for name in (
        "variant_core.py",
        "variant_explorer.py",
        "ensembl_vep.py",
        "clinvar_evidence.py",
        "protein_domains.py",
    ):
        text = (REPO_ROOT / "modules" / name).read_text(encoding="utf-8")
        for pattern in forbidden:
            assert not pattern.search(text), f"{name} matches {pattern.pattern}"


def test_no_module_uses_shell_true_or_unsafe_deserialization():
    offenders = []
    for path in sorted((REPO_ROOT / "modules").glob("*.py")):
        text = path.read_text(encoding="utf-8")
        if re.search(r"[,(]\s*shell\s*=\s*True", text):
            offenders.append(f"{path.name}: shell=True")
        if re.search(r"\bpickle\.loads?\b", text):
            offenders.append(f"{path.name}: pickle.load")
        if re.search(r"\byaml\.load\s*\(", text):
            offenders.append(f"{path.name}: yaml.load")
    assert offenders == []
