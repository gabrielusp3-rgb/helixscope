"""Testes do cliente Ensembl VEP. Offline, com rede injetada.

Os payloads usados sao recortes reduzidos de respostas reais do endpoint
`/vep/human/region` (BRCA1, GRCh38), preservando os nomes de campo originais.
"""

from __future__ import annotations

import io
import json
from email.message import EmailMessage
from urllib.error import HTTPError, URLError

import pytest

from modules import ensembl_vep


class _FakeHandle:
    def __init__(self, body, url: str = "https://rest.ensembl.org/vep/human/region/x") -> None:
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


BRCA1_MISSENSE = [
    {
        "input": "17 43093464 . A T",
        "id": "17_43093464_A/T",
        "assembly_name": "GRCh38",
        "seq_region_name": "17",
        "start": 43093464,
        "end": 43093464,
        "strand": 1,
        "allele_string": "A/T",
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [
            {
                "transcript_id": "ENST00000357654",
                "gene_id": "ENSG00000012048",
                "gene_symbol": "BRCA1",
                "gene_symbol_source": "HGNC",
                "hgnc_id": "HGNC:1100",
                "biotype": "protein_coding",
                "strand": -1,
                "canonical": 1,
                "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
                "cdna_start": 2180,
                "cdna_end": 2180,
                "cds_start": 2067,
                "cds_end": 2067,
                "protein_start": 689,
                "protein_end": 689,
                "protein_id": "ENSP00000350283",
                "amino_acids": "S/R",
                "codons": "agT/agA",
                "hgvsc": "ENST00000357654.9:c.2067T>A",
                "hgvsp": "ENSP00000350283.3:p.Ser689Arg",
                "swissprot": ["P38398.282"],
                "uniprot_isoform": ["P38398-1"],
                "sift_score": 0.82,
                "sift_prediction": "tolerated",
                "polyphen_score": 0.062,
                "polyphen_prediction": "benign",
                "domains": [{"db": "PIRSF", "name": "PIRSF001734"}],
            },
            {
                "transcript_id": "ENST00000352993",
                "gene_id": "ENSG00000012048",
                "gene_symbol": "BRCA1",
                "biotype": "protein_coding",
                "strand": -1,
                "impact": "MODIFIER",
                "consequence_terms": ["intron_variant"],
            },
        ],
        "colocated_variants": [
            {
                "id": "rs80357382",
                "allele_string": "A/T",
                "start": 43093464,
                "end": 43093464,
                "clin_sig": ["uncertain_significance"],
            }
        ],
    }
]

INTERGENIC = [
    {
        "input": "2 30000000 . A G",
        "assembly_name": "GRCh38",
        "seq_region_name": "2",
        "start": 30000000,
        "end": 30000000,
        "strand": 1,
        "allele_string": "A/G",
        "most_severe_consequence": "intergenic_variant",
        "intergenic_consequences": [
            {"impact": "MODIFIER", "consequence_terms": ["intergenic_variant"]}
        ],
    }
]


def _responder(payload, url: str = "https://rest.ensembl.org/vep/human/region/x"):
    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(payload), url)

    return opener


def test_host_is_chosen_by_assembly() -> None:
    assert ensembl_vep.host_for_assembly("GRCh38.p14") == ensembl_vep.VEP_HOST_GRCH38
    assert ensembl_vep.host_for_assembly("GRCh37") == ensembl_vep.VEP_HOST_GRCH37
    assert ensembl_vep.host_for_assembly("hg19") == ensembl_vep.VEP_HOST_GRCH37


def test_missing_or_unknown_assembly_is_refused_not_defaulted() -> None:
    for value in ("", "T2T-CHM13v2.0", "GRCz11"):
        with pytest.raises(ensembl_vep.VepError) as exc:
            ensembl_vep.host_for_assembly(value)
        assert exc.value.category == "INVALID_INPUT"


def test_url_allowlist_blocks_other_hosts() -> None:
    assert ensembl_vep.url_is_allowed("https://rest.ensembl.org/vep/human/region/17:1-1/T")
    assert ensembl_vep.url_is_allowed("https://grch37.rest.ensembl.org/info/ping")
    assert not ensembl_vep.url_is_allowed("http://rest.ensembl.org/vep/human/region/x")
    assert not ensembl_vep.url_is_allowed("https://evil.example.org/vep/human/region/x")
    assert not ensembl_vep.url_is_allowed("https://rest.ensembl.org/admin/secret")
    assert not ensembl_vep.url_is_allowed("https://rest.ensembl.org.evil.tld/vep/")
    assert not ensembl_vep.url_is_allowed("file:///etc/passwd")
    assert not ensembl_vep.url_is_allowed("")


def test_network_disclosure_names_what_leaves_the_machine() -> None:
    disclosure = ensembl_vep.network_disclosure("GRCh38.p14")
    assert disclosure["remote"] is True
    assert disclosure["host"] == ensembl_vep.VEP_HOST_GRCH38
    assert any("position" in item for item in disclosure["data_sent"])
    assert any("reference genome files" in item for item in disclosure["data_not_sent"])


def test_all_transcripts_are_preserved_and_terms_untouched() -> None:
    result = ensembl_vep.annotate_region(
        contig="17",
        position_1based=43093464,
        ref="A",
        alt="T",
        assembly="GRCh38.p14",
        urlopen_fn=_responder(BRCA1_MISSENSE),
    )
    assert result["transcript_count"] == 2
    assert result["transcript_ambiguity"] is True
    assert result["canonical_transcript_count"] == 1
    assert result["most_severe_consequence"] == "missense_variant"
    terms = [row["consequence_terms"] for row in result["transcript_consequences"]]
    assert ["missense_variant"] in terms
    assert ["intron_variant"] in terms
    assert result["consequence_vocabulary"].startswith("Sequence Ontology")
    assert result["evidence_status"] == "RETRIEVED"


def test_protein_fields_are_split_without_being_renamed() -> None:
    result = ensembl_vep.annotate_region(
        contig="17",
        position_1based=43093464,
        ref="A",
        alt="T",
        assembly="GRCh38.p14",
        urlopen_fn=_responder(BRCA1_MISSENSE),
    )
    row = result["transcript_consequences"][0]
    assert row["amino_acids"] == "S/R"
    assert row["reference_amino_acid"] == "S"
    assert row["variant_amino_acid"] == "R"
    assert row["reference_codon"] == "agT"
    assert row["variant_codon"] == "agA"
    assert row["protein_start"] == 689
    assert row["protein_id"] == "ENSP00000350283"
    assert row["swissprot"] == ["P38398.282"]
    assert row["hgvsp"] == "ENSP00000350283.3:p.Ser689Arg"


def test_sift_and_polyphen_are_labelled_as_external_predictions() -> None:
    result = ensembl_vep.annotate_region(
        contig="17",
        position_1based=43093464,
        ref="A",
        alt="T",
        assembly="GRCh38.p14",
        urlopen_fn=_responder(BRCA1_MISSENSE),
    )
    row = result["transcript_consequences"][0]
    assert row["sift_prediction"] == "tolerated"
    assert row["polyphen_prediction"] == "benign"
    note = row["prediction_note"].lower()
    assert "predicted" in note
    assert "not a helixscope pathogenicity call" in note
    assert "not a clinical classification" in note


def test_consequence_classes_never_replace_the_source_term() -> None:
    assert ensembl_vep.classify_consequence("missense_variant") == "missense"
    assert ensembl_vep.classify_consequence("stop_gained") == "stop_gained"
    assert ensembl_vep.classify_consequence("frameshift_variant") == "frameshift"
    assert ensembl_vep.classify_consequence("synonymous_variant") == "synonymous"
    assert ensembl_vep.classify_consequence("intron_variant") == "intronic"
    assert ensembl_vep.classify_consequence("5_prime_UTR_variant") == "UTR"
    assert ensembl_vep.classify_consequence("intergenic_variant") == "intergenic"
    assert ensembl_vep.classify_consequence("splice_acceptor_variant") == "splice_related"
    assert ensembl_vep.classify_consequence("something_new_from_SO") == "unclassified"


def test_intergenic_variant_yields_no_protein_mapping() -> None:
    result = ensembl_vep.annotate_region(
        contig="2",
        position_1based=30000000,
        ref="A",
        alt="G",
        assembly="GRCh38.p14",
        urlopen_fn=_responder(INTERGENIC),
    )
    assert result["transcript_consequences"] == []
    assert result["protein_mapped_transcript_count"] == 0
    assert result["intergenic_consequences"][0]["consequence_terms"] == ["intergenic_variant"]


def test_colocated_variants_are_attributed_to_ensembl() -> None:
    result = ensembl_vep.annotate_region(
        contig="17",
        position_1based=43093464,
        ref="A",
        alt="T",
        assembly="GRCh38.p14",
        urlopen_fn=_responder(BRCA1_MISSENSE),
    )
    colocated = result["colocated_variants"][0]
    assert colocated["id"] == "rs80357382"
    assert colocated["clin_sig"] == ["uncertain_significance"]
    assert "not a HelixScope assertion" in colocated["source"]


def test_empty_response_is_not_found_not_silence() -> None:
    with pytest.raises(ensembl_vep.VepError) as exc:
        ensembl_vep.parse_vep_response([])
    assert exc.value.category == "NOT_FOUND"


def test_malformed_response_is_parsing_error() -> None:
    for payload in ("a string", 42):
        with pytest.raises(ensembl_vep.VepError) as exc:
            ensembl_vep.parse_vep_response(payload)
        assert exc.value.category == "PARSING_ERROR"


def test_http_status_maps_to_distinct_categories() -> None:
    def boom(code: int):
        def opener(request, timeout=None):
            raise HTTPError(
                str(request.full_url), code, "err", hdrs=EmailMessage(), fp=io.BytesIO(b"")
            )

        return opener

    expected = {404: "NOT_FOUND", 400: "INVALID_INPUT"}
    for code, category in expected.items():
        with pytest.raises(ensembl_vep.VepError) as exc:
            ensembl_vep.annotate_region(
                contig="17",
                position_1based=1,
                ref="A",
                alt="T",
                assembly="GRCh38",
                urlopen_fn=boom(code),
            )
        assert exc.value.category == category


def test_rate_limit_and_server_error_are_not_collapsed_to_not_found(monkeypatch) -> None:
    monkeypatch.setattr(ensembl_vep.time, "sleep", lambda *_a: None)
    for code, category in ((429, "RATE_LIMITED"), (503, "SERVICE_UNAVAILABLE")):
        def opener(request, timeout=None, _code=code):
            headers = EmailMessage()
            headers["Retry-After"] = "1"
            raise HTTPError(
                str(request.full_url), _code, "err", hdrs=headers, fp=io.BytesIO(b"")
            )

        with pytest.raises(ensembl_vep.VepError) as exc:
            ensembl_vep.annotate_region(
                contig="17",
                position_1based=1,
                ref="A",
                alt="T",
                assembly="GRCh38",
                urlopen_fn=opener,
            )
        assert exc.value.category == category
        assert exc.value.category != "NOT_FOUND"


def test_timeout_and_network_error_are_distinguished(monkeypatch) -> None:
    monkeypatch.setattr(ensembl_vep.time, "sleep", lambda *_a: None)

    def timeout_opener(request, timeout=None):
        raise TimeoutError("timed out")

    with pytest.raises(ensembl_vep.VepError) as timed:
        ensembl_vep.annotate_region(
            contig="17", position_1based=1, ref="A", alt="T",
            assembly="GRCh38", urlopen_fn=timeout_opener,
        )
    assert timed.value.category == "TIMEOUT"

    def net_opener(request, timeout=None):
        raise URLError("connection refused")

    with pytest.raises(ensembl_vep.VepError) as net:
        ensembl_vep.annotate_region(
            contig="17", position_1based=1, ref="A", alt="T",
            assembly="GRCh38", urlopen_fn=net_opener,
        )
    assert net.value.category == "NETWORK_ERROR"


def test_oversized_response_is_resource_limit() -> None:
    def opener(request, timeout=None):
        return _FakeHandle(b"[" + b"0" * (ensembl_vep.MAX_RESPONSE_BYTES + 1) + b"]")

    with pytest.raises(ensembl_vep.VepError) as exc:
        ensembl_vep.annotate_region(
            contig="17", position_1based=1, ref="A", alt="T",
            assembly="GRCh38", urlopen_fn=opener,
        )
    assert exc.value.category == "RESOURCE_LIMIT"


def test_response_from_unexpected_host_is_refused() -> None:
    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(BRCA1_MISSENSE), "https://evil.example.org/vep/")

    with pytest.raises(ensembl_vep.VepError) as exc:
        ensembl_vep.annotate_region(
            contig="17", position_1based=43093464, ref="A", alt="T",
            assembly="GRCh38", urlopen_fn=opener,
        )
    assert exc.value.category == "INVALID_INPUT"


@pytest.mark.parametrize(
    "contig",
    ["", "17/../admin", "17?x=1", "17 43", "17&y=2", "17#frag", "17%2f"],
)
def test_unsafe_contig_tokens_are_refused(contig: str) -> None:
    with pytest.raises(ensembl_vep.VepError) as exc:
        ensembl_vep.annotate_region(
            contig=contig, position_1based=1, ref="A", alt="T", assembly="GRCh38"
        )
    assert exc.value.category == "INVALID_INPUT"


@pytest.mark.parametrize(
    "notation",
    ["", "NM_007294.4:c.1974G>C and rm -rf /", "a?b", "x/../y", "n" * 300],
)
def test_unsafe_hgvs_notation_is_refused(notation: str) -> None:
    with pytest.raises(ensembl_vep.VepError):
        ensembl_vep.annotate_hgvs(notation=notation, assembly="GRCh38")


def test_hgvs_endpoint_uses_the_service_not_a_local_parser() -> None:
    result = ensembl_vep.annotate_hgvs(
        notation="ENST00000357654.9:c.2067T>A",
        assembly="GRCh38",
        urlopen_fn=_responder(BRCA1_MISSENSE, "https://rest.ensembl.org/vep/homo_sapiens/hgvs/x"),
    )
    assert result["transcript_count"] == 2
    assert "vep/homo_sapiens/hgvs" in result["query"] or "vep/" in result["query"]


def test_identifier_endpoint_validates_the_identifier() -> None:
    with pytest.raises(ensembl_vep.VepError):
        ensembl_vep.annotate_identifier(identifier="rs123; drop", assembly="GRCh38")
    result = ensembl_vep.annotate_identifier(
        identifier="rs80357382",
        assembly="GRCh38",
        urlopen_fn=_responder(BRCA1_MISSENSE, "https://rest.ensembl.org/vep/homo_sapiens/id/x"),
    )
    assert result["transcript_count"] == 2


def test_region_endpoint_requires_an_alternate_allele() -> None:
    with pytest.raises(ensembl_vep.VepError) as exc:
        ensembl_vep.annotate_region(
            contig="17", position_1based=1, ref="A", alt="", assembly="GRCh38"
        )
    assert exc.value.category == "INVALID_INPUT"


def test_service_availability_reports_failure_without_raising() -> None:
    def opener(request, timeout=None):
        raise URLError("down")

    status = ensembl_vep.service_availability(urlopen_fn=opener, assembly="GRCh38")
    assert status["available"] is False
    assert status["category"] == "NETWORK_ERROR"
    assert "not an absence of consequences" in status["reason"]


def test_service_availability_reports_release_when_reachable() -> None:
    calls = {"n": 0}

    def opener(request, timeout=None):
        calls["n"] += 1
        body = {"ping": 1} if "ping" in str(request.full_url) else {"releases": [116]}
        return _FakeHandle(json.dumps(body), "https://rest.ensembl.org/info/ping")

    status = ensembl_vep.service_availability(urlopen_fn=opener, assembly="GRCh38")
    assert status["available"] is True
    assert status["release"] == 116
    assert calls["n"] == 2
