"""Testes de dominios InterPro e anotacao UniProtKB. Offline, com rede injetada.

Os payloads sao recortes reduzidos de respostas reais de
`/interpro/api/entry/interpro/protein/uniprot/P38398` e
`https://rest.uniprot.org/uniprotkb/P38398.json`.
"""

from __future__ import annotations

import io
import json
from email.message import EmailMessage
from urllib.error import HTTPError, URLError

import pytest

from modules import protein_domains


class _FakeHandle:
    def __init__(self, body, url: str = "https://www.ebi.ac.uk/interpro/api/x", status: int = 200) -> None:
        self._body = body.encode("utf-8") if isinstance(body, str) else body
        self._url = url
        self.status = status

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


INTERPRO_PAGE = {
    "count": 2,
    "next": None,
    "results": [
        {
            "metadata": {
                "accession": "IPR001357",
                "name": "BRCT domain",
                "source_database": "interpro",
                "type": "domain",
                "integrated": None,
                "member_databases": {
                    "pfam": {"PF00533": "BRCA1 C Terminus (BRCT) domain"},
                    "smart": {"SM00292": "breast cancer carboxy-terminal domain"},
                },
                "go_terms": None,
            },
            "proteins": [
                {
                    "accession": "p38398",
                    "protein_length": 1863,
                    "source_database": "reviewed",
                    "entry_protein_locations": [
                        {
                            "fragments": [
                                {"start": 1642, "end": 1736, "dc-status": "CONTINUOUS"}
                            ],
                            "representative": False,
                        },
                        {
                            "fragments": [
                                {"start": 1756, "end": 1855, "dc-status": "CONTINUOUS"}
                            ],
                            "representative": False,
                        },
                    ],
                }
            ],
        },
        {
            "metadata": {
                "accession": "IPR001841",
                "name": "Zinc finger, RING-type",
                "source_database": "interpro",
                "type": "domain",
                "member_databases": {"pfam": {"PF13639": "Ring finger domain"}},
            },
            "proteins": [
                {
                    "accession": "p38398",
                    "protein_length": 1863,
                    "entry_protein_locations": [
                        {"fragments": [{"start": 24, "end": 65, "dc-status": "CONTINUOUS"}]}
                    ],
                }
            ],
        },
    ],
}

UNIPROT_ENTRY = {
    "primaryAccession": "P38398",
    "uniProtkbId": "BRCA1_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Breast cancer type 1 susceptibility protein"}
        }
    },
    "genes": [{"geneName": {"value": "BRCA1"}}],
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "sequence": {"length": 1863, "crc64": "F5F0F4F0F4F0F4F0"},
    "proteinExistence": "1: Evidence at protein level",
    "uniProtKBCrossReferences": [
        {"database": "PDB", "id": "1JM7"},
        {"database": "AlphaFoldDB", "id": "P38398"},
    ],
}


def _responder(payload, url: str = "https://www.ebi.ac.uk/interpro/api/x", status: int = 200):
    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(payload), url, status)

    return opener


def test_url_allowlist_covers_interpro_and_uniprot_only() -> None:
    assert protein_domains.url_is_allowed(
        "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398"
    )
    assert protein_domains.url_is_allowed("https://rest.uniprot.org/uniprotkb/P38398.json")
    assert not protein_domains.url_is_allowed("http://www.ebi.ac.uk/interpro/api/x")
    assert not protein_domains.url_is_allowed("https://www.ebi.ac.uk/other/api/x")
    assert not protein_domains.url_is_allowed("https://evil.example.org/interpro/api/x")
    assert not protein_domains.url_is_allowed("https://rest.uniprot.org/../etc/passwd")


@pytest.mark.parametrize("accession", ["P38398", "P38398.282", "p38398", "P38398-1", "Q9Y6K9"])
def test_valid_uniprot_accessions_are_accepted(accession: str) -> None:
    assert protein_domains.validate_uniprot_accession(accession).startswith(
        accession.split(".")[0].split("-")[0].upper()
    )


@pytest.mark.parametrize(
    "accession",
    ["", "   ", "NOTANACC", "../../etc/passwd", "P38398; rm -rf /", "<script>", "P" * 20],
)
def test_invalid_uniprot_accessions_are_refused(accession: str) -> None:
    with pytest.raises(protein_domains.DomainError) as exc:
        protein_domains.validate_uniprot_accession(accession)
    assert exc.value.category == "INVALID_INPUT"


def test_domains_are_parsed_with_source_database_preserved() -> None:
    result = protein_domains.fetch_protein_domains(
        accession="P38398", urlopen_fn=_responder(INTERPRO_PAGE)
    )
    assert result["status"] == "AVAILABLE"
    assert result["count"] == 2
    assert result["positional_count"] == 2
    brct = result["entries"][0]
    assert brct["entry_accession"] == "IPR001357"
    assert brct["source_database"] == "interpro"
    assert brct["type"] == "domain"
    assert brct["protein_length"] == 1863
    databases = {member["source_database"] for member in brct["member_databases"]}
    assert databases == {"pfam", "smart"}
    pfam = [m for m in brct["member_databases"] if m["source_database"] == "pfam"][0]
    assert pfam["accession"] == "PF00533"


def test_discontinuous_fragments_are_kept_separate() -> None:
    result = protein_domains.fetch_protein_domains(
        accession="P38398", urlopen_fn=_responder(INTERPRO_PAGE)
    )
    brct = result["entries"][0]
    assert len(brct["ranges"]) == 2
    spans = [
        (fragment["start_1based"], fragment["end_1based"])
        for interval in brct["ranges"]
        for fragment in interval["fragments"]
    ]
    assert spans == [(1642, 1736), (1756, 1855)]


def test_residue_coverage_is_positional_not_functional() -> None:
    result = protein_domains.fetch_protein_domains(
        accession="P38398", urlopen_fn=_responder(INTERPRO_PAGE)
    )
    inside = protein_domains.domains_covering_residue(result, residue_1based=1699)
    assert inside["covering_count"] == 1
    assert inside["covering"][0]["entry_accession"] == "IPR001357"
    assert "does not demonstrate a functional effect" in inside["note"]
    between = protein_domains.domains_covering_residue(result, residue_1based=1745)
    assert between["covering_count"] == 0
    assert "absence of annotation" in between["note"]
    ring = protein_domains.domains_covering_residue(result, residue_1based=30)
    assert ring["covering"][0]["entry_accession"] == "IPR001841"


def test_residue_coverage_boundaries_are_inclusive() -> None:
    result = protein_domains.fetch_protein_domains(
        accession="P38398", urlopen_fn=_responder(INTERPRO_PAGE)
    )
    for residue in (1642, 1736):
        assert protein_domains.domains_covering_residue(
            result, residue_1based=residue
        )["covering_count"] == 1
    assert protein_domains.domains_covering_residue(
        result, residue_1based=1641
    )["covering_count"] == 0


def test_residue_coverage_rejects_invalid_positions() -> None:
    for residue in (0, -5):
        with pytest.raises(protein_domains.DomainError):
            protein_domains.domains_covering_residue([], residue_1based=residue)
    with pytest.raises(protein_domains.DomainError):
        protein_domains.domains_covering_residue([], residue_1based="x")


def test_empty_interpro_result_is_not_found_not_error() -> None:
    result = protein_domains.fetch_protein_domains(
        accession="P38398", urlopen_fn=_responder({"count": 0, "results": [], "next": None})
    )
    assert result["status"] == "NOT_FOUND"
    assert result["entries"] == []
    assert "not evidence of absence of structure" in result["message"]


def test_http_404_yields_not_found_payload_without_raising() -> None:
    def opener(request, timeout=None):
        raise HTTPError(
            str(request.full_url), 404, "missing", hdrs=EmailMessage(), fp=io.BytesIO(b"")
        )

    result = protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=opener)
    assert result["status"] == "NOT_FOUND"
    assert "not evidence that the region is" in result["message"]


def test_interpro_408_is_timeout_not_absence(monkeypatch) -> None:
    monkeypatch.setattr(protein_domains.time, "sleep", lambda *_a: None)

    def opener(request, timeout=None):
        raise HTTPError(
            str(request.full_url), 408, "background", hdrs=EmailMessage(), fp=io.BytesIO(b"")
        )

    with pytest.raises(protein_domains.DomainError) as exc:
        protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=opener)
    assert exc.value.category == "TIMEOUT"


def test_network_failures_keep_their_category(monkeypatch) -> None:
    monkeypatch.setattr(protein_domains.time, "sleep", lambda *_a: None)

    def timeout_opener(request, timeout=None):
        raise TimeoutError("timed out")

    with pytest.raises(protein_domains.DomainError) as timed:
        protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=timeout_opener)
    assert timed.value.category == "TIMEOUT"
    assert "not an absence of domains" in str(timed.value)

    def net_opener(request, timeout=None):
        raise URLError("refused")

    with pytest.raises(protein_domains.DomainError) as net:
        protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=net_opener)
    assert net.value.category == "NETWORK_ERROR"


def test_pagination_outside_allowlist_is_refused() -> None:
    hostile = dict(INTERPRO_PAGE)
    hostile = json.loads(json.dumps(INTERPRO_PAGE))
    hostile["next"] = "https://evil.example.org/interpro/api/next"

    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(hostile))

    with pytest.raises(protein_domains.DomainError) as exc:
        protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=opener)
    assert exc.value.category == "INVALID_INPUT"


def test_pagination_is_bounded(monkeypatch) -> None:
    monkeypatch.setattr(protein_domains.time, "sleep", lambda *_a: None)
    looping = json.loads(json.dumps(INTERPRO_PAGE))
    looping["next"] = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P38398?page=2"

    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(looping))

    result = protein_domains.fetch_protein_domains(accession="P38398", urlopen_fn=opener)
    assert result["pages_fetched"] == protein_domains.MAX_PAGES
    assert result["truncated"] is True


def test_malformed_interpro_payload_is_parsing_error() -> None:
    with pytest.raises(protein_domains.DomainError) as exc:
        protein_domains.parse_interpro_entries({"nope": 1}, accession="P38398")
    assert exc.value.category == "PARSING_ERROR"
    with pytest.raises(protein_domains.DomainError):
        protein_domains.parse_interpro_entries("a string", accession="P38398")
    with pytest.raises(protein_domains.DomainError):
        protein_domains.parse_interpro_entries({"results": "not a list"}, accession="P38398")


def test_uniprot_entry_fields_are_retrieved_verbatim() -> None:
    entry = protein_domains.fetch_uniprot_protein(
        accession="P38398",
        urlopen_fn=_responder(UNIPROT_ENTRY, "https://rest.uniprot.org/uniprotkb/P38398.json"),
    )
    assert entry["accession"] == "P38398"
    assert entry["entry_name"] == "BRCA1_HUMAN"
    assert entry["protein_name"] == "Breast cancer type 1 susceptibility protein"
    assert entry["gene_names"] == ["BRCA1"]
    assert entry["organism"] == "Homo sapiens"
    assert entry["taxon_id"] == 9606
    assert entry["sequence_length"] == 1863
    assert entry["reviewed"] is True
    assert entry["pdb_cross_references"] == ["1JM7"]
    assert entry["alphafold_cross_references"] == ["P38398"]
    assert entry["evidence_status"] == "RETRIEVED"


def test_uniprot_version_suffix_is_stripped_before_lookup() -> None:
    entry = protein_domains.fetch_uniprot_protein(
        accession="P38398.282",
        urlopen_fn=_responder(UNIPROT_ENTRY, "https://rest.uniprot.org/uniprotkb/P38398.json"),
    )
    assert entry["accession"] == "P38398"


def test_network_disclosure_states_the_scientific_limit() -> None:
    disclosure = protein_domains.network_disclosure()
    assert disclosure["remote"] is True
    assert "not evidence of a functional effect" in disclosure["scientific_limit"]
    hosts = {service["host"] for service in disclosure["services"]}
    assert hosts == {protein_domains.INTERPRO_HOST, protein_domains.UNIPROT_HOST}
