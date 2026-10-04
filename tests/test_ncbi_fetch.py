"""Testes unitarios do modulo modules.ncbi_fetch.

Cobrem a interpretacao de localizacoes de features GenBank, a extracao das regioes
codificadoras anotadas e a validacao de argumentos de fetch_by_accession. Nenhum
teste faz chamadas de rede: os registros usados sao dicionarios sinteticos com o
mesmo formato retornado por fetch_by_accession.
"""

import urllib.error

import pytest

from modules import ncbi_fetch


@pytest.fixture
def annotated_record() -> dict:
    """Registro sintetico com duas CDS validas e uma feature sem localizacao."""
    return {
        "accession": "NC_000000.1",
        "description": "Synthetic record for tests",
        "organism": "Synthetic organism",
        "length": 1000,
        "sequence": "A" * 1000,
        "features": [
            {"type": "source", "location": "[0:1000](+)", "qualifiers": {}},
            {
                "type": "CDS",
                "location": "[10:100](+)",
                "qualifiers": {
                    "gene": ["orf1ab"],
                    "product": ["polyprotein"],
                    "protein_id": ["ABC12345.1"],
                    "codon_start": ["1"],
                },
            },
            {
                "type": "CDS",
                "location": "join{[200:250](-), [300:350](-)}",
                "qualifiers": {"codon_start": ["2"]},
            },
            {"type": "CDS", "location": "", "qualifiers": {}},
        ],
    }


def test_parse_genbank_location_simple():
    segments = ncbi_fetch.parse_genbank_location("[0:1234](+)")
    assert segments == [{"start": 0, "end": 1234, "strand": "+"}]


def test_parse_genbank_location_compound():
    segments = ncbi_fetch.parse_genbank_location(
        "join{[0:100](+), [200:300](+)}"
    )
    assert len(segments) == 2
    assert segments[1]["start"] == 200
    assert segments[1]["end"] == 300


def test_parse_genbank_location_fuzzy_and_reverse():
    segments = ncbi_fetch.parse_genbank_location("[<0:>500](-)")
    assert segments == [{"start": 0, "end": 500, "strand": "-"}]


def test_parse_genbank_location_invalid_returns_empty():
    assert ncbi_fetch.parse_genbank_location("") == []
    assert ncbi_fetch.parse_genbank_location("nonsense") == []
    assert ncbi_fetch.parse_genbank_location("[100:100](+)") == []


def test_extract_cds_regions_skips_non_cds(annotated_record):
    regions = ncbi_fetch.extract_cds_regions(annotated_record)
    assert len(regions) == 2


def test_extract_cds_regions_reads_qualifiers(annotated_record):
    first = ncbi_fetch.extract_cds_regions(annotated_record)[0]
    assert first["gene"] == "orf1ab"
    assert first["product"] == "polyprotein"
    assert first["protein_id"] == "ABC12345.1"
    assert first["strand"] == "+"
    assert first["segments"] == [(10, 100)]
    assert first["length_nt"] == 90
    assert first["codon_start"] == 1


def test_extract_cds_regions_handles_joined_reverse_strand(annotated_record):
    second = ncbi_fetch.extract_cds_regions(annotated_record)[1]
    assert second["strand"] == "-"
    assert second["segments"] == [(200, 250), (300, 350)]
    assert second["start"] == 200
    assert second["end"] == 350
    assert second["length_nt"] == 100
    assert second["codon_start"] == 2


def test_extract_cds_regions_without_features():
    assert ncbi_fetch.extract_cds_regions({}) == []
    assert ncbi_fetch.extract_cds_regions({"features": []}) == []


def test_fetch_by_accession_requires_accession_and_email():
    with pytest.raises(ncbi_fetch.NCBIQueryError) as empty_id:
        ncbi_fetch.fetch_by_accession("", "user@example.com")
    assert empty_id.value.category == "invalid_input"
    with pytest.raises(ValueError):
        ncbi_fetch.fetch_by_accession("NC_000000.1", "  ")


def test_classify_identifier_nucleotide_and_version():
    with_version = ncbi_fetch.classify_identifier("NM_000518.5")
    assert with_version["kind"] == "nucleotide"
    assert with_version["suggested_db"] == "nucleotide"
    assert with_version["has_version"] is True
    without = ncbi_fetch.classify_identifier("NM_000518")
    assert without["kind"] == "nucleotide"
    assert without["has_version"] is False


def test_classify_identifier_protein_gene_wgs_assembly():
    protein = ncbi_fetch.classify_identifier("NP_000509.1")
    assert protein["kind"] == "protein"
    assert protein["suggested_db"] == "protein"
    gene = ncbi_fetch.classify_identifier("3043")
    assert gene["kind"] == "gene_id"
    assert gene["suggested_db"] == "gene"
    wgs = ncbi_fetch.classify_identifier("AGTM000000000.1")
    assert wgs["kind"] == "wgs_master"
    wgs6 = ncbi_fetch.classify_identifier("JBFSEQ000000000")
    assert wgs6["kind"] == "wgs_master"
    assembly = ncbi_fetch.classify_identifier("GCA_000241425.1")
    assert assembly["kind"] == "assembly"


def test_classify_identifier_organism_and_invalid():
    organism = ncbi_fetch.classify_identifier("Daubentonia madagascariensis")
    assert organism["kind"] == "organism_name"
    typo = ncbi_fetch.classify_identifier("Dalbentonia madagascariensis")
    assert typo["kind"] == "organism_name"
    empty = ncbi_fetch.classify_identifier("   ")
    assert empty["kind"] == "empty"
    unknown = ncbi_fetch.classify_identifier("ZZZZ99999.1")
    assert unknown["kind"] in {"unknown", "wgs_master"}


def test_resolve_database_routes_numeric_id_to_gene():
    db, info, warnings = ncbi_fetch.resolve_database("3043", "nucleotide")
    assert db == "gene"
    assert info["kind"] == "gene_id"
    assert warnings


def test_resolve_database_routes_protein_away_from_nucleotide():
    db, info, warnings = ncbi_fetch.resolve_database("NP_000509.1", "nucleotide")
    assert db == "protein"
    assert info["kind"] == "protein"
    assert warnings


def test_resolve_database_rejects_organism_and_assembly():
    with pytest.raises(ncbi_fetch.NCBIQueryError) as organism:
        ncbi_fetch.resolve_database("Daubentonia madagascariensis", "nucleotide")
    assert organism.value.category == "invalid_input"
    with pytest.raises(ncbi_fetch.NCBIQueryError) as assembly:
        ncbi_fetch.resolve_database("GCA_000241425.1", "nucleotide")
    assert assembly.value.category == "wrong_database"


def test_sequence_payload_undefined_is_not_fake_bases():
    class FakeSeq:
        defined = False

        def __str__(self) -> str:
            raise RuntimeError("should not stringify undefined sequence")

    text, available = ncbi_fetch._sequence_payload(FakeSeq(), RuntimeError)
    assert text == ""
    assert available is False


def test_locus_length_from_genbank_text():
    raw = (
        "LOCUS       NM_007294               7094 bp    mRNA    linear   "
        "PRI 13-AUG-2026\n"
    )
    assert ncbi_fetch._locus_length_from_genbank_text(raw) == 7094
    protein = "LOCUS       NP_000509                147 aa            linear\n"
    assert ncbi_fetch._locus_length_from_genbank_text(protein) == 147
    assert ncbi_fetch._locus_length_from_genbank_text("") == 0


def test_http_error_categories():
    not_found = ncbi_fetch._http_error(
        "NP_000509.1", "nucleotide", urllib.error.HTTPError(
            url="https://eutils.ncbi.nlm.nih.gov",
            code=400,
            msg="Bad Request",
            hdrs=None,
            fp=None,
        )
    )
    assert not_found.category == "not_found"
    limited = ncbi_fetch._http_error(
        "NM_000518.5",
        "nucleotide",
        urllib.error.HTTPError(
            url="https://eutils.ncbi.nlm.nih.gov",
            code=429,
            msg="Too Many Requests",
            hdrs=None,
            fp=None,
        ),
    )
    assert limited.category == "source_unavailable"
    assert limited.failure_kind == "rate_limited"


def test_fetch_empty_response_is_source_unavailable(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)
    monkeypatch.setattr(ncbi_fetch, "_fetch_sequence_summary", lambda *_args: None)

    class FakeHandle:
        def read(self) -> str:
            return "   "

        def close(self) -> None:
            return None

    import Bio.Entrez as entrez

    monkeypatch.setattr(entrez, "efetch", lambda **_kwargs: FakeHandle())
    with pytest.raises(ncbi_fetch.NCBIQueryError) as exc:
        ncbi_fetch.fetch_by_accession("NM_000518.5", "user@example.com")
    assert exc.value.category == "source_unavailable"
    assert "vazia" in str(exc.value)


def test_fetch_timeout_is_source_unavailable(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)
    monkeypatch.setattr(ncbi_fetch, "_fetch_sequence_summary", lambda *_args: None)

    import Bio.Entrez as entrez

    def boom(**_kwargs):
        raise TimeoutError("timed out")

    monkeypatch.setattr(entrez, "efetch", boom)
    with pytest.raises(ncbi_fetch.NCBIQueryError) as exc:
        ncbi_fetch.fetch_by_accession("NM_000518.5", "user@example.com")
    assert exc.value.category == "source_unavailable"


def test_classify_wgs_contig_accession():
    info = ncbi_fetch.classify_identifier("JBFSEQ010000002.1")
    assert info["kind"] == "wgs_master"
    assert info["suggested_db"] == "nucleotide"


def test_summary_length_parses_large_contig():
    class IntegerElement(int):
        def __str__(self) -> str:
            return f"IntegerElement({int(self)}, attributes={{}})"

    assert ncbi_fetch._summary_length({"Length": IntegerElement(290592686)}) == 290592686
    assert ncbi_fetch._summary_length({"Length": "290592686"}) == 290592686


def test_large_sequence_does_not_invent_bases(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)
    monkeypatch.setattr(
        ncbi_fetch,
        "_fetch_sequence_summary",
        lambda *_args: {
            "Title": "Daubentonia scaffold_2",
            "AccessionVersion": "JBFSEQ010000002.1",
            "Length": 290592686,
        },
    )
    monkeypatch.setattr(
        ncbi_fetch,
        "_large_record_payload",
        lambda accession, db, summary, reported_length: {
            "accession": accession,
            "description": summary["Title"],
            "organism": "Daubentonia madagascariensis",
            "length": reported_length,
            "sequence": "",
            "features": [],
            "record_kind": "large_sequence",
            "sequence_available": False,
            "warnings": ["too large"],
            "related_records": [],
            "gene": None,
            "genbank_date": "",
        },
    )
    record = ncbi_fetch.fetch_by_accession(
        "JBFSEQ010000002.1", "user@example.com", "nucleotide"
    )
    assert record["sequence"] == ""
    assert record["sequence_available"] is False
    assert record["length"] == 290592686
    assert record["record_kind"] == "large_sequence"


def test_search_records_rejects_empty_query_and_unknown_database():
    with pytest.raises(ncbi_fetch.NCBIQueryError) as empty:
        ncbi_fetch.search_records("   ", "user@example.com")
    assert empty.value.category == "invalid_input"
    with pytest.raises(ncbi_fetch.NCBIQueryError) as bad_db:
        ncbi_fetch.search_records("BRCA1", "user@example.com", db="sra")
    assert bad_db.value.category == "wrong_database"
    with pytest.raises(ValueError):
        ncbi_fetch.search_records("BRCA1", " ")


def test_search_records_empty_id_list_is_not_found(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)

    class Handle:
        def close(self) -> None:
            return None

    def fake_esearch(**_kwargs):
        return Handle()

    def fake_read(_handle):
        return {"Count": "0", "IdList": []}

    from Bio import Entrez

    monkeypatch.setattr(Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(Entrez, "read", fake_read)
    result = ncbi_fetch.search_records("xyzzy_no_such_record", "user@example.com")
    assert result["status"] == "no_records_found"
    assert result["count"] == 0
    assert result["records"] == []
    assert result["source"] == "NCBI Entrez"


def test_search_records_parses_esummary_hits(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)
    calls = {"n": 0}

    class Handle:
        def close(self) -> None:
            return None

    def fake_esearch(**_kwargs):
        return Handle()

    def fake_esummary(**_kwargs):
        return Handle()

    def fake_read(_handle):
        calls["n"] += 1
        if calls["n"] == 1:
            return {"Count": "1", "IdList": ["123"]}
        return [
            {
                "Id": "123",
                "Caption": "NM_000059.4",
                "Title": "Homo sapiens BRCA2",
                "Organism": "Homo sapiens",
                "Length": 11386,
            }
        ]

    from Bio import Entrez

    monkeypatch.setattr(Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(Entrez, "esummary", fake_esummary)
    monkeypatch.setattr(Entrez, "read", fake_read)
    result = ncbi_fetch.search_records("BRCA2", "user@example.com")
    assert result["status"] == "retrieved"
    assert result["count"] == 1
    assert result["records"][0]["accession"] == "NM_000059.4"
    assert result["records"][0]["organism"] == "Homo sapiens"
    assert result["records"][0]["length"] == 11386
    assert result["records"][0]["title"] == "Homo sapiens BRCA2"


def test_homology_search_uses_ncbi_blast_url_api():
    from modules import blast_search

    availability = ncbi_fetch.homology_search_availability()
    assert availability["available"] is True
    assert availability["endpoint"] == blast_search.BLAST_ENDPOINT
    assert availability["backend"] == "NCBI BLAST Common URL API"
    local = availability["local_blast_plus"]
    assert isinstance(local, dict)
    assert "BLAST" in availability["reason"]
    assert "locally" in availability["reason"].lower() or "Blast.cgi" in availability["reason"]


def test_search_records_pagination_retstart(monkeypatch):
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)
    seen = {}

    class Handle:
        def close(self) -> None:
            return None

    def fake_esearch(**kwargs):
        seen.update(kwargs)
        return Handle()

    def fake_read(_handle):
        return {"Count": "40", "IdList": []}

    from Bio import Entrez

    monkeypatch.setattr(Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(Entrez, "read", fake_read)
    result = ncbi_fetch.search_records(
        "BRCA1", "user@example.com", retstart=20, retmax=15
    )
    assert seen["retstart"] == 20
    assert seen["retmax"] == 15
    assert result["retstart"] == 20
    assert result["status"] == "no_records_found"


def test_split_accession_version_and_explorer_omits_empty():
    assert ncbi_fetch.split_accession_version("NM_007294.4") == (
        "NM_007294",
        "NM_007294.4",
    )
    assert ncbi_fetch.split_accession_version("NM_007294") == ("NM_007294", "")
    fields = ncbi_fetch.record_explorer_fields(
        {
            "accession": "NM_007294.4",
            "description": "test",
            "organism": "",
            "taxonomy": [],
            "molecule": "mRNA",
            "length": 10,
            "sequence": "",
            "sequence_available": False,
            "features": [],
            "references": [],
            "related_records": [],
        }
    )
    assert fields["accession"] == "NM_007294"
    assert fields["version"] == "NM_007294.4"
    assert "organism" not in fields
    assert "taxonomy" not in fields
    assert fields["sequence_available"] is False
    assert "sequence" not in fields


def test_extract_annotated_features_marks_unrecoverable_coords(annotated_record):
    annotated_record["features"].append(
        {
            "type": "exon",
            "location": "[0:5000](+)",
            "qualifiers": {"gene": ["too_long"]},
        }
    )
    annotated_record["features"].append(
        {
            "type": "intron",
            "location": "nonsense",
            "qualifiers": {},
        }
    )
    items = ncbi_fetch.extract_annotated_features(annotated_record)
    cds_items = [item for item in items if item["type"] == "CDS"]
    assert any(item["status"] == "AVAILABLE" for item in cds_items)
    assert any(item["status"] == "UNAVAILABLE" for item in cds_items)
    exon = next(item for item in items if item["type"] == "exon")
    intron = next(item for item in items if item["type"] == "intron")
    assert exon["status"] == "UNAVAILABLE"
    assert intron["status"] == "UNAVAILABLE"
    assert exon["region_sequence"] == ""
    join_cds = next(item for item in cds_items if item["is_join"])
    assert join_cds["is_complement"] is True
    assert join_cds["status"] == "AVAILABLE"


def test_sequence_for_workspace_copies_verbatim(annotated_record):
    annotated_record["sequence_available"] = True
    annotated_record["database"] = "nucleotide"
    annotated_record["source"] = "NCBI Entrez"
    annotated_record["retrieved_at_utc"] = "2026-01-01T00:00:00Z"
    original = annotated_record["sequence"]
    payload = ncbi_fetch.sequence_for_workspace(annotated_record)
    assert payload["sequence"] == original
    assert payload["sequence"] is not original or payload["sequence"] == original
    assert payload["hash"]
    assert payload["molecule"] == "DNA"
    assert payload["accession"] == "NC_000000"
    missing = dict(annotated_record)
    missing["sequence"] = ""
    missing["sequence_available"] = False
    with pytest.raises(ValueError):
        ncbi_fetch.sequence_for_workspace(missing)


def test_feature_spans_for_analysis_uses_retrieved_coordinates(annotated_record):
    items = ncbi_fetch.extract_annotated_features(annotated_record)
    simple = next(item for item in items if item["type"] == "CDS" and not item["is_join"])
    spans = ncbi_fetch.feature_spans_for_analysis(simple, annotated_record["sequence"])
    assert spans[0]["sequence"] == annotated_record["sequence"][10:100]
    assert spans[0]["status"] == "RETRIEVED"
    unavailable = {
        "status": "UNAVAILABLE",
        "unavailable_reason": "out of range",
        "segments": [(0, 1)],
        "strand": "+",
        "type": "exon",
        "gene": "x",
    }
    with pytest.raises(ValueError):
        ncbi_fetch.feature_spans_for_analysis(unavailable, annotated_record["sequence"])



