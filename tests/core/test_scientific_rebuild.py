"""Non-3D scientific rebuild: Core domain pipelines and explanations."""

from __future__ import annotations

import pytest

from helixscope_core import analyze_dna, analyze_rna, analyze_protein, explain_metrics
from helixscope_core.crispr import analyze_crispr
from helixscope_core.motif import analyze_motif
from modules import ncbi_fetch, protein_domains
from modules.explain import explain


def test_analyze_dna_domain_pipeline_includes_optional_arrays() -> None:
    result = analyze_dna(
        "ATGC" * 40,
        include_orfs=True,
        include_restriction=True,
        include_profiles=True,
        include_cpg=True,
        include_kmers=True,
        include_santalucia=True,
        profile_window=50,
        profile_step=25,
        kmer_k=3,
    )
    assert result["status"] == "COMPUTED"
    assert result["identity_scope"] == "sequence_derived"
    assert "organism" not in result or result.get("organism") in {None, "", "unknown"}
    assert result["windowed_profiles"]
    assert "kmer_summary" in result
    assert "santalucia_tm" in result
    assert "sequence_feature_tracks" in result or result.get("sequence_features") is not None
    metrics = explain_metrics("dna", result)
    assert "gc_content" in metrics
    assert "does not identify the organism" in metrics["gc_content"]["limitations"].lower()
    assert metrics["gc_content"]["plain_meaning"]


def test_analyze_dna_fasta_header_is_not_validated_as_residues() -> None:
    result = analyze_dna(">desc 9606 numbers\nATGCATGCATGC\n")
    assert result["status"] == "COMPUTED"
    assert result["length"] == 12


def test_analyze_rna_codon_metrics_are_scoped() -> None:
    seq = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA"
    result = analyze_rna(seq, include_codon_metrics=True)
    assert result["status"] == "COMPUTED"
    assert "codon_metrics_scope" in result
    assert "cai" in result
    assert "coding sequence of a known gene" in result["codon_metrics_scope"]
    metrics = explain_metrics("rna", result)
    assert "rscu" in metrics or "cai" in metrics


def test_analyze_protein_does_not_invent_uniprot_identity() -> None:
    result = analyze_protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
    assert result["identity_scope"] == "sequence_derived"
    assert "uniprot" not in str(result.get("protein_name") or "").lower()
    metrics = explain_metrics("protein", result)
    assert "gravy" in metrics
    assert "not a measured solubility" in metrics["gravy"]["interpretation"].lower()


def test_uniprot_extra_fields_are_empty_when_absent_from_payload() -> None:
    entry = protein_domains.fetch_uniprot_protein(
        accession="P38398",
        urlopen_fn=_uniprot_opener(),
    )
    assert entry["evidence_status"] == "RETRIEVED"
    assert entry["function_comments"] == []
    assert entry["retrieved_sequence"] == ""
    assert entry["retrieved_sequence_status"] == "UNAVAILABLE"
    assert entry["identity_scope"] == "database_retrieved"


def test_uniprot_parses_function_and_sequence_when_present() -> None:
    payload = {
        "primaryAccession": "P01308",
        "uniProtkbId": "INS_HUMAN",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Insulin"}}},
        "genes": [{"geneName": {"value": "INS"}}],
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "sequence": {"length": 110, "crc64": "AAAA", "value": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"},
        "comments": [
            {"commentType": "FUNCTION", "texts": [{"value": "Decreases blood glucose."}]}
        ],
        "features": [
            {
                "type": "Signal",
                "description": "Signal peptide",
                "location": {"start": {"value": 1}, "end": {"value": 24}},
            }
        ],
        "references": [
            {
                "citation": {
                    "title": "A paper",
                    "journal": "Nature",
                    "publicationDate": "1980",
                    "citationCrossReferences": [{"database": "PubMed", "id": "123"}],
                }
            }
        ],
        "entryAudit": {"entryVersion": 7, "sequenceVersion": 1},
        "uniProtKBCrossReferences": [],
        "proteinExistence": "1: Evidence at protein level",
    }

    def opener(request, timeout=None):
        return _Handle(payload)

    entry = protein_domains.fetch_uniprot_protein(accession="P01308", urlopen_fn=opener)
    assert entry["function_comments"] == ["Decreases blood glucose."]
    assert entry["retrieved_sequence"].startswith("MALW")
    assert entry["features"][0]["start"] == 1
    assert entry["literature"][0]["pmid"] == "123"
    assert entry["entry_audit"]["entry_version"] == 7


def test_pubmed_resolve_does_not_reroute_pmid_to_gene() -> None:
    db, info, warnings = ncbi_fetch.resolve_database("24144637", "pubmed")
    assert db == "pubmed"
    assert info["kind"] == "pmid"
    assert warnings == []
    gene_db, gene_info, gene_warnings = ncbi_fetch.resolve_database("3043", "nucleotide")
    assert gene_db == "gene"
    assert gene_info["kind"] == "gene_id"
    assert gene_warnings


def test_pubmed_fetch_parses_esummary(monkeypatch) -> None:
    monkeypatch.setattr(ncbi_fetch.time, "sleep", lambda _seconds: None)

    class Handle:
        def close(self) -> None:
            return None

        def read(self) -> str:
            return "Abstract text from Entrez."

    def fake_esummary(**_kwargs):
        return Handle()

    def fake_efetch(**_kwargs):
        return Handle()

    def fake_read(_handle):
        return [
            {
                "Id": "24144637",
                "Title": "A PubMed title",
                "FullJournalName": "Nucleic Acids Research",
                "PubDate": "2013",
                "DOI": "10.1093/nar/gkt111",
                "HasAbstract": "1",
                "AuthorList": ["Doe J", "Roe A"],
            }
        ]

    from Bio import Entrez

    monkeypatch.setattr(Entrez, "esummary", fake_esummary)
    monkeypatch.setattr(Entrez, "efetch", fake_efetch)
    monkeypatch.setattr(Entrez, "read", fake_read)
    record = ncbi_fetch.fetch_by_accession("24144637", "user@example.com", db="pubmed")
    assert record["record_kind"] == "pubmed"
    assert record["title"] == "A PubMed title"
    assert record["authors"] == ["Doe J", "Roe A"]
    assert record["sequence_available"] is False
    assert record["abstract_status"] == "RETRIEVED"
    assert "Abstract text" in record["abstract"]


def test_analyze_motif_and_crispr_domain_results() -> None:
    motif = analyze_motif("ACGTACGTACGT", "ACGT", molecule="DNA")
    assert motif["n_hits"] >= 1
    assert motif["sequence_feature_tracks"]
    assert "not proven biological function" in motif["hit_meaning"]
    crispr = analyze_crispr("ACGTACGTACGTACGTACGTAGG", "SpCas9")
    assert crispr["status"] == "HEURISTIC"
    assert crispr["not_genome_wide"] is True
    assert crispr["spacer_is_not_target"] is True
    assert "PAM-containing" in crispr["input_note"]


def test_metric_explanations_never_turn_missing_into_zero() -> None:
    rec = explain("gravy", {"status": "UNAVAILABLE"})
    assert rec.plain_meaning
    assert "0" not in rec.title or "Unavailable" in rec.title
    assert "Unavailable" in rec.title
    entropy = explain("entropy", {"status": "COMPUTED", "shannon_entropy": float("nan")})
    assert "Unavailable" in entropy.title or "unavailable" in entropy.plain_meaning.lower()


class _Handle:
    def __init__(self, payload: dict) -> None:
        import json

        self._body = json.dumps(payload).encode("utf-8")
        self.status = 200

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return "https://rest.uniprot.org/uniprotkb/P38398.json"

    def __enter__(self) -> "_Handle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _uniprot_opener():
    payload = {
        "primaryAccession": "P38398",
        "uniProtkbId": "BRCA1_HUMAN",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Breast cancer type 1 susceptibility protein"}}
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

    def opener(request, timeout=None):
        return _Handle(payload)

    return opener
