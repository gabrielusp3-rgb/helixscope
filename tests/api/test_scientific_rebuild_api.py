"""API completeness for non-3D domain pipelines. No 3D assertions."""

from __future__ import annotations

from fastapi.testclient import TestClient

from helixscope_core import analyze_dna, analyze_rna, analyze_protein
from helixscope_core.crispr import analyze_crispr
from helixscope_core.motif import analyze_motif
from helixscope_api.serialization import to_transport


INTERNAL_KEYS = frozenset({"__debug__", "stack", "traceback"})


def _lost_fields(core: dict, api: dict) -> set[str]:
    core_keys = {key for key in core if key not in INTERNAL_KEYS}
    api_keys = set(api)
    return core_keys - api_keys


def test_dna_api_does_not_drop_core_domain_fields(client: TestClient) -> None:
    seq = "ATGC" * 40
    kwargs = {
        "include_orfs": True,
        "include_restriction": True,
        "include_profiles": True,
        "include_cpg": True,
        "include_kmers": True,
        "include_santalucia": True,
        "profile_window": 50,
        "profile_step": 25,
        "kmer_k": 3,
    }
    core = to_transport(analyze_dna(seq, **kwargs))
    response = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": seq, "include_explanation": True, **kwargs},
    )
    assert response.status_code == 200, response.text
    api = response.json()["result"]
    lost = _lost_fields(core, api)
    assert lost == set(), lost
    assert api["metric_explanations"]["gc_content"]["plain_meaning"]
    assert api["identity_scope"] == "sequence_derived"


def test_rna_api_keeps_codon_and_profiles(client: TestClient) -> None:
    seq = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA"
    core = to_transport(
        analyze_rna(seq, include_codon_metrics=True, include_profiles=True, include_kmers=True)
    )
    response = client.post(
        "/api/v1/rna/analyze",
        json={
            "sequence": seq,
            "include_codon_metrics": True,
            "include_profiles": True,
            "include_kmers": True,
            "include_explanation": True,
            "fold": False,
        },
    )
    assert response.status_code == 200, response.text
    api = response.json()["result"]
    lost = _lost_fields(core, api) - {"explanation", "metric_explanations"}
    # explanation is API-attached; Core may also omit it
    lost -= {"explanation", "metric_explanations"}
    assert lost == set(), lost
    assert "codon_metrics_scope" in api
    assert "cai" in api


def test_protein_api_keeps_physchem_fields(client: TestClient) -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    core = to_transport(analyze_protein(seq))
    response = client.post(
        "/api/v1/protein/analyze",
        json={"sequence": seq, "include_explanation": True},
    )
    assert response.status_code == 200, response.text
    api = response.json()["result"]
    lost = _lost_fields(core, api) - {"explanation", "metric_explanations"}
    assert lost == set(), lost


def test_protein_uniprot_rejects_url(client: TestClient) -> None:
    response = client.post(
        "/api/v1/protein/uniprot",
        json={"accession": "https://evil.example"},
    )
    assert response.status_code == 400
    assert response.json()["error"]["code"] == "INVALID_INPUT"


def test_motif_and_crispr_api_match_core_hits(client: TestClient) -> None:
    motif_core = analyze_motif("ACGTACGTACGT", "ACGT")
    motif_api = client.post(
        "/api/v1/motif/search",
        json={"sequence": "ACGTACGTACGT", "pattern": "ACGT", "molecule": "DNA", "include_explanation": True},
    )
    assert motif_api.status_code == 200, motif_api.text
    assert motif_api.json()["result"]["hits"] == to_transport(motif_core)["hits"]
    locus = "ACGTACGTACGTACGTACGTAGG"
    core = analyze_crispr(locus, "SpCas9")
    api = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": locus, "cas_system": "SpCas9", "include_explanation": True},
    )
    assert api.status_code == 200, api.text
    body = api.json()["result"]
    assert body["guides"][0]["guide_sequence"] == core["guides"][0]["guide_sequence"]
    assert body["input_note"]
    assert body["metric_explanations"]


def test_capabilities_include_remote_providers_without_secrets(client: TestClient) -> None:
    response = client.get("/api/v1/system/capabilities")
    assert response.status_code == 200
    providers = response.json()["result"]["remote_providers"]
    dumped = str(providers)
    assert "NCBI E-utilities" in dumped
    for row in providers:
        assert row["status"] != "REMOTE_VALIDATED"
        assert row.get("api_key") in {None, False, "", True} or "api_key" not in row
        assert not str(row.get("api_key") or "").startswith("secret")
