"""Prompt 1 reconnection: FastAPI contracts, CRISPR tokens, Compare field pass-through."""

from __future__ import annotations

import time
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from helixscope_api.job_handlers import run_compare_remote_job, run_usalign_job
from helixscope_core.msa import import_alignment
from modules.crispr import CAS_SYSTEMS

PREALIGNED = """>seq_a
ACGTACGTACGTACGTACGT
>seq_b
ACGTACGTACGTACGTTTTT
>seq_c
TTTTACGTACGTACGTACGT
>seq_d
ACGTACGTAAAAACGTACGT
"""
SPCAS9_LOCUS = "ACGTACGTACGTACGTACGTAGG"
ANNOTATED_TABLE = "Humano      [ A G T C ]\nChimpanze   [ A G T - ] <- comment\n"
ORF_RNA = "AUG" + "GCA" * 48 + "UAA"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"

SYSTEM_FIXTURES = {
    "SpCas9": SPCAS9_LOCUS,
    "SaCas9": "ACGTACGTACGTACGTACGTA" + "AAGAGT",
    "AsCas12a": "TTTA" + "ACGTACGTACGTACGTACGTACG",
    "Cas12b": "TTA" + "ACGTACGTACGTACGTACGT" + "CCCC",
    "Cas13": "ACGUACGUACGUACGUACGUACGUACGUAC",
    "CBE (SpCas9)": SPCAS9_LOCUS,
    "ABE (SpCas9)": SPCAS9_LOCUS,
    "Prime Editor (SpCas9 nickase)": SPCAS9_LOCUS,
}


def _poll(client: TestClient, job_id: str, *, attempts: int = 120, delay: float = 0.1) -> dict:
    snapshot: dict = {}
    for _ in range(attempts):
        got = client.get(f"/api/v1/jobs/{job_id}")
        assert got.status_code == 200, got.text
        snapshot = got.json()
        if snapshot.get("execution_status") in {"COMPLETED", "FAILED", "TIMEOUT", "RESOURCE_LIMIT"}:
            return snapshot
        time.sleep(delay)
    return snapshot


def test_crispr_systems_registry_matches_core(client: TestClient) -> None:
    response = client.get("/api/v1/crispr/systems")
    assert response.status_code == 200, response.text
    systems = response.json()["result"]["systems"]
    keys = [row["canonical_key"] for row in systems]
    assert keys == list(CAS_SYSTEMS)
    assert "Cas12a" not in keys
    assert "CBE" not in keys
    cas13 = next(row for row in systems if row["canonical_key"] == "Cas13")
    assert cas13["target_molecule"] == "RNA"


@pytest.mark.parametrize("token", ["Cas12a", "CBE", "ABE", "Prime Editor"])
def test_informal_cas_tokens_are_invalid_input(client: TestClient, token: str) -> None:
    response = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": SPCAS9_LOCUS, "cas_system": token, "include_explanation": False},
    )
    assert response.status_code == 400, response.text
    assert response.json()["error"]["code"] == "INVALID_INPUT"


@pytest.mark.parametrize("system", list(CAS_SYSTEMS))
def test_every_cas_system_token_is_accepted(client: TestClient, system: str) -> None:
    sequence = SYSTEM_FIXTURES[system]
    response = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": sequence, "cas_system": system, "include_explanation": False},
    )
    assert response.status_code == 200, response.text
    result = response.json()["result"]
    assert result["cas_system"] == system
    assert result["cas_system"] != "SpCas9" or system == "SpCas9"
    for guide in result.get("guides") or []:
        assert guide.get("ruleset2_score") is None
        assert guide.get("deephf_score") is None


def test_cas13_dna_is_rejected_and_rna_is_accepted(client: TestClient) -> None:
    dna = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": "ACGTACGTACGTACGTACGTACGTACGTAC", "cas_system": "Cas13"},
    )
    assert dna.status_code == 400
    rna = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": SYSTEM_FIXTURES["Cas13"], "cas_system": "Cas13", "include_explanation": False},
    )
    assert rna.status_code == 200, rna.text
    assert rna.json()["result"]["cas_system"] == "Cas13"
    assert rna.json()["result"]["target_molecule"] == "RNA"


def test_msa_formats_and_annotated_table_error(client: TestClient) -> None:
    fasta = client.post("/api/v1/msa/prealigned", json={"fasta": PREALIGNED, "format": "fasta"})
    assert fasta.status_code == 200, fasta.text
    assert fasta.json()["result"]["n_sequences"] == 4
    table = client.post("/api/v1/msa/prealigned", json={"fasta": ANNOTATED_TABLE, "format": "auto"})
    assert table.status_code == 400, table.text
    body = table.json()
    assert body["error"]["code"] == "INVALID_FORMAT"
    assert "not a recognized alignment format" in body["error"]["message"]
    clustal = """CLUSTAL W (1.83) multiple sequence alignment

seq_a      ACGTACGTACGTACGTACGT
seq_b      ACGTACGTACGTACGTTTTT
seq_c      TTTTACGTACGTACGTACGT
seq_d      ACGTACGTAAAAACGTACGT
"""
    imported = client.post("/api/v1/msa/prealigned", json={"fasta": clustal, "format": "clustal"})
    assert imported.status_code == 200, imported.text
    assert imported.json()["result"]["input_format"] == "clustal"


def test_msa_to_phylogeny_preserves_hash(client: TestClient) -> None:
    msa = client.post("/api/v1/msa/prealigned", json={"fasta": PREALIGNED, "format": "fasta"})
    alignment_hash = msa.json()["result"]["alignment_hash"]
    nj = client.post(
        "/api/v1/phylogeny/infer",
        json={
            "fasta": PREALIGNED,
            "method": "neighbor_joining",
            "distance_model": "p_distance",
            "alignment_hash": alignment_hash,
        },
    )
    assert nj.status_code == 200, nj.text
    result = nj.json()["result"]
    assert result["source_msa_hash"] == alignment_hash
    assert result["alignment_hash"] == alignment_hash
    assert result["method"] == "neighbor_joining"
    assert result["newick"]
    upgma = client.post(
        "/api/v1/phylogeny/infer",
        json={"fasta": PREALIGNED, "method": "upgma", "distance_model": "p_distance", "alignment_hash": alignment_hash},
    )
    assert upgma.status_code == 200, upgma.text
    assert upgma.json()["result"]["method"] == "upgma"
    assert upgma.json()["result"]["source_msa_hash"] == alignment_hash


def test_iqtree_and_fasttree_jobs_or_classified_absence(client: TestClient) -> None:
    alignment_hash = import_alignment(PREALIGNED, fmt="fasta")["alignment_hash"]
    for method in ("iqtree_ml", "fasttree_ml"):
        submitted = client.post(
            "/api/v1/phylogeny/jobs",
            json={"fasta": PREALIGNED, "method": method, "alignment_hash": alignment_hash},
        )
        if submitted.status_code in {503, 409}:
            assert submitted.json()["error"]["code"] in {"ENGINE_NOT_INSTALLED", "RESOURCE_LIMIT"}
            continue
        assert submitted.status_code == 202, submitted.text
        snapshot = _poll(client, submitted.json()["job_id"], attempts=200, delay=0.25)
        status = snapshot.get("execution_status")
        if status == "COMPLETED":
            inner = snapshot["result"]["result"]
            assert method.replace("_ml", "") in str(inner.get("method") or method).lower() or inner.get("method") == method
            assert inner.get("source_msa_hash") == alignment_hash or inner.get("alignment_hash") == alignment_hash
            assert inner.get("newick")
        else:
            assert status in {"FAILED", "TIMEOUT", "RESOURCE_LIMIT"}
            code = str((snapshot.get("result") or {}).get("error", {}).get("code") or snapshot.get("scientific_status") or "")
            assert code or status


def test_taxonomy_endpoint_preserves_topology(client: TestClient, monkeypatch: pytest.MonkeyPatch) -> None:
    tree = client.post(
        "/api/v1/phylogeny/infer",
        json={"fasta": PREALIGNED, "method": "neighbor_joining", "distance_model": "p_distance"},
    ).json()["result"]
    original_newick = tree["newick"]
    original_hash = tree["tree_hash"]

    def fake_attach(tree_result, email, api_key=None):
        updated = dict(tree_result)
        updated["newick"] = "CORRUPTED;"
        updated["tree_hash"] = "nope"
        updated["taxonomy_layer"] = {"status": "COMPUTED"}
        return updated

    monkeypatch.setattr("modules.taxonomy.attach_ncbi_taxonomy", fake_attach)
    response = client.post(
        "/api/v1/phylogeny/taxonomy",
        json={"tree": tree, "email": "helixscope@example.org", "declared_organisms": {"seq_a": "Homo sapiens"}},
    )
    assert response.status_code == 200, response.text
    payload = response.json()["result"]
    assert payload["newick"] == original_newick
    assert payload["tree_hash"] == original_hash


def test_translate_nucleic_alignment_and_rna_composition(client: TestClient) -> None:
    silent = client.post(
        "/api/v1/alignment/pairwise",
        json={"seq1": "ATGCCCTAA", "seq2": "ATGCCCTAG", "mode": "global", "translate_nucleic": False},
    )
    assert silent.status_code == 200, silent.text
    assert silent.json()["result"]["translation"]["applied"] is False
    translated = client.post(
        "/api/v1/alignment/pairwise",
        json={"seq1": "ATGCCCGGGTAA", "seq2": "ATGCCCGGGTAG", "mode": "global", "translate_nucleic": True},
    )
    assert translated.status_code == 200, translated.text
    assert translated.json()["result"]["translation"]["applied"] is True
    rna = client.post("/api/v1/rna/translate-coding", json={"sequence": ORF_RNA})
    assert rna.status_code == 200, rna.text
    assert rna.json()["result"]["status"] in {"COMPUTED", "UNAVAILABLE", "ERROR"}
    random_rna = client.post("/api/v1/rna/translate-coding", json={"sequence": "ACGUACGUACGU"})
    assert random_rna.status_code == 200
    assert random_rna.json()["result"]["status"] in {"UNAVAILABLE", "ERROR"}


def test_protein_structure_capabilities_load_and_sasa(client: TestClient) -> None:
    caps = client.get("/api/v1/protein/structure/capabilities")
    assert caps.status_code == 200, caps.text
    result = caps.json()["result"]
    assert "dssp_local" in result
    assert "dssp_pdb_redo" in result
    assert "sasa" in result
    assert "edtsurf" in result
    loaded = client.post("/api/v1/protein/structure/load", json={"structure_id": "1CRN", "bundled": True})
    assert loaded.status_code == 200, loaded.text
    assert loaded.json()["result"]["structure_id"] in {"1CRN", "1crn"} or "1CRN" in str(loaded.json()["result"])
    analyzed = client.post(
        "/api/v1/protein/structure/analyze",
        json={"structure_id": "1CRN", "run_sasa": True, "run_dssp": False, "run_edtsurf": False},
    )
    assert analyzed.status_code == 200, analyzed.text
    sasa = analyzed.json()["result"]["sasa"]
    assert sasa["status"] in {"COMPUTED", "EXPERIMENTAL", "UNAVAILABLE", "ERROR"}
    search = client.post(
        "/api/v1/protein/structure/search",
        json={"sequence": CRAMBIN, "molecule": "PROTEIN", "pdb_id": "1CRN", "search_pdb_by_sequence": False},
    )
    assert search.status_code == 200, search.text
    hits = search.json()["result"].get("hits") or []
    assert any(str(hit.get("structure_id") or "").upper() == "1CRN" for hit in hits if isinstance(hit, dict))


def test_usalign_schema_rejects_non_bundled_and_handler_keeps_user_ids(
    client: TestClient, monkeypatch: pytest.MonkeyPatch
) -> None:
    rejected = client.post(
        "/api/v1/compare/jobs/usalign",
        json={"reference_entry": "8HSK", "target_entry": "8HSF", "reference_chain": "B", "target_chain": "C"},
    )
    assert rejected.status_code == 422
    captured: dict = {}

    def fake_compare(**kwargs):
        captured.update(kwargs)
        return {
            "mode": "structures",
            "reference_entry": kwargs["reference_entry"],
            "target_entry": kwargs["target_entry"],
            "original_coordinates_preserved": True,
        }

    monkeypatch.setattr("helixscope_api.job_handlers.read_bundled_mmcif", lambda name: "data_X")
    monkeypatch.setattr("helixscope_api.job_handlers.compare_structures_usalign", fake_compare)
    out = run_usalign_job(
        {
            "reference_entry": "1RNA",
            "target_entry": "1CRN",
            "reference_chain": "B",
            "target_chain": "C",
        },
        SimpleNamespace(job_id="t"),
    )
    assert captured["reference_entry"] == "1RNA"
    assert captured["target_entry"] == "1CRN"
    assert captured["reference_chain"] == "B"
    assert captured["target_chain"] == "C"
    assert out["reference_entry"] == "1RNA"


def test_rcsb_compare_uses_selected_chains_and_method(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}

    def fake_compare(**kwargs):
        captured.update(kwargs)
        return {"mode": "structures", "method": kwargs.get("method")}

    monkeypatch.setattr("helixscope_api.job_handlers.compare_structures_remote", fake_compare)
    run_compare_remote_job(
        {
            "reference_entry": "8HSK",
            "target_entry": "8HSF",
            "reference_chain": "B",
            "target_chain": "C",
            "method": "fatcat-rigid",
            "fetch_coordinates": False,
        },
        SimpleNamespace(job_id="t"),
    )
    assert captured["reference_entry"] == "8HSK"
    assert captured["target_entry"] == "8HSF"
    assert captured["reference_chain"] == "B"
    assert captured["target_chain"] == "C"
    assert captured["method"] == "fatcat-rigid"
    assert captured["method"] != "tm-align"


def test_compare_six_modes(client: TestClient) -> None:
    modes = client.get("/api/v1/compare/modes")
    assert modes.status_code == 200, modes.text
    ids = {row["id"] for row in modes.json()["result"]["modes"]}
    assert ids == {"structures", "variants", "proteins", "guides", "evolution", "evidence"}
    methods = client.get("/api/v1/compare/methods")
    assert methods.status_code == 200
    api_names = [row["api_name"] for row in methods.json()["result"]["methods"]]
    assert "tm-align" in api_names
    assert "fatcat-rigid" in api_names
    variants = client.post(
        "/api/v1/compare/variants",
        json={
            "text_a": "17 43093557 C G",
            "assembly_a": "GRCh38.p14",
            "text_b": "17 43093558 A T",
            "assembly_b": "GRCh38.p14",
        },
    )
    assert variants.status_code == 200, variants.text
    assert variants.json()["result"]["ranking"] is None
    proteins = client.post(
        "/api/v1/compare/proteins",
        json={"sequence_a": "MKTAYIAKQR", "sequence_b": "MKTAYIAKQR", "identifier_a": "a", "identifier_b": "b"},
    )
    assert proteins.status_code == 200, proteins.text
    guides = client.post(
        "/api/v1/compare/guides",
        json={
            "guide_a": {"guide_sequence": "ACGTACGTACGTACGTACGT", "pam_sequence": "AGG"},
            "guide_b": {"guide_sequence": "TGCATGCATGCATGCATGCA", "pam_sequence": "TGG"},
        },
    )
    assert guides.status_code == 200, guides.text
    evo = client.post("/api/v1/compare/evolution", json={"fasta": PREALIGNED, "column": 0})
    assert evo.status_code == 200, evo.text
    evidence = client.post(
        "/api/v1/compare/evidence",
        json={
            "text_a": "17 43093557 C G",
            "assembly_a": "GRCh38.p14",
            "text_b": "17 43093558 A T",
            "assembly_b": "GRCh38.p14",
        },
    )
    assert evidence.status_code == 200, evidence.text
    assert evidence.json()["result"]["confidence_score"] is None


def test_references_admission_and_install_confirm(client: TestClient) -> None:
    test_ref = client.post("/api/v1/references/admission", json={"assembly_id": "HELIXSCOPE_TEST_REF"})
    assert test_ref.status_code == 200, test_ref.text
    assert test_ref.json()["result"]["not_a_public_assembly"] is True
    grch = client.post("/api/v1/references/admission", json={"assembly_id": "GRCh38.p14"})
    assert grch.status_code == 200, grch.text
    unconfirmed = client.post("/api/v1/references/install", json={"assembly_id": "GRCh38.p14", "confirm": False})
    assert unconfirmed.status_code == 400
    if grch.json()["result"].get("decision") != "PROCEED":
        blocked = client.post("/api/v1/references/install", json={"assembly_id": "GRCh38.p14", "confirm": True})
        assert blocked.status_code in {409, 400}, blocked.text
        assert blocked.json()["error"]["code"] in {"RESOURCE_LIMIT", "INVALID_INPUT"}
