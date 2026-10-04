"""Direct helixscope_core vs FastAPI transport parity."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from helixscope_core import (
    analyze_dna,
    analyze_protein,
    analyze_rna,
    design_guides,
    evidence_conflict_pack,
    explain_result,
    identify_variant,
    infer_tree,
    pairwise_align,
    search_motifs,
    transcribe_dna,
)
from helixscope_core.compare import parse_alignment_payload
from helixscope_core.crispr import model_availability
from helixscope_core.dna import gc_content
from helixscope_core.evidence import evidence_item
from helixscope_core.msa import import_prealigned_fasta
from helixscope_core.phylogeny import DISTANCE_P, METHOD_NJ
from helixscope_core.rna import codon_adaptation_index_report
from helixscope_core.serialize import provenance_compare_fields
from helixscope_core.structure import classify_mapping_status
from helixscope_api.serialization import assert_strict_json, to_transport

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures"
PREALIGNED = """>seq_a
ACGTACGTACGTACGTACGT
>seq_b
ACGTACGTACGTACGTTTTT
>seq_c
TTTTACGTACGTACGTACGT
>seq_d
ACGTACGTAAAAACGTACGT
"""


def _result(client: TestClient, method: str, path: str, payload: dict | None = None) -> dict:
    if method == "GET":
        response = client.get(path)
    else:
        response = client.post(path, json=payload)
    assert response.status_code == 200, response.text
    body = response.json()
    assert_strict_json(body)
    return body


def test_dna_parity_and_nan_not_zero(client: TestClient) -> None:
    seq = "ATGCATGCATGC"
    core = analyze_dna(seq)
    body = _result(client, "POST", "/api/v1/dna/analyze", {"sequence": seq, "include_explanation": False})
    api = body["result"]
    transported = to_transport(core)
    assert api["status"] == core["status"] == "COMPUTED"
    assert api["length"] == core["length"] == 12
    assert api["sequence_hash"] == core["sequence_hash"]
    assert api["gc_content"] == transported["gc_content"] == pytest.approx(50.0)
    nan_body = _result(
        client,
        "POST",
        "/api/v1/dna/analyze",
        {"sequence": "AT" * 80, "include_profiles": True, "include_explanation": False},
    )
    assert nan_body["result"]["gc_content"] == pytest.approx(0.0)
    profile = nan_body["result"]["windowed_profiles"][0]
    assert profile["gc_percent"] == pytest.approx(0.0)
    assert isinstance(profile["gc_skew"], dict)
    assert profile["gc_skew"]["value"] is None
    assert profile["gc_skew"]["value_state"] == "NAN"
    assert math.isnan(gc_content("NNNN"))
    assert to_transport(gc_content("NNNN")) == {"value": None, "value_state": "NAN"}


def test_dna_50k(client: TestClient) -> None:
    seq = "ATGC" * 12_500
    core = analyze_dna(seq)
    body = _result(client, "POST", "/api/v1/dna/analyze", {"sequence": seq, "include_explanation": False})
    assert body["result"]["length"] == core["length"] == 50_000
    assert body["result"]["gc_content"] == pytest.approx(50.0)


def test_dna_optional_core_fields_and_extras(client: TestClient) -> None:
    seq = "ATGCATGCATGC"
    core = analyze_dna(seq)
    body = _result(
        client,
        "POST",
        "/api/v1/dna/analyze",
        {
            "sequence": seq,
            "include_explanation": False,
            "include_cpg": True,
            "include_kmers": True,
            "include_santalucia": True,
        },
    )
    api = body["result"]
    assert api["reverse_complement"] == core["reverse_complement"]
    assert api["composition"]["A"]["count"] == core["composition"]["A"]["count"]
    assert "cpg_islands" in api
    assert api["kmer_summary"]["k"] == 3
    assert api["santalucia_tm"]["status"] in {"COMPUTED", "UNAVAILABLE", "INVALID_INPUT"}
    default_body = _result(
        client,
        "POST",
        "/api/v1/dna/analyze",
        {"sequence": seq, "include_explanation": False},
    )
    assert "cpg_islands" not in default_body["result"]
    assert "kmer_summary" not in default_body["result"]
    assert "santalucia_tm" not in default_body["result"]


def test_rna_transcribe_and_cai(client: TestClient) -> None:
    trans = _result(client, "POST", "/api/v1/rna/transcribe", {"sequence": "ATGC"})
    assert trans["result"]["rna"] == transcribe_dna("ATGC") == "AUGC"
    coding = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA"
    core = analyze_rna(coding)
    body = _result(
        client,
        "POST",
        "/api/v1/rna/analyze",
        {"sequence": coding, "fold": False, "include_codon_metrics": True, "include_explanation": False},
    )
    assert body["result"]["status"] == core["status"]
    assert body["result"]["length"] == core["length"]
    cai = body["result"]["cai"]
    core_cai = codon_adaptation_index_report(coding)
    assert cai["status"] == core_cai["status"] == "HEURISTIC"
    assert cai["value"] == pytest.approx(core_cai["value"], abs=1e-4)
    enc = body["result"]["enc"]
    if isinstance(enc.get("value"), dict):
        assert enc["value"]["value"] is None
        assert enc["value"]["value_state"] == "NAN"
        assert enc["value"] != 0


def test_protein_parity(client: TestClient) -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    core = analyze_protein(seq)
    body = _result(client, "POST", "/api/v1/protein/analyze", {"sequence": seq, "include_explanation": False})
    api = body["result"]
    assert api["status"] == core["status"]
    assert api["length"] == core["length"]
    assert api["molecular_weight_kda"] == pytest.approx(core["molecular_weight_kda"], abs=1e-6)
    assert api["isoelectric_point"] == pytest.approx(core["isoelectric_point"], abs=1e-2)
    assert api["gravy"] == pytest.approx(core["gravy"], abs=1e-3)
    assert api["ss_prediction"]["available"] is False


def test_alignment_parity(client: TestClient) -> None:
    core_g = pairwise_align("ACGTACGTACGT", "ACGTACGGACGT", mode="global")
    api_g = _result(
        client,
        "POST",
        "/api/v1/alignment/pairwise",
        {"seq1": "ACGTACGTACGT", "seq2": "ACGTACGGACGT", "mode": "global", "include_explanation": False},
    )["result"]
    assert api_g["method"] == core_g["method"] == "Needleman-Wunsch"
    assert api_g["aligned_seq1"] == core_g["aligned_seq1"]
    assert api_g["identity_pct"] == pytest.approx(core_g["identity_pct"])
    core_l = pairwise_align("ACGTACGTACGT", "ACGTACGGACGT", mode="local")
    api_l = _result(
        client,
        "POST",
        "/api/v1/alignment/pairwise",
        {"seq1": "ACGTACGTACGT", "seq2": "ACGTACGGACGT", "mode": "local", "include_explanation": False},
    )["result"]
    assert api_l["method"] == core_l["method"] == "Smith-Waterman"


def test_motif_parity(client: TestClient) -> None:
    core = search_motifs("ACGTACGTACGT", "ACGT", molecule="DNA")
    api = _result(
        client,
        "POST",
        "/api/v1/motif/search",
        {"sequence": "ACGTACGTACGT", "pattern": "ACGT", "molecule": "DNA", "include_explanation": False},
    )["result"]
    assert api["hits"] == to_transport(core)
    assert api["hit_meaning"].startswith("pattern occurrence")


def test_crispr_parity_and_none_scores(client: TestClient) -> None:
    seq = "ACGTACGTACGTACGTACGTAGG"
    core = design_guides(seq, "SpCas9")
    api = _result(
        client,
        "POST",
        "/api/v1/crispr/guides",
        {"sequence": seq, "cas_system": "SpCas9", "include_explanation": False},
    )["result"]
    assert api["guides"][0]["guide_sequence"] == core[0]["guide_sequence"]
    assert api["guides"][0]["position"] == core[0]["position"] == 1
    assert api["guides"][0].get("ruleset2_score") is None
    assert api["guides"][0].get("deephf_score") is None
    assert api["not_genome_wide"] is True
    avail = api["model_availability"]
    core_avail = model_availability()
    assert avail["on_target_azimuth"]["available"] is False
    assert core_avail["on_target_azimuth"]["available"] is False


def test_variant_and_evidence_parity(client: TestClient) -> None:
    core = identify_variant("17 43093557 C G", assembly="GRCh38.p14")
    api = _result(
        client,
        "POST",
        "/api/v1/variants/identify",
        {"text": "17 43093557 C G", "assembly": "GRCh38.p14"},
    )["result"]
    assert api["identity_hash"] == core["identity_hash"]
    assert api["position_0based"] == core["position_0based"]
    assert api["effect_status"] == "NOT_COMPUTED"
    assert api.get("helixscope_pathogenic") not in {True, "true", "True"}
    items = [
        evidence_item(
            field="Variant A most_severe_consequence",
            value="missense_variant",
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
        evidence_item(
            field="Variant B most_severe_consequence",
            value="synonymous_variant",
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
    ]
    core_pack = evidence_conflict_pack(items, kind="variants")
    pack = _result(
        client,
        "POST",
        "/api/v1/evidence/pack",
        {
            "kind": "variants",
            "items": [
                {
                    "field": item["field"],
                    "value": item["value"],
                    "source": item["source"],
                    "evidence_status": item["evidence_status"],
                    "retrieved_at_utc": item["retrieved_at_utc"],
                }
                for item in items
            ],
        },
    )["result"]
    assert pack["confidence_score"] is None
    assert core_pack["confidence_score"] is None
    assert pack["conflicts"]


def test_explain_parity(client: TestClient) -> None:
    values = {"pattern": "GAATTC", "n_hits": 2, "status": "COMPUTED"}
    core = explain_result("motif", values)
    api = _result(client, "POST", "/api/v1/explain", {"kind": "motif", "values": values})["result"]
    assert api["limitations"] == core["limitations"]
    assert "does not prove biological function" in api["limitations"].lower()


def test_msa_and_phylogeny_nj_parity(client: TestClient) -> None:
    core_msa = import_prealigned_fasta(PREALIGNED)
    api_msa = _result(client, "POST", "/api/v1/msa/prealigned", {"fasta": PREALIGNED})["result"]
    assert api_msa["alignment_hash"] == core_msa["alignment_hash"]
    assert api_msa["alignment_hash"] == "3b861c3b34ed5bbb6f01a27588fc50a07eed28fb557c9912a6278fe41086748a"
    core_tree = infer_tree(core_msa, method=METHOD_NJ, distance_model=DISTANCE_P)
    api_tree = _result(
        client,
        "POST",
        "/api/v1/phylogeny/infer",
        {"fasta": PREALIGNED, "method": "neighbor_joining", "distance_model": "p_distance"},
    )["result"]
    assert api_tree["status"] == core_tree["status"] == "COMPUTED"
    assert api_tree["method"] == core_tree["method"]
    assert {leaf["tree_id"] for leaf in api_tree["leaves"]} == {leaf["tree_id"] for leaf in core_tree["leaves"]}
    assert api_tree["not_taxonomic_tree"] is True
    assert api_tree["newick"]
    assert provenance_compare_fields(api_tree).get("status") == "COMPUTED"


def test_compare_and_mapping_uncertain(client: TestClient) -> None:
    payload = json.loads((FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8"))
    core = parse_alignment_payload(payload)
    api = _result(client, "POST", "/api/v1/compare/parse", {"payload": payload})["result"]
    assert api["n_aligned_residue_pairs"] == core["n_aligned_residue_pairs"] == 20
    assert api["rmsd_block0_angstrom"] == pytest.approx(core["rmsd_block0_angstrom"], abs=0.05)
    assert api["rmsd_global_angstrom"] == pytest.approx(core["rmsd_global_angstrom"], abs=0.05)
    assert api["rmsd_block0_angstrom"] != api["rmsd_global_angstrom"]
    core_map = classify_mapping_status(
        sequences_identical=True,
        polymer_index_mode=False,
        observed_coordinate_residues=0,
        polymer_length=0,
    )
    api_map = _result(
        client,
        "POST",
        "/api/v1/structures/mapping/classify",
        {
            "sequences_identical": True,
            "polymer_index_mode": False,
            "observed_coordinate_residues": 0,
            "polymer_length": 0,
        },
    )["result"]
    assert api_map["mapping_status"] == core_map["mapping_status"] == "UNCERTAIN"


def test_structure_fixture_has_coordinates_no_camera(client: TestClient) -> None:
    body = _result(client, "POST", "/api/v1/structures/fixture", {"structure_id": "1CRN"})
    result = body["result"]
    assert result["structure_id"] == "1CRN"
    assert result["kind"].lower() == "experimental"
    viewer = result["viewer"]
    assert viewer.get("coordinates")
    assert "camera" not in viewer
    assert "mediapipe" not in json.dumps(viewer).lower()
    dumped = json.dumps(body)
    assert "C:\\\\Users" not in dumped
    catalog = _result(client, "GET", "/api/v1/structures/catalog/experimental-complexes")
    pointers = json.dumps(catalog["result"])
    assert "4UN3" in pointers
    assert catalog["result"]["mapping_contract"]["structure_mapping"] is None
