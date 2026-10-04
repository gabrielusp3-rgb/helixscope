"""FastAPI pasted-sequence recovery: FASTA headers vs invalid alphabet."""

from __future__ import annotations

from fastapi.testclient import TestClient

from helixscope_core import analyze_dna, analyze_protein


INSULIN_FASTA = (
    ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n"
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n"
)


def test_api_dna_fasta_matches_core(client: TestClient) -> None:
    fasta = ">dna_fixture\nATGCATGCATGC\n"
    core = analyze_dna(fasta)
    response = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": fasta, "include_explanation": True, "include_profiles": True},
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["status"] == "COMPUTED"
    assert body["length"] == core["length"] == 12
    assert body["gc_content"] == core["gc_content"]
    assert body["gc_skew"] == core["gc_skew"]
    assert body["at_skew"] == core["at_skew"]
    assert "explanation" in body


def test_api_dna_invalid_is_error_not_fake_metrics(client: TestClient) -> None:
    response = client.post(
        "/api/v1/dna/analyze",
        json={"sequence": "NNNN", "include_explanation": True},
    )
    assert response.status_code == 200
    body = response.json()["result"]
    assert body["status"] == "ERROR"
    assert "explanation" not in body
    assert body["validation"]["invalid_chars"]


def test_api_protein_fasta_matches_core(client: TestClient) -> None:
    core = analyze_protein(INSULIN_FASTA)
    response = client.post(
        "/api/v1/protein/analyze",
        json={"sequence": INSULIN_FASTA, "include_explanation": False},
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["length"] == core["length"]
    assert body["isoelectric_point"] == core["isoelectric_point"]
    assert "|" not in str(body.get("sequence") or "")


def test_api_protein_invalid_lists_residues_not_fasta_punctuation(client: TestClient) -> None:
    response = client.post(
        "/api/v1/protein/analyze",
        json={"sequence": ">hdr extra=1\nMKTAYIOUAK\n", "include_explanation": False},
    )
    assert response.status_code == 400
    message = response.json()["error"]["message"]
    assert "O" in message and "U" in message
    assert ">" not in message
    assert "=" not in message


def test_api_crispr_short_spacer_zero_guides_not_error(client: TestClient) -> None:
    response = client.post(
        "/api/v1/crispr/guides",
        json={"sequence": "ACGTACGTACGTACGTACGT", "cas_system": "SpCas9", "include_explanation": True},
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["n_guides"] == 0
    assert body["target_length"] == 20
    assert body["status"] == "HEURISTIC"


def test_api_crispr_fasta_target_has_guides(client: TestClient) -> None:
    response = client.post(
        "/api/v1/crispr/guides",
        json={
            "sequence": ">locus\nACGTACGTACGTACGTACGTAGG\n",
            "cas_system": "SpCas9",
            "include_explanation": True,
        },
    )
    assert response.status_code == 200, response.text
    body = response.json()["result"]
    assert body["n_guides"] >= 1
    assert body["guides"][0]["pam_sequence"]
