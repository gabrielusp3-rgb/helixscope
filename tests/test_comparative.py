"""Comparacoes: estruturas (HTTP injetado), variantes, guias, MSA, proteinas."""

from __future__ import annotations

import json
from pathlib import Path
from urllib.request import Request

import pytest

from modules import (
    comparative,
    evidence_workspace,
    msa,
    rcsb_alignment,
    variant_core,
)
from tests.test_rcsb_alignment import _FakeHandle

FIXTURES = Path(__file__).parent / "fixtures"


def test_compare_structures_remote_parses_official_fixture():
    complete = (FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8")
    ticket = "095be615-a8ad-4c33-8e9c-c7612fbf6c9f"

    def opener(request: Request, timeout: float = 0) -> _FakeHandle:
        if request.get_method() == "POST":
            return _FakeHandle(json.dumps({"ticket": ticket}), request.full_url)
        return _FakeHandle(complete, request.full_url)

    result = comparative.compare_structures_remote(
        reference_entry="8HSK",
        target_entry="8HSF",
        method="fatcat-rigid",
        urlopen_fn=opener,
        sleep_fn=lambda _s: None,
    )
    alignment = result["alignment"]
    assert alignment["n_aligned_residue_pairs"] == 20
    assert alignment["rmsd_global_angstrom"] == pytest.approx(0.99)
    assert result["engine_location"] == "remote"
    assert result["export"]["confidence_score"] is None
    assert result["superposition"] is None


def test_compare_variants_does_not_rank():
    va = variant_core.build_variant(
        text="17 43093557 C G", assembly="GRCh38.p14", source=variant_core.SOURCE_USER
    )
    vb = variant_core.build_variant(
        text="17 43094500 T C", assembly="GRCh38.p14", source=variant_core.SOURCE_USER
    )
    fake_a = {
        "variant": va,
        "layers": {
            "consequences": {"status": "AVAILABLE", "data": {"most_severe_consequence": "missense_variant", "transcript_count": 50}},
            "protein_mapping": {"status": "AVAILABLE", "data": {"count": 1, "transcripts": [{"hgvsp": "p.Met658Ile", "amino_acids": "M/I", "protein_start": 658}]}},
            "clinical_evidence": {"status": "AVAILABLE", "data": {"status": "AVAILABLE"}},
            "domains": {"status": "SKIPPED"},
            "residue_properties": {"data": {"reference_amino_acid": "M", "variant_amino_acid": "I", "stability_claim": None}},
        },
    }
    fake_b = {
        "variant": vb,
        "layers": {
            "consequences": {"status": "AVAILABLE", "data": {"most_severe_consequence": "missense_variant", "transcript_count": 12}},
            "protein_mapping": {"status": "AVAILABLE", "data": {"count": 1, "transcripts": [{"hgvsp": "p.Ala344Gly", "amino_acids": "A/G", "protein_start": 344}]}},
            "clinical_evidence": {"status": "AVAILABLE", "data": {"status": "AVAILABLE"}},
            "domains": {"status": "SKIPPED"},
            "residue_properties": {"data": {"reference_amino_acid": "A", "variant_amino_acid": "G", "stability_claim": None}},
        },
    }
    compared = comparative.compare_variants(fake_a, fake_b)
    assert compared["ranking"] is None
    assert compared["side_by_side"]["a"]["hgvsp"] == "p.Met658Ile"
    assert compared["side_by_side"]["b"]["hgvsp"] == "p.Ala344Gly"
    assert compared["mutant_structure"]["status"] == "UNAVAILABLE"
    assert compared["side_by_side"]["a"]["helixscope_effect"] is None


def test_guide_comparison_refuses_genome_wide_and_mixed_scope():
    guides = (
        {"guide_sequence": "ACGTACGTACGTACGTACGT", "pam": "NGG"},
        {"guide_sequence": "TGCATGCATGCATGCATGCA", "pam": "NGG"},
    )
    search_gw = {"status": "COMPLETED", "genome_wide": True, "hits": [1, 2]}
    search_test = {
        "status": "COMPLETED_TEST_REFERENCE",
        "genome_wide": False,
        "hits": [1],
        "assembly_id": "HELIXSCOPE_TEST_REF",
        "guide_specificity": {"mit_sguide": {"score": 31.3, "genome_wide": False}},
    }
    mixed = comparative.compare_guides(*guides, search_a=search_gw, search_b=search_test)
    assert mixed["specificity_comparable"] is False
    assert mixed["genome_wide_comparison"] is False
    same = comparative.compare_guides(*guides, search_a=search_test, search_b=dict(search_test))
    assert same["specificity_comparable"] is True
    assert same["side_by_side"]["a"]["test_only"] is True


def test_msa_column_structure_roundtrip_and_group_pattern():
    fasta = ">leaf_a\nAC-DEF\n>leaf_b\nAC-DEF\n>outgroup\nAG-DKF\n"
    result = msa.import_prealigned_fasta(fasta)
    parsed = {
        "atoms": [
            {
                "group": "ATOM",
                "atom_name": "CA",
                "label_asym_id": "A",
                "label_seq_id": 1,
                "comp_id": "ALA",
                "x": 0,
                "y": 0,
                "z": 0,
            },
            {
                "group": "ATOM",
                "atom_name": "CA",
                "label_asym_id": "A",
                "label_seq_id": 2,
                "comp_id": "CYS",
                "x": 1,
                "y": 0,
                "z": 0,
            },
        ]
    }
    col0 = comparative.msa_column_to_structure(
        msa_result=result,
        column=0,
        member_id="leaf_a",
        parsed_structure=parsed,
        chain="A",
    )
    assert col0["status"] == "AVAILABLE"
    assert col0["msa_residue"] == "A"
    assert col0["conservation_note"].startswith("MSA conservation")
    back = comparative.structure_residue_to_msa_column(
        msa_result=result, member_id="leaf_a", ungapped_index=0
    )
    assert back["column"] == 0
    grouped = comparative.group_associated_positions(
        msa_result=result, group_ids=["leaf_a", "leaf_b"]
    )
    assert grouped["not_dn_ds"] is True
    assert grouped["not_positive_selection"] is True
    columns = {row["column"] for row in grouped["positions"]}
    assert 1 in columns
    assert "not an adaptive mutation" in grouped["statement"].lower()


def test_compare_proteins_does_not_guess_domains():
    seq_a = "MKTAYIAKQRQISFVKSHFSRQLE"
    seq_b = "MKTAYIAKQRQISFVKSHFSRQLE"
    result = comparative.compare_proteins(sequence_a=seq_a, sequence_b=seq_b)
    assert result["identity"]["ungapped_identity_pct"] == pytest.approx(100.0)
    assert result["domains"]["a"]["status"] == "UNAVAILABLE"
    assert "not guessed" in result["domains"]["a"]["reason"].lower()
    with pytest.raises(comparative.ComparativeError):
        comparative.compare_proteins(sequence_a="", sequence_b=seq_b)


def test_load_parsed_rejects_hostile_identifier():
    loaded = comparative.load_parsed_structure_for_entry("https://evil.example/1.cif")
    assert loaded["status"] == "INVALID_INPUT"
    loaded2 = comparative.load_parsed_structure_for_entry("BRCA1")
    assert loaded2["status"] == "INVALID_INPUT"
