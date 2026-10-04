"""Validacao remota opcional da RCSB Alignment API (8HSK vs 8HSF).

Nao e LIVE_VALIDATED local. Skip se a rede falhar. Golden esperado: o JSON
oficial da documentacao RCSB (tests/fixtures/rcsb_alignment_complete.json)
usa fatcat-rigid no mesmo par; a corrida live pode usar tm-align ou
fatcat-rigid e compara metodo, pares, RMSD e TM-score dentro da tolerancia.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from modules import (
    engine_validation,
    protein_structure,
    rcsb_alignment,
    structure_alignment,
    structure_superposition,
)

GOLDEN = Path(__file__).parent / "fixtures" / "rcsb_alignment_complete.json"


def test_rcsb_alignment_remote_8hsk_8hsf(tmp_path: Path):
    try:
        wrapped = rcsb_alignment.align_pairwise(
            reference_entry="8HSK",
            target_entry="8HSF",
            reference_chain="A",
            target_chain="A",
            method=rcsb_alignment.METHOD_FATCAT_RIGID,
        )
    except rcsb_alignment.AlignmentApiError as exc:
        pytest.skip(f"RCSB Alignment API unavailable: {exc.category}: {exc}")
    record = structure_alignment.parse_alignment_payload(
        wrapped["raw"],
        method=wrapped["method"],
        method_kind=wrapped["method_kind"],
        query=wrapped.get("query"),
    )
    golden = structure_alignment.parse_alignment_payload(json.loads(GOLDEN.read_text(encoding="utf-8")))
    assert record["method"] == golden["method"] == "fatcat-rigid"
    assert record["n_aligned_residue_pairs"] == golden["n_aligned_residue_pairs"] == 20
    assert record["atoms_fitted"] == "CA"
    if record["rmsd_block0_angstrom"] is not None and golden["rmsd_block0_angstrom"] is not None:
        delta = abs(float(record["rmsd_block0_angstrom"]) - float(golden["rmsd_block0_angstrom"]))
        assert delta <= 0.15 or delta / max(float(golden["rmsd_block0_angstrom"]), 1e-9) <= 0.1
    live_tm = record["tm_scores"]
    gold_tm = golden["tm_scores"]
    if live_tm and gold_tm:
        assert live_tm[0]["type"] == gold_tm[0]["type"]
        tm_delta = abs(float(live_tm[0]["value"]) - float(gold_tm[0]["value"]))
        assert tm_delta <= 0.05 or tm_delta / max(abs(float(gold_tm[0]["value"])), 1e-9) <= 0.15
    mapped = structure_alignment.lookup_correspondence(
        record, side="reference", asym_id="A", label_seq_id=2
    )
    gold_map = structure_alignment.lookup_correspondence(
        golden, side="reference", asym_id="A", label_seq_id=2
    )
    assert mapped["status"] == gold_map["status"] == "AVAILABLE"
    superposition_proof = {"status": "UNAVAILABLE"}
    try:
        ref_parsed = protein_structure.load_deposited_macromolecule(
            source="RCSB PDB", structure_id="8HSK"
        )
        tgt_parsed = protein_structure.load_deposited_macromolecule(
            source="RCSB PDB", structure_id="8HSF"
        )
        original_x = float((tgt_parsed.get("atoms") or [{}])[0].get("x"))
        bundle = structure_superposition.superposition_bundle(
            alignment=record,
            reference_parsed=ref_parsed,
            target_parsed=tgt_parsed,
            block_index=0,
        )
        after_x = float((tgt_parsed.get("atoms") or [{}])[0].get("x"))
        assert after_x == pytest.approx(original_x)
        assert bundle["original_target_untouched"] is True
        check_global = structure_superposition.rmsd_agrees(
            bundle["independent_rmsd"],
            record.get("rmsd_global_angstrom"),
        )
        check_block = structure_superposition.rmsd_agrees(
            bundle["independent_rmsd"],
            record.get("rmsd_block0_angstrom"),
        )
        assert check_global.get("comparable") is True
        assert check_global.get("agree") is True
        superposition_proof = {
            "status": "COMPUTED",
            "visual_status": bundle.get("visual_status"),
            "independent_rmsd_angstrom": bundle["independent_rmsd"].get("rmsd_angstrom"),
            "n_pairs_used": bundle["independent_rmsd"].get("n_pairs_used"),
            "api_rmsd_block0": record.get("rmsd_block0_angstrom"),
            "api_rmsd_global": record.get("rmsd_global_angstrom"),
            "agreement_vs_global": check_global,
            "agreement_vs_block0": check_block,
            "note": (
                "Independent Euclidean CA RMSD after the RCSB 4x4 is compared "
                "to the global summary RMSD. Block RMSD is a separate API "
                "value and is not fused with the global score."
            ),
            "original_target_untouched": True,
        }
    except protein_structure.StructureError as exc:
        superposition_proof = {"status": "UNAVAILABLE", "reason": str(exc)}
    engine_validation.record_remote_validation(
        "RCSB Alignment API",
        ok=True,
        version="api/v1",
        details={
            "pair": "8HSK vs 8HSF",
            "method": record["method"],
            "n_pairs": record["n_aligned_residue_pairs"],
            "rmsd_block0": record["rmsd_block0_angstrom"],
            "rmsd_global": record["rmsd_global_angstrom"],
            "ticket": record.get("ticket"),
        },
    )
    evidence_dir = tmp_path / "phase17_evidence"
    evidence_dir.mkdir()
    (evidence_dir / "live_rcsb_alignment.json").write_text(
        json.dumps(
            {
                "engine_location": "remote",
                "status": "REMOTE_VALIDATED",
                "method": record["method"],
                "n_aligned_residue_pairs": record["n_aligned_residue_pairs"],
                "rmsd_block0_angstrom": record["rmsd_block0_angstrom"],
                "rmsd_global_angstrom": record["rmsd_global_angstrom"],
                "tm_scores": record["tm_scores"],
                "aln_coverage_percent": record["aln_coverage_percent"],
                "mapping_A2": mapped,
                "ticket": record.get("ticket"),
                "identity_hash": record.get("identity_hash"),
                "superposition_proof": superposition_proof,
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )


def test_tool_registry_lists_rcsb_as_remote():
    from modules import tool_registry

    snap = tool_registry.collect_tool_snapshot(refresh=True)
    row = snap["tools"]["RCSB Alignment API"]
    assert row["engine_location"] == "remote"
    assert "tm-align" in row["methods"]
    assert row["live"] is None or row["live"].get("status") in {
        engine_validation.STATUS_REMOTE_VALIDATED,
        engine_validation.STATUS_BROKEN,
        None,
    }
    us_row = snap["tools"]["US-align"]
    assert us_row["name"] == "US-align"
