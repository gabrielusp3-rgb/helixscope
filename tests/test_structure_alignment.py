"""Parser do JSON oficial RCSB Alignment API. Golden = fixture da documentacao."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from modules import structure_alignment

FIXTURES = Path(__file__).parent / "fixtures"


def _official_complete() -> dict:
    return json.loads((FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8"))


def test_official_8hsk_8hsf_scores_are_not_fused():
    record = structure_alignment.parse_alignment_payload(_official_complete())
    assert record["method"] == "fatcat-rigid"
    assert record["method_display"] == "jFATCAT-rigid"
    assert record["method_kind"] == "rigid"
    assert record["atoms_fitted"] == "CA"
    assert record["n_aligned_residue_pairs"] == 20
    assert record["n_blocks"] == 1
    assert record["rmsd_block0_angstrom"] == pytest.approx(0.8)
    assert record["rmsd_global_angstrom"] == pytest.approx(0.99)
    assert record["rmsd_block0_angstrom"] != record["rmsd_global_angstrom"]
    tm = record["tm_scores"]
    assert len(tm) == 1
    assert tm[0]["type"] == "TM-score"
    assert tm[0]["value"] == pytest.approx(0.23)
    assert record["sequence_identity"] == pytest.approx(1)
    assert record["aln_coverage_percent"] == [95, 100]
    assert record["reference"]["entry_id"] == "8HSK"
    assert record["reference"]["kind_label"] == "EXPERIMENTAL"
    assert record["target"]["entry_id"] == "8HSF"
    assert record["superposition_visual"] == "AVAILABLE"
    assert record["mapping_reference_to_target"]["A:2"] == "A:2"
    unmapped = structure_alignment.lookup_correspondence(
        record, side="reference", asym_id="A", label_seq_id=1
    )
    assert unmapped["status"] == "UNMAPPED"
    mapped = structure_alignment.lookup_correspondence(
        record, side="reference", asym_id="A", label_seq_id=2
    )
    assert mapped["status"] == "AVAILABLE"
    assert mapped["corresponding"] == "A:2"
    reverse = structure_alignment.lookup_correspondence(
        record, side="target", asym_id="A", label_seq_id=2
    )
    assert reverse["corresponding"] == "A:2"


def test_predicted_kind_label_for_alphafold_csm():
    payload = _official_complete()
    payload["results"][0]["structures"][1]["entry_id"] = "AF_AFP01308F1"
    record = structure_alignment.parse_alignment_payload(payload)
    assert record["target"]["kind_label"] == "PREDICTED"
    assert record["reference"]["kind_label"] == "EXPERIMENTAL"


def test_flexible_multi_block_is_partial_visual_not_one_transform():
    payload = json.loads(
        (FIXTURES / "rcsb_alignment_flexible.json").read_text(encoding="utf-8")
    )
    record = structure_alignment.parse_alignment_payload(payload)
    assert record["method_kind"] == "flexible"
    assert record["n_blocks"] == 2
    assert record["superposition_visual"] == "PARTIAL"
    assert record["blocks"][0]["rmsd_angstrom"] == pytest.approx(1.2)
    assert record["blocks"][1]["rmsd_angstrom"] == pytest.approx(2.4)
    assert record["rmsd_global_angstrom"] == pytest.approx(2.1)
    assert record["n_aligned_residue_pairs"] == 7
    assert record["mapping_reference_to_target"]["A:10"] == "A:20"
    assert record["mapping_reference_to_target"]["A:80"] == "A:90"


def test_incomplete_payload_is_parsing_error():
    with pytest.raises(structure_alignment.StructureAlignmentError) as exc:
        structure_alignment.parse_alignment_payload({"info": {"status": "RUNNING"}})
    assert exc.value.category in {"INVALID_INPUT", "PARSING_ERROR"}


def test_summary_rows_show_na_not_zero_when_tm_missing():
    payload = _official_complete()
    payload["results"][0]["summary"]["scores"] = [
        {"value": 0.5, "type": "RMSD"},
    ]
    record = structure_alignment.parse_alignment_payload(payload)
    assert record["tm_scores"] == []
    rows = dict(structure_alignment.comparison_summary_rows(record))
    assert rows["TM-score(s)"] == "N/A"
    assert "0" != rows["TM-score(s)"]
