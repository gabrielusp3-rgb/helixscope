"""Evidence workspace: proveniencia, conflitos sem voto, sem score inventado."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from modules import evidence_workspace, structure_alignment

FIXTURES = Path(__file__).parent / "fixtures"


def test_evidence_item_rejects_unknown_status():
    with pytest.raises(evidence_workspace.EvidenceError):
        evidence_workspace.evidence_item(
            field="x",
            value=1,
            source="test",
            evidence_status="PRETTY_SURE",
        )


def test_alignment_evidence_and_export_have_no_confidence_percent():
    payload = json.loads((FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8"))
    record = structure_alignment.parse_alignment_payload(payload)
    items = evidence_workspace.collect_from_alignment(record)
    fields = {item["field"] for item in items}
    assert "rmsd_global_angstrom" in fields
    assert "n_aligned_residue_pairs" in fields
    pack = evidence_workspace.export_pack(kind="structures", inputs={"a": "8HSK"}, items=items)
    assert pack["confidence_score"] is None
    assert "percentage" in pack["confidence_note"].lower() or "N/A" in pack["confidence_note"]


def test_source_conflicts_are_listed_not_voted():
    items = [
        evidence_workspace.evidence_item(
            field="Variant A most_severe_consequence",
            value="missense_variant",
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
        ),
        evidence_workspace.evidence_item(
            field="Variant B most_severe_consequence",
            value="synonymous_variant",
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
        ),
    ]
    conflicts = evidence_workspace.source_conflicts(items)
    assert len(conflicts) == 1
    assert conflicts[0]["field"] == "most_severe_consequence"
    assert "majority" in conflicts[0]["policy"].lower()
