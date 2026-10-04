"""Superposicao: transformacao 4x4, RMSD independente, coordenadas originais."""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path

import pytest

from modules import structure_alignment, structure_superposition
from ui import structure_viewer

FIXTURES = Path(__file__).parent / "fixtures"


def _atom(chain: str, seq: int, x: float, y: float, z: float) -> dict:
    return {
        "group": "ATOM",
        "atom_name": "CA",
        "label_asym_id": chain,
        "auth_asym_id": chain,
        "label_seq_id": seq,
        "auth_seq_id": seq,
        "comp_id": "ALA",
        "x": x,
        "y": y,
        "z": z,
    }


def test_column_major_translation_matches_hand_computation():
    matrix = [
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        -1.0, -2.0, -3.0, 1.0,
    ]
    xp, yp, zp = structure_superposition.apply_column_major_4x4(1.0, 2.0, 3.0, matrix)
    assert (xp, yp, zp) == pytest.approx((0.0, 0.0, 0.0))
    rotated = [
        0.0, 1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ]
    xr, yr, zr = structure_superposition.apply_column_major_4x4(1.0, 0.0, 0.0, rotated)
    assert (xr, yr, zr) == pytest.approx((0.0, 1.0, 0.0))


def test_independent_rmsd_zero_after_known_translation_not_api_echo():
    reference = {"atoms": [_atom("A", 1, 0.0, 0.0, 0.0), _atom("A", 2, 3.0, 0.0, 0.0)]}
    target = {"atoms": [_atom("A", 1, 1.0, 2.0, 3.0), _atom("A", 2, 4.0, 2.0, 3.0)]}
    original_target = copy.deepcopy(target)
    matrix = [
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        -1.0, -2.0, -3.0, 1.0,
    ]
    pairs = [
        {
            "reference_asym_id": "A",
            "target_asym_id": "A",
            "reference_label_seq_id": 1,
            "target_label_seq_id": 1,
        },
        {
            "reference_asym_id": "A",
            "target_asym_id": "A",
            "reference_label_seq_id": 2,
            "target_label_seq_id": 2,
        },
    ]
    untransformed = math.sqrt(((1.0) ** 2 + (2.0) ** 2 + (3.0) ** 2))
    assert untransformed == pytest.approx(math.sqrt(14.0))
    result = structure_superposition.independent_aligned_ca_rmsd(
        reference_parsed=reference,
        target_parsed=target,
        pairs=pairs,
        target_matrix=matrix,
    )
    assert result["status"] == "COMPUTED"
    assert result["n_pairs_used"] == 2
    assert result["rmsd_angstrom"] == pytest.approx(0.0, abs=1e-9)
    assert target == original_target
    copies = structure_superposition.transform_atom_copies(target["atoms"], matrix)
    assert copies[0]["x"] == pytest.approx(0.0)
    assert copies[0]["original_x"] == pytest.approx(1.0)
    assert target["atoms"][0]["x"] == pytest.approx(1.0)


def test_missing_ca_is_unmapped_not_zero_rmsd():
    reference = {"atoms": [_atom("A", 1, 0.0, 0.0, 0.0)]}
    target = {"atoms": []}
    identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    result = structure_superposition.independent_aligned_ca_rmsd(
        reference_parsed=reference,
        target_parsed=target,
        pairs=[
            {
                "reference_asym_id": "A",
                "target_asym_id": "A",
                "reference_label_seq_id": 1,
                "target_label_seq_id": 1,
            }
        ],
        target_matrix=identity,
    )
    assert result["status"] == "UNMAPPED"
    assert result["rmsd_angstrom"] is None


def test_superposition_bundle_does_not_mutate_inputs_and_flexible_is_partial():
    payload = json.loads(
        (FIXTURES / "rcsb_alignment_flexible.json").read_text(encoding="utf-8")
    )
    alignment = structure_alignment.parse_alignment_payload(payload)
    reference = {
        "atoms": [
            _atom("A", 10, 0, 0, 0),
            _atom("A", 11, 1, 0, 0),
            _atom("A", 12, 2, 0, 0),
            _atom("A", 13, 3, 0, 0),
        ]
    }
    target = {
        "atoms": [
            _atom("A", 20, 5, 0, 0),
            _atom("A", 21, 6, 0, 0),
            _atom("A", 22, 7, 0, 0),
            _atom("A", 23, 8, 0, 0),
        ]
    }
    snapshot = copy.deepcopy(target)
    bundle = structure_superposition.superposition_bundle(
        alignment=alignment,
        reference_parsed=reference,
        target_parsed=target,
        block_index=0,
    )
    assert bundle["visual_status"] == "PARTIAL"
    assert bundle["original_target_untouched"] is True
    assert target == snapshot
    assert any(atom.get("transformed") for atom in bundle["transformed_target_atoms"])
    figure = structure_viewer.figure_from_superposition(bundle)
    assert len(figure.data) == 2
    check = structure_superposition.rmsd_agrees(
        bundle["independent_rmsd"],
        alignment["blocks"][0]["rmsd_angstrom"],
        abs_tolerance_angstrom=0.15,
        rel_tolerance=0.1,
    )
    assert "abs_tolerance_angstrom" in check
