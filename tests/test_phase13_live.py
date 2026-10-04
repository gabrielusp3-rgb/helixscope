"""Testes live opcionais. Skip se o binario oficial nao estiver instalado.

MOCK VALIDATION != LIVE VALIDATION. Estes testes chamam o executavel real.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from modules import msa, phylogeny, protein_structure, provenance

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN_CIF = FIXTURES / "1CRN.cif"


def _tiny_dna_msa() -> dict:
    return msa.import_prealigned_fasta(
        ">seq_a\nACGTACGTACGTACGTACGT\n"
        ">seq_b\nACGTACGTACGTACGTTTTT\n"
        ">seq_c\nTTTTACGTACGTACGTACGT\n"
        ">seq_d\nACGTACGTAAAAACGTACGT\n"
    )


def test_live_iqtree_version_and_ml_tree():
    info = phylogeny.detect_iqtree()
    if not info.get("available"):
        pytest.skip("IQ-TREE official binary not installed")
    assert str(info.get("version") or "").strip()
    result = phylogeny.infer_phylogeny(
        _tiny_dna_msa(),
        method=phylogeny.METHOD_IQTREE,
        distance_model=phylogeny.DISTANCE_JC69,
        rooting=phylogeny.ROOTING_UNROOTED,
        iqtree_model_mode=phylogeny.IQTREE_MODEL_USER,
    )
    assert result["status"] == "COMPUTED"
    assert result["n_leaves"] == 4
    assert result["tool"] == "IQ-TREE"
    assert str(result.get("tool_version") or "").strip()
    assert result["newick"].count("seq_") == 4
    assert result["branch_lengths_present"] is True


def test_live_iqtree_modelfinder_ufboot_and_alrt():
    info = phylogeny.detect_iqtree()
    if not info.get("available"):
        pytest.skip("IQ-TREE official binary not installed")
    result = phylogeny.infer_phylogeny(
        _tiny_dna_msa(),
        method=phylogeny.METHOD_IQTREE,
        distance_model=phylogeny.DISTANCE_JC69,
        rooting=phylogeny.ROOTING_UNROOTED,
        iqtree_model_mode=phylogeny.IQTREE_MODEL_MFP,
        ufboot_replicates=1000,
        alrt_replicates=1000,
    )
    assert result["status"] == "COMPUTED"
    assert result["iqtree_model_mode"] == phylogeny.IQTREE_MODEL_MFP
    methods = [str(item.get("method")) for item in list((result.get("support") or {}).get("methods") or [])]
    assert phylogeny.SUPPORT_UFBOOT in methods
    assert phylogeny.SUPPORT_SH_ALRT in methods
    assert phylogeny.SUPPORT_FELSENSTEIN not in methods
    support = result.get("support") or {}
    assert "ULTRAFAST" in str(support.get("method") or support.get("label") or "").upper() or phylogeny.SUPPORT_UFBOOT in methods


def test_live_iqtree_timeout():
    info = phylogeny.detect_iqtree()
    if not info.get("available"):
        pytest.skip("IQ-TREE official binary not installed")
    original = phylogeny.ML_TIMEOUT_S
    phylogeny.ML_TIMEOUT_S = 0.01
    try:
        with pytest.raises(phylogeny.PhylogenyError) as exc:
            phylogeny.infer_phylogeny(
                _tiny_dna_msa(),
                method=phylogeny.METHOD_IQTREE,
                distance_model=phylogeny.DISTANCE_JC69,
                rooting=phylogeny.ROOTING_UNROOTED,
            )
        assert exc.value.category == "TIMEOUT"
    finally:
        phylogeny.ML_TIMEOUT_S = original


def test_live_fasttree_approximate_ml():
    info = phylogeny.detect_fasttree()
    if not info.get("available"):
        pytest.skip("FastTree official binary not installed")
    result = phylogeny.infer_phylogeny(
        _tiny_dna_msa(),
        method=phylogeny.METHOD_FASTTREE,
        rooting=phylogeny.ROOTING_UNROOTED,
    )
    assert result["status"] == "COMPUTED"
    assert "approximate" in str(result.get("method_label") or result.get("algorithm") or "").lower() or result["tool"] == "FastTree"
    assert result["n_leaves"] == 4


def test_live_dssp_on_1crn():
    text = CRAMBIN_CIF.read_text(encoding="utf-8")
    local = protein_structure.dssp_availability()
    if local.get("available"):
        assignment = protein_structure.assign_secondary_structure_dssp(text)
        assert assignment.get("engine_location") != "remote"
    else:
        assignment = protein_structure.assign_secondary_structure_dssp_remote(text)
        if str(assignment.get("status") or "") not in {"EXPERIMENTAL", "COMPUTED"}:
            pytest.skip(
                "local mkdssp is not installed and the PDB-REDO DSSP API did not "
                f"return an assignment ({assignment.get('category')}: "
                f"{assignment.get('reason')})"
            )
        assert assignment.get("engine_location") == "remote"
        assert assignment.get("not_local_mkdssp") is True
        assert assignment.get("not_sequence_prediction") is True
    assert assignment["status"] in {"EXPERIMENTAL", "COMPUTED"}
    assert assignment.get("available") is True
    dssp_residues = list(assignment.get("residues") or [])
    assert dssp_residues
    parsed = protein_structure.parse_mmcif(text)
    structure_residues = []
    for chain in list(parsed.get("chains") or []):
        chain_id = str(chain.get("chain_id") or "A")
        for residue in list(chain.get("residues") or []):
            row = dict(residue)
            row.setdefault("chain_id", chain_id)
            structure_residues.append(row)
    if not structure_residues:
        pytest.skip("1CRN fixture did not expose residue records for mapping")
    mapped = protein_structure.map_dssp_assignments(dssp_residues, structure_residues)
    assert int(mapped.get("n_mapped") or 0) >= 1


def test_live_edtsurf_on_1crn():
    info = protein_structure.surface_availability()
    if not info.get("available") or "edtsurf" not in str(info.get("tool") or "").lower():
        pytest.skip("official EDTSurf Windows executable not installed")
    text = CRAMBIN_CIF.read_text(encoding="utf-8")
    parsed = protein_structure.parse_mmcif(text)
    mesh = protein_structure.generate_molecular_surface(
        kind="experimental",
        structure_hash=str(parsed.get("content_hash") or "1CRN"),
        parsed=parsed,
        surface_code=3,
        probe_radius=1.4,
    )
    assert mesh["status"] == "EXPERIMENTAL"
    assert mesh.get("backend") == "EDTSurf"
    assert mesh.get("surface_definition") == "MS"
    assert mesh.get("not_shrake_rupley") is True
    assert int(mesh.get("n_vertices") or 0) >= 3
    assert int(mesh.get("n_faces") or 0) >= 1
    protein_structure.validate_surface_mesh(list(mesh.get("vertices") or []))


def test_live_msms_on_1crn():
    info = protein_structure.msms_availability()
    if not info.get("available"):
        pytest.skip("official MSMS Windows executable not installed")
    text = CRAMBIN_CIF.read_text(encoding="utf-8")
    parsed = protein_structure.parse_mmcif(text)
    mesh = protein_structure._run_msms_surface(
        list(parsed.get("atoms") or []),
        kind="experimental",
        structure_hash=str(parsed.get("content_hash") or "1CRN"),
        probe_radius=1.4,
    )
    assert mesh["status"] == "EXPERIMENTAL"
    assert mesh.get("backend") == "MSMS"
    assert mesh.get("surface_definition") == "SES"
    assert mesh.get("not_shrake_rupley") is True
    assert int(mesh.get("n_vertices") or 0) >= 3
    assert int(mesh.get("n_faces") or 0) >= 1
    protein_structure.validate_surface_mesh(list(mesh.get("vertices") or []))
