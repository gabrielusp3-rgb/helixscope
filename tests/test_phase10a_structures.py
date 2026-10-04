"""Fase 10A: DNA/RNA experimentais reais, mapping, DSSP/STRIDE/surface, complexo Cas.

1BNA.cif e 1RNA.cif sao copias publicas RCSB PDB (nao geometria inventada).
parser_nucleic_mini.cif continua a NAO ser uma deposicao. 4UN3 mmCIF completo
excede o antigo teto de 2 MB; identidade de cadeias vem da Data API
polymer_entity. Nenhum teste monta Cas9+guide+DNA procedural.
"""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path

import pytest

from modules import (
    complex_structure,
    dna_3d,
    molecule_selection,
    na_structure_catalog,
    protein_structure,
    provenance,
    region_nav,
    rna_3d,
    scientific_checks,
    structure_scene,
)
from ui import structure_viewer
from ui.tokens import STRUCTURE_VIEWPORT_HEIGHT

FIXTURES = Path(__file__).parent / "fixtures"
DICKERSON = "CGCGAATTCGCG"
RNA_DUPLEX = "UUAUAUAUAUAUAA"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _fail_network(*_args, **_kwargs):
    raise AssertionError("network")


def _load_json(name: str) -> dict:
    return json.loads(_load(name))


def _1bna_result(*, sequence: str = DICKERSON, chain_id: str = "A") -> dict:
    meta = _load_json("rcsb_entry_1BNA.json")
    return protein_structure.load_protein_structure(
        sequence,
        source="RCSB PDB",
        structure_id="1BNA",
        structure_text=_load("1BNA.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
        molecule="DNA",
        chain_id=chain_id,
    )


def _1rna_result(*, sequence: str = RNA_DUPLEX) -> dict:
    meta = _load_json("rcsb_entry_1RNA.json")
    return protein_structure.load_protein_structure(
        sequence,
        source="RCSB PDB",
        structure_id="1RNA",
        structure_text=_load("1RNA.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
        molecule="RNA",
    )


def _1crn_result() -> dict:
    meta = protein_structure.parse_rcsb_entry_metadata(_load_json("rcsb_entry_1CRN.json"))
    return protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
    )


def _4un3_entities() -> list[dict]:
    return [
        _load_json("rcsb_polymer_entity_4UN3_1.json"),
        _load_json("rcsb_polymer_entity_4UN3_2.json"),
        _load_json("rcsb_polymer_entity_4UN3_3.json"),
        _load_json("rcsb_polymer_entity_4UN3_4.json"),
    ]


def test_twelve_evidence_statuses_are_reachable():
    reached: dict[str, str] = {}
    retrieved = _1bna_result()
    assert retrieved["status"] == "RETRIEVED"
    reached["RETRIEVED"] = "load_protein_structure 1BNA"
    from modules import dna_analysis

    _ = dna_analysis.gc_content(DICKERSON)
    computed = provenance.analysis_envelope(
        module="TEST",
        payload={"gc": True},
        status="COMPUTED",
        algorithm="gc_content",
    )
    reached[computed["status"]] = "analysis_envelope COMPUTED"
    pred = provenance.analysis_envelope(
        module="TEST",
        payload={"kind": "predicted"},
        status="PREDICTED",
        algorithm="AlphaFold DB metadata",
        source="AlphaFold DB",
    )
    reached[pred["status"]] = "AlphaFold metadata envelope"
    heur = provenance.analysis_envelope(
        module="TEST",
        payload={},
        status="HEURISTIC",
        algorithm="test heuristic",
    )
    reached[heur["status"]] = "analysis_envelope HEURISTIC"
    exp = dna_3d.representation_for_sequence(
        DICKERSON, deposited=protein_structure.structure_3d_input(retrieved)
    )
    assert exp["status"] == "EXPERIMENTAL"
    reached["EXPERIMENTAL"] = "dna_3d 1BNA"
    ill = dna_3d.representation_for_sequence(DICKERSON, prefer="illustrative")
    assert ill["status"] == "ILLUSTRATIVE"
    reached["ILLUSTRATIVE"] = "dna_3d default helix"
    dssp = protein_structure.assign_secondary_structure_dssp("")
    assert dssp["status"] == "UNAVAILABLE"
    reached["UNAVAILABLE"] = "DSSP absent"
    bad = dna_3d.representation_for_sequence("ZZZZ")
    assert bad["status"] == "ERROR"
    reached["ERROR"] = "invalid DNA"
    partial = protein_structure.map_query_to_chain(
        DICKERSON + "AAAA", retrieved["selected_chain"], molecule="DNA"
    )
    assert partial["mapping_status"] == "PARTIAL"
    reached["PARTIAL"] = "1BNA plus unmatched 3' bases"
    unmapped = protein_structure.map_query_to_chain(
        "WWWWWWWWWW", _1crn_result()["selected_chain"], molecule="PROTEIN"
    )
    assert unmapped["mapping_status"] == "UNMAPPED"
    reached["UNMAPPED"] = "poly-W vs crambin"
    envelope = protein_structure.structure_3d_input(retrieved)
    stale = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=provenance.sequence_digest("AT"),
        representation="backbone",
    )
    assert stale["status"] == "STALE"
    reached["STALE"] = "structure_scene hash mismatch"
    previous = protein_structure.MAX_STRUCTURE_BYTES
    try:
        protein_structure.MAX_STRUCTURE_BYTES = 20
        with pytest.raises(protein_structure.StructureError) as cap:
            protein_structure.parse_mmcif(_load("1BNA.cif"))
        assert cap.value.category == "RESOURCE_LIMIT"
    finally:
        protein_structure.MAX_STRUCTURE_BYTES = previous
    reached["RESOURCE_LIMIT"] = "parse_mmcif oversized cap"
    missing = set(provenance.EVIDENCE_STATUSES) - set(reached)
    assert not missing, missing


def test_1bna_experimental_flow_matches_deposition():
    parsed = protein_structure.parse_mmcif(_load("1BNA.cif"))
    assert parsed["entry_id"] == "1BNA"
    assert parsed["method"] == "X-RAY DIFFRACTION"
    assert parsed["resolution_angstrom"] == pytest.approx(1.9)
    assert parsed["n_atoms"] == 566
    chain_a = next(item for item in parsed["chains"] if item["chain_id"] == "A")
    assert chain_a["sequence"] == DICKERSON
    assert chain_a["polymer_type"] == "polydeoxyribonucleotide"
    protein_structure.check_structure_integrity(
        parsed, molecule="DNA", chain_id="A", expected_sequence=DICKERSON
    )
    c1 = next(
        atom
        for atom in parsed["atoms"]
        if atom.get("auth_asym_id") == "A"
        and atom.get("label_seq_id") == 1
        and atom.get("atom_name") == "C1'"
    )
    mapped = protein_structure.map_query_to_chain(DICKERSON, chain_a, molecule="DNA")
    assert mapped["mapping_status"] == "EXACT"
    row0 = mapped["mapping"][0]
    assert row0["label_seq_id"] == 1
    assert row0["c1"]["x"] == pytest.approx(c1["x"])
    assert row0["c1"]["y"] == pytest.approx(c1["y"])
    assert row0["c1"]["z"] == pytest.approx(c1["z"])
    result = _1bna_result()
    prov = protein_structure.experimental_structure_provenance(result)
    assert prov["pdb_id"] == "1BNA"
    assert prov["source"] == "RCSB PDB"
    assert prov["chain"] == "A"
    assert prov["method"] == "X-RAY DIFFRACTION"
    assert prov["resolution_angstrom"] == pytest.approx(1.9)
    assert prov["sequence_hash"] == provenance.sequence_digest(DICKERSON)
    assert prov["structure_hash"] == result["content_hash"]
    assert prov["timestamp"]
    envelope = protein_structure.structure_3d_input(result)
    auto = dna_3d.representation_for_sequence(DICKERSON, deposited=envelope)
    ill = dna_3d.representation_for_sequence(DICKERSON, deposited=envelope, prefer="illustrative")
    assert auto["kind"] == dna_3d.KIND_EXPERIMENTAL
    assert ill["kind"] == dna_3d.KIND_ILLUSTRATIVE


def test_1bna_chain_b_auth_and_label_diverge():
    result = _1bna_result(chain_id="B")
    row = result["residue_mapping"][0]
    assert row["label_seq_id"] == 1
    assert row["auth_seq_id"] == 13
    assert row["query_index_0based"] == 0
    last = result["residue_mapping"][-1]
    assert last["label_seq_id"] == 12
    assert last["auth_seq_id"] == 24


def test_1bna_partial_and_unmapped_and_validation_fields():
    parsed = protein_structure.parse_mmcif(_load("1BNA.cif"))
    chain_a = next(item for item in parsed["chains"] if item["chain_id"] == "A")
    partial = protein_structure.map_query_to_chain(DICKERSON + "AAAA", chain_a, molecule="DNA")
    assert partial["mapping_status"] == "PARTIAL"
    assert partial["mapped_query_indices"]
    assert partial["unmapped_query_indices"]
    assert max(partial["unmapped_query_indices"]) >= 12
    unmapped = protein_structure.map_query_to_chain(
        "WWWWWWWWWW", _1crn_result()["selected_chain"], molecule="PROTEIN"
    )
    assert unmapped["mapping_status"] == "UNMAPPED"
    assert unmapped["n_letter_matches"] == 0
    with pytest.raises(protein_structure.StructureError) as missing:
        protein_structure.check_structure_integrity(parsed, molecule="DNA", chain_id="Z")
    assert missing.value.category == "INVALID_INPUT"
    with pytest.raises(protein_structure.StructureError) as polymer:
        protein_structure.check_structure_integrity(parsed, molecule="PROTEIN", chain_id="A")
    assert polymer.value.category == "INVALID_INPUT"
    broken = copy.deepcopy(parsed)
    broken["chains"][0]["residues"][0]["label_seq_id"] = 99
    broken["chains"][0]["residues"][1]["label_seq_id"] = 1
    with pytest.raises(protein_structure.StructureError) as order:
        protein_structure.check_structure_integrity(broken, molecule="DNA", chain_id="A")
    assert order.value.category == "PARSING_ERROR"
    length_broken = copy.deepcopy(parsed)
    length_broken["chains"][0]["length"] = 3
    with pytest.raises(protein_structure.StructureError):
        protein_structure.check_structure_integrity(length_broken, molecule="DNA", chain_id="A")
    seq_broken = copy.deepcopy(parsed)
    seq_broken["chains"][0]["sequence"] = "CGCGAATTCGCN"
    with pytest.raises(protein_structure.StructureError):
        protein_structure.check_structure_integrity(seq_broken, molecule="DNA", chain_id="A")
    nan_atoms = copy.deepcopy(parsed)
    nan_atoms["atoms"][0]["x"] = float("nan")
    with pytest.raises(protein_structure.StructureError) as nan_exc:
        protein_structure.check_structure_integrity(nan_atoms, molecule="DNA", chain_id="A")
    assert nan_exc.value.category == "PARSING_ERROR"
    mismatch = dna_3d.representation_for_sequence(
        "ATGCATGCATGC", deposited=protein_structure.structure_3d_input(_1bna_result())
    )
    assert mismatch["status"] == "UNMAPPED"


def test_1bna_motif_and_selection_use_deposited_coords():
    result = _1bna_result()
    envelope = protein_structure.structure_3d_input(result)
    motif = dna_3d.highlight_from_motif(DICKERSON, "GAATTC")
    assert motif["n_hits"] >= 1
    assert any(DICKERSON[span["start"] : span["end"]] == "GAATTC" for span in motif["spans"])
    assert motif["kind"] == "sequence_span"
    rep = dna_3d.representation_for_sequence(DICKERSON, deposited=envelope)
    mapped = dna_3d.position_to_3d(rep, 4)
    assert mapped["status"] == "READY"
    assert mapped["mapping"]["label_seq_id"] == 5
    assert mapped["mapping"]["has_coordinates"] is True
    scene = dna_3d.scene_from_representation(
        rep, selected_query_index=4, highlight_query_indices=motif["indices"]
    )
    assert scene["status"] == "READY"
    assert 4 in scene["selection_ids"]


def test_1rna_experimental_and_pairing_geometry():
    result = _1rna_result()
    assert result["kind"] == "experimental"
    assert result["structure_id"] == "1RNA"
    envelope = protein_structure.structure_3d_input(result)
    rep = rna_3d.representation_for_sequence(RNA_DUPLEX, deposited=envelope)
    assert rep["kind"] == rna_3d.KIND_EXPERIMENTAL
    n = len(RNA_DUPLEX)
    sequential = rna_3d.pairing_geometry_consistency(
        rep, [{"position_0based": 0, "paired_position_0based": 1}]
    )
    assert sequential["n_within"] >= 1
    parsed = protein_structure.parse_mmcif(_load("1RNA.cif"))
    chain_a = next(item for item in parsed["chains"] if item["chain_id"] == "A")
    chain_b = next(item for item in parsed["chains"] if item["chain_id"] == "B")
    pa = chain_a["residues"][0]["c1"]
    pb = chain_b["residues"][n - 1]["c1"]
    dx = float(pa["x"]) - float(pb["x"])
    dy = float(pa["y"]) - float(pb["y"])
    dz = float(pa["z"]) - float(pb["z"])
    terminal_pair = (dx * dx + dy * dy + dz * dz) ** 0.5
    assert terminal_pair < 25.0
    mfe_as_3d = rna_3d.pairing_sync(
        rep,
        {
            "status": "PREDICTED",
            "base_pairs": [{"position_0based": 0, "paired_position_0based": 1}],
            "sequence_hash": result["sequence_hash"],
        },
    )
    assert mfe_as_3d["geometry_kind"] == rna_3d.KIND_EXPERIMENTAL
    rnac = na_structure_catalog.rnacentral_is_appropriate_for_3d()
    assert rnac["appropriate_for_3d_coordinates"] is False
    predictors = na_structure_catalog.rna_3d_predictor_research()
    assert len(predictors) >= 3
    assert all(item["decision"] == "not integrated" for item in predictors)


def test_1crn_deposited_ss_is_not_dssp_prediction():
    parsed = protein_structure.parse_mmcif(_load("1CRN.cif"))
    ss = parsed["deposited_secondary_structure"]
    assert ss["not_dssp"] is True
    assert ss["not_sequence_prediction"] is True
    helix_spans = {(item["beg_label_seq_id"], item["end_label_seq_id"]) for item in ss["helices"]}
    assert (7, 19) in helix_spans
    assert (23, 30) in helix_spans
    sheet_spans = {(item["beg_label_seq_id"], item["end_label_seq_id"]) for item in ss["sheets"]}
    assert (1, 4) in sheet_spans
    assert (32, 35) in sheet_spans
    dssp_absent = protein_structure.assign_secondary_structure_dssp(_load("1CRN.cif"))
    assert dssp_absent["status"] == "UNAVAILABLE"
    malformed = protein_structure.assign_secondary_structure_dssp(output_text="returncode zero but empty")
    assert malformed["status"] == "ERROR"
    valid = protein_structure.assign_secondary_structure_dssp(output_text=_load("dssp_legacy_min.txt"))
    assert valid["status"] == "EXPERIMENTAL"
    assert valid["n_helix"] == 1
    assert valid["n_sheet"] == 1
    assert valid["n_coil"] == 1
    stride = protein_structure.stride_availability()
    assert stride["available"] is False
    source = Path(protein_structure.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    assert "assign_stride" not in names
    assert "house_dssp" not in names


def test_surface_refuses_illustrative_and_validates_vertices():
    with pytest.raises(protein_structure.StructureError):
        protein_structure.generate_molecular_surface(kind="illustrative")
    missing = protein_structure.generate_molecular_surface(kind="experimental", structure_hash="abc")
    assert missing["status"] == "UNAVAILABLE"
    ok = protein_structure.validate_surface_mesh(
        [
            {"x": 1.0, "y": 2.0, "z": 3.0},
            {"x": 4.0, "y": 5.0, "z": 6.0},
            {"x": 7.0, "y": 8.0, "z": 9.0},
        ]
    )
    assert ok["finite"] is True
    with pytest.raises(protein_structure.StructureError):
        protein_structure.validate_surface_mesh([{"x": float("nan"), "y": 0.0, "z": 0.0}])
    research = protein_structure.surface_tool_research()
    assert {item["tool"] for item in research} >= {"MSMS", "EDTSurf", "Mol*"}


def test_msms_xyzr_and_parsers_skip_headers():
    atoms = [
        {"x": 0.0, "y": 0.0, "z": 0.0, "element": "C"},
        {"x": 1.5, "y": 0.0, "z": 0.0, "element": "N"},
        {"x": 0.0, "y": 1.5, "z": 0.0, "element": "O"},
    ]
    xyzr = protein_structure.atoms_to_xyzr(atoms)
    rows = [line for line in xyzr.splitlines() if line.strip()]
    assert len(rows) == 3
    assert "1.700" in xyzr or "1.55" in xyzr or "1.52" in xyzr
    vert = (
        "# MSMS 2.6.1\n"
        "     3      1      1.000      1.400\n"
        "  0.0 0.0 0.0  1.0 0.0 0.0  1  1\n"
        "  1.0 0.0 0.0  0.0 1.0 0.0  1  1\n"
        "  0.0 1.0 0.0  0.0 0.0 1.0  1  1\n"
    )
    vertices = protein_structure.parse_msms_vert(vert)
    assert len(vertices) == 3
    assert vertices[0]["x"] == pytest.approx(0.0)
    assert vertices[0]["nx"] == pytest.approx(1.0)
    face = (
        "# faces\n"
        "     1      1      1.000      1.400\n"
        "     1     2     3      0\n"
    )
    faces = protein_structure.parse_msms_face(face, n_vertices=3)
    assert faces == [[0, 1, 2]]
    edt = protein_structure.edtsurf_availability()
    msms = protein_structure.msms_availability()
    assert "available" in edt and "available" in msms
    if edt.get("available") and msms.get("available"):
        surface = protein_structure.surface_availability()
        assert "edtsurf" in str(surface.get("tool") or "").lower()


def test_4un3_polymer_entities_identify_cas_guide_target_nontarget():
    meta = {
        "method": "X-RAY DIFFRACTION",
        "resolution_angstrom": 2.593,
        "kind": "experimental",
        "structure_id": "4UN3",
    }
    shuffled = list(reversed(_4un3_entities()))
    annotated = complex_structure.annotate_from_polymer_entities(
        shuffled,
        structure_id="4UN3",
        metadata=meta,
        guide_sequence="GGAUAACUCAAUUUGUAAAA",
        target_sequence="CAATACCATTTTTTACAAATTGAGTTAT",
        pam="NGG",
    )
    assert annotated["synthetic_complex"] is False
    roles = {item["role"]: item for item in annotated["chains"]}
    assert complex_structure.ROLE_CAS_PROTEIN in roles
    assert complex_structure.ROLE_GUIDE_RNA in roles
    assert complex_structure.ROLE_TARGET_DNA in roles
    assert complex_structure.ROLE_NONTARGET_DNA in roles
    assert roles[complex_structure.ROLE_CAS_PROTEIN]["chain_id"] == "B"
    assert roles[complex_structure.ROLE_GUIDE_RNA]["chain_id"] == "A"
    assert roles[complex_structure.ROLE_TARGET_DNA]["chain_id"] == "C"
    assert roles[complex_structure.ROLE_NONTARGET_DNA]["chain_id"] == "D"
    assert annotated["guide_mapping"]["status"] == "EXACT"
    assert annotated["target_mapping"]["status"] == "EXACT"
    assert annotated["pam_mapping"]["status"] == "EXACT"
    wrong_guide = complex_structure.map_guide_to_rna_chain("AAAAAAAAAAAAAAAAAAAA", annotated["chains"])
    assert wrong_guide["status"] == "UNMAPPED"
    wrong_target = complex_structure.map_target_to_dna_chain(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", annotated["chains"]
    )
    assert wrong_target["status"] == "UNMAPPED"
    wrong_pam = complex_structure.map_pam_on_nontarget("CCC", annotated["chains"])
    assert wrong_pam["status"] == "UNMAPPED"
    prov = complex_structure.complex_provenance(annotated)
    assert prov["pdb"] == "4UN3"
    assert prov["method"] == "X-RAY DIFFRACTION"
    assert prov["resolution_angstrom"] == pytest.approx(2.593)
    assert prov["chains"]
    assert prov["sequence_identities"]
    assert prov["mapping_status"]["guide"] == "EXACT"
    source = Path(complex_structure.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    assert "place_cas_guide_dna" not in names
    assert "synthetic_cas_complex" not in names


def test_cache_keys_distinguish_molecule_types():
    dna_key = protein_structure.cache_key(
        source="RCSB PDB", identifier="1BNA", chain="A", parameters={"molecule": "DNA"}
    )
    prot_key = protein_structure.cache_key(
        source="RCSB PDB", identifier="1BNA", chain="A", parameters={"molecule": "PROTEIN"}
    )
    rna_key = protein_structure.cache_key(
        source="RCSB PDB", identifier="1BNA", chain="A", parameters={"molecule": "RNA"}
    )
    complex_key = protein_structure.cache_key(
        source="RCSB PDB", identifier="1BNA", chain="A", parameters={"cache_molecule": "COMPLEX"}
    )
    assert len({dna_key, prot_key, rna_key, complex_key}) == 4
    again = protein_structure.cache_key(
        source="RCSB PDB", identifier="1BNA", chain="A", parameters={"molecule": "DNA"}
    )
    assert again == dna_key


def test_tree_msa_and_hash_use_1bna():
    result = _1bna_result()
    digest = result["sequence_hash"]
    available = molecule_selection.structure_availability_for_hash(
        digest, nucleic_result=protein_structure.structure_3d_input(result)
    )
    assert available["status"] == "AVAILABLE"
    missing = molecule_selection.structure_availability_for_hash(
        provenance.sequence_digest("ATGCATGCATGC"),
        nucleic_result=protein_structure.structure_3d_input(result),
    )
    assert missing["status"] == "UNAVAILABLE"
    assert missing["reason"] == "No validated structure available."
    col = protein_structure.map_msa_column_to_structure_residue([None, 0, 1], 0, result)
    assert col is None
    mapped_col = protein_structure.map_msa_column_to_structure_residue([None, 0, 1], 2, result)
    assert mapped_col is not None
    assert mapped_col["query_index_0based"] == 1
    h1 = result["content_hash"]
    h2 = protein_structure.parse_mmcif(_load("1BNA.cif"))["content_hash"]
    assert h1 == h2
    scene = structure_scene.build_scene(
        protein_structure.structure_3d_input(result),
        expected_sequence_hash="deadbeef",
        representation="backbone",
    )
    assert scene["status"] == "STALE"
    assert scene["n_atoms_rendered"] == 0


def test_lod_step_is_not_stride_algorithm_and_counts_decrease():
    assert region_nav.backbone_stride("high") == 1
    assert region_nav.backbone_stride("medium") == 2
    assert region_nav.backbone_stride("low") == 4
    source = Path(region_nav.__file__).read_text(encoding="utf-8")
    assert "LOD_BACKBONE_STEP" in source
    assert "nao e o programa Frishman/Argos" in source
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    high = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        lod="high",
    )
    low = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
        lod="low",
    )
    assert low["n_backbone_points"] < high["n_backbone_points"]
    assert structure_scene.scientific_snapshot(envelope) == structure_scene.scientific_snapshot(envelope)


def test_mapping_statuses_include_unmapped_and_partial():
    assert scientific_checks.mapping_status_is_declared("UNMAPPED")
    assert scientific_checks.mapping_status_is_declared("PARTIAL")
    assert set(protein_structure.MAPPING_STATUSES) == {
        "EXACT",
        "ALIGNED",
        "PARTIAL",
        "BEST_EFFORT",
        "UNMAPPED",
    }


def test_no_gpu_claim_in_new_modules():
    blob = Path(na_structure_catalog.__file__).read_text(encoding="utf-8").lower()
    assert "gpu optimized" not in blob
    assert "real-time" not in blob


def test_illustrative_dna_window_builds_plotly_scatter3d():
    seq = "ATGC" * 700
    representation = dna_3d.representation_for_sequence(
        seq,
        region_start=0,
        region_end=10,
        lod=region_nav.LOD_HIGH,
    )
    assert representation["kind"] == dna_3d.KIND_ILLUSTRATIVE
    envelope = representation["envelope"]
    assert isinstance(envelope, dict)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=str(representation["sequence_hash"]),
        representation="backbone",
        selected_query_index=0,
        lod=region_nav.LOD_HIGH,
    )
    assert scene["status"] == "READY"
    assert int(scene["n_backbone_points"]) >= 2
    figure = structure_viewer.figure_from_scene(scene)
    assert int(figure.layout.height or 0) == STRUCTURE_VIEWPORT_HEIGHT
    assert figure.layout.autosize is False
    assert len(figure.data) >= 1
    xs = list(figure.data[0].x or [])
    assert len(xs) >= 2
    assert all(item is not None for item in xs)
    hover = " ".join(str(item) for item in (figure.data[0].text or []))
    assert "kind illustrative" in hover.lower() or "illustrative" in hover.lower()


def test_bundled_1bna_is_experimental_not_illustrative():
    bundled = dna_3d.bundled_experimental_for_sequence(DICKERSON, "DNA")
    assert bundled is not None
    assert bundled["kind"] == dna_3d.KIND_EXPERIMENTAL
    assert bundled["structure_id"] == "1BNA"
    assert "1BNA.cif" in str(bundled.get("bundled_fixture") or "")
    meta = bundled.get("metadata") or {}
    assert meta.get("method") == "X-RAY DIFFRACTION"
    assert meta.get("resolution_angstrom") == pytest.approx(1.9)
    auto = dna_3d.representation_for_sequence(DICKERSON)
    assert auto["kind"] == dna_3d.KIND_EXPERIMENTAL
    assert auto["status"] == "EXPERIMENTAL"
    scene = structure_scene.build_scene(
        auto["envelope"],
        expected_sequence_hash=str(auto["sequence_hash"]),
        representation="backbone",
    )
    assert scene["status"] == "READY"
    assert int(scene["n_backbone_points"]) >= 12
    chains = {str(trace.get("chain_id") or "") for trace in scene["objects"]["backbone_traces"]}
    assert "A" in chains and "B" in chains
    figure = structure_viewer.figure_from_scene(scene)
    finite = 0
    for trace in figure.data:
        finite += structure_viewer.count_finite_xyz(list(trace.x), list(trace.y), list(trace.z))
    assert finite > 0
    assert int(figure.layout.height or 0) == STRUCTURE_VIEWPORT_HEIGHT
    other = dna_3d.bundled_experimental_for_sequence("ATGCATGCATGC", "DNA")
    assert other is None


def test_bundled_1rna_is_experimental_and_mfe_stays_separate():
    bundled = dna_3d.bundled_experimental_for_sequence(RNA_DUPLEX, "RNA")
    assert bundled is not None
    assert bundled["kind"] == rna_3d.KIND_EXPERIMENTAL
    auto = rna_3d.representation_for_sequence(RNA_DUPLEX)
    assert auto["kind"] == rna_3d.KIND_EXPERIMENTAL
    pairing = rna_3d.pairing_sync(auto, None)
    assert pairing["geometry_kind"] == rna_3d.KIND_EXPERIMENTAL
    assert pairing["pairing_status"] == "UNAVAILABLE"
    assert "MFE" in str(pairing.get("note") or "")
    scene = structure_scene.build_scene(
        auto["envelope"],
        expected_sequence_hash=str(auto["sequence_hash"]),
    )
    assert scene["status"] == "READY"
    figure = structure_viewer.figure_from_scene(scene)
    assert len(figure.data) >= 1


def test_figure_rejects_mismatched_xyz_arrays():
    with pytest.raises(ValueError, match="Insufficient data"):
        structure_viewer.count_finite_xyz([1.0, 2.0], [1.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="Insufficient data"):
        structure_viewer.count_finite_xyz([float("nan")], [0.0], [0.0])
