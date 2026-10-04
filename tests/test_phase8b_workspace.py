"""Fase 8B: DNA/RNA ilustrativo, complexo depositado, MSA FASTA, ligacoes.

Nenhum teste inventa coordenadas experimentais. parser_nucleic_mini.cif nao
e uma deposicao PDB. 1CRN.cif permanece a estrutura experimental de crambin.
"""

from __future__ import annotations

from pathlib import Path

from modules import (
    complex_structure,
    dna_analysis,
    molecule_selection,
    msa,
    nucleic_geometry,
    protein_structure,
    provenance,
    structure_scene,
)

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _fail_network(*_args, **_kwargs):
    raise AssertionError("network")


def test_parse_unaligned_fasta_keeps_protein_multifasta():
    text = (
        ">seq1_humano\nMVLSPADKTNVKAAW\n"
        ">seq2_chimpunze\nMVLSPADKTNVKAAW\n"
        ">seq3_gorila\nMVLSPADKTNVKAAW\n"
    )
    records = msa.parse_unaligned_fasta(text)
    assert len(records) == 3
    assert records[0]["identifier"] == "seq1_humano"
    assert records[2]["identifier"] == "seq3_gorila"
    assert records[1]["sequence"] == "MVLSPADKTNVKAAW"


def test_parse_unaligned_fasta_survives_missing_dna_helper(monkeypatch):
    monkeypatch.setattr(msa.dna_analysis, "parse_fasta_records", None, raising=False)
    monkeypatch.setattr(msa.importlib, "reload", lambda module: module)
    records = msa.parse_unaligned_fasta(
        ">seq1_humano\nMVLSPADKTNVKAAW\n>seq2_chimpunze\nMVLSPADKTNVKAAW\n"
    )
    assert len(records) == 2
    assert records[0]["identifier"] == "seq1_humano"


def test_illustrative_bdna_is_not_experimental():
    model = nucleic_geometry.illustrative_bdna_model("ATGCATGCAT")
    assert model["kind"] == nucleic_geometry.KIND_ILLUSTRATIVE
    assert model["status"] == "ILLUSTRATIVE"
    assert "experimental" not in str(model["kind"]).lower()
    assert model["residues"][0]["base"] == "A"
    assert model["residues"][1]["base"] == "T"
    envelope = nucleic_geometry.structure_3d_input_from_illustrative(model)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=model["sequence_hash"],
        representation="backbone",
    )
    assert scene["status"] == "READY"
    assert envelope["kind"] == "illustrative"
    snapshot = structure_scene.scientific_snapshot(envelope)
    again = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=model["sequence_hash"],
        representation="backbone",
    )
    assert structure_scene.scientific_snapshot(envelope) == snapshot
    assert again["status"] == "READY"


def test_illustrative_arna_does_not_use_mfe():
    model = nucleic_geometry.illustrative_arna_model("AUGCAUGC")
    assert "mfe" not in str(model.get("algorithm") or "").lower()
    assert model["kind"] == "illustrative"
    source = Path(nucleic_geometry.__file__).read_text(encoding="utf-8")
    assert "ViennaRNA" in source
    assert "nao converte dot-bracket" in source.lower() or "not" in model["disclaimer"].lower()


def test_complex_roles_from_deposited_descriptions():
    parsed = protein_structure.parse_structure_text(_load("parser_nucleic_mini.cif"), hint="mmcif")
    annotated = complex_structure.annotate_parsed_complex(
        parsed,
        structure_id="NA01",
        source="fixture",
        metadata={"kind": "experimental"},
        guide_sequence="GG",
        pam="NGG",
    )
    roles = {item["chain_id"]: item["role"] for item in annotated["chains"]}
    assert roles["A"] == complex_structure.ROLE_GUIDE_RNA
    assert roles["B"] == complex_structure.ROLE_CAS_PROTEIN
    assert roles["C"] == complex_structure.ROLE_TARGET_DNA
    assert roles["D"] == complex_structure.ROLE_NONTARGET_DNA
    assert annotated["synthetic_complex"] is False
    assert annotated["guide_mapping"]["status"] == "EXACT"
    assert annotated["kind"] == "experimental"


def test_guide_unmapped_when_spacer_absent():
    parsed = protein_structure.parse_structure_text(_load("parser_nucleic_mini.cif"), hint="mmcif")
    annotated = complex_structure.annotate_parsed_complex(
        parsed,
        structure_id="NA01",
        source="fixture",
        metadata={"kind": "experimental"},
        guide_sequence="AAAAAAAAAA",
        pam="NGG",
    )
    assert annotated["guide_mapping"]["status"] == "UNMAPPED"
    assert annotated["guide_mapping"]["start_0based"] is None


def test_nucleic_backbone_uses_p_not_invented_ca():
    parsed = protein_structure.parse_structure_text(_load("parser_nucleic_mini.cif"), hint="mmcif")
    dna_chain = next(item for item in parsed["chains"] if item["chain_id"] == "C")
    mapped = protein_structure.map_query_to_chain("AT", dna_chain, molecule="DNA")
    assert mapped["mapping"][0]["p"] is not None
    assert mapped["mapping"][0]["ca"] is None
    result = protein_structure.load_protein_structure(
        "AT",
        source="RCSB PDB",
        structure_id="NA01",
        structure_text=_load("parser_nucleic_mini.cif"),
        metadata={"kind": "experimental", "method": "X-RAY DIFFRACTION"},
        urlopen_fn=_fail_network,
        molecule="DNA",
        chain_id="C",
    )
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
    )
    assert scene["status"] == "READY"


def test_ca_contact_map_uses_deposited_coordinates():
    import json

    meta = protein_structure.parse_rcsb_entry_metadata(
        json.loads(_load("rcsb_entry_1CRN.json"))
    )
    result = protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
    )
    contacts = protein_structure.ca_contact_map(result, threshold_angstrom=8.0)
    assert contacts["status"] == "COMPUTED"
    assert contacts["threshold_angstrom"] == 8.0
    assert contacts["n_pairs"] >= 1
    assert "proprietary" not in contacts["method"].lower()
    assert "euclidean" in contacts["method"].lower()


def test_dssp_absent_is_unavailable(monkeypatch):
    monkeypatch.setattr(
        protein_structure.tool_detection,
        "resolve_allowlisted_executable",
        lambda *_args, **_kwargs: None,
    )
    info = protein_structure.dssp_availability()
    assert info["available"] is False
    assert "not installed" in info["reason"].lower()


def test_tree_structure_availability_requires_matching_hash():
    missing = molecule_selection.structure_availability_for_hash(
        "abc",
        protein_result={"sequence_hash": "zzz", "kind": "experimental", "structure_id": "1CRN"},
    )
    assert missing["status"] == "UNAVAILABLE"
    assert missing["reason"] == "No validated structure available."
    present = molecule_selection.structure_availability_for_hash(
        "abc",
        protein_result={
            "sequence_hash": "abc",
            "kind": "experimental",
            "structure_id": "1CRN",
            "source": "RCSB PDB",
        },
    )
    assert present["status"] == "AVAILABLE"
    assert present["structure_id"] == "1CRN"


def test_complex_scene_does_not_fabricate_atoms():
    parsed = protein_structure.parse_structure_text(_load("parser_nucleic_mini.cif"), hint="mmcif")
    annotated = complex_structure.annotate_parsed_complex(
        parsed,
        structure_id="NA01",
        source="fixture",
        metadata={"kind": "experimental"},
    )
    envelope = complex_structure.deposited_to_structure_3d_input(annotated)
    assert envelope["synthetic_complex"] is False
    assert envelope["n_atoms"] == parsed["n_atoms"]
    scene = structure_scene.build_scene(envelope, representation="backbone")
    assert scene["status"] == "READY"


def test_1crn_protein_scene_still_ready():
    import json

    meta = protein_structure.parse_rcsb_entry_metadata(
        json.loads(_load("rcsb_entry_1CRN.json"))
    )
    result = protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
    )
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="backbone",
    )
    assert scene["status"] == "READY"
    assert envelope["kind"] == "experimental"


def test_parse_fasta_records_still_exported():
    parsed = dna_analysis.parse_fasta_records(">one\nAAAA\n>two\nCCCC\n")
    assert parsed["record_count"] == 2
    assert provenance.HELIXSCOPE_VERSION.startswith("0.")


def test_msa_tab_add_protein_multifasta_does_not_crash():
    from tests.helix_apptest import open_module

    fasta = (
        ">seq1_humano\nMVLSPADKTNVKAAW\n"
        ">seq2_chimpunze\nMVLSPADKTNVKAAW\n"
        ">seq3_gorila\nMVLSPADKTNVKAAW\n"
    )
    app = open_module("msa", timeout=60)
    assert not app.exception
    app.text_area(key="msa_fasta_text").set_value(fasta).run()
    app.button(key="msa_add_fasta").click().run()
    assert not app.exception
    collection = app.session_state["msa_collection"]
    assert len(collection) == 3
    assert collection[0]["identifier"] == "seq1_humano"
    assert collection[2]["identifier"] == "seq3_gorila"
    errors = [str(item.value) for item in list(app.error) or []]
    assert errors == []
