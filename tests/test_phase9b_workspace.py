"""Fase 9B: caminho DNA 3D estavel, RNA 3D honesto, mappings, escala, LOD.

Nenhum teste inventa coordenadas experimentais, arvores, genomas ou complexos
Cas. parser_nucleic_mini.cif nao e uma deposicao PDB. 1CRN.cif e crambin
publico. O registro NCBI dos testes de feature e sintetico (NC_000000).
"""

from __future__ import annotations

import ast
import inspect
import json
import math
import time
from pathlib import Path
from types import MappingProxyType

import pytest

from modules import (
    complex_structure,
    dna_3d,
    dna_analysis,
    molecule_selection,
    msa,
    ncbi_fetch,
    nucleic_geometry,
    protein_structure,
    provenance,
    region_nav,
    rna_3d,
    rna_folding,
    scale_profile,
    structure_scene,
)
from ui import charts, components, structure_viewer

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
HAIRPIN_RNA = "GGGAAACCC"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _fail_network(*_args, **_kwargs):
    raise AssertionError("network")


def _seq_50k() -> str:
    return "ACGT" * 12_500


def _1crn_result() -> dict:
    meta = protein_structure.parse_rcsb_entry_metadata(json.loads(_load("rcsb_entry_1CRN.json")))
    return protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=_fail_network,
    )


def _dna_experimental_at() -> dict:
    return protein_structure.load_protein_structure(
        "AT",
        source="parser fixture",
        structure_id="helixscope_na_mini",
        structure_text=_load("parser_nucleic_mini.cif"),
        metadata={"kind": "experimental", "method": "X-RAY DIFFRACTION"},
        urlopen_fn=_fail_network,
        molecule="DNA",
        chain_id="C",
    )


def test_default_dna_3d_not_unavailable_small_medium_50k():
    samples = {
        "small": "ATGCATGCAT",
        "medium": "ATGC" * 80,
        "50k": _seq_50k(),
    }
    for _label, seq in samples.items():
        rep = dna_3d.representation_for_sequence(seq)
        assert rep["kind"] != dna_3d.KIND_UNAVAILABLE
        assert rep["status"] == "ILLUSTRATIVE"
        assert rep["envelope"] is not None
        scene = dna_3d.scene_from_representation(rep)
        assert scene["status"] == "READY"


def test_ten_valid_dna_sequences_never_unavailable():
    seqs = [
        "ATGCATGCAT",
        "G" * 20,
        "ACGT" * 25,
        "TTTTAAAA" * 8,
        "GAATTC" + "C" * 20,
        "ATG" + "C" * 147,
        "A" * 200,
        "CG" * 150,
        "ATGC" * 200,
        _seq_50k(),
    ]
    assert len(seqs) == 10
    for seq in seqs:
        info = dna_analysis.validate_for_molecule(seq, "DNA")
        assert info["is_valid"] is True
        rep = dna_3d.representation_for_sequence(seq)
        assert rep["kind"] != dna_3d.KIND_UNAVAILABLE, seq[:20]


def test_four_dna_3d_kinds_are_reachable():
    seq = "ATGCATGCAT"
    illustrative = dna_3d.representation_for_sequence(seq)
    assert illustrative["kind"] == dna_3d.KIND_ILLUSTRATIVE
    predicted_missing = dna_3d.representation_for_sequence(seq, prefer="predicted")
    assert predicted_missing["kind"] == dna_3d.KIND_UNAVAILABLE
    assert predicted_missing["status"] == "UNAVAILABLE"
    experimental_missing = dna_3d.representation_for_sequence(seq, prefer="experimental")
    assert experimental_missing["kind"] == dna_3d.KIND_UNAVAILABLE
    loaded = _dna_experimental_at()
    experimental = dna_3d.representation_for_sequence(
        "AT", deposited=loaded, prefer="experimental"
    )
    assert experimental["kind"] == dna_3d.KIND_EXPERIMENTAL
    predicted_attached = dict(protein_structure.structure_3d_input(loaded))
    predicted_attached["kind"] = dna_3d.KIND_PREDICTED
    predicted_attached["sequence_hash"] = provenance.sequence_digest("AT")
    predicted = dna_3d.representation_for_sequence(
        "AT", deposited=predicted_attached, prefer="predicted"
    )
    assert predicted["kind"] == dna_3d.KIND_PREDICTED
    assert predicted["status"] == "PREDICTED"
    assert set(dna_3d.DNA_3D_KINDS) == {
        dna_3d.KIND_EXPERIMENTAL,
        dna_3d.KIND_PREDICTED,
        dna_3d.KIND_ILLUSTRATIVE,
        dna_3d.KIND_UNAVAILABLE,
    }


def test_illustrative_bases_match_input_and_window_contiguous():
    seq = "ATGCGAATTCGG"
    rep = dna_3d.representation_for_sequence(seq, region_start=4, region_end=10)
    check = dna_3d.window_integrity(rep, seq)
    assert check["ok"] is True
    assert check["fragment"] == seq[4:10]
    assert check["fragment"] == "GAATTC"
    assert check["n"] == check["end"] - check["start"]
    assert rep["kind"] == dna_3d.KIND_ILLUSTRATIVE
    assert "experimental" not in str(rep["kind"]).lower()


def test_helix_has_two_strands_and_pairs():
    seq = "ATGCATGC"
    rep = dna_3d.representation_for_sequence(seq)
    assert rep["n_strands"] == 2
    assert rep["n_pairs"] == len(seq)
    scene = dna_3d.scene_from_representation(rep)
    assert scene["n_pair_traces"] == len(seq)
    chains = list(rep["envelope"]["chains"])
    assert {item["chain_id"] for item in chains} == {"A", "B"}


def test_region_edges_defined_behavior():
    seq = "ATGC" * 10
    start_zero = dna_3d.representation_for_sequence(seq, region_start=0, region_end=8)
    assert start_zero["kind"] == dna_3d.KIND_ILLUSTRATIVE
    assert start_zero["view_start"] == 0
    last = dna_3d.representation_for_sequence(seq, region_start=len(seq) - 4, region_end=len(seq))
    assert last["fragment"] == seq[-4:]
    size_one = dna_3d.representation_for_sequence(seq, region_start=0, region_end=1)
    assert size_one["status"] == "ERROR"
    full = dna_3d.representation_for_sequence(seq, region_start=0, region_end=len(seq))
    assert full["fragment"] == seq
    with pytest.raises(dna_3d.Dna3dError):
        dna_3d.navigate(len(seq), action="region", start=10, end=4)
    inverted = dna_3d.representation_for_sequence(seq, region_start=10, region_end=4)
    assert inverted["status"] == "ERROR"


def test_50k_overview_object_count_is_capped():
    seq = _seq_50k()
    budget = dna_3d.overview_object_budget(len(seq), region_nav.LOD_HIGH)
    assert budget["scales_linearly_with_full_length"] is False
    assert budget["max_nt_atomic_detail"] == region_nav.MAX_HELIX_NT_HIGH
    rep = dna_3d.representation_for_sequence(seq, lod=region_nav.LOD_HIGH)
    assert rep["n_scene_objects"] <= budget["max_atoms"] + budget["max_pairs"]
    assert rep["n_scene_objects"] < 5_000
    assert rep["n_residues_drawn"] == region_nav.MAX_HELIX_NT_HIGH
    assert rep["partial_view"] is True


def test_navigate_50k_overview_jump_region_zoom_reset():
    length = 50_000
    overview = dna_3d.navigate(length, action="overview")
    assert overview["partial_view"] is True
    assert overview["n"] == region_nav.MAX_HELIX_NT_HIGH
    jump = dna_3d.navigate(length, action="jump", position=25_000)
    assert jump["start"] <= 25_000 < jump["end"]
    region = dna_3d.navigate(length, action="region", start=20_000, end=20_400)
    assert region["start"] == 20_000
    assert region["end"] - region["start"] == region["n"]
    zoom = dna_3d.navigate(
        length, action="zoom", start=region["start"], end=region["end"], lod=region_nav.LOD_LOW
    )
    assert zoom["n"] <= region_nav.MAX_HELIX_NT_LOW
    seq = _seq_50k()
    local = dna_3d.representation_for_sequence(
        seq, region_start=zoom["start"], region_end=zoom["end"], lod=region_nav.LOD_LOW
    )
    assert local["kind"] == dna_3d.KIND_ILLUSTRATIVE
    reset = dna_3d.navigate(length, action="reset")
    assert reset["start"] == overview["start"]
    assert reset["end"] == overview["end"]


def test_ecori_motif_highlight_exact_span():
    with_site = "CCCCGAATTCCCCC"
    without = "CCCCCCCCCCCCCC"
    hits = dna_3d.highlight_from_motif(with_site, "GAATTC")
    assert hits["n_hits"] >= 1
    assert with_site[hits["spans"][0]["start"] : hits["spans"][0]["end"]] == "GAATTC"
    empty = dna_3d.highlight_from_motif(without, "GAATTC")
    assert empty["indices"] == []
    assert empty["n_hits"] == 0


def test_crispr_highlight_requires_real_spacer_and_adjacent_pam():
    spacer = "GAGTCCGAGCAGAAGAAGAA"
    seq = "ATGC" + spacer + "GGGTT"
    guide = {
        "guide_sequence": spacer,
        "start_0based": 4,
        "pam_sequence": "NGG",
        "pam_start_0based": 24,
        "strand": "+",
    }
    mapped = dna_3d.highlight_from_guide(seq, guide)
    assert mapped["status"] == "EXACT"
    assert seq[4:24] == spacer
    assert seq[24:27] == "GGG"
    missing_pam = dict(guide)
    missing_pam["pam_start_0based"] = 24
    bad_seq = "ATGC" + spacer + "AAATT"
    unmapped = dna_3d.highlight_from_guide(bad_seq, missing_pam)
    assert unmapped["status"] == "UNMAPPED"
    assert unmapped["indices"] == []


def test_ncbi_feature_spans_map_to_highlight_synthetic_record():
    record = {
        "accession": "NC_000000.1",
        "description": "Synthetic record for tests",
        "organism": "Synthetic organism",
        "length": 1000,
        "sequence": "A" * 1000,
        "features": [
            {
                "type": "CDS",
                "location": "[10:100](+)",
                "qualifiers": {"gene": ["orf1ab"]},
            }
        ],
    }
    items = ncbi_fetch.extract_annotated_features(record)
    simple = next(item for item in items if item["type"] == "CDS")
    spans = ncbi_fetch.feature_spans_for_analysis(simple, record["sequence"])
    mapped = dna_3d.highlight_from_ncbi_spans(
        record["sequence"], spans, accession=record["accession"]
    )
    assert 10 in mapped["indices"]
    assert 99 in mapped["indices"]
    assert 9 not in mapped["indices"]
    assert 100 not in mapped["indices"]


def test_msa_gap_column_is_unmapped():
    members = [
        msa.build_collection_member(sequence="ACGT", identifier="seq_a", molecule="DNA"),
        msa.build_collection_member(sequence="ACCGT", identifier="seq_b", molecule="DNA"),
    ]
    result = msa.build_msa_result(
        raw_alignment=_load("msa_aligned_gaps.fa"),
        members=members,
        tool="fixture",
        tool_version="test",
        method="prealigned fixture",
        parameters={},
    )
    gap = dna_3d.highlight_from_msa_column(result, 2, "ACGT")
    assert gap["status"] == "UNMAPPED"
    assert gap["indices"] == []
    mapped = dna_3d.highlight_from_msa_column(result, 3, "ACGT")
    assert mapped["indices"] == [2]
    assert "ACGT"[2] == "G"


def test_graph_window_boundary_is_unambiguous():
    seq = "A" * 200
    profile = dna_analysis.gc_sliding_window(seq, window=100, step=100)
    assert profile[0]["end"] == 100
    first = dna_3d.highlight_from_window(len(seq), profile[0]["start"], profile[0]["end"])
    second = dna_3d.highlight_from_window(len(seq), profile[1]["start"], profile[1]["end"])
    assert 99 in first["indices"]
    assert 100 not in first["indices"]
    assert 100 in second["indices"]
    assert 99 not in second["indices"]


def test_export_unavailable_is_not_success():
    missing = dna_3d.representation_for_sequence("ATGCATGC", prefer="experimental")
    exported = dna_3d.export_representation(missing)
    assert exported["status"] in {"UNAVAILABLE", "ERROR"}
    assert exported["successful_result"] is False
    for field in ("source", "method", "model", "version", "hashes", "status", "limitations"):
        assert field in exported


def test_export_stale_after_sequence_change():
    seq = "ATGCATGCAT"
    rep = dna_3d.representation_for_sequence(seq)
    other = provenance.sequence_digest("GGGGGGGGGG")
    exported = dna_3d.export_representation(rep, expected_sequence_hash=other)
    assert exported["status"] == "STALE"
    assert exported["successful_result"] is False


def test_lod_omits_points_without_changing_remaining_xyz_dna():
    seq = "ATGC" * 50
    rep = dna_3d.representation_for_sequence(seq)
    envelope = rep["envelope"]
    snapshot = structure_scene.scientific_snapshot(envelope)
    high = structure_scene.build_scene(
        envelope, expected_sequence_hash=rep["sequence_hash"], lod="high"
    )
    low = structure_scene.build_scene(
        envelope, expected_sequence_hash=rep["sequence_hash"], lod="low"
    )
    hx = high["objects"]["backbone_traces"][0]["x"]
    lx = low["objects"]["backbone_traces"][0]["x"]
    assert len(lx) < len(hx)
    high_set = {
        (round(float(x), 6), round(float(y), 6), round(float(z), 6))
        for x, y, z in zip(
            high["objects"]["backbone_traces"][0]["x"],
            high["objects"]["backbone_traces"][0]["y"],
            high["objects"]["backbone_traces"][0]["z"],
        )
    }
    for x, y, z in zip(
        low["objects"]["backbone_traces"][0]["x"],
        low["objects"]["backbone_traces"][0]["y"],
        low["objects"]["backbone_traces"][0]["z"],
    ):
        assert (round(float(x), 6), round(float(y), 6), round(float(z), 6)) in high_set
    assert structure_scene.scientific_snapshot(envelope) == snapshot


def test_bidirectional_sequence_3d_selection():
    seq = "ATGCGAATTC"
    rep = dna_3d.representation_for_sequence(seq)
    forward = dna_3d.position_to_3d(rep, 3)
    assert forward["status"] == "READY"
    assert forward["mapping"]["query_residue"] == seq[3]
    pick = [3, "A", 1, 4, 4, "", "P"]
    reverse = dna_3d.pick_to_position(pick)
    assert reverse["status"] == "READY"
    assert reverse["query_index_0based"] == 3
    assert dna_3d.position_to_3d(rep, 99)["status"] == "UNMAPPED"


def test_detailed_mode_highlights():
    seq = "CCCCGAATTCCCCC"
    motif = dna_3d.highlight_from_motif(seq, "GAATTC")
    rep = dna_3d.representation_for_sequence(seq, highlight_indices=motif["indices"])
    residues = list(rep["model"]["residues"])
    highlighted = [item["index_0based"] for item in residues if item.get("highlighted")]
    assert highlighted == motif["indices"]
    scene = dna_3d.scene_from_representation(
        rep, selected_query_index=4, highlight_query_indices=motif["indices"]
    )
    assert scene["selected_query_index"] == 4
    assert 4 in scene["selection_ids"]
    labels = dna_3d.base_labels_for_view(rep)
    assert {item["base"] for item in labels} <= set(seq)
    assert {item["chain_id"] for item in rep["envelope"]["chains"]} == {"A", "B"}
    assert all(item.get("kind") == dna_3d.KIND_ILLUSTRATIVE for item in rep["model"]["pairs"])


def test_base_labels_capped_independent_of_length():
    long_seq = "ATGC" * 200
    rep = dna_3d.representation_for_sequence(long_seq)
    labels = dna_3d.base_labels_for_view(rep)
    assert len(labels) <= dna_3d.MAX_BASE_LABELS
    assert len(rep["model"]["residues"]) > dna_3d.MAX_BASE_LABELS


def test_round_trip_window_concat_matches_sequence():
    seq = "ATGC" * 10
    chunks = []
    step = 8
    for start in range(0, len(seq), step):
        end = min(len(seq), start + step)
        if end - start < 2:
            break
        rep = dna_3d.representation_for_sequence(seq, region_start=start, region_end=end)
        chunks.append(rep["fragment"])
    rebuilt = "".join(chunks)
    assert rebuilt == seq[: len(rebuilt)]
    assert rebuilt == seq


def test_rna_four_kinds_and_mfe_not_experimental():
    seq = "AUGCAUGC"
    illustrative = rna_3d.representation_for_sequence(seq)
    assert illustrative["kind"] == rna_3d.KIND_ILLUSTRATIVE
    assert illustrative["mfe_does_not_yield_3d"] is True
    predicted_missing = rna_3d.representation_for_sequence(seq, prefer="predicted")
    assert predicted_missing["status"] == "UNAVAILABLE"
    assert predicted_missing["kind"] == rna_3d.KIND_UNAVAILABLE
    experimental_missing = rna_3d.representation_for_sequence(seq, prefer="experimental")
    assert experimental_missing["kind"] == rna_3d.KIND_UNAVAILABLE
    deposited = dict(illustrative["envelope"])
    deposited["kind"] = rna_3d.KIND_PREDICTED
    deposited["sequence_hash"] = illustrative["sequence_hash"]
    predicted = rna_3d.representation_for_sequence(seq, deposited=deposited, prefer="predicted")
    assert predicted["kind"] == rna_3d.KIND_PREDICTED
    exp_dep = dict(illustrative["envelope"])
    exp_dep["kind"] = rna_3d.KIND_EXPERIMENTAL
    exp_dep["sequence_hash"] = illustrative["sequence_hash"]
    experimental = rna_3d.representation_for_sequence(seq, deposited=exp_dep, prefer="experimental")
    assert experimental["kind"] == rna_3d.KIND_EXPERIMENTAL
    source = Path(rna_3d.__file__).read_text(encoding="utf-8")
    assert "ViennaRNA MFE" in source
    assert "mfe_does_not_yield_3d" in source


def test_rna_pairing_origin_propagates():
    seq = HAIRPIN_RNA
    geom = rna_3d.representation_for_sequence(seq)
    none = rna_3d.pairing_sync(geom, None)
    assert none["pairing_status"] == "UNAVAILABLE"
    assert none["synced"] is False
    assert none["geometry_kind"] == rna_3d.KIND_ILLUSTRATIVE
    body = _load("rnafold_output_hairpin.txt")
    parsed = rna_folding.parse_rnafold_output(body, seq)
    pairs = rna_folding.parse_dot_bracket(parsed["structure"])
    mfe = {
        "status": "PREDICTED",
        "sequence_hash": provenance.sequence_digest(seq),
        "tool": "RNAfold fixture",
        "base_pairs": [
            {"position_0based": left, "paired_position_0based": right}
            for left, right in pairs
        ],
    }
    synced = rna_3d.pairing_sync(geom, mfe)
    assert synced["pairing_status"] == "PREDICTED"
    assert synced["geometry_kind"] == rna_3d.KIND_ILLUSTRATIVE
    assert synced["synced"] is True
    illustrative_mfe = dict(mfe)
    illustrative_mfe["status"] = "ILLUSTRATIVE"
    blocked = rna_3d.pairing_sync(geom, illustrative_mfe)
    assert blocked["pairing_status"] == "UNAVAILABLE"


def test_1crn_contact_euclidean_spot_check():
    result = _1crn_result()
    contacts = protein_structure.ca_contact_map(result, threshold_angstrom=8.0)
    assert contacts["threshold_angstrom"] == 8.0
    assert contacts["status"] == "COMPUTED"
    pair = contacts["pairs"][0]
    by_index = {
        int(item["query_index_0based"]): item["ca"]
        for item in result["residue_mapping"]
        if item.get("ca")
    }
    a = by_index[pair["i"]]
    b = by_index[pair["j"]]
    dist = math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)
    assert dist == pytest.approx(pair["distance_angstrom"])
    assert dist <= 8.0


def test_ca_representation_one_point_per_residue():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    ca = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="ca",
    )
    n_ca = sum(1 for item in result["residue_mapping"] if item.get("ca"))
    assert len(ca["objects"]["atom_points"]["x"]) == n_ca
    atoms = structure_scene.build_scene(
        envelope,
        expected_sequence_hash=result["sequence_hash"],
        representation="atoms",
    )
    assert len(atoms["objects"]["atom_points"]["x"]) == len(result["atoms"])


def test_no_synthetic_cas_complex():
    parsed = protein_structure.parse_structure_text(_load("parser_nucleic_mini.cif"), hint="mmcif")
    annotated = complex_structure.annotate_parsed_complex(
        parsed,
        structure_id="NA01",
        source="fixture",
        metadata={"kind": "experimental"},
        guide_sequence="GG",
        pam="NGG",
    )
    assert annotated["synthetic_complex"] is False
    source = Path(complex_structure.__file__).read_text(encoding="utf-8")
    assert "synthetic_complex" in source
    tree = ast.parse(source)
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    assert "place_cas_guide_dna" not in names
    assert "synthetic_cas_complex" not in names


def test_freeze_structure_3d_input_rejects_mutation():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    frozen = protein_structure.freeze_structure_3d_input(envelope)
    assert isinstance(frozen, MappingProxyType)
    with pytest.raises(TypeError):
        frozen["kind"] = "illustrative"
    envelope["kind"] = "illustrative"
    assert protein_structure.structure_3d_input(result)["kind"] == "experimental"


def test_selection_hash_stable_and_sensitive():
    digest = provenance.sequence_digest("ATGC")
    left = molecule_selection.make_selection(
        molecule="DNA", sequence_hash=digest, query_indices=[1, 2], source="sequence"
    )
    right = molecule_selection.make_selection(
        molecule="DNA", sequence_hash=digest, query_indices=[2, 1], source="sequence"
    )
    assert left["selection_hash"] == right["selection_hash"]
    changed = molecule_selection.make_selection(
        molecule="DNA", sequence_hash=digest, query_indices=[1], source="sequence"
    )
    assert changed["selection_hash"] != left["selection_hash"]


def test_tree_unmapped_leaf_is_unavailable():
    missing = molecule_selection.structure_availability_for_hash(
        provenance.sequence_digest("ATGCATGC"),
        protein_result=None,
        nucleic_result=None,
        complex_result=None,
    )
    assert missing["status"] == "UNAVAILABLE"
    assert missing["reason"] == "No validated structure available."


def test_protein_1crn_and_alphafold_p01308_kinds():
    result = _1crn_result()
    envelope = protein_structure.structure_3d_input(result)
    scene = structure_scene.build_scene(
        envelope, expected_sequence_hash=result["sequence_hash"]
    )
    assert result["kind"] == "experimental"
    assert scene["status"] == "READY"
    payload = json.loads(_load("alphafold_prediction_P01308.json"))
    parsed = protein_structure.parse_alphafold_prediction(payload)
    assert parsed["kind"] == "predicted"
    assert parsed["uniprot"] == "P01308"
    assert parsed["source"] == "AlphaFold DB"


def test_experimental_dna_coordinates_match_parser_cif():
    result = _dna_experimental_at()
    chain_c = [
        atom
        for atom in result["atoms"]
        if str(atom.get("auth_asym_id") or atom.get("label_asym_id")) == "C"
    ]
    assert chain_c[0]["x"] == pytest.approx(0.0)
    assert chain_c[0]["y"] == pytest.approx(8.0)
    assert chain_c[1]["x"] == pytest.approx(3.0)
    cif = _load("parser_nucleic_mini.cif")
    assert "0.000 8.000 0.000" in cif
    assert "3.000 8.000 0.000" in cif
    assert "Not a PDB deposition" in cif or "not a PDB deposition" in cif.lower()
    rep = dna_3d.representation_for_sequence("AT", deposited=result, prefer="experimental")
    assert rep["kind"] == dna_3d.KIND_EXPERIMENTAL


def test_sequence_display_position_45000():
    seq = _seq_50k()
    html = components.sequence_display(seq, "DNA", view_start=45_000, view_end=45_060)
    assert "Showing residues 45,000" in html
    assert seq[45_000] in html
    with pytest.raises(ValueError):
        components.sequence_display(seq, "DNA", view_start=50_000)
    with pytest.raises(ValueError):
        components.sequence_display(seq, "DNA", view_start=10, view_end=10)


def test_partial_view_badge_html():
    html = components.status_badge("PARTIAL_VIEW")
    assert "PARTIAL VIEW" in html
    seq = _seq_50k()
    rep = dna_3d.representation_for_sequence(seq)
    assert rep["partial_view"] is True


def test_infographics_still_render():
    seq = "ATGCATGCATGC"
    composition = dna_analysis.nucleotide_composition(seq)
    bar = charts.nucleotide_bar_chart(composition)
    assert bar.data
    profile = dna_analysis.gc_sliding_window(seq, window=4, step=2)
    gc_fig = charts.gc_sliding_window_plot(profile)
    assert len(gc_fig.data) >= 1
    hist = charts.orf_length_histogram([150, 300, 450])
    assert hist.data
    exported = list(profile)
    assert len(exported) == len(profile)


def test_renderer_has_no_network_or_science():
    viewer = Path(structure_viewer.__file__).read_text(encoding="utf-8")
    for banned in ("urllib", "requests", "httpx", "socket", "gc_content", "find_guides"):
        assert banned not in viewer
    scene_src = Path(structure_scene.__file__).read_text(encoding="utf-8")
    assert "from . import dna_analysis" not in scene_src
    assert "from . import crispr" not in scene_src
    assert "from . import msa" not in scene_src
    assert "from . import phylogeny" not in scene_src


def test_stale_triggers_are_isolated():
    seq = "ATGCATGCAT"
    dna_rep = dna_3d.representation_for_sequence(seq)
    other = "GGGGGGGGGG"
    stale_export = dna_3d.export_representation(
        dna_rep, expected_sequence_hash=provenance.sequence_digest(other)
    )
    assert stale_export["status"] == "STALE"
    members = [
        msa.build_collection_member(sequence="ACGT", identifier="seq_a", molecule="DNA"),
        msa.build_collection_member(sequence="ACCGT", identifier="seq_b", molecule="DNA"),
    ]
    msa_result = msa.build_msa_result(
        raw_alignment=_load("msa_aligned_gaps.fa"),
        members=members,
        tool="fixture",
        tool_version="test",
        method="prealigned fixture",
        parameters={},
    )
    stale_msa = dna_3d.highlight_from_msa_column(msa_result, 0, seq)
    assert stale_msa["status"] == "STALE"
    guide = {
        "guide_sequence": "AAAAAAAAAA",
        "start_0based": 0,
        "pam_sequence": "NGG",
        "pam_start_0based": 10,
        "strand": "+",
    }
    unmapped = dna_3d.highlight_from_guide(seq, guide)
    assert unmapped["status"] == "UNMAPPED"
    loaded = _dna_experimental_at()
    mismatch = dna_3d.representation_for_sequence(seq, deposited=loaded, prefer="experimental")
    assert mismatch["status"] in {"STALE", "UNAVAILABLE", "ERROR", "UNMAPPED"}


def test_50k_benchmark_validation_analysis_navigation_setup():
    rows = {}
    for n in (10_000, 25_000, 50_000):
        seq = ("ACGT" * ((n // 4) + 1))[:n]
        t0 = time.perf_counter()
        info = dna_analysis.validate_for_molecule(seq, "DNA")
        validation_ms = (time.perf_counter() - t0) * 1000.0
        t1 = time.perf_counter()
        dna_analysis.gc_content(seq)
        dna_analysis.nucleotide_composition(seq)
        analysis_ms = (time.perf_counter() - t1) * 1000.0
        t2 = time.perf_counter()
        dna_3d.navigate(n, action="jump", position=n // 2)
        navigation_ms = (time.perf_counter() - t2) * 1000.0
        t3 = time.perf_counter()
        dna_3d.representation_for_sequence(seq)
        viz_ms = (time.perf_counter() - t3) * 1000.0
        assert info["is_valid"] is True
        rows[n] = {
            "validation_ms": validation_ms,
            "analysis_ms": analysis_ms,
            "navigation_ms": navigation_ms,
            "visualization_setup_ms": viz_ms,
        }
        assert validation_ms < 15_000
        assert analysis_ms < 15_000
        assert viz_ms < 15_000
    assert set(rows) == {10_000, 25_000, 50_000}


def test_b_dna_parameters_are_documented():
    assert nucleic_geometry.B_DNA_RISE_A == pytest.approx(3.38)
    assert nucleic_geometry.B_DNA_TWIST_DEG == pytest.approx(36.0)
    source = inspect.getsource(nucleic_geometry)
    assert "3.38" in source
    assert "heuristica simplificada" in source.lower() or "heuristic" in source.lower()


def test_apptest_dna_3d_auto_path_sets_representation():
    from tests.helix_apptest import open_module

    app = open_module("dna", timeout=90)
    assert not app.exception
    app.text_area(key="dna_text").set_value("ATGGCATTACGTACGTACGTATGC").run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    app.toggle(key="dna_3d").set_value(True).run()
    assert not app.exception
    assert "dna_3d_representation" in app.session_state
    rep = app.session_state["dna_3d_representation"]
    assert isinstance(rep, dict)
    assert rep.get("kind") == dna_3d.KIND_ILLUSTRATIVE
    app.toggle(key="dna_illust_view_3d").set_value(True).run()
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("markdown", "caption", "info", "error")
        for item in list(getattr(app, bucket, []) or [])
    )
    assert "ILLUSTRATIVE" in combined or "illustrative" in combined.lower()
    assert "3D viewport (Plotly Scatter3d)" in combined
    assert "Insufficient data" not in combined
