"""Selecao molecular compartilhada: mapping, stale, gap e tipos distintos.

Nenhum teste inventa coordenadas. MSA de fixture e alinhamento de parser,
nao Clustal ao vivo. parser_two_chains.cif nao e uma deposicao.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from modules import molecule_selection, msa, protein_structure, provenance, scientific_checks

FIXTURES = Path(__file__).parent / "fixtures"


def _msa_protein_with_gap() -> tuple[dict, str]:
    sequence = "AK"
    other = "AGK"
    raw = ">seq_1\nA-K\n>seq_2\nAGK\n"
    members = [
        msa.build_collection_member(
            sequence=sequence, identifier="seq_1", source="user input", molecule="PROTEIN"
        ),
        msa.build_collection_member(
            sequence=other, identifier="seq_2", source="user input", molecule="PROTEIN"
        ),
    ]
    result = msa.build_msa_result(
        raw_alignment=raw,
        members=members,
        tool="fixture parser",
        tool_version="",
        method="aligned FASTA fixture, not a live Clustal Omega job",
        parameters={"fixture": True},
        source="test fixture",
    )
    return result, sequence


def test_molecules_are_not_interchangeable():
    assert molecule_selection.molecules_are_compatible("PROTEIN", "PROTEIN") is True
    assert molecule_selection.molecules_are_compatible("DNA", "PROTEIN") is False
    assert molecule_selection.molecules_are_compatible("DNA", "RNA") is False
    dna = molecule_selection.selection_from_sequence_position(
        molecule="DNA", sequence="ATGC", position=1
    )
    protein_hash = provenance.sequence_digest("AK")
    assert (
        molecule_selection.selection_applies_to(
            dna, sequence_hash=protein_hash, molecule="PROTEIN"
        )
        is False
    )


def test_empty_hashes_never_match_and_unknown_source_fails():
    selection = molecule_selection.make_selection(
        molecule="PROTEIN",
        sequence_hash="",
        query_indices=[0],
        source="sequence",
    )
    assert (
        molecule_selection.selection_applies_to(
            selection, sequence_hash="", molecule="PROTEIN"
        )
        is False
    )
    with pytest.raises(molecule_selection.SelectionError) as exc:
        molecule_selection.make_selection(
            molecule="PROTEIN",
            sequence_hash="abc",
            query_indices=[0],
            source="invented",
        )
    assert exc.value.category == "INVALID_INPUT"


def test_structure_pick_unmapped_does_not_invent_index():
    digest = provenance.sequence_digest("AK")
    unmapped = molecule_selection.selection_from_structure_pick(
        sequence_hash=digest,
        structure_id="1CRN",
        query_index=None,
        chain_id="A",
    )
    assert unmapped["status"] == "UNMAPPED"
    assert unmapped["query_indices"] == []
    assert "No mapped structural residue" in unmapped["message"]
    picked = molecule_selection.selection_from_structure_pick(
        sequence_hash=digest,
        structure_id="1CRN",
        structure_hash="hash-a",
        query_index=0,
        chain_id="A",
        atom_name="CA",
        label_seq_id=1,
        auth_seq_id=10,
    )
    assert picked["query_indices"] == [0]
    assert picked["source"] == "structure"
    other = molecule_selection.selection_applies_to(
        picked,
        sequence_hash=digest,
        molecule="PROTEIN",
        structure_id="9ZZZ",
    )
    assert other is False
    mismatch_hash = molecule_selection.selection_applies_to(
        picked,
        sequence_hash=digest,
        molecule="PROTEIN",
        structure_id="1CRN",
        structure_hash="hash-b",
    )
    assert mismatch_hash is False


def test_msa_gap_does_not_create_a_structure_residue():
    msa_result, sequence = _msa_protein_with_gap()
    gap = molecule_selection.selection_from_msa_column(
        msa_result, 1, sequence=sequence, molecule="PROTEIN"
    )
    assert gap["status"] == "UNMAPPED"
    assert gap["query_indices"] == []
    assert gap["message"] == molecule_selection.GAP_MESSAGE
    mapped = molecule_selection.selection_from_msa_column(
        msa_result, 0, sequence=sequence, molecule="PROTEIN"
    )
    assert mapped["status"] == "READY"
    assert mapped["query_indices"] == [0]
    disclaimer = str(mapped.get("conservation_disclaimer") or "").lower()
    assert "among analyzed sequences" in disclaimer
    assert "functional residue annotation" in disclaimer
    stale = molecule_selection.selection_from_msa_column(
        msa_result, 0, sequence="KR", molecule="PROTEIN"
    )
    assert stale["status"] == "STALE"


def test_motif_and_ncbi_spans_use_real_intervals():
    sequence = "AKDE"
    motif = molecule_selection.selection_from_motif_hit(
        {"start": 1, "end": 3, "match": "KD", "pattern": "KD"},
        molecule="PROTEIN",
        sequence=sequence,
    )
    assert motif["query_indices"] == [1, 2]
    assert motif["source"] == "motif"
    span = scientific_checks.make_span(
        start=0,
        end=2,
        sequence=sequence,
        source="NCBI",
        label="Region",
        kind="Region",
    )
    ncbi = molecule_selection.selection_from_ncbi_spans(
        [span], sequence=sequence, molecule="PROTEIN", accession="P00000"
    )
    assert ncbi["query_indices"] == [0, 1]
    assert ncbi.get("feature_source") == "NCBI"
    empty = molecule_selection.selection_from_ncbi_spans(
        [], sequence=sequence, molecule="PROTEIN"
    )
    assert empty["status"] == "UNAVAILABLE"


def test_selection_does_not_mutate_scientific_hashes():
    digest = provenance.sequence_digest("AK")
    selection = molecule_selection.make_selection(
        molecule="PROTEIN",
        sequence_hash=digest,
        query_indices=[0, 1],
        source="sequence",
    )
    selection["query_indices"].append(99)
    again = molecule_selection.make_selection(
        molecule="PROTEIN",
        sequence_hash=digest,
        query_indices=[0, 1],
        source="sequence",
    )
    assert again["query_indices"] == [0, 1]
    assert again["sequence_hash"] == digest


def test_two_chain_fixture_is_not_concatenated():
    text = (FIXTURES / "parser_two_chains.cif").read_text(encoding="utf-8")
    result = protein_structure.load_protein_structure(
        "AK",
        source="parser fixture",
        structure_id="parser_two_chains",
        structure_text=text,
        metadata={
            "kind": "experimental",
            "method": "X-ray",
            "structure_id": "parser_two_chains",
        },
        urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
    )
    chain_ids = {str(item.get("chain_id")) for item in result["chains"]}
    assert chain_ids == {"A", "B"}
    assert result["selected_chain"]["chain_id"] == "A"
    assert result["sequence"] == "AK"
    assert result["structure_sequence"] == "AK"
    mapped = protein_structure.query_position_to_residue(result, 0)
    assert mapped["chain_id"] == "A"
    b_comp = {
        str(atom.get("comp_id"))
        for atom in result["atoms"]
        if str(atom.get("auth_asym_id") or atom.get("label_asym_id")) == "B"
    }
    assert b_comp == {"GLY", "LEU"}
    assert result["mapping_status"] == "EXACT"


def test_tree_leaf_selection_uses_sequence_hash():
    from tests.test_phylogeny import _msa_from_aligned
    from modules import phylogeny

    msa_result = _msa_from_aligned(
        ["a", "b"],
        ["ACGTACGTAA", "ACGTTCGTAA"],
    )
    tree = phylogeny.infer_phylogeny(
        msa_result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    leaf = tree["leaves"][0]
    selection = molecule_selection.selection_from_tree_leaf(tree, leaf["tree_id"])
    assert selection["source"] == "tree"
    assert selection["sequence_hash"] == leaf["sequence_hash"]
    assert selection["alignment_hash"] == tree["alignment_hash"]
    stale = molecule_selection.selection_from_tree_leaf(
        tree, leaf["tree_id"], sequence="GGGGGGGGGG", molecule="DNA"
    )
    assert stale["status"] == "STALE"
    with pytest.raises(molecule_selection.SelectionError):
        molecule_selection.selection_from_tree_leaf(tree, "missing_leaf")
