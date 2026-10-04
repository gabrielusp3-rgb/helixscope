"""Testes da camada taxonomica, isolada da topologia filogenetica.

Nenhum teste inventa taxon, nome cientifico ou linhagem. Fetch NCBI e mockado.
"""

from __future__ import annotations

import pytest

from modules import ncbi_fetch, phylogeny, taxonomy
from tests.test_phylogeny import _msa_from_aligned


def _tree_with_organisms() -> dict:
    msa_result = _msa_from_aligned(
        ["human", "mouse", "user"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
        organisms=["Homo sapiens", "Mus musculus", ""],
        accessions=["NM_000001.1", "NM_000002.1", ""],
    )
    return phylogeny.infer_phylogeny(
        msa_result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )


def test_declared_organism_is_not_invented_for_user_sequence():
    tree = _tree_with_organisms()
    user = next(leaf for leaf in tree["leaves"] if leaf["original_id"] == "user")
    assert user["organism_display"] == "unknown / user-provided"
    assert user["taxonomy"] is None
    human = next(leaf for leaf in tree["leaves"] if leaf["original_id"] == "human")
    assert human["organism_declared"] == "Homo sapiens"
    assert human["accession"] == "NM_000001.1"


def test_duplicate_species_are_not_collapsed():
    msa_result = _msa_from_aligned(
        ["a", "b"],
        ["AAAAAAAAAA", "AAAAAAAATT"],
        organisms=["Homo sapiens", "Homo sapiens"],
        accessions=["NM_1.1", "NM_2.1"],
    )
    tree = phylogeny.infer_phylogeny(
        msa_result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    assert tree["n_leaves"] == 2
    assert [leaf["organism_declared"] for leaf in tree["leaves"]] == [
        "Homo sapiens",
        "Homo sapiens",
    ]
    assert tree["leaves"][0]["sequence_hash"] != tree["leaves"][1]["sequence_hash"]


def test_taxonomy_fetch_failure_does_not_destroy_tree():
    tree = _tree_with_organisms()
    original_hash = tree["tree_hash"]
    original_newick = tree["newick"]

    def boom(query, email, api_key=None):
        raise ncbi_fetch.NCBIQueryError("NCBI down", "source_unavailable")

    updated = taxonomy.attach_ncbi_taxonomy(
        tree,
        email="user@example.com",
        fetch_fn=boom,
    )
    assert updated["tree_hash"] == original_hash
    assert updated["newick"] == original_newick
    assert updated["n_leaves"] == tree["n_leaves"]
    layer = updated["taxonomy_layer"]
    assert layer["not_phylogenetic_tree"] is True
    assert layer["n_failed"] >= 1
    assert layer["tree_hash_unchanged"] == original_hash


def test_ncbi_taxonomy_attach_uses_real_payload_fields():
    tree = _tree_with_organisms()

    def fake_fetch(query, email, api_key=None):
        if query == "Homo sapiens":
            return {
                "taxon_id": "9606",
                "scientific_name": "Homo sapiens",
                "rank": "species",
                "lineage": [
                    {"taxon_id": "33208", "scientific_name": "Metazoa", "rank": "kingdom"},
                    {"taxon_id": "40674", "scientific_name": "Mammalia", "rank": "class"},
                    {"taxon_id": "9604", "scientific_name": "Hominidae", "rank": "family"},
                    {"taxon_id": "9605", "scientific_name": "Homo", "rank": "genus"},
                ],
                "common_name": "human",
                "source": "NCBI Taxonomy",
                "status": "RETRIEVED",
            }
        if query == "Mus musculus":
            return {
                "taxon_id": "10090",
                "scientific_name": "Mus musculus",
                "rank": "species",
                "lineage": [
                    {"taxon_id": "40674", "scientific_name": "Mammalia", "rank": "class"},
                    {"taxon_id": "10066", "scientific_name": "Muridae", "rank": "family"},
                    {"taxon_id": "10088", "scientific_name": "Mus", "rank": "genus"},
                ],
                "common_name": "house mouse",
                "source": "NCBI Taxonomy",
                "status": "RETRIEVED",
            }
        raise ncbi_fetch.NCBIQueryError("not found", "not_found")

    cache: dict = {}
    updated = taxonomy.attach_ncbi_taxonomy(
        tree,
        email="user@example.com",
        cache=cache,
        fetch_fn=fake_fetch,
    )
    human = next(leaf for leaf in updated["leaves"] if leaf["original_id"] == "human")
    assert human["taxonomy"]["taxon_id"] == "9606"
    assert human["taxonomy"]["scientific_name"] == "Homo sapiens"
    assert human["taxonomy"]["common_name"] == "human"
    assert human["organism_display"] == "Homo sapiens"
    user = next(leaf for leaf in updated["leaves"] if leaf["original_id"] == "user")
    assert user["organism_display"] == "unknown / user-provided"
    assert user["taxonomy"] is None
    groups = taxonomy.color_groups_by_rank(updated, "genus")
    assert groups["groups"][human["tree_id"]] == "Homo"
    assert "Unknown" in groups["legend"]
    assert "not phylogenetic distance" in groups["disclaimer"].lower()
    again = taxonomy.attach_ncbi_taxonomy(
        tree,
        email="user@example.com",
        cache=cache,
        fetch_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    human2 = next(leaf for leaf in again["leaves"] if leaf["original_id"] == "human")
    assert human2["taxonomy"]["taxon_id"] == "9606"


def test_taxonomy_cache_key_is_independent_of_tree_hash():
    key_a = taxonomy.cache_key("Homo sapiens")
    key_b = taxonomy.cache_key("Mus musculus")
    assert key_a != key_b
    tree_key = phylogeny.cache_key(
        alignment_hash="abc",
        method="neighbor_joining",
        distance_model="p_distance",
        rooting="unrooted",
        outgroup_tree_id="",
        bootstrap_replicates=0,
        bootstrap_seed=0,
        tool_version="1.85",
    )
    assert key_a != tree_key


def test_fetch_taxonomy_rejects_urls_and_empty_email():
    with pytest.raises(ValueError):
        ncbi_fetch.fetch_taxonomy("9606", "")
    with pytest.raises(ncbi_fetch.NCBIQueryError) as url:
        ncbi_fetch.fetch_taxonomy("https://evil.example/taxonomy", "user@example.com")
    assert url.value.category == "invalid_input"
    with pytest.raises(ncbi_fetch.NCBIQueryError) as empty:
        ncbi_fetch.fetch_taxonomy("  ", "user@example.com")
    assert empty.value.category == "invalid_input"


def test_parse_taxonomy_record_requires_taxid_and_name():
    parsed = ncbi_fetch._parse_taxonomy_record(
        {
            "TaxId": "9606",
            "ScientificName": "Homo sapiens",
            "Rank": "species",
            "LineageEx": [
                {"TaxId": "9605", "ScientificName": "Homo", "Rank": "genus"}
            ],
            "OtherNames": {"GenbankCommonName": "human"},
        }
    )
    assert parsed["taxon_id"] == "9606"
    assert parsed["common_name"] == "human"
    assert parsed["lineage"][0]["scientific_name"] == "Homo"
    with pytest.raises(ncbi_fetch.NCBIQueryError):
        ncbi_fetch._parse_taxonomy_record({"TaxId": "9606"})
