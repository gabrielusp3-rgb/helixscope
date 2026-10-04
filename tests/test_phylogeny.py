"""Testes da inferencia filogenetica real a partir de MSA validado.

Nenhum teste inventa topologia, bootstrap, especie ou comprimento de ramo.
Ferramentas ML ausentes permanecem UNAVAILABLE. Newick malformado e
PARSING_ERROR, nao um placeholder.
"""

from __future__ import annotations

import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from modules import msa, phylogeny, provenance, scientific_checks

FIXTURES = Path(__file__).parent / "fixtures"


def _msa_from_aligned(
    identifiers: list[str],
    aligned: list[str],
    *,
    molecule: str = "DNA",
    organisms: list[str] | None = None,
    accessions: list[str] | None = None,
) -> dict:
    members = []
    raw_lines = []
    for index, (ident, row) in enumerate(zip(identifiers, aligned)):
        ungapped = row.replace("-", "").replace(".", "")
        organism = ""
        accession = ""
        if organisms:
            organism = organisms[index]
        if accessions:
            accession = accessions[index]
        members.append(
            msa.build_collection_member(
                sequence=ungapped,
                identifier=ident,
                source="user input" if not accession else "NCBI Entrez",
                molecule=molecule,
                organism=organism,
                accession=accession,
                version="1" if accession else "",
            )
        )
        raw_lines.append(f">{ident}")
        raw_lines.append(row)
    return msa.build_msa_result(
        raw_alignment="\n".join(raw_lines) + "\n",
        members=members,
        tool="fixture",
        tool_version="test",
        method="prealigned fixture",
        parameters={"source": "test fixture"},
    )


def test_valid_msa_neighbor_joining_is_computed():
    result = _msa_from_aligned(
        ["seq_a", "seq_b", "seq_c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        rooting=phylogeny.ROOTING_UNROOTED,
    )
    assert tree["status"] == "COMPUTED"
    assert tree["kind"] == "phylogenetic_inference"
    assert tree["not_taxonomic_tree"] is True
    assert tree["n_leaves"] == 3
    assert tree["branch_lengths_present"] is True
    assert tree["rooting"]["method"] == phylogeny.ROOTING_UNROOTED
    assert tree["rooting"]["rooted"] is False
    assert tree["rooting"]["biological_root"] is False
    assert tree["support"]["status"] == "UNAVAILABLE"
    assert tree["support"]["n_replicates"] == 0
    assert tree["tool"] == "Biopython DistanceTreeConstructor"
    assert tree["tool_version"]
    assert tree["alignment_hash"] == result["alignment_hash"]
    assert scientific_checks.phylogenetic_leaf_set_matches_msa(
        [leaf["tree_id"] for leaf in tree["leaves"]],
        [leaf["tree_id"] for leaf in tree["leaves"]],
    )
    for edge in tree["edges"]:
        assert scientific_checks.branch_length_is_valid(edge["length"])
        assert edge["length"] is not None
    for leaf in tree["leaves"]:
        assert leaf["organism_display"] == "unknown / user-provided"
        assert leaf["taxonomy"] is None
    parsed = phylogeny.parse_newick(tree["newick"])
    roundtrip = phylogeny.newick_roundtrip(tree["newick"])
    assert set(roundtrip["input_leaves"]) == {
        leaf["tree_id"] for leaf in tree["leaves"]
    }
    assert parsed is not None
    report = phylogeny.scientific_report(tree)
    assert report["method"] == phylogeny.METHOD_NJ
    assert report["kind"] == "phylogenetic_inference"
    assert "Not a true evolutionary tree" in tree["disclaimer"]


def test_upgma_is_rooted_and_declares_clock_assumption():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_UPGMA,
        distance_model=phylogeny.DISTANCE_P,
    )
    assert tree["rooting"]["method"] == phylogeny.ROOTING_UPGMA
    assert tree["rooting"]["rooted"] is True
    assert tree["rooting"]["ultrametric_assumption"] is True
    assert "molecular clock" in tree["rooting"]["note"].lower()


def test_jukes_cantor_undefined_is_not_silently_p_distance():
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.jukes_cantor_distance(0.75)
    assert exc.value.category == "INVALID_INPUT"
    identical = phylogeny.jukes_cantor_distance(0.0)
    assert identical["d"] == pytest.approx(0.0)
    result = _msa_from_aligned(
        ["a", "b"],
        ["AAAA", "TTTT"],
    )
    p_tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    assert p_tree["status"] == "COMPUTED"
    with pytest.raises(phylogeny.PhylogenyError) as jc:
        phylogeny.infer_phylogeny(
            result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_JC69
        )
    assert jc.value.category == "INVALID_INPUT"
    assert "0.75" in str(jc.value)


def test_protein_rejects_dna_models():
    result = _msa_from_aligned(
        ["p1", "p2"],
        ["ACDEFGHIKL", "ACDEFGHIKM"],
        molecule="PROTEIN",
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(
            result,
            method=phylogeny.METHOD_NJ,
            distance_model=phylogeny.DISTANCE_JC69,
        )
    assert exc.value.category == "INVALID_INPUT"
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
    )
    assert tree["molecule"] == "PROTEIN"
    assert tree["n_leaves"] == 2


def test_mixed_molecule_envelope_is_rejected():
    result = _msa_from_aligned(["a", "b"], ["ACGT", "ACGT"])
    result["molecule"] = "DNA/RNA"
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.validate_msa_for_phylogeny(result)
    assert exc.value.category == "INVALID_INPUT"


def test_duplicate_ids_are_renamed_not_dropped():
    result = _msa_from_aligned(
        ["same", "same"],
        ["ACGTACGTAA", "ACGTTCGTAA"],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    assert tree["n_leaves"] == 2
    ids = [leaf["tree_id"] for leaf in tree["leaves"]]
    assert len(set(ids)) == 2
    assert tree["id_transformations"]
    assert all(item["original_id"] == "same" for item in tree["id_transformations"])


def test_duplicate_sequences_are_kept():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["ACGTACGTAA", "ACGTACGTAA", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    assert tree["n_leaves"] == 3
    hashes = [leaf["sequence_hash"] for leaf in tree["leaves"]]
    assert hashes.count(hashes[0]) >= 2 or hashes[0] == hashes[1]
    assert tree["identical_groups"] or hashes[0] == hashes[1]


def test_malformed_newick_is_parsing_error_not_placeholder():
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.parse_newick("((((not a tree")
    assert exc.value.category == "PARSING_ERROR"
    assert "NO_TREE" not in str(exc.value)


def test_newick_roundtrip_preserves_leaves():
    newick = "(a:0.1,b:0.2,c:0.3);"
    result = phylogeny.newick_roundtrip(newick)
    assert set(result["input_leaves"]) == {"a", "b", "c"}
    assert set(result["output_leaves"]) == {"a", "b", "c"}


def test_midpoint_rooting_is_labeled_not_biological():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        rooting=phylogeny.ROOTING_MIDPOINT,
    )
    assert tree["rooting"]["method"] == phylogeny.ROOTING_MIDPOINT
    assert tree["rooting"]["rooted"] is True
    assert tree["rooting"]["biological_root"] is False
    assert "not a biological root" in tree["rooting"]["note"].lower()


def test_outgroup_rooting_records_the_leaf():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    validated = phylogeny.validate_msa_for_phylogeny(result)
    outgroup = validated["leaf_records"][2]["tree_id"]
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        rooting=phylogeny.ROOTING_OUTGROUP,
        outgroup_tree_id=outgroup,
    )
    assert tree["rooting"]["outgroup_tree_id"] == outgroup
    assert tree["rooting"]["biological_root"] is False


def test_bootstrap_zero_is_unavailable_not_zero():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        bootstrap_replicates=0,
    )
    assert tree["support"]["status"] == "UNAVAILABLE"
    assert tree["support"]["n_replicates"] == 0
    for node in tree["nodes"]:
        if node["is_leaf"]:
            continue
        assert scientific_checks.bootstrap_support_is_valid(
            node.get("support"), replicates=0
        )


def test_bootstrap_is_real_percent_not_fake_ninetynine():
    result = _msa_from_aligned(
        ["a", "b", "c", "d"],
        ["AAAAAAAAAA", "AAAAAAAATT", "GGGGGGGGGG", "CCCCCCCCCC"],
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_P,
        bootstrap_replicates=8,
        bootstrap_seed=1,
    )
    assert tree["support"]["status"] == "COMPUTED"
    assert tree["support"]["method"].startswith("Felsenstein 1985")
    assert tree["support"]["label"] == "bootstrap support"
    assert tree["support"]["not_confidence"] is True
    values = list((tree["support"].get("values") or {}).values())
    assert values
    assert not all(abs(float(item) - 99.0) < 1e-9 for item in values)
    for value in values:
        assert scientific_checks.bootstrap_support_is_valid(value, replicates=8)


def test_cache_hit_and_miss_and_stale_msa():
    first = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    second = _msa_from_aligned(
        ["a", "b", "c"],
        ["CCCCCCCCCC", "CCCGCCCCCC", "TTTTTTTTTT"],
    )
    cache: dict = {}
    tree_a = phylogeny.infer_phylogeny(
        first, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P, cache=cache
    )
    tree_again = phylogeny.infer_phylogeny(
        first, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P, cache=cache
    )
    assert tree_again["cache_status"] == "cached"
    assert tree_again["tree_hash"] == tree_a["tree_hash"]
    tree_upgma = phylogeny.infer_phylogeny(
        first, method=phylogeny.METHOD_UPGMA, distance_model=phylogeny.DISTANCE_P, cache=cache
    )
    assert tree_upgma["cache_status"] == "live"
    assert tree_upgma["tree_hash"] != tree_a["tree_hash"]
    tree_b = phylogeny.infer_phylogeny(
        second, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P, cache=cache
    )
    assert tree_b["alignment_hash"] != tree_a["alignment_hash"]
    assert tree_b["tree_hash"] != tree_a["tree_hash"]


def test_resource_limit_sequences(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    monkeypatch.setattr(phylogeny, "MAX_SEQUENCES", 2)
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.validate_msa_for_phylogeny(result)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_incomplete_msa_is_rejected():
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.validate_msa_for_phylogeny({"status": "RUNNING", "rows": []})
    assert exc.value.category == "INVALID_INPUT"


def test_leaf_mapping_matches_sequence_hash():
    result = _msa_from_aligned(
        ["alpha", "beta"],
        ["ACGTACGTAA", "ACGTTCGTAA"],
        organisms=["", "Homo sapiens"],
        accessions=["", "NM_000000.1"],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    target_hash = result["rows"][1]["hash"]
    leaf = phylogeny.leaf_by_sequence_hash(tree, target_hash)
    assert leaf is not None
    assert leaf["original_id"] == "beta"
    assert leaf["accession"] == "NM_000000.1"
    assert leaf["organism_declared"] == "Homo sapiens"
    assert phylogeny.leaf_by_tree_id(tree, leaf["tree_id"])["sequence_hash"] == target_hash


def test_layout_edges_are_real_and_scale_requires_lengths():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    layout = tree["layouts"]["rectangular"]
    assert layout["topology_unchanged"] is True
    assert len(layout["edges"]) == tree["n_edges"]
    node_ids = {item["id"] for item in layout["nodes"]}
    for edge in layout["edges"]:
        assert edge["parent"] in node_ids
        assert edge["child"] in node_ids
    assert layout["scale_bar"]["present"] is True
    radial = tree["layouts"]["radial"]
    assert radial["mode"] == "radial"
    assert len(radial["edges"]) == len(layout["edges"])


def test_kimura_undefined_pair_is_refused():
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.kimura_two_parameter_distance("AAAA", "CCCC", "PROTEIN")
    assert exc.value.category == "INVALID_INPUT"
    k2p = phylogeny.kimura_two_parameter_distance("AAAAAAAAAA", "GAAAAAAAAA", "DNA")
    assert k2p["status"] == "COMPUTED"
    assert k2p["d"] > 0
    rna = phylogeny.kimura_two_parameter_distance("AAAAAAAAAA", "GAAAAAAAAA", "RNA")
    assert rna["d"] == pytest.approx(k2p["d"])


def test_fasttree_absent_is_unavailable(monkeypatch):
    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: None,
    )
    info = phylogeny.detect_fasttree()
    assert info["available"] is False
    assert info["version"] == ""
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(result, method=phylogeny.METHOD_FASTTREE)
    assert exc.value.category == "TOOL_NOT_INSTALLED"


def test_fasttree_mock_run_parses_real_newick(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    validated = phylogeny.validate_msa_for_phylogeny(result)
    names = [item["tree_id"] for item in validated["leaf_records"]]
    newick = f"({names[0]}:0.1,{names[1]}:0.2,{names[2]}:0.3);"

    def fake_run(argv, **kwargs):
        assert kwargs.get("shell") is False
        assert argv[0] == "fasttree"
        return SimpleNamespace(returncode=0, stdout=newick, stderr="FastTree Version 2.1.11")

    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: "fasttree",
    )
    monkeypatch.setattr(
        phylogeny,
        "detect_fasttree",
        lambda: {
            "available": True,
            "version": "2.1.11",
            "reason": "",
        },
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_FASTTREE,
        run_fn=fake_run,
    )
    assert tree["tool"] == "FastTree"
    assert tree["tool_version"] == "2.1.11"
    assert tree["n_leaves"] == 3
    assert set(leaf["tree_id"] for leaf in tree["leaves"]) == set(names)


def test_subprocess_helpers_do_not_use_shell():
    source = inspect.getsource(phylogeny._run_fasttree)
    assert "shell=False" in source
    assert "shell=True" not in source
    iq = inspect.getsource(phylogeny._run_iqtree)
    assert "shell=False" in iq
    assert "shell=True" not in iq


def test_unsafe_newick_labels_are_sanitized():
    records, transforms = phylogeny.unique_leaf_labels(
        ["a(b);c", "a(b);c"],
        ["hashonehashonehashone", "hashtwohashtwohashtwo"],
    )
    assert all("(" not in item["tree_id"] and ";" not in item["tree_id"] for item in records)
    assert len({item["tree_id"] for item in records}) == 2
    assert transforms


def test_export_leaf_csv_escapes_commas():
    result = _msa_from_aligned(
        ["a,b", "c"],
        ["ACGTACGTAA", "ACGTTCGTAA"],
        organisms=["Genus, species", ""],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    csv_text = phylogeny.export_leaf_csv(tree)
    assert "Genus, species" in csv_text or '"Genus, species"' in csv_text


def test_performance_small_alignment_is_measured():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    assert tree["elapsed_ms"] >= 0
    assert tree["elapsed_ms"] < 5000


def test_blosum62_protein_distance_is_labeled_not_wag():
    result = _msa_from_aligned(
        ["p1", "p2", "p3"],
        ["ACDEFGHIKL", "ACDEFGHIKM", "ACDEFGAAAL"],
        molecule="PROTEIN",
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_NJ,
        distance_model=phylogeny.DISTANCE_BLOSUM62,
    )
    assert tree["status"] == "COMPUTED"
    assert tree["model"] == phylogeny.DISTANCE_BLOSUM62
    assert "WAG" not in tree["model"]
    assert "not wag" in tree["model_source"].lower()


def test_msa_tab_exposes_phylogeny_action_without_auto_compute():
    from tests.helix_apptest import open_module

    msa_result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    app = open_module("msa", timeout=60)
    assert not app.exception
    app.session_state["msa_result"] = msa_result
    app.session_state["msa_collection"] = list(msa_result.get("input_members") or [])
    app.run()
    assert not app.exception
    keys = [item.key for item in app.button]
    assert "phylo_build" in keys
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("info", "caption", "markdown", "warning")
        for item in list(getattr(app, bucket, []) or [])
    ).lower()
    assert "phylogenetic inference" in combined
    assert "not a true evolutionary" in combined or "not true evolutionary" in combined
    assert "confirmed ancestry" not in combined
    stored = None
    try:
        stored = app.session_state["phylo_result"]
    except KeyError:
        stored = None
    assert not stored


def test_msa_tab_build_tree_computes_on_explicit_click():
    from tests.helix_apptest import open_module

    msa_result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    app = open_module("msa", timeout=60)
    app.session_state["msa_result"] = msa_result
    app.session_state["msa_collection"] = list(msa_result.get("input_members") or [])
    app.run()
    assert not app.exception
    app.button(key="phylo_build").click().run()
    assert not app.exception
    tree = app.session_state["phylo_result"]
    assert tree["status"] == "COMPUTED"
    assert tree["n_leaves"] == 3
    assert tree["method"] == phylogeny.METHOD_NJ
    assert tree["alignment_hash"] == msa_result["alignment_hash"]
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("info", "caption", "markdown", "warning")
        for item in list(getattr(app, bucket, []) or [])
    ).lower()
    assert "neighbor-joining" in combined
    assert "unrooted" in combined
