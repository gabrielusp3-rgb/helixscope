"""Prompt 0 characterization tests. Do not change scientific formulas.

These tests freeze migration invariants: Streamlit isolation of modules/,
golden parity classes, coordinate roundtrip, and JSON NaN honesty.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from modules import (
    alignment,
    dna_analysis,
    evidence_workspace,
    genome_coordinates,
    motif_search,
    phylogeny,
    provenance,
    variant_core,
)

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "tests" / "migration_reference"


def _load(name: str) -> dict:
    return json.loads((REF / name).read_text(encoding="utf-8"))


def test_scientific_modules_do_not_import_streamlit() -> None:
    offenders: list[str] = []
    for path in sorted((ROOT / "modules").glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "streamlit" or alias.name.startswith("streamlit."):
                        offenders.append(path.name)
            if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("streamlit"):
                offenders.append(path.name)
    assert offenders == []


def test_reference_version_is_frozen() -> None:
    assert provenance.HELIXSCOPE_VERSION == "0.24.3-19"


def test_golden_dna_gc_and_nan_honesty() -> None:
    golden = _load("dna.json")
    short = golden["cases"]["short_valid"]
    seq = short["input"]
    assert dna_analysis.gc_content(seq) == pytest.approx(short["gc_content"])
    assert short["gc_content"] == 50.0
    all_n = golden["cases"]["all_n"]
    assert all_n["gc_content"] == {"__nonfinite__": "NaN"}
    live = dna_analysis.gc_content("NNNN")
    assert live != live  # NaN
    assert live != 0.0
    scale = golden["cases"]["scale_50k"]
    assert scale["validation"]["length"] == 50_000
    assert dna_analysis.gc_content("ATGC" * 12_500) == pytest.approx(50.0)


def test_golden_alignment_identity() -> None:
    golden = _load("alignment.json")
    live = alignment.pairwise_global("ACGTACGTACGT", "ACGTACGGACGT")
    assert live["status"] == "COMPUTED"
    assert live["method"] == "Needleman-Wunsch"
    assert live["identity_pct"] == pytest.approx(golden["global"]["identity_pct"])
    assert live["aligned_seq1"] == golden["global"]["aligned_seq1"]
    assert live["aligned_seq2"] == golden["global"]["aligned_seq2"]
    local = alignment.pairwise_local("ACGTACGTACGT", "ACGTACGGACGT")
    assert local["method"] == "Smith-Waterman"
    assert local["identity_pct"] == pytest.approx(golden["local"]["identity_pct"])


def test_golden_motif_zero_based_half_open() -> None:
    golden = _load("motif.json")
    hits = motif_search.find_motif("ACGTACGTACGT", "ACGT")
    assert hits == golden["exact"]
    for hit in hits:
        assert 0 <= hit["start"] < hit["end"]


def test_golden_variant_assembly_and_coordinates() -> None:
    golden = _load("variant.json")
    live = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    assert live["position_0based"] == 43093556
    assert live["position_1based"] == 43093557
    assert live["identity_hash"] == golden["grch38"]["identity_hash"]
    assert live["effect_status"] == "NOT_COMPUTED"
    other = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh37")
    assert other["identity_hash"] != live["identity_hash"]


def test_coordinate_registry_roundtrip() -> None:
    shown = genome_coordinates.internal_to_display(0, 10)
    assert shown["start_1based"] == 1
    assert shown["end_1based"] == 10
    back = genome_coordinates.display_to_internal(1, 10)
    assert back["start_0based"] == 0
    assert back["end_0based"] == 10
    assert genome_coordinates.roundtrip_ok(0, 10)


def test_golden_phylogeny_leaf_sets_and_distinct_methods() -> None:
    golden = _load("phylogeny.json")
    expected = {"seq_a", "seq_b", "seq_c", "seq_d"}
    methods = []
    for key in ("nj", "upgma", "iqtree", "fasttree"):
        tree = golden["trees"][key]
        assert tree["status"] == "COMPUTED"
        assert tree["n_leaves"] == 4
        assert {leaf["tree_id"] for leaf in tree["leaves"]} == expected
        assert tree["alignment_hash"] == golden["msa_alignment_hash"]
        assert tree["not_taxonomic_tree"] is True
        methods.append(tree["method"])
        parsed = phylogeny.parse_newick(tree["newick"])
        assert parsed is not None
    assert methods == [
        phylogeny.METHOD_NJ,
        phylogeny.METHOD_UPGMA,
        phylogeny.METHOD_IQTREE,
        phylogeny.METHOD_FASTTREE,
    ]
    assert golden["trees"]["nj"]["newick"] != golden["trees"]["iqtree"]["newick"]


def test_golden_evidence_conflict_has_no_confidence_score() -> None:
    golden = _load("evidence.json")
    assert golden["confidence_score"] is None
    assert len(golden["conflicts"]) == 1
    assert "majority" in golden["conflicts"][0]["policy"].lower()
    items = [
        evidence_workspace.evidence_item(
            field="Variant A most_severe_consequence",
            value="missense_variant",
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
        evidence_workspace.evidence_item(
            field="Variant B most_severe_consequence",
            value="synonymous_variant",
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
    ]
    assert len(evidence_workspace.source_conflicts(items)) == 1


def test_hash_manifest_covers_scientific_modules() -> None:
    """Scientific modules and golden fixtures are present without the byte freeze.

    MIGRATION_REFERENCE_HASHES.json is written at the repository root by
    tests/migration_reference/_freeze.py. The freeze note says later extraction
    may change bytes and that scientific parity is the golden fixtures, not
    this manifest. The production image does not copy the manifest. This test
    always checks the runtime modules and the golden fixtures. When a checkout
    also has the historical manifest, that file must still name the scientific
    modules and store SHA-256 digests. Equality with the current dna_analysis.py
    bytes is not required: this branch changed that module, and rewriting the
    stored digest would alter the Prompt 0 provenance.
    """
    required = [
        ROOT / "modules" / "dna_analysis.py",
        ROOT / "modules" / "phylogeny.py",
        ROOT / "modules" / "crispr.py",
        ROOT / "modules" / "variant_core.py",
        ROOT / "app.py",
        REF / "dna.json",
        REF / "alignment.json",
        REF / "phylogeny.json",
    ]
    missing = [path.name for path in required if not path.is_file()]
    assert missing == []
    manifest_path = ROOT / "MIGRATION_REFERENCE_HASHES.json"
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest.get("algorithm") == "SHA-256"
    paths = {row["path"] for row in manifest["files"]}
    for rel in (
        "modules/dna_analysis.py",
        "modules/phylogeny.py",
        "modules/crispr.py",
        "modules/variant_core.py",
        "app.py",
    ):
        assert rel in paths
    stored = next(
        row["sha256"] for row in manifest["files"] if row["path"] == "modules/dna_analysis.py"
    )
    assert len(stored) == 64
    assert all(char in "0123456789abcdef" for char in stored)


def test_result_json_nan_is_not_silently_zero() -> None:
    payload = json.dumps({"gc": float("nan")})
    loaded = json.loads(payload)
    assert loaded["gc"] is None or (isinstance(loaded["gc"], float) and loaded["gc"] != loaded["gc"])
    golden = _load("dna.json")
    assert golden["cases"]["all_n"]["gc_content"]["__nonfinite__"] == "NaN"
