"""Golden scientific parity through helixscope_core. Same numbers as Prompt 0."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import pytest

from helixscope_core.alignment import pairwise_align
from helixscope_core.compare import parse_alignment_payload
from helixscope_core.coordinates import internal_to_display, roundtrip_ok
from helixscope_core.crispr import design_guides, model_availability
from helixscope_core.dna import analyze_dna, gc_content
from helixscope_core.errors import FoldingError, PhylogenyError
from helixscope_core.evidence import evidence_conflict_pack, evidence_item
from helixscope_core.explain import DIAGNOSIS_MARKERS, explain_result
from helixscope_core.motif import search_motifs
from helixscope_core.msa import build_collection_member, build_msa_result
from helixscope_core.phylogeny import (
    DISTANCE_P,
    METHOD_FASTTREE,
    METHOD_IQTREE,
    METHOD_NJ,
    METHOD_UPGMA,
    infer_tree,
    parse_newick,
)
from helixscope_core.protein import analyze_protein
from helixscope_core.rna import analyze_rna, codon_adaptation_index_report, fold_rna, transcribe_dna
from helixscope_core.structure import classify_mapping_status, mapping_contract
from helixscope_core.variants import identify_variant
from modules import provenance

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "tests" / "migration_reference"
FIXTURES = ROOT / "tests" / "fixtures"


def _load(name: str) -> dict:
    return json.loads((REF / name).read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def test_product_version_unchanged() -> None:
    assert provenance.HELIXSCOPE_VERSION == "0.24.3-19"


def test_dna_gc_and_nan_honesty() -> None:
    golden = _load("dna.json")
    short = golden["cases"]["short_valid"]
    seq = short["input"]
    assert gc_content(seq) == pytest.approx(short["gc_content"])
    live = gc_content("NNNN")
    assert math.isnan(live)
    assert live != 0.0
    bundled = analyze_dna(seq)
    assert bundled["status"] == "COMPUTED"
    assert bundled["gc_content"] == pytest.approx(50.0)
    assert bundled["sequence_hash"] == short["sequence_hash"]


def test_dna_50k_length_and_gc() -> None:
    seq = "ATGC" * 12_500
    assert gc_content(seq) == pytest.approx(50.0)
    info = analyze_dna(seq)
    assert info["length"] == 50_000


def test_rna_transcription_and_cai() -> None:
    golden = _load("rna.json")
    assert transcribe_dna("ATGC") == "AUGC"
    cai = codon_adaptation_index_report(golden["codon_fixture"]["coding_rna"])
    assert cai["status"] == "HEURISTIC"
    assert cai["value"] == pytest.approx(golden["codon_fixture"]["cai_report"]["value"], abs=1e-4)


def test_rna_viennarna_mfe_when_available() -> None:
    golden = _load("rna.json")["fold_python_viennarna"]
    try:
        folded = fold_rna(golden["sequence"], identifier="prompt0_hairpin")
    except FoldingError as exc:
        if getattr(exc, "category", "") == "TOOL_NOT_INSTALLED":
            pytest.skip("ViennaRNA not installed in this process")
        raise
    if folded.get("status") == "UNAVAILABLE":
        pytest.skip("ViennaRNA not installed in this process")
    assert folded["status"] == "PREDICTED"
    assert folded["backend"] == "viennarna_python"
    assert folded["dot_bracket"] == golden["dot_bracket"]
    assert folded["mfe_kcal_mol"] == pytest.approx(golden["mfe_kcal_mol"], abs=1e-4)
    rna = analyze_rna(golden["sequence"], fold=True)
    assert rna["status"] == "COMPUTED"
    assert rna["fold"]["status"] == "PREDICTED"


def test_protein_physicochemical_report() -> None:
    golden = _load("protein.json")
    seq = golden["sequence"]
    live = analyze_protein(seq)
    report = golden["physicochemical_report"]
    assert live["status"] == "COMPUTED"
    assert live["length"] == report["length"]
    assert live["molecular_weight_kda"] == pytest.approx(report["molecular_weight_kda"], abs=1e-6)
    assert live["isoelectric_point"] == pytest.approx(report["isoelectric_point"], abs=1e-2)
    assert live["gravy"] == pytest.approx(report["gravy"], abs=1e-3)
    assert live["ss_prediction"]["available"] is False


def test_pairwise_alignment_identity() -> None:
    golden = _load("alignment.json")
    global_aln = pairwise_align("ACGTACGTACGT", "ACGTACGGACGT", mode="global")
    assert global_aln["method"] == "Needleman-Wunsch"
    assert global_aln["identity_pct"] == pytest.approx(golden["global"]["identity_pct"])
    assert global_aln["aligned_seq1"] == golden["global"]["aligned_seq1"]
    local_aln = pairwise_align("ACGTACGTACGT", "ACGTACGGACGT", mode="local")
    assert local_aln["method"] == "Smith-Waterman"
    assert local_aln["identity_pct"] == pytest.approx(golden["local"]["identity_pct"])


def test_motif_zero_based_half_open() -> None:
    golden = _load("motif.json")
    hits = search_motifs("ACGTACGTACGT", "ACGT", molecule="DNA")
    assert hits == golden["exact"]
    for hit in hits:
        assert 0 <= hit["start"] < hit["end"]


def test_msa_prealigned_hash_and_phylogeny_leaf_sets() -> None:
    golden = _load("phylogeny.json")
    identifiers = ["seq_a", "seq_b", "seq_c", "seq_d"]
    aligned = [
        "ACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTTTTT",
        "TTTTACGTACGTACGTACGT",
        "ACGTACGTAAAAACGTACGT",
    ]
    members = [
        build_collection_member(
            sequence=row.replace("-", ""),
            identifier=ident,
            source="user input",
            molecule="DNA",
        )
        for ident, row in zip(identifiers, aligned)
    ]
    raw = "".join(f">{i}\n{row}\n" for i, row in zip(identifiers, aligned))
    msa_result = build_msa_result(
        raw_alignment=raw,
        members=members,
        tool="fixture",
        tool_version="prompt0",
        method="prealigned fixture",
        parameters={"source": "tests/phase13_streamlit_flow.PREALIGNED"},
    )
    assert msa_result["alignment_hash"] == golden["msa_alignment_hash"]
    expected = {"seq_a", "seq_b", "seq_c", "seq_d"}
    specs = [
        ("nj", METHOD_NJ, {"distance_model": DISTANCE_P}),
        ("upgma", METHOD_UPGMA, {"distance_model": DISTANCE_P}),
        ("iqtree", METHOD_IQTREE, {}),
        ("fasttree", METHOD_FASTTREE, {}),
    ]
    computed: list[str] = []
    for key, method, extra in specs:
        try:
            tree = infer_tree(msa_result, method=method, **extra)
        except PhylogenyError as exc:
            if method in (METHOD_IQTREE, METHOD_FASTTREE) and exc.category == "TOOL_NOT_INSTALLED":
                continue
            raise
        frozen = golden["trees"][key]
        assert tree["status"] == frozen["status"] == "COMPUTED"
        assert tree["method"] == frozen["method"]
        assert tree["n_leaves"] == 4
        assert {leaf["tree_id"] for leaf in tree["leaves"]} == expected
        assert tree["alignment_hash"] == golden["msa_alignment_hash"]
        assert tree["not_taxonomic_tree"] is True
        parsed = parse_newick(tree["newick"])
        assert parsed is not None
        computed.append(key)
    assert "nj" in computed
    assert "upgma" in computed


def test_crispr_pasted_spcas9_and_unavailable_models() -> None:
    golden = _load("crispr.json")
    guides = design_guides("ACGTACGTACGTACGTACGTAGG", "SpCas9")
    assert guides
    first = guides[0]
    scan = golden["guide_scan"][0]
    assert first["guide_sequence"] == scan["guide_sequence"]
    assert first["position"] == 1
    assert first["start_0based"] == 0
    assert first["end_0based"] == 20
    avail = model_availability()
    assert avail["on_target_azimuth"]["available"] is False
    assert avail["on_target_deephf"]["available"] is False
    assert avail["on_target_doench_ruleset2"]["available"] is False
    assert avail["cas9_guide_dna_3d_complex"]["available"] is False
    assert first.get("ruleset2_score") is None
    assert first.get("deephf_score") is None


def test_variant_identity_and_coordinates() -> None:
    golden = _load("variant.json")
    live = identify_variant("17 43093557 C G", assembly="GRCh38.p14")
    assert live["position_0based"] == 43093556
    assert live["position_1based"] == 43093557
    assert live["identity_hash"] == golden["grch38"]["identity_hash"]
    assert live["effect_status"] == "NOT_COMPUTED"
    other = identify_variant("17 43093557 C G", assembly="GRCh37")
    assert other["identity_hash"] != live["identity_hash"]
    shown = internal_to_display(0, 10)
    assert shown["start_1based"] == 1
    assert shown["end_1based"] == 10
    assert roundtrip_ok(0, 10)


def test_evidence_conflict_has_no_confidence_score() -> None:
    golden = _load("evidence.json")
    items = [
        evidence_item(
            field="Variant A most_severe_consequence",
            value="missense_variant",
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
        evidence_item(
            field="Variant B most_severe_consequence",
            value="synonymous_variant",
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
    ]
    pack = evidence_conflict_pack(items, kind="variants")
    assert pack["confidence_score"] is None
    assert golden["confidence_score"] is None
    assert pack["conflicts"]
    assert "majority" in pack["conflicts"][0]["policy"].lower()


def test_explainability_negative_claims() -> None:
    predicted = explain_result(
        "rna_fold",
        {
            "status": "PREDICTED",
            "mfe_kcal_mol": -1.2,
            "backend": "viennarna_python",
            "tool": "ViennaRNA RNAlib",
            "method": "RNA.fold MFE",
        },
    )
    blob = (predicted["plain_meaning"] + predicted["limitations"]).lower()
    assert predicted["status"] == "PREDICTED"
    assert "predicted" in blob
    assert "not an experimental" in blob or "not converted into rna 3d" in blob
    variant = explain_result(
        "variant",
        {
            "hgvs": "17:43093557 C>G",
            "assembly": "GRCh38",
            "most_severe_consequence": "missense_variant",
            "clinvar_significance": "Pathogenic",
        },
    )
    text = (variant["plain_meaning"] + variant["interpretation"] + variant["limitations"]).lower()
    for banned in DIAGNOSIS_MARKERS:
        assert banned not in text
    assert "not a helixscope diagnosis" in variant["interpretation"].lower()
    motif = explain_result("motif", {"pattern": "GAATTC", "n_hits": 2, "status": "COMPUTED"})
    assert "does not prove biological function" in motif["limitations"].lower()
    evo = explain_result("evolution", {"conservation_label": "HIGH CONSERVATION", "status": "COMPUTED"})
    assert "does not prove" in evo["interpretation"].lower()


def test_structure_fixture_hashes_and_4un3_uncertain() -> None:
    golden = _load("structures.json")
    for name, meta in golden["fixtures"].items():
        path = FIXTURES / name
        assert path.is_file()
        assert _sha256(path) == meta["sha256"]
    uncertain = classify_mapping_status(
        sequences_identical=True,
        polymer_index_mode=False,
        observed_coordinate_residues=0,
        polymer_length=0,
    )
    assert uncertain["mapping_status"] == "UNCERTAIN"
    contract = mapping_contract()
    assert contract["structure_mapping"] is None
    assert "4UN3" in contract["experimental_complex_pointers"]
    assert golden["4UN3_cache"] is None


def test_compare_8hsk_8hsf_fixture_parser() -> None:
    golden = _load("compare.json")
    payload = json.loads((FIXTURES / "rcsb_alignment_complete.json").read_text(encoding="utf-8"))
    record = parse_alignment_payload(payload)
    assert record["method"] == "fatcat-rigid"
    assert record["n_aligned_residue_pairs"] == golden["expected_semantics"]["n_aligned_ca_pairs"] == 20
    assert record["rmsd_block0_angstrom"] == pytest.approx(
        golden["expected_semantics"]["block_rmsd_angstrom"], abs=0.05
    )
    assert record["rmsd_global_angstrom"] == pytest.approx(
        golden["expected_semantics"]["global_rmsd_angstrom"], abs=0.05
    )
    tm = record["tm_scores"][0]["value"]
    assert tm == pytest.approx(golden["expected_semantics"]["tm_score"], abs=0.05)
    assert record["reference"]["entry_id"] == "8HSK"
    assert record["target"]["entry_id"] == "8HSF"
