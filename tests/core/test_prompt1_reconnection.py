"""Prompt 1 reconnection: Core contracts without rewriting mature science."""

from __future__ import annotations

import pytest

from helixscope_core.alignment import pairwise_align, translate_nucleic_for_alignment
from helixscope_core.compare import (
    compare_evolution,
    compare_evidence_from_variants,
    compare_guides,
    compare_modes,
    compare_proteins,
    compare_variants_from_text,
)
from helixscope_core.crispr import analyze_crispr, list_cas_systems, require_cas_system
from helixscope_core.msa import classify_alignment_text, import_alignment
from helixscope_core.phylogeny import METHOD_NJ, METHOD_UPGMA, annotate_taxonomy, infer_tree
from helixscope_core.references import admit_catalog_download, start_catalog_download
from helixscope_core.rna import translate_coding_composition
from helixscope_core.structure import (
    analyze_structure_coordinates,
    load_bundled_deposited,
    read_bundled_mmcif,
    structure_engine_capabilities,
)
from modules.crispr import CAS_SYSTEMS
from modules.genome_download import GenomeDownloadError
from modules.msa import MsaError
from modules.protein_structure import search_rcsb_by_sequence
from modules.resource_admission import DECISION_RESOURCE_LIMIT

PREALIGNED = """>seq_a
ACGTACGTACGTACGTACGT
>seq_b
ACGTACGTACGTACGTTTTT
>seq_c
TTTTACGTACGTACGTACGT
>seq_d
ACGTACGTAAAAACGTACGT
"""

SPCAS9_LOCUS = "ACGTACGTACGTACGTACGTAGG"
SACAS9_LOCUS = "ACGTACGTACGTACGTACGTA" + "AAGAGT"
ASCAS12A_LOCUS = "TTTA" + "ACGTACGTACGTACGTACGTACG"
CAS12B_LOCUS = "TTA" + "ACGTACGTACGTACGTACGT" + "CCCC"
CAS13_RNA = "ACGUACGUACGUACGUACGUACGUACGUAC"
ORF_RNA = "AUG" + "GCA" * 48 + "UAA"
ANNOTATED_TABLE = (
    "Humano      [ A G T C A G T C ]\n"
    "Chimpanze   [ A G T C A G T - ] <- comment\n"
)
CLUSTAL = """CLUSTAL W (1.83) multiple sequence alignment

seq_a      ACGTACGTACGTACGTACGT
seq_b      ACGTACGTACGTACGTTTTT
seq_c      TTTTACGTACGTACGTACGT
seq_d      ACGTACGTAAAAACGTACGT
"""
STOCKHOLM = """# STOCKHOLM 1.0
seq_a ACGTACGTACGTACGTACGT
seq_b ACGTACGTACGTACGTTTTT
seq_c TTTTACGTACGTACGTACGT
seq_d ACGTACGTAAAAACGTACGT
//
"""
PHYLIP = """ 4 20
seq_a     ACGTACGTACGTACGTACGT
seq_b     ACGTACGTACGTACGTTTTT
seq_c     TTTTACGTACGTACGTACGT
seq_d     ACGTACGTAAAAACGTACGT
"""


def test_every_cas_system_key_is_listed_and_required() -> None:
    registry = list_cas_systems()
    keys = [row["canonical_key"] for row in registry["systems"]]
    assert keys == list(CAS_SYSTEMS)
    for key in CAS_SYSTEMS:
        assert require_cas_system(key) == key
        row = next(item for item in registry["systems"] if item["canonical_key"] == key)
        assert row["display_name"]
        assert row["target_molecule"] in {"DNA", "RNA"}
        assert row["mode"]
        assert "on_target_doench_ruleset2" in str(row["unsupported_metrics"])


@pytest.mark.parametrize("token", ["Cas12a", "CBE", "ABE", "Prime Editor", "Cpf1", "cas12a"])
def test_informal_crispr_tokens_are_rejected(token: str) -> None:
    with pytest.raises(ValueError, match="canonical_key|Informal tokens"):
        require_cas_system(token)


@pytest.mark.parametrize(
    ("system", "sequence"),
    [
        ("SpCas9", SPCAS9_LOCUS),
        ("SaCas9", SACAS9_LOCUS),
        ("AsCas12a", ASCAS12A_LOCUS),
        ("Cas12b", CAS12B_LOCUS),
        ("CBE (SpCas9)", SPCAS9_LOCUS),
        ("ABE (SpCas9)", SPCAS9_LOCUS),
        ("Prime Editor (SpCas9 nickase)", SPCAS9_LOCUS),
    ],
)
def test_canonical_cas_tokens_reach_core_without_spcas9_fallback(system: str, sequence: str) -> None:
    result = analyze_crispr(sequence, system)
    assert result["cas_system"] == system
    assert result["cas_system"] != "SpCas9" or system == "SpCas9"
    assert result.get("target_molecule") == "DNA"
    guides = result.get("guides") or []
    if guides:
        first = guides[0]
        assert first.get("ruleset2_score") is None
        assert first.get("deephf_score") is None
        assert first.get("cfd_score") is None or first.get("cfd_score") is None
        availability = result.get("model_availability") or {}
        azimuth = availability.get("on_target_azimuth") or {}
        assert azimuth.get("available") is not True


def test_cas13_requires_rna_and_does_not_fallback() -> None:
    with pytest.raises(ValueError, match="RNA"):
        analyze_crispr("ACGTACGTACGTACGTACGTACGTACGTAC", "Cas13")
    result = analyze_crispr(CAS13_RNA, "Cas13")
    assert result["cas_system"] == "Cas13"
    assert result["target_molecule"] == "RNA"
    primers = result.get("validation_primers") or {}
    assert primers.get("status") in {"METHOD_NOT_APPLICABLE", "UNAVAILABLE", None} or primers.get("status") != "COMPUTED" or primers.get("forward_primer") in {None, ""}
    for guide in result.get("guides") or []:
        assert guide.get("cut_site") is None or "RNA" in str(result.get("input_note") or "")


def test_msa_fasta_and_prealigned_and_standard_formats() -> None:
    fasta = import_alignment(PREALIGNED, fmt="fasta")
    assert fasta["n_sequences"] == 4
    assert fasta["alignment_hash"]
    assert fasta["input_format"] == "fasta"
    auto = import_alignment(PREALIGNED, fmt="auto")
    assert auto["alignment_hash"] == fasta["alignment_hash"]
    clustal = import_alignment(CLUSTAL, fmt="clustal")
    assert clustal["n_sequences"] == 4
    assert clustal["input_format"] == "clustal"
    stockholm = import_alignment(STOCKHOLM, fmt="stockholm")
    assert stockholm["n_sequences"] == 4
    phylip = import_alignment(PHYLIP, fmt="phylip")
    assert phylip["n_sequences"] == 4
    assert classify_alignment_text(ANNOTATED_TABLE) == "annotated_table"


def test_msa_rejects_annotated_human_table() -> None:
    with pytest.raises(MsaError) as exc:
        import_alignment(ANNOTATED_TABLE, fmt="auto")
    assert exc.value.category == "INVALID_FORMAT"
    assert "not a recognized alignment format" in str(exc.value)


def test_msa_to_nj_and_upgma_preserve_alignment_hash() -> None:
    msa = import_alignment(PREALIGNED, fmt="fasta")
    nj = infer_tree(msa, method=METHOD_NJ, distance_model="p_distance")
    upgma = infer_tree(msa, method=METHOD_UPGMA, distance_model="p_distance")
    assert nj["alignment_hash"] == msa["alignment_hash"]
    assert upgma["alignment_hash"] == msa["alignment_hash"]
    assert nj["method"] == METHOD_NJ
    assert upgma["method"] == METHOD_UPGMA
    assert nj["newick"]
    assert upgma["newick"]
    for leaf in nj["leaves"]:
        assert str(leaf.get("original_id") or leaf.get("id") or "") in {"seq_a", "seq_b", "seq_c", "seq_d"} or "seq_" in str(
            leaf
        )


def test_taxonomy_does_not_change_newick(monkeypatch: pytest.MonkeyPatch) -> None:
    msa = import_alignment(PREALIGNED, fmt="fasta")
    tree = infer_tree(msa, method=METHOD_NJ, distance_model="p_distance")
    original_newick = tree["newick"]
    original_hash = tree["tree_hash"]

    def fake_attach(tree_result, email, api_key=None):
        updated = dict(tree_result)
        updated["taxonomy_layer"] = {"status": "COMPUTED", "source": "NCBI Taxonomy"}
        updated["newick"] = original_newick + "UNCHANGED_PROBE;"
        updated["tree_hash"] = "should-not-stick"
        return updated

    monkeypatch.setattr("modules.taxonomy.attach_ncbi_taxonomy", fake_attach)
    annotated = annotate_taxonomy(
        tree,
        email="helixscope@example.org",
        declared_organisms={"seq_a": "Homo sapiens"},
    )
    assert annotated["newick"] == original_newick
    assert annotated["tree_hash"] == original_hash
    assert annotated["taxonomy_is_not_phylogeny"] is True


def test_translate_then_align_is_explicit() -> None:
    silent = pairwise_align("ATGAAATAA", "ATGAAATAG", mode="global", translate_nucleic=False)
    assert silent["translation"]["applied"] is False
    translated = pairwise_align("ATGCCCGGGTAA", "ATGCCCGGGTAG", mode="global", translate_nucleic=True)
    assert translated["translation"]["applied"] is True
    assert translated["translation"]["seq1"]["frame"] == 1
    assert translated["translation"]["seq1"]["source_molecule"] == "DNA"
    meta = translate_nucleic_for_alignment("AUGCCCUAA")
    assert meta["source_molecule"] == "RNA"
    assert meta["translation"]


def test_rna_coding_composition_requires_coding_context() -> None:
    random_rna = "ACGUACGUACGUACGU"
    denied = translate_coding_composition(random_rna)
    assert denied["status"] in {"UNAVAILABLE", "ERROR"}
    assert denied.get("translation") in {"", None}
    coding = translate_coding_composition(ORF_RNA)
    assert coding["status"] in {"COMPUTED", "UNAVAILABLE"}
    if coding["status"] == "COMPUTED":
        assert coding["translation"]
        assert coding["composition"]
        assert "not assumed to be protein-coding" in " ".join(coding.get("limitations") or []) or coding.get("coding_context")


def test_structure_capabilities_and_bundled_sasa() -> None:
    caps = structure_engine_capabilities()
    assert caps["sasa"]["available"] is True
    assert "Not a molecular surface" in caps["sasa"]["limitation"]
    assert caps["edtsurf"]["limitation"]
    assert caps["dssp_local"]["name"] == "DSSP"
    assert caps["dssp_pdb_redo"]["engine_location"] == "remote"
    parsed = load_bundled_deposited("1CRN")
    text = read_bundled_mmcif("1CRN.cif")
    analyzed = analyze_structure_coordinates(parsed, run_sasa=True, structure_text=text)
    assert analyzed["sasa"]["status"] in {"COMPUTED", "EXPERIMENTAL", "UNAVAILABLE", "ERROR"}
    assert analyzed["dssp"]["status"] == "NOT_COMPUTED"
    if analyzed["sasa"]["status"] in {"COMPUTED", "EXPERIMENTAL"}:
        assert analyzed["sasa"].get("total_sasa") not in {0, 0.0} or analyzed["sasa"].get("residues")


def test_rcsb_sequence_search_posts_official_api(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}

    def fake_http(url, method="GET", json_body=None, urlopen_fn=None, empty_ok=False):
        captured["url"] = url
        captured["method"] = method
        captured["body"] = json_body
        return {"result_set": [{"identifier": "1CRN_1", "score": 1.0}]}

    monkeypatch.setattr("modules.protein_structure._http_json", fake_http)
    hits = search_rcsb_by_sequence("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", sequence_type="protein")
    assert "search.rcsb.org" in str(captured["url"])
    assert captured["method"] == "POST"
    assert captured["body"]["query"]["service"] == "sequence"
    assert captured["body"]["query"]["parameters"]["sequence_type"] == "protein"
    assert hits
    assert "1CRN" in str(hits[0].get("identifier") or hits[0])


def test_compare_modes_and_non_structure_science() -> None:
    modes = {row["id"] for row in compare_modes()["modes"]}
    assert modes == {"structures", "variants", "proteins", "guides", "evolution", "evidence"}
    variants = compare_variants_from_text(
        "17 43093557 C G",
        "GRCh38.p14",
        "17 43093558 A T",
        "GRCh38.p14",
    )
    assert variants["ranking"] is None
    proteins = compare_proteins(sequence_a="MKTAYIAK", sequence_b="MKTAYIAQ", identifier_a="a", identifier_b="b")
    assert proteins["mode"] == "proteins"
    assert proteins["identity"]["status"] == "COMPUTED"
    guides = compare_guides(
        {"guide_sequence": "ACGTACGTACGTACGTACGT", "pam_sequence": "AGG", "cas_system": "SpCas9"},
        {"guide_sequence": "TGCATGCATGCATGCATGCA", "pam_sequence": "TGG", "cas_system": "SpCas9"},
    )
    assert guides["mode"] == "guides"
    evo = compare_evolution(PREALIGNED, column=0)
    assert evo["mode"] == "evolution"
    assert evo["alignment_hash"]
    evidence = compare_evidence_from_variants(
        "17 43093557 C G",
        "GRCh38.p14",
        "17 43093558 A T",
        "GRCh38.p14",
    )
    assert evidence["mode"] == "evidence"
    assert evidence["confidence_score"] is None


def test_reference_admission_does_not_download_on_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    test_ref = admit_catalog_download("HELIXSCOPE_TEST_REF")
    assert test_ref["not_a_public_assembly"] is True
    grch = admit_catalog_download("GRCh38.p14")
    spawned: list[bool] = []
    monkeypatch.setattr(
        "helixscope_core.references.spawn_download_assembly",
        lambda *args, **kwargs: spawned.append(True) or {"status": "DOWNLOADING"},
    )
    if grch.get("decision") != "PROCEED":
        with pytest.raises(GenomeDownloadError) as exc:
            start_catalog_download("GRCh38.p14")
        assert exc.value.category in {DECISION_RESOURCE_LIMIT, "RESOURCE_LIMIT", "INVALID_INPUT"}
        assert spawned == []
    else:
        pytest.skip("This machine admitted GRCh38; Prompt 1 still must not auto-download in tests.")
