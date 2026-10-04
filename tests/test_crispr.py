"""Testes unitarios do modulo modules.crispr.

Cobrem a deteccao de guias por PAM NGG (posicao e fita), o intervalo da
pontuacao heuristica, o comportamento da varredura de off-targets (alvo perfeito
excluido e sitios divergentes ausentes), a faixa de Tm dos primers de validacao,
o registro de sistemas Cas, a declaracao de indisponibilidade de modelos reais,
as metricas por guia (poli-T, auto-complementaridade, especificidade local,
risco), o ranking combinado, a janela de editores de base e as notas de prime
editing. Nenhum teste faz chamadas de rede.
"""

import pytest

from modules import crispr
from tests.helix_apptest import open_module


@pytest.fixture
def guide() -> str:
    """Protoespacador de 20 nt com GC equilibrado, sem homopolimeros."""
    return "ACGTACGTACGTACGTACGT"


@pytest.fixture
def target_with_pam(guide) -> str:
    """Sequencia alvo com o guia seguido de um PAM NGG (AGG)."""
    return guide + "AGG"


@pytest.fixture
def flanking_sequence() -> str:
    """Sequencia flanqueadora balanceada (GC 50%) que contem o guia."""
    return "ACGT" * 150


def test_find_guides_position_and_strand(guide, target_with_pam):
    guides = crispr.find_guides(target_with_pam, "SpCas9")
    match = [
        g
        for g in guides
        if g["guide_sequence"] == guide and g["strand"] == "+"
    ]
    assert match, "esperado ao menos um guia na fita sense"
    hit = match[0]
    assert hit["position"] == 1
    assert hit["pam_sequence"] == "AGG"


def test_score_guide_doench_within_unit_interval(guide):
    doc = " ".join((crispr.score_guide_doench.__doc__ or "").split())
    assert "Implementacao heuristica simplificada inspirada em" in doc
    assert "nao reproduz o modelo/algoritmo original publicado" in doc
    for candidate in [guide, "AAAAAAAAAAAAAAAAAAAA", "GGGGCCCCGGGGCCCCGGGG", ""]:
        score = crispr.score_guide_doench(candidate)
        assert 0.0 <= score <= 1.0


def test_find_off_targets_identical_sequence_excludes_on_target(guide):
    assert crispr.find_off_targets(guide, guide, max_mismatches=3) == []


def test_find_off_targets_distant_sequence_absent(guide):
    distant = "TGCATGCATGCATGCATGCA"
    assert crispr.find_off_targets(guide, distant, max_mismatches=1) == []


def test_design_primer_pair_tm_in_range(guide, flanking_sequence):
    primers = crispr.design_primer_pair(guide, flanking_sequence)
    assert 55.0 <= primers["forward_tm"] <= 65.0
    assert 55.0 <= primers["reverse_tm"] <= 65.0


def test_list_supported_systems_includes_editors():
    systems = crispr.list_supported_systems()
    for expected in ["SpCas9", "SaCas9", "AsCas12a", "Cas12b", "Cas13"]:
        assert expected in systems
    assert any("CBE" in name for name in systems)
    assert any("ABE" in name for name in systems)
    assert any("Prime Editor" in name for name in systems)
    assert crispr.list_interface_systems() == ["SpCas9"]


def test_system_info_rejects_unknown_system():
    with pytest.raises(ValueError):
        crispr.system_info("NotACas")


def test_system_info_cas13_targets_rna():
    info = crispr.system_info("Cas13")
    assert info["target_molecule"] == "RNA"
    assert info["pam_system"] is None


def test_model_availability_keeps_published_models_unavailable():
    availability = crispr.model_availability()
    for name in (
        "on_target_doench_ruleset2",
        "on_target_azimuth",
        "on_target_deephf",
        "specificity_mit",
        "indel_frameshift_profile",
        "cas9_guide_dna_3d_complex",
        "distributed_workers",
    ):
        assert availability[name]["available"] is False
        assert availability[name]["requires"]
    gw = availability["genome_wide_off_targets"]
    assert gw["requires"]
    if gw["available"]:
        assert "COMPLETED_FULL_REFERENCE" in gw["reason"]
    else:
        assert "genome-wide" in gw["reason"].lower()
    assert availability["hsu_2013_single_hit"]["available"] is True
    assert availability["provided_reference_offtarget_search"]["available"] is True
    assert availability["specificity_cfd"]["available"] is True
    assert availability["specificity_cfd_guide_aggregate"]["available"] is True
    assert availability["specificity_mit_aggregate_over_search"]["available"] is True
    assert "Rule Set 2" in availability["on_target_doench_ruleset2"]["reason"]
    assert "silent" in availability["on_target_doench_ruleset2"]["reason"].lower() or (
        "not used as a silent" in availability["on_target_doench_ruleset2"]["reason"]
    )
    assert "pickle" not in availability["specificity_cfd"]["reason"].lower() or (
        "no pickle" in availability["specificity_cfd"]["reason"].lower()
    )


def test_find_guides_cas12b_uses_upstream_pam():
    target = "TTA" + "ACGTACGTACGTACGTACGT" + "CCCC"
    guides = crispr.find_guides(target, "Cas12b")
    assert any(g["guide_sequence"] == "ACGTACGTACGTACGTACGT" for g in guides)


def test_find_guides_cas13_has_no_pam_and_spacer():
    rna_target = "ACGUACGUACGUACGUACGUACGUACGUAC"
    guides = crispr.find_guides(rna_target, "Cas13")
    assert guides
    assert all(g["pam_sequence"] == "" for g in guides)
    assert all("U" in g["guide_sequence"] or "A" in g["guide_sequence"] for g in guides)
    assert guides[0]["pfs"] == "A"


def test_find_guides_cas13_rejects_dna():
    with pytest.raises(ValueError, match="RNA"):
        crispr.find_guides("ACGTACGTACGTACGTACGTACGTACGTAC", "Cas13")


def test_find_guides_spcas9_rejects_rna():
    with pytest.raises(ValueError, match="DNA"):
        crispr.find_guides("ACGUACGUACGUACGUACGUAGG", "SpCas9")
    with pytest.raises(ValueError):
        crispr.find_guides("ACGTACGT", "Unknown")


def test_poly_t_signal_detects_run():
    assert crispr.poly_t_signal("ACGTTTTACG")["has_signal"] is True
    clean = crispr.poly_t_signal("ACGTACGTAC")
    assert clean["has_signal"] is False
    assert clean["max_run"] == 1


def test_self_complementarity_flags_hairpin():
    hairpin = crispr.self_complementarity("GGGGGGAAAACCCCCC")
    assert hairpin["max_stem"] >= 6
    assert hairpin["risk"] == "high"
    flat = crispr.self_complementarity("AAAAAAAAAAAA")
    assert flat["max_stem"] == 0
    assert flat["risk"] == "low"


def test_off_target_mismatch_summary_buckets():
    off_targets = [
        {"mismatches": 1, "risk_score": 0.5},
        {"mismatches": 1, "risk_score": 0.4},
        {"mismatches": 3, "risk_score": 0.1},
    ]
    summary = crispr.off_target_mismatch_summary(off_targets, max_mismatches=3)
    assert summary == {0: 0, 1: 2, 2: 0, 3: 1}


def test_local_specificity_proxy_bounds():
    assert crispr.local_specificity_proxy([]) == 1.0
    proxy = crispr.local_specificity_proxy(
        [{"risk_score": 0.5}, {"risk_score": 0.5}]
    )
    assert 0.0 < proxy < 1.0


def test_classify_guide_risk_levels():
    high = crispr.classify_guide_risk({0: 0, 1: 1, 2: 0, 3: 0}, 0.9)
    moderate = crispr.classify_guide_risk({0: 0, 1: 0, 2: 1, 3: 0}, 0.9)
    low = crispr.classify_guide_risk({0: 0, 1: 0, 2: 0, 3: 2}, 0.95)
    assert high["level"] == "High"
    assert moderate["level"] == "Moderate"
    assert low["level"] == "Low"


def test_composite_score_none_without_specificity():
    assert crispr.composite_score(0.8, None) is None
    combined = crispr.composite_score(0.8, 0.6)
    assert combined == pytest.approx(0.7)


def test_base_editing_window_detects_target_base():
    guide = {"guide_sequence": "AAACAAAAAAAAAAAAAAAA"}
    window = crispr.base_editing_window(guide, "CBE (SpCas9)")
    assert window["editor"] == "C>T"
    assert 4 in window["target_positions"]
    assert window["editable"] is True


def test_base_editing_window_rejects_non_base_editor():
    with pytest.raises(ValueError):
        crispr.base_editing_window({"guide_sequence": "ACGT"}, "SpCas9")


def test_prime_editing_notes_declares_incomplete_design():
    notes = crispr.prime_editing_notes(
        {"guide_sequence": "ACGTACGTACGTACGTACGT", "position": 1, "strand": "+"},
        "Prime Editor (SpCas9 nickase)",
    )
    assert notes["designed"] is False
    assert notes["recommended_tools"]


def test_evaluate_guides_ranks_and_enriches(guide):
    target = guide + "AGG" + "TTTT" + "GCGCGCGCGCGCGCGCGCGC" + "TGG"
    found = crispr.find_guides(target, "SpCas9")
    evaluated = crispr.evaluate_guides(found, target, "SpCas9", run_off_target=False)
    assert evaluated
    assert evaluated[0]["rank"] == 1
    assert "poly_t" in evaluated[0]
    assert "self_complementarity" in evaluated[0]
    assert evaluated[0]["specificity_proxy"] is None


def test_evaluate_guides_with_off_target_sets_specificity(guide):
    target = guide + "AGG" + "ACGTACGTGGGGCCCCATToff".replace("off", "")
    found = crispr.find_guides(target, "SpCas9")
    evaluated = crispr.evaluate_guides(found, target, "SpCas9", run_off_target=True)
    scored = [g for g in evaluated if g["specificity_proxy"] is not None]
    assert scored
    assert all(0.0 <= g["specificity_proxy"] <= 1.0 for g in scored)
    assert all(g["risk"] is not None for g in scored)


def test_evaluate_guides_defers_seed_beyond_detail_cap():
    target = ("ACGTACGTGG") * 80
    found = crispr.find_guides(target, "SpCas9")
    assert len(found) > crispr.MAX_GUIDES_DETAILED
    evaluated = crispr.evaluate_guides(found, target, "SpCas9", run_off_target=False)
    detailed = [g for g in evaluated if g.get("seed") is not None]
    skipped = [g for g in evaluated if g.get("seed") is None]
    assert len(detailed) == crispr.MAX_GUIDES_DETAILED
    assert skipped
    assert all(item.get("doench_score") is not None for item in evaluated)


def test_cut_site_position_none_for_cas13():
    guide = {"guide_sequence": "ACGUACGUACGUACGUACGUACGUACGU", "position": 1, "strand": "+"}
    assert crispr.cut_site_position(guide, "Cas13") is None


def test_generate_report_labels_efficiency_as_heuristic(guide):
    report = crispr.generate_report(
        [{"rank": 1, "guide_sequence": guide, "pam_sequence": "AGG", "position": 1, "strand": "+", "gc_content": 50.0, "doench_score": 0.7}],
        "Target",
    )
    assert "Doench_Score" not in report.columns
    assert "Efficiency_heuristic" in report.columns
    assert "Efficiency_method" in report.columns
    assert report.iloc[0]["Efficiency_heuristic"] == pytest.approx(0.7)
    assert "Doench" not in str(report.iloc[0]["Efficiency_method"]) or "not" in str(
        report.iloc[0]["Efficiency_method"]
    ).lower()


def test_example_emx1_is_published_and_has_spcas9_site():
    assert crispr.EXAMPLE_EMX1_SPCAS9 == "GAGTCCGAGCAGAAGAAGAAGGG"
    guides = crispr.find_guides(crispr.EXAMPLE_EMX1_SPCAS9, "SpCas9")
    spacers = [g["guide_sequence"] for g in guides if g["strand"] == "+"]
    assert "GAGTCCGAGCAGAAGAAGAA" in spacers
    assert any(g["pam_sequence"] == "GGG" for g in guides)


def test_diagnose_spcas9_explains_missing_ngg_and_short_input():
    agct = crispr.diagnose_spcas9_target("AGCT" * 20)
    assert agct["ok"] is False
    assert agct["ngg_total"] == 0
    assert any("NGG" in message for message in agct["messages"])
    short = crispr.diagnose_spcas9_target("ACGTACGT")
    assert short["too_short"] is True
    emx = crispr.diagnose_spcas9_target(crispr.EXAMPLE_EMX1_SPCAS9)
    assert emx["ok"] is True
    assert emx["guides_possible"] >= 1


def test_off_targets_require_adjacent_pam():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    decoy = "ACGT" * 5
    without_pam = guide + decoy
    assert crispr.find_off_targets(guide, without_pam, max_mismatches=3) == []
    with_pam = guide + "AGG" + decoy
    hits = crispr.find_off_targets(guide, with_pam, max_mismatches=3)
    assert all(h.get("pam_sequence") in {"AGG", "TGG", "CGG", "GGG", "AAG", "TAG", "CAG", "GAG"} for h in hits)


def test_seed_restriction_microhomology_and_order():
    guide = {"guide_sequence": "GAGTCCGAGCAGAAGAAGAA", "pam_sequence": "GGG"}
    seed = crispr.seed_occurrences(guide["guide_sequence"], crispr.EXAMPLE_EMX1_SPCAS9)
    assert seed["count"] >= 1
    assert "T" not in crispr.order_ready_oligos(guide)["spacer_rna"]
    assert "U" in crispr.order_ready_oligos(guide)["spacer_rna"]
    assert "EcoRI" in crispr.restriction_overlaps("AAAAGAATTCGGG")
    seq = "AAAAAACGT" + "ACGT" + "T" * 9
    mh = crispr.microhomology_at_cut(seq, 10, max_k=4)
    assert mh["length"] == 4
    assert mh["motif"] == "ACGT"


def test_hsu_single_hit_perfect_match_is_100():
    spacer = "GAGTCCGAGCAGAAGAAGAA"
    result = crispr.hsu_single_hit_score(spacer, spacer)
    assert result["status"] == "COMPUTED"
    assert result["score"] == pytest.approx(100.0)
    assert "single-hit" in result["method"]
    assert "genome-wide" not in result["method"].lower() or "not" in result["method"].lower()


def test_hsu_single_hit_one_pam_proximal_mismatch():
    spacer = "G" * 20
    off = "G" * 19 + "A"
    result = crispr.hsu_single_hit_score(spacer, off)
    expected = (1.0 - 0.583) * 100.0
    assert result["status"] == "COMPUTED"
    assert result["score"] == pytest.approx(expected)


def test_hsu_single_hit_rejects_non_20nt():
    result = crispr.hsu_single_hit_score("ACGTACGTACGTACGTA", "ACGTACGTACGTACGTA")
    assert result["status"] == "UNAVAILABLE"
    assert result["score"] is None


def test_validate_spcas9_guide_accepts_emx1_and_rejects_bad_pam():
    target = crispr.EXAMPLE_EMX1_SPCAS9
    guides = crispr.find_guides(target, "SpCas9")
    hit = next(g for g in guides if g["guide_sequence"] == "GAGTCCGAGCAGAAGAAGAA")
    ok = crispr.validate_spcas9_guide(hit, target_sequence=target)
    assert ok["valid"] is True
    bad = dict(hit)
    bad["pam_sequence"] = "TTC"
    failed = crispr.validate_spcas9_guide(bad, target_sequence=target)
    assert failed["valid"] is False
    short = {"guide_sequence": "ACGTACGTAC", "pam_sequence": "AGG", "position": 1, "strand": "+"}
    assert crispr.validate_spcas9_guide(short)["valid"] is False


def test_validate_spcas9_guide_rejects_wrong_strand_symbol():
    guide = {
        "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
        "pam_sequence": "GGG",
        "position": 1,
        "strand": "sense",
    }
    assert crispr.validate_spcas9_guide(guide)["valid"] is False


def test_find_off_targets_include_perfect_keeps_on_target_with_pam():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    target = guide + "GGG"
    assert crispr.find_off_targets(guide, target, max_mismatches=3) == []
    hits = crispr.find_off_targets(
        guide, target, max_mismatches=3, include_perfect=True
    )
    assert hits
    assert hits[0]["mismatches"] == 0
    assert hits[0]["pam_class"] == "NGG"
    assert hits[0]["hsu_hit_score"] == pytest.approx(100.0)


def test_find_guides_records_internal_zero_based_coordinates():
    target = "GAGTCCGAGCAGAAGAAGAAGGG"
    guides = crispr.find_guides(target, "SpCas9")
    hit = next(g for g in guides if g["strand"] == "+" and g["guide_sequence"].startswith("GAGT"))
    assert hit["start_0based"] == 0
    assert hit["end_0based"] == 20
    assert hit["coordinate_system"] == "internal_sense"


def test_offtarget_cache_key_changes_with_assembly_and_parameters():
    kwargs = dict(
        guide_sequence="GAGTCCGAGCAGAAGAAGAA",
        target_hash="aaa",
        reference_hash="bbb",
        assembly_declared="asmA",
        algorithm=crispr.ALGORITHM_SPCAS9_PAM_INDEX,
        model="Hsu2013_single_hit",
        model_version="CRISPOR-hitScoreM-20nt",
        max_mismatches=3,
        include_perfect=True,
        pam_mode="ngg_or_nag",
    )
    key_a = crispr.offtarget_cache_key(**kwargs)
    kwargs["assembly_declared"] = "asmB"
    key_b = crispr.offtarget_cache_key(**kwargs)
    kwargs["assembly_declared"] = "asmA"
    kwargs["max_mismatches"] = 2
    key_c = crispr.offtarget_cache_key(**kwargs)
    assert key_a != key_b
    assert key_a != key_c


def test_scientific_report_does_not_emit_zero_hits_when_not_run():
    guide = {
        "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
        "pam_sequence": "GGG",
        "position": 1,
        "strand": "+",
        "doench_score": 0.5,
        "efficiency_method": crispr.HEURISTIC_EFFICIENCY_METHOD,
    }
    report = crispr.scientific_crispr_report(
        target_sequence=crispr.EXAMPLE_EMX1_SPCAS9,
        guide=guide,
        search=None,
    )
    assert report["search_status"] == "NOT_RUN"
    assert report["verified_hit_count"] is None
    assert report["genome_wide"] is False
    assert any("not genome-wide" in item.lower() or "Genome-wide" in item or "genome-wide" in item for item in report["limitations"])


def test_crispr_ui_designs_emx1_without_autorunning_reference_search():
    app = open_module("crispr", timeout=60)
    assert not app.exception
    app.text_area(key="crispr_text").set_value(crispr.EXAMPLE_EMX1_SPCAS9).run()
    assert not app.exception
    app.button(key="crispr_button").click().run()
    assert not app.exception
    assert "crispr_reference_search" not in app.session_state
    keys = [item.key for item in app.button]
    assert "crispr_button" in keys
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("info", "error", "markdown", "caption", "warning")
        for item in list(getattr(app, bucket, []) or [])
    ).lower()
    assert "rule set 2" in combined
    assert "heuristic" in combined
    if "genome-wide" in combined:
        assert "not genome-wide" in combined or "not a public" in combined
    app.toggle(key="crispr_reference_panel").set_value(True).run()
    assert not app.exception
    keys_after = [item.key for item in app.button]
    assert "crispr_ref_search_button" in keys_after
    assert "crispr_reference_search" not in app.session_state

