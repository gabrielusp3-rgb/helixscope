"""Testes da busca de off-target verificavel em referencia fornecida."""

import pytest

from modules import crispr, crispr_offtarget, crispr_reference


EMX1_SPACER = "GAGTCCGAGCAGAAGAAGAA"
EMX1_SITE = crispr.EXAMPLE_EMX1_SPCAS9


def _emx1_guide() -> dict:
    guides = crispr.find_guides(EMX1_SITE, "SpCas9")
    return next(g for g in guides if g["guide_sequence"] == EMX1_SPACER)


def test_reference_search_finds_published_emx1_perfect_site():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", organism="Homo sapiens", assembly=""
    )
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=3, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    assert result["genome_wide"] is False
    assert "genome-wide" not in str(result.get("scope_label") or "").lower() or "not" in str(
        result.get("scope_label") or ""
    ).lower()
    assert result["verified_hit_count"] >= 1
    perfect = [h for h in result["hits"] if h["mismatches"] == 0 and h["PAM"] == "GGG"]
    assert perfect
    hit = perfect[0]
    assert hit["target_sequence"] == EMX1_SPACER
    assert hit["chromosome_or_contig"] == "emx1"
    assert hit["score"] == pytest.approx(100.0)
    assert "Hsu" in hit["score_method"]
    assert crispr_offtarget.validate_hit_against_reference(hit, reference)


def test_similar_sequence_without_pam_is_not_a_hit():
    decoy = EMX1_SPACER + "ACGTACGTACGT"
    fasta = ">decoy\n" + decoy + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=3, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    assert result["verified_hit_count"] == 0
    assert result["hits"] == []
    assert "not a safety claim" in result["reason"].lower() or "0 verified" in result["reason"]


def test_one_and_multiple_mismatches_are_counted_not_estimated():
    # one mismatch at PAM-distal base, still NGG
    off = "A" + EMX1_SPACER[1:] + "GGG"
    fasta = ">mm1\n" + off + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=3, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    ones = [h for h in result["hits"] if h["mismatches"] == 1]
    assert ones
    assert ones[0]["mismatch_positions"] == [1]
    assert ones[0]["PAM"] == "GGG"


def test_reverse_strand_hit_is_extracted_from_real_contig():
    from modules.crispr import _reverse_complement

    spacer_pam = EMX1_SPACER + "GGG"
    sense = "TTTT" + _reverse_complement(spacer_pam) + "AAAA"
    fasta = ">rev\n" + sense + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=0, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    minus = [h for h in result["hits"] if h["strand"] == "-" and h["mismatches"] == 0]
    assert minus
    hit = minus[0]
    assert hit["target_sequence"] == EMX1_SPACER
    assert hit["PAM"] == "GGG"
    assert crispr_offtarget.validate_hit_against_reference(hit, reference)


def test_invalid_pam_window_never_becomes_verified_hit():
    fasta = ">x\n" + EMX1_SPACER + "TTC" + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=3, include_perfect=True
    )
    assert result["verified_hit_count"] == 0


def test_verification_failure_is_error_not_zero_hits(monkeypatch):
    fasta = ">x\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")

    def boom(**kwargs):
        raise crispr.CrisprError(
            "Index spacer/PAM does not match contig extraction.",
            "PARSING_ERROR",
        )

    monkeypatch.setattr(crispr_offtarget, "_verified_hit", boom)
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=0, include_perfect=True
    )
    assert result["status"] == "ERROR"
    assert result["verified_hit_count"] is None
    assert result["hits"] == []
    assert result["unverified_count"] >= 1
    assert "zero off-targets" in result["reason"].lower() or "re-extraction" in result["reason"].lower()


def test_wrong_length_guide_is_invalid_input_not_zero_hits():
    guide = {
        "guide_sequence": "ACGTACGTAC",
        "pam_sequence": "AGG",
        "position": 1,
        "strand": "+",
    }
    fasta = ">x\nGAGTCCGAGCAGAAGAAGAAGGG\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(guide, reference)
    assert result["status"] == "INVALID_INPUT"
    assert result["verified_hit_count"] is None
    assert result["hits"] == []


def test_export_refuses_not_run_error_timeout_as_zero_hits():
    with pytest.raises(crispr.CrisprError) as exc:
        crispr_offtarget.export_offtarget_rows(crispr_offtarget.not_run_envelope())
    assert "NOT_RUN" in str(exc.value)
    timeout = {
        "status": "TIMEOUT",
        "hits": [],
        "verified_hit_count": None,
    }
    with pytest.raises(crispr.CrisprError):
        crispr_offtarget.export_offtarget_rows(timeout)
    failed = {"status": "ERROR", "hits": [], "verified_hit_count": None}
    with pytest.raises(crispr.CrisprError):
        crispr_offtarget.export_offtarget_rows(failed)


def test_export_completed_zero_hits_is_empty_list_not_a_fake_row():
    fasta = ">x\nACGTACGTACGTACGTACGTACGT\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=1, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    assert result["verified_hit_count"] == 0
    rows = crispr_offtarget.export_offtarget_rows(result)
    assert rows == []


def test_hit_cap_is_resource_limit_not_a_complete_catalogue():
    copies = (EMX1_SITE + "AAAA") * 3
    fasta = ">multi\n" + copies + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(),
        reference,
        max_mismatches=3,
        include_perfect=True,
        max_hits=1,
    )
    assert result["status"] == "RESOURCE_LIMIT"
    assert result["truncated"] is True
    assert result["verified_hit_count"] is None
    with pytest.raises(crispr.CrisprError) as exc:
        crispr_offtarget.export_offtarget_rows(result)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_timeout_is_not_converted_to_zero_hits():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, timeout_s=0.0, include_perfect=True
    )
    assert result["status"] == "TIMEOUT"
    assert result["verified_hit_count"] is None
    assert result["hits"] == []
    assert result["truncated"] is True


def test_result_identity_changes_when_assembly_label_changes():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    ref_a = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", assembly="plasmid-v1"
    )
    ref_b = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", assembly="plasmid-v2"
    )
    a = crispr_offtarget.search_off_targets_in_reference(_emx1_guide(), ref_a)
    b = crispr_offtarget.search_off_targets_in_reference(_emx1_guide(), ref_b)
    assert a["result_identity"] != b["result_identity"]
    assert a["reference_identity_hash"] == b["reference_identity_hash"]


def test_parameter_change_changes_result_identity():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    a = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=1
    )
    b = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=3
    )
    assert a["result_identity"] != b["result_identity"]


def test_ncbi_feature_overlap_does_not_claim_gene_disrupted():
    hit = {"start_0based": 10, "end_0based": 30}
    features = [
        {
            "status": "AVAILABLE",
            "type": "gene",
            "gene": "EMX1",
            "product": "",
            "start": 0,
            "end": 50,
        }
    ]
    rows = crispr_offtarget.overlapping_ncbi_features(hit, features, contig_length=80)
    assert rows
    assert "not a demonstration" in rows[0]["claim"].lower()
    assert "disrupted" in rows[0]["claim"].lower()


def test_search_performance_small_reference_is_measured():
    fasta = ">emx1\n" + (EMX1_SITE + "ACGT") * 20 + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, include_perfect=True
    )
    assert result["status"] in {"COMPLETED", "RESOURCE_LIMIT"}
    assert result["elapsed_ms"] is not None
    assert result["elapsed_ms"] < 5000.0
    assert result["index_build_ms"] is not None
    assert result["search_ms"] is not None


def test_offtarget_chart_requires_real_hits():
    from ui import charts

    with pytest.raises(ValueError):
        charts.off_target_contig_histogram([])
    fig = charts.off_target_contig_histogram(
        [{"chromosome_or_contig": "c1"}, {"chromosome_or_contig": "c1"}]
    )
    assert list(fig.data[0].y) == [2]
    with pytest.raises(ValueError):
        charts.off_target_pam_histogram({})
    pam_fig = charts.off_target_pam_histogram({"NGG": 2, "NAG_reduced_activity": 1})
    assert list(pam_fig.data[0].y) == [1, 2]
    with pytest.raises(ValueError):
        charts.off_target_score_histogram([], title="x", x_title="y")
    score_fig = charts.off_target_score_histogram(
        [100.0, 50.0], title="Hsu 2013 single-hit scores (0-100)", x_title="Hsu"
    )
    assert score_fig.layout.title.text.startswith("Hsu")
