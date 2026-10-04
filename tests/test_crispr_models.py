"""Testes de CFD, agregado MIT/CFD, catalogo de assemblies e Cas-OFFinder."""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import pytest

from modules import (
    crispr,
    crispr_assemblies,
    crispr_casoffinder,
    crispr_cfd,
    crispr_offtarget,
    crispr_reference,
    crispr_specificity,
    crispr_structure_catalog,
    provenance,
)


EMX1_SPACER = "GAGTCCGAGCAGAAGAAGAA"
EMX1_SITE = crispr.EXAMPLE_EMX1_SPCAS9


def _emx1_guide() -> dict:
    guides = crispr.find_guides(EMX1_SITE, "SpCas9")
    return next(g for g in guides if g["guide_sequence"] == EMX1_SPACER)


def test_cfd_crispor_doctest_vectors():
    perfect = crispr_cfd.cfd_pair_score("G" * 23, "G" * 23)
    assert perfect["status"] == "COMPUTED"
    assert perfect["score"] == pytest.approx(1.0)
    assert perfect["scale"] == "0-1"
    assert perfect["model"] == "Doench2016_CFD"
    mm = crispr_cfd.cfd_pair_score("G" * 23, "G" * 17 + "AAA" + "GGG")
    assert mm["score"] == pytest.approx(0.4635989007074176)
    assert crispr_cfd.cfd_pair_score(
        "ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG"
    )["score"] == pytest.approx(0.5)
    assert crispr_cfd.cfd_pair_score(
        "ATGGTCGGACTCCCTGCCAGAGG", "ATGATCCAAATCCCTGCCAGAGG"
    )["score"] == pytest.approx(0.53625000020625)
    assert crispr_cfd.cfd_pair_score(
        "ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG"
    )["score"] == pytest.approx(0.384615385)


def test_cfd_wrong_length_is_unavailable_not_zero():
    result = crispr_cfd.cfd_pair_score("ACGT", "ACGT")
    assert result["status"] == "UNAVAILABLE"
    assert result["score"] is None


def test_cfd_non_dna_is_invalid_input_not_zero():
    result = crispr_cfd.cfd_pair_score("N" * 23, "G" * 23)
    assert result["status"] == "INVALID_INPUT"
    assert result["score"] is None


def test_cfd_uses_embedded_tables_not_pickle():
    record = crispr_cfd.cfd_model_record()
    assert record["pickle_used"] is False
    assert len(crispr_cfd.CFD_MM_SCORES) == 240
    assert len(crispr_cfd.CFD_PAM_SCORES) == 16
    assert crispr_cfd.CFD_PAM_SCORES["GG"] == pytest.approx(1.0)
    source = inspect.getsource(crispr_cfd)
    assert "pickle" in source.lower()
    assert "no pickle" in source.lower() or "pickle_used" in source


def test_mit_sguide_formula_and_crispor_integer():
    empty = crispr_specificity.mit_sguide_from_hit_sum(0.0)
    assert empty["score"] == pytest.approx(100.0)
    assert empty["score_crispor_integer"] == 100
    one_perfect = crispr_specificity.mit_sguide_from_hit_sum(100.0)
    assert one_perfect["score"] == pytest.approx(50.0)
    assert one_perfect["score_crispor_integer"] == 50
    assert "genome-wide" in one_perfect["method"].lower()
    assert "not" in one_perfect["method"].lower()


def test_cfd_sguide_formula():
    none_ot = crispr_specificity.cfd_sguide_from_cfd_sum(0.0)
    assert none_ot["score"] == pytest.approx(100.0)
    one = crispr_specificity.cfd_sguide_from_cfd_sum(1.0)
    assert one["score"] == pytest.approx(50.0)
    assert one["scale"] == "0-100"


def test_nag_hsu_is_downweighted_for_mit_not_for_cfd():
    hit = {"score": 100.0, "PAM": "GAG", "cfd_score": 0.25925925899999996}
    weighted = crispr_specificity.mit_weighted_hsu(hit)
    assert weighted == pytest.approx(20.0)
    gg = crispr_specificity.mit_weighted_hsu({"score": 100.0, "PAM": "GGG"})
    assert gg == pytest.approx(100.0)


def test_completed_search_excludes_intended_and_keeps_extra_exact():
    fasta = ">dup\n" + EMX1_SITE + "AAAA" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, max_mismatches=0, include_perfect=True
    )
    assert result["status"] == "COMPLETED"
    roles = [hit["site_role"] for hit in result["hits"]]
    assert roles.count("intended_on_target") == 1
    assert roles.count("additional_exact_site") == 1
    mit = result["guide_specificity"]["mit_sguide"]
    assert mit["status"] == "COMPUTED"
    assert mit["n_intended_excluded"] == 1
    assert mit["n_additional_exact"] == 1
    assert mit["score"] == pytest.approx(50.0)
    assert mit["genome_wide"] is False
    cfd = result["guide_specificity"]["cfd_sguide"]
    assert cfd["status"] == "COMPUTED"
    assert cfd["score"] == pytest.approx(50.0)
    assert result["hits"][0]["cfd_score"] == pytest.approx(1.0)


def test_not_run_and_timeout_do_not_become_specificity_zero():
    empty = crispr_specificity.scores_from_search({"status": "NOT_RUN", "hits": []})
    assert empty["mit_sguide"]["score"] is None
    assert empty["cfd_sguide"]["score"] is None
    assert empty["mit_sguide"]["status"] == "NOT_RUN"
    timeout = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(),
        crispr_reference.load_reference_from_text(
            ">emx1\n" + EMX1_SITE + "\n", source="user FASTA"
        ),
        timeout_s=0.0,
        include_perfect=True,
    )
    assert timeout["status"] == "TIMEOUT"
    assert timeout["guide_specificity"]["mit_sguide"]["score"] is None
    assert timeout["guide_specificity"]["cfd_sguide"]["score"] is None


def test_resource_limit_does_not_emit_complete_sguide():
    copies = (EMX1_SITE + "AAAA") * 3
    reference = crispr_reference.load_reference_from_text(
        ">multi\n" + copies + "\n", source="user FASTA"
    )
    result = crispr_offtarget.search_off_targets_in_reference(
        _emx1_guide(), reference, include_perfect=True, max_hits=1
    )
    assert result["status"] == "RESOURCE_LIMIT"
    assert result["guide_specificity"]["mit_sguide"]["score"] is None


def test_hg38_declaration_is_not_verified_assembly():
    fasta = ">chr1\n" + EMX1_SITE + "\n"
    loaded = crispr_reference.load_reference_from_text(
        fasta, source="user FASTA", organism="Homo sapiens", assembly="hg38"
    )
    assert loaded["assembly_declared"] == "hg38"
    assert loaded["filename_trusted_as_assembly"] is False
    assert loaded["verified_assembly"] == ""
    assert loaded["assembly_identity_status"] == "nickname_matches_catalog_not_verified"
    assert loaded["assembly_catalog"]["catalog_match"]["accession"] == "GCF_000001405.40"
    assert loaded["assembly_catalog"]["checksum"] == ""
    assert loaded["genome_wide"] is False


def test_matching_accession_sets_verified_assembly_name_not_fasta_proof():
    lookup = crispr_assemblies.lookup_declaration(
        "", accession="GCF_000001405.40"
    )
    assert lookup["verified_assembly"] == "GRCh38.p14"
    assert lookup["assembly_identity_status"] == "accession_matches_catalog_not_fasta"
    assert lookup["downloaded"] is False


def test_cas_offinder_missing_is_unavailable_not_zero_hits():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")
    result = crispr_casoffinder.search_with_cas_offinder(
        _emx1_guide(),
        reference,
        which_info={"available": False, "reason": "not installed", "path": ""},
    )
    assert result["status"] == "UNAVAILABLE"
    assert result["verified_hit_count"] is None
    assert result["hits"] == []
    assert result["genome_wide"] is False


def test_cas_offinder_mocked_output_is_reextracted():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")

    def run_fn(args, **_kwargs):
        output_path = args[3]
        with open(args[1], encoding="utf-8") as handle:
            body = handle.read()
        if "NGG" in body.splitlines()[1]:
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "GAGTCCGAGCAGAAGAAGAAGGG emx1 0 "
                    "GAGTCCGAGCAGAAGAAGAAGGG + 0\n"
                )
        else:
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write("")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    result = crispr_casoffinder.search_with_cas_offinder(
        _emx1_guide(),
        reference,
        include_perfect=True,
        max_mismatches=0,
        run_fn=run_fn,
        which_info={
            "available": True,
            "path": "cas-offinder",
            "version": "test-mock",
            "report_path": "cas-offinder",
        },
    )
    assert result["status"] == "COMPLETED"
    assert result["algorithm"] == crispr_casoffinder.ALGORITHM_CAS_OFFINDER
    assert result["verified_hit_count"] == 1
    hit = result["hits"][0]
    assert hit["target_sequence"] == EMX1_SPACER
    assert hit["PAM"] == "GGG"
    assert hit["cfd_score"] == pytest.approx(1.0)
    assert result["genome_wide"] is False
    assert result["tool"]["bulges_supported"] is False


def test_cas_offinder_subprocess_does_not_use_shell():
    source = inspect.getsource(crispr_casoffinder._run_queries)
    assert "shell=False" in source
    assert "shell=True" not in source
    assert "DEVICE_CPU" in source or '"C"' in source


def test_cas_offinder_malformed_output_is_parsing_error_not_zero():
    fasta = ">emx1\n" + EMX1_SITE + "\n"
    reference = crispr_reference.load_reference_from_text(fasta, source="user FASTA")

    def run_fn(args, **_kwargs):
        with open(args[3], "w", encoding="utf-8") as handle:
            handle.write("not-a-cas-offinder-line\n")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    result = crispr_casoffinder.search_with_cas_offinder(
        _emx1_guide(),
        reference,
        run_fn=run_fn,
        which_info={"available": True, "path": "cas-offinder", "report_path": "cas-offinder"},
    )
    assert result["status"] == "PARSING_ERROR"
    assert result["verified_hit_count"] is None


def test_structure_catalog_does_not_invent_coordinates():
    items = crispr_structure_catalog.list_experimental_complexes()
    assert {item["pdb_id"] for item in items} >= {"4UN3", "4OO8"}
    contract = crispr_structure_catalog.mapping_contract(guide=_emx1_guide())
    assert contract["status"] == "UNAVAILABLE"
    assert contract["structure_mapping"] is None
    assert "xyz" not in contract
    assert contract["helixscope_3d"] == "FETCH_ON_DEMAND"


def test_include_perfect_defaults_unchanged_for_local_and_reference():
    local_params = inspect.signature(crispr.find_off_targets).parameters
    assert local_params["include_perfect"].default is False
    ref_params = inspect.signature(
        crispr_offtarget.search_off_targets_in_reference
    ).parameters
    assert ref_params["include_perfect"].default is True


def test_software_version_is_phase_8b():
    assert provenance.HELIXSCOPE_VERSION.startswith("0.")
