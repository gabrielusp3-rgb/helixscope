"""Testes de identidade, normalizacao e validacao de REF de variantes.

Todos offline. Usam um FASTA sintetico pequeno para exercitar left-alignment e
REFERENCE_MISMATCH sem depender de uma assembly publica instalada.
"""

from __future__ import annotations

import os

import pytest

from modules import genome_fasta, variant_core

CONTIG = "HS_TEST_1"
# Homopolimero de A em 1-based 11-16 para testar left-alignment de indel.
SEQUENCE = "GATTACAGAT" + "AAAAAA" + "CGTACGTACG" + "TTTTGGGGCC"


@pytest.fixture()
def reference(tmp_path) -> tuple[str, str]:
    """Cria um FASTA sintetico indexado e devolve (fasta_path, fai_path)."""
    fasta_path = os.path.join(str(tmp_path), "tiny.fna")
    with open(fasta_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{CONTIG} synthetic test contig\n")
        for offset in range(0, len(SEQUENCE), 10):
            handle.write(SEQUENCE[offset : offset + 10] + "\n")
    fai_path = fasta_path + ".fai"
    genome_fasta.build_fai_for_fasta(fasta_path, fai_path=fai_path)
    return fasta_path, fai_path


def test_assembly_is_mandatory_for_identity() -> None:
    with pytest.raises(variant_core.VariantInputError) as exc:
        variant_core.build_variant(text="17 43093557 C G", assembly="")
    assert exc.value.category == "INVALID_INPUT"
    assert "assembly" in str(exc.value).lower()


def test_same_coordinate_in_two_assemblies_has_different_identity() -> None:
    grch38 = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    grch37 = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh37")
    assert grch38["identity_hash"] != grch37["identity_hash"]


def test_snv_identity_and_kind() -> None:
    variant = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    assert variant["status"] == "NORMALIZED"
    assert variant["variant_kind"] == variant_core.KIND_SNV
    assert variant["contig"] == "17"
    assert variant["position_1based"] == 43093557
    assert variant["position_0based"] == 43093556
    assert variant["ref"] == "C"
    assert variant["alt"] == "G"
    assert variant["effect"] is None
    assert variant["effect_status"] == "NOT_COMPUTED"
    assert variant["ref_validation"] == variant_core.REF_NOT_CHECKED


@pytest.mark.parametrize(
    "text,expected_format",
    [
        ("17 43093557 C G", "vcf_like"),
        ("17-43093557-C-G", "vcf_like"),
        ("NC_000017.11:43093556:C:G", "spdi"),
        ("NC_000017.11:g.43093557C>G", "hgvs_g_substitution"),
    ],
)
def test_supported_input_formats_agree_on_position(text: str, expected_format: str) -> None:
    parsed = variant_core.parse_variant_input(text)
    assert parsed["format"] == expected_format
    assert parsed["position_0based"] == 43093556
    assert parsed["ref"] == "C"
    assert parsed["alt"] == "G"


def test_rsid_and_vcv_are_identifier_only() -> None:
    for text, fmt in (("rs80357382", "rsid"), ("VCV000054425", "clinvar_vcv")):
        variant = variant_core.build_variant(text=text, assembly="GRCh38.p14")
        assert variant["status"] == "IDENTIFIER_ONLY"
        assert variant["format"] == fmt
        assert variant["position_0based"] is None
        assert variant["contig"] == ""


def test_hgvs_c_notation_is_refused_as_unsupported() -> None:
    with pytest.raises(variant_core.VariantInputError) as exc:
        variant_core.parse_variant_input("NM_007294.4:c.1974G>C")
    assert exc.value.category == "UNSUPPORTED"
    assert "house HGVS parser" in str(exc.value)


def test_trim_alleles_reduces_to_minimal_spdi_form() -> None:
    # CTT -> CT, TT -> T and T -> '' all describe deleting one T. The minimal
    # form is the only unique one, so it is what identity uses.
    trimmed = variant_core.trim_alleles(position_0based=100, ref="CTT", alt="CT")
    assert (trimmed["ref"], trimmed["alt"]) == ("T", "")
    # C at 0-based 100, T at 101, T at 102. Trimming right then left leaves the
    # T at 101 as the deleted base, which is the left-aligned choice.
    assert trimmed["position_0based"] == 101
    same = variant_core.trim_alleles(position_0based=100, ref="AGGG", alt="AGG")
    assert (same["ref"], same["alt"]) == ("G", "")


def test_anchored_and_minimal_input_reach_the_same_identity() -> None:
    anchored = variant_core.build_variant(text="17 100 CT C", assembly="GRCh38.p14")
    minimal = variant_core.build_variant(text="17 101 T .", assembly="GRCh38.p14")
    assert anchored["identity_hash"] == minimal["identity_hash"]
    assert anchored["variant_kind"] == variant_core.KIND_DELETION


def test_trim_refuses_identical_alleles() -> None:
    with pytest.raises(variant_core.VariantInputError):
        variant_core.trim_alleles(position_0based=5, ref="A", alt="A")


def test_classify_variant_kind_covers_indels() -> None:
    assert variant_core.classify_variant_kind("A", "G") == variant_core.KIND_SNV
    assert variant_core.classify_variant_kind("AT", "GC") == variant_core.KIND_MNV
    assert variant_core.classify_variant_kind("A", "AGG") == variant_core.KIND_INSERTION
    assert variant_core.classify_variant_kind("AGG", "A") == variant_core.KIND_DELETION
    assert variant_core.classify_variant_kind("", "GG") == variant_core.KIND_INSERTION
    assert variant_core.classify_variant_kind("GG", "") == variant_core.KIND_DELETION


def test_normalize_allele_rejects_ambiguous_bases() -> None:
    assert variant_core.normalize_allele("acgu") == "ACGT"
    assert variant_core.normalize_allele("-") == ""
    with pytest.raises(variant_core.VariantInputError):
        variant_core.normalize_allele("ACGN")
    with pytest.raises(variant_core.VariantInputError) as big:
        variant_core.normalize_allele("A" * (variant_core.MAX_ALLELE_NT + 1))
    assert big.value.category == "RESOURCE_LIMIT"


def test_ref_validation_accepts_matching_base(reference) -> None:
    fasta_path, fai_path = reference
    outcome = variant_core.validate_ref_against_reference(
        fasta_path=fasta_path, fai_path=fai_path, contig=CONTIG, position_0based=0, ref="G"
    )
    assert outcome["status"] == variant_core.REF_VALIDATED
    assert outcome["observed"] == "G"


def test_ref_mismatch_is_reported_not_silently_accepted(reference) -> None:
    fasta_path, fai_path = reference
    outcome = variant_core.validate_ref_against_reference(
        fasta_path=fasta_path, fai_path=fai_path, contig=CONTIG, position_0based=0, ref="T"
    )
    assert outcome["status"] == variant_core.REF_MISMATCH
    assert outcome["observed"] == "G"
    assert "wrong assembly" in outcome["reason"]


def test_build_variant_raises_reference_mismatch(reference) -> None:
    fasta_path, fai_path = reference
    with pytest.raises(variant_core.VariantInputError) as exc:
        variant_core.build_variant(
            text=f"{CONTIG} 1 T A",
            assembly="HELIXSCOPE_TEST",
            fasta_path=fasta_path,
            fai_path=fai_path,
        )
    assert exc.value.category == "REFERENCE_MISMATCH"


def test_build_variant_validates_ref_when_reference_present(reference) -> None:
    fasta_path, fai_path = reference
    variant = variant_core.build_variant(
        text=f"{CONTIG} 1 G A",
        assembly="HELIXSCOPE_TEST",
        fasta_path=fasta_path,
        fai_path=fai_path,
    )
    assert variant["ref_validation"] == variant_core.REF_VALIDATED
    assert variant["normalization"] == variant_core.NORMALIZATION_TRIMMED
    assert variant["ref_observed"] == "G"


def test_left_alignment_shifts_indel_in_homopolymer(reference) -> None:
    fasta_path, fai_path = reference
    # SEQUENCE 1-based 11..16 is AAAAAA. Deleting the last A must left-align to 11.
    variant = variant_core.build_variant(
        text=f"{CONTIG} 15 AA A",
        assembly="HELIXSCOPE_TEST",
        fasta_path=fasta_path,
        fai_path=fai_path,
    )
    assert variant["variant_kind"] == variant_core.KIND_DELETION
    assert variant["left_shifted_nt"] > 0
    assert variant["position_1based"] < 15
    assert variant["ref_validation"] == variant_core.REF_VALIDATED


def test_left_align_is_noop_for_substitution(reference) -> None:
    fasta_path, fai_path = reference
    outcome = variant_core.left_align_indel(
        fasta_path=fasta_path,
        fai_path=fai_path,
        contig=CONTIG,
        position_0based=0,
        ref="G",
        alt="A",
    )
    assert outcome["shifted_nt"] == 0
    assert outcome["position_0based"] == 0


def test_unknown_contig_is_not_silently_mapped(reference) -> None:
    fasta_path, fai_path = reference
    with pytest.raises(variant_core.VariantInputError) as exc:
        variant_core.validate_ref_against_reference(
            fasta_path=fasta_path,
            fai_path=fai_path,
            contig="chr17",
            position_0based=0,
            ref="G",
        )
    assert exc.value.category == "INVALID_INPUT"
    assert "does not silently map" in str(exc.value)


def test_position_outside_contig_is_invalid(reference) -> None:
    fasta_path, fai_path = reference
    with pytest.raises(variant_core.VariantInputError):
        variant_core.validate_ref_against_reference(
            fasta_path=fasta_path,
            fai_path=fai_path,
            contig=CONTIG,
            position_0based=10_000,
            ref="G",
        )


def test_variant_region_string_for_substitution_and_insertion() -> None:
    snv = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    assert variant_core.variant_region_string(snv) == "17:43093557-43093557:1"
    insertion = variant_core.build_variant(text="17 43093557 . GG", assembly="GRCh38.p14")
    assert insertion["variant_kind"] == variant_core.KIND_INSERTION
    assert variant_core.variant_region_string(insertion) == "17:43093557-43093556:1"


def test_variant_region_string_refuses_identifier_only() -> None:
    identifier = variant_core.build_variant(text="rs80357382", assembly="GRCh38.p14")
    with pytest.raises(variant_core.VariantInputError):
        variant_core.variant_region_string(identifier)


def test_ensembl_contig_name_only_strips_canonical_chr() -> None:
    assert variant_core.ensembl_contig_name("chr17") == "17"
    assert variant_core.ensembl_contig_name("17") == "17"
    assert variant_core.ensembl_contig_name("chrX") == "X"
    assert variant_core.ensembl_contig_name("chrM") == "MT"
    assert variant_core.ensembl_contig_name("NC_000017.11") == "NC_000017.11"
    assert variant_core.ensembl_contig_name("chrUn_KI270302v1") == "chrUn_KI270302v1"


def test_identity_hash_is_stable_and_allele_sensitive() -> None:
    base = dict(assembly="GRCh38", contig="17", position_0based=1, ref="C", alt="G")
    first = variant_core.variant_identity_hash(**base)
    assert first == variant_core.variant_identity_hash(**base)
    other = variant_core.variant_identity_hash(**{**base, "alt": "T"})
    assert first != other


@pytest.mark.parametrize(
    "payload",
    [
        "",
        "   ",
        "17 43093557",
        "17 43093557 C",
        "17 notanumber C G",
        "17 43093557 C G extra field here",
        "../../etc/passwd 1 A T",
        "17 43093557 C G; rm -rf /",
        "<script>alert(1)</script>",
        "17\t43093557\tC\tG\nchr1 1 A T",
    ],
)
def test_malformed_and_hostile_input_is_refused(payload: str) -> None:
    with pytest.raises(variant_core.VariantInputError):
        variant_core.build_variant(text=payload, assembly="GRCh38.p14")


def test_oversized_input_hits_resource_limit() -> None:
    with pytest.raises(variant_core.VariantInputError) as exc:
        variant_core.parse_variant_input("A" * (variant_core.MAX_INPUT_CHARS + 1))
    assert exc.value.category == "RESOURCE_LIMIT"


def test_contig_with_path_separators_is_refused() -> None:
    with pytest.raises(variant_core.VariantInputError):
        variant_core.build_variant(
            text="17 43093557 C G", assembly="GRCh38.p14", contig_override="../../secret"
        )


def test_variant_from_mapping_still_available_for_legacy_callers() -> None:
    from modules import genomic_variant

    variant = genomic_variant.variant_from_mapping(
        {"assembly": "GRCh38.p14", "contig": "17", "position_0based": 10, "ref": "A", "alt": "G"}
    )
    assert variant["effect"] is None
