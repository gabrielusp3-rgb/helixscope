"""Variant identity only. No effect prediction."""

from __future__ import annotations

import pytest

from modules import genomic_variant


def test_genomic_variant_identity_no_effect():
    row = genomic_variant.genomic_variant(
        assembly="GRCh38.p14",
        accession="GCF_000001405.40",
        contig="NC_000007.14",
        position_0based=100,
        ref="A",
        alt="G",
        assembly_sha256="abc",
    )
    assert row["effect"] is None
    assert row["effect_status"] == "NOT_COMPUTED"
    assert row["position_1based"] == 101
    assert row["ref"] == "A"
    assert row["alt"] == "G"


def test_variant_rejects_path_contig_and_bad_alleles():
    with pytest.raises(genomic_variant.VariantError):
        genomic_variant.genomic_variant(
            assembly="GRCh38.p14",
            accession="GCF_000001405.40",
            contig="../etc/passwd",
            position_0based=0,
            ref="A",
            alt="T",
        )
    with pytest.raises(genomic_variant.VariantError):
        genomic_variant.genomic_variant(
            assembly="GRCh38.p14",
            accession="GCF_000001405.40",
            contig="NC_000001.11",
            position_0based=0,
            ref="X",
            alt="A",
        )
