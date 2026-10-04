"""Testes unitarios do modulo modules.rna_analysis.

Cobrem a transcricao de DNA em RNA, a deteccao de numero de acesso, as notas sobre
genomas com mecanismos especiais de traducao, a conversao de ORFs em regioes, a
extracao de codons a partir de regioes codificadoras, o pipeline completo de
resolucao de CDS, a consistencia da tabela de uso de codons e o tratamento de
organismo invalido no Codon Adaptation Index. Nenhum teste faz chamadas de rede.
"""

import math

import pytest

from modules import rna_analysis


@pytest.fixture
def coding_dna() -> str:
    """Fita codificante de DNA reutilizavel nos testes."""
    return "ATGAAATTTGGG"


@pytest.fixture
def mrna_sequence() -> str:
    """Sequencia de mRNA com numero inteiro de codons."""
    return "AUGAAAUUUGGG"


@pytest.fixture
def single_orf_genome() -> str:
    """Sequencia de 156 nt contendo exatamente uma ORF completa."""
    return "ATG" + "AAA" * 50 + "TAA"


@pytest.fixture
def nested_orf_genome() -> str:
    """Sequencia de 159 nt com duas ORFs aninhadas no mesmo frame e stop."""
    return "ATG" + "AAA" * 20 + "ATG" + "AAA" * 30 + "TAA"


@pytest.fixture
def non_coding_sequence() -> str:
    """Sequencia de 400 nt sem nenhum codon de iniciacao ATG."""
    return "ACGT" * 100


def test_transcribe_dna_to_rna(coding_dna):
    assert rna_analysis.transcribe(coding_dna) == "AUGAAAUUUGGG"


def test_transcribe_rejects_rna_input():
    with pytest.raises(ValueError):
        rna_analysis.transcribe("AUGCAUGC")


def test_codon_usage_count_matches_total_codons(mrna_sequence):
    table = rna_analysis.codon_usage_table(mrna_sequence)
    expected_codons = len(mrna_sequence) // 3
    assert int(table["Count"].sum()) == expected_codons


def test_codon_usage_table_exposes_unambiguous_columns(mrna_sequence):
    table = rna_analysis.codon_usage_table(mrna_sequence)
    assert list(table.columns) == [
        "Codon",
        "AminoAcid",
        "Count",
        "Synonymous_pct",
        "Absolute_pct",
        "RSCU",
    ]
    assert "RSCU_pct" not in table.columns


def test_synonymous_family_sizes_matches_genetic_code():
    sizes = rna_analysis.synonymous_family_sizes()
    assert sizes["L"] == 6
    assert sizes["R"] == 6
    assert sizes["S"] == 6
    assert sizes["I"] == 3
    assert sizes["A"] == 4
    assert sizes["M"] == 1
    assert sizes["W"] == 1
    assert sizes["*"] == 3
    assert sum(sizes.values()) == 64


def test_rscu_maximum_for_exclusive_codon():
    table = rna_analysis.relative_synonymous_codon_usage("CUG" * 4)
    leucine = table[table["AminoAcid"] == "L"]
    assert len(leucine) == 6
    assert float(leucine[leucine["Codon"] == "CUG"]["RSCU"].iloc[0]) == 6.0
    assert float(leucine["RSCU"].sum()) == pytest.approx(6.0)


def test_rscu_equals_one_without_bias():
    table = rna_analysis.relative_synonymous_codon_usage("AAAAAG")
    lysine = table[table["AminoAcid"] == "K"]
    assert len(lysine) == 2
    assert all(float(value) == 1.0 for value in lysine["RSCU"])


def test_rscu_excludes_stop_codons():
    table = rna_analysis.relative_synonymous_codon_usage("AAAUAA")
    assert "*" not in set(table["AminoAcid"])


def test_rscu_rejects_stop_only_sequence():
    with pytest.raises(ValueError):
        rna_analysis.relative_synonymous_codon_usage("UAAUAG")


def test_effective_number_of_codons_maximum_bias():
    biased = "AAA" * 5 + "GCU" * 5 + "CUG" * 5
    assert rna_analysis.effective_number_of_codons(biased) == rna_analysis.ENC_MIN


def test_effective_number_of_codons_stays_in_range():
    mixed = ("AAAAAG" * 4) + ("GCUGCCGCAGCG" * 3) + ("CUGCUCUUAUUG" * 3)
    enc = rna_analysis.effective_number_of_codons(mixed)
    assert rna_analysis.ENC_MIN <= enc <= rna_analysis.ENC_MAX


def test_effective_number_of_codons_requires_enough_families():
    with pytest.raises(ValueError):
        rna_analysis.effective_number_of_codons("AUGAUGAUG")


def test_codon_adaptation_index_invalid_organism(mrna_sequence):
    with pytest.raises(ValueError):
        rna_analysis.codon_adaptation_index(mrna_sequence, organism="martian")


def test_detect_accession_from_fasta_header():
    header = ">NC_045512.2 Severe acute respiratory syndrome coronavirus 2\nATGCATGC"
    assert rna_analysis.detect_accession(header) == "NC_045512.2"


def test_detect_accession_from_bare_identifier():
    assert rna_analysis.detect_accession("MN908947.3") == "MN908947.3"


def test_detect_accession_ignores_plain_sequence(non_coding_sequence):
    assert rna_analysis.detect_accession(non_coding_sequence) is None
    assert rna_analysis.detect_accession("ATGCATGCATGC") is None
    assert rna_analysis.detect_accession("") is None


def test_special_genome_notes_by_accession():
    notes = rna_analysis.special_genome_notes(accession="NC_045512.2")
    assert notes
    assert any("frameshift" in note.lower() for note in notes)


def test_special_genome_notes_by_organism_and_product():
    organism_notes = rna_analysis.special_genome_notes(
        organism="Human immunodeficiency virus 1"
    )
    product_notes = rna_analysis.special_genome_notes(
        product_names=["ORF1ab polyprotein"]
    )
    assert organism_notes
    assert product_notes
    assert rna_analysis.special_genome_notes() == []


def test_orfs_to_regions_projects_negative_frame():
    orfs = [
        {"frame": 1, "start": 0, "end": 30, "length_bp": 30},
        {"frame": -1, "start": 0, "end": 30, "length_bp": 30},
    ]
    regions = rna_analysis.orfs_to_regions(orfs, seq_length=100)
    assert regions[0]["segments"] == [(0, 30)]
    assert regions[0]["strand"] == "+"
    assert regions[1]["segments"] == [(70, 100)]
    assert regions[1]["strand"] == "-"


def test_orfs_to_regions_rejects_invalid_length():
    with pytest.raises(ValueError):
        rna_analysis.orfs_to_regions([], seq_length=0)


def test_codons_from_regions_forward_and_reverse():
    forward = rna_analysis.codons_from_regions(
        "ATGAAATTTTAG", [{"segments": [(0, 12)], "strand": "+"}]
    )
    reverse = rna_analysis.codons_from_regions(
        "ATGAAA", [{"segments": [(0, 6)], "strand": "-"}]
    )
    assert forward == "AUGAAAUUUUAG"
    assert reverse == "UUUCAU"


def test_codons_from_regions_respects_codon_start():
    codons = rna_analysis.codons_from_regions(
        "ATGAAATTTTAG",
        [{"segments": [(0, 12)], "strand": "+", "codon_start": 2}],
    )
    assert len(codons) == 9
    assert codons == "UGAAAUUUU"


def test_codons_from_regions_without_regions():
    assert rna_analysis.codons_from_regions("ATGAAA", []) == ""


def test_resolve_coding_codons_prediction(single_orf_genome):
    resolution = rna_analysis.resolve_coding_codons(
        single_orf_genome, source="prediction", min_cds_length=150
    )
    assert resolution["method"] == rna_analysis.CDS_SOURCE_PREDICTION
    assert resolution["n_cds"] == 1
    assert resolution["n_codons"] == 52
    assert resolution["n_codons"] * 3 <= resolution["total_input_nt"]


def test_resolve_coding_codons_never_exceeds_sequence_when_filtering(
    nested_orf_genome,
):
    filtered = rna_analysis.resolve_coding_codons(
        nested_orf_genome,
        source="prediction",
        min_cds_length=90,
        drop_overlapping=True,
    )
    assert filtered["n_cds"] == 1
    assert filtered["n_codons"] * 3 <= filtered["total_input_nt"]

    unfiltered = rna_analysis.resolve_coding_codons(
        nested_orf_genome,
        source="prediction",
        min_cds_length=90,
        drop_overlapping=False,
    )
    assert unfiltered["n_codons"] * 3 > unfiltered["total_input_nt"]
    assert any("sobrepostas" in message for message in unfiltered["warnings"])


def test_resolve_coding_codons_refuses_silent_fallback(non_coding_sequence):
    with pytest.raises(ValueError):
        rna_analysis.resolve_coding_codons(
            non_coding_sequence, source="prediction", min_cds_length=150
        )


def test_resolve_coding_codons_whole_sequence_is_explicit(non_coding_sequence):
    resolution = rna_analysis.resolve_coding_codons(
        non_coding_sequence, source="whole"
    )
    assert resolution["method"] == rna_analysis.CDS_SOURCE_WHOLE
    assert resolution["n_codons"] == len(non_coding_sequence) // 3
    assert resolution["warnings"]


def test_resolve_coding_codons_prefers_annotation(single_orf_genome):
    annotation = [
        {
            "segments": [(0, 30)],
            "strand": "+",
            "codon_start": 1,
            "length_nt": 30,
            "product": "test product",
        }
    ]
    resolution = rna_analysis.resolve_coding_codons(
        single_orf_genome, cds_regions=annotation, source="auto"
    )
    assert resolution["method"] == rna_analysis.CDS_SOURCE_ANNOTATION
    assert resolution["n_cds"] == 1
    assert resolution["n_codons"] == 10


def test_resolve_coding_codons_reports_naive_reference(single_orf_genome):
    resolution = rna_analysis.resolve_coding_codons(
        single_orf_genome, source="prediction", min_cds_length=150
    )
    assert resolution["naive_codon_count"] == len(single_orf_genome) // 3
    assert resolution["method_label"] in rna_analysis.CDS_SOURCE_LABELS.values()


def test_resolve_coding_codons_validates_arguments(single_orf_genome):
    with pytest.raises(ValueError):
        rna_analysis.resolve_coding_codons("")
    with pytest.raises(ValueError):
        rna_analysis.resolve_coding_codons(single_orf_genome, source="magic")
    with pytest.raises(ValueError):
        rna_analysis.resolve_coding_codons(single_orf_genome, min_cds_length=0)


def test_rna_shannon_entropy_and_kmers():
    assert rna_analysis.shannon_entropy("AAAA") == pytest.approx(0.0)
    assert rna_analysis.shannon_entropy("ACGU") == pytest.approx(2.0)
    assert math.isnan(rna_analysis.shannon_entropy("NNNN"))
    assert math.isnan(rna_analysis.shannon_entropy(""))
    counts = rna_analysis.kmer_counts("AAA", k=2)
    assert counts["AA"]["count"] == 2
    with pytest.raises(ValueError):
        rna_analysis.kmer_counts("ACGU", k=0)


def test_rna_secondary_structure_availability_is_honest():
    availability = rna_analysis.secondary_structure_availability()
    assert availability["recommended_tools"]
    assert "ViennaRNA" in availability["reason"] or "RNAfold" in availability["reason"]
    if availability["available"]:
        assert availability["experimental"] is False
    else:
        assert availability["available"] is False
        assert (
            "unavailable" in availability["reason"].lower()
            or "not detected" in availability["reason"].lower()
            or "not installed" in availability["reason"].lower()
        )


def test_codon_adaptation_index_report_uses_embedded_table():
    report = rna_analysis.codon_adaptation_index_report("CUGCUGCUGAAA", organism="human")
    assert report["organism"] == "human"
    assert "Kazusa" in report["table_source"]
    assert report["status"] == "HEURISTIC"
    assert 0.0 <= report["value"] <= 1.0
    methionine_only = rna_analysis.codon_adaptation_index("AUG" * 4)
    assert math.isnan(methionine_only)
    empty_report = rna_analysis.codon_adaptation_index_report("AUG" * 4)
    assert empty_report["status"] == "UNAVAILABLE"
    assert empty_report["n_informative"] == 0


def test_effective_number_of_codons_report_unavailable_for_short():
    report = rna_analysis.effective_number_of_codons_report("AUGAUGAUG")
    assert report["status"] == "UNAVAILABLE"
    assert math.isnan(report["value"])
    biased = "AAA" * 5 + "GCU" * 5 + "CUG" * 5
    ok = rna_analysis.effective_number_of_codons_report(biased)
    assert ok["status"] == "COMPUTED"
    assert ok["value"] == rna_analysis.ENC_MIN


def test_rna_kmer_rejects_explosive_k():
    with pytest.raises(ValueError, match="memoria|memory"):
        rna_analysis.kmer_counts("ACGU", k=1_000_000)
