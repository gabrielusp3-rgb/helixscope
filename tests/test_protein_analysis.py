"""Testes unitarios do modulo modules.protein_analysis.

Cobrem a traducao com parada no codon stop, a normalizacao da composicao de
aminoacidos, o tratamento de janela maior que a proteina no perfil de
hidrofobicidade, as propriedades fisico-quimicas derivadas (GRAVY, carga liquida,
indice alifatico), o agrupamento de residuos por categoria e a declaracao de
indisponibilidade da predicao de estrutura secundaria. Nenhum teste faz chamadas
de rede.
"""

import pytest

from modules import protein_analysis


@pytest.fixture
def mrna_with_stop() -> str:
    """mRNA que codifica dois aminoacidos seguidos de um codon de parada."""
    return "AUGAAAUAA"


@pytest.fixture
def protein_sequence() -> str:
    """Sequencia curta de proteina reutilizavel nos testes."""
    return "ACDEFGHIKL"


def test_translate_stops_at_stop_codon(mrna_with_stop):
    assert protein_analysis.translate(mrna_with_stop) == "MK"


def test_amino_acid_composition_frequencies_sum_to_100(protein_sequence):
    composition = protein_analysis.amino_acid_composition(protein_sequence)
    total = sum(entry["frequency"] for entry in composition.values())
    assert total == pytest.approx(100.0, abs=0.01)


def test_hydrophobicity_profile_window_larger_than_protein():
    with pytest.raises(ValueError):
        protein_analysis.hydrophobicity_profile("AC", window=9)


def test_gravy_index_matches_scale_extremes():
    assert protein_analysis.gravy_index("IIII") == pytest.approx(4.5)
    assert protein_analysis.gravy_index("RRRR") == pytest.approx(-4.5)


def test_amino_acid_composition_rejects_nonstandard_residues():
    with pytest.raises(ValueError):
        protein_analysis.amino_acid_composition("MKVLX")
    with pytest.raises(ValueError):
        protein_analysis.gravy_index("AXA")


def test_aliphatic_index_known_values():
    assert protein_analysis.aliphatic_index("AAAA") == pytest.approx(100.0)
    assert protein_analysis.aliphatic_index("LLLL") == pytest.approx(390.0)
    assert protein_analysis.aliphatic_index("GGGG") == pytest.approx(0.0)


def test_net_charge_sign_follows_residue_type():
    assert protein_analysis.net_charge("KKKK", ph=7.0) > 0.0
    assert protein_analysis.net_charge("DDDD", ph=7.0) < 0.0


def test_net_charge_rejects_invalid_ph():
    with pytest.raises(ValueError):
        protein_analysis.net_charge("KKKK", ph=20.0)


def test_amino_acid_categories_exclusive_classes_sum_to_100(protein_sequence):
    categories = protein_analysis.amino_acid_categories(protein_sequence)
    exclusive = sum(
        categories[name]["frequency"]
        for name in protein_analysis.AMINO_ACID_CATEGORIES
    )
    assert exclusive == pytest.approx(100.0, abs=0.05)


def test_amino_acid_categories_counts_charged_residues():
    categories = protein_analysis.amino_acid_categories("KRHDE")
    assert categories["positively_charged"]["count"] == 3
    assert categories["negatively_charged"]["count"] == 2
    assert categories["charged_total"]["count"] == 5
    assert categories["charged_total"]["frequency"] == pytest.approx(100.0)


def test_amino_acid_categories_rejects_empty_protein():
    with pytest.raises(ValueError):
        protein_analysis.amino_acid_categories("123")


def test_secondary_structure_is_declared_unavailable():
    availability = protein_analysis.secondary_structure_availability()
    assert availability["available"] is False
    assert availability["reason"]
    assert availability["recommended_tools"]


def test_aromaticity_matches_protparam_definition():
    assert protein_analysis.aromaticity("FFFF") == pytest.approx(1.0)
    assert protein_analysis.aromaticity("AAAA") == pytest.approx(0.0)
    assert protein_analysis.aromaticity("AF") == pytest.approx(0.5)


def test_charge_count_profile_counts_formal_side_chains():
    profile = protein_analysis.charge_count_profile("KKKKEEEE", window=4)
    assert profile[0] == pytest.approx(1.0)
    assert profile[-1] == pytest.approx(-1.0)
    with pytest.raises(ValueError):
        protein_analysis.charge_count_profile("KK", window=9)


def test_composition_is_consistent_and_zero_count_is_not_na():
    seq = "ACDEFGHIKL"
    composition = protein_analysis.amino_acid_composition(seq)
    total = sum(int(item["count"]) for item in composition.values())
    assert total == len(seq)
    assert protein_analysis.composition_is_consistent(seq)
    assert composition["W"]["count"] == 0
    assert composition["A"]["count"] == 1


def test_physicochemical_report_declares_theoretical_methods():
    report = protein_analysis.physicochemical_report("ACDEFGHIKL", ph=7.0)
    assert report["status"] == "COMPUTED"
    assert report["isoelectric_point_kind"] == "theoretical"
    assert "Bjellqvist" in report["isoelectric_point_method"]
    assert report["gravy_method"] == "Kyte-Doolittle"
    assert report["net_charge_ph"] == pytest.approx(7.0)
    assert "Guruprasad" in report["instability_index_method"]
    assert "Ikai" in report["aliphatic_index_method"]
    records = protein_analysis.hydrophobicity_profile_records("ACDEFGHIKL", window=5)
    assert records[0]["start"] == 0
    assert records[0]["end"] == 5
    assert records[0]["method"] == "Kyte-Doolittle"
    with pytest.raises(ValueError):
        protein_analysis.hydrophobicity_profile_records("AC", window=9)
    availability = protein_analysis.feature_prediction_availability()
    assert availability["domains"]["available"] is False
    assert availability["signal_peptide"]["available"] is False
    assert availability["transmembrane"]["available"] is False
    assert availability["conservation"]["available"] is False


def test_per_residue_formal_charge_matches_window_criterion():
    values = protein_analysis.per_residue_formal_charge("KRHDE")
    assert values == [1.0, 1.0, 0.0, -1.0, -1.0]
    with pytest.raises(ValueError):
        protein_analysis.per_residue_formal_charge("AX")


def test_extinction_matches_protparam_and_seg_stays_unavailable():
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    sequence = "ACDEFGHIKL"
    reduced, cystine = ProteinAnalysis(sequence).molar_extinction_coefficient()
    report = protein_analysis.extinction_coefficients(sequence)
    assert report["reduced_cysteines"] == int(reduced)
    assert report["cystine_bridges"] == int(cystine)
    assert report["wavelength_nm"] == 280
    packed = protein_analysis.physicochemical_report(sequence)
    assert packed["extinction_coefficient_280_reduced"] == int(reduced)
    assert packed["seg_low_complexity"] == "UNAVAILABLE"
    assert "not the SEG" in packed["sequence_complexity_method"]
    assert protein_analysis.sequence_complexity_bits("AAAA") == 0.0
    curve = protein_analysis.charge_curve(sequence, ph_step=1.0)
    assert curve[0]["ph"] == 0.0
    assert curve[-1]["ph"] == 14.0
    assert len(curve) == 15

