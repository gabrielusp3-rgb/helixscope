"""Testes da termodinamica nearest-neighbor SantaLucia & Hicks 2004.

Os valores esperados sao calculados a partir da Tabela 1 publicada
(Annu. Rev. Biophys. Biomol. Struct. 33:415-440, 2004), nao copiando a
funcao sob teste. Wallace/salt-adjusted permanecem metodos separados.
"""

from __future__ import annotations

import math

import pytest

from modules import dna_analysis, provenance, thermodynamics

# SantaLucia & Hicks 2004 Table 1 (kcal/mol, cal K^-1 mol^-1). Independent copy.
NN = {
    "init": (0.2, -5.7),
    "init_A/T": (2.2, 6.9),
    "init_G/C": (0.0, 0.0),
    "sym": (0.0, -1.4),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
}
R = 1.987
SALT_K = 0.368


def _hand_tm_cgcgcgcg(*, na_m: float, oligo_nm: float) -> float:
    """Tm de 5'-CGCGCGCG-3' (auto-complementar) a partir da tabela publicada."""
    # 4 x CG/GC + 3 x GC/CG
    dh = 4 * NN["CG/GC"][0] + 3 * NN["GC/CG"][0] + NN["init"][0] + NN["sym"][0]
    ds = 4 * NN["CG/GC"][1] + 3 * NN["GC/CG"][1] + NN["init"][1] + NN["sym"][1]
    # ends C and G: two G/C termini
    ds = ds + SALT_K * (8 - 1) * math.log(na_m)
    k = oligo_nm * 1e-9
    return (1000.0 * dh) / (ds + R * math.log(k)) - 273.15


def test_santalucia_self_complementary_matches_published_table():
    seq = "CGCGCGCG"
    expected = _hand_tm_cgcgcgcg(na_m=1.0, oligo_nm=250.0)
    result = thermodynamics.santalucia_nearest_neighbor_tm(
        seq, na_m=1.0, oligo_nm=250.0, mg_m=0.0
    )
    assert result["status"] == "COMPUTED"
    assert result["parameter_set"] == "SantaLucia_Hicks_2004_DNA_NN"
    assert result["self_complementary"] is True
    assert result["value_c"] == pytest.approx(round(expected, 2), abs=0.02)
    assert result["value_c"] != 0


def test_santalucia_salt_and_concentration_change_tm():
    seq = "ATGCATGCATGC"
    low_salt = thermodynamics.santalucia_nearest_neighbor_tm(seq, na_m=0.05, oligo_nm=250.0)
    high_salt = thermodynamics.santalucia_nearest_neighbor_tm(seq, na_m=1.0, oligo_nm=250.0)
    dilute = thermodynamics.santalucia_nearest_neighbor_tm(seq, na_m=0.05, oligo_nm=25.0)
    assert low_salt["status"] == high_salt["status"] == dilute["status"] == "COMPUTED"
    assert high_salt["value_c"] > low_salt["value_c"]
    assert dilute["cache_key"] != low_salt["cache_key"]
    assert dilute["value_c"] != low_salt["value_c"]


def test_santalucia_out_of_domain_is_na_not_zero():
    too_short = thermodynamics.santalucia_nearest_neighbor_tm("ATGCAT")
    assert too_short["status"] == "UNAVAILABLE"
    assert too_short["value_c"] is None
    too_long = thermodynamics.santalucia_nearest_neighbor_tm("A" * 61)
    assert too_long["status"] == "UNAVAILABLE"
    assert too_long["value_c"] is None
    ambiguous = thermodynamics.santalucia_nearest_neighbor_tm("ATGCATGCNTGCATGC")
    assert ambiguous["status"] == "UNAVAILABLE"
    assert ambiguous["value_c"] is None
    magnesium = thermodynamics.santalucia_nearest_neighbor_tm(
        "ATGCATGCATGC", mg_m=0.002
    )
    assert magnesium["status"] == "UNAVAILABLE"
    assert magnesium["value_c"] is None
    assert "Owczarzy" in magnesium["reason"]


def test_santalucia_does_not_replace_wallace():
    wallace = dna_analysis.melting_temperature_report("ATGC")
    nn = thermodynamics.santalucia_nearest_neighbor_tm("ATGC")
    assert wallace["method"] == "Wallace"
    assert wallace["status"] == "COMPUTED"
    assert nn["status"] == "UNAVAILABLE"
    salt = dna_analysis.melting_temperature_report("ACGTACGTACGTACGT")
    nn_long = thermodynamics.santalucia_nearest_neighbor_tm("ACGTACGTACGTACGT")
    assert salt["method"] == "salt-adjusted"
    assert nn_long["status"] == "COMPUTED"
    assert salt["value_c"] != nn_long["value_c"]


def test_santalucia_invalid_conditions_raise():
    with pytest.raises(thermodynamics.ThermodynamicsError):
        thermodynamics.santalucia_nearest_neighbor_tm("ATGCATGC", na_m=0)
    with pytest.raises(thermodynamics.ThermodynamicsError):
        thermodynamics.santalucia_nearest_neighbor_tm("ATGCATGC", oligo_nm=-1)


def test_santalucia_helpers_and_parameter_record():
    assert thermodynamics.reverse_complement_dna("ATGC") == "GCAT"
    assert thermodynamics.is_self_complementary("CGCGCGCG") is True
    assert thermodynamics.is_self_complementary("ATGCATGC") is False
    record = thermodynamics.nn_parameter_record()
    assert record["owczarzy"] is False
    assert record["mismatches"] is False
    key_a = thermodynamics.cache_key(
        sequence_hash="abc", na_m=0.05, oligo_nm=250.0, selfcomp=False
    )
    key_b = thermodynamics.cache_key(
        sequence_hash="abc", na_m=1.0, oligo_nm=250.0, selfcomp=False
    )
    assert key_a != key_b
    result = thermodynamics.santalucia_nearest_neighbor_tm("ATGCATGCATGC")
    assert result["sequence_hash"] == provenance.sequence_digest("ATGCATGCATGC")
    assert result["software_version"] == provenance.HELIXSCOPE_VERSION
