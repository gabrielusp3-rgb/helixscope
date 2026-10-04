"""Termodinamica nearest-neighbor de oligos de DNA (SantaLucia).

Implementa o modelo published de duplex DNA/DNA perfeito com o conjunto de
parametros de SantaLucia e Hicks (2004) e a correcao de entropia de sodio de
SantaLucia (1998). Nao e Wallace, nao e a formula salt-adjusted empirica, e nao
mistura Owczarzy 2004/2008.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import math
from typing import Dict, Mapping, Optional, Tuple

from . import dna_analysis, provenance

GAS_CONSTANT_CAL: float = 1.987
"""Constante dos gases em cal K^-1 mol^-1 (SantaLucia 1998 / 2004)."""

SALT_ENTROPY_COEFFICIENT: float = 0.368
"""Coeficiente 0.368 x (N-1) x ln[Na+] sobre deltaS (SantaLucia 1998 PNAS)."""

MIN_OLIGO_NT: int = 8
"""Limite inferior do dominio de oligo para este parametro set calorimetrico."""

MAX_OLIGO_NT: int = 60
"""Limite superior declarado; o parametro set nao e um modelo de DNA genomico."""

DEFAULT_NA_M: float = 0.05
"""Sodio 50 mM (0.05 M)."""

DEFAULT_OLIGO_NM: float = 250.0
"""Concentracao total de oligo 250 nM (0.25 uM)."""

PARAMETER_SET: str = "SantaLucia_Hicks_2004_DNA_NN"
PARAMETER_CITATION: str = (
    "SantaLucia J Jr, Hicks D. Annu Rev Biophys Biomol Struct 33:415-440 (2004). "
    "Sodium entropy correction: SantaLucia J Jr. Proc Natl Acad Sci USA "
    "95:1460-1465 (1998)."
)

COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}

# deltaH kcal/mol, deltaS cal K^-1 mol^-1. SantaLucia & Hicks 2004 Table 1 /
# Biopython DNA_NN4 (same numbers). Keys are 5'-NN-3'/5'-complement-3'.
DNA_NN_SANTALUCIA_HICKS_2004: Dict[str, Tuple[float, float]] = {
    "init": (0.2, -5.7),
    "init_A/T": (2.2, 6.9),
    "init_G/C": (0.0, 0.0),
    "sym": (0.0, -1.4),
    "AA/TT": (-7.6, -21.3),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}


class ThermodynamicsError(Exception):
    """Entrada invalida para o calculo nearest-neighbor.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def nn_parameter_record() -> dict:
    """Metadados do conjunto de parametros realmente usado.

    Args:
        Nenhum.

    Returns:
        Dict name, citation, salt_correction, mismatches, mg.

    Raises:
        Nenhum.
    """
    return {
        "name": PARAMETER_SET,
        "citation": PARAMETER_CITATION,
        "salt_correction": (
            "SantaLucia 1998 entropy: deltaS += 0.368 * (N-1) * ln([Na+] / 1 M)"
        ),
        "mismatches": False,
        "dangling_ends": False,
        "mg_correction": False,
        "owczarzy": False,
        "wallace": False,
        "salt_adjusted_empirical": False,
    }


def reverse_complement_dna(seq: str) -> str:
    """Complemento reverso canonico A/C/G/T.

    Args:
        seq: DNA canonico.

    Returns:
        Complemento reverso.

    Raises:
        Nenhum.
    """
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def is_self_complementary(seq: str) -> bool:
    """True se a sequencia for identica ao seu complemento reverso.

    Args:
        seq: DNA canonico.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    return seq == reverse_complement_dna(seq)


def cache_key(
    *,
    sequence_hash: str,
    na_m: float,
    oligo_nm: float,
    selfcomp: bool,
) -> str:
    """Chave de cache do Tm nearest-neighbor.

    Args:
        sequence_hash: SHA-256 da sequencia canonica.
        na_m: Sodio em M.
        oligo_nm: Concentracao total de oligo em nM.
        selfcomp: Duplex auto-complementar.

    Returns:
        Hex SHA-256.

    Raises:
        Nenhum.
    """
    import hashlib

    payload = "|".join(
        [
            str(sequence_hash),
            PARAMETER_SET,
            f"{float(na_m):.12g}",
            f"{float(oligo_nm):.12g}",
            "self" if selfcomp else "nonself",
            provenance.HELIXSCOPE_VERSION,
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _lookup_dimer(dimer: str) -> Tuple[float, float]:
    complement = "".join(COMPLEMENT[base] for base in dimer)
    key = f"{dimer}/{complement}"
    if key in DNA_NN_SANTALUCIA_HICKS_2004:
        return DNA_NN_SANTALUCIA_HICKS_2004[key]
    reversed_key = key[::-1]
    if reversed_key in DNA_NN_SANTALUCIA_HICKS_2004:
        return DNA_NN_SANTALUCIA_HICKS_2004[reversed_key]
    raise ThermodynamicsError(
        f"No SantaLucia & Hicks 2004 nearest-neighbor pair for {key}.",
        "INVALID_INPUT",
    )


def _canonical_dna(seq: str) -> str:
    cleaned = dna_analysis._clean(seq)
    if any(base in cleaned for base in "U"):
        raise ThermodynamicsError(
            "SantaLucia DNA/DNA nearest-neighbor does not accept RNA (U).",
            "INVALID_INPUT",
        )
    return cleaned


def santalucia_nearest_neighbor_tm(
    seq: str,
    *,
    na_m: float = DEFAULT_NA_M,
    oligo_nm: float = DEFAULT_OLIGO_NM,
    mg_m: float = 0.0,
) -> dict:
    """Tm nearest-neighbor DNA/DNA (SantaLucia & Hicks 2004 + sal 1998).

    Apenas duplex perfeito Watson-Crick. Sem mismatches, dangling ends,
    Owczarzy nem Mg2+. Fora do dominio o valor e N/A (None), nunca 0.

    Args:
        seq: Oligo de DNA (sera normalizado).
        na_m: Concentracao de sodio em mol/L. Deve ser positiva.
        oligo_nm: Concentracao total de fita em nanomolar. Deve ser positiva.
        mg_m: Magnesio em mol/L. Qualquer valor > 0 devolve N/A: a correcao
            de Mg nao faz parte deste parametro set.

    Returns:
        Envelope com value_c (float ou None), status COMPUTED ou UNAVAILABLE,
        delta_h_kcal, delta_s_cal, method, parameter_set, conditions e reason.

    Raises:
        ThermodynamicsError: INVALID_INPUT se Na, oligo ou Mg forem nao finitos
            ou se Na/oligo nao forem positivos.

    Nota biologica:
        Tm de oligo em solucao nao e Tm genomica nem Tm de qPCR com polimerase.
        Implementacao do modelo nearest-neighbor publicado (parametros
        SantaLucia & Hicks 2004; correcao de entropia SantaLucia 1998).
        Nao reproduz calculadoras comerciais (IDT, Primer3 defaults) se as
        concentracoes ou a correcao de sal forem diferentes.
    """
    if not isinstance(na_m, (int, float)) or not math.isfinite(float(na_m)) or float(na_m) <= 0:
        raise ThermodynamicsError("na_m must be a finite positive molar concentration.")
    if (
        not isinstance(oligo_nm, (int, float))
        or not math.isfinite(float(oligo_nm))
        or float(oligo_nm) <= 0
    ):
        raise ThermodynamicsError("oligo_nm must be a finite positive nanomolar concentration.")
    if not isinstance(mg_m, (int, float)) or not math.isfinite(float(mg_m)) or float(mg_m) < 0:
        raise ThermodynamicsError("mg_m must be a finite concentration >= 0.")

    na = float(na_m)
    oligo = float(oligo_nm)
    mg = float(mg_m)
    cleaned = _canonical_dna(seq)
    canonical = "".join(base for base in cleaned if base in "ACGT")
    n = len(canonical)
    seq_hash = provenance.sequence_digest(canonical)
    selfcomp = is_self_complementary(canonical) if n >= 2 else False
    key = cache_key(
        sequence_hash=seq_hash,
        na_m=na,
        oligo_nm=oligo,
        selfcomp=selfcomp,
    )
    base_payload = {
        "value_c": None,
        "status": "UNAVAILABLE",
        "method": "SantaLucia & Hicks 2004 nearest-neighbor DNA/DNA",
        "parameter_set": PARAMETER_SET,
        "citation": PARAMETER_CITATION,
        "n_canonical": n,
        "sequence_hash": seq_hash,
        "self_complementary": selfcomp,
        "na_m": na,
        "mg_m": mg,
        "oligo_nm": oligo,
        "delta_h_kcal_mol": None,
        "delta_s_cal_k_mol": None,
        "salt_correction": "SantaLucia 1998 entropy",
        "cache_key": key,
        "reason": "",
    }
    if mg > 0:
        base_payload["reason"] = (
            "Mg2+ is outside the SantaLucia 1998 sodium-only entropy correction. "
            "Owczarzy 2008 is not mixed in. Tm is N/A, not 0."
        )
        return _envelope(base_payload)
    if n == 0:
        base_payload["reason"] = (
            "No canonical A/C/G/T bases. SantaLucia nearest-neighbor Tm is N/A."
        )
        return _envelope(base_payload)
    if n < MIN_OLIGO_NT or n > MAX_OLIGO_NT:
        base_payload["reason"] = (
            f"SantaLucia nearest-neighbor Tm applies to DNA oligos of "
            f"{MIN_OLIGO_NT}-{MAX_OLIGO_NT} canonical bases. Length {n} is N/A."
        )
        return _envelope(base_payload)
    if any(base not in "ACGT" for base in cleaned):
        base_payload["reason"] = (
            "Ambiguous or non-canonical bases are not padded or averaged. Tm is N/A."
        )
        return _envelope(base_payload)

    delta_h, delta_s = _duplex_thermo(canonical, selfcomp=selfcomp)
    delta_s = delta_s + SALT_ENTROPY_COEFFICIENT * (n - 1) * math.log(na)
    if selfcomp:
        k_molar = oligo * 1e-9
    else:
        dnac1 = oligo / 2.0
        dnac2 = oligo / 2.0
        k_molar = (dnac1 - (dnac2 / 2.0)) * 1e-9
    if k_molar <= 0:
        base_payload["reason"] = "Strand concentration produced a non-positive k. Tm is N/A."
        return _envelope(base_payload)
    denom = delta_s + GAS_CONSTANT_CAL * math.log(k_molar)
    if denom == 0:
        base_payload["reason"] = "Nearest-neighbor denominator is zero. Tm is N/A."
        return _envelope(base_payload)
    tm_c = (1000.0 * delta_h) / denom - 273.15
    if not math.isfinite(tm_c):
        base_payload["reason"] = "Nearest-neighbor Tm is not finite. Value is N/A, not 0."
        return _envelope(base_payload)
    base_payload.update(
        {
            "value_c": round(tm_c, 2),
            "status": "COMPUTED",
            "delta_h_kcal_mol": round(delta_h, 4),
            "delta_s_cal_k_mol": round(delta_s, 4),
            "k_molar": k_molar,
            "reason": "",
        }
    )
    return _envelope(base_payload)


def _duplex_thermo(seq: str, *, selfcomp: bool) -> Tuple[float, float]:
    init_h, init_s = DNA_NN_SANTALUCIA_HICKS_2004["init"]
    delta_h = init_h
    delta_s = init_s
    ends = seq[0] + seq[-1]
    at_ends = ends.count("A") + ends.count("T")
    gc_ends = ends.count("G") + ends.count("C")
    at_h, at_s = DNA_NN_SANTALUCIA_HICKS_2004["init_A/T"]
    gc_h, gc_s = DNA_NN_SANTALUCIA_HICKS_2004["init_G/C"]
    delta_h += at_h * at_ends + gc_h * gc_ends
    delta_s += at_s * at_ends + gc_s * gc_ends
    if selfcomp:
        sym_h, sym_s = DNA_NN_SANTALUCIA_HICKS_2004["sym"]
        delta_h += sym_h
        delta_s += sym_s
    for index in range(len(seq) - 1):
        pair_h, pair_s = _lookup_dimer(seq[index : index + 2])
        delta_h += pair_h
        delta_s += pair_s
    return delta_h, delta_s


def _envelope(payload: Mapping[str, object]) -> dict:
    status = str(payload.get("status") or "UNAVAILABLE")
    record = provenance.analysis_envelope(
        module="DNA thermodynamics",
        payload=dict(payload),
        status=status,
        algorithm=str(payload.get("method") or PARAMETER_SET),
        parameters={
            "parameter_set": PARAMETER_SET,
            "na_m": payload.get("na_m"),
            "mg_m": payload.get("mg_m"),
            "oligo_nm": payload.get("oligo_nm"),
            "n_canonical": payload.get("n_canonical"),
            "self_complementary": payload.get("self_complementary"),
        },
        source="HelixScope SantaLucia nearest-neighbor",
        input_identifier=str(payload.get("sequence_hash") or ""),
    )
    record.update(payload)
    record["status"] = status
    return record
