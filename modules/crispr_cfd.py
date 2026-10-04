"""Cutting Frequency Determination (CFD) scores for SpCas9 off-targets.

Implements the published Doench et al. 2016 per-site CFD algorithm as used by
CRISPOR (`calcCfdScore` / `cfd-score-calculator.py`). Mismatch and PAM weights
are the published matrices, embedded as constants (GuideMaker USDA dump of the
same CRISPOR/Doench pickle tables). No pickle is loaded. No invented matrix.

Honestidade cientifica: este modulo NAO e Rule Set 2 nem Azimuth. CFD descreve
a frequencia relativa de corte num sitio off-target (escala 0-1). Agregacao
ao nivel do guia e feita em crispr_specificity.py e so e valida para o escopo
da busca realmente executada.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Mapping

from . import provenance

CFD_MODEL: str = "Doench2016_CFD"
CFD_VERSION: str = "CRISPOR-calcCfdScore-Doench2016-matrices"
CFD_METHOD: str = (
    "Doench et al. 2016 Cutting Frequency Determination; CRISPOR calcCfdScore "
    "(Haeussler et al. 2016). Per-site scale 0-1. Not Rule Set 2. Not Azimuth."
)
CFD_SOURCE: str = (
    "Doench JG et al. Nature Biotechnology 34:184-191 (2016); matrices as "
    "distributed with CRISPOR cfd-score-calculator.py / GuideMaker "
    "guidemaker/data/cfd_data.json (USDA-ARS-GBRU). Numeric tables embedded; "
    "no pickle deserialization."
)
CFD_SCALE: str = "0-1"
CFD_SPACER_NT: int = 20
CFD_CONTEXT_NT: int = 23

_REVCOM: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
}

CFD_MM_SCORES: dict[str, float] = {"rA:dA,1": 1.0, "rA:dA,10": 0.882352941, "rA:dA,11": 0.307692308, "rA:dA,12": 0.333333333, "rA:dA,13": 0.3, "rA:dA,14": 0.533333333, "rA:dA,15": 0.2, "rA:dA,16": 0.0, "rA:dA,17": 0.133333333, "rA:dA,18": 0.5, "rA:dA,19": 0.538461538, "rA:dA,2": 0.727272727, "rA:dA,20": 0.6, "rA:dA,3": 0.705882353, "rA:dA,4": 0.636363636, "rA:dA,5": 0.363636364, "rA:dA,6": 0.7142857140000001, "rA:dA,7": 0.4375, "rA:dA,8": 0.428571429, "rA:dA,9": 0.6, "rA:dC,1": 1.0, "rA:dC,10": 0.5555555560000001, "rA:dC,11": 0.65, "rA:dC,12": 0.7222222220000001, "rA:dC,13": 0.6521739129999999, "rA:dC,14": 0.46666666700000003, "rA:dC,15": 0.65, "rA:dC,16": 0.192307692, "rA:dC,17": 0.176470588, "rA:dC,18": 0.4, "rA:dC,19": 0.375, "rA:dC,2": 0.8, "rA:dC,20": 0.764705882, "rA:dC,3": 0.611111111, "rA:dC,4": 0.625, "rA:dC,5": 0.72, "rA:dC,6": 0.7142857140000001, "rA:dC,7": 0.705882353, "rA:dC,8": 0.7333333329999999, "rA:dC,9": 0.666666667, "rA:dG,1": 0.857142857, "rA:dG,10": 0.333333333, "rA:dG,11": 0.4, "rA:dG,12": 0.263157895, "rA:dG,13": 0.21052631600000002, "rA:dG,14": 0.214285714, "rA:dG,15": 0.272727273, "rA:dG,16": 0.0, "rA:dG,17": 0.176470588, "rA:dG,18": 0.19047619, "rA:dG,19": 0.20689655199999998, "rA:dG,2": 0.7857142859999999, "rA:dG,20": 0.22727272699999998, "rA:dG,3": 0.428571429, "rA:dG,4": 0.352941176, "rA:dG,5": 0.5, "rA:dG,6": 0.454545455, "rA:dG,7": 0.4375, "rA:dG,8": 0.428571429, "rA:dG,9": 0.571428571, "rC:dA,1": 1.0, "rC:dA,10": 0.9411764709999999, "rC:dA,11": 0.307692308, "rC:dA,12": 0.538461538, "rC:dA,13": 0.7, "rC:dA,14": 0.7333333329999999, "rC:dA,15": 0.066666667, "rC:dA,16": 0.307692308, "rC:dA,17": 0.46666666700000003, "rC:dA,18": 0.642857143, "rC:dA,19": 0.46153846200000004, "rC:dA,2": 0.9090909090000001, "rC:dA,20": 0.3, "rC:dA,3": 0.6875, "rC:dA,4": 0.8, "rC:dA,5": 0.636363636, "rC:dA,6": 0.9285714290000001, "rC:dA,7": 0.8125, "rC:dA,8": 0.875, "rC:dA,9": 0.875, "rC:dC,1": 0.913043478, "rC:dC,10": 0.38888888899999996, "rC:dC,11": 0.25, "rC:dC,12": 0.444444444, "rC:dC,13": 0.13636363599999998, "rC:dC,14": 0.0, "rC:dC,15": 0.05, "rC:dC,16": 0.153846154, "rC:dC,17": 0.058823529000000006, "rC:dC,18": 0.133333333, "rC:dC,19": 0.125, "rC:dC,2": 0.695652174, "rC:dC,20": 0.058823529000000006, "rC:dC,3": 0.5, "rC:dC,4": 0.5, "rC:dC,5": 0.6, "rC:dC,6": 0.5, "rC:dC,7": 0.470588235, "rC:dC,8": 0.642857143, "rC:dC,9": 0.6190476189999999, "rC:dT,1": 1.0, "rC:dT,10": 0.8666666670000001, "rC:dT,11": 0.75, "rC:dT,12": 0.7142857140000001, "rC:dT,13": 0.384615385, "rC:dT,14": 0.35, "rC:dT,15": 0.222222222, "rC:dT,16": 1.0, "rC:dT,17": 0.46666666700000003, "rC:dT,18": 0.538461538, "rC:dT,19": 0.428571429, "rC:dT,2": 0.727272727, "rC:dT,20": 0.5, "rC:dT,3": 0.8666666670000001, "rC:dT,4": 0.842105263, "rC:dT,5": 0.571428571, "rC:dT,6": 0.9285714290000001, "rC:dT,7": 0.75, "rC:dT,8": 0.65, "rC:dT,9": 0.857142857, "rG:dA,1": 1.0, "rG:dA,10": 0.8125, "rG:dA,11": 0.384615385, "rG:dA,12": 0.384615385, "rG:dA,13": 0.3, "rG:dA,14": 0.26666666699999997, "rG:dA,15": 0.14285714300000002, "rG:dA,16": 0.0, "rG:dA,17": 0.25, "rG:dA,18": 0.666666667, "rG:dA,19": 0.666666667, "rG:dA,2": 0.636363636, "rG:dA,20": 0.7, "rG:dA,3": 0.5, "rG:dA,4": 0.363636364, "rG:dA,5": 0.3, "rG:dA,6": 0.666666667, "rG:dA,7": 0.571428571, "rG:dA,8": 0.625, "rG:dA,9": 0.533333333, "rG:dG,1": 0.7142857140000001, "rG:dG,10": 0.4, "rG:dG,11": 0.428571429, "rG:dG,12": 0.529411765, "rG:dG,13": 0.42105263200000004, "rG:dG,14": 0.428571429, "rG:dG,15": 0.272727273, "rG:dG,16": 0.0, "rG:dG,17": 0.235294118, "rG:dG,18": 0.47619047600000003, "rG:dG,19": 0.448275862, "rG:dG,2": 0.692307692, "rG:dG,20": 0.428571429, "rG:dG,3": 0.384615385, "rG:dG,4": 0.529411765, "rG:dG,5": 0.7857142859999999, "rG:dG,6": 0.681818182, "rG:dG,7": 0.6875, "rG:dG,8": 0.615384615, "rG:dG,9": 0.538461538, "rG:dT,1": 0.9, "rG:dT,10": 0.933333333, "rG:dT,11": 1.0, "rG:dT,12": 0.933333333, "rG:dT,13": 0.923076923, "rG:dT,14": 0.75, "rG:dT,15": 0.9411764709999999, "rG:dT,16": 1.0, "rG:dT,17": 0.933333333, "rG:dT,18": 0.692307692, "rG:dT,19": 0.7142857140000001, "rG:dT,2": 0.846153846, "rG:dT,20": 0.9375, "rG:dT,3": 0.75, "rG:dT,4": 0.9, "rG:dT,5": 0.8666666670000001, "rG:dT,6": 1.0, "rG:dT,7": 1.0, "rG:dT,8": 1.0, "rG:dT,9": 0.642857143, "rU:dC,1": 0.956521739, "rU:dC,10": 0.5, "rU:dC,11": 0.4, "rU:dC,12": 0.5, "rU:dC,13": 0.260869565, "rU:dC,14": 0.0, "rU:dC,15": 0.05, "rU:dC,16": 0.346153846, "rU:dC,17": 0.117647059, "rU:dC,18": 0.333333333, "rU:dC,19": 0.25, "rU:dC,2": 0.84, "rU:dC,20": 0.176470588, "rU:dC,3": 0.5, "rU:dC,4": 0.625, "rU:dC,5": 0.64, "rU:dC,6": 0.571428571, "rU:dC,7": 0.588235294, "rU:dC,8": 0.7333333329999999, "rU:dC,9": 0.6190476189999999, "rU:dG,1": 0.857142857, "rU:dG,10": 0.533333333, "rU:dG,11": 0.666666667, "rU:dG,12": 0.947368421, "rU:dG,13": 0.7894736840000001, "rU:dG,14": 0.28571428600000004, "rU:dG,15": 0.272727273, "rU:dG,16": 0.666666667, "rU:dG,17": 0.705882353, "rU:dG,18": 0.428571429, "rU:dG,19": 0.275862069, "rU:dG,2": 0.857142857, "rU:dG,20": 0.090909091, "rU:dG,3": 0.428571429, "rU:dG,4": 0.647058824, "rU:dG,5": 1.0, "rU:dG,6": 0.9090909090000001, "rU:dG,7": 0.6875, "rU:dG,8": 1.0, "rU:dG,9": 0.923076923, "rU:dT,1": 1.0, "rU:dT,10": 0.857142857, "rU:dT,11": 0.75, "rU:dT,12": 0.8, "rU:dT,13": 0.692307692, "rU:dT,14": 0.6190476189999999, "rU:dT,15": 0.578947368, "rU:dT,16": 0.9090909090000001, "rU:dT,17": 0.533333333, "rU:dT,18": 0.666666667, "rU:dT,19": 0.28571428600000004, "rU:dT,2": 0.846153846, "rU:dT,20": 0.5625, "rU:dT,3": 0.7142857140000001, "rU:dT,4": 0.47619047600000003, "rU:dT,5": 0.5, "rU:dT,6": 0.8666666670000001, "rU:dT,7": 0.875, "rU:dT,8": 0.8, "rU:dT,9": 0.9285714290000001}
"""Pesos mismatch tipo x posicao (Doench 2016 / CRISPOR)."""

CFD_PAM_SCORES: dict[str, float] = {"AA": 0.0, "AC": 0.0, "AG": 0.25925925899999996, "AT": 0.0, "CA": 0.0, "CC": 0.0, "CG": 0.107142857, "CT": 0.0, "GA": 0.06944444400000001, "GC": 0.022222222000000003, "GG": 1.0, "GT": 0.016129031999999998, "TA": 0.0, "TC": 0.0, "TG": 0.038961038999999996, "TT": 0.0}
"""Pesos dos dinucleotidos do PAM (GG=1.0, AG ~ 0.259)."""


def _revcom_base(base: str) -> str:
    """Complemento de uma base para a chave CFD (CRISPOR revcom)."""
    return _REVCOM.get(base, "")


def cfd_model_record() -> dict:
    """Metadados do modelo CFD embutido.

    Args:
        Nenhum.

    Returns:
        Dict com model, method, version, source, scale, status e n_weights.

    Raises:
        Nenhum.
    """
    return {
        "model": CFD_MODEL,
        "method": CFD_METHOD,
        "version": CFD_VERSION,
        "source": CFD_SOURCE,
        "scale": CFD_SCALE,
        "status": "COMPUTED" if CFD_MM_SCORES and CFD_PAM_SCORES else "UNAVAILABLE",
        "n_mismatch_weights": len(CFD_MM_SCORES),
        "n_pam_weights": len(CFD_PAM_SCORES),
        "spacer_nt": CFD_SPACER_NT,
        "context_nt": CFD_CONTEXT_NT,
        "bulges": False,
        "pickle_used": False,
    }


def cfd_pair_score(guide_23: str, off_target_23: str) -> dict:
    """CFD entre um par 23-mer (20 nt spacer + 3 nt PAM).

    Implementacao alinhada a CRISPOR calcCfdScore: PAM = dois nucleotidos
    finais do off-target; mismatches no spacer de 20 nt com chave
    r{wt}:d{revcom(off)},{pos 1-20}; produto dos pesos; vezes peso do PAM.

    Args:
        guide_23: On-target 23-mer DNA (spacer + PAM), A/C/G/T.
        off_target_23: Off-target 23-mer DNA (spacer + PAM), A/C/G/T.

    Returns:
        Dict com score (float|None), status, method, model, version, scale.

    Raises:
        Nenhum.

    Nota biologica:
        CFD 1.0 e o on-target com PAM GG. NAG usa pam_scores['AG'] ~ 0.259.
        Nao e especificidade genoma-wide. Sem bulges. Sequencias com letras
        nao-ATCG devolvem INVALID_INPUT, nunca 0.0 fabricado. Comprimento
        diferente de 23 devolve UNAVAILABLE.
    """
    record = {
        "score": None,
        "status": "UNAVAILABLE",
        "method": CFD_METHOD,
        "model": CFD_MODEL,
        "version": CFD_VERSION,
        "scale": CFD_SCALE,
        "source": CFD_SOURCE,
        "parameters": {
            "spacer_nt": CFD_SPACER_NT,
            "context_nt": CFD_CONTEXT_NT,
            "bulges": False,
            "pickle_used": False,
        },
        "reason": "",
        "pam_dinucleotide": "",
        "mismatch_keys": [],
    }
    if not CFD_MM_SCORES or not CFD_PAM_SCORES:
        record["reason"] = "CFD matrices are not embedded."
        return record
    wt = str(guide_23 or "").strip().upper().replace("U", "T")
    off = str(off_target_23 or "").strip().upper().replace("U", "T")
    if len(wt) != CFD_CONTEXT_NT or len(off) != CFD_CONTEXT_NT:
        record["reason"] = (
            "CFD is defined here for 23-mer spacer+PAM pairs only "
            f"(got {len(wt)} and {len(off)} nt)."
        )
        return record
    allowed = set("ACGT")
    if set(wt) - allowed or set(off) - allowed:
        record["status"] = "INVALID_INPUT"
        record["reason"] = "CFD 23-mers must be canonical DNA A/C/G/T."
        return record
    pam = off[-2:]
    pam_weight = CFD_PAM_SCORES.get(pam)
    if pam_weight is None:
        record["status"] = "INVALID_INPUT"
        record["reason"] = f"No published CFD PAM weight for dinucleotide {pam}."
        return record
    wt_rna = wt.replace("T", "U")
    off_rna = off.replace("T", "U")
    wt_spacer = wt_rna[:CFD_SPACER_NT]
    off_spacer = off_rna[:CFD_SPACER_NT]
    score = 1.0
    keys: list[str] = []
    for index, rna_base in enumerate(wt_spacer):
        dna_off = off_spacer[index]
        if rna_base == dna_off:
            continue
        rev = _revcom_base(dna_off)
        key = f"r{rna_base}:d{rev},{index + 1}"
        weight = CFD_MM_SCORES.get(key)
        if weight is None:
            record["status"] = "UNAVAILABLE"
            record["reason"] = f"Missing published CFD weight for key {key}."
            record["score"] = None
            return record
        score *= float(weight)
        keys.append(key)
    score *= float(pam_weight)
    record["score"] = float(score)
    record["status"] = "COMPUTED"
    record["pam_dinucleotide"] = pam
    record["pam_weight"] = float(pam_weight)
    record["mismatch_keys"] = keys
    record["reason"] = ""
    record["software_version"] = provenance.HELIXSCOPE_VERSION
    return record


def cfd_from_spacer_and_pam(
    guide_spacer: str,
    guide_pam: str,
    off_spacer: str,
    off_pam: str,
) -> dict:
    """CFD a partir de spacers de 20 nt e PAMs de 3 nt.

    Args:
        guide_spacer: Spacer on-target (20 nt).
        guide_pam: PAM on-target (3 nt, tipicamente NGG).
        off_spacer: Spacer do sitio candidato (20 nt).
        off_pam: PAM extraido do contig (3 nt).

    Returns:
        Saida de cfd_pair_score, ou UNAVAILABLE/INVALID_INPUT.

    Raises:
        Nenhum.
    """
    wt = (str(guide_spacer or "").strip().upper().replace("U", "T")
          + str(guide_pam or "").strip().upper().replace("U", "T"))
    off = (str(off_spacer or "").strip().upper().replace("U", "T")
           + str(off_pam or "").strip().upper().replace("U", "T"))
    return cfd_pair_score(wt, off)


def attach_cfd_to_hit(hit: Mapping[str, object], guide_pam: str) -> dict:
    """Acrescenta campos CFD a um hit ja verificado, sem alterar Hsu.

    Args:
        hit: Hit com guide_sequence, target_sequence e PAM.
        guide_pam: PAM NGG do guia on-target.

    Returns:
        Novo dict do hit com cfd_score e metadados. Hsu permanece em score.

    Raises:
        Nenhum.
    """
    enriched = dict(hit)
    result = cfd_from_spacer_and_pam(
        str(hit.get("guide_sequence") or ""),
        str(guide_pam or ""),
        str(hit.get("target_sequence") or ""),
        str(hit.get("PAM") or hit.get("pam_sequence") or ""),
    )
    enriched["cfd_score"] = result.get("score")
    enriched["cfd_status"] = result.get("status")
    enriched["cfd_method"] = result.get("method")
    enriched["cfd_model"] = result.get("model")
    enriched["cfd_version"] = result.get("version")
    enriched["cfd_scale"] = result.get("scale")
    enriched["cfd_reason"] = result.get("reason") or ""
    enriched["cfd_pam_dinucleotide"] = result.get("pam_dinucleotide") or ""
    return enriched
