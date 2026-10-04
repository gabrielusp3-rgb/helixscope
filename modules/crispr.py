"""Desenho de guias CRISPR, varredura de off-targets e primers de validacao.

Reune utilidades educacionais para projetar RNAs guia (gRNAs) para diferentes
sistemas Cas, pontuar a eficiencia esperada com uma heuristica posicional,
varrer off-targets por janela deslizante e desenhar pares de primers para
validacao por sequenciamento.

A pontuacao de eficiencia (score_guide_doench) e inspirada no espirito do
Rule Set 2 de Doench et al. (2016), "Optimized sgRNA design to maximize activity
and minimize off-target effects of CRISPR-Cas9", Nature Biotechnology 34:184-191.

Honestidade cientifica: este modulo NAO reproduz o modelo publicado de Doench et
al. (2016) nem qualquer modelo treinado em dados experimentais. As constantes
DOENCH_SCORES sao pesos posicionais puramente ilustrativos, definidos a mao para
fins didaticos. Da mesma forma, a varredura de off-targets e uma heuristica
simplificada e nao substitui ferramentas validadas como Cas-OFFinder ou CRISPOR
em uso laboratorial real. Os resultados servem apenas para exploracao e ensino.
"""

from __future__ import annotations

import hashlib
import math
from typing import Any, Callable, Dict, List, Mapping, Optional

import pandas as pd

from . import provenance

PAM_PATTERNS: Dict[str, str] = {
    "SpCas9": "NGG",       # downstream, 3' do protoespacador
    "SaCas9": "NNGRRT",    # downstream
    "AsCas12a": "TTTV",    # upstream, 5' do protoespacador
    "CjCas9": "NNNNRYAC",  # downstream
    "Cas12b": "TTN",       # upstream, 5' do protoespacador
}
"""Padroes de reconhecimento de PAM (codigo IUPAC) por nuclease que reconhece DNA
por PAM. Cas13 nao aparece aqui porque reconhece RNA e nao depende de PAM."""

PAM_ORIENTATION: Dict[str, str] = {
    "SpCas9": "downstream",
    "SaCas9": "downstream",
    "AsCas12a": "upstream",
    "CjCas9": "downstream",
    "Cas12b": "upstream",
}
"""Posicao do PAM em relacao ao protoespacador (downstream=3', upstream=5')."""

GUIDE_LENGTH_BY_CAS: Dict[str, int] = {
    "SpCas9": 20,
    "SaCas9": 21,
    "AsCas12a": 23,
    "CjCas9": 22,
    "Cas12b": 20,
    "Cas13": 28,
}
"""Comprimento do protoespacador (ou espacador, no caso do Cas13) por sistema."""

CAS_SYSTEMS: Dict[str, dict] = {
    "SpCas9": {
        "label": "SpCas9",
        "pam_system": "SpCas9",
        "pam_display": "NGG (3')",
        "target_molecule": "DNA",
        "editor": "nuclease",
        "cut_type": "blunt",
        "cut_offset_from_pam": 3,
        "notes": "Corte cego cerca de 3 pb a montante do PAM NGG.",
    },
    "SaCas9": {
        "label": "SaCas9",
        "pam_system": "SaCas9",
        "pam_display": "NNGRRT (3')",
        "target_molecule": "DNA",
        "editor": "nuclease",
        "cut_type": "blunt",
        "cut_offset_from_pam": 3,
        "notes": "Nuclease compacta; PAM mais longo reduz a densidade de sitios.",
    },
    "AsCas12a": {
        "label": "AsCas12a (Cpf1)",
        "pam_system": "AsCas12a",
        "pam_display": "TTTV (5')",
        "target_molecule": "DNA",
        "editor": "nuclease",
        "cut_type": "staggered",
        "cut_offset_from_pam": 18,
        "notes": "Corte escalonado distal ao PAM, gerando extremidades coesivas.",
    },
    "Cas12b": {
        "label": "Cas12b (C2c1)",
        "pam_system": "Cas12b",
        "pam_display": "TTN (5')",
        "target_molecule": "DNA",
        "editor": "nuclease",
        "cut_type": "staggered",
        "cut_offset_from_pam": 18,
        "notes": "Nuclease termoestavel; corte escalonado distal ao PAM.",
    },
    "Cas13": {
        "label": "Cas13 (RNA)",
        "pam_system": None,
        "pam_display": "sem PAM; PFS 3' nao-G (LwaCas13a)",
        "target_molecule": "RNA",
        "editor": "nuclease_rna",
        "cut_type": "collateral",
        "cut_offset_from_pam": None,
        "notes": (
            "Alveja RNA de fita simples; nao produz corte sitio-especifico no DNA "
            "e exibe clivagem colateral de RNAs vizinhos apos ativacao."
        ),
    },
    "CBE (SpCas9)": {
        "label": "Cytosine Base Editor (SpCas9)",
        "pam_system": "SpCas9",
        "pam_display": "NGG (3')",
        "target_molecule": "DNA",
        "editor": "base_editor",
        "edit_from": "C",
        "edit_to": "T",
        "edit_window": (4, 8),
        "cut_type": "nick",
        "cut_offset_from_pam": 3,
        "notes": (
            "Deamina citosinas para uracila na janela de edicao, resultando em "
            "conversao C para T; nao gera quebra de fita dupla."
        ),
    },
    "ABE (SpCas9)": {
        "label": "Adenine Base Editor (SpCas9)",
        "pam_system": "SpCas9",
        "pam_display": "NGG (3')",
        "target_molecule": "DNA",
        "editor": "base_editor",
        "edit_from": "A",
        "edit_to": "G",
        "edit_window": (4, 8),
        "cut_type": "nick",
        "cut_offset_from_pam": 3,
        "notes": (
            "Deamina adeninas na janela de edicao, resultando em conversao A para "
            "G; nao gera quebra de fita dupla."
        ),
    },
    "Prime Editor (SpCas9 nickase)": {
        "label": "Prime Editor (SpCas9 H840A nickase)",
        "pam_system": "SpCas9",
        "pam_display": "NGG (3')",
        "target_molecule": "DNA",
        "editor": "prime_editor",
        "cut_type": "nick",
        "cut_offset_from_pam": 3,
        "notes": (
            "Cria uma incisao na fita que contem o PAM; a edicao e codificada em "
            "uma pegRNA com sequencia de ligacao ao iniciador (PBS) e molde de "
            "transcricao reversa (RTT)."
        ),
    },
}
"""Registro dos sistemas oferecidos na interface. Editores de base e prime editor
reutilizam o reconhecimento por PAM do SpCas9; suas metricas especificas sao
adaptadas nas funcoes dedicadas."""

PFS_FORBIDDEN_CAS13: frozenset[str] = frozenset({"G"})
"""Base do sitio flanqueador do protoespacador (PFS) 3' desfavoravel ao LwaCas13a;
um G nessa posicao reduz a atividade."""

DOENCH_SCORES: Dict[int, Dict[str, float]] = {
    1: {"A": 0.02, "C": 0.01, "G": 0.00, "T": -0.01},
    2: {"A": 0.01, "C": 0.02, "G": 0.01, "T": -0.02},
    3: {"A": 0.00, "C": 0.02, "G": 0.03, "T": -0.02},
    4: {"A": -0.01, "C": 0.01, "G": 0.02, "T": 0.00},
    5: {"A": 0.00, "C": 0.02, "G": 0.02, "T": -0.01},
    6: {"A": 0.01, "C": 0.01, "G": 0.02, "T": -0.01},
    7: {"A": 0.00, "C": 0.02, "G": 0.03, "T": -0.02},
    8: {"A": -0.01, "C": 0.02, "G": 0.02, "T": -0.01},
    9: {"A": 0.00, "C": 0.01, "G": 0.03, "T": -0.02},
    10: {"A": 0.01, "C": 0.02, "G": 0.03, "T": -0.02},
    11: {"A": 0.00, "C": 0.02, "G": 0.03, "T": -0.02},
    12: {"A": 0.00, "C": 0.02, "G": 0.03, "T": -0.02},
    13: {"A": -0.01, "C": 0.02, "G": 0.04, "T": -0.02},
    14: {"A": 0.00, "C": 0.03, "G": 0.04, "T": -0.03},
    15: {"A": 0.00, "C": 0.03, "G": 0.04, "T": -0.03},
    16: {"A": -0.01, "C": 0.03, "G": 0.05, "T": -0.03},
    17: {"A": 0.00, "C": 0.03, "G": 0.05, "T": -0.03},
    18: {"A": -0.01, "C": 0.04, "G": 0.06, "T": -0.04},
    19: {"A": 0.00, "C": 0.04, "G": 0.06, "T": -0.04},
    20: {"A": -0.02, "C": 0.04, "G": 0.08, "T": -0.05},
}
"""Pesos posicionais ilustrativos (posicoes 1-20) inspirados no espirito do
Rule Set 2 (Doench et al., 2016). Valores fixados a mao, sem ajuste a dados
experimentais; nao reproduzem o modelo publicado."""

_IUPAC_SETS: Dict[str, set] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}
"""Conjuntos de bases permitidas por simbolo IUPAC, para casar PAMs ambiguos."""

_COMPLEMENT: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
}
"""Complemento de cada base, incluindo codigos de ambiguidade IUPAC."""

MAX_OFF_TARGET_SCAN_NT: int = 8_000
"""Teto tecnico da varredura local de off-target (janela deslizante nas duas fitas)."""

MAX_GUIDE_SCAN_NT: int = 50_000
"""Teto da varredura de PAMs. Genes longos cabem; cromossomos devem ser recortados."""

MAX_GUIDES_DETAILED: int = 50
"""Numero maximo de guias que recebem varredura de semente no alvo inteiro.

A contagem de copias da semente e O(guias x comprimento). PAM sites alem deste
teto continuam listados com eficiencia heuristica; Seed_copies fica vazio
(nao calculado), nunca zero inventado.
"""

INTERFACE_CAS_SYSTEMS: tuple[str, ...] = ("SpCas9",)
"""Unico sistema oferecido na interface: SpCas9 e o mais usado, com PAM NGG denso."""

EXAMPLE_EMX1_SPCAS9: str = "GAGTCCGAGCAGAAGAAGAAGGG"
"""Sitio EMX1 humano publicado (spacer 20 nt + PAM GGG). Sequencia publica de
Cong et al., Science 339:819-823 (2013) e Mali et al., Science 339:823-826
(2013); nao e uma sequencia inventada."""

SPCAS9_MIN_NT: int = 23
"""Minimo de bases para um sitio SpCas9 completo: 20 nt de spacer + PAM NGG."""

SPCAS9_SEED_NT: int = 12
"""Bases PAM-proximais do spacer usadas para contar repeticoes da semente
(posicoes 9-20 do protoespacador de 20 nt). Orientacao: adjacentes ao PAM,
no 3' do protoespacador de 20 nt (indices 8:20 em 0-based)."""

HSU_2013_HIT_WEIGHTS: tuple[float, ...] = (
    0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583,
)
"""Pesos posicionais do score de um unico hit de Hsu et al., Nature
Biotechnology 31:827-832 (2013), na ordem 5' para 3' do spacer (posicao 1
mais distal do PAM). Vetor reproduzido da implementacao publica do CRISPOR
(Haeussler et al., Genome Biology 17:148, 2016; crispor.py hitScoreM).
Nao e o MIT Specificity Score agregado de todo o genoma."""

HSU_2013_HIT_METHOD: str = (
    "Hsu 2013 MIT single-hit score (crispr.mit.edu / CRISPOR calcHitScore); "
    "scale 0-100; 20 nt spacer; no bulges"
)
"""Rotulo obrigatorio do score Hsu por sitio. Nao usar como 'MIT score'."""

HEURISTIC_EFFICIENCY_METHOD: str = (
    "HelixScope positional efficiency heuristic inspired by Rule Set 2; "
    "not Doench Rule Set 2, not Azimuth"
)
"""Rotulo da eficiencia heuristica. Nunca apresentar como Doench/Azimuth."""

LOCAL_SPECIFICITY_METHOD: str = (
    "Proxy local from the pasted sequence; not CFD and not MIT Specificity Score"
)
"""Rotulo do proxy local de especificidade."""

MAX_REFERENCE_GUIDES: int = 10
"""Teto de guias por busca em referencia fornecida (limite tecnico)."""

MAX_REFERENCE_HITS: int = 200
"""Teto de hits verificados retornados por busca; excesso e RESOURCE_LIMIT."""

REFERENCE_SEARCH_TIMEOUT_S: float = 20.0
"""Timeout da busca em referencia fornecida, em segundos."""

ALGORITHM_SPCAS9_PAM_INDEX: str = "spcas9_pam_index_ngg_nag_no_bulge"
"""Algoritmo real da busca em referencia: indice de PAM NGG/NAG, sem bulges."""


class CrisprError(ValueError):
    """Erro classificado da camada CRISPR.

    Attributes:
        category: INVALID_INPUT, RESOURCE_LIMIT, TIMEOUT, ERROR, UNAVAILABLE
            ou PARSING_ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"

RESTRICTION_IN_GUIDE: tuple[tuple[str, str], ...] = (
    ("EcoRI", "GAATTC"),
    ("BamHI", "GGATCC"),
    ("HindIII", "AAGCTT"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("XhoI", "CTCGAG"),
)
"""Sitios das mesmas enzimas do modulo de DNA, para marcar sobreposicao no
spacer+PAM. Padroes canonicos palindromicos, nao cortam inventados."""


def _clean(seq: str) -> str:
    """Normaliza uma sequencia para maiusculas sem espacos (uso interno)."""
    return "".join(seq.split()).upper()


def _reverse_complement(seq: str) -> str:
    """Retorna o complemento reverso de uma sequencia (uso interno)."""
    normalized = seq.replace("U", "T")
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(normalized))


def _gc_content(seq: str) -> float:
    """Calcula o conteudo GC (%) entre bases canonicas (uso interno).

    Devolve NaN quando nao ha A/C/G/T/U: isso nao e GC 0%.
    """
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T") + seq.count("U")
    canonical = a + c + g + t
    if canonical == 0:
        return float("nan")
    return round(((g + c) / canonical) * 100.0, 2)


def _iupac_match(segment: str, pattern: str) -> bool:
    """Verifica se um segmento casa um padrao IUPAC (uso interno)."""
    if len(segment) != len(pattern):
        return False
    for base, symbol in zip(segment, pattern):
        allowed = _IUPAC_SETS.get(symbol)
        if allowed is None or base not in allowed:
            return False
    return True


def _wallace_tm(primer: str) -> float:
    """Calcula a Tm de um primer pela Regra de Wallace (uso interno)."""
    a = primer.count("A")
    t = primer.count("T")
    g = primer.count("G")
    c = primer.count("C")
    return float(2 * (a + t) + 4 * (g + c))


def find_guides(
    sequence: str,
    cas_system: str = "SpCas9",
    guide_length: int = 20,
) -> List[dict]:
    """Encontra candidatos a RNA guia varrendo as duas fitas de uma sequencia.

    Para cada sistema Cas, localiza o PAM correspondente e extrai o
    protoespacador no comprimento e na orientacao adequados: SpCas9 (PAM NGG,
    20 nt a montante), SaCas9 (PAM NNGRRT, 21 nt a montante), CjCas9
    (PAM NNNNRYAC, 22 nt a montante) e AsCas12a (PAM TTTV, 23 nt a jusante).
    A varredura cobre a fita informada e o seu complemento reverso. Para SpCas9
    o comprimento do guia pode ser ajustado por guide_length; para os demais
    sistemas usa-se o comprimento canonico.

    Cada guia tem o conteudo GC calculado para permitir avaliar os criterios de
    validade (ausencia de TTTT, que atua como terminador da Pol III em promotores
    U6, e GC entre 40% e 80%). Candidatos com bases nao canonicas (por exemplo N)
    no protoespacador sao descartados; os demais sao retornados sem filtragem por
    GC ou TTTT, para que a aplicacao possa distinguir o total de candidatos dos
    guias validos.

    Args:
        sequence: Sequencia de DNA alvo (sera normalizada).
        cas_system: Sistema Cas: "SpCas9", "SaCas9", "AsCas12a" ou "CjCas9".
        guide_length: Comprimento do protoespacador para SpCas9 (padrao 20).

    Returns:
        Lista de dicionarios, um por candidato, com "guide_sequence" (str),
        "pam_sequence" (str), "position" (int, 1-based em coordenadas da fita
        sense), "strand" (str, "+" ou "-"), "gc_content" (float, %),
        "cas_system" (str) e "full_target" (str, protoespacador combinado ao PAM
        na orientacao correta).

    Raises:
        ValueError: Se cas_system for invalido, guide_length nao for positivo,
            ou se a sequencia exceder MAX_GUIDE_SCAN_NT.

    Nota biologica:
        O PAM e indispensavel para o reconhecimento pela nuclease; como o DNA e
        de fita dupla, varrer ambas as fitas revela todos os sitios passiveis de
        clivagem.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    if guide_length <= 0:
        raise ValueError("guide_length deve ser positivo.")

    info = CAS_SYSTEMS[cas_system]
    seq = _clean(sequence)
    if len(seq) > MAX_GUIDE_SCAN_NT:
        raise ValueError(
            f"CRISPR PAM scan is limited to {MAX_GUIDE_SCAN_NT:,} nt. "
            "Paste a gene or locus, not a chromosome-scale sequence."
        )
    if info["target_molecule"] == "RNA":
        extra = set(seq) - {"A", "C", "G", "U"}
        if "T" in seq or extra:
            raise ValueError(
                "Cas13 alveja RNA. Forneca uma sequencia de RNA com A, U, C e G; "
                "DNA (T) nao e transcrito automaticamente."
            )
        return _find_rna_guides(seq, cas_system)

    extra = set(seq) - {"A", "C", "G", "T"}
    if "U" in seq or extra:
        raise ValueError(
            "Este sistema alveja DNA. Forneca uma sequencia de DNA com A, T, C e G; "
            "RNA (U) nao e aceito."
        )
    backbone = info["pam_system"]
    pam = PAM_PATTERNS[backbone]
    lp = len(pam)
    orientation = PAM_ORIENTATION[backbone]
    glen = guide_length if backbone == "SpCas9" else GUIDE_LENGTH_BY_CAS[backbone]

    guides: List[dict] = []
    for strand, strand_seq in (("+", seq), ("-", _reverse_complement(seq))):
        n = len(strand_seq)
        for i in range(0, n - lp + 1):
            segment = strand_seq[i : i + lp]
            if not _iupac_match(segment, pam):
                continue
            if orientation == "downstream":
                g_start = i - glen
                if g_start < 0:
                    continue
                guide = strand_seq[g_start : i]
                full_target = guide + segment
            else:
                g_start = i + lp
                g_end = g_start + glen
                if g_end > n:
                    continue
                guide = strand_seq[g_start:g_end]
                full_target = segment + guide
            if len(guide) < glen or set(guide) - {"A", "C", "G", "T"}:
                continue
            if strand == "+":
                position = g_start + 1
                start_0based = g_start
            else:
                position = n - (g_start + glen) + 1
                start_0based = n - (g_start + glen)
            end_0based = start_0based + glen
            guides.append(
                {
                    "guide_sequence": guide,
                    "pam_sequence": segment,
                    "position": position,
                    "start_0based": start_0based,
                    "end_0based": end_0based,
                    "coordinate_system": "internal_sense",
                    "strand": strand,
                    "gc_content": _gc_content(guide),
                    "cas_system": cas_system,
                    "full_target": full_target,
                    "guide_length": glen,
                    "canonical_spcas9_length": glen == GUIDE_LENGTH_BY_CAS.get(
                        "SpCas9", 20
                    )
                    if backbone == "SpCas9"
                    else None,
                }
            )
    return guides


def _find_rna_guides(sequence: str, cas_system: str = "Cas13") -> List[dict]:
    """Encontra espacadores para um sistema Cas que alveja RNA (uso interno).

    Diferente das nucleases dependentes de PAM, o Cas13 reconhece RNA de fita
    simples sem PAM; a atividade e modulada pelo sitio flanqueador do
    protoespacador (PFS) na extremidade 3' do alvo. Esta funcao varre apenas a
    fita informada (tratada como o RNA alvo), extrai janelas do comprimento do
    espacador e avalia o PFS.

    Args:
        sequence: Sequencia alvo de RNA ou DNA (T e interpretado como U).
        cas_system: Nome do sistema, por padrao "Cas13".

    Returns:
        Lista de dicionarios, um por espacador candidato, com "guide_sequence"
        (str, o espacador complementar ao alvo, em alfabeto de RNA),
        "pam_sequence" (str vazia; nao ha PAM), "position" (int, 1-based no alvo),
        "strand" (str, sempre "+"), "gc_content" (float), "cas_system" (str),
        "full_target" (str, o trecho do alvo) e "pfs" (str, a base do PFS 3').

    Raises:
        Nenhum.

    Nota biologica:
        O espacador da crRNA e complementar ao RNA alvo; a clivagem do Cas13 e
        seguida de atividade colateral sobre RNAs vizinhos, o que inviabiliza a
        nocao de sitio de corte pontual usada nas nucleases de DNA.
    """
    glen = GUIDE_LENGTH_BY_CAS["Cas13"]
    target = _clean(sequence).replace("U", "T")
    complement = {"A": "U", "T": "A", "C": "G", "G": "C", "N": "N"}

    guides: List[dict] = []
    n = len(target)
    for i in range(0, n - glen + 1):
        window = target[i : i + glen]
        if set(window) - {"A", "C", "G", "T"}:
            continue
        spacer = "".join(complement.get(base, "N") for base in reversed(window))
        pfs_base = target[i + glen] if i + glen < n else ""
        pfs_rna = pfs_base.replace("T", "U") if pfs_base else ""
        guides.append(
            {
                "guide_sequence": spacer,
                "pam_sequence": "",
                "position": i + 1,
                "strand": "+",
                "gc_content": _gc_content(window),
                "cas_system": cas_system,
                "full_target": window.replace("T", "U"),
                "pfs": pfs_rna,
            }
        )
    return guides


def score_guide_doench(guide_sequence: str) -> float:
    """Estima a eficiencia de um guia com uma heuristica posicional (0.0 a 1.0).

    Honestidade cientifica: Implementacao heuristica simplificada inspirada em
    Doench et al. 2016 Rule Set 2 (Nature Biotechnology 34:184-191); nao reproduz
    o modelo/algoritmo original publicado. Os pesos em DOENCH_SCORES foram
    definidos a mao, sem treinamento em dados experimentais. Use somente para
    fins educacionais.

    A pontuacao combina:
    - uma componente posicional, somando os pesos de DOENCH_SCORES para o
      nucleotideo presente em cada posicao (1-20), transformada por uma funcao
      logistica para a faixa (0, 1);
    - penalidade quando o conteudo GC e menor que 40% ou maior que 70%;
    - penalidade por homopolimeros (AAAA, TTTT, CCCC, GGGG);
    - bonus quando ha G na posicao 20 (adjacente ao PAM em SpCas9).

    Args:
        guide_sequence: Sequencia do protoespacador (sera normalizada).

    Returns:
        Pontuacao heuristica de eficiencia no intervalo [0.0, 1.0], com 4 casas
        decimais. Retorna 0.0 para entrada vazia.

    Raises:
        Nenhum.

    Nota biologica:
        A escolha do nucleotideo proximo ao PAM e a composicao GC influenciam a
        atividade do sgRNA; os pesos aqui apenas imitam essas tendencias gerais.
    """
    guide = _clean(guide_sequence).replace("U", "T")
    if not guide:
        return 0.0

    positional = 0.0
    for index, base in enumerate(guide[:20], start=1):
        positional += DOENCH_SCORES.get(index, {}).get(base, 0.0)

    score = 1.0 / (1.0 + math.exp(-positional * 4.0))

    gc = _gc_content(guide)
    if not (isinstance(gc, float) and math.isnan(gc)):
        if gc < 40.0 or gc > 70.0:
            score -= 0.20

    for homopolymer in ("AAAA", "TTTT", "CCCC", "GGGG"):
        if homopolymer in guide:
            score -= 0.10

    if len(guide) >= 20 and guide[19] == "G":
        score += 0.10

    return round(max(0.0, min(1.0, score)), 4)


def find_off_targets(
    guide_sequence: str,
    genome_sequence: str,
    max_mismatches: int = 3,
    cas_system: str = "SpCas9",
    require_pam: bool = True,
    include_perfect: bool = False,
    on_target_pam: str = "",
) -> List[dict]:
    """Varre off-targets potenciais por janela deslizante nas duas fitas.

    Desliza uma janela do tamanho do guia ao longo da sequencia fornecida (fita
    informada e complemento reverso), contando mismatches posicao a posicao.
    Sitios com correspondencia perfeita (0 mismatch) sao o alvo pretendido e
    ficam de fora por omissao (include_perfect=False). Para o risco local, os
    mismatches recebem pesos: posicoes 12-20 (regiao seed, proximas ao PAM)
    pesam o dobro (2.0), posicoes 9-11 pesam 1.0 e as posicoes 1-8 (PAM-distal)
    pesam 0.5.

    Honestidade cientifica: trata-se de uma varredura simplificada para fins
    educacionais. Nao modela energias de hibridizacao, bulges nem efeitos
    epigeneticos, e NAO substitui ferramentas validadas como Cas-OFFinder ou
    CRISPOR em uso laboratorial real. O escopo e a sequencia passada, nunca um
    genoma de referencia publico, a menos que essa sequencia seja explicitamente
    a referencia fornecida.

    Args:
        guide_sequence: Sequencia do protoespacador (sera normalizada).
        genome_sequence: Sequencia onde procurar off-targets (sera normalizada).
        max_mismatches: Numero maximo de mismatches aceitos por sitio.
        cas_system: Sistema Cas; SpCas9 exige PAM NGG (NAG e anotado a parte).
        require_pam: Se verdadeiro, so conta janelas seguidas de NGG ou NAG na
            mesma fita. Sem PAM adjacente o sitio nao e um alvo Cas9.
        include_perfect: Se verdadeiro, inclui sitios com 0 mismatch e PAM
            valido. Use na busca em referencia fornecida. Na varredura local do
            alvo de desenho o padrao continua a excluir o on-target.
        on_target_pam: PAM NGG do guia desenhado, necessario para o 23-mer CFD
            on-target. Vazio => CFD por sitio fica UNAVAILABLE, nunca 0.0.

    Returns:
        Lista de dicionarios ordenada por "risk_score" decrescente, cada um com
        "position" (int, 1-based em coordenadas da fita sense), "start_0based",
        "end_0based", "strand" (str), "off_target_sequence" (str), "mismatches"
        (int), "mismatch_positions" (list[int], 1-based), "risk_score" (float),
        "risk_score_method" (str), "pam_sequence", "pam_class", e quando o
        spacer tem 20 nt, "hsu_hit_score" com metodo declarado.

    Raises:
        ValueError: Se max_mismatches for negativo ou o guia for vazio.

    Nota biologica:
        Mismatches longe do PAM tendem a ser melhor tolerados pela nuclease.
        SpCas9 praticamente exige NGG; NAG e um PAM alternativo de atividade
        reduzida (Hsu et al., 2013) e e rotulado, nao equiparado ao NGG.
        Sequencia semelhante sem PAM valido nao e um hit.
    """
    guide = _clean(guide_sequence).replace("U", "T")
    if not guide:
        raise ValueError("O guia nao pode ser vazio.")
    if max_mismatches < 0:
        raise ValueError("max_mismatches nao pode ser negativo.")

    genome = _clean(genome_sequence).replace("U", "T")
    if len(genome) > MAX_OFF_TARGET_SCAN_NT:
        raise ValueError(
            f"A varredura local de off-target nao e executada acima de "
            f"{MAX_OFF_TARGET_SCAN_NT:,} nt. Desligue a analise de off-target "
            "ou cole um alvo menor."
        )
    length = len(guide)
    results: List[dict] = []

    for strand, strand_seq in (("+", genome), ("-", _reverse_complement(genome))):
        n = len(strand_seq)
        for i in range(0, n - length + 1):
            window = strand_seq[i : i + length]
            mismatch_positions: List[int] = []
            weighted = 0.0
            for pos, (a, b) in enumerate(zip(guide, window), start=1):
                if a != b:
                    mismatch_positions.append(pos)
                    if pos <= 8:
                        weighted += 0.5
                    elif 9 <= pos <= 11:
                        weighted += 1.0
                    else:
                        weighted += 2.0
            mismatches = len(mismatch_positions)
            if mismatches > max_mismatches:
                continue
            if mismatches == 0 and not include_perfect:
                continue
            pam_sequence = ""
            pam_class = ""
            if require_pam and cas_system == "SpCas9":
                pam_end = i + length + 3
                if pam_end > n:
                    continue
                pam_sequence = strand_seq[i + length : pam_end]
                if _iupac_match(pam_sequence, "NGG"):
                    pam_class = "NGG"
                elif _iupac_match(pam_sequence, "NAG"):
                    pam_class = "NAG_reduced_activity"
                else:
                    continue
            if strand == "+":
                position = i + 1
                start_0based = i
                end_0based = i + length
                pam_start_0based = i + length if pam_sequence else None
                pam_end_0based = i + length + 3 if pam_sequence else None
            else:
                position = n - (i + length) + 1
                start_0based = n - (i + length)
                end_0based = n - i
                pam_start_0based = start_0based - 3 if pam_sequence else None
                pam_end_0based = start_0based if pam_sequence else None
            hsu = hsu_single_hit_score(guide, window)
            cfd = {
                "score": None,
                "status": "UNAVAILABLE",
                "method": "",
                "model": "",
                "version": "",
                "scale": "",
                "reason": (
                    "On-target PAM was not supplied; CFD 23-mer was not built. "
                    "Not reported as 0.0."
                ),
            }
            if on_target_pam and pam_sequence:
                from . import crispr_cfd

                cfd = crispr_cfd.cfd_from_spacer_and_pam(
                    guide, on_target_pam, window, pam_sequence
                )
            results.append(
                {
                    "position": position,
                    "start_0based": start_0based,
                    "end_0based": end_0based,
                    "coordinate_system": "internal_sense",
                    "strand": strand,
                    "off_target_sequence": window,
                    "mismatches": mismatches,
                    "mismatch_positions": mismatch_positions,
                    "risk_score": round(1.0 / (1.0 + weighted), 4),
                    "risk_score_method": (
                        "local weighted mismatch heuristic (not CFD, not MIT, "
                        "not a safety claim)"
                    ),
                    "pam_sequence": pam_sequence,
                    "pam_class": pam_class,
                    "pam_start_0based": pam_start_0based,
                    "pam_end_0based": pam_end_0based,
                    "hsu_hit_score": hsu["score"],
                    "hsu_hit_score_method": hsu["method"],
                    "hsu_hit_score_status": hsu["status"],
                    "cfd_score": cfd.get("score"),
                    "cfd_status": cfd.get("status"),
                    "cfd_method": cfd.get("method"),
                    "cfd_model": cfd.get("model"),
                    "cfd_version": cfd.get("version"),
                    "cfd_scale": cfd.get("scale"),
                    "bulges_supported": False,
                    "search_scope": "pasted_sequence",
                }
            )

    results.sort(key=lambda item: item["risk_score"], reverse=True)
    return results


def _ruleset2_for_guide(guide: Mapping[str, Any]) -> dict:
    from . import crispr_ontarget

    spacer = str(guide.get("guide_sequence") or "")
    return crispr_ontarget.score_doench_ruleset2(
        spacer,
        nuclease=str(guide.get("cas_system") or "SpCas9"),
    )


def _deephf_for_guide(guide: Mapping[str, Any]) -> dict:
    from . import crispr_ontarget

    spacer = str(guide.get("guide_sequence") or "")
    pam = str(guide.get("pam_sequence") or "")
    return crispr_ontarget.score_deephf(
        spacer + pam,
        nuclease=str(guide.get("cas_system") or "SpCas9"),
    )


def rank_guides(guides: List[dict]) -> List[dict]:
    """Pontua e classifica uma lista de guias por eficiencia heuristica.

    Calcula score_guide_doench para cada guia, adiciona as chaves "doench_score",
    "efficiency_method" e "rank" e ordena a lista por "doench_score" decrescente.
    O nome da funcao score_guide_doench permanece; o valor nao e Rule Set 2.

    Args:
        guides: Lista de guias, cada um com ao menos a chave "guide_sequence".

    Returns:
        Nova lista de dicionarios (copias dos originais) com "doench_score"
        (float) e "rank" (int, 1-based) adicionados, ordenada por pontuacao
        decrescente.

    Raises:
        Nenhum.

    Nota biologica:
        Ordenar por eficiencia esperada ajuda a priorizar guias mais promissores
        para validacao experimental.
    """
    scored: List[dict] = []
    for guide in guides:
        enriched = dict(guide)
        enriched["doench_score"] = score_guide_doench(enriched.get("guide_sequence", ""))
        enriched["efficiency_score"] = enriched["doench_score"]
        enriched["efficiency_method"] = HEURISTIC_EFFICIENCY_METHOD
        enriched["efficiency_model"] = "helixscope_positional_heuristic"
        rs2 = _ruleset2_for_guide(enriched)
        deephf = _deephf_for_guide(enriched)
        enriched["ruleset2_score"] = rs2.get("score")
        enriched["ruleset2_status"] = rs2.get("status")
        enriched["ruleset2_reason"] = rs2.get("reason")
        enriched["deephf_score"] = deephf.get("score")
        enriched["deephf_status"] = deephf.get("status")
        enriched["deephf_reason"] = deephf.get("reason")
        scored.append(enriched)

    scored.sort(key=lambda item: item["doench_score"], reverse=True)
    for rank, guide in enumerate(scored, start=1):
        guide["rank"] = rank
    return scored


def _pick_primer(provider: Callable[[int], Optional[str]]) -> tuple:
    """Escolhe o melhor primer entre comprimentos de 18 a 28 bp (uso interno).

    Args:
        provider: Funcao que recebe um comprimento e retorna a sequencia do
            primer candidato (ou None se nao houver bases suficientes).

    Returns:
        Tupla (sequencia, Tm) do primer escolhido, priorizando Tm dentro de
        55-65 graus Celsius (Wallace) e, em seguida, a Tm mais proxima de 60.

    Raises:
        ValueError: Se nenhum primer puder ser construido.
    """
    best = None
    for length in range(18, 29):
        candidate = provider(length)
        if candidate is None or len(candidate) < length:
            continue
        tm = _wallace_tm(candidate)
        in_range = 55.0 <= tm <= 65.0
        key = (0 if in_range else 1, abs(tm - 60.0))
        if best is None or key < best[0]:
            best = (key, candidate, tm)
    if best is None:
        raise ValueError("Nao ha bases suficientes para desenhar o primer.")
    return best[1], best[2]


def design_primer_pair(
    guide_sequence: str,
    flanking_sequence: str,
    amplicon_size: int = 500,
) -> dict:
    """Desenha um par de primers flanqueando o sitio de corte do guia.

    Localiza o guia na sequencia flanqueadora (em qualquer das duas fitas),
    estima o sitio de corte tipo SpCas9 (corte cego cerca de 3 pb a montante do
    PAM, ou seja, entre as posicoes 17 e 18 do protoespacador) e desenha primers
    direto e reverso de 18 a 28 pb, escolhendo o comprimento cuja Tm (Regra de
    Wallace) fique entre 55 e 65 graus Celsius; quando nenhum candidato atinge a
    faixa, escolhe o de Tm mais proxima de 60 graus.

    Args:
        guide_sequence: Sequencia do protoespacador (sera normalizada).
        flanking_sequence: Sequencia genomica que contem o sitio alvo.
        amplicon_size: Tamanho alvo do amplicon, em pares de base.

    Returns:
        Dicionario com "forward_primer" (str), "reverse_primer" (str),
        "forward_tm" (float, graus Celsius), "reverse_tm" (float, graus Celsius),
        "amplicon_size" (int, tamanho efetivo do amplicon) e "cut_site_position"
        (int, 1-based na fita usada).

    Raises:
        ValueError: Se o guia ou a sequencia flanqueadora forem vazios, se
            amplicon_size nao for positivo ou se o guia nao for encontrado em
            nenhuma das fitas da sequencia flanqueadora.

    Nota biologica:
        Primers flanqueando o sitio de corte permitem amplificar a regiao
        editada para validacao por sequenciamento de Sanger e analises de
        eficiencia como TIDE ou ICE.
    """
    guide = _clean(guide_sequence)
    flank = _clean(flanking_sequence)
    if not guide:
        raise ValueError("O guia nao pode ser vazio.")
    if not flank:
        raise ValueError("A sequencia flanqueadora nao pode ser vazia.")
    if amplicon_size <= 0:
        raise ValueError("amplicon_size deve ser positivo.")

    work = flank
    idx = work.find(guide)
    if idx == -1:
        work = _reverse_complement(flank)
        idx = work.find(guide)
        if idx == -1:
            raise ValueError(
                "guide_sequence nao encontrado na sequencia flanqueadora (em "
                "nenhuma das fitas)."
            )

    cut_site = idx + len(guide) - 3
    half = amplicon_size // 2
    left = max(0, cut_site - half)
    right = min(len(work), cut_site + half)

    forward_primer, forward_tm = _pick_primer(
        lambda length: work[left : left + length]
        if left + length <= len(work)
        else None
    )
    reverse_primer, reverse_tm = _pick_primer(
        lambda length: _reverse_complement(work[right - length : right])
        if right - length >= 0
        else None
    )

    return {
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer,
        "forward_tm": round(forward_tm, 2),
        "reverse_tm": round(reverse_tm, 2),
        "amplicon_size": right - left,
        "cut_site_position": cut_site + 1,
    }


def generate_report(guides: List[dict], sequence_name: str = "Target") -> pd.DataFrame:
    """Monta um relatorio tabular dos guias para exportacao.

    Para cada guia produz uma linha com a classificacao, a sequencia, o PAM, a
    posicao, a fita, o conteudo GC, a pontuacao heuristica e os indicadores de
    validade (ausencia de TTTT e GC entre 40% e 80%). Se um guia nao tiver
    "doench_score" ou "rank", esses valores sao calculados ou atribuidos.

    Args:
        guides: Lista de guias (idealmente ja processada por rank_guides).
        sequence_name: Nome do alvo, usado como nome do indice do DataFrame.

    Returns:
        DataFrame com as colunas "Rank", "Guide_Sequence", "PAM", "Position",
        "Strand", "GC_Content_pct", "Efficiency_heuristic", "TTTT_Free",
        "Valid_GC", "Seed_copies", "Microhomology_bp", "Restriction" e
        "Spacer_RNA", ordenado por "Rank".

    Raises:
        Nenhum.

    Nota biologica:
        Consolidar as metricas em uma tabela facilita a comparacao e a selecao
        de guias para validacao experimental.
    """
    rows: List[dict] = []
    for offset, guide in enumerate(guides, start=1):
        sequence = guide.get("guide_sequence", "")
        gc = guide.get("gc_content", _gc_content(_clean(sequence)))
        doench = guide.get("doench_score")
        if doench is None:
            doench = score_guide_doench(sequence)
        rows.append(
            {
                "Rank": guide.get("rank", offset),
                "Guide_Sequence": sequence,
                "PAM": guide.get("pam_sequence", ""),
                "Position": guide.get("position", 0),
                "Strand": guide.get("strand", ""),
                "GC_Content_pct": gc,
                "Efficiency_heuristic": doench,
                "Efficiency_method": guide.get(
                    "efficiency_method", HEURISTIC_EFFICIENCY_METHOD
                ),
                "TTTT_Free": "TTTT" not in sequence,
                "Valid_GC": 40.0 <= gc <= 80.0,
                "Seed_copies": (guide.get("seed") or {}).get("count", ""),
                "Microhomology_bp": (guide.get("microhomology") or {}).get("length", ""),
                "Restriction": ",".join(guide.get("restriction_sites") or []),
                "Spacer_RNA": (guide.get("order") or {}).get("spacer_rna", ""),
            }
        )

    columns = [
        "Rank", "Guide_Sequence", "PAM", "Position", "Strand",
        "GC_Content_pct", "Efficiency_heuristic", "Efficiency_method",
        "TTTT_Free", "Valid_GC",
        "Seed_copies", "Microhomology_bp", "Restriction", "Spacer_RNA",
    ]
    frame = pd.DataFrame(rows, columns=columns)
    if not frame.empty:
        frame = frame.sort_values(by="Rank").reset_index(drop=True)
    frame.index.name = _safe_report_index_name(sequence_name)
    return frame


def _safe_report_index_name(sequence_name: str) -> str:
    """Evita nomes de indice que o Excel interpretaria como formula (interno)."""
    text = str(sequence_name or "Target").strip() or "Target"
    if text[:1] in {"=", "+", "-", "@"}:
        text = "_" + text[1:]
    return text[:80]


def list_supported_systems() -> List[str]:
    """Lista os nomes dos sistemas Cas suportados pela interface.

    Args:
        Nenhum.

    Returns:
        Lista dos nomes de sistema aceitos por find_guides e system_info, na
        ordem de definicao de CAS_SYSTEMS.

    Raises:
        Nenhum.

    Nota biologica:
        A escolha do sistema define o PAM exigido, o comprimento do espacador e o
        mecanismo de edicao, que por sua vez determinam quais metricas fazem
        sentido calcular.
    """
    return list(CAS_SYSTEMS.keys())


def list_interface_systems() -> List[str]:
    """Lista os sistemas Cas oferecidos na interface Streamlit.

    Args:
        Nenhum.

    Returns:
        Somente ("SpCas9",). Os demais sistemas continuam no modulo para testes
        e uso programatico, mas nao aparecem na aba CRISPR.

    Raises:
        Nenhum.

    Nota biologica:
        SpCas9 (PAM NGG) e o sistema com maior densidade de sitios e maior
        ecossistema experimental; PAMs mais longos (SaCas9 NNGRRT) deixam
        sequencias curtas sem nenhum guia.
    """
    return list(INTERFACE_CAS_SYSTEMS)


def system_info(cas_system: str) -> dict:
    """Retorna os metadados de um sistema Cas registrado.

    Args:
        cas_system: Nome do sistema, conforme list_supported_systems().

    Returns:
        Copia do dicionario de metadados do sistema, incluindo PAM exibido,
        molecula alvo, tipo de editor e notas.

    Raises:
        ValueError: Se o sistema nao estiver registrado.

    Nota biologica:
        Expor os metadados torna explicito, por exemplo, que um editor de base nao
        gera quebra de fita dupla e que o Cas13 age sobre RNA.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    return dict(CAS_SYSTEMS[cas_system])


def model_availability() -> Dict[str, dict]:
    """Declara quais modelos e recursos CRISPR estao ou nao disponiveis.

    CFD por sitio e o agregado MIT/CFD sobre uma busca COMPLETED passaram a
    ser executados com matrizes e formulas publicadas. Rule Set 2, Azimuth,
    DeepHF, MIT genoma-wide, Cas-OFFinder ausente, inDelphi e o complexo 3D
    permanecem UNAVAILABLE. A heuristica local nao e relabeled como Rule Set 2.

    Args:
        Nenhum.

    Returns:
        Dicionario nome -> {available, reason, requires}.

    Raises:
        Nenhum.

    Nota biologica:
        Agregar Hsu/CFD sobre um FASTA de locus nao e o MIT Specificity Score
        de um assembly humano. Cas-OFFinder so e AVAILABLE se o binario real
        estiver instalado.
    """
    from . import crispr_cfd, crispr_casoffinder, crispr_ontarget

    cfd_info = crispr_cfd.cfd_model_record()
    cof = crispr_casoffinder.detect_cas_offinder()
    from . import genome_store, engine_validation

    public_ready = genome_store.any_public_assembly_ready()
    live = engine_validation.live_record("Cas-OFFinder")
    live_ok = bool(
        live
        and live.get("ok")
        and str(live.get("status") or "") == engine_validation.STATUS_LIVE_VALIDATED
    )
    genome_wide_ready = bool(cof.get("available")) and public_ready and live_ok
    rs2 = crispr_ontarget.ruleset2_availability()
    deephf = crispr_ontarget.deephf_availability()
    return {
        "on_target_doench_ruleset2": {
            "available": False,
            "reason": rs2["reason"],
            "requires": rs2["requires"],
            "paper": rs2["paper"],
            "identical_to_azimuth": True,
        },
        "on_target_azimuth": {
            "available": False,
            "reason": (
                "Azimuth 2 is Rule Set 2 (Doench et al. 2016), not a second on-target "
                "model. The Microsoft Research repository is archived (2024) and is "
                "not a supported runtime here."
            ),
            "requires": rs2["requires"],
            "identical_to_ruleset2": True,
        },
        "on_target_deephf": {
            "available": False,
            "reason": deephf["reason"],
            "requires": deephf["requires"],
            "paper": deephf["paper"],
        },
        "specificity_cfd": {
            "available": True,
            "reason": (
                "Per-site CFD (Doench et al. 2016) with embedded published "
                "mismatch and PAM matrices (CRISPOR calcCfdScore / GuideMaker "
                "dump). Scale 0-1. Not Rule Set 2. Not a genome-wide score."
            ),
            "requires": "verified 23-mer spacer+PAM pair (20 nt spacer + 3 nt PAM)",
            "model": cfd_info["model"],
            "version": cfd_info["version"],
            "scale": cfd_info["scale"],
        },
        "specificity_cfd_guide_aggregate": {
            "available": True,
            "reason": (
                "Guide-level CFD specificity (GuideScan / CRISPOR May 2019) "
                "100/(100+sum) on CFD hits of a COMPLETED search. Scope-limited. "
                "NOT_RUN/TIMEOUT/RESOURCE_LIMIT never become 0. Not genome-wide "
                "unless that search processed a complete declared assembly."
            ),
            "requires": "completed off-target search with verified CFD per hit",
            "model": "Doench2016_CFD_guide_specificity",
            "version": "GuideScan-CRISPOR-100-over-100-plus-sum",
            "scale": "0-100",
        },
        "specificity_mit": {
            "available": False,
            "reason": (
                "The genome-wide MIT Specificity Score Sguide aggregates Hsu "
                "hit scores over a complete public assembly. This app does not "
                "download or process hg38/mm39. Summing hits on a plasmid or "
                "locus is a scope-limited aggregate (see "
                "specificity_mit_aggregate_over_search), not MIT genome-wide."
            ),
            "requires": "complete declared assembly search plus Hsu aggregate",
        },
        "specificity_mit_aggregate_over_search": {
            "available": True,
            "reason": (
                "CRISPOR calcMitGuideScore = 100/(100+sum of Hsu 2013 "
                "single-hit scores) over a COMPLETED search, excluding the "
                "intended on-target, including additional exact copies. "
                "Alt-PAM Hsu is multiplied by 0.2 as in CRISPOR. Not "
                "genome-wide unless the search was."
            ),
            "requires": "completed off-target search with Hsu single-hit scores",
            "model": "Hsu2013_MIT_Sguide_aggregate",
            "version": "CRISPOR-calcMitGuideScore",
            "scale": "0-100",
        },
        "hsu_2013_single_hit": {
            "available": True,
            "reason": (
                "Per-site Hsu 2013 score for a verified 20 nt spacer pair "
                "(CRISPOR calcHitScore / crispr.mit.edu). Scale 0-100. Not the "
                "genome-wide MIT Specificity Score. No bulge term."
            ),
            "requires": "20 nt spacer and a verified site with computed mismatches",
        },
        "provided_reference_offtarget_search": {
            "available": True,
            "reason": (
                "PAM-indexed scan of a user FASTA, workspace DNA, or NCBI "
                "nucleotide record already in this session, within size limits. "
                "Not genome-wide. No automatic assembly download."
            ),
            "requires": "validated reference FASTA in this session",
        },
        "cas_offinder": {
            "available": bool(cof.get("available")),
            "reason": str(cof.get("reason") or ""),
            "requires": (
                "allowlisted cas-offinder executable (env HELIXSCOPE_CAS_OFFINDER "
                "or PATH) and a working OpenCL ICD. Not a Python reimplementation."
            ),
            "version": str(cof.get("version") or ""),
        },
        "genome_wide_off_targets": {
            "available": genome_wide_ready,
            "reason": (
                "Cas-OFFinder and a checksum-verified READY public assembly "
                "are present. A search is labelled genome-wide only after "
                "status COMPLETED_FULL_REFERENCE (full assembly, not truncated). "
                "The HelixScope TEST REFERENCE is never genome-wide."
                if genome_wide_ready
                else (
                    "Genome-wide means a checksum-verified READY public assembly "
                    "(GRCh38.p14, T2T-CHM13v2.0, or GRCm39) was searched to "
                    "completion by Cas-OFFinder. Opening CRISPR does not download "
                    "genomes. A provided FASTA or TEST REFERENCE is not genome-wide. "
                    "Partial, timed-out, cancelled, or capped jobs do not produce "
                    "a global specificity score."
                )
            ),
            "requires": (
                "Cas-OFFinder binary, READY assembly with NCBI MD5 verification, "
                "and COMPLETED_FULL_REFERENCE. Auth, quotas and public workers "
                "remain UNAVAILABLE."
            ),
        },
        "indel_frameshift_profile": {
            "available": False,
            "reason": (
                "Indel profiles need trained models such as inDelphi or ForeCasT."
            ),
            "requires": "inDelphi or ForeCasT",
        },
        "cas9_guide_dna_3d_complex": {
            "available": False,
            "reason": (
                "No synthetic Cas9+guide+DNA structure is generated. PDB "
                "pointers (4UN3, 4OO8) are catalogued for a future mapping "
                "workflow. Coordinates are not loaded."
            ),
            "requires": "deposited Cas-guide-DNA structure and validated mapping",
        },
        "distributed_workers": {
            "available": False,
            "reason": (
                "Local subprocess workers exist for Cas-OFFinder genome jobs. "
                "They are not cloud workers, not multi-tenant, and not a "
                "public job queue. Authentication and per-user quotas remain "
                "UNAVAILABLE."
            ),
            "requires": "authentication, quotas and isolated public workers",
        },
    }


def poly_t_signal(guide_sequence: str, run_length: int = 4) -> dict:
    """Detecta sinais de terminacao poli-T no espacador.

    A RNA polimerase III, usada em promotores U6/H1 para expressar o sgRNA,
    termina a transcricao ao encontrar uma corrida de timinas (uridinas no RNA).
    Uma corrida de quatro ou mais T no espacador pode truncar o guia.

    Args:
        guide_sequence: Sequencia do espacador (sera normalizada; U vira T).
        run_length: Tamanho minimo da corrida de T considerada um terminador.

    Returns:
        Dicionario com "has_signal" (bool), "max_run" (int, maior corrida de T) e
        "run_length" (int, o limiar usado).

    Raises:
        ValueError: Se run_length nao for positivo.

    Nota biologica:
        Guias com poli-T interno costumam ser expressos de forma incompleta, o que
        reduz ou elimina a atividade de edicao.
    """
    if run_length <= 0:
        raise ValueError("run_length deve ser positivo.")
    guide = _clean(guide_sequence).replace("U", "T")
    max_run = 0
    current = 0
    for base in guide:
        if base == "T":
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return {
        "has_signal": max_run >= run_length,
        "max_run": max_run,
        "run_length": run_length,
    }


def self_complementarity(guide_sequence: str, min_loop: int = 3) -> dict:
    """Estima o risco de auto-complementaridade (grampo) de um espacador.

    Procura o maior grampo (stem-loop) possivel pareando o espacador consigo
    mesmo: para cada par de posicoes separado por ao menos min_loop bases, estende
    o pareamento Watson-Crick para fora e mede o comprimento do braco. O resultado
    e classificado em risco baixo, moderado ou alto pelo comprimento do braco.

    Args:
        guide_sequence: Sequencia do espacador (sera normalizada; U vira T).
        min_loop: Numero minimo de bases no laco entre os dois bracos do grampo.

    Returns:
        Dicionario com "max_stem" (int, maior braco pareado), "risk" (str: "low",
        "moderate" ou "high") e "min_loop" (int, o valor usado).

    Raises:
        ValueError: Se min_loop for negativo.

    Nota biologica:
        Estruturas secundarias estaveis no espacador competem com o dobramento do
        scaffold do sgRNA e podem impedir a formacao do complexo ribonucleoproteico
        ativo. Esta e uma estimativa puramente baseada em sequencia, nao um calculo
        de energia livre.

        Implementacao heuristica simplificada inspirada em regras de desenho de
        sgRNA sobre auto-complementaridade; nao reproduz o modelo/algoritmo
        original publicado.
    """
    if min_loop < 0:
        raise ValueError("min_loop nao pode ser negativo.")
    guide = _clean(guide_sequence).replace("U", "T")
    pair = {"A": "T", "T": "A", "C": "G", "G": "C"}
    n = len(guide)
    max_stem = 0
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            stem = 0
            left = i
            right = j
            while (
                left >= 0
                and right < n
                and right - left > min_loop
                and pair.get(guide[left]) == guide[right]
            ):
                stem += 1
                left -= 1
                right += 1
            max_stem = max(max_stem, stem)

    if max_stem >= 6:
        risk = "high"
    elif max_stem >= 4:
        risk = "moderate"
    else:
        risk = "low"
    return {"max_stem": max_stem, "risk": risk, "min_loop": min_loop}


def off_target_mismatch_summary(
    off_targets: List[dict], max_mismatches: int = 3
) -> Dict[int, int]:
    """Conta os off-targets locais por numero de mismatches.

    Args:
        off_targets: Lista retornada por find_off_targets.
        max_mismatches: Maior numero de mismatches a incluir no resumo.

    Returns:
        Dicionario que mapeia cada contagem de mismatches, de 0 ate
        max_mismatches, para o numero de sitios com aquele numero de diferencas.
        A chave 0 representa sitios identicos adicionais na sequencia fornecida.

    Raises:
        ValueError: Se max_mismatches for negativo.

    Nota biologica:
        Sitios com poucos mismatches, sobretudo perto do PAM, sao os que a nuclease
        tem maior chance de clivar; distribuir os off-targets por numero de
        mismatches e a leitura direta do risco relativo.
    """
    if max_mismatches < 0:
        raise ValueError("max_mismatches nao pode ser negativo.")
    summary: Dict[int, int] = {i: 0 for i in range(0, max_mismatches + 1)}
    for off in off_targets:
        mm = int(off.get("mismatches", -1))
        if 0 <= mm <= max_mismatches:
            summary[mm] += 1
    return summary


def local_specificity_proxy(off_targets: List[dict]) -> float:
    """Calcula um proxy local de especificidade a partir da varredura embutida.

    O proxy e 1 / (1 + soma dos risk_score dos off-targets encontrados na
    sequencia fornecida). Sem off-targets, vale 1.0. NAO e o CFD nem o MIT score:
    considera apenas a sequencia colada, nunca o genoma inteiro.

    Args:
        off_targets: Lista retornada por find_off_targets para o guia analisado.

    Returns:
        Valor em (0, 1], com 4 casas decimais; quanto maior, mais especifico
        dentro do escopo local analisado.

    Raises:
        Nenhum.

    Nota biologica:
        A especificidade real depende de todos os sitios do genoma e das energias
        de hibridizacao; este proxy so ordena guias entre si dentro da regiao
        fornecida e nao deve ser lido como probabilidade absoluta de off-target.

        Implementacao heuristica simplificada inspirada no conceito de agregacao de
        risco de off-targets do MIT specificity score; nao reproduz o
        modelo/algoritmo original publicado.
    """
    total = sum(float(off.get("risk_score", 0.0)) for off in off_targets)
    return round(1.0 / (1.0 + total), 4)


def classify_guide_risk(
    off_summary: Dict[int, int], specificity_proxy: float
) -> dict:
    """Classifica o risco global de off-target de um guia em Low, Moderate ou High.

    Criterio, aplicado apenas a sequencia fornecida: risco alto se houver qualquer
    sitio adicional com 0 ou 1 mismatch, ou se o proxy local de especificidade for
    menor que 0.5; risco moderado se houver sitios com 2 mismatches ou se o proxy
    ficar entre 0.5 e 0.8; risco baixo caso contrario.

    Args:
        off_summary: Resumo de off_target_mismatch_summary.
        specificity_proxy: Valor de local_specificity_proxy.

    Returns:
        Dicionario com "level" (str: "Low", "Moderate" ou "High") e "criterion"
        (str, a regra aplicada).

    Raises:
        Nenhum.

    Nota biologica:
        Mismatches proximos do PAM sao mal tolerados pela nuclease, entao sitios
        com pouquissimas diferencas concentram o risco real de clivagem
        indesejada.

        Implementacao heuristica simplificada inspirada nas faixas de risco de
        especificidade do CRISPOR; nao reproduz o modelo/algoritmo original
        publicado.
    """
    near_perfect = off_summary.get(0, 0) + off_summary.get(1, 0)
    two_mismatch = off_summary.get(2, 0)
    criterion = (
        "Escopo local apenas. Alto: algum sitio com 0 ou 1 mismatch, ou proxy < "
        "0.5. Moderado: sitios com 2 mismatches, ou proxy entre 0.5 e 0.8. Baixo: "
        "caso contrario."
    )
    if near_perfect > 0 or specificity_proxy < 0.5:
        level = "High"
    elif two_mismatch > 0 or specificity_proxy < 0.8:
        level = "Moderate"
    else:
        level = "Low"
    return {"level": level, "criterion": criterion}


def composite_score(
    efficiency: float, specificity_proxy: Optional[float], weight_efficiency: float = 0.5
) -> Optional[float]:
    """Combina eficiencia heuristica e proxy de especificidade em um unico valor.

    Quando o proxy de especificidade nao esta disponivel (varredura de off-target
    nao executada), retorna None em vez de fingir uma especificidade.

    Args:
        efficiency: Eficiencia heuristica no intervalo [0, 1].
        specificity_proxy: Proxy local de especificidade em [0, 1], ou None.
        weight_efficiency: Peso da eficiencia na media ponderada; o complemento
            pondera a especificidade. Deve estar entre 0 e 1.

    Returns:
        Media ponderada em [0, 1] com 4 casas decimais, ou None quando o proxy de
        especificidade nao foi calculado.

    Raises:
        ValueError: Se weight_efficiency estiver fora de [0, 1].

    Nota biologica:
        Um guia util precisa cortar bem o alvo e pouco fora dele; combinar as duas
        dimensoes evita priorizar guias eficientes porem pouco especificos.

        Implementacao heuristica simplificada inspirada na combinacao de eficiencia
        e especificidade adotada por ferramentas como o CRISPOR; nao reproduz o
        modelo/algoritmo original publicado.
    """
    if not 0.0 <= weight_efficiency <= 1.0:
        raise ValueError("weight_efficiency deve estar entre 0 e 1.")
    if specificity_proxy is None:
        return None
    combined = weight_efficiency * efficiency + (1.0 - weight_efficiency) * specificity_proxy
    return round(combined, 4)


def count_spcas9_pams(sequence: str) -> dict:
    """Conta motivos NGG e NAG nas duas fitas e quantos comportam spacer de 20 nt.

    Args:
        sequence: DNA alvo (sera normalizado).

    Returns:
        Dicionario com "length", "ngg_plus", "ngg_minus", "ngg_total",
        "nag_total", "guides_possible_plus", "guides_possible_minus" e
        "guides_possible".

    Raises:
        Nenhum.

    Nota biologica:
        SpCas9 reconhece NGG. NAG e um PAM alternativo de atividade reduzida e
        nao e usado para desenhar guias aqui; so e contado para transparencia.
    """
    seq = _clean(sequence).replace("U", "T")
    rc = _reverse_complement(seq)

    def _scan(strand_seq: str) -> tuple[int, int, int]:
        ngg = 0
        nag = 0
        usable = 0
        n = len(strand_seq)
        for i in range(0, n - 2):
            motif = strand_seq[i : i + 3]
            if _iupac_match(motif, "NGG"):
                ngg += 1
                if i >= 20:
                    usable += 1
            elif _iupac_match(motif, "NAG"):
                nag += 1
        return ngg, nag, usable

    ngg_plus, nag_plus, usable_plus = _scan(seq)
    ngg_minus, nag_minus, usable_minus = _scan(rc)
    return {
        "length": len(seq),
        "ngg_plus": ngg_plus,
        "ngg_minus": ngg_minus,
        "ngg_total": ngg_plus + ngg_minus,
        "nag_total": nag_plus + nag_minus,
        "guides_possible_plus": usable_plus,
        "guides_possible_minus": usable_minus,
        "guides_possible": usable_plus + usable_minus,
    }


def diagnose_spcas9_target(sequence: str) -> dict:
    """Explica, com contagens reais, por que a varredura SpCas9 pode falhar.

    Args:
        sequence: Texto colado pelo usuario.

    Returns:
        Dicionario com as chaves de count_spcas9_pams, mais "too_short" (bool),
        "has_u" (bool), "ok" (bool) e "messages" (list[str] em ingles, para a UI).

    Raises:
        Nenhum.

    Nota biologica:
        Sem NGG com 20 nt a montante nao existe sitio SpCas9. Uma fita AGCT
        repetida, por exemplo, pode nao ter nenhum GG e portanto nenhum guia.
    """
    raw = sequence or ""
    seq = _clean(raw)
    counts = count_spcas9_pams(seq)
    has_u = "U" in seq
    too_short = counts["length"] < SPCAS9_MIN_NT
    messages: List[str] = []
    if not seq:
        messages.append("Paste a DNA sequence (A, T, C, G).")
    elif has_u and "T" not in seq:
        messages.append(
            "This looks like RNA (contains U). SpCas9 cuts DNA. Paste DNA with T, "
            "not U."
        )
    elif too_short:
        messages.append(
            f"SpCas9 needs at least {SPCAS9_MIN_NT} bp (20 nt spacer + NGG PAM). "
            f"This input has {counts['length']} bp after cleaning."
        )
    if seq and not has_u and counts["ngg_total"] == 0:
        messages.append(
            "No NGG motif on either strand. SpCas9 cannot bind without NGG. "
            "A GG dinucleotide is required. Repeats such as AGCT often have none."
        )
    elif seq and counts["ngg_total"] > 0 and counts["guides_possible"] == 0:
        messages.append(
            f"Found {counts['ngg_total']} NGG motif(s), but none have 20 nt of "
            "DNA upstream on that strand, so no complete spacer can be extracted."
        )
    ok = (
        bool(seq)
        and not has_u
        and not too_short
        and counts["guides_possible"] > 0
    )
    if ok:
        messages.append(
            f"{counts['guides_possible']} complete SpCas9 site(s) "
            f"(NGG with 20 nt spacer) in {counts['length']} bp."
        )
    counts.update(
        {
            "too_short": too_short,
            "has_u": has_u,
            "ok": ok,
            "messages": messages,
        }
    )
    return counts


def seed_occurrences(guide_sequence: str, sequence: str, seed_nt: int = SPCAS9_SEED_NT) -> dict:
    """Conta ocorrencias exatas da semente PAM-proximal nas duas fitas.

    Args:
        guide_sequence: Spacer de 20 nt.
        sequence: DNA alvo.
        seed_nt: Comprimento da semente a partir da extremidade 3' do spacer.

    Returns:
        Dicionario com "seed" (str), "count" (int, incluindo o alvo) e
        "unique" (bool, True quando count == 1).

    Raises:
        ValueError: Se seed_nt nao for positivo.

    Nota biologica:
        Mismatches na semente (proximos ao PAM) bloqueiam Cas9 com mais forca
        do que mismatches distais. Contar copias identicas da semente na
        sequencia colada e uma observacao, nao um score genomico.
    """
    if seed_nt <= 0:
        raise ValueError("seed_nt deve ser positivo.")
    guide = _clean(guide_sequence).replace("U", "T")
    genome = _clean(sequence).replace("U", "T")
    seed = guide[-seed_nt:] if len(guide) >= seed_nt else guide
    if not seed:
        return {"seed": "", "count": 0, "unique": False}
    rc = _reverse_complement(genome)
    count = 0
    for text in (genome, rc):
        start = 0
        while True:
            found = text.find(seed, start)
            if found < 0:
                break
            count += 1
            start = found + 1
    return {"seed": seed, "count": count, "unique": count == 1}


def longest_homopolymer(guide_sequence: str) -> dict:
    """Mede a maior corrida de uma mesma base no spacer.

    Args:
        guide_sequence: Spacer (sera normalizado).

    Returns:
        Dicionario com "base" (str), "length" (int) e "has_run4" (bool).

    Raises:
        Nenhum.
    """
    guide = _clean(guide_sequence).replace("U", "T")
    if not guide:
        return {"base": "", "length": 0, "has_run4": False}
    best_base = guide[0]
    best_len = 1
    current_base = guide[0]
    current_len = 1
    for base in guide[1:]:
        if base == current_base:
            current_len += 1
        else:
            current_base = base
            current_len = 1
        if current_len > best_len:
            best_len = current_len
            best_base = current_base
    return {"base": best_base, "length": best_len, "has_run4": best_len >= 4}


def restriction_overlaps(full_target: str) -> List[str]:
    """Lista enzimas cujos sitios cabem inteiros no spacer+PAM.

    Args:
        full_target: Protoespacador concatenado ao PAM (23 nt para SpCas9).

    Returns:
        Nomes das enzimas com sitio palindromico presente, na ordem de
        RESTRICTION_IN_GUIDE.

    Raises:
        Nenhum.

    Nota biologica:
        Um sitio de restricao no alvo pode ser usado para genotipar o corte, ou
        pode ser destruido pela edicao. A busca e correspondencia exata do
        padrao, nas duas orientacoes.
    """
    seq = _clean(full_target).replace("U", "T")
    hits: List[str] = []
    for name, pattern in RESTRICTION_IN_GUIDE:
        rc = _reverse_complement(pattern)
        if pattern in seq or rc in seq:
            hits.append(name)
    return hits


def microhomology_at_cut(sequence: str, cut_site: Optional[int], max_k: int = 8) -> dict:
    """Observa micro-homologia imediata ao redor do corte cego SpCas9.

    Compara o k-mer terminal a esquerda do corte com o k-mer inicial a direita
    (k de 2 a 8) na fita sense. Nao prediz eficiencia de indel nem frameshift.

    Args:
        sequence: DNA alvo.
        cut_site: Posicao 1-based do corte (cut_site_position), ou None.
        max_k: Maior k-mer examinado (2 a 8).

    Returns:
        Dicionario com "length" (int), "motif" (str) e "method" (str).

    Raises:
        ValueError: Se max_k for menor que 2.

    Nota biologica:
        Micro-homologias flanqueando uma DSB podem favorecer delecoes por MMEJ.
        Isto e uma leitura da sequencia colada, nao o preditor inDelphi/ForeCasT.

        Implementacao heuristica simplificada inspirada em Bae et al., Nature
        Methods 11:705-706 (2014), Microhomology-Predictor; nao reproduz o
        modelo/algoritmo original publicado.
    """
    if max_k < 2:
        raise ValueError("max_k deve ser pelo menos 2.")
    empty = {
        "length": 0,
        "motif": "",
        "method": (
            "Longest identical k-mer (2-8 nt) immediately flanking the blunt cut "
            "on the sense strand. Observation only; not inDelphi or ForeCasT."
        ),
    }
    if cut_site is None:
        return empty
    seq = _clean(sequence).replace("U", "T")
    cut = int(cut_site) - 1
    if cut < 0 or cut > len(seq):
        return empty
    left = seq[max(0, cut - 25) : cut]
    right = seq[cut : cut + 25]
    for k in range(min(max_k, len(left), len(right)), 1, -1):
        motif = left[-k:]
        if motif == right[:k]:
            empty = dict(empty)
            empty["length"] = k
            empty["motif"] = motif
            return empty
    return empty


def order_ready_oligos(guide: dict) -> dict:
    """Monta as formas DNA e RNA do spacer para encomenda, sem bases extras.

    Args:
        guide: Guia com "guide_sequence" e "pam_sequence".

    Returns:
        Dicionario com "spacer_dna", "spacer_rna", "pam" e "full_target_dna".
        Nao adiciona scaffold de tracrRNA; isso depende do fornecedor.

    Raises:
        Nenhum.

    Nota biologica:
        O protoespacador de 20 nt e o que a Cas9 usa para reconhecer o DNA. O
        sgRNA completo inclui o scaffold constante, que nao e gerado aqui.
    """
    spacer = _clean(guide.get("guide_sequence", "")).replace("U", "T")
    pam = _clean(guide.get("pam_sequence", "")).replace("U", "T")
    return {
        "spacer_dna": spacer,
        "spacer_rna": spacer.replace("T", "U"),
        "pam": pam,
        "full_target_dna": spacer + pam,
    }


def cut_site_position(guide: dict, cas_system: str) -> Optional[int]:
    """Estima a posicao do sitio de corte de um guia na sequencia fornecida.

    Para nucleases de corte cego tipo SpCas9, o corte ocorre cerca de 3 pb a
    montante do PAM. Para nucleases de corte escalonado tipo Cas12a/Cas12b, usa-se
    um deslocamento representativo distal ao PAM. Editores de base e prime editor
    produzem incisao na mesma posicao de referencia da nuclease de origem. Cas13
    nao possui sitio de corte pontual e retorna None.

    Args:
        guide: Dicionario de guia com "position", "strand" e o comprimento do
            espacador implicito em "guide_sequence".
        cas_system: Nome do sistema Cas.

    Returns:
        Posicao 1-based do sitio de corte em coordenadas da fita sense, ou None
        quando o conceito nao se aplica (Cas13).

    Raises:
        ValueError: Se o sistema nao estiver registrado.

    Nota biologica:
        Conhecer o sitio de corte orienta o desenho de primers de validacao e a
        previsao de qual regiao do gene sera afetada pela edicao.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    info = CAS_SYSTEMS[cas_system]
    if info["target_molecule"] == "RNA":
        return None
    offset = info.get("cut_offset_from_pam")
    if offset is None:
        return None

    glen = len(_clean(guide.get("guide_sequence", "")))
    position = int(guide.get("position", 1))
    strand = guide.get("strand", "+")
    backbone = info["pam_system"]
    orientation = PAM_ORIENTATION[backbone]

    if orientation == "downstream":
        distance_from_5p = glen - offset
    else:
        distance_from_5p = offset

    if strand == "+":
        return position + distance_from_5p
    return position + glen - distance_from_5p


def base_editing_window(guide: dict, cas_system: str) -> dict:
    """Avalia a janela de edicao de um editor de base para um guia.

    Args:
        guide: Dicionario de guia com "guide_sequence".
        cas_system: Nome do sistema; deve ser um editor de base.

    Returns:
        Dicionario com "editor" (str, ex.: "C>T"), "window" (tuple[int, int],
        posicoes 1-based no espacador), "target_positions" (list[int], posicoes da
        base editavel dentro da janela), "editable" (bool) e "product" (str,
        descricao do resultado esperado).

    Raises:
        ValueError: Se o sistema nao for um editor de base.

    Nota biologica:
        Editores de base so convertem a base alvo quando ela cai na janela de
        atividade da deaminase; a ausencia da base na janela significa que aquele
        guia nao produz a edicao pretendida, ainda que corte bem o alvo.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    info = CAS_SYSTEMS[cas_system]
    if info.get("editor") != "base_editor":
        raise ValueError(f"{cas_system} nao e um editor de base.")

    guide_seq = _clean(guide.get("guide_sequence", ""))
    start, end = info["edit_window"]
    edit_from = info["edit_from"]
    edit_to = info["edit_to"]
    target_positions = [
        pos
        for pos in range(start, end + 1)
        if 1 <= pos <= len(guide_seq) and guide_seq[pos - 1] == edit_from
    ]
    return {
        "editor": f"{edit_from}>{edit_to}",
        "window": (start, end),
        "target_positions": target_positions,
        "editable": bool(target_positions),
        "product": (
            f"Converte {edit_from} para {edit_to} nas posicoes {target_positions} "
            f"do espacador"
            if target_positions
            else f"Nenhum {edit_from} na janela {start}-{end}; sem edicao esperada"
        ),
    }


def prime_editing_notes(guide: dict, cas_system: str) -> dict:
    """Descreve o que e possivel e o que nao e possivel afirmar para prime editing.

    Args:
        guide: Dicionario de guia com "position" e "strand".
        cas_system: Nome do sistema; deve ser um prime editor.

    Returns:
        Dicionario com "nick_position" (int ou None, sitio da incisao na sequencia
        fornecida), "designed" (bool, sempre False para a pegRNA completa),
        "reason" (str) e "recommended_tools" (list[str]).

    Raises:
        ValueError: Se o sistema nao for um prime editor.

    Nota biologica:
        O prime editor introduz uma incisao e usa uma pegRNA com PBS e RTT para
        escrever a edicao por transcricao reversa. Desenhar a pegRNA exige otimizar
        comprimento de PBS, comprimento de RTT e a posicao da edicao, o que depende
        de regras e, na pratica, de ferramentas dedicadas.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    info = CAS_SYSTEMS[cas_system]
    if info.get("editor") != "prime_editor":
        raise ValueError(f"{cas_system} nao e um prime editor.")
    return {
        "nick_position": cut_site_position(guide, cas_system),
        "designed": False,
        "reason": (
            "O desenho completo da pegRNA (comprimento de PBS, comprimento de RTT e "
            "codificacao da edicao) nao e calculado aqui; exige regras dedicadas e "
            "a edicao pretendida do usuario."
        ),
        "recommended_tools": ["PrimeDesign", "pegFinder", "PE-Designer"],
    }


def evaluate_guides(
    guides: List[dict],
    sequence: str,
    cas_system: str,
    max_mismatches: int = 3,
    run_off_target: bool = False,
    specificity_limit: int = 25,
) -> List[dict]:
    """Enriquece guias com metricas por guia e ordena por qualidade combinada.

    Para cada guia calcula eficiencia heuristica, poli-T, auto-complementaridade e
    o sitio de corte. A varredura de copias da semente no alvo (O(n) por guia)
    corre apenas nos MAX_GUIDES_DETAILED melhores por eficiencia, para nao
    esgotar CPU em regioes densas em NGG. Quando run_off_target e verdadeiro,
    varre off-targets na sequencia fornecida para os guias mais eficientes (ate
    specificity_limit), resume por numero de mismatches, calcula o proxy local de
    especificidade e classifica o risco global. A ordenacao usa o score combinado
    quando ha especificidade e a eficiencia heuristica caso contrario.

    Args:
        guides: Lista de guias de find_guides.
        sequence: Sequencia alvo fornecida (escopo da varredura de off-target).
        cas_system: Nome do sistema Cas.
        max_mismatches: Maior numero de mismatches na varredura de off-target.
        run_off_target: Se verdadeiro, calcula especificidade local e risco.
        specificity_limit: Numero maximo de guias, por eficiencia, para os quais a
            varredura de off-target e executada; limita o custo computacional.

    Returns:
        Nova lista de dicionarios enriquecidos, com "doench_score",
        "efficiency_method", "poly_t", "self_complementarity", "cut_site",
        "off_target_summary" (ou None), "specificity_proxy" (ou None),
        "specificity_method", "risk" (ou None), "composite_score" (ou None) e
        "rank". Ordenada por qualidade decrescente.

    Raises:
        ValueError: Se o sistema nao estiver registrado ou specificity_limit for
            negativo.

    Nota biologica:
        Priorizar por eficiencia e especificidade ao mesmo tempo aproxima a escolha
        do que um pesquisador faria na bancada: um guia que corta muito o alvo mas
        tambem fora dele nao e util.

        Implementacao heuristica simplificada inspirada no fluxo de priorizacao de
        guias do CRISPOR; nao reproduz o modelo/algoritmo original publicado.
    """
    if cas_system not in CAS_SYSTEMS:
        suportados = ", ".join(sorted(CAS_SYSTEMS))
        raise ValueError(f"Sistema Cas invalido. Use um de: {suportados}.")
    if specificity_limit < 0:
        raise ValueError("specificity_limit nao pode ser negativo.")

    info = CAS_SYSTEMS[cas_system]
    is_rna = info["target_molecule"] == "RNA"
    enriched: List[dict] = []
    for guide in guides:
        item = dict(guide)
        guide_seq = item.get("guide_sequence", "")
        item["doench_score"] = score_guide_doench(guide_seq)
        item["efficiency_score"] = item["doench_score"]
        item["efficiency_method"] = HEURISTIC_EFFICIENCY_METHOD
        item["efficiency_model"] = "helixscope_positional_heuristic"
        item["guide_id"] = _guide_id(item)
        item["sequence_hash"] = provenance.sequence_digest(str(guide_seq))
        if item.get("start_0based") is None and item.get("position"):
            try:
                start = int(item["position"]) - 1
                item["start_0based"] = start
                item["end_0based"] = start + len(str(guide_seq))
                item["coordinate_system"] = "internal_sense"
            except (TypeError, ValueError):
                pass
        item["poly_t"] = poly_t_signal(guide_seq)
        item["self_complementarity"] = self_complementarity(guide_seq)
        item["cut_site"] = None if is_rna else cut_site_position(item, cas_system)
        item["off_target_summary"] = None
        item["specificity_proxy"] = None
        item["specificity_method"] = (
            "Nao calculado (varredura de off-target desligada); NOT_RUN"
        )
        item["risk"] = None
        item["composite_score"] = None
        item["seed"] = None
        if not is_rna:
            item["homopolymer"] = longest_homopolymer(guide_seq)
            item["restriction_sites"] = restriction_overlaps(
                str(item.get("full_target") or guide_seq)
            )
            item["microhomology"] = microhomology_at_cut(sequence, item["cut_site"])
            item["order"] = order_ready_oligos(item)
        enriched.append(item)

    enriched.sort(key=lambda g: g["doench_score"], reverse=True)
    if not is_rna:
        for item in enriched[:MAX_GUIDES_DETAILED]:
            item["seed"] = seed_occurrences(item.get("guide_sequence", ""), sequence)

    if run_off_target and not is_rna:
        genome = _clean(sequence).replace("U", "T")
        if len(genome) > MAX_OFF_TARGET_SCAN_NT:
            raise ValueError(
                f"A varredura local de off-target nao e executada acima de "
                f"{MAX_OFF_TARGET_SCAN_NT:,} nt. Desligue a analise de off-target "
                "ou cole um alvo menor."
            )
        for item in enriched[: specificity_limit]:
            off_targets = find_off_targets(
                item["guide_sequence"],
                sequence,
                max_mismatches,
                cas_system=cas_system,
                require_pam=True,
                on_target_pam=str(item.get("pam_sequence") or ""),
            )
            summary = off_target_mismatch_summary(off_targets, max_mismatches)
            proxy = local_specificity_proxy(off_targets)
            item["off_target_summary"] = summary
            item["off_targets"] = off_targets
            item["n_off_targets"] = len(off_targets)
            item["specificity_proxy"] = proxy
            item["specificity_method"] = LOCAL_SPECIFICITY_METHOD
            item["risk"] = classify_guide_risk(summary, proxy)
            item["composite_score"] = composite_score(item["doench_score"], proxy)

    has_composite = any(item["composite_score"] is not None for item in enriched)
    if has_composite:
        enriched.sort(
            key=lambda g: (
                g["composite_score"] if g["composite_score"] is not None else -1.0,
                g["doench_score"],
            ),
            reverse=True,
        )
    else:
        enriched.sort(key=lambda g: g["doench_score"], reverse=True)

    for rank, item in enumerate(enriched, start=1):
        item["rank"] = rank
    return enriched


def hsu_single_hit_score(guide_sequence: str, off_target_sequence: str) -> dict:
    """Score de um unico sitio segundo Hsu et al. 2013 (MIT hit, nao Sguide).

    Implementa a formula documentada em http://crispr.mit.edu/about e no
    CRISPOR (calcHitScore): produto (1 - peso) nos mismatches, termo de
    distancia media entre mismatches quando ha dois ou mais, e 1/mmCount^2,
    multiplicado por 100. Pesos: HSU_2013_HIT_WEIGHTS.

    Honestidade cientifica: isto NAO e o MIT Specificity Score do guia, que
    agrega hits de todo o genoma. NAO e CFD. So e definido para spacers de
    20 nt sem bulges. Comprimentos diferentes devolvem status UNAVAILABLE,
    nunca um numero inventado.

    Args:
        guide_sequence: Spacer on-target (DNA, T).
        off_target_sequence: Spacer do sitio candidato, mesmo comprimento.

    Returns:
        Dict com score (float|None), status, method, version e parameters.

    Raises:
        Nenhum.

    Nota biologica:
        Hsu et al., Nat Biotechnol 31:827-832 (2013). A implementacao segue
        CRISPOR/Haeussler 2016 para os casos especiais mmCount < 2 (score2=1)
        e mmCount==0 (score3=1, score=100).
    """
    guide = _clean(guide_sequence).replace("U", "T")
    off = _clean(off_target_sequence).replace("U", "T")
    parameters = {
        "weights": list(HSU_2013_HIT_WEIGHTS),
        "scale": "0-100",
        "spacer_nt": 20,
        "bulges": False,
        "source": (
            "Hsu 2013; CRISPOR hitScoreM / calcHitScore "
            "(Haeussler et al., 2016)"
        ),
    }
    if len(guide) != 20 or len(off) != 20:
        return {
            "score": None,
            "status": "UNAVAILABLE",
            "method": HSU_2013_HIT_METHOD,
            "model": "Hsu2013_single_hit",
            "version": "CRISPOR-hitScoreM-20nt",
            "parameters": parameters,
            "reason": (
                "Hsu 2013 single-hit score is defined here only for 20 nt "
                "SpCas9 spacers without length fudging."
            ),
        }
    if set(guide) - {"A", "C", "G", "T"} or set(off) - {"A", "C", "G", "T"}:
        return {
            "score": None,
            "status": "INVALID_INPUT",
            "method": HSU_2013_HIT_METHOD,
            "model": "Hsu2013_single_hit",
            "version": "CRISPOR-hitScoreM-20nt",
            "parameters": parameters,
            "reason": "Spacer pair must be canonical DNA A/C/G/T.",
        }

    dists: List[int] = []
    mm_count = 0
    last_mm: Optional[int] = None
    score1 = 1.0
    max_dist = 19
    for pos in range(20):
        if guide[pos] != off[pos]:
            mm_count += 1
            if last_mm is not None:
                dists.append(pos - last_mm)
            score1 *= 1.0 - HSU_2013_HIT_WEIGHTS[pos]
            last_mm = pos
    if mm_count < 2:
        score2 = 1.0
    else:
        avg_dist = sum(dists) / len(dists)
        score2 = 1.0 / (((max_dist - avg_dist) / float(max_dist)) * 4 + 1)
    if mm_count == 0:
        score3 = 1.0
    else:
        score3 = 1.0 / (mm_count ** 2)
    score = score1 * score2 * score3 * 100.0
    return {
        "score": round(score, 10),
        "status": "COMPUTED",
        "method": HSU_2013_HIT_METHOD,
        "model": "Hsu2013_single_hit",
        "version": "CRISPOR-hitScoreM-20nt",
        "parameters": parameters,
        "reason": "",
    }


def validate_spcas9_guide(
    guide: Mapping[str, Any],
    *,
    target_sequence: str = "",
) -> dict:
    """Valida alfabeto, comprimento 20, PAM NGG, fita e coordenadas internas.

    Args:
        guide: Dict com guide_sequence, pam_sequence, strand, position.
        target_sequence: DNA alvo opcional para conferir o protoespacador.

    Returns:
        Dict com valid (bool), errors (list[str]) e warnings (list[str]).

    Raises:
        Nenhum.

    Nota biologica:
        Um spacer que nao tem NGG real nao e um guia SpCas9. Comprimentos
        diferentes de 20 nt nao sao pontuados pelos metodos SpCas9 desta fase.
    """
    errors: List[str] = []
    warnings: List[str] = []
    spacer = _clean(str(guide.get("guide_sequence") or "")).replace("U", "T")
    pam = _clean(str(guide.get("pam_sequence") or "")).replace("U", "T")
    strand = str(guide.get("strand") or "")
    if not spacer:
        errors.append("guide_sequence is empty.")
    elif set(spacer) - {"A", "C", "G", "T"}:
        errors.append("guide_sequence must be canonical DNA A/C/G/T.")
    if len(spacer) != GUIDE_LENGTH_BY_CAS["SpCas9"]:
        errors.append(
            f"SpCas9 spacer length must be {GUIDE_LENGTH_BY_CAS['SpCas9']} nt "
            f"for the methods implemented here; got {len(spacer)}."
        )
    if not pam:
        errors.append("PAM is missing; SpCas9 requires NGG.")
    elif not _iupac_match(pam, "NGG"):
        errors.append(f"PAM '{pam}' is not NGG on the reported strand.")
    if strand not in {"+", "-"}:
        errors.append('strand must be "+" or "-".')
    position = guide.get("position")
    try:
        pos_i = int(position)
        if pos_i < 1:
            errors.append("position must be a 1-based positive integer.")
    except (TypeError, ValueError):
        errors.append("position is not an integer.")
        pos_i = None
    target = _clean(target_sequence).replace("U", "T") if target_sequence else ""
    if target and spacer and pos_i is not None:
        from . import scientific_checks

        if not scientific_checks.crispr_guide_in_sequence(target, dict(guide)):
            errors.append(
                "Protospacer does not match the target at the reported "
                "internal coordinates."
            )
    if target and pam and pos_i is not None and strand in {"+", "-"} and len(spacer) == 20:
        if not _pam_exists_on_target(target, spacer, pam, strand, pos_i):
            errors.append(
                "Reported PAM does not exist adjacent to the spacer on that strand."
            )
    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "guide_sequence": spacer,
        "pam_sequence": pam,
        "strand": strand,
    }


def _pam_exists_on_target(
    target: str, spacer: str, pam: str, strand: str, position_1based: int
) -> bool:
    """Confere NGG/PAM adjacente nas coordenadas internas sense."""
    start = position_1based - 1
    end = start + len(spacer)
    if start < 0 or end > len(target):
        return False
    if strand == "+":
        return target[start:end] == spacer and target[end : end + len(pam)] == pam
    window = target[start:end]
    if _reverse_complement(window) != spacer:
        return False
    if start < len(pam):
        return False
    pam_sense = target[start - len(pam) : start]
    return _reverse_complement(pam_sense) == pam


def guide_target_span(guide: Mapping[str, Any], target_sequence: str) -> dict:
    """Span do protoespacador na sequencia alvo para mapeamento 3D futuro.

    Args:
        guide: Guia com position 1-based e strand.
        target_sequence: DNA onde o guia foi desenhado.

    Returns:
        Dict com start/end 0-based, strand, subsequence sense, e
        structure_mapping UNAVAILABLE. Nao fabrica complexo Cas.

    Raises:
        CrisprError: INVALID_INPUT se o intervalo nao couber na sequencia.
    """
    from . import scientific_checks

    seq = _clean(target_sequence).replace("U", "T")
    spacer = _clean(str(guide.get("guide_sequence") or ""))
    try:
        start = int(guide.get("start_0based", int(guide.get("position")) - 1))
    except (TypeError, ValueError) as exc:
        raise CrisprError(
            "Guide is missing valid internal coordinates.",
            "INVALID_INPUT",
        ) from exc
    end = start + len(spacer)
    strand = str(guide.get("strand") or "+")
    try:
        span = scientific_checks.make_span(
            start=start,
            end=end,
            sequence=seq,
            strand=strand,
            source="CRISPR protospacer on pasted DNA",
            label=str(guide.get("guide_id") or "guide"),
            kind="guide",
            status="COMPUTED",
        )
    except ValueError as exc:
        raise CrisprError(str(exc), "INVALID_INPUT") from exc
    span["structure_mapping"] = "UNAVAILABLE"
    span["structure_mapping_reason"] = (
        "No deposited Cas-guide-DNA complex is linked. HelixScope does not "
        "fabricate a 3D Cas9 structure."
    )
    span["pam_sequence"] = str(guide.get("pam_sequence") or "")
    return span


def offtarget_cache_key(
    *,
    guide_sequence: str,
    target_hash: str,
    reference_hash: str,
    assembly_declared: str,
    algorithm: str,
    model: str,
    model_version: str,
    max_mismatches: int,
    include_perfect: bool,
    pam_mode: str,
) -> str:
    """Identidade deterministica de um resultado de off-target.

    Args:
        guide_sequence: Spacer.
        target_hash: Hash da sequencia de desenho.
        reference_hash: identity_hash da referencia.
        assembly_declared: Rotulo de assembly (faz parte da identidade).
        algorithm: Nome do algoritmo de busca.
        model: Nome do modelo de score.
        model_version: Versao do modelo.
        max_mismatches: Parametro da busca.
        include_perfect: Se sitios 0-mismatch entram.
        pam_mode: Por exemplo ngg_or_nag.

    Returns:
        SHA-256 hexadecimal.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [
            _clean(guide_sequence).replace("U", "T"),
            str(target_hash or ""),
            str(reference_hash or ""),
            str(assembly_declared or ""),
            str(algorithm or ""),
            str(model or ""),
            str(model_version or ""),
            str(int(max_mismatches)),
            "1" if include_perfect else "0",
            str(pam_mode or ""),
            provenance.HELIXSCOPE_VERSION,
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def scientific_crispr_report(
    *,
    target_sequence: str,
    guide: Mapping[str, Any],
    reference: Optional[Mapping[str, Any]] = None,
    search: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Relatorio cientifico JSON-safe de um guia e, se houver, da busca.

    Args:
        target_sequence: DNA de desenho.
        guide: Guia avaliado.
        reference: Referencia fornecida, se a busca correu.
        search: Envelope da busca de off-target, se existir.

    Returns:
        Dict com sequencia, guia, PAM, fita, coordenadas, modelos, escopo,
        contagem, limitacoes e proveniencia. Nunca converte NOT_RUN em 0 hits.

    Raises:
        Nenhum.
    """
    search_status = "NOT_RUN"
    hit_count: Any = None
    if isinstance(search, Mapping):
        search_status = str(search.get("status") or "NOT_RUN")
        if search_status in {"COMPLETED"}:
            hit_count = int(search.get("verified_hit_count") or 0)
        else:
            hit_count = None
    limitations = [
        "Computational result. Not experimental validation.",
        "Heuristic efficiency is not Doench Rule Set 2 or Azimuth.",
        "Hsu single-hit is not the genome-wide MIT Specificity Score.",
        "Per-site CFD is Doench 2016 / CRISPOR calcCfdScore (scale 0-1).",
        "MIT/CFD guide aggregates are scope-limited to a COMPLETED search.",
        "An off-target overlapping a gene annotation is not 'gene disrupted'.",
        "BLAST similarity is not CRISPR specificity.",
        "MSA conservation is not guide efficiency.",
        "This deploy has no authentication, per-user quotas or distributed workers.",
        "Bulges are not modelled.",
    ]
    if not (isinstance(search, Mapping) and search.get("genome_wide")):
        limitations.append(
            "Search scope is not genome-wide unless a complete declared "
            "assembly was fully processed (it was not)."
        )
    payload = {
        "target_length": len(_clean(target_sequence)),
        "target_hash": provenance.sequence_digest(_clean(target_sequence))
        if target_sequence
        else "",
        "guide_id": guide.get("guide_id"),
        "guide_sequence": guide.get("guide_sequence"),
        "guide_hash": guide.get("sequence_hash")
        or provenance.sequence_digest(str(guide.get("guide_sequence") or "")),
        "PAM": guide.get("pam_sequence"),
        "strand": guide.get("strand"),
        "position_1based_sense": guide.get("position"),
        "start_0based": guide.get("start_0based"),
        "end_0based": guide.get("end_0based"),
        "coordinate_system": guide.get("coordinate_system") or "internal_sense",
        "efficiency_score": guide.get("efficiency_score", guide.get("doench_score")),
        "efficiency_method": guide.get("efficiency_method", HEURISTIC_EFFICIENCY_METHOD),
        "efficiency_model": guide.get("efficiency_model", "helixscope_positional_heuristic"),
        "specificity_score": guide.get("specificity_proxy"),
        "specificity_method": guide.get("specificity_method"),
        "reference_source": (reference or {}).get("source") if reference else "",
        "organism_declared": (reference or {}).get("organism_declared") if reference else "",
        "assembly_declared": (reference or {}).get("assembly_declared") if reference else "",
        "assembly_version_declared": (reference or {}).get("version_declared")
        if reference
        else "",
        "reference_identity_hash": (reference or {}).get("identity_hash")
        if reference
        else "",
        "search_scope": (search or {}).get("scope") if search else "NOT_RUN",
        "search_algorithm": (search or {}).get("algorithm") if search else "",
        "search_status": search_status,
        "verified_hit_count": hit_count,
        "truncated": bool((search or {}).get("truncated")) if search else False,
        "genome_wide": False,
        "guide_specificity": (search or {}).get("guide_specificity") if search else None,
        "verified_assembly": (reference or {}).get("verified_assembly") if reference else "",
        "assembly_identity_status": (
            (reference or {}).get("assembly_identity_status") if reference else ""
        ),
        "limitations": limitations,
        "provenance": {
            "timestamp": provenance.utc_now(),
            "software_version": provenance.HELIXSCOPE_VERSION,
            "kind": "PREDICTED",
        },
    }
    return provenance.json_safe(payload)


def search_off_targets_in_reference(*args: Any, **kwargs: Any) -> dict:
    """Encaminha para o motor de busca em referencia fornecida.

    Ver modules.crispr_offtarget.search_off_targets_in_reference.
    """
    from . import crispr_offtarget

    return crispr_offtarget.search_off_targets_in_reference(*args, **kwargs)


def export_offtarget_rows(search: Mapping[str, Any]) -> List[dict]:
    """Encaminha a exportacao tabular que recusa 0 hits em buscas nao completas."""
    from . import crispr_offtarget

    return crispr_offtarget.export_offtarget_rows(search)


def _guide_id(guide: Mapping[str, Any]) -> str:
    """Identificador estavel: hash curto + posicao + fita."""
    spacer = str(guide.get("guide_sequence") or "")
    digest = provenance.sequence_digest(spacer)[:12]
    position = guide.get("position", "")
    strand = guide.get("strand", "")
    return f"{digest}:{position}:{strand}"

