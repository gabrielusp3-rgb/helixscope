"""Analise de proteinas: traducao e propriedades fisico-quimicas.

Reune traducao de mRNA, composicao de aminoacidos, ponto isoeletrico, massa
molecular, perfil de hidrofobicidade (Kyte-Doolittle) e indice de instabilidade.

As propriedades fisico-quimicas usam o modulo ProtParam do Biopython; a escala
de hidrofobicidade esta embutida como dicionario constante. Nenhuma funcao aqui
importa Streamlit.
"""

from __future__ import annotations

import math
from typing import Dict, List

STANDARD_AMINO_ACIDS: str = "ACDEFGHIKLMNPQRSTVWY"
"""Os 20 aminoacidos padrao em codigo de uma letra."""

KYTE_DOOLITTLE: Dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
"""Escala de hidropatia de Kyte-Doolittle (1982) por aminoacido."""

AMINO_ACID_CATEGORIES: Dict[str, str] = {
    "nonpolar_aliphatic": "GAPVLIM",
    "aromatic": "FYW",
    "polar_uncharged": "STCNQ",
    "positively_charged": "KRH",
    "negatively_charged": "DE",
}
"""Classificacao mutuamente exclusiva dos 20 aminoacidos padrao segundo o esquema
de Lehninger (cadeia lateral apolar alifatica, aromatica, polar sem carga, com
carga positiva e com carga negativa)."""

HYDROPHOBIC_KYTE_DOOLITTLE: str = "AVLIMFC"
"""Aminoacidos com valor positivo na escala de Kyte-Doolittle, isto e,
hidrofobicos por esse criterio especifico."""

ALIPHATIC_INDEX_COEFFICIENTS: Dict[str, float] = {"A": 1.0, "V": 2.9, "I": 3.9, "L": 3.9}
"""Coeficientes do indice alifatico de Ikai (1980), aplicados ao percentual molar
de alanina, valina, isoleucina e leucina."""

SECONDARY_STRUCTURE_TOOLS: List[str] = [
    "PSIPRED (predicao baseada em perfis de PSI-BLAST)",
    "JPred4 (consenso de multiplos preditores)",
    "AlphaFold DB (estrutura tridimensional, da qual a estrutura secundaria e derivada)",
    "DSSP (atribuicao a partir de estrutura experimental resolvida)",
]
"""Ferramentas adequadas para estrutura secundaria, nenhuma delas integrada ao
HelixScope."""


def _require_standard_protein(protein: str) -> str:
    """Exige os 20 aminoacidos padrao e rejeita qualquer outro simbolo.

    Args:
        protein: Sequencia de proteina bruta.

    Returns:
        Sequencia em maiusculas contendo somente os 20 aminoacidos padrao.

    Raises:
        ValueError: Se a sequencia estiver vazia ou contiver simbolos fora dos
            20 aminoacidos padrao (incluindo X, B, Z, U, O, J e parada).

    Nota biologica:
        ProtParam e a escala de Kyte-Doolittle nao definem parametros para
        residuos desconhecidos ou nao canonicos; remove-los em silencio distorce
        pI, GRAVY e massa.
    """
    cleaned = "".join(protein.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    invalid = sorted({residue for residue in cleaned if residue not in STANDARD_AMINO_ACIDS})
    if invalid:
        raise ValueError(
            "Caracteres fora dos 20 aminoacidos padrao: " + ", ".join(invalid) + "."
        )
    return cleaned


def translate(mrna_seq: str) -> str:
    """Traduz uma sequencia de mRNA em proteina pelo codigo genetico padrao.

    Usa Bio.Seq.Seq.translate(table=1, to_stop=True), interrompendo no primeiro
    codon de parada.

    Args:
        mrna_seq: Sequencia de mRNA (sera normalizada para maiusculas).

    Returns:
        Sequencia de aminoacidos (codigo de uma letra), sem o codon de parada.

    Raises:
        ValueError: Se a sequencia for vazia, contiver timina (T) ou nao
            contiver uracila (U), ou seja, se nao for RNA.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        O ribossomo le o mRNA em trincas a partir do codon de iniciacao; a
        traducao termina ao encontrar um codon de parada.
    """
    cleaned = "".join(mrna_seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia de mRNA esta vazia.")
    if "T" in cleaned:
        raise ValueError("A sequencia contem T; forneca RNA (com U), nao DNA.")
    if "U" not in cleaned:
        raise ValueError("A sequencia nao contem U; nao parece ser RNA.")

    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc

    return str(Seq(cleaned).translate(table=1, to_stop=True))


def amino_acid_composition(protein: str) -> dict:
    """Calcula a composicao dos 20 aminoacidos padrao de uma proteina.

    Args:
        protein: Sequencia de proteina (sera normalizada para maiusculas).

    Returns:
        Dicionario que mapeia cada um dos 20 aminoacidos padrao para um
        sub-dicionario com "count" (int) e "frequency" (float, % com 2 casas
        relativa ao comprimento da sequencia).

    Raises:
        ValueError: Se a sequencia estiver vazia ou contiver simbolos fora dos
            20 aminoacidos padrao.

    Nota biologica:
        A composicao de aminoacidos influencia carga, hidrofobicidade,
        estabilidade e funcao da proteina. Residuos desconhecidos (X) nao sao
        contados como aminoacido padrao.
    """
    sanitized = _require_standard_protein(protein)
    total = len(sanitized)
    result: dict = {}
    for symbol in STANDARD_AMINO_ACIDS:
        count = sanitized.count(symbol)
        result[symbol] = {
            "count": count,
            "frequency": round((count / total) * 100.0, 2),
        }
    return result


def isoelectric_point(protein: str) -> float:
    """Calcula o ponto isoeletrico (pI) de uma proteina via Biopython ProtParam.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        O ponto isoeletrico estimado, com 2 casas decimais.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        No pI a carga liquida da proteina e zero; nesse pH a solubilidade tende a
        ser minima e a mobilidade eletroforetica tambem.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    return round(ProteinAnalysis(sanitized).isoelectric_point(), 2)


def molecular_weight_protein(protein: str) -> float:
    """Calcula a massa molecular de uma proteina em kDa, via Biopython ProtParam.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        Massa molecular em quilodaltons (kDa, valor em Daltons dividido por 1000),
        com 3 casas decimais.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        A massa molecular e essencial para interpretar geis SDS-PAGE,
        espectrometria de massa e calculos de concentracao molar.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    weight_da = ProteinAnalysis(sanitized).molecular_weight()
    return round(weight_da / 1000.0, 3)


def hydrophobicity_profile(protein: str, window: int = 9) -> List[float]:
    """Calcula o perfil de hidrofobicidade pela escala de Kyte-Doolittle.

    Aplica uma janela deslizante centrada, calculando a media dos valores de
    hidropatia dos residuos de cada janela. Os valores resultantes sao limitados
    ao intervalo [-4.5, 4.5] da escala.

    Args:
        protein: Sequencia de proteina (sera normalizada para maiusculas).
        window: Tamanho da janela em residuos; deve ser positivo e nao maior que
            o comprimento da proteina.

    Returns:
        Lista de valores medios de hidropatia (float, 3 casas decimais), um por
        janela, com comprimento len(protein) - window + 1.

    Raises:
        ValueError: Se window nao for positivo ou se len(protein) for menor que
            window.

    Nota biologica:
        Perfis de hidrofobicidade ajudam a prever regioes transmembrana (picos
        positivos prolongados) e segmentos expostos ao solvente (vales).
    """
    if window <= 0:
        raise ValueError("window deve ser positivo.")
    cleaned = _require_standard_protein(protein)
    if len(cleaned) < window:
        raise ValueError("len(protein) nao pode ser menor que window.")

    profile: List[float] = []
    for start in range(0, len(cleaned) - window + 1):
        segment = cleaned[start : start + window]
        score = sum(KYTE_DOOLITTLE[residue] for residue in segment) / window
        score = max(-4.5, min(4.5, score))
        profile.append(round(score, 3))
    return profile


def gravy_index(protein: str) -> float:
    """Calcula o GRAVY (Grand Average of Hydropathy) de uma proteina.

    O GRAVY e a media aritmetica dos valores de hidropatia de Kyte-Doolittle de
    todos os residuos, sem janela deslizante.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        Valor medio de hidropatia com 3 casas decimais, no intervalo teorico
        [-4.5, 4.5].

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.

    Nota biologica:
        GRAVY negativo indica proteina globalmente hidrofilica, compativel com
        proteina solavel ou citoplasmatica; GRAVY positivo indica caracter
        hidrofobico global, comum em proteinas de membrana. Ao contrario do perfil
        por janela, o GRAVY nao localiza dominios transmembrana, apenas resume a
        tendencia media da molecula inteira.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    total = sum(KYTE_DOOLITTLE[residue] for residue in sanitized)
    return round(total / len(sanitized), 3)


def net_charge(protein: str, ph: float = 7.0) -> float:
    """Calcula a carga liquida de uma proteina em um pH definido.

    Usa o modulo ProtParam do Biopython, que aplica a equacao de
    Henderson-Hasselbalch aos grupos ionizaveis das cadeias laterais e das
    extremidades N e C terminais.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).
        ph: Valor de pH no qual calcular a carga; deve estar entre 0 e 14.

    Returns:
        Carga liquida em unidades de carga elementar, com 2 casas decimais.
        Valores positivos indicam proteina basica no pH informado.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao ou se ph
            estiver fora do intervalo [0, 14].
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        A carga liquida em pH fisiologico determina interacoes eletrostaticas,
        ligacao a acidos nucleicos e comportamento em cromatografia de troca
        ionica. Proteinas que ligam DNA tendem a carga liquida positiva, ja que o
        esqueleto de fosfato do DNA e negativo.
    """
    if not 0.0 <= ph <= 14.0:
        raise ValueError("ph deve estar entre 0 e 14.")
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    return round(ProteinAnalysis(sanitized).charge_at_pH(ph), 2)


def aliphatic_index(protein: str) -> float:
    """Calcula o indice alifatico de Ikai (1980) de uma proteina.

    O indice e o volume relativo ocupado por cadeias laterais alifaticas e e dado
    por AI = X(Ala) + 2.9 * X(Val) + 3.9 * (X(Ile) + X(Leu)), onde X e o
    percentual molar de cada residuo.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        Indice alifatico com 2 casas decimais. Valores tipicos ficam entre 60 e
        110 para proteinas globulares.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.

    Nota biologica:
        Indices alifaticos altos correlacionam-se com maior termoestabilidade,
        motivo pelo qual proteinas de organismos termofilicos costumam apresentar
        valores elevados. O indice mede volume de cadeias laterais alifaticas, nao
        hidrofobicidade, e por isso e complementar ao GRAVY.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    total = len(sanitized)
    index = 0.0
    for residue, coefficient in ALIPHATIC_INDEX_COEFFICIENTS.items():
        mole_percent = (sanitized.count(residue) / total) * 100.0
        index += coefficient * mole_percent
    return round(index, 2)


def amino_acid_categories(protein: str) -> dict:
    """Agrupa os residuos de uma proteina por natureza da cadeia lateral.

    Aplica a classificacao mutuamente exclusiva de AMINO_ACID_CATEGORIES e
    acrescenta duas categorias derivadas e sobrepostas: o total de residuos
    carregados e os residuos hidrofobicos pelo criterio de Kyte-Doolittle.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        Dicionario que mapeia o nome de cada categoria para um sub-dicionario com
        "count" (int), "frequency" (float, % com 2 casas do total de residuos
        padrao) e "residues" (str, os aminoacidos que compoem a categoria). Alem
        das cinco categorias exclusivas, inclui "charged_total" e
        "hydrophobic_kyte_doolittle", que se sobrepoem as anteriores.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.

    Nota biologica:
        A proporcao entre residuos carregados, polares e apolares determina
        solubilidade, dobramento e localizacao celular. Excesso de residuos
        apolares sugere dominio transmembrana ou nucleo hidrofobico grande;
        excesso de carregados sugere proteina intrinsecamente desordenada ou de
        superficie.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    total = len(sanitized)

    result: dict = {}
    for name, residues in AMINO_ACID_CATEGORIES.items():
        count = sum(sanitized.count(residue) for residue in residues)
        result[name] = {
            "count": count,
            "frequency": round((count / total) * 100.0, 2),
            "residues": residues,
        }

    charged_residues = (
        AMINO_ACID_CATEGORIES["positively_charged"]
        + AMINO_ACID_CATEGORIES["negatively_charged"]
    )
    charged_count = sum(sanitized.count(residue) for residue in charged_residues)
    result["charged_total"] = {
        "count": charged_count,
        "frequency": round((charged_count / total) * 100.0, 2),
        "residues": charged_residues,
    }

    hydrophobic_count = sum(
        sanitized.count(residue) for residue in HYDROPHOBIC_KYTE_DOOLITTLE
    )
    result["hydrophobic_kyte_doolittle"] = {
        "count": hydrophobic_count,
        "frequency": round((hydrophobic_count / total) * 100.0, 2),
        "residues": HYDROPHOBIC_KYTE_DOOLITTLE,
    }
    return result


def aromaticity(protein: str) -> float:
    """Fracao de residuos aromaticos (Phe, Tyr, Trp) segundo ProtParam.

    Args:
        protein: Sequencia de proteina (somente os 20 aminoacidos padrao).

    Returns:
        Fracao no intervalo [0, 1], com 4 casas, igual a (F + Y + W) / n.
        Este e o aromaticity do Biopython ProtParam, nao um score de funcao.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        Aromaticidade descreve a proporcao de cadeias laterais aromaticas. Nao
        identifica sitios ativos nem ligantes.
    """
    sanitized = _require_standard_protein(protein)
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    return round(float(ProteinAnalysis(sanitized).aromaticity()), 4)


def charge_count_profile(protein: str, window: int = 9) -> List[float]:
    """Perfil de carga formal de cadeias laterais em janela deslizante.

    Para cada janela calcula (n(K)+n(R) - n(D)-n(E)) / window. Histidina nao
    entra na conta porque seu pKa (~6.0) deixa a protonacao parcial em pH 7;
    a carga liquida da cadeia inteira continua vindo de ProtParam em
    net_charge().

    Args:
        protein: Sequencia de proteina.
        window: Tamanho da janela; deve ser positivo e nao maior que a proteina.

    Returns:
        Lista de valores medios de carga formal por janela, com 3 casas.

    Raises:
        ValueError: Se window nao for positivo ou for maior que a proteina.

    Nota biologica:
        Este perfil e uma contagem de residuos ionizaveis, nao a carga de
        Henderson-Hasselbalch. Nao prediz transmembrana nem dominio.
    """
    if window <= 0:
        raise ValueError("window deve ser positivo.")
    cleaned = _require_standard_protein(protein)
    if len(cleaned) < window:
        raise ValueError("len(protein) nao pode ser menor que window.")
    profile: List[float] = []
    for start in range(0, len(cleaned) - window + 1):
        segment = cleaned[start : start + window]
        formal = (
            segment.count("K")
            + segment.count("R")
            - segment.count("D")
            - segment.count("E")
        )
        profile.append(round(formal / window, 3))
    return profile


def per_residue_formal_charge(protein: str) -> List[float]:
    """Carga formal por residuo, o mesmo criterio do perfil por janela.

    K e R = +1; D e E = -1; histidina = 0. Nao e Henderson-Hasselbalch e nao
    substitui net_charge(protein, ph).

    Args:
        protein: Sequencia de proteina.

    Returns:
        Lista com um valor por residuo.

    Raises:
        ValueError: Se a sequencia nao for proteina padrao.

    Nota biologica:
        Em pH 7 a histidina fica parcialmente protonada; por isso este valor
        por residuo nao a conta, alinhado a charge_count_profile. A carga
        liquida da cadeia inteira continua a ser net_charge no pH informado.
    """
    cleaned = _require_standard_protein(protein)
    values: List[float] = []
    for residue in cleaned:
        if residue in {"K", "R"}:
            values.append(1.0)
        elif residue in {"D", "E"}:
            values.append(-1.0)
        else:
            values.append(0.0)
    return values


def secondary_structure_availability() -> dict:
    """Declara que o HelixScope nao preve estrutura secundaria com confianca.

    Args:
        Nenhum.

    Returns:
        Dicionario com "available" (bool, sempre False), "reason" (str,
        justificativa) e "recommended_tools" (list[str], ferramentas adequadas).

    Raises:
        Nenhum.

    Nota biologica:
        Predicao confiavel de estrutura secundaria depende de informacao evolutiva
        (perfis de alinhamento multiplo ou acoplamento entre residuos), que este
        projeto nao computa. Metodos baseados apenas em propensao local de
        aminoacidos, como Chou-Fasman e GOR, alcancam cerca de 60 por cento de
        acerto por residuo, insuficiente para sustentar qualquer conclusao
        biologica. O modulo ProtParam do Biopython expoe uma fracao de estrutura
        secundaria por propensao local, deliberadamente nao usada aqui para nao
        apresentar estimativa fraca como resultado.
    """
    return {
        "available": False,
        "reason": (
            "Predicao de estrutura secundaria exige informacao evolutiva de "
            "alinhamentos multiplos, ausente neste projeto. Metodos por propensao "
            "local acertam cerca de 60 por cento dos residuos, o que nao sustenta "
            "conclusao biologica."
        ),
        "recommended_tools": list(SECONDARY_STRUCTURE_TOOLS),
    }


def instability_index(protein: str) -> float:
    """Calcula o indice de instabilidade de uma proteina via Biopython ProtParam.

    Args:
        protein: Sequencia de proteina (residuos nao canonicos sao removidos).

    Returns:
        Indice de instabilidade com 2 casas decimais. Por convencao, a proteina
        e considerada estavel se o indice for menor que 40 e instavel se for
        maior ou igual a 40.

    Raises:
        ValueError: Se a sequencia nao contiver nenhum aminoacido padrao.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        O indice estima a estabilidade in vitro a partir de dipeptideos
        (Guruprasad et al., 1990, via Biopython ProtParam). Nao e uma medida
        experimental de estabilidade da proteina na celula.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    return round(ProteinAnalysis(sanitized).instability_index(), 2)


def composition_is_consistent(protein: str) -> bool:
    """Confere se a soma das contagens de aminoacidos iguala o comprimento.

    Args:
        protein: Sequencia de proteina dos 20 aminoacidos padrao.

    Returns:
        True quando as somas coincidem.

    Raises:
        ValueError: Se a sequencia for invalida.
    """
    sanitized = _require_standard_protein(protein)
    composition = amino_acid_composition(sanitized)
    total = sum(int(item["count"]) for item in composition.values())
    return total == len(sanitized)


def hydrophobicity_profile_records(protein: str, window: int = 9) -> List[dict]:
    """Perfil Kyte-Doolittle com coordenadas por janela.

    Args:
        protein: Sequencia de proteina.
        window: Tamanho da janela.

    Returns:
        Lista de dicts com "start" (0-based), "end" (exclusivo), "midpoint",
        "score", "window" e "method".

    Raises:
        ValueError: Se window for invalido ou maior que a proteina.
    """
    scores = hydrophobicity_profile(protein, window=window)
    records: List[dict] = []
    for index, score in enumerate(scores):
        records.append(
            {
                "start": index,
                "end": index + window,
                "midpoint": index + window // 2,
                "score": score,
                "window": window,
                "method": "Kyte-Doolittle",
            }
        )
    return records


def physicochemical_report(protein: str, ph: float = 7.0) -> dict:
    """Pacote de propriedades fisico-quimicas com metodo e parametros.

    Args:
        protein: Sequencia dos 20 aminoacidos padrao.
        ph: pH da carga liquida; deve estar em [0, 14].

    Returns:
        Dicionario com pI teorico, massa, GRAVY, carga, indice alifatico,
        instabilidade, metodos e status COMPUTED.

    Raises:
        ValueError / RuntimeError: Propagados das funcoes subjacentes.

    Nota biologica:
        pI e carga sao teoricos (ProtParam, escala de pKa de Bjellqvist). GRAVY
        usa Kyte-Doolittle. Instabilidade e o indice de Guruprasad. Nenhum
        desses valores e uma medicao experimental.
    """
    sanitized = _require_standard_protein(protein)
    mw_kda = molecular_weight_protein(sanitized)
    extinction = extinction_coefficients(sanitized)
    return {
        "length": len(sanitized),
        "molecular_weight_kda": mw_kda,
        "molecular_weight_da": round(mw_kda * 1000.0, 2),
        "molecular_weight_method": (
            "Biopython ProtParam average isotopic residue masses minus water "
            "for peptide bonds"
        ),
        "isoelectric_point": isoelectric_point(sanitized),
        "isoelectric_point_kind": "theoretical",
        "isoelectric_point_method": (
            "Biopython ProtParam / Bjellqvist pKa scale; not an experimental pI"
        ),
        "net_charge": net_charge(sanitized, ph=ph),
        "net_charge_ph": ph,
        "net_charge_method": (
            "Henderson-Hasselbalch on ionizable side chains and termini via ProtParam"
        ),
        "gravy": gravy_index(sanitized),
        "gravy_method": "Kyte-Doolittle",
        "aliphatic_index": aliphatic_index(sanitized),
        "aliphatic_index_method": "Ikai (1980)",
        "instability_index": instability_index(sanitized),
        "instability_index_method": (
            "Guruprasad et al. (1990) dipeptide instability, via ProtParam; "
            "computational index, not experimental stability"
        ),
        "aromaticity": aromaticity(sanitized),
        "composition_consistent": composition_is_consistent(sanitized),
        "extinction_coefficient_280_reduced": extinction["reduced_cysteines"],
        "extinction_coefficient_280_cystine": extinction["cystine_bridges"],
        "extinction_unit": "M^-1 cm^-1",
        "extinction_method": (
            "Biopython ProtParam molar extinction at 280 nm. The reduced value "
            "assumes cysteines are reduced. The cystine value assumes every pair "
            "of cysteines forms a cystine. Not an experimental absorbance."
        ),
        "sequence_complexity_bits": sequence_complexity_bits(sanitized),
        "sequence_complexity_method": (
            "Shannon entropy over the 20 standard amino acids, in bits. "
            "This is not the SEG low-complexity algorithm."
        ),
        "seg_low_complexity": "UNAVAILABLE",
        "status": "COMPUTED",
    }


def extinction_coefficients(protein: str) -> dict:
    """Coeficientes de extincao molar a 280 nm via ProtParam.

    Args:
        protein: Sequencia dos 20 aminoacidos padrao.

    Returns:
        Dict com reduced_cysteines, cystine_bridges, wavelength_nm, unit,
        method e status COMPUTED. Os dois coeficientes sao inteiros de ProtParam.

    Raises:
        ValueError: Se a sequencia nao tiver aminoacido padrao.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        Trp, Tyr e Cys dominam a absorbancia a 280 nm. O valor com cistina
        assume pontes dissulfeto; nao e uma medida no espectrofotometro.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    reduced, cystine = ProteinAnalysis(sanitized).molar_extinction_coefficient()
    return {
        "reduced_cysteines": int(reduced),
        "cystine_bridges": int(cystine),
        "wavelength_nm": 280,
        "unit": "M^-1 cm^-1",
        "method": (
            "Biopython ProtParam molar extinction at 280 nm. Reduced cysteines "
            "ignore cystine. The second value assumes every cysteine pair forms "
            "a cystine."
        ),
        "status": "COMPUTED",
    }


def charge_curve(protein: str, *, ph_step: float = 0.5) -> List[dict]:
    """Carga liquida teorica de pH 0 a 14 no passo pedido.

    Args:
        protein: Sequencia dos 20 aminoacidos padrao.
        ph_step: Passo de pH. Deve ser positivo e produzir no maximo 57 pontos.

    Returns:
        Lista de dicts com ph, net_charge, method e status COMPUTED.

    Raises:
        ValueError: Se o passo for invalido ou a sequencia for vazia.

    Nota biologica:
        Cada ponto usa a mesma equacao de Henderson-Hasselbalch de net_charge.
        A curva nao e uma titulacao experimental.
    """
    if not isinstance(ph_step, (int, float)) or isinstance(ph_step, bool):
        raise ValueError("ph_step must be a positive number.")
    step = float(ph_step)
    if step <= 0 or step > 14:
        raise ValueError("ph_step must be in (0, 14].")
    points = int(round(14.0 / step)) + 1
    if points > 57:
        raise ValueError(
            "Charge curve would exceed 57 pH points. Increase ph_step. "
            "This is not a truncated protein."
        )
    rows: List[dict] = []
    ph = 0.0
    while ph <= 14.0 + (step / 10.0):
        clamped = 14.0 if ph > 14.0 else round(ph, 4)
        rows.append(
            {
                "ph": clamped,
                "net_charge": net_charge(protein, ph=clamped),
                "method": "ProtParam charge_at_pH / Henderson-Hasselbalch",
                "status": "COMPUTED",
            }
        )
        if clamped >= 14.0:
            break
        ph += step
    return rows


def sequence_complexity_bits(protein: str) -> float:
    """Entropia de Shannon da composicao de aminoacidos, em bits.

    Args:
        protein: Sequencia dos 20 aminoacidos padrao.

    Returns:
        Entropia com 4 casas. 0.0 e real para um homopolimero. float('nan')
        nao e usado: sequencia vazia levanta ValueError.

    Raises:
        ValueError: Se nao houver aminoacido padrao.

    Nota biologica:
        Entropia baixa indica composicao enviesada. Nao e o algoritmo SEG de
        Wootton e Federhen e nao marca regioes de baixa complexidade por si so.
    """
    sanitized = _require_standard_protein(protein)
    if not sanitized:
        raise ValueError("A sequencia nao contem aminoacidos padrao.")
    total = len(sanitized)
    entropy = 0.0
    for residue in STANDARD_AMINO_ACIDS:
        count = sanitized.count(residue)
        if count:
            probability = count / total
            entropy -= probability * math.log2(probability)
    return round(entropy, 4)


def feature_prediction_availability() -> dict:
    """Declara preditores de dominio, sinal e transmembrana como unavailable.

    Args:
        Nenhum.

    Returns:
        Mapa de capacidade -> {available, reason, recommended_tools}.

    Raises:
        Nenhum.

    Nota biologica:
        Dominios Pfam, peptideo sinal e helices transmembrana exigem modelos
        treinados ou bancos externos. Inferi-los so da composicao seria
        fabricacao.
    """
    return {
        "domains": {
            "available": False,
            "reason": (
                "Protein domain assignment is not implemented. HelixScope does "
                "not query Pfam, InterPro or HMMER."
            ),
            "recommended_tools": ["InterProScan", "Pfam", "HMMER"],
        },
        "signal_peptide": {
            "available": False,
            "reason": (
                "Signal peptide prediction is not implemented. A hydropathy peak "
                "is not a SignalP result."
            ),
            "recommended_tools": ["SignalP 6.0"],
        },
        "transmembrane": {
            "available": False,
            "reason": (
                "Transmembrane helix prediction is not implemented. Kyte-Doolittle "
                "windows are hydropathy, not TMHMM."
            ),
            "recommended_tools": ["DeepTMHMM", "TMHMM"],
        },
        "conservation": {
            "available": False,
            "reason": "No multiple-sequence alignment or phylogenetic model is run.",
            "recommended_tools": ["ConSurf", "Jalview"],
        },
    }
