"""Nucleo de analise de sequencias de DNA do HelixScope.

Reune validacao e deteccao de tipo de sequencia, composicao de nucleotideos,
conteudo GC, temperatura de melting, massa molecular, complemento reverso com
suporte IUPAC, GC skew em janelas, busca de sitios de restricao, deteccao de
ORFs nos seis reading frames e identificacao de ilhas CpG.

Os dados externos (codigos de ambiguidade IUPAC, pesos moleculares e tabela de
codons NCBI Standard #1) estao embutidos como constantes neste modulo. Nenhuma
funcao aqui importa Streamlit; o Biopython e importado de forma tardia apenas na
funcao que usa RestrictionBatch.
"""

from __future__ import annotations

import io
import math
from typing import Dict, List, Optional

from . import scale_profile

IUPAC_DNA_AMBIGUITY: Dict[str, str] = {
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
"""Codigos de ambiguidade IUPAC do DNA mapeados para as bases que representam."""

CANONICAL_DNA_ALPHABET: frozenset[str] = frozenset("ACGT")
"""Alfabeto canonico de DNA aceito na entrada dos modulos (uma fita, A/T/C/G)."""

CANONICAL_RNA_ALPHABET: frozenset[str] = frozenset("ACGU")
"""Alfabeto canonico de RNA aceito na entrada dos modulos (fita simples, A/U/C/G)."""

CANONICAL_PROTEIN_ALPHABET: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
"""Os 20 aminoacidos padrao em codigo de uma letra, unicos aceitos como proteina."""

DNA_ALPHABET: frozenset[str] = CANONICAL_DNA_ALPHABET | frozenset(IUPAC_DNA_AMBIGUITY)
"""Alfabeto de DNA com IUPAC, usado em complemento reverso e busca de motivos."""

RNA_ALPHABET: frozenset[str] = CANONICAL_RNA_ALPHABET | frozenset(IUPAC_DNA_AMBIGUITY)
"""Alfabeto de RNA com IUPAC, usado apenas em funcoes internas de degeneracao."""

PROTEIN_ALPHABET: frozenset[str] = CANONICAL_PROTEIN_ALPHABET
"""Alfabeto de proteina na entrada: somente os 20 aminoacidos padrao."""

PROTEIN_ONLY_SYMBOLS: frozenset[str] = frozenset("DEFHIKLMPQVWY")
"""Residuos que existem no alfabeto proteico e nao pertencem a ACGT/ACGU.
A, C, G, T e U sozinhos nao bastam para classificar uma fita como proteina,
porque coincidem com nucleotideos."""

DNA_COMPLEMENT_IUPAC: Dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
"""Mapa de complementaridade de bases com suporte completo a codigos IUPAC."""

DNA_WEIGHTS: Dict[str, float] = {"A": 313.21, "T": 304.19, "G": 329.21, "C": 289.18}
"""Massas medias (Da, forma sodio) dos nucleotideos monofosfato de DNA."""

RNA_WEIGHTS: Dict[str, float] = {"A": 329.21, "U": 306.17, "G": 345.21, "C": 305.18}
"""Massas medias (Da, forma sodio) dos nucleotideos monofosfato de RNA."""

PHOSPHODIESTER_WATER: float = 18.02
"""Massa de agua (Da) liberada por ligacao fosfodiester formada."""

CODON_TABLE_STANDARD: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
"""Codigo genetico NCBI Standard #1 (DNA): codons para aminoacidos; o asterisco
indica codon de parada."""

START_CODON: str = "ATG"
"""Codon de iniciacao canonico no codigo genetico padrao."""

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
"""Codons de parada no codigo genetico padrao."""

RESTRICTION_ENZYME_NAMES: List[str] = [
    "EcoRI", "BamHI", "HindIII", "NcoI", "NdeI", "XhoI",
]
"""Enzimas de restricao pesquisadas por find_restriction_sites."""

DINUCLEOTIDES: List[str] = [
    first + second for first in "ACGT" for second in "ACGT"
]
"""Os 16 dinucleotideos canonicos, em ordem alfabetica por primeira base."""

MIN_CDS_LENGTH_NT: int = 150
"""Comprimento minimo (nt) adotado para aceitar uma ORF prevista como candidata a
regiao codificadora. Trincas curtas ocorrem por acaso em qualquer sequencia e,
se contadas, inflam artificialmente estatisticas de uso de codons."""

MAX_ORF_OVERLAP_FRACTION: float = 0.5
"""Fracao maxima do proprio comprimento que uma ORF pode sobrepor a uma ORF maior
ja selecionada antes de ser considerada redundante."""

MAX_RAW_INPUT_CHARS: int = 250_000
"""Limite tecnico do texto bruto (incluindo espacos e cabecalhos FASTA) aceito
na validacao. Impede alocar dezenas de megabytes antes da limpeza."""

MAX_INPUT_RESIDUES: int = 200_000
"""Limite tecnico de residuos apos a limpeza, alinhado ao teto de download NCBI.
Nao e um limite biologico: cromossomos e scaffolds maiores devem ser recortados
antes da analise na interface."""

MAX_KMER_K: int = 5
"""Maior k aceito em kmer_counts. 4**5 = 1024 chaves; k maior explode memoria."""

IUPAC_DNA_SYMBOLS: tuple[str, ...] = tuple(
    sorted(CANONICAL_DNA_ALPHABET | frozenset(IUPAC_DNA_AMBIGUITY))
)
"""Simbolos de DNA efetivamente suportados (canonicos + ambiguidades listadas)."""


def _clean(seq: str) -> str:
    """Normaliza uma sequencia removendo espacos e convertendo para maiusculas.

    Args:
        seq: Sequencia bruta de entrada.

    Returns:
        A sequencia sem espacos em branco e em letras maiusculas.

    Nota biologica:
        A limpeza padroniza a entrada para que a contagem de bases independa de
        capitalizacao ou formatacao de texto.
    """
    return "".join(seq.split()).upper()


def _size_rejection_reason(seq: str) -> str:
    """Devolve a mensagem de rejeicao se a entrada exceder os limites tecnicos.

    Args:
        seq: Texto bruto ainda nao normalizado.

    Returns:
        Motivo em ingles quando o texto ou a sequencia limpa e grande demais;
        string vazia quando o tamanho e aceitavel.

    Raises:
        Nenhum.
    """
    if len(seq) > MAX_RAW_INPUT_CHARS:
        return (
            f"Input exceeds the technical limit of {MAX_RAW_INPUT_CHARS:,} "
            f"characters. HelixScope analyzes at most {MAX_INPUT_RESIDUES:,} "
            "residues. Paste a gene, transcript or contig region, not a "
            "chromosome-scale sequence."
        )
    cleaned = _clean(seq)
    if len(cleaned) > MAX_INPUT_RESIDUES:
        return (
            f"Sequence length {len(cleaned):,} exceeds the technical limit of "
            f"{MAX_INPUT_RESIDUES:,} residues. Paste a smaller region of interest."
        )
    return ""


def _empty_validation(cleaned: str = "") -> dict:
    """Monta o dicionario de validacao para sequencia vazia (uso interno)."""
    return {
        "sequence": cleaned,
        "type": "UNKNOWN",
        "length": len(cleaned),
        "is_valid": False,
        "invalid_chars": [],
        "expected": "",
        "rejection_reason": "The sequence is empty.",
    }


def validate_sequence(seq: str) -> dict:
    """Valida uma sequencia e detecta automaticamente seu tipo biologico.

    Remove espacos e quebras de linha, converte para maiusculas e classifica a
    sequencia como DNA, RNA, PROTEIN ou UNKNOWN. DNA valido usa somente A, T, C
    e G; RNA valido usa somente A, U, C e G; proteina valida usa somente os 20
    aminoacidos padrao. Codigos IUPAC, N, B, Z, X e o codon de parada nao sao
    aceitos como entrada valida.

    Args:
        seq: Sequencia bruta de nucleotideos ou aminoacidos.

    Returns:
        Dicionario com "sequence", "type" ("DNA", "RNA", "PROTEIN" ou
        "UNKNOWN"), "length", "is_valid", "invalid_chars", "expected" (str
        vazia neste detector automatico) e "rejection_reason" (str).

    Raises:
        TypeError: Se seq nao for uma string.

    Nota biologica:
        A, C, G e T existem tanto no DNA quanto no codigo de aminoacidos
        (Ala, Cys, Gly, Thr). Uma fita feita so desses quatro simbolos e
        classificada como DNA, nao como proteina. Uracila distingue RNA.
        Residuos como E, F, L ou K forcam a classificacao como proteina.
    """
    if not isinstance(seq, str):
        raise TypeError("A sequencia deve ser uma string.")
    oversized = _size_rejection_reason(seq)
    if oversized:
        return {
            "sequence": "",
            "type": "UNKNOWN",
            "length": 0,
            "is_valid": False,
            "invalid_chars": [],
            "expected": "",
            "rejection_reason": oversized,
        }
    cleaned = _clean(seq)
    chars = set(cleaned)

    if not cleaned:
        return _empty_validation()

    if "T" in chars and "U" in chars:
        mixed = sorted(chars)
        return {
            "sequence": cleaned,
            "type": "UNKNOWN",
            "length": len(cleaned),
            "is_valid": False,
            "invalid_chars": mixed,
            "expected": "",
            "rejection_reason": (
                "The sequence mixes T and U, so it is neither DNA nor RNA."
            ),
        }

    if chars <= CANONICAL_DNA_ALPHABET:
        seq_type = "DNA"
        invalid: List[str] = []
    elif chars <= CANONICAL_RNA_ALPHABET:
        seq_type = "RNA"
        invalid = []
    elif chars <= CANONICAL_PROTEIN_ALPHABET:
        seq_type = "PROTEIN"
        invalid = []
    elif "U" in chars:
        seq_type = "RNA"
        invalid = sorted(chars - CANONICAL_RNA_ALPHABET)
    elif chars & PROTEIN_ONLY_SYMBOLS:
        seq_type = "PROTEIN"
        invalid = sorted(chars - CANONICAL_PROTEIN_ALPHABET)
    else:
        seq_type = "DNA"
        invalid = sorted(chars - CANONICAL_DNA_ALPHABET)
        if invalid and not (chars <= (CANONICAL_DNA_ALPHABET | CANONICAL_RNA_ALPHABET | CANONICAL_PROTEIN_ALPHABET | set("NXBZ*"))):
            seq_type = "UNKNOWN"

    reason = ""
    if invalid or seq_type == "UNKNOWN":
        extra = ", ".join(invalid) if invalid else "none"
        reason = (
            f"Detected type {seq_type} is not a valid canonical sequence. "
            f"Invalid characters: {extra}."
        )

    return {
        "sequence": cleaned,
        "type": seq_type,
        "length": len(cleaned),
        "is_valid": seq_type != "UNKNOWN" and not invalid,
        "invalid_chars": invalid,
        "expected": "",
        "rejection_reason": reason,
    }


def validate_for_molecule(seq: str, molecule: str) -> dict:
    """Valida a sequencia contra o alfabeto de um unico tipo molecular.

    Cada modulo deve chamar esta funcao com o tipo que ele analisa. DNA so
    aceita A, T, C e G; RNA so aceita A, U, C e G; proteina so aceita os 20
    aminoacidos padrao. Uma fita de DNA colada no modulo de RNA (ou o inverso)
    e rejeitada com motivo explicito. Nao ha conversao T/U nem remocao silenciosa
    de simbolos.

    Args:
        seq: Sequencia bruta.
        molecule: "DNA", "RNA" ou "PROTEIN".

    Returns:
        O mesmo dicionario de validate_sequence(), com "expected" igual a
        molecule e "is_valid" verdadeiro somente quando o alfabeto casa com
        esse tipo.

    Raises:
        TypeError: Se seq nao for uma string.
        ValueError: Se molecule nao for DNA, RNA ou PROTEIN.

    Nota biologica:
        O DNA de entrada e uma fita no alfabeto A/T/C/G; a fita complementar e
        calculada quando a analise precisa dela. O RNA e fita simples com U no
        lugar de T. Peptideos compostos so por A, C, G e T sao indistinguiveis
        de DNA e sao rejeitados no modulo de proteina.
    """
    if not isinstance(seq, str):
        raise TypeError("A sequencia deve ser uma string.")
    expected = molecule.strip().upper()
    if expected not in {"DNA", "RNA", "PROTEIN"}:
        raise ValueError("molecule deve ser 'DNA', 'RNA' ou 'PROTEIN'.")

    oversized = _size_rejection_reason(seq)
    if oversized:
        return {
            "sequence": "",
            "length": 0,
            "expected": expected,
            "is_valid": False,
            "invalid_chars": [],
            "type": "UNKNOWN",
            "rejection_reason": oversized,
        }

    cleaned = _clean(seq)
    chars = set(cleaned)
    base = {
        "sequence": cleaned,
        "length": len(cleaned),
        "expected": expected,
        "is_valid": False,
        "invalid_chars": [],
        "type": "UNKNOWN",
        "rejection_reason": "",
    }
    if not cleaned:
        base["rejection_reason"] = f"Provide a non-empty {expected} sequence."
        return base

    alphabets = {
        "DNA": CANONICAL_DNA_ALPHABET,
        "RNA": CANONICAL_RNA_ALPHABET,
        "PROTEIN": CANONICAL_PROTEIN_ALPHABET,
    }
    allowed = alphabets[expected]
    extra = sorted(chars - allowed)

    looks_dna = chars <= CANONICAL_DNA_ALPHABET
    looks_rna = chars <= CANONICAL_RNA_ALPHABET and "U" in chars
    looks_protein = (
        chars <= CANONICAL_PROTEIN_ALPHABET
        and not looks_dna
        and not (chars <= CANONICAL_RNA_ALPHABET)
    )

    if expected == "DNA":
        if looks_rna:
            base["type"] = "RNA"
            base["invalid_chars"] = ["U"]
            base["rejection_reason"] = (
                "This sequence is RNA (contains U). The DNA module accepts only "
                "A, T, C and G."
            )
            return base
        if looks_protein:
            base["type"] = "PROTEIN"
            base["invalid_chars"] = extra
            base["rejection_reason"] = (
                "This sequence is protein. The DNA module accepts only A, T, C "
                "and G."
            )
            return base
        if extra:
            base["type"] = "DNA"
            base["invalid_chars"] = extra
            base["rejection_reason"] = (
                "Invalid characters for DNA (A, T, C, G only): "
                + ", ".join(extra)
                + "."
            )
            return base
        base["type"] = "DNA"
        base["is_valid"] = True
        return base

    if expected == "RNA":
        if chars <= CANONICAL_DNA_ALPHABET and "T" in chars:
            base["type"] = "DNA"
            base["invalid_chars"] = ["T"]
            base["rejection_reason"] = (
                "This sequence is DNA (contains T). The RNA module accepts only "
                "A, U, C and G and does not transcribe DNA."
            )
            return base
        if "T" in chars and "U" in chars:
            base["type"] = "UNKNOWN"
            base["invalid_chars"] = ["T", "U"]
            base["rejection_reason"] = (
                "The sequence mixes T and U. The RNA module accepts only A, U, C "
                "and G."
            )
            return base
        if looks_protein:
            base["type"] = "PROTEIN"
            base["invalid_chars"] = extra
            base["rejection_reason"] = (
                "This sequence is protein. The RNA module accepts only A, U, C "
                "and G."
            )
            return base
        if extra:
            base["type"] = "RNA"
            base["invalid_chars"] = extra
            base["rejection_reason"] = (
                "Invalid characters for RNA (A, U, C, G only): "
                + ", ".join(extra)
                + "."
            )
            return base
        base["type"] = "RNA"
        base["is_valid"] = True
        return base

    if looks_dna:
        base["type"] = "DNA"
        base["rejection_reason"] = (
            "This sequence is DNA (A, T, C, G only). The protein module accepts "
            "the 20 standard amino acids and rejects nucleic-acid input."
        )
        return base
    if looks_rna:
        base["type"] = "RNA"
        base["rejection_reason"] = (
            "This sequence is RNA (contains U). The protein module accepts the "
            "20 standard amino acids."
        )
        return base
    if extra:
        base["type"] = "PROTEIN"
        base["invalid_chars"] = extra
        base["rejection_reason"] = (
            "Invalid characters for protein (20 standard amino acids only): "
            + ", ".join(extra)
            + "."
        )
        return base
    base["type"] = "PROTEIN"
    base["is_valid"] = True
    return base


def nucleotide_composition(seq: str) -> dict:
    """Calcula a composicao de bases de uma sequencia, incluindo ambiguidades.

    Conta A, T, C, G, N e todos os codigos de ambiguidade IUPAC, retornando a
    contagem absoluta e a frequencia relativa (em porcentagem) de cada simbolo.

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        Dicionario que mapeia cada simbolo (A, T, C, G, N e ambiguidades IUPAC)
        para um sub-dicionario com "count" (int) e "frequency" (float, % com 2
        casas decimais, relativa ao comprimento total da sequencia).

    Raises:
        Nenhum.

    Nota biologica:
        Incluir os codigos de ambiguidade permite avaliar vies composicional,
        controle de qualidade e artefatos de sequenciamento sem descartar
        posicoes incertas.
    """
    cleaned = _clean(seq)
    total = len(cleaned)
    symbols = ["A", "T", "U", "C", "G", "N"] + sorted(IUPAC_DNA_AMBIGUITY.keys() - {"N"})
    result: dict = {}
    for symbol in symbols:
        count = cleaned.count(symbol)
        frequency = round((count / total) * 100.0, 2) if total else 0.0
        result[symbol] = {"count": count, "frequency": frequency}
    return result


def gc_content(seq: str) -> float:
    """Calcula o conteudo GC de uma sequencia, em porcentagem.

    Considera apenas as bases canonicas A, T, C e G no denominador, excluindo N
    e codigos de ambiguidade.

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        Percentual de bases G e C entre as bases canonicas, com 2 casas
        decimais. Retorna float('nan') quando nao ha bases canonicas (desconhecido,
        nao zero).

    Raises:
        Nenhum.

    Nota biologica:
        Excluir ambiguidades do denominador torna o valor adequado a calculos
        termodinamicos como a temperatura de melting e o desenho de sondas.
    """
    cleaned = _clean(seq).replace("U", "T")
    a = cleaned.count("A")
    t = cleaned.count("T")
    g = cleaned.count("G")
    c = cleaned.count("C")
    canonical = a + t + g + c
    if canonical == 0:
        return float("nan")
    return round(((g + c) / canonical) * 100.0, 2)


def at_content(seq: str) -> float:
    """Calcula o conteudo AT de uma sequencia, em porcentagem.

    Considera apenas as bases canonicas A, T, C e G no denominador, excluindo N e
    codigos de ambiguidade, de forma simetrica a gc_content().

    Args:
        seq: Sequencia de DNA ou RNA (sera normalizada; U e tratado como T).

    Returns:
        Percentual de bases A e T entre as bases canonicas, com 2 casas decimais.
        Retorna float('nan') quando nao ha bases canonicas (desconhecido, nao zero).

    Raises:
        Nenhum.

    Nota biologica:
        Regioes ricas em AT desnaturam com mais facilidade porque pares A-T fazem
        duas ligacoes de hidrogenio contra tres de G-C; origens de replicacao,
        promotores e regioes de curvatura do DNA costumam ser AT-ricas.
    """
    cleaned = _clean(seq).replace("U", "T")
    a = cleaned.count("A")
    t = cleaned.count("T")
    g = cleaned.count("G")
    c = cleaned.count("C")
    canonical = a + t + g + c
    if canonical == 0:
        return float("nan")
    return round(((a + t) / canonical) * 100.0, 2)


def dinucleotide_frequencies(seq: str) -> dict:
    """Calcula frequencias dos 16 dinucleotideos e sua razao observado/esperado.

    Percorre a sequencia em janelas sobrepostas de duas bases, conta cada
    dinucleotideo canonico e compara a frequencia observada com a esperada sob
    independencia entre posicoes adjacentes. A razao observado/esperado (O/E) de
    um dinucleotideo XY e obtida por
    (freq_XY) / (freq_X * freq_Y), onde freq_XY usa o total de dinucleotideos
    contados e freq_X e freq_Y usam o total de bases canonicas.

    Args:
        seq: Sequencia de DNA ou RNA (sera normalizada; U e tratado como T).

    Returns:
        Dicionario que mapeia cada um dos 16 dinucleotideos canonicos para um
        sub-dicionario com "count" (int), "frequency" (float, % com 2 casas do
        total de dinucleotideos) e "observed_expected" (float com 3 casas, ou
        NaN quando o esperado e zero: indefinido, nao zero).

    Raises:
        Nenhum.

    Nota biologica:
        Sob independencia, O/E fica proximo de 1.0. O caso classico de desvio e o
        CpG em genomas de vertebrados, onde O/E cai para cerca de 0.2 porque
        citosinas metiladas em contexto CpG desaminam para timina ao longo da
        evolucao, depletando o dinucleotideo. Valores de CpG O/E proximos de 1.0
        em regiao genomica sugerem ilha CpG hipometilada, tipicamente promotora.
    """
    cleaned = _clean(seq).replace("U", "T")
    canonical = "".join(base for base in cleaned if base in "ACGT")
    total_bases = len(canonical)

    pair_counts: Dict[str, int] = {pair: 0 for pair in DINUCLEOTIDES}
    total_pairs = 0
    for index in range(len(cleaned) - 1):
        pair = cleaned[index : index + 2]
        if pair in pair_counts:
            pair_counts[pair] += 1
            total_pairs += 1

    base_counts = {base: canonical.count(base) for base in "ACGT"}

    result: dict = {}
    for pair in DINUCLEOTIDES:
        count = pair_counts[pair]
        frequency = round((count / total_pairs) * 100.0, 2) if total_pairs else 0.0
        first_freq = base_counts[pair[0]] / total_bases if total_bases else 0.0
        second_freq = base_counts[pair[1]] / total_bases if total_bases else 0.0
        expected = first_freq * second_freq
        observed = count / total_pairs if total_pairs else 0.0
        odds = round(observed / expected, 3) if expected > 0 else float("nan")
        result[pair] = {
            "count": count,
            "frequency": frequency,
            "observed_expected": odds,
        }
    return result


def gc_sliding_window(seq: str, window: int = 100, step: int = 50) -> List[dict]:
    """Calcula o conteudo GC ao longo da sequencia em janela deslizante.

    Desliza uma janela de tamanho fixo com passo definido e calcula o percentual
    de G e C de cada janela, permitindo localizar variacao composicional interna.

    Args:
        seq: Sequencia de DNA ou RNA (sera normalizada; U e tratado como T).
        window: Tamanho da janela em bases; deve ser positivo e nao maior que o
            comprimento da sequencia.
        step: Passo de deslocamento em bases; deve ser positivo.

    Returns:
        Lista de dicionarios, um por janela, com "start" (int, base 0), "end"
        (int, exclusivo), "midpoint" (int) e "gc_percent" (float com 2 casas,
        ou NaN se a janela nao tiver A/C/G/T).

    Raises:
        ValueError: Se window ou step nao forem positivos, ou se a sequencia for
            menor que window.

    Nota biologica:
        O perfil de GC ao longo do genoma revela isocoros, ilhas CpG, regioes
        transferidas horizontalmente (que costumam divergir do GC medio do
        hospedeiro) e limites entre regiao codificante e nao codificante, ja que
        exons tendem a ser mais ricos em GC que introns.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window e step devem ser positivos.")
    cleaned = _clean(seq).replace("U", "T")
    if len(cleaned) < window:
        raise ValueError("len(seq) nao pode ser menor que window.")
    scale_profile.enforce_window_budget(len(cleaned), window, step)

    profile: List[dict] = []
    for start in range(0, len(cleaned) - window + 1, step):
        end = start + window
        chunk = cleaned[start:end]
        g = chunk.count("G")
        c = chunk.count("C")
        a = chunk.count("A")
        t = chunk.count("T")
        canonical = a + t + g + c
        gc_percent = (
            round(((g + c) / canonical) * 100.0, 2) if canonical else float("nan")
        )
        profile.append(
            {
                "start": start,
                "end": end,
                "midpoint": start + window // 2,
                "gc_percent": gc_percent,
            }
        )
    return profile


def melting_temperature(seq: str, na_concentration: float = 0.05) -> Optional[float]:
    """Estima a temperatura de melting (Tm) de uma sequencia, em graus Celsius.

    Para sequencias com menos de 14 bases canonicas usa a Regra de Wallace
    (Tm = 2*(A+T) + 4*(G+C)). Para 14 bases ou mais usa a formula salt-adjusted
    (Tm = 81.5 + 16.6*log10([Na+]) + 0.41*GC% - 675/N). Acima de 200 bases
    canônicas a funcao devolve None: essas formulas de oligo nao se aplicam.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        na_concentration: Concentracao de sodio em Molar; deve ser positiva.

    Returns:
        Temperatura de melting estimada em graus Celsius, com 2 casas decimais,
        ou None quando N > 200.

    Raises:
        ValueError: Se na_concentration nao for positiva ou se a sequencia nao
            contiver bases canonicas.

    Nota biologica:
        Acima de 14 bp os efeitos ionicos do sodio estabilizam o duplex de forma
        relevante, por isso a formula longa incorpora log10([Na+]) e correcao
        por comprimento, ao contrario da Regra de Wallace para oligos curtos.
        Esta funcao nao implementa o metodo nearest-neighbor de SantaLucia;
        esse calculo vive em modules.thermodynamics e nao substitui silenciosamente
        Wallace nem a formula salt-adjusted. Calculadoras de primers modernas
        (por exemplo IDT) usam nearest-neighbor e nao devem ser comparadas
        numericamente com o valor exibido aqui.
    """
    if na_concentration <= 0:
        raise ValueError("na_concentration deve ser positiva (em Molar).")
    cleaned = _clean(seq)
    a = cleaned.count("A")
    t = cleaned.count("T")
    g = cleaned.count("G")
    c = cleaned.count("C")
    n = a + t + g + c
    if n == 0:
        raise ValueError("A sequencia nao contem bases canonicas para o calculo.")
    if n > 200:
        return None
    if n < 14:
        return round(2.0 * (a + t) + 4.0 * (g + c), 2)
    gc_percent = ((g + c) / n) * 100.0
    tm = 81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc_percent - 675.0 / n
    return round(tm, 2)


def melting_temperature_report(
    seq: str, na_concentration: float = 0.05
) -> dict:
    """Rotula o metodo de Tm realmente aplicavel, sem calcular SantaLucia.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        na_concentration: Concentracao de sodio em Molar; deve ser positiva.

    Returns:
        Dicionario com "value_c" (float ou None), "method" (str: "Wallace",
        "salt-adjusted" ou ""), "status" ("COMPUTED" ou "UNAVAILABLE"),
        "n_canonical" (int), "na_concentration_m" (float) e "reason" (str).

    Raises:
        ValueError: Se na_concentration nao for positiva.

    Nota biologica:
        Wallace e a formula salt-adjusted sao estimativas de oligo. Nenhuma
        delas e o modelo nearest-neighbor de SantaLucia. Acima de 200 bases
        canônicas o valor fica UNAVAILABLE, nao zero.
    """
    if na_concentration <= 0:
        raise ValueError("na_concentration deve ser positiva (em Molar).")
    cleaned = _clean(seq)
    a = cleaned.count("A")
    t = cleaned.count("T")
    g = cleaned.count("G")
    c = cleaned.count("C")
    n = a + t + g + c
    base = {
        "value_c": None,
        "method": "",
        "status": "UNAVAILABLE",
        "n_canonical": n,
        "na_concentration_m": na_concentration,
        "reason": "",
    }
    if n == 0:
        base["reason"] = "No canonical A/C/G/T bases; Tm is not applicable."
        return base
    if n > 200:
        base["reason"] = (
            "Unavailable for sequences with more than 200 canonical bases. "
            "Wallace and salt-adjusted oligo formulas do not apply. "
            "SantaLucia nearest-neighbor Tm is a separate method with its own length domain."
        )
        return base
    try:
        value = melting_temperature(seq, na_concentration=na_concentration)
    except ValueError as exc:
        base["reason"] = str(exc)
        return base
    method = "Wallace" if n < 14 else "salt-adjusted"
    return {
        "value_c": value,
        "method": method,
        "status": "COMPUTED",
        "n_canonical": n,
        "na_concentration_m": na_concentration,
        "reason": "",
    }


def molecular_weight(seq: str, mol_type: str = "DNA") -> float:
    """Calcula a massa molecular media de um acido nucleico, em Daltons.

    Soma as massas dos nucleotideos canonicos e subtrai a agua liberada em cada
    ligacao fosfodiester (uma molecula de agua por ligacao, ou seja, N-1
    moleculas para N nucleotideos).

    Args:
        seq: Sequencia de DNA ou RNA (sera normalizada).
        mol_type: Tipo de molecula, "DNA" ou "RNA".

    Returns:
        Massa molecular media em Daltons, com 2 casas decimais. Retorna
        float('nan') quando nao ha nucleotideos canonicos (desconhecido, nao
        zero).

    Raises:
        ValueError: Se mol_type nao for "DNA" nem "RNA".

    Nota biologica:
        A formacao de cada ligacao fosfodiester entre nucleotideos libera uma
        molecula de agua; por isso a massa do polimero e menor que a soma das
        massas dos monomeros isolados.
    """
    normalized_type = mol_type.strip().upper()
    if normalized_type not in {"DNA", "RNA"}:
        raise ValueError("mol_type deve ser 'DNA' ou 'RNA'.")
    cleaned = _clean(seq)
    weights = DNA_WEIGHTS if normalized_type == "DNA" else RNA_WEIGHTS

    total = 0.0
    count = 0
    for base in cleaned:
        if base in weights:
            total += weights[base]
            count += 1
    if count == 0:
        return float("nan")
    total -= (count - 1) * PHOSPHODIESTER_WATER
    return round(total, 2)


def reverse_complement(seq: str) -> str:
    """Retorna o complemento reverso de uma sequencia, com suporte IUPAC.

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        A sequencia complementar lida no sentido oposto (5'->3' da fita oposta),
        com complementaridade definida tambem para os codigos de ambiguidade
        IUPAC (R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N).

    Raises:
        ValueError: Se a sequencia contiver simbolos sem complemento definido.

    Nota biologica:
        O complemento reverso representa a fita antisentido lida no sentido
        biologico 5'->3'; o suporte a IUPAC preserva a degeneracao das posicoes
        ambiguas.
    """
    cleaned = _clean(seq)
    complemented: List[str] = []
    for base in reversed(cleaned):
        if base not in DNA_COMPLEMENT_IUPAC:
            raise ValueError(f"Simbolo sem complemento definido: {base}.")
        complemented.append(DNA_COMPLEMENT_IUPAC[base])
    return "".join(complemented)


def gc_skew(seq: str, window: int = 100) -> List[float]:
    """Calcula o GC skew em janelas consecutivas ao longo da sequencia.

    Para cada janela de tamanho fixo calcula (G - C) / (G + C). Janelas sem G
    nem C recebem float('nan'). A ultima janela pode ser parcial.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        window: Tamanho da janela em bases; deve ser positivo.

    Returns:
        Lista de valores de GC skew (float), um por janela, com 4 casas
        decimais; janelas sem G e C retornam float('nan').

    Raises:
        ValueError: Se window nao for positivo ou se len(seq) for menor que
            window.

    Nota biologica:
        O GC skew ajuda a localizar origens e terminos de replicacao em genomas,
        pois as fitas lider e tardia acumulam vieses distintos de G e C.
    """
    if window <= 0:
        raise ValueError("window deve ser positivo.")
    cleaned = _clean(seq)
    if len(cleaned) < window:
        raise ValueError("len(seq) nao pode ser menor que window.")
    n_windows = (len(cleaned) + window - 1) // window
    if n_windows > scale_profile.MAX_PLOT_POINTS:
        raise ValueError(
            f"GC skew would produce {n_windows:,} non-overlapping windows. "
            f"Technical plot limit is {scale_profile.MAX_PLOT_POINTS:,}. "
            "Increase the window size. The (G-C)/(G+C) formula is unchanged."
        )

    skews: List[float] = []
    for start in range(0, len(cleaned), window):
        chunk = cleaned[start : start + window]
        g = chunk.count("G")
        c = chunk.count("C")
        if g + c == 0:
            skews.append(float("nan"))
        else:
            skews.append(round((g - c) / (g + c), 4))
    return skews


def at_skew_windows(seq: str, window: int = 100) -> List[float]:
    """Calcula o AT skew em janelas consecutivas ao longo da sequencia.

    Para cada janela de tamanho fixo calcula (A - T) / (A + T). Janelas sem A
    nem T recebem float('nan'). A ultima janela pode ser parcial.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        window: Tamanho da janela em bases; deve ser positivo.

    Returns:
        Lista de valores de AT skew (float), um por janela, com 4 casas
        decimais; janelas sem A e T retornam float('nan').

    Raises:
        ValueError: Se window nao for positivo ou se len(seq) for menor que
            window.

    Nota biologica:
        O AT skew descreve vies composicional complementar ao GC skew; nao
        localiza origens de replicacao por si so.
    """
    if window <= 0:
        raise ValueError("window deve ser positivo.")
    cleaned = _clean(seq)
    if len(cleaned) < window:
        raise ValueError("len(seq) nao pode ser menor que window.")
    n_windows = (len(cleaned) + window - 1) // window
    if n_windows > scale_profile.MAX_PLOT_POINTS:
        raise ValueError(
            f"AT skew would produce {n_windows:,} non-overlapping windows. "
            f"Technical plot limit is {scale_profile.MAX_PLOT_POINTS:,}. "
            "Increase the window size. The (A-T)/(A+T) formula is unchanged."
        )

    skews: List[float] = []
    for start in range(0, len(cleaned), window):
        chunk = cleaned[start : start + window]
        a = chunk.count("A")
        t = chunk.count("T")
        if a + t == 0:
            skews.append(float("nan"))
        else:
            skews.append(round((a - t) / (a + t), 4))
    return skews


def orf_coordinates_on_input(orf: dict, sequence_length: int) -> tuple[int, int]:
    """Converte coordenadas de uma ORF para a fita de entrada (base 0, fim exclusivo).

    Args:
        orf: Dicionario retornado por find_orfs(), com "frame", "start" e "end".
        sequence_length: Comprimento da sequencia de entrada.

    Returns:
        Tupla (start, end) em coordenadas da fita informada, base 0, intervalo
        semiaberto.

    Raises:
        ValueError: Se sequence_length for nao positivo ou se as coordenadas
            resultantes forem inconsistentes.

    Nota biologica:
        Frames negativos sao encontrados no complemento reverso; plotar start/end
        crus desses frames sobre a fita direta deslocaria o ORF para a posicao
        errada.
    """
    if sequence_length <= 0:
        raise ValueError("sequence_length deve ser positivo.")
    start = int(orf["start"])
    end = int(orf["end"])
    frame = int(orf["frame"])
    if frame > 0:
        mapped_start, mapped_end = start, end
    else:
        mapped_start = sequence_length - end
        mapped_end = sequence_length - start
    if mapped_start < 0 or mapped_end > sequence_length or mapped_start >= mapped_end:
        raise ValueError("ORF coordinates are inconsistent with the sequence length.")
    return mapped_start, mapped_end


def find_restriction_sites(seq: str) -> dict:
    """Localiza sitios de enzimas de restricao usando o Biopython.

    Obtem os padroes de reconhecimento das enzimas EcoRI, BamHI, HindIII, NcoI,
    NdeI e XhoI a partir de um RestrictionBatch do Biopython e procura cada
    padrao e tambem o seu complemento reverso na sequencia informada, retornando
    posicoes em um unico sistema de coordenadas (base 1).

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        Dicionario que mapeia o nome de cada enzima para um sub-dicionario com
        "sites" (list[int] de posicoes iniciais do sitio na fita informada,
        base 1 e ordenadas), "count" (int), "pattern" (str da sequencia de
        reconhecimento) e "hits" (lista de dicts com "enzyme", "site",
        "position" base 1 e "strand"). As seis enzimas listadas sao palindromos;
        a fita e reportada como "+" nas coordenadas da entrada.

    Raises:
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        Enzimas de restricao reconhecem sequencias especificas em DNA de fita
        dupla, em geral palindromicas; procurar o padrao e o seu complemento
        reverso cobre sitios presentes em qualquer das duas fitas sem duplicar
        sitios palindromicos.
    """
    try:
        from Bio.Restriction import RestrictionBatch
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc

    cleaned = _clean(seq)
    batch = RestrictionBatch(RESTRICTION_ENZYME_NAMES)

    result: dict = {}
    for enzyme in batch:
        pattern = str(enzyme.site).upper()
        targets = {pattern, reverse_complement(pattern)}
        positions: set[int] = set()
        for target in targets:
            start = cleaned.find(target)
            while start != -1:
                positions.add(start + 1)
                start = cleaned.find(target, start + 1)
        ordered = sorted(positions)
        hits = [
            {
                "enzyme": str(enzyme),
                "site": pattern,
                "position": position,
                "strand": "+",
            }
            for position in ordered
        ]
        result[str(enzyme)] = {
            "sites": ordered,
            "count": len(ordered),
            "pattern": pattern,
            "hits": hits,
        }
    return result


def _translate(seq: str) -> str:
    """Traduz uma sequencia de DNA em proteina ate o primeiro codon de parada.

    Args:
        seq: Sequencia de DNA ja normalizada, em quadro de leitura 0.

    Returns:
        Sequencia de aminoacidos (codigo de uma letra), sem incluir o codon de
        parada; codons com bases ambiguas sao traduzidos como "X".

    Nota biologica:
        A traducao para no primeiro codon de parada encontrado, reproduzindo o
        termino da sintese proteica pelo ribossomo.
    """
    residues: List[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        amino = CODON_TABLE_STANDARD.get(codon, "X")
        if amino == "*":
            break
        residues.append(amino)
    return "".join(residues)


def find_orfs(seq: str, min_length: int = 100) -> List[dict]:
    """Identifica ORFs nos seis reading frames de uma sequencia de DNA.

    Varre os tres frames da fita sense (+1, +2, +3) e os tres frames do
    complemento reverso (-1, -2, -3), usando a tabela de codons NCBI Standard #1
    (ATG inicia; TAA, TAG e TGA terminam). Cada ORF vai de um ATG ao primeiro
    codon de parada no mesmo frame.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        min_length: Comprimento minimo da ORF em pares de base, incluindo o
            codon de parada; deve ser positivo. ORFs menores sao descartadas.

    Returns:
        Lista de dicionarios ordenada por "length_bp" decrescente, cada um com
        "frame" (int: 1, 2, 3, -1, -2, -3), "strand" ("+" ou "-"), "start" e
        "end" na fita analisada (base 0, fim exclusivo), "input_start" e
        "input_end" na fita de entrada, "length_bp" (int), "start_codon" (str),
        "stop_codon" (str), "protein" (str) e "status" ("PREDICTED"). Frames
        negativos usam o complemento reverso para start/end; input_* ja estao
        mapeados para a fita colada.

    Raises:
        ValueError: Se min_length nao for positivo.

    Nota biologica:
        Como o DNA de fita dupla e bidirecional, um gene pode estar em qualquer
        das duas fitas e em qualquer um dos tres frames; por isso os seis frames
        sao varridos para encontrar todas as ORFs possiveis.
    """
    if min_length <= 0:
        raise ValueError("min_length deve ser positivo.")
    cleaned = _clean(seq)
    rc = reverse_complement(cleaned)

    orfs: List[dict] = []
    strands = [(1, cleaned), (-1, rc)]
    for sign, strand_seq in strands:
        for offset in range(3):
            frame = sign * (offset + 1)
            i = offset
            limit = len(strand_seq) - 2
            while i < limit:
                if strand_seq[i : i + 3] == START_CODON:
                    j = i
                    found_stop = False
                    while j < limit:
                        codon = strand_seq[j : j + 3]
                        if codon in STOP_CODONS:
                            end = j + 3
                            length_bp = end - i
                            if length_bp >= min_length:
                                orfs.append(
                                    {
                                        "frame": frame,
                                        "strand": "+" if sign > 0 else "-",
                                        "start": i,
                                        "end": end,
                                        "input_start": (
                                            i if sign > 0 else len(cleaned) - end
                                        ),
                                        "input_end": (
                                            end if sign > 0 else len(cleaned) - i
                                        ),
                                        "length_bp": length_bp,
                                        "start_codon": START_CODON,
                                        "stop_codon": codon,
                                        "protein": _translate(strand_seq[i:end]),
                                        "status": "PREDICTED",
                                    }
                                )
                            found_stop = True
                            break
                        j += 3
                    i += 3
                else:
                    i += 3
    orfs.sort(key=lambda orf: orf["length_bp"], reverse=True)
    return orfs


def select_non_overlapping_orfs(
    orfs: List[dict],
    min_length: int = MIN_CDS_LENGTH_NT,
    max_overlap_fraction: float = MAX_ORF_OVERLAP_FRACTION,
) -> List[dict]:
    """Seleciona ORFs nao redundantes, descartando as contidas em ORFs maiores.

    Aplica um corte de comprimento minimo e, em seguida, uma selecao gulosa da
    maior para a menor: uma ORF e descartada quando sobrepoe uma ORF maior ja
    aceita em mais de max_overlap_fraction do seu proprio comprimento. A
    comparacao ocorre dentro da mesma fita, porque frames positivos e negativos
    usam sistemas de coordenadas distintos em find_orfs.

    Args:
        orfs: Lista de ORFs no formato retornado por find_orfs(), com as chaves
            "frame", "start", "end" e "length_bp".
        min_length: Comprimento minimo em nucleotideos para manter uma ORF; deve
            ser positivo.
        max_overlap_fraction: Fracao maxima de sobreposicao tolerada, entre 0.0 e
            1.0.

    Returns:
        Lista de ORFs mantidas, ordenada por "length_bp" decrescente. Os
        dicionarios originais sao reaproveitados sem copia.

    Raises:
        ValueError: Se min_length nao for positivo ou se max_overlap_fraction
            estiver fora do intervalo [0.0, 1.0].

    Nota biologica:
        find_orfs() reporta uma ORF para cada ATG a montante de um mesmo codon de
        parada, de modo que um unico gene gera muitas ORFs aninhadas. Contar
        todas elas duplica os mesmos nucleotideos e distorce qualquer estatistica
        derivada, como uso de codons ou CAI. Manter apenas a ORF mais longa de
        cada regiao aproxima a nocao de um CDS por locus.

        Implementacao heuristica simplificada inspirada nos criterios de selecao
        de ORFs do ORFfinder (NCBI) e do getorf (EMBOSS); nao reproduz o
        modelo/algoritmo original publicado.
    """
    if min_length <= 0:
        raise ValueError("min_length deve ser positivo.")
    if not 0.0 <= max_overlap_fraction <= 1.0:
        raise ValueError("max_overlap_fraction deve estar entre 0.0 e 1.0.")

    candidates = [
        orf for orf in orfs if int(orf.get("length_bp", 0)) >= min_length
    ]
    candidates.sort(
        key=lambda orf: (
            -int(orf["length_bp"]),
            int(orf["start"]),
            int(orf["frame"]),
        )
    )

    kept: List[dict] = []
    kept_spans: List[tuple[int, int, int]] = []
    for orf in candidates:
        start = int(orf["start"])
        end = int(orf["end"])
        strand = 1 if int(orf["frame"]) > 0 else -1
        span = end - start
        if span <= 0:
            continue
        redundant = False
        for other_start, other_end, other_strand in kept_spans:
            if other_strand != strand:
                continue
            overlap = min(end, other_end) - max(start, other_start)
            if overlap > 0 and (overlap / span) > max_overlap_fraction:
                redundant = True
                break
        if not redundant:
            kept.append(orf)
            kept_spans.append((start, end, strand))

    kept.sort(key=lambda orf: int(orf["length_bp"]), reverse=True)
    return kept


def cpg_islands(
    seq: str,
    window: int = 200,
    step: int = 50,
    gc_threshold: float = 0.50,
    oe_threshold: float = 0.60,
) -> List[dict]:
    """Identifica ilhas CpG pelos criterios de Gardiner-Garden e Frommer (1987).

    Desliza uma janela ao longo da sequencia e seleciona as que satisfazem,
    simultaneamente, conteudo GC acima de gc_threshold e razao CpG observado/
    esperado (O/E) acima de oe_threshold. A razao e calculada como
    (CpG_count / (C_count * G_count)) * comprimento_da_janela. Janelas
    aprovadas que se sobrepoem sao mescladas em uma unica regiao.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        window: Tamanho da janela em bases; deve ser positivo.
        step: Passo de deslocamento da janela em bases; deve ser positivo.
        gc_threshold: Fracao minima de conteudo GC (0 a 1) para aprovar a janela.
        oe_threshold: Razao CpG O/E minima para aprovar a janela.

    Returns:
        Lista de dicionarios, um por regiao mesclada, com "start" (int), "end"
        (int), "gc_percent" (float, % com 2 casas) e "cpg_oe" (float com 2
        casas). Lista vazia quando nenhuma janela e aprovada ou quando a
        sequencia e menor que a janela.

    Raises:
        ValueError: Se window ou step nao forem positivos.

    Nota biologica:
        Ilhas CpG sao regioes ricas em GC com dinucleotideos CpG perto do
        esperado, frequentemente associadas a promotores; a baixa frequencia de
        CpG no restante do genoma decorre da metilacao e desaminacao de
        citosinas metiladas.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window e step devem ser positivos.")
    cleaned = _clean(seq)
    if len(cleaned) < window:
        return []

    passing: List[tuple[int, int]] = []
    for start in range(0, len(cleaned) - window + 1, step):
        end = start + window
        chunk = cleaned[start:end]
        c = chunk.count("C")
        g = chunk.count("G")
        gc_fraction = (c + g) / window
        if c == 0 or g == 0:
            continue
        cpg = chunk.count("CG")
        oe = (cpg / (c * g)) * window
        if gc_fraction >= gc_threshold and oe >= oe_threshold:
            passing.append((start, end))

    if not passing:
        return []

    merged: List[tuple[int, int]] = []
    current_start, current_end = passing[0]
    for start, end in passing[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))

    islands: List[dict] = []
    for start, end in merged:
        region = cleaned[start:end]
        c = region.count("C")
        g = region.count("G")
        length = end - start
        gc_percent = round(((c + g) / length) * 100.0, 2) if length else 0.0
        cpg = region.count("CG")
        oe = round((cpg / (c * g)) * length, 2) if c and g else 0.0
        islands.append(
            {
                "start": start,
                "end": end,
                "gc_percent": gc_percent,
                "cpg_oe": oe,
                "cpg_count": cpg,
                "length_bp": length,
            }
        )
    return islands


def cpg_dinucleotide_positions(seq: str) -> List[int]:
    """Lista posicoes 0-based de dinucleotideos CG na fita informada.

    Args:
        seq: Sequencia de DNA (sera normalizada; U e tratado como T).

    Returns:
        Lista crescente de indices i tais que a sequencia limpa tem "CG" em
        [i:i+2]. Toda posicao satisfaz 0 <= i < len(seq)-1.

    Raises:
        Nenhum.

    Nota biologica:
        Esta e a contagem direta de CpG na fita colada, nao uma predicao de
        ilha CpG. Densidade = len(positions) / max(len-1, 1) sobre a fita.
    """
    cleaned = _clean(seq).replace("U", "T")
    positions: List[int] = []
    for index in range(len(cleaned) - 1):
        if cleaned[index : index + 2] == "CG":
            positions.append(index)
    return positions


def shannon_entropy(seq: str) -> float:
    """Entropia de Shannon (bits) das bases canonicas A, C, G e T.

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        Entropia em bits com 4 casas. float('nan') quando nao ha A/C/G/T
        (desconhecido, nao zero). Maximo teorico para 4 simbolos
        equiprovaveis: 2.0. Uma sequencia so de A tem entropia 0.0 real.

    Raises:
        Nenhum.

    Nota biologica:
        Entropia baixa indica composicao enviesada ou baixa complexidade; nao e
        um detector de repeats publicado (RepeatMasker/Dust).
    """
    cleaned = _clean(seq)
    counts = [cleaned.count(base) for base in "ACGT"]
    total = sum(counts)
    if total == 0:
        return float("nan")
    entropy = 0.0
    for count in counts:
        if count:
            p = count / total
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def at_skew(seq: str) -> float:
    """Calcula o AT skew global (A - T) / (A + T).

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        AT skew com 4 casas, ou float('nan') se A+T for zero.

    Raises:
        Nenhum.

    Nota biologica:
        Complementa o GC skew na descricao de vies composicional entre fitas.
    """
    cleaned = _clean(seq)
    a = cleaned.count("A")
    t = cleaned.count("T")
    if a + t == 0:
        return float("nan")
    return round((a - t) / (a + t), 4)


def gc_skew_global(seq: str) -> float:
    """Calcula o GC skew global (G - C) / (G + C).

    Args:
        seq: Sequencia de DNA (sera normalizada).

    Returns:
        GC skew com 4 casas, ou float('nan') se G+C for zero.

    Raises:
        Nenhum.

    Nota biologica:
        Denominador zero nao e vies nulo: uma sequencia so de A/T nao define
        GC skew.
    """
    cleaned = _clean(seq)
    g = cleaned.count("G")
    c = cleaned.count("C")
    if g + c == 0:
        return float("nan")
    return round((g - c) / (g + c), 4)


def cumulative_skew_summary(seq: str, plot_step: Optional[int] = None) -> dict:
    """Soma cumulativa de GC e AT no sentido de Lobry, sobre a sequencia inteira.

    Cada G soma +1 ao cumulativo GC e cada C soma -1. Cada A soma +1 ao
    cumulativo AT e cada T soma -1. As outras letras somam 0 e continuam no
    eixo. O valor final usa todas as bases. A serie devolvida para grafico
    pode ter passo maior que 1; esse passo nao descarta bases do somatorio.

    Args:
        seq: Sequencia de DNA (sera normalizada).
        plot_step: Passo dos pontos visualizados. None escolhe o menor passo
            que cabe no teto tecnico de plotagem, incluindo o ultimo residuo.

    Returns:
        Dict com input_length, analyzed_length, visualized_length, plot_step,
        visualization, method, final_gc, final_at e series. final_gc e final_at
        sao inteiros. series traz position (1-based) e os dois cumulativos.

    Raises:
        ValueError: Se a sequencia ficar vazia, se plot_step nao for positivo,
            ou se o passo pedido exceder o teto de pontos.

    Nota biologica:
        Este cumulativo e a soma de Lobry (1996), nao o quociente (G-C)/(G+C).
        O quociente global continua em gc_skew_global e at_skew. Origem de
        replicacao nao e declarada a partir deste numero.
    """
    cleaned = _clean(seq)
    length = len(cleaned)
    if length == 0:
        raise ValueError("Sequence is empty; cumulative skew is not defined as 0.")
    cap = scale_profile.MAX_PLOT_POINTS
    if plot_step is None:
        step = 1 if length <= cap else (length + cap - 1) // cap
        automatic = True
    else:
        if isinstance(plot_step, bool) or not isinstance(plot_step, int) or plot_step < 1:
            raise ValueError("plot_step must be a positive integer.")
        step = plot_step
        automatic = False
    estimated = (length + step - 1) // step
    if length % step != 0:
        estimated += 0
    if estimated > cap + 1:
        raise ValueError(
            f"Cumulative skew plot would draw {estimated:,} points "
            f"(step={step}). Technical plot limit is {cap:,}. "
            "Increase plot_step. The running sum itself is not changed and "
            "is not replaced by a shorter sequence."
        )
    gc_total = 0
    at_total = 0
    series: List[dict] = []
    for index, base in enumerate(cleaned, start=1):
        if base == "G":
            gc_total += 1
        elif base == "C":
            gc_total -= 1
        elif base == "A":
            at_total += 1
        elif base == "T":
            at_total -= 1
        if index % step == 0 or index == length:
            if not series or int(series[-1]["position"]) != index:
                series.append(
                    {
                        "position": index,
                        "cumulative_gc": gc_total,
                        "cumulative_at": at_total,
                    }
                )
    if len(series) > cap + 1:
        raise ValueError(
            f"Cumulative skew plot produced {len(series):,} points. "
            f"Technical plot limit is {cap:,}."
        )
    from . import scale_contract

    contract = scale_contract.length_contract(
        input_length=length,
        analyzed_length=length,
        visualized_length=len(series),
        visualization="FULL" if step == 1 and len(series) == length else "SAMPLED",
        full_scale=True,
    )
    return {
        **contract,
        "plot_step": step,
        "plot_step_automatic": automatic,
        "method": (
            "Lobry cumulative skew: running sum of +1 for G and -1 for C, "
            "and +1 for A and -1 for T. Other symbols add 0. "
            "This is not (G-C)/(G+C)."
        ),
        "final_gc": gc_total,
        "final_at": at_total,
        "global_gc_skew_ratio": gc_skew_global(cleaned),
        "global_at_skew_ratio": at_skew(cleaned),
        "series": series,
    }


def kmer_counts(seq: str, k: int = 3) -> dict:
    """Conta k-mers canonicos (alfabeto ACGT) com deslocamento de 1.

    Args:
        seq: Sequencia de DNA.
        k: Comprimento do k-mer; 1 a 5 inclusive.

    Returns:
        Dicionario k-mer -> {"count": int, "frequency": float %} relativo ao
        numero de janelas (len - k + 1) que so contem ACGT. Janelas com N ou
        outro simbolo sao ignoradas.

    Raises:
        ValueError: Se k estiver fora de 1 a 5.

    Nota biologica:
        Frequencias de k-mers sao contagens diretas da sequencia, nao um modelo
        de selecao.
    """
    if k < 1 or k > MAX_KMER_K:
        raise ValueError(
            f"k deve estar entre 1 e {MAX_KMER_K}. k={k} excederia o limite de "
            "memoria desta analise."
        )
    cleaned = _clean(seq)
    windows = max(0, len(cleaned) - k + 1)
    counts: Dict[str, int] = {}
    valid = 0
    for i in range(windows):
        mer = cleaned[i : i + k]
        if set(mer) - CANONICAL_DNA_ALPHABET:
            continue
        counts[mer] = counts.get(mer, 0) + 1
        valid += 1
    result: dict = {}
    for mer, count in sorted(counts.items()):
        result[mer] = {
            "count": count,
            "frequency": round((count / valid) * 100.0, 4) if valid else 0.0,
        }
    return result


def kmer_summary(seq: str, k: int = 3, top_n: int = 20) -> dict:
    """Resume k-mers: ranking, diversidade observada e janelas validas.

    Args:
        seq: Sequencia de DNA.
        k: Comprimento do k-mer; 1 a MAX_KMER_K.
        top_n: Quantos k-mers mais abundantes listar; deve ser positivo.

    Returns:
        Dicionario com "k" (int), "valid_windows" (int), "unique_kmers" (int),
        "alphabet_size" (4), "possible_kmers" (int, 4**k), "diversity"
        (unique / min(valid, 4**k), ou NaN se valid=0), "top" (lista de
        {"kmer", "count", "frequency"} ordenada por count) e "status".

    Raises:
        ValueError: Se k ou top_n forem invalidos.

    Nota biologica:
        Diversidade proxima de 1.0 significa que os k-mers observados ocupam
        quase todo o alfabeto possivel dado o comprimento. Nao e complexidade
        Kolmogorov nem um detector de repeats.
    """
    if top_n <= 0:
        raise ValueError("top_n deve ser positivo.")
    counts = kmer_counts(seq, k)
    valid = sum(int(item["count"]) for item in counts.values())
    unique = len(counts)
    possible = 4 ** k
    if valid == 0:
        diversity = float("nan")
        status = "UNAVAILABLE"
    else:
        diversity = round(unique / min(valid, possible), 4)
        status = "COMPUTED"
    ranked = sorted(
        counts.items(), key=lambda item: int(item[1]["count"]), reverse=True
    )[:top_n]
    return {
        "k": k,
        "valid_windows": valid,
        "unique_kmers": unique,
        "alphabet_size": 4,
        "possible_kmers": possible,
        "diversity": diversity,
        "top": [
            {
                "kmer": mer,
                "count": int(data["count"]),
                "frequency": float(data["frequency"]),
            }
            for mer, data in ranked
        ],
        "status": status,
        "window_size": k,
        "step_size": 1,
        "alphabet": "ACGT",
    }


def entropy_sliding_window(seq: str, window: int = 100, step: int = 50) -> List[dict]:
    """Entropia de Shannon em janela deslizante.

    Args:
        seq: Sequencia de DNA.
        window: Tamanho da janela; deve ser positivo.
        step: Passo; deve ser positivo.

    Returns:
        Lista de dicts com "start" (1-based), "end" (inclusivo, 1-based),
        "midpoint" (1-based) e "entropy" (float). Lista vazia se a sequencia
        for menor que a janela.

    Raises:
        ValueError: Se window ou step nao forem positivos.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window e step devem ser positivos.")
    cleaned = _clean(seq)
    if len(cleaned) < window:
        return []
    scale_profile.enforce_window_budget(len(cleaned), window, step)
    profile: List[dict] = []
    for start in range(0, len(cleaned) - window + 1, step):
        chunk = cleaned[start : start + window]
        profile.append(
            {
                "start": start + 1,
                "end": start + window,
                "midpoint": start + (window // 2) + 1,
                "entropy": shannon_entropy(chunk),
            }
        )
    return profile


def composition_is_consistent(seq: str) -> bool:
    """Verifica se a soma das contagens de nucleotide_composition iguala o comprimento.

    Args:
        seq: Sequencia de DNA.

    Returns:
        True quando a soma dos counts bate com len(_clean(seq)).

    Raises:
        Nenhum.
    """
    cleaned = _clean(seq)
    composition = nucleotide_composition(cleaned)
    total = sum(int(item["count"]) for item in composition.values())
    return total == len(cleaned)


def supported_iupac_dna_symbols() -> Dict[str, str]:
    """Lista os simbolos IUPAC de DNA que o HelixScope realmente trata.

    Args:
        Nenhum.

    Returns:
        Mapa simbolo -> bases representadas. Inclui A, C, G, T e as
        ambiguidades de IUPAC_DNA_AMBIGUITY. Nao inclui U como simbolo de DNA
        autonomo (U e mapeado para T em algumas formulas).

    Raises:
        Nenhum.

    Nota:
        Isto nao e o alfabeto IUPAC completo de aminoacidos nem de RNA.
    """
    mapping: Dict[str, str] = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
    }
    mapping.update(IUPAC_DNA_AMBIGUITY)
    return mapping


def windowed_profiles(
    seq: str,
    window: int = 100,
    step: int = 50,
) -> List[dict]:
    """Perfis GC, AT, GC skew, AT skew e entropia na mesma grade de janelas.

    Args:
        seq: Sequencia de DNA.
        window: Tamanho da janela; deve ser positivo.
        step: Passo; deve ser positivo.

    Returns:
        Lista de dicts com "start" (0-based), "end" (exclusivo), "midpoint",
        "window_size", "step_size", "alphabet", "gc_percent", "at_percent",
        "gc_skew", "at_skew" e "entropy". Metricas indefinidas sao NaN.
        Lista vazia se a sequencia for menor que a janela.

    Raises:
        ValueError: Se window ou step nao forem positivos.

    Nota biologica:
        Todas as metricas desta grade usam o alfabeto canonico ACGT. Janelas
        so de N nao produzem 0% de GC.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window e step devem ser positivos.")
    cleaned = _clean(seq).replace("U", "T")
    if len(cleaned) < window:
        return []
    scale_profile.enforce_window_budget(len(cleaned), window, step)
    profile: List[dict] = []
    for start in range(0, len(cleaned) - window + 1, step):
        end = start + window
        chunk = cleaned[start:end]
        a = chunk.count("A")
        t = chunk.count("T")
        g = chunk.count("G")
        c = chunk.count("C")
        canonical = a + t + g + c
        gc_percent = (
            round(((g + c) / canonical) * 100.0, 2) if canonical else float("nan")
        )
        at_percent = (
            round(((a + t) / canonical) * 100.0, 2) if canonical else float("nan")
        )
        gc_skew = round((g - c) / (g + c), 4) if (g + c) else float("nan")
        at_skew_value = round((a - t) / (a + t), 4) if (a + t) else float("nan")
        profile.append(
            {
                "start": start,
                "end": end,
                "midpoint": start + window // 2,
                "window_size": window,
                "step_size": step,
                "alphabet": "ACGT",
                "gc_percent": gc_percent,
                "at_percent": at_percent,
                "gc_skew": gc_skew,
                "at_skew": at_skew_value,
                "entropy": shannon_entropy(chunk),
            }
        )
    return profile


def parse_sequence_payload(text: str) -> dict:
    """Interpreta texto colado ou conteudo de arquivo como uma unica sequencia.

    FASTA com um registro devolve so os residuos. FASTA com varios registros e
    recusado: o HelixScope nao escolhe silenciosamente o primeiro.

    Args:
        text: Texto bruto (FASTA ou sequencia continua).

    Returns:
        Dicionario com "sequence" (str), "record_count" (int), "format"
        ("empty", "fasta" ou "raw") e "identifier" (str do cabecalho FASTA).

    Raises:
        ValueError: Se o texto exceder os limites tecnicos, se o FASTA tiver
            cabecalho sem registros, ou se houver mais de um registro FASTA.

    Nota biologica:
        Concatenar registros FASTA distintos inventaria uma molecula que nao
        existe. Analises devem recair sobre uma sequencia identificavel.
    """
    payload = text or ""
    size_reason = _size_rejection_reason(payload)
    if size_reason:
        raise ValueError(size_reason)
    stripped = payload.strip()
    if not stripped:
        return {
            "sequence": "",
            "record_count": 0,
            "format": "empty",
            "identifier": "",
        }
    if stripped.startswith(">"):
        try:
            from Bio import SeqIO
        except ImportError as exc:
            raise RuntimeError(
                "Biopython nao esta instalado; instale a dependencia 'biopython'."
            ) from exc
        records = list(SeqIO.parse(io.StringIO(payload), "fasta"))
        if not records:
            raise ValueError(
                "FASTA header found but no sequence records could be parsed."
            )
        if len(records) > 1:
            raise ValueError(
                f"FASTA contains {len(records)} records. HelixScope analyzes one "
                "sequence at a time. Submit a single-record FASTA or paste one "
                "sequence. The first record was not selected automatically."
            )
        record = records[0]
        return {
            "sequence": str(record.seq),
            "record_count": 1,
            "format": "fasta",
            "identifier": str(record.id or ""),
        }
    return {
        "sequence": stripped,
        "record_count": 1,
        "format": "raw",
        "identifier": "",
    }


def parse_fasta_records(text: str, *, max_records: int = 50) -> dict:
    """Interpreta FASTA com um ou varios registros, sem escolher o primeiro.

    Analises de uma sequencia continuam usando parse_sequence_payload, que
    recusa multi-FASTA. Esta funcao existe para colecoes (MSA).

    Args:
        text: FASTA ou uma sequencia crua sem cabecalho.
        max_records: Teto de registros aceitos (limite tecnico, nao biologico).

    Returns:
        Dict com "format", "record_count" e "records" (lista de identifier,
        sequence, length). FASTA vazio devolve record_count 0.

    Raises:
        TypeError: Se text nao for string.
        ValueError: Limite tecnico, FASTA sem registros, ou excesso de records.

    Nota biologica:
        Cada registro FASTA e uma sequencia distinta. Concatenar ou descartar
        silenciosamente o restante inventaria ou omitiria moleculas.
    """
    if not isinstance(text, str):
        raise TypeError("A sequencia deve ser uma string.")
    try:
        limit = int(max_records)
    except (TypeError, ValueError) as exc:
        raise ValueError("max_records must be an integer.") from exc
    if limit < 1:
        raise ValueError("max_records must be at least 1.")
    size_reason = _size_rejection_reason(text or "")
    if size_reason:
        raise ValueError(size_reason)
    stripped = (text or "").strip()
    if not stripped:
        return {"format": "empty", "record_count": 0, "records": []}
    if not stripped.startswith(">"):
        cleaned = _clean(stripped)
        return {
            "format": "raw",
            "record_count": 1,
            "records": [
                {
                    "identifier": "",
                    "sequence": cleaned,
                    "length": len(cleaned),
                }
            ],
        }
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not records:
        raise ValueError(
            "FASTA header found but no sequence records could be parsed."
        )
    if len(records) > limit:
        raise ValueError(
            f"FASTA contains {len(records)} records; the limit is {limit}. "
            "The first record was not selected automatically."
        )
    parsed: List[dict] = []
    for record in records:
        seq = str(record.seq or "")
        parsed.append(
            {
                "identifier": str(record.id or ""),
                "sequence": seq,
                "length": len(seq),
            }
        )
    return {
        "format": "fasta",
        "record_count": len(parsed),
        "records": parsed,
    }
