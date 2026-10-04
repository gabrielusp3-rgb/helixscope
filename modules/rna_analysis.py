"""Analise de RNA: transcricao, uso de codons e indice de adaptacao de codons.

Inclui transcricao de DNA em mRNA, resolucao das regioes codificadoras a partir
de anotacao oficial ou de predicao de ORFs, tabela de uso de codons e o calculo do
Codon Adaptation Index (CAI) para humano, E. coli e levedura.

Toda contagem de codons deve partir de regioes codificadoras resolvidas por
resolve_coding_codons(); fatiar a sequencia bruta de tres em tres nucleotideos
produz numeros sem significado biologico e infla as estatisticas.

Os dados externos (codigo genetico NCBI Standard #1 em RNA, tabelas de uso de
codons por organismo em frequencia por mil, padroes de numero de acesso e notas
sobre genomas com mecanismos especiais de traducao) estao embutidos como
dicionarios constantes neste modulo. Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import math
import re
from typing import Dict, List, Optional

import pandas as pd

from . import dna_analysis

RNA_CODON_TABLE: Dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
"""Codigo genetico NCBI Standard #1 em codons de RNA; o asterisco indica parada."""

CODON_USAGE_PER_1000: Dict[str, Dict[str, float]] = {
    "human": {
        "UUU": 17.6, "UUC": 20.3, "UUA": 7.7, "UUG": 12.9,
        "CUU": 13.2, "CUC": 19.6, "CUA": 7.2, "CUG": 39.6,
        "AUU": 16.0, "AUC": 20.8, "AUA": 7.5, "AUG": 22.0,
        "GUU": 11.0, "GUC": 14.5, "GUA": 7.1, "GUG": 28.1,
        "UCU": 15.2, "UCC": 17.7, "UCA": 12.2, "UCG": 4.4,
        "CCU": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
        "ACU": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
        "GCU": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
        "UAU": 12.2, "UAC": 15.3, "UAA": 1.0, "UAG": 0.8,
        "CAU": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
        "AAU": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
        "GAU": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
        "UGU": 10.6, "UGC": 12.6, "UGA": 1.6, "UGG": 13.2,
        "CGU": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
        "AGU": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
        "GGU": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
    },
    "ecoli": {
        "UUU": 22.4, "UUC": 16.6, "UUA": 13.9, "UUG": 13.7,
        "CUU": 11.0, "CUC": 11.0, "CUA": 3.9, "CUG": 52.6,
        "AUU": 30.5, "AUC": 25.1, "AUA": 4.4, "AUG": 27.9,
        "GUU": 18.3, "GUC": 15.3, "GUA": 10.9, "GUG": 26.4,
        "UCU": 8.5, "UCC": 8.6, "UCA": 7.2, "UCG": 8.9,
        "CCU": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.2,
        "ACU": 9.0, "ACC": 23.4, "ACA": 7.1, "ACG": 14.4,
        "GCU": 15.3, "GCC": 25.5, "GCA": 20.1, "GCG": 33.6,
        "UAU": 16.2, "UAC": 12.2, "UAA": 2.0, "UAG": 0.3,
        "CAU": 12.9, "CAC": 9.7, "CAA": 15.3, "CAG": 28.8,
        "AAU": 17.7, "AAC": 21.7, "AAA": 33.6, "AAG": 10.3,
        "GAU": 32.1, "GAC": 19.1, "GAA": 39.4, "GAG": 17.8,
        "UGU": 5.2, "UGC": 6.4, "UGA": 1.0, "UGG": 15.2,
        "CGU": 20.9, "CGC": 22.0, "CGA": 3.6, "CGG": 5.4,
        "AGU": 8.8, "AGC": 16.1, "AGA": 2.1, "AGG": 1.2,
        "GGU": 24.8, "GGC": 29.6, "GGA": 8.0, "GGG": 11.1,
    },
    "yeast": {
        "UUU": 26.1, "UUC": 18.4, "UUA": 26.2, "UUG": 27.2,
        "CUU": 12.3, "CUC": 5.4, "CUA": 13.4, "CUG": 10.5,
        "AUU": 30.1, "AUC": 17.2, "AUA": 17.8, "AUG": 20.9,
        "GUU": 22.1, "GUC": 11.8, "GUA": 11.8, "GUG": 10.8,
        "UCU": 23.5, "UCC": 14.2, "UCA": 18.7, "UCG": 8.6,
        "CCU": 13.5, "CCC": 6.8, "CCA": 18.3, "CCG": 5.3,
        "ACU": 20.3, "ACC": 12.7, "ACA": 17.8, "ACG": 8.0,
        "GCU": 21.2, "GCC": 12.6, "GCA": 16.2, "GCG": 6.2,
        "UAU": 18.8, "UAC": 14.8, "UAA": 1.1, "UAG": 0.5,
        "CAU": 13.6, "CAC": 7.8, "CAA": 27.3, "CAG": 12.1,
        "AAU": 35.7, "AAC": 24.8, "AAA": 41.9, "AAG": 30.8,
        "GAU": 37.6, "GAC": 20.2, "GAA": 45.6, "GAG": 19.2,
        "UGU": 8.1, "UGC": 4.8, "UGA": 0.7, "UGG": 10.4,
        "CGU": 6.4, "CGC": 2.6, "CGA": 3.0, "CGG": 1.7,
        "AGU": 14.2, "AGC": 9.8, "AGA": 21.3, "AGG": 9.2,
        "GGU": 23.9, "GGC": 9.8, "GGA": 10.9, "GGG": 6.0,
    },
}
"""Tabelas de uso de codons (frequencia por mil) para humano, E. coli e levedura.
Valores aproximados derivados de tabelas tipo Kazusa, usados como referencia
relativa para o CAI."""

CAI_TABLE_SOURCE: Dict[str, str] = {
    "human": (
        "Approximate Kazusa-style codon usage (frequency per thousand) for "
        "Homo sapiens, embedded as a HelixScope constant. Not a live Kazusa query."
    ),
    "ecoli": (
        "Approximate Kazusa-style codon usage (frequency per thousand) for "
        "Escherichia coli, embedded as a HelixScope constant. Not a live Kazusa query."
    ),
    "yeast": (
        "Approximate Kazusa-style codon usage (frequency per thousand) for "
        "Saccharomyces cerevisiae, embedded as a HelixScope constant. Not a live "
        "Kazusa query."
    ),
}
"""Origem declarada das tabelas de referencia do CAI."""

MAX_KMER_K: int = 5
"""Maior k aceito em kmer_counts de RNA."""

SINGLE_CODON_AMINO_ACIDS: frozenset[str] = frozenset({"M", "W"})
"""Aminoacidos codificados por um unico codon no codigo padrao; nao possuem
sinonimos e por isso sao excluidos de metricas de vies de codons."""

WRIGHT_FAMILY_WEIGHTS: Dict[int, int] = {2: 9, 3: 1, 4: 5, 6: 3}
"""Numero de aminoacidos em cada classe de familia sinonima usado na formula do
Numero Efetivo de Codons de Wright (1990): 9 familias de 2 codons, 1 de 3, 5 de 4
e 3 de 6."""

ENC_MIN: float = 20.0
"""Valor minimo teorico do Numero Efetivo de Codons: um unico codon por
aminoacido, ou seja, vies maximo."""

ENC_MAX: float = 61.0
"""Valor maximo teorico do Numero Efetivo de Codons: todos os codons sinonimos
usados com igual frequencia, ou seja, ausencia de vies."""

CDS_SOURCE_ANNOTATION: str = "ncbi_annotation"
"""Identificador do metodo em que as regioes codificadoras vem de anotacao oficial."""

CDS_SOURCE_PREDICTION: str = "orf_prediction"
"""Identificador do metodo em que as regioes codificadoras sao previstas por ORFs."""

CDS_SOURCE_WHOLE: str = "whole_sequence"
"""Identificador do metodo bruto, sem filtro, em que a sequencia inteira e lida no
frame +1; nao corresponde a regioes codificadoras reais."""

CDS_SOURCE_LABELS: Dict[str, str] = {
    CDS_SOURCE_ANNOTATION: "Anotacao oficial (NCBI)",
    CDS_SOURCE_PREDICTION: "Predicao de ORFs",
    CDS_SOURCE_WHOLE: "Sequencia inteira no frame +1 (sem filtro)",
}
"""Rotulo legivel de cada metodo de obtencao das regioes codificadoras."""

ACCESSION_PATTERNS: tuple = (
    re.compile(r"\b([A-Z]{2}_\d{6,9}(?:\.\d+)?)\b"),
    re.compile(r"\b([A-Z]{4}\d{8,10}(?:\.\d+)?)\b"),
    re.compile(r"\b([A-Z]{1,2}\d{5,6}(?:\.\d+)?)\b"),
)
"""Padroes de numero de acesso reconhecidos, em ordem de prioridade: RefSeq com
prefixo sublinhado, projetos WGS e acessos GenBank classicos."""

SEQUENCE_LIKE_THRESHOLD: float = 0.9
"""Fracao minima de simbolos de acido nucleico para considerar que uma linha e
sequencia, e nao um cabecalho ou identificador."""

SPECIAL_GENOME_NOTES: Dict[str, str] = {
    "NC_045512": (
        "SARS-CoV-2: ORF1ab e traduzida como poliproteina unica e depende de "
        "frameshift ribossomal -1 no sitio de escorregamento entre ORF1a e ORF1b. "
        "As proteinas maduras nsp1 a nsp16 resultam de clivagem proteolitica, "
        "portanto um ORF previsto nao equivale a uma proteina madura."
    ),
    "MN908947": (
        "SARS-CoV-2 (isolado Wuhan-Hu-1): ORF1ab depende de frameshift ribossomal "
        "-1 e gera poliproteina processada por proteases virais; a regra de um ORF "
        "por proteina nao se aplica."
    ),
    "NC_004718": (
        "SARS-CoV: ORF1ab depende de frameshift ribossomal -1 e produz poliproteina "
        "processada proteoliticamente."
    ),
    "NC_019843": (
        "MERS-CoV: ORF1ab depende de frameshift ribossomal -1 e produz poliproteina "
        "processada proteoliticamente."
    ),
    "NC_001802": (
        "HIV-1: gag-pol depende de frameshift ribossomal -1, e tat, rev e env "
        "resultam de transcritos com splicing alternativo. A predicao de ORFs em "
        "fase unica subestima o repertorio proteico real."
    ),
    "NC_002549": (
        "Zaire ebolavirus: o gene GP gera GP1,2 por edicao transcricional (insercao "
        "de adenina pela polimerase viral), de modo que mais de uma proteina deriva "
        "da mesma regiao codificadora."
    ),
    "NC_001547": (
        "Sindbis virus: a poliproteina nao estrutural depende de readthrough do "
        "codon de parada em nsP3, gerando produtos de tamanhos distintos."
    ),
}
"""Notas sobre genomas cujo mecanismo de traducao invalida a regra de um ORF por
proteina, indexadas pelo numero de acesso sem versao."""

SPECIAL_ORGANISM_NOTES: Dict[str, str] = {
    "severe acute respiratory syndrome coronavirus": (
        "Coronavirus do grupo SARS: ORF1ab depende de frameshift ribossomal -1 e de "
        "processamento proteolitico da poliproteina."
    ),
    "middle east respiratory syndrome": (
        "MERS-CoV: ORF1ab depende de frameshift ribossomal -1 e de processamento "
        "proteolitico da poliproteina."
    ),
    "human immunodeficiency virus": (
        "HIV: gag-pol depende de frameshift ribossomal -1 e varios produtos vem de "
        "splicing alternativo."
    ),
    "influenza a virus": (
        "Influenza A: genoma segmentado; os segmentos M e NS geram proteinas "
        "adicionais por splicing e PA-X depende de frameshift +1."
    ),
    "ebolavirus": (
        "Ebolavirus: o gene GP sofre edicao transcricional e origina mais de um "
        "produto proteico."
    ),
    "hepatitis b virus": (
        "Virus da hepatite B: os quadros de leitura sao sobrepostos e o genoma "
        "compacto codifica proteinas em fases distintas da mesma regiao."
    ),
}
"""Notas sobre mecanismos especiais de traducao indexadas por trecho do nome do
organismo, em minusculas."""

POLYPROTEIN_KEYWORDS: tuple = (
    "orf1ab",
    "gag-pol",
    "polyprotein",
    "poliproteina",
    "readthrough",
)
"""Termos que, presentes em nomes de gene ou produto, indicam poliproteina ou
recodificacao traducional e exigem interpretacao cuidadosa."""


def transcribe(dna_seq: str) -> str:
    """Transcreve uma fita codificante de DNA em mRNA, substituindo T por U.

    Args:
        dna_seq: Sequencia de DNA da fita codificante (sera normalizada).

    Returns:
        Sequencia de mRNA correspondente, com uracila no lugar de timina.

    Raises:
        ValueError: Se a sequencia for vazia, contiver uracila (ou seja, ja for
            RNA) ou apresentar caracteres fora do alfabeto de DNA.

    Nota biologica:
        A RNA polimerase sintetiza um transcrito identico a fita codificante,
        exceto pela troca de T por U; a fita molde e a complementar reversa.
    """
    cleaned = "".join(dna_seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia de DNA esta vazia.")
    if "U" in cleaned:
        raise ValueError("O input ja contem U; forneca DNA, nao RNA.")
    valid = set("ACGTN") | set("RYSWKMBDHV")
    invalid = set(cleaned) - valid
    if invalid:
        encontrados = ", ".join(sorted(invalid))
        raise ValueError(f"Caracteres invalidos em sequencia de DNA: {encontrados}.")
    return cleaned.replace("T", "U")


def _looks_like_sequence(text: str) -> bool:
    """Indica se uma linha parece sequencia de acido nucleico (uso interno).

    Args:
        text: Linha de texto a avaliar.

    Returns:
        True quando ao menos SEQUENCE_LIKE_THRESHOLD dos caracteres pertencem ao
        alfabeto A, C, G, T, U e N.

    Nota biologica:
        Cabecalhos FASTA e numeros de acesso contem digitos, sublinhados e pontos,
        praticamente ausentes em sequencias, o que permite separar identificadores
        de dados de sequencia sem parsear o formato inteiro.
    """
    stripped = text.strip().upper()
    if not stripped:
        return False
    nucleotides = sum(1 for char in stripped if char in "ACGTUN")
    return (nucleotides / len(stripped)) >= SEQUENCE_LIKE_THRESHOLD


def _reverse_complement_rna(seq: str) -> str:
    """Calcula o complemento reverso de uma sequencia em alfabeto de RNA (interno).

    Args:
        seq: Sequencia normalizada em maiusculas, usando U no lugar de T.

    Returns:
        Complemento reverso; simbolos desconhecidos tornam-se N.

    Nota biologica:
        Um CDS na fita negativa e transcrito a partir da fita molde, portanto sua
        leitura em codons exige inverter e complementar as coordenadas da fita
        direta.
    """
    complement = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq))


def detect_accession(text: str) -> Optional[str]:
    """Detecta um numero de acesso GenBank ou RefSeq em um texto de entrada.

    Procura apenas em cabecalhos FASTA e em linhas curtas que nao parecem
    sequencia, evitando falsos positivos dentro dos nucleotideos.

    Args:
        text: Texto bruto colado pelo usuario ou conteudo de um arquivo FASTA.

    Returns:
        O primeiro numero de acesso reconhecido, em maiusculas e preservando a
        versao quando presente, ou None quando nenhum e encontrado.

    Raises:
        Nenhum.

    Nota biologica:
        Um numero de acesso identifica de forma estavel e versionada um registro
        depositado, permitindo recuperar a anotacao oficial de CDS em vez de
        prever regioes codificadoras estatisticamente.
    """
    if not text:
        return None

    candidates: List[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            candidates.append(stripped[1:])
        elif len(stripped) <= 32 and not _looks_like_sequence(stripped):
            candidates.append(stripped)

    for candidate in candidates:
        upper = candidate.upper()
        for pattern in ACCESSION_PATTERNS:
            match = pattern.search(upper)
            if match:
                return match.group(1)
    return None


def special_genome_notes(
    accession: Optional[str] = None,
    organism: Optional[str] = None,
    product_names: Optional[List[str]] = None,
) -> List[str]:
    """Reune avisos sobre genomas com mecanismos especiais de traducao.

    Combina tres fontes de evidencia: o numero de acesso, o nome do organismo e os
    nomes de gene ou produto das regioes codificadoras.

    Args:
        accession: Numero de acesso, com ou sem versao; opcional.
        organism: Nome do organismo tal como anotado no registro; opcional.
        product_names: Nomes de gene ou produto a inspecionar; opcional.

    Returns:
        Lista de avisos sem repeticoes, na ordem em que foram identificados. Lista
        vazia quando nenhum mecanismo especial e reconhecido.

    Raises:
        Nenhum.

    Nota biologica:
        Frameshift ribossomal, readthrough de codon de parada, edicao
        transcricional e poliproteinas quebram a correspondencia entre um quadro de
        leitura e uma proteina madura. Sinalizar esses casos evita interpretar
        contagens de ORFs como contagem de genes.
    """
    notes: List[str] = []

    if accession:
        base = accession.strip().upper().split(".")[0]
        note = SPECIAL_GENOME_NOTES.get(base)
        if note:
            notes.append(note)

    if organism:
        lowered = organism.strip().lower()
        for fragment, note in SPECIAL_ORGANISM_NOTES.items():
            if fragment in lowered and note not in notes:
                notes.append(note)

    for name in product_names or []:
        lowered_name = str(name).strip().lower()
        if not lowered_name:
            continue
        for keyword in POLYPROTEIN_KEYWORDS:
            if keyword in lowered_name:
                note = (
                    f"A anotacao inclui '{name}', indicativo de poliproteina ou "
                    "recodificacao traducional: uma unica regiao codificadora "
                    "origina multiplas proteinas maduras."
                )
                if note not in notes:
                    notes.append(note)
                break

    return notes


def orfs_to_regions(orfs: List[dict], seq_length: int) -> List[dict]:
    """Converte ORFs previstas em regioes com coordenadas da fita direta.

    find_orfs() reporta frames negativos em coordenadas do complemento reverso;
    esta funcao projeta essas coordenadas de volta para a fita direta, de modo que
    todas as regioes fiquem no mesmo sistema usado pela anotacao GenBank.

    Args:
        orfs: Lista de ORFs no formato de dna_analysis.find_orfs(), com as chaves
            "frame", "start", "end" e "length_bp".
        seq_length: Comprimento total da sequencia analisada, em nucleotideos;
            deve ser positivo.

    Returns:
        Lista de regioes com "segments" (lista de tuplas (start, end) na fita
        direta), "strand" (str), "codon_start" (int), "length_nt" (int),
        "frame" (int) e "label" (str).

    Raises:
        ValueError: Se seq_length nao for positivo.

    Nota biologica:
        Unificar o sistema de coordenadas e necessario para comparar predicao e
        anotacao oficial e para somar corretamente a extensao codificante sem
        contar o mesmo nucleotideo em duas convencoes diferentes.
    """
    if seq_length <= 0:
        raise ValueError("seq_length deve ser positivo.")

    regions: List[dict] = []
    for orf in orfs:
        frame = int(orf["frame"])
        start = int(orf["start"])
        end = int(orf["end"])
        if frame > 0:
            plus_start, plus_end = start, end
            strand = "+"
        else:
            plus_start, plus_end = seq_length - end, seq_length - start
            strand = "-"
        plus_start = max(0, plus_start)
        plus_end = min(seq_length, plus_end)
        if plus_end <= plus_start:
            continue
        regions.append(
            {
                "segments": [(plus_start, plus_end)],
                "strand": strand,
                "codon_start": 1,
                "length_nt": plus_end - plus_start,
                "frame": frame,
                "label": f"ORF frame {frame:+d} ({plus_end - plus_start} nt)",
            }
        )
    return regions


def codons_from_regions(sequence: str, regions: List[dict]) -> str:
    """Concatena os codons completos das regioes codificadoras informadas.

    Para cada regiao, une os segmentos em coordenadas da fita direta, aplica o
    complemento reverso quando a fita e negativa, respeita o deslocamento
    codon_start e descarta a trinca final incompleta.

    Args:
        sequence: Sequencia completa de DNA ou RNA; T e convertido em U.
        regions: Regioes no formato produzido por orfs_to_regions() ou por
            ncbi_fetch.extract_cds_regions(), com "segments", "strand" e,
            opcionalmente, "codon_start".

    Returns:
        String com os codons concatenados em alfabeto de RNA. String vazia quando
        nenhuma regiao produz uma trinca completa.

    Raises:
        Nenhum.

    Nota biologica:
        Somente nucleotideos pertencentes a regioes codificadoras entram na
        contagem, e sempre no quadro de leitura correto. Fatiar a sequencia inteira
        de tres em tres misturaria regioes intergenicas, introns e frames errados,
        gerando um perfil de codons que nao existe na celula.
    """
    cleaned = "".join(sequence.split()).upper().replace("T", "U")
    if not cleaned or not regions:
        return ""

    codons: List[str] = []
    for region in regions:
        segments = region.get("segments") or []
        pieces: List[str] = []
        for span in segments:
            start = max(0, int(span[0]))
            end = min(len(cleaned), int(span[1]))
            if end > start:
                pieces.append(cleaned[start:end])
        joined = "".join(pieces)
        if not joined:
            continue
        if str(region.get("strand", "+")) == "-":
            joined = _reverse_complement_rna(joined)
        offset = int(region.get("codon_start", 1)) - 1
        if not 0 <= offset <= 2:
            offset = 0
        joined = joined[offset:]
        for i in range(0, len(joined) - 2, 3):
            codons.append(joined[i : i + 3])
    return "".join(codons)


def resolve_coding_codons(
    sequence: str,
    cds_regions: Optional[List[dict]] = None,
    source: str = "auto",
    min_cds_length: int = dna_analysis.MIN_CDS_LENGTH_NT,
    drop_overlapping: bool = True,
    allow_whole_sequence: bool = False,
) -> Dict[str, object]:
    """Resolve as regioes codificadoras de uma sequencia e extrai seus codons.

    A ordem de prioridade e: anotacao oficial fornecida em cds_regions; predicao
    de ORFs com criterios minimos (inicio em AUG, terminacao em codon de parada,
    comprimento minimo e remocao de sobreposicoes redundantes); e, apenas quando
    explicitamente autorizado, leitura bruta da sequencia inteira no frame +1.
    Nunca ha recurso silencioso a sequencia inteira.

    Args:
        sequence: Sequencia de DNA ou RNA a analisar.
        cds_regions: Regioes codificadoras oficiais, no formato de
            ncbi_fetch.extract_cds_regions(); opcional.
        source: Estrategia desejada: "auto", "annotation", "prediction" ou
            "whole".
        min_cds_length: Comprimento minimo, em nucleotideos, de uma ORF prevista.
        drop_overlapping: Remove ORFs previstas redundantes por sobreposicao.
        allow_whole_sequence: Autoriza a leitura bruta da sequencia inteira quando
            nenhuma regiao codificadora e encontrada.

    Returns:
        Dicionario com "method" (str), "method_label" (str), "codons" (str),
        "n_cds" (int), "n_codons" (int), "total_input_nt" (int), "coding_nt"
        (int), "coverage_pct" (float), "naive_codon_count" (int),
        "min_cds_length" (int), "warnings" (list[str]) e "regions" (list[dict]).

    Raises:
        ValueError: Se a sequencia for vazia, se source for invalido, se
            min_cds_length nao for positivo, se nenhuma regiao codificadora for
            encontrada sem autorizacao para a sequencia inteira, ou se as regioes
            encontradas nao contiverem nenhum codon completo.

    Nota biologica:
        O numero de codons deve corresponder a soma das regioes codificadoras
        reais, nunca ao comprimento do genoma dividido por tres. Um genoma pequeno
        que resulta em centenas de milhares de codons e sinal de que o filtro
        falhou e que regioes nao codificantes ou frames sobrepostos foram
        contados.

        Implementacao heuristica simplificada inspirada nos criterios de predicao
        de ORFs do ORFfinder (NCBI) e do getorf (EMBOSS); nao reproduz o
        modelo/algoritmo original publicado.
    """
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia esta vazia.")
    if min_cds_length <= 0:
        raise ValueError("min_cds_length deve ser positivo.")

    normalized_source = (source or "auto").strip().lower()
    if normalized_source not in {"auto", "annotation", "prediction", "whole"}:
        raise ValueError(
            "source deve ser 'auto', 'annotation', 'prediction' ou 'whole'."
        )

    total_nt = len(cleaned)
    dna_view = cleaned.replace("U", "T")
    warnings: List[str] = []
    regions: List[dict] = []
    method = ""

    if normalized_source in {"auto", "annotation"} and cds_regions:
        annotated = [
            region
            for region in cds_regions
            if sum(int(span[1]) - int(span[0]) for span in (region.get("segments") or []))
            > 0
        ]
        if annotated:
            regions = annotated
            method = CDS_SOURCE_ANNOTATION

    if not regions and normalized_source in {"auto", "prediction"}:
        predicted = dna_analysis.find_orfs(dna_view, min_length=min_cds_length)
        raw_count = len(predicted)
        if drop_overlapping:
            predicted = dna_analysis.select_non_overlapping_orfs(
                predicted, min_length=min_cds_length
            )
        if predicted:
            regions = orfs_to_regions(predicted, total_nt)
            method = CDS_SOURCE_PREDICTION
            if raw_count > len(predicted):
                warnings.append(
                    f"{raw_count} ORFs brutas foram reduzidas a {len(predicted)} "
                    "apos descartar quadros aninhados e sobrepostos; sem esse "
                    "filtro os mesmos nucleotideos seriam contados varias vezes."
                )

    if not regions:
        if normalized_source == "whole" or allow_whole_sequence:
            method = CDS_SOURCE_WHOLE
            regions = [
                {
                    "segments": [(0, total_nt)],
                    "strand": "+",
                    "codon_start": 1,
                    "length_nt": total_nt,
                    "label": "Sequencia inteira no frame +1",
                }
            ]
            warnings.append(
                "Leitura bruta da sequencia inteira no frame +1: os codons NAO "
                "correspondem a regioes codificadoras reais e servem apenas para "
                "inspecao exploratoria."
            )
        else:
            raise ValueError(
                "Nenhuma regiao codificadora foi identificada com o comprimento "
                f"minimo de {min_cds_length} nt. Forneca um numero de acesso com "
                "anotacao oficial, reduza o comprimento minimo ou autorize "
                "explicitamente a leitura da sequencia inteira."
            )

    codons = codons_from_regions(cleaned, regions)
    n_codons = len(codons) // 3
    if n_codons == 0:
        raise ValueError(
            "As regioes identificadas nao contem nenhum codon completo."
        )

    coding_nt = 0
    for region in regions:
        coding_nt += sum(
            int(span[1]) - int(span[0]) for span in (region.get("segments") or [])
        )
    coverage_pct = (coding_nt / total_nt) * 100.0 if total_nt else 0.0

    if n_codons * 3 > total_nt:
        warnings.append(
            f"O total de {n_codons} codons excede os {total_nt} nucleotideos "
            "disponiveis: ha regioes sobrepostas sendo contadas mais de uma vez."
        )

    if method == CDS_SOURCE_PREDICTION:
        if coverage_pct < 5.0:
            warnings.append(
                f"Apenas {coverage_pct:.2f}% da sequencia foi classificada como "
                "codificante; o resultado pode subestimar o repertorio real ou "
                "indicar que a entrada nao e uma regiao genica."
            )
        elif coverage_pct > 95.0:
            warnings.append(
                f"{coverage_pct:.2f}% da sequencia foi classificada como "
                "codificante; verifique se houve superpredicao de ORFs."
            )
        density = (len(regions) / total_nt) * 1000.0 if total_nt else 0.0
        if density > 3.0:
            warnings.append(
                f"Densidade de {density:.2f} regioes codificadoras por kb, acima "
                "do esperado mesmo para genomas compactos; revise o comprimento "
                "minimo adotado."
            )

    return {
        "method": method,
        "method_label": CDS_SOURCE_LABELS[method],
        "codons": codons,
        "n_cds": len(regions),
        "n_codons": n_codons,
        "total_input_nt": total_nt,
        "coding_nt": coding_nt,
        "coverage_pct": round(coverage_pct, 2),
        "naive_codon_count": total_nt // 3,
        "min_cds_length": min_cds_length,
        "warnings": warnings,
        "regions": regions,
    }


def _count_codons(mrna_seq: str) -> Dict[str, int]:
    """Conta os codons de uma sequencia em trincas consecutivas (uso interno).

    Args:
        mrna_seq: Sequencia de codons concatenados, normalmente a saida de
            resolve_coding_codons().

    Returns:
        Dicionario de codon para contagem, ignorando a trinca final incompleta.

    Raises:
        ValueError: Se a sequencia for vazia ou nao contiver nenhum codon
            completo.

    Nota biologica:
        A entrada precisa ser a concatenacao de regioes codificadoras em quadro de
        leitura correto; trincas lidas fora de fase nao correspondem a codons.
    """
    cleaned = "".join(mrna_seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia de mRNA esta vazia.")
    counts: Dict[str, int] = {}
    for i in range(0, len(cleaned) - 2, 3):
        codon = cleaned[i : i + 3]
        counts[codon] = counts.get(codon, 0) + 1
    if not counts:
        raise ValueError("A sequencia nao contem nenhum codon completo.")
    return counts


def synonymous_family_sizes() -> Dict[str, int]:
    """Retorna o numero de codons sinonimos de cada aminoacido no codigo padrao.

    Args:
        Nenhum.

    Returns:
        Dicionario de aminoacido (codigo de uma letra, incluindo o asterisco de
        parada) para o numero de codons que o especificam no codigo genetico
        NCBI Standard #1.

    Raises:
        Nenhum.

    Nota biologica:
        O tamanho da familia sinonima determina o valor maximo possivel de RSCU
        para os codons daquele aminoacido e define as classes de familia usadas na
        formula do Numero Efetivo de Codons.
    """
    sizes: Dict[str, int] = {}
    for amino in RNA_CODON_TABLE.values():
        sizes[amino] = sizes.get(amino, 0) + 1
    return sizes


def relative_synonymous_codon_usage(mrna_seq: str) -> pd.DataFrame:
    """Calcula o RSCU (Relative Synonymous Codon Usage) de cada codon.

    O RSCU de um codon e a sua contagem observada dividida pela contagem que seria
    esperada se todos os codons sinonimos do mesmo aminoacido fossem usados com
    igual frequencia: RSCU = n * X_i / soma(X_j), onde n e o numero de codons
    sinonimos do aminoacido. Codons de parada sao excluidos; aminoacidos com um
    unico codon (Met e Trp) recebem RSCU igual a 1.0 por definicao.

    Args:
        mrna_seq: Codons concatenados das regioes codificadoras, tal como
            produzido por resolve_coding_codons().

    Returns:
        DataFrame com as colunas "Codon", "AminoAcid", "Count", "RSCU" (float com
        3 casas) e "FamilySize" (int), ordenado por "AminoAcid" e "RSCU"
        decrescente. Inclui os codons nao observados dos aminoacidos presentes,
        com RSCU igual a 0.0.

    Raises:
        ValueError: Se a sequencia for vazia, nao contiver nenhum codon completo
            ou nao contiver nenhum aminoacido nao terminador.

    Nota biologica:
        RSCU igual a 1.0 indica ausencia de vies: o codon e usado exatamente na
        proporcao esperada entre seus sinonimos. Valores acima de 1.0 indicam
        codon preferido e abaixo de 1.0 codon evitado. Como a metrica e normalizada
        dentro de cada aminoacido, ela separa o vies de uso de codons da
        composicao de aminoacidos da proteina, ao contrario da frequencia
        absoluta.
    """
    counts = _count_codons(mrna_seq)
    sizes = synonymous_family_sizes()

    amino_totals: Dict[str, int] = {}
    for codon, count in counts.items():
        amino = RNA_CODON_TABLE.get(codon)
        if amino is None or amino == "*":
            continue
        amino_totals[amino] = amino_totals.get(amino, 0) + count

    if not amino_totals:
        raise ValueError(
            "A sequencia nao contem codons de aminoacidos nao terminadores."
        )

    rows = []
    for codon, amino in RNA_CODON_TABLE.items():
        if amino == "*" or amino not in amino_totals:
            continue
        family_size = sizes[amino]
        total = amino_totals[amino]
        count = counts.get(codon, 0)
        rscu = (family_size * count) / total if total else 0.0
        rows.append(
            {
                "Codon": codon,
                "AminoAcid": amino,
                "Count": count,
                "RSCU": round(rscu, 3),
                "FamilySize": family_size,
            }
        )

    frame = pd.DataFrame(rows, columns=["Codon", "AminoAcid", "Count", "RSCU", "FamilySize"])
    return frame.sort_values(
        by=["AminoAcid", "RSCU"], ascending=[True, False]
    ).reset_index(drop=True)


def effective_number_of_codons(mrna_seq: str) -> float:
    """Calcula o Numero Efetivo de Codons (ENC ou Nc) de Wright (1990).

    Para cada aminoacido com pelo menos dois codons observados, calcula a
    homozigosidade estimada F = (n * soma(p_i^2) - 1) / (n - 1), onde n e o total
    de codons observados para aquele aminoacido e p_i a proporcao de cada codon
    sinonimo. As estimativas sao mediadas por classe de familia sinonima e
    combinadas em Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6.

    Args:
        mrna_seq: Codons concatenados das regioes codificadoras, tal como
            produzido por resolve_coding_codons().

    Returns:
        Valor de Nc limitado ao intervalo [20.0, 61.0], com 2 casas decimais.

    Raises:
        ValueError: Se a sequencia for vazia, nao contiver nenhum codon completo,
            ou se as familias de 2, 4 ou 6 codons nao tiverem dados suficientes
            para estimar F de forma confiavel.

    Nota biologica:
        Nc igual a 61 indica que todos os codons sinonimos sao usados com igual
        frequencia, ou seja, ausencia de vies; Nc igual a 20 indica que apenas um
        codon e usado por aminoacido, ou seja, vies maximo. Genes altamente
        expressos tendem a valores baixos (cerca de 30 ou menos) porque a selecao
        traducional favorece codons com tRNAs abundantes, enquanto genes pouco
        expressos ficam proximos do esperado pela composicao de bases. Quando Nc
        acompanha o conteudo GC na terceira posicao, o vies e provavelmente
        mutacional e nao seletivo.

        A familia de 3 codons (isoleucina) e a unica de sua classe; quando ela nao
        tem dados, F3 e interpolado como a media de F2 e F4, conforme a convencao
        do artigo original.
    """
    counts = _count_codons(mrna_seq)
    sizes = synonymous_family_sizes()

    per_amino_counts: Dict[str, Dict[str, int]] = {}
    for codon, count in counts.items():
        amino = RNA_CODON_TABLE.get(codon)
        if amino is None or amino == "*" or amino in SINGLE_CODON_AMINO_ACIDS:
            continue
        per_amino_counts.setdefault(amino, {})[codon] = count

    family_estimates: Dict[int, List[float]] = {2: [], 3: [], 4: [], 6: []}
    for amino, codon_counts in per_amino_counts.items():
        family_size = sizes[amino]
        if family_size not in family_estimates:
            continue
        total = sum(codon_counts.values())
        if total < 2:
            continue
        homozygosity = sum((count / total) ** 2 for count in codon_counts.values())
        estimate = ((total * homozygosity) - 1.0) / (total - 1.0)
        if estimate <= 0.0:
            continue
        family_estimates[family_size].append(estimate)

    averages: Dict[int, float] = {}
    for family_size, estimates in family_estimates.items():
        if estimates:
            averages[family_size] = sum(estimates) / len(estimates)

    missing = [size for size in (2, 4, 6) if size not in averages]
    if missing:
        faltando = ", ".join(f"{size} codons" for size in missing)
        raise ValueError(
            "Nao ha dados suficientes para estimar o ENC de forma confiavel: as "
            f"familias sinonimas de {faltando} nao tem codons observados em "
            "quantidade minima. Analise uma regiao codificadora mais longa."
        )

    if 3 not in averages:
        averages[3] = (averages[2] + averages[4]) / 2.0

    enc = 2.0
    for family_size, weight in WRIGHT_FAMILY_WEIGHTS.items():
        enc += weight / averages[family_size]
    return round(max(ENC_MIN, min(ENC_MAX, enc)), 2)


def codon_usage_table(mrna_seq: str) -> pd.DataFrame:
    """Constroi a tabela de uso de codons de um mRNA no frame de leitura +1.

    Conta os codons da sequencia (trincas a partir da primeira base, ignorando a
    trinca final incompleta) e calcula, para cada codon, a contagem absoluta, a
    porcentagem dentro do conjunto de codons sinonimos do mesmo aminoacido
    (Synonymous_pct), a porcentagem em relacao ao total de codons (Absolute_pct) e
    o RSCU.

    Args:
        mrna_seq: Sequencia de mRNA (sera normalizada para maiusculas).

    Returns:
        DataFrame com as colunas "Codon", "AminoAcid", "Count", "Synonymous_pct",
        "Absolute_pct" e "RSCU", ordenado por "AminoAcid" e "Count" decrescente.
        Codons de parada sao movidos para o final do DataFrame e recebem RSCU
        calculado dentro da propria familia de terminadores.

    Raises:
        ValueError: Se a sequencia for vazia ou nao contiver nenhum codon
            completo.

    Nota biologica:
        A distribuicao do uso de codons sinonimos reflete o vies de codons do
        transcrito, relevante para eficiencia de traducao e expressao heterologa.
        A entrada deve ser a concatenacao de regioes codificadoras produzida por
        resolve_coding_codons(); passar um genoma bruto faz a funcao ler regioes
        nao codificantes em frame arbitrario e produz um perfil sem significado
        biologico.
    """
    counts = _count_codons(mrna_seq)
    sizes = synonymous_family_sizes()

    amino_totals: Dict[str, int] = {}
    total_codons = sum(counts.values())
    for codon, count in counts.items():
        amino = RNA_CODON_TABLE.get(codon, "X")
        amino_totals[amino] = amino_totals.get(amino, 0) + count

    columns = [
        "Codon",
        "AminoAcid",
        "Count",
        "Synonymous_pct",
        "Absolute_pct",
        "RSCU",
    ]
    rows = []
    stop_rows = []
    for codon, count in counts.items():
        amino = RNA_CODON_TABLE.get(codon, "X")
        total_for_amino = amino_totals[amino]
        family_size = sizes.get(amino, 1)
        synonymous_pct = (
            round((count / total_for_amino) * 100.0, 2) if total_for_amino else 0.0
        )
        abs_frequency = (
            round((count / total_codons) * 100.0, 2) if total_codons else 0.0
        )
        rscu = (
            round((family_size * count) / total_for_amino, 3)
            if total_for_amino
            else 0.0
        )

        row = {
            "Codon": codon,
            "AminoAcid": amino,
            "Count": count,
            "Synonymous_pct": synonymous_pct,
            "Absolute_pct": abs_frequency,
            "RSCU": rscu,
        }
        if amino == "*":
            stop_rows.append(row)
        else:
            rows.append(row)

    frame = pd.DataFrame(rows, columns=columns)
    frame = frame.sort_values(
        by=["AminoAcid", "Count"], ascending=[True, False]
    ).reset_index(drop=True)

    if stop_rows:
        stop_frame = pd.DataFrame(stop_rows, columns=columns)
        stop_frame = stop_frame.sort_values(by="Count", ascending=False)
        frame = pd.concat([frame, stop_frame], ignore_index=True)

    return frame


def codon_adaptation_index(mrna_seq: str, organism: str = "human", exclude_stops: bool = True) -> float:
    """Calcula o Codon Adaptation Index (CAI) de um mRNA para um organismo.

    O CAI e a media geometrica da adaptabilidade relativa (w) dos codons do gene,
    onde w de cada codon e a sua frequencia dividida pela frequencia do codon
    sinonimo mais usado no organismo de referencia. Codons de aminoacidos com um 
    unico codon (Met e Trp) sao desconsiderados.

    Args:
        mrna_seq: Sequencia de mRNA (sera normalizada para maiusculas).
        organism: Organismo de referencia: "human", "ecoli" ou "yeast".
        exclude_stops: Ignora codons de parada (True por padrao).

    Returns:
        CAI no intervalo [0.0, 1.0], com 4 casas decimais. Retorna float('nan')
        quando nao ha codons informativos (desconhecido, nao zero).

    Raises:
        ValueError: Se o organismo for invalido ou a sequencia for vazia.

    Nota biologica:
        O CAI estima o quao bem o uso de codons de um gene se ajusta ao de genes
        altamente expressos do hospedeiro, indicando potencial de expressao.
        Referencia: Sharp & Li (1987), Nucleic Acids Research.

        Implementacao heuristica simplificada inspirada em Sharp & Li (1987),
        Nucleic Acids Research; nao reproduz o modelo/algoritmo original
        publicado.
    """
    key = organism.strip().lower()
    if key not in CODON_USAGE_PER_1000:
        suportados = ", ".join(sorted(CODON_USAGE_PER_1000))
        raise ValueError(f"Organismo invalido. Use um de: {suportados}.")
    cleaned = "".join(mrna_seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia de mRNA esta vazia.")

    usage = CODON_USAGE_PER_1000[key]

    max_synonymous: Dict[str, float] = {}
    for codon, amino in RNA_CODON_TABLE.items():
        freq = usage.get(codon, 0.0)
        if amino not in max_synonymous or freq > max_synonymous[amino]:
            max_synonymous[amino] = freq

    single_codon_aminos = {"M", "W"}
    log_sum = 0.0
    informative = 0
    for i in range(0, len(cleaned) - 2, 3):
        codon = cleaned[i : i + 3]
        amino = RNA_CODON_TABLE.get(codon)
        if amino is None or amino in single_codon_aminos:
            continue
        if exclude_stops and amino == "*":
            continue
        reference_max = max_synonymous.get(amino, 0.0)
        freq = usage.get(codon, 0.0)
        if reference_max <= 0.0 or freq <= 0.0:
            continue
        log_sum += math.log(freq / reference_max)
        informative += 1

    if informative == 0:
        return float("nan")
    cai = math.exp(log_sum / informative)
    return round(max(0.0, min(1.0, cai)), 4)


def codon_adaptation_index_report(
    mrna_seq: str, organism: str = "human", exclude_stops: bool = True
) -> dict:
    """CAI com a tabela de referencia, o metodo e o numero de codons usados.

    Args:
        mrna_seq: Codons concatenados das regioes codificadoras.
        organism: "human", "ecoli" ou "yeast".
        exclude_stops: Ignora codons de parada.

    Returns:
        Dicionario com "value", "status", "organism", "table_source", "method",
        "n_informative" e "exclude_stops". value e None/NaN quando inaplicavel.

    Raises:
        ValueError: Se o organismo for invalido ou a sequencia for vazia.
    """
    key = organism.strip().lower()
    if key not in CODON_USAGE_PER_1000:
        suportados = ", ".join(sorted(CODON_USAGE_PER_1000))
        raise ValueError(f"Organismo invalido. Use um de: {suportados}.")
    value = codon_adaptation_index(
        mrna_seq, organism=key, exclude_stops=exclude_stops
    )
    unavailable = isinstance(value, float) and math.isnan(value)
    cleaned = "".join(mrna_seq.split()).upper()
    informative = 0
    if not unavailable:
        sizes_skip = {"M", "W"}
        for i in range(0, len(cleaned) - 2, 3):
            codon = cleaned[i : i + 3]
            amino = RNA_CODON_TABLE.get(codon)
            if amino is None or amino in sizes_skip:
                continue
            if exclude_stops and amino == "*":
                continue
            informative += 1
    return {
        "value": value,
        "status": "UNAVAILABLE" if unavailable else "HEURISTIC",
        "organism": key,
        "table_source": CAI_TABLE_SOURCE[key],
        "method": (
            "Geometric mean of relative adaptiveness w = f(codon)/f(best synonym) "
            "from the embedded frequency-per-thousand table. Heuristic inspired by "
            "Sharp and Li (1987); not the original published software."
        ),
        "n_informative": 0 if unavailable else informative,
        "exclude_stops": exclude_stops,
    }


def effective_number_of_codons_report(mrna_seq: str) -> dict:
    """ENC de Wright (1990) ou UNAVAILABLE quando as familias nao permitem F.

    Args:
        mrna_seq: Codons concatenados das regioes codificadoras.

    Returns:
        Dicionario com "value" (float ou NaN), "status", "method" e "reason".

    Raises:
        Nenhum. Sequencia vazia ou familias insuficientes viram UNAVAILABLE.
    """
    method = (
        "Wright (1990) Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, with F3 interpolated "
        "from F2 and F4 when isoleucine is missing."
    )
    try:
        value = effective_number_of_codons(mrna_seq)
    except ValueError as exc:
        return {
            "value": float("nan"),
            "status": "UNAVAILABLE",
            "method": method,
            "reason": str(exc),
        }
    return {
        "value": value,
        "status": "COMPUTED",
        "method": method,
        "reason": "",
    }


def shannon_entropy(seq: str) -> float:
    """Entropia de Shannon (bits) das bases canonicas A, C, G e U.

    Args:
        seq: Sequencia de RNA (sera normalizada).

    Returns:
        Entropia em bits com 4 casas. float('nan') quando nao ha A/C/G/U
        (desconhecido, nao zero). Maximo teorico para 4 simbolos
        equiprovaveis: 2.0.

    Raises:
        Nenhum.

    Nota biologica:
        Entropia baixa indica composicao enviesada; nao e predicao de estrutura
        secundaria.
    """
    cleaned = "".join(seq.split()).upper()
    counts = [cleaned.count(base) for base in "ACGU"]
    total = sum(counts)
    if total == 0:
        return float("nan")
    entropy = 0.0
    for count in counts:
        if count:
            p = count / total
            entropy -= p * math.log2(p)
    return round(entropy, 4)


def kmer_counts(seq: str, k: int = 3) -> dict:
    """Conta k-mers canonicos de RNA (alfabeto ACGU) com deslocamento de 1.

    Args:
        seq: Sequencia de RNA.
        k: Comprimento do k-mer; 1 a 5 inclusive.

    Returns:
        Dicionario k-mer -> {"count": int, "frequency": float %} relativo ao
        numero de janelas que so contem ACGU.

    Raises:
        ValueError: Se k estiver fora de 1 a 5.
    """
    if k < 1 or k > MAX_KMER_K:
        raise ValueError(
            f"k deve estar entre 1 e {MAX_KMER_K}. k={k} excederia o limite de "
            "memoria desta analise."
        )
    cleaned = "".join(seq.split()).upper()
    windows = max(0, len(cleaned) - k + 1)
    counts: Dict[str, int] = {}
    valid = 0
    canonical = set("ACGU")
    for i in range(windows):
        mer = cleaned[i : i + k]
        if set(mer) - canonical:
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


RNA_STRUCTURE_TOOLS: List[str] = [
    "RNAfold (ViennaRNA, energia livre minima)",
    "mfold / UNAfold",
    "RNAstructure",
]
"""Ferramentas de folding. ViennaRNA e usada quando detectada; as outras nao estao integradas."""


def secondary_structure_availability() -> dict:
    """Disponibilidade de ViennaRNA; pairwise/CAI nao sao folding.

    Args:
        Nenhum.

    Returns:
        Dict de rna_folding.tool_availability, com recommended_tools preservadas.

    Raises:
        Nenhum.

    Nota biologica:
        Predicao de folding de RNA exige ViennaRNA (ou equivalente). Sem o
        motor, o estado permanece UNAVAILABLE. Desenhar uma estrutura
        ilustrativa seria fabricacao.
    """
    from . import rna_folding

    info = rna_folding.tool_availability()
    tools = list(info.get("recommended_tools") or RNA_STRUCTURE_TOOLS)
    info["recommended_tools"] = tools
    return info
