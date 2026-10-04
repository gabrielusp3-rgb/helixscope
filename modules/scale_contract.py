"""Contrato de comprimento para analises nucleotidicas de escala completa.

Separa o que entrou, o que foi calculado e o que foi desenhado. Nao trunca
sequencia e nao calcula biologia por conta propria: os numeros cientificos
continuam nas funcoes de DNA e RNA. Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import hashlib
from typing import Mapping

FULL_SCALE_LENGTHS: tuple[int, ...] = (100_000, 125_000, 150_000)
"""Comprimentos sinteticos obrigatorios. Nenhum deles e um organismo real."""

FULL_SCALE_TARGET_NT: int = 150_000
"""Meta de entrada para analises lineares. Nao autoriza folding global."""

_REPEAT: str = "ACGT"
"""Bloco deterministico. Nao representa genoma, gene ou isolado."""


def deterministic_nucleotide_fixture(length: int) -> str:
    """Gera DNA sintetico ACGT repetido, com comprimento exato.

    Args:
        length: Numero de residuos. Deve ser positivo e nao maior que o teto
            tecnico de dna_analysis.MAX_INPUT_RESIDUES.

    Returns:
        Sequencia A, C, G, T em ciclo, sem cabecalho e sem organismo.

    Raises:
        ValueError: Se length nao for um inteiro positivo dentro do teto.

    Nota biologica:
        A repeticao ACGT e um controle de composicao, nao uma sequencia
        biologica. Quando o comprimento e multiplo de 4, cada base aparece
        o mesmo numero de vezes.
    """
    from . import dna_analysis

    if isinstance(length, bool) or not isinstance(length, int):
        raise ValueError("length must be a positive integer.")
    if length < 1 or length > dna_analysis.MAX_INPUT_RESIDUES:
        raise ValueError(
            f"length must be between 1 and {dna_analysis.MAX_INPUT_RESIDUES}."
        )
    copies = (length + len(_REPEAT) - 1) // len(_REPEAT)
    return (_REPEAT * copies)[:length]


def fixture_digest(sequence: str) -> str:
    """SHA-256 hexadecimal da sequencia ja normalizada.

    Args:
        sequence: Residuos sem espaco.

    Returns:
        Digest hexadecimal de 64 caracteres.

    Raises:
        TypeError: Se sequence nao for str.

    Nota biologica:
        O digest identifica a entrada analisada. Nao e uma anotacao biologica.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string.")
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def length_contract(
    *,
    input_length: int,
    analyzed_length: int,
    visualized_length: int,
    visualization: str,
    full_scale: bool,
) -> dict:
    """Registra as tres grandezas do contrato de escala.

    Args:
        input_length: Residuos recebidos apos normalizacao.
        analyzed_length: Residuos usados no calculo.
        visualized_length: Pontos ou residuos efetivamente desenhados.
        visualization: Nome do modo visual (FULL, WINDOW, PREVIEW, SAMPLED).
        full_scale: Quando verdadeiro, analyzed_length tem de igualar input_length.

    Returns:
        Dict com input_length, analyzed_length, visualized_length, visualization
        e full_scale.

    Raises:
        ValueError: Se algum comprimento for negativo, se a visualizacao exceder
            a entrada, ou se uma analise full-scale nao cobrir a entrada inteira.

    Nota biologica:
        Um grafico com menos pontos do que a sequencia nao encurta a molecula
        analisada. A diferenca precisa estar neste registro.
    """
    for name, value in (
        ("input_length", input_length),
        ("analyzed_length", analyzed_length),
        ("visualized_length", visualized_length),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a non-negative integer.")
    if analyzed_length > input_length or visualized_length > input_length:
        raise ValueError(
            "analyzed_length and visualized_length cannot exceed input_length."
        )
    if full_scale and analyzed_length != input_length:
        raise ValueError(
            "Full-scale analysis requires analyzed_length == input_length. "
            f"Got analyzed_length={analyzed_length}, input_length={input_length}."
        )
    mode = str(visualization or "").strip().upper()
    if mode not in {"FULL", "WINDOW", "PREVIEW", "SAMPLED", "LOCAL_WINDOW"}:
        raise ValueError(
            "visualization must be FULL, WINDOW, PREVIEW, SAMPLED, or LOCAL_WINDOW."
        )
    return {
        "input_length": input_length,
        "analyzed_length": analyzed_length,
        "visualized_length": visualized_length,
        "visualization": mode,
        "full_scale": bool(full_scale),
    }


def assert_reported_length(sequence: str, reported_length: int) -> None:
    """Confere que o comprimento reportado e o da sequencia normalizada.

    Args:
        sequence: Sequencia ja limpa.
        reported_length: Comprimento declarado pelo resultado.

    Returns:
        None.

    Raises:
        ValueError: Se os comprimentos divergirem, incluindo corte de um residuo.

    Nota biologica:
        149999 ou 150001 residuos numa entrada de 150000 e falha de integridade,
        nao um arredondamento aceitavel.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string.")
    if len(sequence) != int(reported_length):
        raise ValueError(
            f"Reported length {reported_length} does not match sequence length "
            f"{len(sequence)}."
        )


def independent_base_counts(sequence: str) -> Mapping[str, int]:
    """Conta A, C, G e T com um laco proprio deste modulo.

    Args:
        sequence: DNA ja em maiusculas.

    Returns:
        Mapa com as quatro bases. Simbolos fora de ACGT nao entram na conta.

    Raises:
        TypeError: Se sequence nao for str.

    Nota biologica:
        Esta contagem existe para o contrato de fixture. A analise de producao
        continua em dna_analysis.nucleotide_composition.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string.")
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for base in sequence:
        if base in counts:
            counts[base] += 1
    return counts
