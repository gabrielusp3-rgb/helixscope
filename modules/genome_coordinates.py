"""Conversao explicita entre coordenadas internas e de exibicao.

Interno: 0-based half-open [start, end).
Exibicao genomica: 1-based inclusive [start, end].

Nenhuma funcao inventa um contig. Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional


class CoordinateError(ValueError):
    """Intervalo ou origem invalidos.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT").strip().upper() or "INVALID_INPUT"


INTERNAL_CONVENTION: str = "0-based half-open"
DISPLAY_CONVENTION: str = "1-based inclusive"


def internal_to_display(start_0based: int, end_0based: int) -> dict:
    """Converte [start, end) interno para exibicao 1-based inclusive.

    Args:
        start_0based: Inicio interno (incluido).
        end_0based: Fim interno (excluido).

    Returns:
        Dict start_1based, end_1based, empty, convention.

    Raises:
        CoordinateError: INVALID_INPUT.
    """
    try:
        start = int(start_0based)
        end = int(end_0based)
    except (TypeError, ValueError) as exc:
        raise CoordinateError("Coordinates must be integers.") from exc
    if start < 0 or end < start:
        raise CoordinateError(
            "Internal interval must be 0-based half-open with end >= start >= 0."
        )
    if start == end:
        return {
            "start_1based": None,
            "end_1based": None,
            "empty": True,
            "internal_start_0based": start,
            "internal_end_0based": end,
            "internal_convention": INTERNAL_CONVENTION,
            "display_convention": DISPLAY_CONVENTION,
        }
    return {
        "start_1based": start + 1,
        "end_1based": end,
        "empty": False,
        "internal_start_0based": start,
        "internal_end_0based": end,
        "internal_convention": INTERNAL_CONVENTION,
        "display_convention": DISPLAY_CONVENTION,
    }


def display_to_internal(start_1based: int, end_1based: int) -> dict:
    """Converte exibicao 1-based inclusive para [start, end) interno.

    Args:
        start_1based: Primeira base visivel (>= 1).
        end_1based: Ultima base visivel (>= start).

    Returns:
        Dict start_0based, end_0based, convention.

    Raises:
        CoordinateError: INVALID_INPUT.
    """
    try:
        start = int(start_1based)
        end = int(end_1based)
    except (TypeError, ValueError) as exc:
        raise CoordinateError("Display coordinates must be integers.") from exc
    if start < 1 or end < start:
        raise CoordinateError(
            "Display interval must be 1-based inclusive with end >= start >= 1."
        )
    return {
        "start_0based": start - 1,
        "end_0based": end,
        "empty": False,
        "internal_convention": INTERNAL_CONVENTION,
        "display_convention": DISPLAY_CONVENTION,
    }


def roundtrip_ok(start_0based: int, end_0based: int) -> bool:
    """True se internal -> display -> internal reproduz o intervalo.

    Args:
        start_0based: Inicio interno.
        end_0based: Fim interno.

    Returns:
        bool. Intervalos vazios nao tem roundtrip de exibicao.

    Raises:
        Nenhum.
    """
    try:
        shown = internal_to_display(start_0based, end_0based)
    except CoordinateError:
        return False
    if shown.get("empty"):
        return start_0based == end_0based
    back = display_to_internal(int(shown["start_1based"]), int(shown["end_1based"]))
    return back["start_0based"] == int(start_0based) and back["end_0based"] == int(
        end_0based
    )


def pam_interval_on_strand(
    spacer_start_0based: int,
    spacer_end_0based: int,
    *,
    strand: str,
    pam_nt: int = 3,
) -> dict:
    """Intervalo interno do PAM dado o spacer e a fita.

    SpCas9: PAM 3' do spacer na fita do guia. Na fita + o PAM e imediatamente
    apos o spacer. Na fita - o PAM e imediatamente antes do spacer no
    contig sense.

    Args:
        spacer_start_0based: Inicio do spacer no contig sense.
        spacer_end_0based: Fim exclusivo do spacer no contig sense.
        strand: '+' ou '-'.
        pam_nt: Comprimento do PAM (3 para NGG).

    Returns:
        Dict pam_start_0based, pam_end_0based, strand.

    Raises:
        CoordinateError: INVALID_INPUT.
    """
    if strand not in {"+", "-"}:
        raise CoordinateError("Strand must be '+' or '-'.")
    try:
        start = int(spacer_start_0based)
        end = int(spacer_end_0based)
        width = int(pam_nt)
    except (TypeError, ValueError) as exc:
        raise CoordinateError("PAM interval inputs must be integers.") from exc
    if start < 0 or end <= start or width < 1:
        raise CoordinateError("Spacer interval or PAM width is invalid.")
    if strand == "+":
        pam_start = end
        pam_end = end + width
    else:
        pam_start = start - width
        pam_end = start
    if pam_start < 0:
        raise CoordinateError("PAM would start before contig coordinate 0.")
    return {
        "pam_start_0based": pam_start,
        "pam_end_0based": pam_end,
        "strand": strand,
        "internal_convention": INTERNAL_CONVENTION,
    }


def hit_display_fields(hit: Mapping[str, Any]) -> dict:
    """Campos de exibicao 1-based para um hit com coordenadas internas.

    Args:
        hit: Off-target com start_0based e end_0based.

    Returns:
        Dict para tabela. Nao altera o hit.

    Raises:
        CoordinateError: se as coordenadas internas forem invalidas.
    """
    shown = internal_to_display(
        int(hit.get("start_0based") or 0),
        int(hit.get("end_0based") or 0),
    )
    pam_start = hit.get("pam_start_0based")
    pam_end = hit.get("pam_end_0based")
    pam_shown: Optional[dict] = None
    if pam_start is not None and pam_end is not None:
        pam_shown = internal_to_display(int(pam_start), int(pam_end))
    return {
        "contig": str(hit.get("chromosome_or_contig") or ""),
        "strand": str(hit.get("strand") or ""),
        "start_1based": shown.get("start_1based"),
        "end_1based": shown.get("end_1based"),
        "pam_start_1based": None if pam_shown is None else pam_shown.get("start_1based"),
        "pam_end_1based": None if pam_shown is None else pam_shown.get("end_1based"),
        "internal_convention": INTERNAL_CONVENTION,
        "display_convention": DISPLAY_CONVENTION,
    }
