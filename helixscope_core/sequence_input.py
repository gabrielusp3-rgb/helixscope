"""Authoritative pasted-sequence contract for one-molecule modules.

Matches Streamlit ``_resolve_sequence`` for the text path: parse FASTA first,
then let each molecule validator clean whitespace. Does not silently strip
digits or punctuation from non-FASTA garbage.
"""

from __future__ import annotations

from typing import Any

from modules.dna_analysis import parse_sequence_payload

__all__ = ("resolve_pasted_sequence",)


def resolve_pasted_sequence(text: str) -> dict[str, Any]:
    """Parse a single pasted payload the same way Streamlit does for text input.

    FASTA with one record returns only residues. FASTA with several records is
    refused (the first record is never chosen). Non-FASTA text is returned as
    raw payload; molecule validators still remove ordinary whitespace later.

    Args:
        text: Browser textarea or equivalent pasted string.

    Returns:
        Dict with ``sequence``, ``record_count``, ``format`` (``empty``,
        ``fasta``, or ``raw``) and ``identifier``.

    Raises:
        TypeError: If ``text`` is not a string.
        ValueError: Multi-FASTA, empty FASTA parse, or size rejection from
            ``parse_sequence_payload``.

    Nota biologica:
        Concatenar registros FASTA distintos inventaria uma molecula que nao
        existe. Um cabecalho FASTA nao e sequencia.
    """
    if not isinstance(text, str):
        raise TypeError("A sequencia deve ser uma string.")
    return parse_sequence_payload(text)
