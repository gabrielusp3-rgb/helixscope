"""Pairwise alignment public surface. Smith-Waterman is not BLAST."""

from __future__ import annotations

from typing import Any

from helixscope_core.sequence_input import resolve_pasted_sequence
from modules import alignment
from modules.alignment import (
    MAX_ALIGN_LENGTH,
    MAX_DOTPLOT_LENGTH,
    MATCH_SCORE,
    MISMATCH_SCORE,
    OPEN_GAP_SCORE,
    EXTEND_GAP_SCORE,
    dotplot_matrix,
    pairwise_global,
    pairwise_local,
)

__all__ = (
    "EXTEND_GAP_SCORE",
    "MATCH_SCORE",
    "MAX_ALIGN_LENGTH",
    "MAX_DOTPLOT_LENGTH",
    "MISMATCH_SCORE",
    "OPEN_GAP_SCORE",
    "dotplot_matrix",
    "pairwise_align",
    "pairwise_global",
    "pairwise_local",
    "translate_nucleic_for_alignment",
)


def translate_nucleic_for_alignment(sequence: str) -> dict[str, Any]:
    """Translate one nucleic sequence, frame +1, stop at first stop.

    Args:
        sequence: DNA or RNA (or single-record FASTA).

    Returns:
        translation, source_molecule, frame, stop_semantics, warnings.

    Raises:
        ValueError: Empty translation or invalid alphabet for this workflow.
    """
    from modules.protein_analysis import translate as translate_rna

    residues = str(resolve_pasted_sequence(sequence).get("sequence") or "")
    cleaned = "".join(residues.split()).upper()
    if not cleaned:
        raise ValueError("Nucleic sequence is empty after FASTA parsing.")
    has_u = "U" in cleaned
    has_t = "T" in cleaned
    if has_u and has_t:
        raise ValueError("Mixed T and U is not translated. Provide DNA or RNA.")
    warnings: list[str] = []
    if has_u and not has_t:
        molecule = "RNA"
        protein = translate_rna(cleaned)
    else:
        molecule = "DNA"
        try:
            from Bio.Seq import Seq
        except ImportError as exc:
            raise RuntimeError("Biopython is required for DNA translation.") from exc
        protein = str(Seq(cleaned).translate(table=1, to_stop=True))
        warnings.append("DNA was translated in frame +1 with NCBI standard code 1.")
    if not protein:
        raise ValueError("Translation produced an empty peptide (immediate stop or empty CDS).")
    return {
        "status": "COMPUTED",
        "source_molecule": molecule,
        "frame": 1,
        "stop_semantics": "to_first_stop",
        "genetic_code": "NCBI standard 1",
        "translation": protein,
        "warnings": warnings,
    }


def pairwise_align(
    seq1: str,
    seq2: str,
    *,
    mode: str = "global",
    translate_nucleic: bool = False,
) -> dict[str, Any]:
    """Needleman-Wunsch (global) or Smith-Waterman (local).

    Args:
        seq1: First sequence.
        seq2: Second sequence.
        mode: \"global\" or \"local\".
        translate_nucleic: If True, translate each nucleic input (frame +1,
            first stop) then align as protein. Never silent.

    Returns:
        Alignment dict from modules.alignment. Translation metadata when used.

    Raises:
        ValueError: Invalid mode or sequences.
    """
    chosen = str(mode or "global").strip().lower()
    first = str(resolve_pasted_sequence(seq1).get("sequence") or "")
    second = str(resolve_pasted_sequence(seq2).get("sequence") or "")
    translation_meta: dict[str, Any] | None = None
    if translate_nucleic:
        left = translate_nucleic_for_alignment(first)
        right = translate_nucleic_for_alignment(second)
        first = str(left["translation"])
        second = str(right["translation"])
        translation_meta = {
            "status": "COMPUTED",
            "applied": True,
            "seq1": left,
            "seq2": right,
            "note": "Nucleic sequences were translated before protein alignment. This is not six-frame translation.",
        }
    if chosen == "local":
        payload = pairwise_local(first, second)
    elif chosen == "global":
        payload = pairwise_global(first, second)
    else:
        raise ValueError("Alignment mode must be 'global' or 'local'.")
    aligned_a = str(payload.get("aligned_seq1") or "")
    aligned_b = str(payload.get("aligned_seq2") or "")
    if aligned_a and aligned_b and len(aligned_a) == len(aligned_b):
        payload["column_classes"] = alignment.classify_alignment_columns(aligned_a, aligned_b)
    if translation_meta is not None:
        payload["translation"] = translation_meta
    else:
        payload["translation"] = {"applied": False, "status": "NOT_COMPUTED"}
    return payload
