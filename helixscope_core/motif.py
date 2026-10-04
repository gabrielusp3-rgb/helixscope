"""Motif search public surface. Hits are pattern matches, not function."""

from __future__ import annotations

from typing import Any, List

from helixscope_core.sequence_input import resolve_pasted_sequence
from modules.motif_search import find_motif, find_motif_protein, find_motif_rna, find_motifs

__all__ = (
    "analyze_motif",
    "find_motif",
    "find_motif_protein",
    "find_motif_rna",
    "find_motifs",
    "search_motifs",
)


def search_motifs(sequence: str, pattern: str, *, molecule: str = "DNA") -> List[dict[str, Any]]:
    """IUPAC/exact search on one pattern.

    Args:
        sequence: Target sequence.
        pattern: Motif including IUPAC codes.
        molecule: DNA, RNA, or PROTEIN.

    Returns:
        Hit list with 0-based half-open coordinates.
    """
    kind = str(molecule or "DNA").strip().upper()
    target = str(resolve_pasted_sequence(sequence).get("sequence") or "")
    if kind == "RNA":
        return find_motif_rna(target, pattern)
    if kind == "PROTEIN":
        return find_motif_protein(target, pattern)
    return find_motif(target, pattern)


def analyze_motif(sequence: str, pattern: str, *, molecule: str = "DNA") -> dict[str, Any]:
    """Complete motif search result. Hits are pattern matches, not function.

    Args:
        sequence: Target sequence or single-record FASTA.
        pattern: Motif including IUPAC codes.
        molecule: DNA, RNA, or PROTEIN.

    Returns:
        Domain dict with hits, coordinates, n_hits, hit_meaning, feature track.

    Raises:
        TypeError / ValueError: From parsers or motif search.
    """
    parsed = resolve_pasted_sequence(sequence)
    hits = search_motifs(sequence, pattern, molecule=molecule)
    features = []
    for hit in hits:
        if not isinstance(hit, dict):
            continue
        features.append(
            {
                "start": hit.get("start"),
                "end": hit.get("end"),
                "label": str(hit.get("match") or pattern),
                "kind": "motif",
                "strand": hit.get("strand"),
            }
        )
    return {
        "status": "COMPUTED",
        "hits_status": "NO_HITS" if not hits else "COMPUTED",
        "pattern": pattern,
        "molecule": str(molecule or "DNA").strip().upper(),
        "sequence": str(parsed.get("sequence") or ""),
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
        "hits": hits,
        "n_hits": len(hits),
        "hit_meaning": "pattern occurrence, not proven biological function",
        "sequence_feature_tracks": (
            [{"name": "Motif hits", "features": features}] if features else []
        ),
    }
