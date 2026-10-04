"""MSA facade. Pairwise alignment is not MSA."""

from __future__ import annotations

import re
from io import StringIO
from typing import Any

from modules.msa import (
    MAX_SEQUENCES,
    MsaError,
    build_collection_member,
    build_msa_result,
    import_prealigned_fasta,
    parse_aligned_fasta,
    parse_unaligned_fasta,
    poll_msa,
    retrieve_msa,
    submit_msa,
    tool_availability,
    validate_collection,
)

SUPPORTED_ALIGNMENT_FORMATS: tuple[str, ...] = (
    "auto",
    "fasta",
    "clustal",
    "stockholm",
    "phylip",
)

_ANNOTATED_TABLE_EXAMPLE = (
    ">seq1\nATGCATGCATGC\n>seq2\nATGTATGCATGC\n"
)

__all__ = (
    "MAX_SEQUENCES",
    "SUPPORTED_ALIGNMENT_FORMATS",
    "build_collection_member",
    "build_msa_result",
    "classify_alignment_text",
    "import_alignment",
    "import_prealigned_fasta",
    "parse_aligned_fasta",
    "parse_unaligned_fasta",
    "poll_msa",
    "retrieve_msa",
    "submit_msa",
    "tool_availability",
    "validate_collection",
)


def classify_alignment_text(text: str) -> str:
    """Classify pasted text as a supported alignment format or reject tables.

    Args:
        text: User paste.

    Returns:
        One of fasta, clustal, stockholm, phylip, annotated_table, empty, unknown.

    Raises:
        None.
    """
    raw = str(text or "")
    stripped = raw.strip()
    if not stripped:
        return "empty"
    head = stripped.lstrip()
    first = head.splitlines()[0].strip() if head.splitlines() else ""
    if _looks_like_annotated_table(stripped):
        return "annotated_table"
    if first.upper().startswith("CLUSTAL"):
        return "clustal"
    if first.startswith("# STOCKHOLM") or first.startswith("#=GF"):
        return "stockholm"
    if re.match(r"^\s*\d+\s+\d+\s*$", first):
        return "phylip"
    if first.startswith(">"):
        return "fasta"
    return "unknown"


def _bracket_holds_bases(text: str) -> bool:
    """True when one bracket pair holds at least three base characters.

    Args:
        text: Candidate alignment body.

    Returns:
        True only for a closed bracket whose interior has three base
        characters. The scan is linear. A long unclosed bracket is not a match.

    Raises:
        None.
    """
    bases = frozenset("ACGTUacgtu ")
    body = text[:200_000]
    start = 0
    while start < len(body):
        open_at = body.find("[", start)
        if open_at < 0:
            return False
        close_at = body.find("]", open_at + 1)
        if close_at < 0:
            return False
        window = body[open_at + 1:close_at]
        if len(window) <= 2_000:
            count = 0
            for char in window:
                if char in bases:
                    count += 1
                    if count >= 3:
                        return True
        start = open_at + 1
    return False


def _looks_like_annotated_table(text: str) -> bool:
    """True for markdown/human tables that are not FASTA/Clustal/Stockholm."""
    if text.lstrip().startswith(">"):
        return False
    first = text.lstrip().splitlines()[0].strip() if text.strip() else ""
    if first.upper().startswith("CLUSTAL") or first.startswith("# STOCKHOLM"):
        return False
    if "```" in text:
        return True
    if "|" in text and re.search(r"\|[ \t]*-{2,}", text):
        return True
    if any(mark in text for mark in ("←", "→", "->", "<-")):
        return True
    if _bracket_holds_bases(text):
        return True
    if re.search(r"(?im)^(humano|chimpanz[eé]|gorilla|sequence\s+\d+)\b", text):
        return True
    numbered = bool(re.search(r"(?m)^\s*\d{1,4}(\s+\d{1,4}){3,}\s*$", text))
    if numbered and not first.startswith(">"):
        return True
    return False


def _alignio_to_aligned_fasta(text: str, fmt: str) -> str:
    """Convert a Biopython AlignIO format into aligned FASTA.

    Args:
        text: Alignment body.
        fmt: clustal, stockholm, or phylip.

    Returns:
        FASTA with gap columns preserved.

    Raises:
        MsaError: INVALID_FORMAT or PARSING_ERROR.
    """
    try:
        from Bio import AlignIO
    except ImportError as exc:
        raise MsaError("Biopython AlignIO is required to parse this format.", "PARSING_ERROR") from exc
    try:
        alignment = AlignIO.read(StringIO(text), fmt)
    except Exception as exc:
        raise MsaError(
            "This input is not a recognized alignment format. "
            "Use FASTA, aligned FASTA, or another supported alignment format "
            f"(Clustal, Stockholm, PHYLIP). Parser: {type(exc).__name__}.",
            "INVALID_FORMAT",
        ) from exc
    if alignment is None or len(alignment) < 2:
        raise MsaError(
            "This input is not a recognized alignment format. "
            "Use FASTA, aligned FASTA, or another supported alignment format.",
            "INVALID_FORMAT",
        )
    lines: list[str] = []
    for record in alignment:
        ident = str(record.id or "").strip() or "seq"
        lines.append(f">{ident}")
        lines.append(str(record.seq or ""))
    return "\n".join(lines) + "\n"


def import_alignment(text: str, *, fmt: str = "auto") -> dict[str, Any]:
    """Import a user-declared alignment. Does not run Clustal/MAFFT/MUSCLE.

    Args:
        text: Alignment body.
        fmt: auto, fasta, clustal, stockholm, or phylip.

    Returns:
        MSA COMPLETED envelope from import_prealigned_fasta plus input_format.

    Raises:
        MsaError: INVALID_FORMAT, PARSING_ERROR, INVALID_INPUT, RESOURCE_LIMIT.
    """
    chosen = str(fmt or "auto").strip().lower()
    if chosen not in SUPPORTED_ALIGNMENT_FORMATS:
        raise MsaError(
            "Input format must be auto, fasta, clustal, stockholm, or phylip.",
            "INVALID_INPUT",
        )
    stripped = str(text or "").strip()
    if not stripped:
        raise MsaError(
            "This input is not a recognized alignment format. "
            "Use FASTA, aligned FASTA, or another supported alignment format.",
            "INVALID_FORMAT",
        )
    detected = classify_alignment_text(stripped)
    if detected == "annotated_table" and chosen in {"auto", "fasta"}:
        raise MsaError(
            "This input is not a recognized alignment format. "
            "Use FASTA, aligned FASTA, or another supported alignment format. "
            f"Example:\n{_ANNOTATED_TABLE_EXAMPLE}",
            "INVALID_FORMAT",
        )
    if chosen == "auto":
        if detected == "unknown":
            raise MsaError(
                "This input is not a recognized alignment format. "
                "Use FASTA, aligned FASTA, or another supported alignment format. "
                f"Example:\n{_ANNOTATED_TABLE_EXAMPLE}",
                "INVALID_FORMAT",
            )
        chosen = detected if detected != "empty" else "fasta"
    if chosen == "fasta":
        try:
            payload = import_prealigned_fasta(stripped)
        except MsaError as exc:
            if exc.category == "PARSING_ERROR" and not stripped.lstrip().startswith(">"):
                raise MsaError(
                    "This input is not a recognized alignment format. "
                    "Use FASTA, aligned FASTA, or another supported alignment format. "
                    f"Example:\n{_ANNOTATED_TABLE_EXAMPLE}",
                    "INVALID_FORMAT",
                ) from exc
            raise
    else:
        fasta = _alignio_to_aligned_fasta(stripped, chosen)
        payload = import_prealigned_fasta(fasta)
    payload["input_format"] = chosen
    payload["format_declared"] = str(fmt or "auto")
    return payload

