"""Identidade de uma variante genomica. Sem predicao de efeito.

Objeto comum para fases posteriores: assembly, contig, posicao, ref, alt.
Nao classifica patogenia, nao escolhe transcrito, nao mapeia para proteina.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping

from . import genome_coordinates, provenance

DNA_LETTERS: frozenset[str] = frozenset("ACGT")


class VariantError(ValueError):
    """Variante mal formada.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT").strip().upper() or "INVALID_INPUT"


def normalize_allele(text: str) -> str:
    """DNA A/C/G/T em maiusculas. U vira T. Vazio e insercao/delecao explicita.

    Args:
        text: Alelo.

    Returns:
        String normalizada (pode ser vazia para indel simbolico).

    Raises:
        VariantError: bases nao DNA.
    """
    seq = str(text or "").strip().upper().replace("U", "T")
    if seq in {".", "-", "N/A"}:
        return ""
    if set(seq) - DNA_LETTERS:
        raise VariantError(
            "Ref/alt must be A/C/G/T (empty string allowed for symbolic indel)."
        )
    return seq


def genomic_variant(
    *,
    assembly: str,
    accession: str,
    contig: str,
    position_0based: int,
    ref: str,
    alt: str,
    assembly_sha256: str = "",
) -> dict:
    """Constroi um objeto de identidade. Nao prediz efeito.

    Args:
        assembly: id de catalogo (ex. GRCh38.p14) ou vazio se desconhecido.
        accession: accession da assembly quando conhecido.
        contig: identificador do contig/cromossoma no FASTA (sem sinonimo chr).
        position_0based: primeira base afetada, 0-based.
        ref: alelo de referencia no contig.
        alt: alelo alternativo.
        assembly_sha256: hash do FASTA READY, se houver.

    Returns:
        Dict identity-only. effect = None sempre.

    Raises:
        VariantError: INVALID_INPUT.
    """
    ident = str(contig or "").strip()
    if not ident or any(ch in ident for ch in "/\\:*?\"<>|") or ident in {".", ".."}:
        raise VariantError("Contig identifier is missing or not a safe token.")
    try:
        pos = int(position_0based)
    except (TypeError, ValueError) as exc:
        raise VariantError("position_0based must be an integer.") from exc
    if pos < 0:
        raise VariantError("position_0based must be >= 0.")
    ref_n = normalize_allele(ref)
    alt_n = normalize_allele(alt)
    if not ref_n and not alt_n:
        raise VariantError("Ref and alt cannot both be empty.")
    display = genome_coordinates.internal_to_display(pos, pos + max(len(ref_n), 1))
    return {
        "kind": "genomic_variant_identity",
        "assembly": str(assembly or ""),
        "accession": str(accession or ""),
        "contig": ident,
        "position_0based": pos,
        "position_1based": display.get("start_1based"),
        "ref": ref_n,
        "alt": alt_n,
        "assembly_sha256": str(assembly_sha256 or ""),
        "effect": None,
        "effect_status": "NOT_COMPUTED",
        "effect_note": (
            "HelixScope records variant identity only in this phase. "
            "No pathogenic, deleterious or functional claim is attached."
        ),
        "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
        "display_convention": genome_coordinates.DISPLAY_CONVENTION,
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def variant_from_mapping(payload: Mapping[str, Any]) -> dict:
    """Valida um dict ja parcialmente preenchido.

    Args:
        payload: Campos assembly, contig, position_0based, ref, alt.

    Returns:
        Objeto normalizado.

    Raises:
        VariantError.
    """
    return genomic_variant(
        assembly=str(payload.get("assembly") or ""),
        accession=str(payload.get("accession") or ""),
        contig=str(payload.get("contig") or payload.get("chromosome_or_contig") or ""),
        position_0based=int(payload.get("position_0based")),
        ref=str(payload.get("ref") or ""),
        alt=str(payload.get("alt") or ""),
        assembly_sha256=str(payload.get("assembly_sha256") or ""),
    )
