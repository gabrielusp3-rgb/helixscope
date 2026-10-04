"""Mapeamento genoma -> CDS -> codon -> residuo. Sem predicao de efeito.

So corre quando a anotacao e o transcrito sao validos. Mapeamentos ambiguos
(varios transcritos) sao todos devolvidos. Nao escolhe transcrito principal.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from . import provenance

CODON_NT: int = 3


class CdsMappingError(ValueError):
    """Mapeamento recusado.

    Attributes:
        category: INVALID_INPUT ou UNMAPPED.
    """

    def __init__(self, message: str, category: str = "UNMAPPED") -> None:
        super().__init__(message)
        self.category = str(category or "UNMAPPED").strip().upper() or "UNMAPPED"


def _cds_blocks_for_transcript(
    features: Sequence[Mapping[str, Any]],
    transcript_id: str,
) -> list[dict]:
    blocks = [
        dict(item)
        for item in features
        if str(item.get("type") or "") == "cds"
        and str(item.get("transcript_id") or "") == transcript_id
    ]
    if not blocks:
        return []
    strands = {str(item.get("strand") or "") for item in blocks}
    if len(strands) != 1 or "." in strands:
        return []
    strand = next(iter(strands))
    if strand == "+":
        blocks.sort(key=lambda item: int(item["start_0based"]))
    else:
        blocks.sort(key=lambda item: int(item["start_0based"]), reverse=True)
    return blocks


def map_genomic_base_to_codon(
    features: Sequence[Mapping[str, Any]],
    *,
    contig: str,
    position_0based: int,
) -> dict:
    """Mapeia uma base genomica a codon/residuo em cada transcrito CDS.

    Args:
        features: Features GFF do contig (ou de toda a anotacao).
        contig: seqid.
        position_0based: base 0-based no contig sense.

    Returns:
        Dict mappings (lista), unmapped se vazia. effect permanece None.

    Raises:
        CdsMappingError: INVALID_INPUT.
    """
    ident = str(contig or "").strip()
    if not ident:
        raise CdsMappingError("Contig is required.", "INVALID_INPUT")
    try:
        pos = int(position_0based)
    except (TypeError, ValueError) as exc:
        raise CdsMappingError("position_0based must be an integer.", "INVALID_INPUT") from exc
    if pos < 0:
        raise CdsMappingError("position_0based must be >= 0.", "INVALID_INPUT")
    relevant = [
        item
        for item in features
        if str(item.get("seqid") or "") == ident and str(item.get("type") or "") == "cds"
    ]
    transcript_ids: list[str] = []
    for item in relevant:
        tx = str(item.get("transcript_id") or "")
        if tx and tx not in transcript_ids:
            transcript_ids.append(tx)
    mappings: list[dict] = []
    for tx in transcript_ids:
        blocks = _cds_blocks_for_transcript(relevant, tx)
        if not blocks:
            continue
        strand = str(blocks[0].get("strand") or "")
        offset = 0
        hit = None
        for block in blocks:
            start = int(block["start_0based"])
            end = int(block["end_0based"])
            length = end - start
            if start <= pos < end:
                if strand == "+":
                    local = pos - start
                else:
                    local = (end - 1) - pos
                nt_index = offset + local
                codon_index = nt_index // CODON_NT
                position_in_codon = nt_index % CODON_NT
                residue_number = codon_index + 1
                hit = {
                    "transcript_id": tx,
                    "gene_id": str(block.get("gene_id") or ""),
                    "gene_name": str(block.get("gene_name") or ""),
                    "strand": strand,
                    "cds_nt_index_0based": nt_index,
                    "codon_index_0based": codon_index,
                    "position_in_codon_0based": position_in_codon,
                    "protein_residue_number_1based": residue_number,
                    "contig": ident,
                    "genomic_position_0based": pos,
                }
                break
            offset += length
        if hit is not None:
            mappings.append(hit)
    return {
        "contig": ident,
        "genomic_position_0based": pos,
        "n_mappings": len(mappings),
        "mappings": mappings,
        "unmapped": not mappings,
        "transcript_policy": "all_overlapping_cds_transcripts",
        "effect": None,
        "effect_note": (
            "Codon/residue identity only when CDS coordinates are valid. "
            "No missense, nonsense or splicing prediction."
        ),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def map_hit_to_protein_residue(
    hit: Mapping[str, Any],
    index: Mapping[str, Any],
) -> dict:
    """Atalho a partir de um hit CRISPR (primeira base do spacer).

    Args:
        hit: Hit com chromosome_or_contig e start_0based.
        index: Indice de anotacao.

    Returns:
        Saida de map_genomic_base_to_codon.

    Raises:
        CdsMappingError.
    """
    contig = str(hit.get("chromosome_or_contig") or "")
    features: list[dict] = []
    by_contig = index.get("by_contig") or {}
    features.extend(list(by_contig.get(contig) or []))
    return map_genomic_base_to_codon(
        features,
        contig=contig,
        position_0based=int(hit.get("start_0based") or 0),
    )
