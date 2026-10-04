"""Vista regional local (nao o genoma inteiro no DOM).

Mostra eixo de coordenadas, hit/PAM e features que intersectam a janela.
Zoom e deslocamento sao recortes do contig, nunca um browser de cromossoma.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from . import genome_coordinates, provenance

DEFAULT_FLANK_NT: int = 80
MAX_WINDOW_NT: int = 4000


class RegionViewError(ValueError):
    """Janela regional invalida.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT").strip().upper() or "INVALID_INPUT"


def clip_window(
    contig_length: int,
    center_start_0based: int,
    center_end_0based: int,
    *,
    flank_nt: int = DEFAULT_FLANK_NT,
) -> dict:
    """Janela [start, end) no contig, limitada ao comprimento.

    Args:
        contig_length: Bases do contig.
        center_start_0based: Inicio do hit.
        center_end_0based: Fim exclusivo do hit.
        flank_nt: Flancos pedidos.

    Returns:
        Dict start_0based, end_0based, clipped.

    Raises:
        RegionViewError.
    """
    try:
        length = int(contig_length)
        start = int(center_start_0based)
        end = int(center_end_0based)
        flank = int(flank_nt)
    except (TypeError, ValueError) as exc:
        raise RegionViewError("Window bounds must be integers.") from exc
    if length <= 0 or start < 0 or end < start or end > length:
        raise RegionViewError("Hit interval is outside the contig.")
    if flank < 0:
        raise RegionViewError("flank_nt must be >= 0.")
    win_start = max(0, start - flank)
    win_end = min(length, end + flank)
    if win_end - win_start > MAX_WINDOW_NT:
        raise RegionViewError(
            f"Region window exceeds {MAX_WINDOW_NT} nt. HelixScope does not "
            "render a full-genome track."
        )
    return {
        "start_0based": win_start,
        "end_0based": win_end,
        "length": win_end - win_start,
        "contig_length": length,
        "clipped": win_start > 0 or win_end < length,
        "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
    }


def local_region_view(
    *,
    contig: str,
    contig_length: int,
    hit: Mapping[str, Any],
    features: Sequence[Mapping[str, Any]] = (),
    flank_nt: int = DEFAULT_FLANK_NT,
    sequence: str = "",
) -> dict:
    """Monta a vista local. Nao e um browser de genoma.

    Args:
        contig: ID do contig.
        contig_length: Comprimento.
        hit: Off-target com start/end/PAM/strand.
        features: Features GFF ja recortadas ou brutas (serao filtradas).
        flank_nt: Flancos.
        sequence: Subsequencia da janela no sentido + se ja extraida.

    Returns:
        Dict axis, hit, pam, features, sequence.

    Raises:
        RegionViewError.
    """
    ident = str(contig or "").strip()
    if not ident:
        raise RegionViewError("Contig is required.")
    hit_start = int(hit.get("start_0based") or 0)
    hit_end = int(hit.get("end_0based") or 0)
    window = clip_window(contig_length, hit_start, hit_end, flank_nt=flank_nt)
    display_axis = genome_coordinates.internal_to_display(
        window["start_0based"], window["end_0based"]
    )
    pam = None
    if hit.get("pam_start_0based") is not None and hit.get("pam_end_0based") is not None:
        pam = {
            "start_0based": int(hit["pam_start_0based"]),
            "end_0based": int(hit["pam_end_0based"]),
            "sequence": str(hit.get("PAM") or hit.get("pam_sequence") or ""),
            "strand": str(hit.get("strand") or ""),
            "display": genome_coordinates.internal_to_display(
                int(hit["pam_start_0based"]), int(hit["pam_end_0based"])
            ),
        }
    visible: list[dict] = []
    for feat in features:
        f_start = int(feat.get("start_0based") or 0)
        f_end = int(feat.get("end_0based") or 0)
        if f_end <= window["start_0based"] or f_start >= window["end_0based"]:
            continue
        visible.append(
            {
                "seqid": str(feat.get("seqid") or ident),
                "type": str(feat.get("type") or ""),
                "start_0based": f_start,
                "end_0based": f_end,
                "strand": str(feat.get("strand") or "."),
                "gene_id": str(feat.get("gene_id") or ""),
                "gene_name": str(feat.get("gene_name") or ""),
                "display": genome_coordinates.internal_to_display(f_start, f_end),
            }
        )
    seq = str(sequence or "")
    expected = int(window["length"])
    if seq and len(seq) != expected:
        raise RegionViewError(
            "Provided window sequence length does not match the region window."
        )
    return {
        "kind": "local_region_view",
        "full_genome_dom": False,
        "contig": ident,
        "axis": {
            **window,
            "display": display_axis,
        },
        "hit": {
            "start_0based": hit_start,
            "end_0based": hit_end,
            "strand": str(hit.get("strand") or ""),
            "sequence": str(hit.get("target_sequence") or hit.get("sequence") or ""),
            "display": genome_coordinates.internal_to_display(hit_start, hit_end),
        },
        "pam": pam,
        "features": visible,
        "sequence_plus_strand": seq,
        "note": (
            "Local window only. Not a full-chromosome or full-genome track. "
            "Feature overlap is coordinate identity, not a functional effect."
        ),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def shift_window(
    view: Mapping[str, Any],
    *,
    delta_nt: int,
    contig_length: int,
) -> dict:
    """Desloca a janela sem sair do contig.

    Args:
        view: Vista anterior.
        delta_nt: Deslocamento (negativo = upstream no contig +).
        contig_length: Comprimento.

    Returns:
        Novos bounds 0-based.

    Raises:
        RegionViewError.
    """
    axis = dict(view.get("axis") or {})
    try:
        start = int(axis.get("start_0based"))
        end = int(axis.get("end_0based"))
        delta = int(delta_nt)
        length = int(contig_length)
    except (TypeError, ValueError) as exc:
        raise RegionViewError("Shift inputs must be integers.") from exc
    width = end - start
    if width <= 0:
        raise RegionViewError("Current window is empty.")
    start = min(max(0, start + delta), max(0, length - width))
    end = min(length, start + width)
    return {
        "start_0based": start,
        "end_0based": end,
        "length": end - start,
        "contig_length": length,
        "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
    }
