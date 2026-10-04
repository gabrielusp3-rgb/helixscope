"""Caminho RNA sequence -> representacao 3D. MFE nao vira EXPERIMENTAL.

Camadas: secondary (ViennaRNA, se existir) vs 3D experimental/predicted/
illustrative. Pares MFE sao PREDICTED; geometria da helice A-RNA e
ILLUSTRATIVE. Os rotulos nao se promovem.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from . import dna_3d, dna_analysis, nucleic_geometry, provenance, region_nav, structure_scene

KIND_EXPERIMENTAL = dna_3d.KIND_EXPERIMENTAL
KIND_PREDICTED = dna_3d.KIND_PREDICTED
KIND_ILLUSTRATIVE = dna_3d.KIND_ILLUSTRATIVE
KIND_UNAVAILABLE = dna_3d.KIND_UNAVAILABLE


class Rna3dError(Exception):
    """Falha classificada RNA 3D.

    Attributes:
        category: INVALID_INPUT, UNAVAILABLE, ERROR.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def predicted_availability() -> dict:
    """Nenhum motor 3D de RNA (FARFAR/RNAComposer/SimRNA) esta integrado.

    Args:
        Nenhum.

    Returns:
        UNAVAILABLE.

    Raises:
        Nenhum.
    """
    return {
        "status": "UNAVAILABLE",
        "kind": KIND_UNAVAILABLE,
        "reason": (
            "No RNA 3D prediction engine is integrated. ViennaRNA MFE is "
            "secondary structure only and is not converted into coordinates. "
            "RhoFold+, trRosettaRNA2 and FARFAR2 were researched and are not "
            "shipped (see na_structure_catalog.rna_3d_predictor_research)."
        ),
        "tool": None,
        "mfe_does_not_yield_3d": True,
    }


def representation_for_sequence(
    sequence: str,
    *,
    deposited: Optional[Mapping[str, Any]] = None,
    prefer: Optional[str] = None,
    region_start: int = 0,
    region_end: Optional[int] = None,
    lod: str = region_nav.LOD_HIGH,
    center: Optional[int] = None,
    highlight_indices: Optional[Sequence[int]] = None,
) -> dict:
    """RNA validado -> 3D. Default ILLUSTRATIVE A-RNA se nao houver deposito.

    Args:
        sequence: RNA 5'->3'.
        deposited: Envelope ja validado. Sem fetch.
        prefer: experimental, predicted ou illustrative.
        region_start: Janela.
        region_end: Fim exclusivo.
        lod: LOD.
        center: Centro.
        highlight_indices: Destaques.

    Returns:
        Dict status/kind/envelope analogo a dna_3d.

    Raises:
        Nenhum.
    """
    checked = dna_analysis.validate_for_molecule(sequence, "RNA")
    if not checked["is_valid"]:
        return {
            "status": "ERROR",
            "kind": KIND_UNAVAILABLE,
            "reason": str(checked.get("rejection_reason") or "Sequence is not RNA."),
            "envelope": None,
            "model": None,
            "n_scene_objects": 0,
        }
    seq = str(checked["sequence"])
    digest = provenance.sequence_digest(seq)
    mode = str(prefer or "").strip().lower() or None
    if deposited is None and mode in {None, "experimental"}:
        deposited = dna_3d.bundled_experimental_for_sequence(seq, "RNA")
    if mode == "predicted":
        payload = predicted_availability()
        payload["sequence_hash"] = digest
        payload["envelope"] = None
        payload["model"] = None
        payload["n_scene_objects"] = 0
        if isinstance(deposited, Mapping) and str(deposited.get("kind") or "") == KIND_PREDICTED:
            if str(deposited.get("sequence_hash") or "") == digest:
                payload = {
                    "status": "PREDICTED",
                    "kind": KIND_PREDICTED,
                    "sequence_hash": digest,
                    "envelope": deposited if deposited.get("mappings") is not None else deposited,
                    "model": None,
                    "n_scene_objects": len(list(deposited.get("atoms") or [])),
                }
        return payload
    if isinstance(deposited, Mapping) and str(deposited.get("sequence_hash") or "") == digest:
        dep_kind = str(deposited.get("kind") or "").strip().lower()
        if dep_kind == KIND_EXPERIMENTAL and mode in {None, "experimental"}:
            return {
                "status": "EXPERIMENTAL",
                "kind": KIND_EXPERIMENTAL,
                "sequence_hash": digest,
                "envelope": deposited,
                "model": None,
                "partial_view": bool(deposited.get("partial_view")),
                "n_scene_objects": len(list(deposited.get("atoms") or [])),
            }
    if isinstance(deposited, Mapping) and str(deposited.get("sequence_hash") or "") not in {"", digest}:
        if mode != "illustrative":
            return {
                "status": "UNMAPPED",
                "kind": KIND_UNAVAILABLE,
                "reason": (
                    "Deposited RNA structure sequence hash does not match. "
                    "ILLUSTRATIVE is not substituted automatically."
                ),
                "sequence_hash": digest,
                "envelope": None,
                "model": None,
                "n_scene_objects": 0,
            }
    if mode == "experimental":
        return {
            "status": "UNAVAILABLE",
            "kind": KIND_UNAVAILABLE,
            "reason": "No deposited experimental RNA structure is attached.",
            "sequence_hash": digest,
            "envelope": None,
            "model": None,
            "n_scene_objects": 0,
        }
    try:
        model = nucleic_geometry.illustrative_arna_model(
            seq,
            highlight_indices=list(highlight_indices or []),
            region_start=region_start,
            region_end=region_end,
            lod=lod,
            center=center,
        )
        envelope = nucleic_geometry.structure_3d_input_from_illustrative(model)
    except nucleic_geometry.NucleicGeometryError as exc:
        return {
            "status": "ERROR",
            "kind": KIND_UNAVAILABLE,
            "reason": str(exc),
            "sequence_hash": digest,
            "envelope": None,
            "model": None,
            "n_scene_objects": 0,
        }
    return {
        "status": "ILLUSTRATIVE",
        "kind": KIND_ILLUSTRATIVE,
        "sequence_hash": digest,
        "envelope": envelope,
        "model": model,
        "partial_view": bool(model.get("partial_view")),
        "n_strands": 2,
        "n_pairs": int(model.get("n_pairs") or 0),
        "fragment": model.get("fragment"),
        "n_scene_objects": int(model.get("n_atoms") or 0) + int(model.get("n_pairs") or 0),
        "disclaimer": model.get("disclaimer"),
        "mfe_does_not_yield_3d": True,
    }


def pairing_sync(
    representation: Mapping[str, Any],
    mfe_result: Optional[Mapping[str, Any]],
) -> dict:
    """Pares MFE (PREDICTED) sobre geometria ILLUSTRATIVE. Sem promocao.

    Args:
        representation: Saida de representation_for_sequence.
        mfe_result: Envelope ViennaRNA ou None.

    Returns:
        Dict geometry_kind, pairing_status, pairs. Se MFE ausente, pairing
        UNAVAILABLE; pares da helice permanecem ILLUSTRATIVE (Watson-Crick
        canonico, nao MFE).

    Raises:
        Nenhum.
    """
    geometry_kind = str(representation.get("kind") or KIND_UNAVAILABLE)
    helix_pairs = list((representation.get("model") or {}).get("pairs") or [])
    if not isinstance(mfe_result, Mapping):
        return {
            "geometry_kind": geometry_kind,
            "pairing_status": "UNAVAILABLE",
            "pairing_source": None,
            "mfe_pairs": [],
            "illustrative_helix_pairs": helix_pairs,
            "synced": False,
            "note": (
                "No secondary-structure computation is attached. Canonical "
                "helix pairs are ILLUSTRATIVE Watson-Crick complements, not MFE."
            ),
        }
    mfe_status = str(mfe_result.get("status") or "")
    mfe_hash = str(mfe_result.get("sequence_hash") or "")
    geom_hash = str(representation.get("sequence_hash") or "")
    if mfe_hash and geom_hash and mfe_hash != geom_hash:
        return {
            "geometry_kind": geometry_kind,
            "pairing_status": "STALE",
            "pairing_source": None,
            "mfe_pairs": [],
            "illustrative_helix_pairs": helix_pairs,
            "synced": False,
            "note": "stale/incompatible secondary structure",
        }
    raw_pairs = list(mfe_result.get("base_pairs") or mfe_result.get("pairs") or [])
    parsed: list = []
    for pair in raw_pairs:
        if isinstance(pair, Mapping):
            left = pair.get("position_0based", pair.get("i"))
            right = pair.get("paired_position_0based", pair.get("j"))
        elif isinstance(pair, (list, tuple)) and len(pair) >= 2:
            left, right = pair[0], pair[1]
        else:
            continue
        if left is None or right is None:
            continue
        parsed.append((int(left), int(right)))
    if mfe_status not in {"PREDICTED", "COMPUTED"} or not parsed:
        return {
            "geometry_kind": geometry_kind,
            "pairing_status": "UNAVAILABLE",
            "pairing_source": mfe_result.get("tool"),
            "mfe_pairs": [],
            "illustrative_helix_pairs": helix_pairs,
            "synced": False,
            "note": "MFE pairs are not available. Helix pairs stay ILLUSTRATIVE.",
        }
    view_start = int(
        representation.get("view_start")
        or (representation.get("model") or {}).get("view_start")
        or 0
    )
    view_end = representation.get("view_end") or (representation.get("model") or {}).get("view_end")
    window_pairs = []
    if view_end is None:
        view_end = view_start + len(helix_pairs)
    for left, right in parsed:
        if view_start <= left < int(view_end) and view_start <= right < int(view_end):
            window_pairs.append({"i": left, "j": right, "source": "ViennaRNA MFE"})
    return {
        "geometry_kind": geometry_kind,
        "pairing_status": "PREDICTED",
        "pairing_source": mfe_result.get("tool") or "ViennaRNA",
        "mfe_pairs": window_pairs,
        "illustrative_helix_pairs": helix_pairs,
        "synced": True,
        "note": (
            "MFE pair indices are PREDICTED secondary structure. They are not "
            "EXPERIMENTAL coordinates and do not relabel the ILLUSTRATIVE helix."
        ),
    }


def pairing_geometry_consistency(
    representation: Mapping[str, Any],
    pairs: Sequence[Mapping[str, Any]],
    *,
    max_distance_angstrom: float = 25.0,
) -> dict:
    """Compara pares secundarios com distancias 3D depositadas.

    Args:
        representation: Representacao EXPERIMENTAL com envelope mapeado.
        pairs: Lista com position_0based e paired_position_0based.
        max_distance_angstrom: Teto de sanidade C1'/P-P (nao e um corte de
            ligacao de hidrogenio publicado).

    Returns:
        Dict n_checked, n_within, n_beyond, distances. Sem inventar pares.

    Raises:
        Rna3dError: INVALID_INPUT se a geometria nao for experimental.

    Nota biologica:
        Distancia C1'-C1' (ou P-P) entre indices de um par ViennaRNA/duplex.
        25 A e um teto de consistencia grosseiro, nao o criterio de Leontis-Westhof.
    """
    if str(representation.get("kind") or "") != KIND_EXPERIMENTAL:
        raise Rna3dError("Pairing geometry check requires an EXPERIMENTAL RNA 3D representation.")
    envelope = representation.get("envelope")
    if not isinstance(envelope, Mapping):
        raise Rna3dError("Experimental representation has no envelope.")
    rows = list(envelope.get("mappings") or envelope.get("residue_mapping") or [])
    by_query = {
        int(item["query_index_0based"]): item
        for item in rows
        if item.get("query_index_0based") is not None
    }
    distances: list = []
    n_within = 0
    n_beyond = 0
    n_unmapped = 0
    for pair in pairs:
        if not isinstance(pair, Mapping):
            continue
        left = pair.get("position_0based", pair.get("i"))
        right = pair.get("paired_position_0based", pair.get("j"))
        if left is None or right is None:
            continue
        a = by_query.get(int(left))
        b = by_query.get(int(right))
        if a is None or b is None:
            n_unmapped += 1
            distances.append({"i": int(left), "j": int(right), "distance_angstrom": None, "status": "UNMAPPED"})
            continue
        pa = structure_scene._backbone_xyz(a)
        pb = structure_scene._backbone_xyz(b)
        if pa is None or pb is None:
            n_unmapped += 1
            distances.append({"i": int(left), "j": int(right), "distance_angstrom": None, "status": "UNMAPPED"})
            continue
        dx = float(pa["x"]) - float(pb["x"])
        dy = float(pa["y"]) - float(pb["y"])
        dz = float(pa["z"]) - float(pb["z"])
        dist = (dx * dx + dy * dy + dz * dz) ** 0.5
        within = dist <= float(max_distance_angstrom)
        if within:
            n_within += 1
        else:
            n_beyond += 1
        distances.append(
            {
                "i": int(left),
                "j": int(right),
                "distance_angstrom": dist,
                "status": "WITHIN" if within else "BEYOND",
            }
        )
    return {
        "status": "COMPUTED",
        "geometry_kind": KIND_EXPERIMENTAL,
        "max_distance_angstrom": float(max_distance_angstrom),
        "n_checked": n_within + n_beyond,
        "n_within": n_within,
        "n_beyond": n_beyond,
        "n_unmapped": n_unmapped,
        "distances": distances,
        "note": (
            "Euclidean C1'/P distance between secondary-structure pair indices "
            "on deposited coordinates. Not a hydrogen-bond assignment."
        ),
    }
