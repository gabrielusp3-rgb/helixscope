"""Scene model for macromolecular 3D visualization.

Transforms a validated structure_3d_input envelope into visual objects.
Does not fetch data, parse mmCIF/PDB, invent coordinates, create mapping,
or mutate the scientific envelope.

Backbone traces connect deposited CA (protein) or P/C1' (nucleic acid)
coordinates in polymer order. This is not a DSSP/STRIDE cartoon, not a
molecular surface, and not a Mol*/PyMOL renderer.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from . import protein_structure, provenance, region_nav, scientific_checks

ENGINE: str = "plotly"
REPRESENTATIONS: Tuple[str, ...] = ("backbone", "ca", "atoms")
COLOR_MODES: Tuple[str, ...] = (
    "chain",
    "residue_type",
    "selection",
    "plddt",
    "hydropathy",
    "conservation",
    "charge",
)
STRUCTURE_EVENT_KINDS: Tuple[str, ...] = (
    "residueSelected",
    "sequencePositionSelected",
    "chainSelected",
    "modelSelected",
    "structureFeatureSelected",
)

CHAIN_PALETTE: Tuple[str, ...] = (
    "#6893D0",
    "#34D399",
    "#F472B6",
    "#FBBF24",
    "#91B4E4",
    "#F87171",
    "#A78BFA",
    "#AFCBF1",
)
NUCLEOTIDE_COLORS: Dict[str, str] = {
    "A": "#91B4E4",
    "T": "#F472B6",
    "U": "#818CF8",
    "G": "#34D399",
    "C": "#FBBF24",
    "N": "#94A3B8",
}
RESIDUE_TYPE_COLORS: Dict[str, str] = {
    "G": "#34D399",
    "A": "#34D399",
    "P": "#34D399",
    "V": "#34D399",
    "L": "#34D399",
    "I": "#34D399",
    "M": "#34D399",
    "F": "#A78BFA",
    "Y": "#A78BFA",
    "W": "#A78BFA",
    "S": "#FBBF24",
    "T": "#FBBF24",
    "C": "#FBBF24",
    "N": "#FBBF24",
    "Q": "#FBBF24",
    "K": "#6893D0",
    "R": "#6893D0",
    "H": "#6893D0",
    "D": "#F472B6",
    "E": "#F472B6",
    "X": "#64748B",
}
DEFAULT_COLOR: str = "#8AABD6"
SELECTION_COLOR: str = "#D9E4FC"
HIGHLIGHT_COLOR: str = "#AFCBF1"
UNMAPPED_COLOR: str = "#475569"
BLOCKED_KEY_PARTS: Tuple[str, ...] = (
    "url",
    "token",
    "api_key",
    "password",
    "secret",
    "http",
)


class SceneError(Exception):
    """Falha classificada do scene model. Nunca deve virar geometria falsa.

    Attributes:
        category: UNAVAILABLE, ERROR, PARSING_ERROR, RESOURCE_LIMIT, STALE,
            INVALID_INPUT.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "ERROR")


def webgl_capability_note() -> str:
    """Texto honesto: o processo Python nao detecta WebGL pelo user-agent.

    Args:
        Nenhum.

    Returns:
        Frase em ingles.

    Raises:
        Nenhum.
    """
    return (
        "This view uses Plotly WebGL. HelixScope does not detect GPU capability "
        "from the Python process or from the user-agent string. If the canvas "
        "is blank, 3D is unavailable in this browser or environment. Mapping "
        "tables and coordinates remain available."
    )


def visual_color_for_plddt(value: Optional[float]) -> str:
    """Bins de cor AlphaFold pLDDT. Visual apenas; nao e accuracy nem probability.

    Args:
        value: pLDDT depositado (B-factor) ou None.

    Returns:
        Hex visual. Cinza se o valor estiver ausente ou nao finito.

    Raises:
        Nenhum.

    Nota biologica:
        Os intervalos 50/70/90 seguem a convencao visual da AlphaFold DB.
        Nao sao uma probabilidade de que a estrutura esteja correta.
    """
    number = _finite_float(value)
    if number is None:
        return DEFAULT_COLOR
    if number < 50.0:
        return "#EF4444"
    if number < 70.0:
        return "#F59E0B"
    if number < 90.0:
        return "#FACC15"
        return "#6893D0"


def visual_color_from_scale(
    value: Optional[float],
    *,
    vmin: float,
    vmax: float,
    low: str,
    high: str,
) -> str:
    """Interpola uma cor visual. Nao e uma metrica cientifica.

    Args:
        value: Numero ja calculado noutro modulo, ou None.
        vmin: Minimo da escala visual.
        vmax: Maximo da escala visual.
        low: Hex do extremo baixo.
        high: Hex do extremo alto.

    Returns:
        Hex interpolado, ou cinza se value nao for finito.
    """
    number = _finite_float(value)
    if number is None or vmax == vmin:
        return DEFAULT_COLOR
    t = (number - vmin) / (vmax - vmin)
    t = 0.0 if t < 0.0 else 1.0 if t > 1.0 else t
    r1, g1, b1 = _hex_rgb(low)
    r2, g2, b2 = _hex_rgb(high)
    return (
        f"#{int(r1 + (r2 - r1) * t):02X}"
        f"{int(g1 + (g2 - g1) * t):02X}"
        f"{int(b1 + (b2 - b1) * t):02X}"
    )


def make_structure_event(
    kind: str,
    *,
    sequence_hash: str,
    structure_id: object,
    structure_hash: object = None,
    chain_id: object = None,
    model: object = None,
    query_index: object = None,
    label_seq_id: object = None,
    auth_seq_id: object = None,
    insertion_code: object = "",
    atom_name: object = None,
) -> dict:
    """Evento interno com IDs inequívocos. Nao dispara fetch.

    Args:
        kind: Um de STRUCTURE_EVENT_KINDS.
        sequence_hash: Hash da sequencia de analise atual.
        structure_id: Identificador estrutural.
        structure_hash: Hash do ficheiro, quando existir.
        chain_id: Chain.
        model: Modelo NMR.
        query_index: Posicao 0-based da query, ou None.
        label_seq_id: label_seq_id depositado.
        auth_seq_id: auth_seq_id depositado.
        insertion_code: Insercao depositada.
        atom_name: Nome do atomo, quando a selecao for atomica.

    Returns:
        Dict do evento.

    Raises:
        SceneError: INVALID_INPUT se kind for desconhecido.
    """
    if kind not in STRUCTURE_EVENT_KINDS:
        raise SceneError(f"Unknown structure event kind: {kind}.", "INVALID_INPUT")
    return {
        "kind": kind,
        "sequence_hash": str(sequence_hash or ""),
        "structure_id": structure_id,
        "structure_hash": structure_hash,
        "chain_id": chain_id,
        "model": model,
        "query_index_0based": query_index,
        "label_seq_id": label_seq_id,
        "auth_seq_id": auth_seq_id,
        "insertion_code": insertion_code or "",
        "atom_name": atom_name,
    }


def parse_selection_payload(payload: object) -> Optional[dict]:
    """Interpreta customdata de um ponto Plotly. Nao adivinha mapping.

    Args:
        payload: Lista/tupla customdata ou None.

    Returns:
        Dict de selecao ou None.
    """
    if payload is None:
        return None
    values = payload
    if isinstance(payload, (list, tuple)) and payload and isinstance(payload[0], (list, tuple)):
        values = payload[0]
    if not isinstance(values, (list, tuple)) or len(values) < 7:
        return None
    return {
        "query_index_0based": _optional_int(values[0]),
        "chain_id": None if values[1] in (None, "") else str(values[1]),
        "model": _optional_int(values[2]),
        "label_seq_id": _optional_int(values[3]),
        "auth_seq_id": _optional_int(values[4]),
        "insertion_code": "" if values[5] is None else str(values[5]),
        "atom_name": None if values[6] in (None, "") else str(values[6]),
    }


def _backbone_xyz(item: Mapping[str, Any]) -> Optional[dict]:
    """CA de proteina, senao P ou C1' de acido nucleico. Nao inventa XYZ."""
    for key in ("ca", "p", "c1"):
        point = item.get(key) or {}
        if point.get("x") is None or point.get("y") is None or point.get("z") is None:
            continue
        if scientific_checks.atom_coordinate_is_finite(point.get("x"), point.get("y"), point.get("z")):
            return {
                "x": float(point["x"]),
                "y": float(point["y"]),
                "z": float(point["z"]),
                "atom": key.upper() if key != "c1" else "C1'",
            }
    return None


def query_to_structure_residue(envelope: Mapping[str, Any], position: int) -> Optional[dict]:
    """sequence position -> structure residue, quando o mapping existir.

    Args:
        envelope: structure_3d_input.
        position: Indice 0-based.

    Returns:
        Linha de mapping ou None.

    Raises:
        StructureError: INVALID_INPUT se fora do intervalo.
    """
    return protein_structure.query_position_to_residue(envelope, position)


def structure_residue_to_query(
    envelope: Mapping[str, Any],
    *,
    chain_id: str,
    label_seq_id: Optional[int] = None,
    auth_seq_id: Optional[int] = None,
    insertion_code: str = "",
    model: Optional[int] = None,
) -> Optional[dict]:
    """structure residue -> sequence position, quando valido.

    Args:
        envelope: structure_3d_input.
        chain_id: Chain.
        label_seq_id: label_seq_id, se conhecido.
        auth_seq_id: auth_seq_id, se conhecido.
        insertion_code: Insercao.
        model: Modelo.

    Returns:
        Linha de mapping ou None.

    Raises:
        StructureError.
    """
    return protein_structure.structure_residue_to_query(
        envelope,
        chain_id=chain_id,
        label_seq_id=label_seq_id,
        auth_seq_id=auth_seq_id,
        insertion_code=insertion_code,
        model=model,
    )


def highlight_query_indices_for_span(start: int, end: int) -> List[int]:
    """Indices 0-based de um intervalo [start, end) para overlay visual.

    Args:
        start: Inclusivo.
        end: Exclusivo.

    Returns:
        Lista de indices.

    Raises:
        SceneError: INVALID_INPUT se o intervalo for vazio.
    """
    if start < 0 or end <= start:
        raise SceneError("Highlight span is invalid.", "INVALID_INPUT")
    return list(range(start, end))


def scientific_snapshot(envelope: Mapping[str, Any]) -> dict:
    """Campos cientificos que o scene model nao pode alterar.

    Args:
        envelope: structure_3d_input.

    Returns:
        Dict com hashes, coordenadas, chains e proveniencia.

    Raises:
        Nenhum.
    """
    atoms = list(envelope.get("atoms") or [])
    mappings = list(envelope.get("mappings") or [])
    return {
        "sequence_hash": envelope.get("sequence_hash"),
        "structure_hash": envelope.get("structure_hash"),
        "chain_id": envelope.get("chain_id"),
        "model_number": envelope.get("model_number"),
        "mapping_status": envelope.get("mapping_status"),
        "n_atoms": len(atoms),
        "n_mappings": len(mappings),
        "first_xyz": (
            (atoms[0].get("x"), atoms[0].get("y"), atoms[0].get("z")) if atoms else None
        ),
        "provenance": copy.deepcopy(envelope.get("provenance")),
        "first_mapping": copy.deepcopy(mappings[0]) if mappings else None,
    }


def scene_fingerprint(scene: Mapping[str, Any]) -> str:
    """Hash deterministico dos objectos visuais (sem timestamps).

    Args:
        scene: Scene model.

    Returns:
        Hex SHA-256.

    Raises:
        Nenhum.
    """
    payload = {
        "status": scene.get("status"),
        "representation": scene.get("representation"),
        "color_mode": scene.get("color_mode"),
        "objects": scene.get("objects"),
        "view_scope": scene.get("view_scope"),
        "n_atoms_rendered": scene.get("n_atoms_rendered"),
    }
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def structure_3d_report(scene: Mapping[str, Any], envelope: Mapping[str, Any]) -> dict:
    """Resumo 3D sem metricas inventadas.

    Args:
        scene: Scene model.
        envelope: structure_3d_input original.

    Returns:
        Dict JSON-safe.

    Raises:
        Nenhum.
    """
    coverage = envelope.get("coverage_query_with_coordinates")
    payload = {
        "structure_id": envelope.get("structure_id"),
        "source": envelope.get("source"),
        "kind": envelope.get("kind"),
        "kind_label": envelope.get("kind_label") or protein_structure.structure_kind_label(envelope.get("kind")),
        "chain": scene.get("chain_id") or envelope.get("chain_id"),
        "model": scene.get("model_number") or envelope.get("model_number"),
        "mapped_residues": envelope.get("n_with_coordinates"),
        "coordinate_coverage": coverage,
        "coverage_formula": envelope.get("coverage_formula"),
        "mapping_status": envelope.get("mapping_status"),
        "mapping_status_note": envelope.get("mapping_status_note"),
        "view_scope": scene.get("view_scope"),
        "visible_chains": scene.get("visible_chains"),
        "partial_view": scene.get("partial_view"),
        "representation": scene.get("representation"),
        "color_mode": scene.get("color_mode"),
        "n_atoms_source": scene.get("n_atoms_source"),
        "n_atoms_rendered": scene.get("n_atoms_rendered"),
        "engine": ENGINE,
        "sequence_hash": envelope.get("sequence_hash"),
        "structure_hash": envelope.get("structure_hash"),
        "scene_fingerprint": (
            scene_fingerprint(scene) if str(scene.get("status") or "") == "READY" else ""
        ),
        "timings_ms": scene.get("timings_ms"),
        "status": scene.get("status"),
    }
    return provenance.json_safe(payload)


def renderer_log(
    event: str,
    *,
    status: str = "",
    n_atoms: Optional[int] = None,
    representation: str = "",
    color_mode: str = "",
    error_category: str = "",
) -> dict:
    """Evento de diagnostico sem secrets, URLs ou coordenadas.

    Args:
        event: Nome curto (renderer_init, scene_build, error).
        status: Estado do scene.
        n_atoms: Tamanho da estrutura.
        representation: Modo visual.
        color_mode: Modo de cor.
        error_category: Categoria, se houver.

    Returns:
        Dict redigido.

    Raises:
        Nenhum.
    """
    event = provenance.diagnostic_event(
        module="structure_scene",
        operation=event,
        status=status or "info",
        exception=error_category,
    )
    event["n_atoms"] = n_atoms
    event["representation"] = representation
    event["color_mode"] = color_mode
    event["engine"] = ENGINE
    event["browser_capability"] = "not_detected_in_python"
    return provenance.redact_secrets(event)


def build_scene(
    envelope: Mapping[str, Any],
    *,
    expected_sequence_hash: str = "",
    representation: str = "backbone",
    color_mode: str = "chain",
    selected_query_index: Optional[int] = None,
    highlight_query_indices: Optional[Sequence[int]] = None,
    chain_filter: str = "",
    model_filter: Optional[int] = None,
    overlay: Optional[Mapping[int, float]] = None,
    visible_chains: Optional[Sequence[str]] = None,
    lod: str = region_nav.LOD_HIGH,
) -> dict:
    """Constroi o scene model visual a partir do envelope validado.

    Args:
        envelope: Copia isolada de structure_3d_input.
        expected_sequence_hash: Hash da sequencia atual; vazio desliga o gate.
        representation: backbone, ca ou atoms.
        color_mode: chain, residue_type, selection, plddt, hydropathy,
            conservation, charge.
        selected_query_index: Residuo destacado (0-based).
        highlight_query_indices: Outros indices a destacar.
        chain_filter: Chain de mapping/foco; nao concatena chains.
        model_filter: Se definido, so esse modelo NMR.
        overlay: Valores por query index para hydropathy/conservation/charge.
        visible_chains: Chains a desenhar. Subconjunto implica Partial structure view.
        lod: low/medium/high. Visual stride only; omitted residues are not interpolated.

    Returns:
        Scene model. status READY so quando ha coordenadas finitas para desenhar.

    Raises:
        Nenhum. Estados cientificos voltam no campo status.

    Nota biologica:
        A transformacao e visual. sequence_hash, structure_hash, coordenadas,
        chain IDs, mapping e proveniencia do envelope nao sao reescritos.
        Residuos sem coordenadas nao recebem atomos inventados.
        LOD baixo omite pontos do desenho; nao inventa XYZ.
    """
    source = copy.deepcopy(dict(envelope or {}))
    snapshot = scientific_snapshot(source)
    started = time.perf_counter()
    try:
        scene = _build_scene_checked(
            source,
            expected_sequence_hash=expected_sequence_hash,
            representation=representation,
            color_mode=color_mode,
            selected_query_index=selected_query_index,
            highlight_query_indices=highlight_query_indices,
            chain_filter=chain_filter,
            model_filter=model_filter,
            overlay=overlay,
            visible_chains=visible_chains,
            snapshot=snapshot,
            lod=lod,
        )
    except SceneError as exc:
        empty = _empty_scene(source, status=exc.category, message=str(exc), snapshot=snapshot)
        empty["diagnostics"] = renderer_log(
            "scene_build",
            status=exc.category,
            n_atoms=len(list(source.get("atoms") or [])),
            representation=representation,
            color_mode=color_mode,
            error_category=exc.category,
        )
        empty["timings_ms"] = {
            "scene_construction": round((time.perf_counter() - started) * 1000.0, 3)
        }
        return empty
    scene["diagnostics"] = renderer_log(
        "scene_build",
        status=str(scene.get("status") or ""),
        n_atoms=scene.get("n_atoms_source"),
        representation=str(scene.get("representation") or ""),
        color_mode=str(scene.get("color_mode") or ""),
    )
    scene["timings_ms"] = {
        "scene_construction": round((time.perf_counter() - started) * 1000.0, 3)
    }
    return scene


def _build_scene_checked(
    source: dict,
    *,
    expected_sequence_hash: str,
    representation: str,
    color_mode: str,
    selected_query_index: Optional[int],
    highlight_query_indices: Optional[Sequence[int]],
    chain_filter: str,
    model_filter: Optional[int],
    overlay: Optional[Mapping[int, float]],
    visible_chains: Optional[Sequence[str]],
    snapshot: dict,
    lod: str = region_nav.LOD_HIGH,
) -> dict:
    atoms = list(source.get("atoms") or [])
    residues = list(source.get("residues") or [])
    if len(atoms) > protein_structure.MAX_ATOMS:
        raise SceneError(
            f"Structure exceeds {protein_structure.MAX_ATOMS:,} atoms (RESOURCE_LIMIT).",
            "RESOURCE_LIMIT",
        )
    if len(residues) > protein_structure.MAX_RESIDUES:
        raise SceneError(
            f"Structure exceeds {protein_structure.MAX_RESIDUES:,} residues (RESOURCE_LIMIT).",
            "RESOURCE_LIMIT",
        )
    _assert_no_network_fields(source)
    if expected_sequence_hash and str(source.get("sequence_hash") or "") != str(expected_sequence_hash):
        raise SceneError(
            "stale/incompatible structure: current sequence hash does not match "
            "the structure sequence hash.",
            "STALE",
        )
    rep = str(representation or "backbone").strip().lower()
    if rep not in REPRESENTATIONS:
        raise SceneError(f"Unknown representation: {representation}.", "INVALID_INPUT")
    requested_color = str(color_mode or "chain").strip().lower()
    if requested_color not in COLOR_MODES:
        raise SceneError(f"Unknown color mode: {color_mode}.", "INVALID_INPUT")

    chain_ids = sorted(
        {
            str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
            for atom in atoms
            if str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        }
    )
    model_ids = sorted({int(atom.get("model") or 1) for atom in atoms} if atoms else [])
    if len(chain_ids) > protein_structure.MAX_RESIDUES:
        raise SceneError("Chain count exceeds the structural resource limit.", "RESOURCE_LIMIT")
    if len(model_ids) > protein_structure.MAX_RESIDUES:
        raise SceneError("Model count exceeds the structural resource limit.", "RESOURCE_LIMIT")

    for atom in atoms:
        if not scientific_checks.atom_coordinate_is_finite(atom.get("x"), atom.get("y"), atom.get("z")):
            raise SceneError("Structure coordinates are not finite numbers.", "PARSING_ERROR")

    molecule = str(source.get("molecule_type") or "").strip().upper()
    selected_model = int(model_filter) if model_filter is not None else int(source.get("model_number") or 1)
    if chain_filter:
        selected_chain = str(chain_filter)
    elif molecule in {"DNA", "RNA", "COMPLEX"}:
        selected_chain = ""
    else:
        selected_chain = str(source.get("chain_id") or "")
    if visible_chains:
        allowed_chains = {str(item) for item in visible_chains if str(item)}
    elif selected_chain:
        allowed_chains = {selected_chain}
    else:
        allowed_chains = None
    filtered = [
        atom
        for atom in atoms
        if str(atom.get("group") or "ATOM") == "ATOM"
        and (
            allowed_chains is None
            or str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "") in allowed_chains
        )
        and int(atom.get("model") or 1) == selected_model
    ]
    if not filtered:
        raise SceneError("Coordinates unavailable for this chain or model.", "UNAVAILABLE")

    mappings = list(source.get("mappings") or [])
    if molecule in {"DNA", "RNA", "COMPLEX"}:
        mappings = mappings + _unmapped_nucleic_chain_rows(
            atoms,
            mappings,
            model_number=selected_model,
        )
    query_by_residue = _query_index_lookup(mappings)
    highlights = {int(index) for index in (highlight_query_indices or [])}
    if selected_query_index is not None:
        highlights.add(int(selected_query_index))

    kind = str(source.get("kind") or "")
    applied_color, color_note = _resolve_color_mode(
        requested_color,
        kind=kind,
        overlay=overlay,
        mappings=mappings,
    )

    backbone = _backbone_traces(
        mappings,
        chain_id=selected_chain,
        visible_chains=allowed_chains,
        model_number=selected_model,
        color_mode=applied_color,
        highlights=highlights,
        selected_query_index=selected_query_index,
        overlay=overlay,
        kind=kind,
    )
    try:
        visual_lod = region_nav.normalize_lod(lod)
        stride = region_nav.backbone_stride(visual_lod)
    except region_nav.RegionError as exc:
        raise SceneError(str(exc), "INVALID_INPUT") from exc
    lod_partial = False
    if stride > 1:
        backbone = [_stride_backbone_trace(trace, stride) for trace in backbone]
        lod_partial = True
    pair_traces = _pair_traces(list(source.get("pairs") or []))
    n_pairs_drawn = 0
    if pair_traces:
        n_pairs_drawn = int(pair_traces[0].get("n_pairs") or len(pair_traces))
    termini = _terminus_markers(
        mappings,
        visible_chains=allowed_chains,
        model_number=selected_model,
        kind=kind,
    )
    atom_points = None
    if rep in {"ca", "atoms"}:
        pool = filtered if rep == "atoms" else [atom for atom in filtered if str(atom.get("atom_name") or "") == "CA"]
        atom_points = _atom_points(
            pool,
            query_by_residue=query_by_residue,
            mappings=mappings,
            color_mode=applied_color,
            highlights=highlights,
            selected_query_index=selected_query_index,
            overlay=overlay,
            kind=kind,
        )

    selection_marker = _selection_marker(
        mappings,
        chain_id=selected_chain,
        model_number=selected_model,
        selected_query_index=selected_query_index,
        visible_chains=allowed_chains,
    )

    rendered_xyz = _collect_xyz(backbone, atom_points)
    if not rendered_xyz:
        raise SceneError("Coordinates unavailable for this residue set.", "UNAVAILABLE")

    source_chains_model = sorted(
        {
            str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
            for atom in atoms
            if int(atom.get("model") or 1) == selected_model
            and str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        }
    )
    partial = False
    if allowed_chains is not None and set(source_chains_model) - allowed_chains:
        partial = True
    if selected_chain and len(source_chains_model) > 1 and (
        allowed_chains is None or allowed_chains != set(source_chains_model)
    ):
        partial = True
    if model_filter is not None and len(model_ids) > 1:
        partial = True
    if lod_partial:
        partial = True
    if bool(source.get("partial_view")):
        partial = True
    view_scope = "partial_structure_view" if partial else "selected_chain"
    n_atoms_rendered = 0 if atom_points is None else len(atom_points.get("x") or [])
    n_ca = sum(len(trace.get("x") or []) for trace in backbone)
    missing = [
        int(item.get("query_index_0based"))
        for item in mappings
        if not item.get("has_coordinates")
    ]
    camera = _camera_from_points(rendered_xyz, selected_query_index, mappings, selected_chain, selected_model)
    mapping_status = source.get("mapping_status")
    if mapping_status and not scientific_checks.mapping_status_is_declared(mapping_status):
        mapping_status = "UNCERTAIN"
    return {
        "engine": ENGINE,
        "status": "READY",
        "message": "",
        "source": source.get("source"),
        "kind": kind,
        "kind_label": source.get("kind_label") or protein_structure.structure_kind_label(kind),
        "structure_id": source.get("structure_id"),
        "sequence_hash": source.get("sequence_hash"),
        "structure_hash": source.get("structure_hash"),
        "chain_id": selected_chain,
        "model_number": selected_model,
        "models": list(source.get("models") or model_ids),
        "mapping_status": mapping_status,
        "mapping_status_note": source.get("mapping_status_note"),
        "coverage_query_with_coordinates": source.get("coverage_query_with_coordinates"),
        "representation": rep,
        "representation_note": _representation_note(rep),
        "color_mode": applied_color,
        "color_note": color_note,
        "requested_color_mode": requested_color,
        "view_scope": view_scope,
        "view_label": "Partial structure view" if partial else "Selected chain",
        "partial_view": partial,
        "visual_lod": visual_lod,
        "visual_stride": stride,
        "objects": {
            "backbone_traces": backbone,
            "pair_traces": pair_traces,
            "atom_points": atom_points,
            "selection_marker": selection_marker,
            "terminus_markers": termini,
        },
        "materials": {
            "marker_size": _visual_marker_size(rep, visual_lod, molecule),
            "line_width": _visual_line_width(visual_lod, molecule),
            "opacity": 0.98 if molecule in {"DNA", "RNA", "COMPLEX"} else 0.95,
            "note": "Materials are visual encodings, not experimental measurements. LOD changes marker/line size and stride only; deposited XYZ are unchanged.",
        },
        "selection_ids": highlights,
        "selected_query_index": selected_query_index,
        "residues_without_coordinates": missing,
        "legend": _legend(applied_color, backbone, chain_ids),
        "camera": camera,
        "n_atoms_source": len(atoms),
        "n_atoms_rendered": n_atoms_rendered,
        "n_backbone_points": n_ca,
        "n_pair_traces": n_pairs_drawn,
        "n_chains_in_file": len(source_chains_model),
        "n_models_in_file": len(model_ids),
        "visible_chains": sorted(allowed_chains) if allowed_chains is not None else list(source_chains_model),
        "lod_note": (
            f"Visual LOD {visual_lod} draws every {stride} mapped backbone "
            "point(s). Omitted residues keep deposited coordinates in the "
            "scientific envelope and are not interpolated."
            if stride > 1
            else (
                f"Atom representation of {len(atoms):,} deposited ATOM records. "
                "Switch to backbone for a lighter view. Coordinates are not approximated."
                if rep == "atoms" and len(atoms) >= 2000
                else ""
            )
        ),
        "alt_location_policy": (source.get("metadata") or {}).get("alt_location_policy"),
        "model_note": (source.get("metadata") or {}).get("model_note"),
        "disclaimer": (source.get("metadata") or {}).get("disclaimer"),
        "webgl_note": webgl_capability_note(),
        "source_snapshot": snapshot,
        "renderer": ENGINE,
    }


def _empty_scene(source: Mapping[str, Any], *, status: str, message: str, snapshot: dict) -> dict:
    return {
        "engine": ENGINE,
        "status": status,
        "message": message,
        "source": source.get("source"),
        "kind": source.get("kind"),
        "kind_label": source.get("kind_label") or protein_structure.structure_kind_label(source.get("kind")),
        "structure_id": source.get("structure_id"),
        "sequence_hash": source.get("sequence_hash"),
        "structure_hash": source.get("structure_hash"),
        "chain_id": source.get("chain_id"),
        "model_number": source.get("model_number"),
        "mapping_status": source.get("mapping_status"),
        "objects": {"backbone_traces": [], "pair_traces": [], "atom_points": None, "selection_marker": None},
        "view_scope": "none",
        "view_label": "Coordinates Unavailable" if status == "UNAVAILABLE" else status,
        "partial_view": False,
        "visual_lod": region_nav.LOD_HIGH,
        "visual_stride": 1,
        "n_atoms_source": len(list(source.get("atoms") or [])),
        "n_atoms_rendered": 0,
        "n_backbone_points": 0,
        "residues_without_coordinates": list(source.get("residues_without_coordinates") or []),
        "source_snapshot": snapshot,
        "webgl_note": webgl_capability_note(),
        "renderer": ENGINE,
    }


def _assert_no_network_fields(value: Any, *, path: str = "") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            low = str(key).lower()
            if any(part in low for part in BLOCKED_KEY_PARTS):
                raise SceneError("Structure envelope contains blocked network fields.", "ERROR")
            _assert_no_network_fields(item, path=f"{path}.{key}")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _assert_no_network_fields(item, path=f"{path}[{index}]")
        return
    if isinstance(value, str) and value.lower().startswith(("http://", "https://")):
        raise SceneError("Structure envelope contains a URL value.", "ERROR")


def _resolve_color_mode(
    requested: str,
    *,
    kind: str,
    overlay: Optional[Mapping[int, float]],
    mappings: Sequence[Mapping[str, Any]],
) -> Tuple[str, str]:
    if requested == "plddt":
        if kind != "predicted":
            return "chain", "pLDDT coloring applies only to predicted models. Chain color is shown."
        if not any(_finite_float(item.get("ca_b_iso")) is not None for item in mappings):
            return "chain", "pLDDT values are not present on this envelope. Chain color is shown."
        return "plddt", "Color uses deposited pLDDT (B-factor convention). Not a probability or accuracy percentage."
    if requested == "hydropathy" or requested == "conservation" or requested == "charge":
        if not overlay:
            return "chain", f"{requested} overlay was not provided. Chain color is shown."
        if requested == "hydropathy":
            return "hydropathy", "Color uses Kyte-Doolittle hydropathy from the Protein engine. Not a structural measurement."
        if requested == "charge":
            return (
                "charge",
                "Formal side-chain charge from the Protein engine (K/R +1, D/E -1, H excluded). "
                "Not net_charge and not Henderson-Hasselbalch at a stated pH.",
            )
        return "conservation", "Color uses MSA conservation among analyzed sequences. Not a functional annotation."
    if requested == "selection":
        return "selection", "Unselected residues are muted. Selection color is visual, not a scientific score."
    if requested == "residue_type":
        return "residue_type", "Residue-type colors follow side-chain chemistry groups. Visual only."
    return "chain", "Each chain has a distinct visual color. Color is not experimental vs predicted."


def _query_index_lookup(mappings: Sequence[Mapping[str, Any]]) -> Dict[tuple, int]:
    lookup: Dict[tuple, int] = {}
    for item in mappings:
        query_index = item.get("query_index_0based")
        if query_index is None:
            continue
        key = (
            str(item.get("chain_id") or ""),
            int(item.get("model") or 1),
            item.get("label_seq_id"),
            item.get("auth_seq_id"),
            str(item.get("insertion_code") or ""),
        )
        lookup[key] = int(query_index)
    return lookup


def _backbone_traces(
    mappings: Sequence[Mapping[str, Any]],
    *,
    chain_id: str,
    model_number: int,
    color_mode: str,
    highlights: set,
    selected_query_index: Optional[int],
    overlay: Optional[Mapping[int, float]],
    kind: str,
    visible_chains: Optional[set] = None,
) -> List[dict]:
    rows = [
        item
        for item in mappings
        if (not chain_id or str(item.get("chain_id") or "") == chain_id)
        and (visible_chains is None or str(item.get("chain_id") or "") in visible_chains)
        and (item.get("model") is None or int(item.get("model") or 1) == model_number)
    ]
    rows.sort(
        key=lambda item: (
            str(item.get("chain_id") or ""),
            int(item.get("model") or 1),
            int(item.get("chain_index_0based") if item.get("chain_index_0based") is not None else 10**9),
        )
    )
    traces: List[dict] = []
    current: List[Mapping[str, Any]] = []
    last_index: Optional[int] = None
    last_chain = ""
    for item in rows:
        point = _backbone_xyz(item)
        if not item.get("has_coordinates") or point is None:
            if current:
                traces.append(_flush_backbone(current, color_mode, highlights, selected_query_index, overlay, kind))
                current = []
            last_index = None
            continue
        chain = str(item.get("chain_id") or "")
        polymer_index = item.get("chain_index_0based")
        if current and (chain != last_chain or last_index is None or polymer_index != last_index + 1):
            traces.append(_flush_backbone(current, color_mode, highlights, selected_query_index, overlay, kind))
            current = []
        current.append(item)
        last_index = polymer_index if polymer_index is None else int(polymer_index)
        last_chain = chain
    if current:
        traces.append(_flush_backbone(current, color_mode, highlights, selected_query_index, overlay, kind))
    return traces


def _stride_backbone_trace(trace: Mapping[str, Any], stride: int) -> dict:
    """Omite pontos do desenho. Nao interpola XYZ."""
    xs = list(trace.get("x") or [])
    keep = region_nav.downsample_indices(len(xs), stride)
    copied = dict(trace)
    for key in ("x", "y", "z", "colors", "customdata", "hovertext"):
        values = list(trace.get(key) or [])
        if values and len(values) == len(xs):
            copied[key] = [values[i] for i in keep]
    copied["visual_stride"] = stride
    return copied


def _pair_traces(pairs: Sequence[Mapping[str, Any]]) -> List[dict]:
    """Uma trace visual com cortes None; XYZ so dos pares ja presentes."""
    xs: List[Optional[float]] = []
    ys: List[Optional[float]] = []
    zs: List[Optional[float]] = []
    hover: List[str] = []
    n_pairs = 0
    for item in pairs:
        plus = item.get("plus") or {}
        minus = item.get("minus") or {}
        if not scientific_checks.atom_coordinate_is_finite(plus.get("x"), plus.get("y"), plus.get("z")):
            continue
        if not scientific_checks.atom_coordinate_is_finite(minus.get("x"), minus.get("y"), minus.get("z")):
            continue
        if n_pairs:
            xs.append(None)
            ys.append(None)
            zs.append(None)
            hover.append("")
        xs.extend([float(plus["x"]), float(minus["x"])])
        ys.extend([float(plus["y"]), float(minus["y"])])
        zs.extend([float(plus["z"]), float(minus["z"])])
        label = (
            f"base pair {item.get('plus_base')}-{item.get('minus_base')} "
            f"query_index {item.get('plus_index_0based')} "
            f"(identities from the analysis sequence, not a new structure)"
        )
        hover.extend([label, label])
        n_pairs += 1
    if not n_pairs:
        return []
    return [
        {
            "x": xs,
            "y": ys,
            "z": zs,
            "color": "#94A3B8",
            "name": "base pairs",
            "n_pairs": n_pairs,
            "hovertext": hover,
        }
    ]


def _unmapped_nucleic_chain_rows(
    atoms: Sequence[Mapping[str, Any]],
    mappings: Sequence[Mapping[str, Any]],
    *,
    model_number: int,
) -> List[dict]:
    """Filas de backbone para cadeias nucleicas depositadas sem query_index.

    Usa atomos P (ou C1') ja presentes. Nao inventa XYZ nem identidades.
    """
    mapped_chains = {
        str(item.get("chain_id") or "")
        for item in mappings
        if str(item.get("chain_id") or "")
    }
    grouped: Dict[str, List[Mapping[str, Any]]] = {}
    for atom in atoms:
        if int(atom.get("model") or 1) != int(model_number):
            continue
        chain = str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        if not chain or chain in mapped_chains:
            continue
        name = str(atom.get("atom_name") or "")
        if name not in {"P", "C1'"}:
            continue
        grouped.setdefault(chain, []).append(atom)
    extra: List[dict] = []
    for chain, pool in grouped.items():
        preferred = [atom for atom in pool if str(atom.get("atom_name") or "") == "P"]
        use = preferred or pool
        use.sort(
            key=lambda atom: (
                int(atom.get("auth_seq_id") or atom.get("label_seq_id") or 10**9),
                str(atom.get("atom_name") or ""),
            )
        )
        seen: set = set()
        local = 0
        for atom in use:
            seq_id = atom.get("auth_seq_id")
            key = (chain, seq_id, atom.get("insertion_code") or "")
            if key in seen:
                continue
            seen.add(key)
            if str(atom.get("atom_name") or "") != "P" and preferred:
                continue
            extra.append(
                {
                    "query_index_0based": None,
                    "query_residue": None,
                    "chain_index_0based": local,
                    "chain_residue": _residue_letter({"comp_id": atom.get("comp_id")}),
                    "chain_id": chain,
                    "label_seq_id": atom.get("label_seq_id"),
                    "auth_seq_id": atom.get("auth_seq_id"),
                    "insertion_code": atom.get("insertion_code") or "",
                    "comp_id": atom.get("comp_id"),
                    "has_coordinates": True,
                    "coordinate_status": "available",
                    "ca": None,
                    "p": {"x": atom.get("x"), "y": atom.get("y"), "z": atom.get("z"), "atom": "P"},
                    "c1": None,
                    "model": int(atom.get("model") or 1),
                    "molecule_type": "DNA",
                }
            )
            local += 1
    return extra


def _residue_letter(item: Mapping[str, Any]) -> str:
    """Letra da sequencia/mapping. Nao infere residuo a partir so do CCD D*."""
    for key in ("query_residue", "chain_residue"):
        raw = str(item.get(key) or "").strip().upper()
        if len(raw) == 1 and raw.isalpha():
            return raw
    comp = str(item.get("comp_id") or "").strip().upper()
    if len(comp) == 1 and comp.isalpha():
        return comp
    if len(comp) == 2 and comp[0] in {"D", "R"} and comp[1] in "ACGTU":
        return comp[1]
    return ""


def _identity_color(item: Mapping[str, Any]) -> str:
    letter = _residue_letter(item)
    if letter in NUCLEOTIDE_COLORS:
        return NUCLEOTIDE_COLORS[letter]
    if letter in RESIDUE_TYPE_COLORS:
        return RESIDUE_TYPE_COLORS[letter]
    return DEFAULT_COLOR


def _terminus_markers(
    mappings: Sequence[Mapping[str, Any]],
    *,
    visible_chains: Optional[set],
    model_number: int,
    kind: str,
) -> List[dict]:
    """Marcadores 5'/3' nos extremos reais de cada cadeia com coordenadas."""
    grouped: Dict[str, List[Mapping[str, Any]]] = {}
    for item in mappings:
        if not item.get("has_coordinates"):
            continue
        if item.get("model") is not None and int(item.get("model") or 1) != int(model_number):
            continue
        chain = str(item.get("chain_id") or "")
        if visible_chains is not None and chain not in visible_chains:
            continue
        point = _backbone_xyz(item)
        if point is None:
            continue
        grouped.setdefault(chain, []).append(item)
    traces: List[dict] = []
    for chain, rows in grouped.items():
        rows.sort(
            key=lambda item: int(
                item.get("chain_index_0based") if item.get("chain_index_0based") is not None else 10**9
            )
        )
        ends = (("5-prime", rows[0]), ("3-prime", rows[-1])) if len(rows) >= 1 else ()
        if len(rows) == 1:
            ends = (("terminus", rows[0]),)
        xs: List[float] = []
        ys: List[float] = []
        zs: List[float] = []
        hover: List[str] = []
        for label, item in ends:
            point = _backbone_xyz(item)
            if point is None:
                continue
            xs.append(float(point["x"]))
            ys.append(float(point["y"]))
            zs.append(float(point["z"]))
            hover.append(
                f"{label} | chain {chain} | residue {_residue_letter(item) or item.get('comp_id') or 'N/A'} "
                f"| query_index {item.get('query_index_0based')} | kind {kind}"
            )
        if xs:
            traces.append(
                {
                    "x": xs,
                    "y": ys,
                    "z": zs,
                    "hovertext": hover,
                    "name": f"termini {chain}",
                    "chain_id": chain,
                }
            )
    return traces


def _flush_backbone(
    rows: Sequence[Mapping[str, Any]],
    color_mode: str,
    highlights: set,
    selected_query_index: Optional[int],
    overlay: Optional[Mapping[int, float]],
    kind: str,
) -> dict:
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    colors: List[str] = []
    custom: List[list] = []
    hover: List[str] = []
    chain = str(rows[0].get("chain_id") or "")
    model = rows[0].get("model")
    chain_color_index = sum(ord(char) for char in chain) % len(CHAIN_PALETTE) if chain else 0
    line_color = CHAIN_PALETTE[chain_color_index] if chain else DEFAULT_COLOR
    for item in rows:
        point = _backbone_xyz(item)
        if point is None:
            continue
        query_index = item.get("query_index_0based")
        xs.append(float(point["x"]))
        ys.append(float(point["y"]))
        zs.append(float(point["z"]))
        marker = _color_for_residue(
            item,
            color_mode=color_mode,
            highlights=highlights,
            selected_query_index=selected_query_index,
            overlay=overlay,
            kind=kind,
        )
        if color_mode == "chain":
            marker = _identity_color(item) if _residue_letter(item) else line_color
            if query_index is not None and int(query_index) in highlights:
                if selected_query_index is not None and int(query_index) == int(selected_query_index):
                    marker = HIGHLIGHT_COLOR
                elif color_mode != "selection":
                    marker = HIGHLIGHT_COLOR
        colors.append(marker)
        atom_name = str(point.get("atom") or "CA")
        custom.append(
            [
                query_index,
                chain,
                model,
                item.get("label_seq_id"),
                item.get("auth_seq_id"),
                item.get("insertion_code") or "",
                atom_name,
            ]
        )
        hover.append(_hover_text(item, atom_name=atom_name, xyz=point, kind=kind))
    role = str(rows[0].get("role") or "").strip()
    if role:
        name = f"{role.replace('_', ' ')} | chain {chain} model {model}"
    else:
        name = f"chain {chain} model {model}"
    return {
        "chain_id": chain,
        "model": model,
        "role": role or None,
        "x": xs,
        "y": ys,
        "z": zs,
        "color": line_color,
        "colors": colors,
        "customdata": custom,
        "hovertext": hover,
        "name": name,
    }


def _selection_marker(
    mappings: Sequence[Mapping[str, Any]],
    *,
    chain_id: str,
    model_number: int,
    selected_query_index: Optional[int],
    visible_chains: Optional[set] = None,
) -> Optional[dict]:
    if selected_query_index is None:
        return None
    for item in mappings:
        if item.get("query_index_0based") is None or int(item.get("query_index_0based")) != int(selected_query_index):
            continue
        if chain_id and str(item.get("chain_id") or "") != chain_id:
            continue
        if visible_chains is not None and str(item.get("chain_id") or "") not in visible_chains:
            continue
        if item.get("model") is not None and int(item.get("model") or 1) != model_number:
            continue
        point = _backbone_xyz(item)
        if not item.get("has_coordinates") or point is None:
            return None
        return {
            "x": [float(point["x"])],
            "y": [float(point["y"])],
            "z": [float(point["z"])],
            "customdata": [
                [
                    selected_query_index,
                    item.get("chain_id"),
                    item.get("model"),
                    item.get("label_seq_id"),
                    item.get("auth_seq_id"),
                    item.get("insertion_code") or "",
                    str(point.get("atom") or "CA"),
                ]
            ],
            "hovertext": [_hover_text(item, atom_name=str(point.get("atom") or "CA"), xyz=point)],
            "name": "selection",
        }
    return None


def _atom_points(
    atoms: Sequence[Mapping[str, Any]],
    *,
    query_by_residue: Mapping[tuple, int],
    mappings: Sequence[Mapping[str, Any]],
    color_mode: str,
    highlights: set,
    selected_query_index: Optional[int],
    overlay: Optional[Mapping[int, float]],
    kind: str,
) -> dict:
    mapping_by_query = {
        int(item.get("query_index_0based")): item
        for item in mappings
        if item.get("query_index_0based") is not None
    }
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    colors: List[str] = []
    sizes: List[float] = []
    custom: List[list] = []
    hover: List[str] = []
    for atom in atoms:
        chain = str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        model = int(atom.get("model") or 1)
        key = (
            chain,
            model,
            atom.get("label_seq_id"),
            atom.get("auth_seq_id"),
            str(atom.get("insertion_code") or ""),
        )
        query_index = query_by_residue.get(key)
        mapped = mapping_by_query.get(int(query_index)) if query_index is not None else None
        residue_row = mapped or {
            "query_index_0based": query_index,
            "query_residue": None,
            "chain_id": chain,
            "label_seq_id": atom.get("label_seq_id"),
            "auth_seq_id": atom.get("auth_seq_id"),
            "insertion_code": atom.get("insertion_code") or "",
            "comp_id": atom.get("comp_id"),
            "ca_b_iso": atom.get("b_iso"),
            "model": model,
            "has_coordinates": True,
        }
        xs.append(float(atom["x"]))
        ys.append(float(atom["y"]))
        zs.append(float(atom["z"]))
        selected = query_index is not None and int(query_index) == selected_query_index
        colors.append(
            _color_for_residue(
                residue_row,
                color_mode=color_mode,
                highlights=highlights,
                selected_query_index=selected_query_index,
                overlay=overlay,
                kind=kind,
            )
        )
        sizes.append(8.0 if selected else 4.0)
        custom.append(
            [
                query_index,
                chain,
                model,
                atom.get("label_seq_id"),
                atom.get("auth_seq_id"),
                atom.get("insertion_code") or "",
                atom.get("atom_name"),
            ]
        )
        hover.append(
            _hover_text(
                residue_row,
                atom_name=str(atom.get("atom_name") or ""),
                xyz=atom,
                element=str(atom.get("element") or ""),
                occupancy=atom.get("occupancy"),
            )
        )
    return {
        "x": xs,
        "y": ys,
        "z": zs,
        "colors": colors,
        "sizes": sizes,
        "customdata": custom,
        "hovertext": hover,
        "name": "atoms",
    }


def _color_for_residue(
    item: Mapping[str, Any],
    *,
    color_mode: str,
    highlights: set,
    selected_query_index: Optional[int],
    overlay: Optional[Mapping[int, float]],
    kind: str,
) -> str:
    query_index = item.get("query_index_0based")
    if query_index is not None and int(query_index) in highlights and color_mode != "selection":
        if selected_query_index is not None and int(query_index) == int(selected_query_index):
            return HIGHLIGHT_COLOR
    if color_mode == "selection":
        if query_index is not None and int(query_index) in highlights:
            return HIGHLIGHT_COLOR
        return UNMAPPED_COLOR
    if color_mode == "residue_type":
        letter = _residue_letter(item)
        if letter in NUCLEOTIDE_COLORS:
            return NUCLEOTIDE_COLORS[letter]
        return RESIDUE_TYPE_COLORS.get(letter, DEFAULT_COLOR)
    if color_mode == "plddt" and kind == "predicted":
        return visual_color_for_plddt(_finite_float(item.get("ca_b_iso")))
    if color_mode in {"hydropathy", "conservation", "charge"} and overlay is not None and query_index is not None:
        value = overlay.get(int(query_index))
        if color_mode == "hydropathy":
            return visual_color_from_scale(value, vmin=-4.5, vmax=4.5, low="#6893D0", high="#F97316")
        if color_mode == "charge":
            return visual_color_from_scale(value, vmin=-1.0, vmax=1.0, low="#F472B6", high="#91B4E4")
        return visual_color_from_scale(value, vmin=0.0, vmax=1.0, low="#1E293B", high="#34D399")
    chain = str(item.get("chain_id") or "")
    if not chain:
        return DEFAULT_COLOR
    index = sum(ord(char) for char in chain) % len(CHAIN_PALETTE)
    return CHAIN_PALETTE[index]


def _hover_text(
    item: Mapping[str, Any],
    *,
    atom_name: str,
    xyz: Mapping[str, Any],
    element: str = "",
    occupancy: object = None,
    kind: str = "",
) -> str:
    query_index = item.get("query_index_0based")
    label = item.get("label_seq_id")
    auth = item.get("auth_seq_id")
    occ = occupancy if occupancy is not None else None
    letter = _residue_letter(item)
    parts = [
        f"atom {atom_name or 'N/A'}",
        f"element {element or 'N/A'}",
        f"residue {letter or item.get('comp_id') or 'N/A'}",
        f"chain {item.get('chain_id') or 'N/A'}",
        f"role {item.get('role') or 'N/A'}",
        f"label_seq_id {label if label is not None else 'N/A'}",
        f"auth_seq_id {auth if auth is not None else 'N/A'}",
        f"model {item.get('model') if item.get('model') is not None else 'N/A'}",
        f"xyz {_fmt_coord(xyz.get('x'))}, {_fmt_coord(xyz.get('y'))}, {_fmt_coord(xyz.get('z'))}",
        f"occupancy {occ if occ is not None else 'N/A'}",
        f"query_index {query_index if query_index is not None else 'N/A'}",
        f"kind {kind or 'N/A'}",
    ]
    return " | ".join(parts)


def _visual_marker_size(representation: str, visual_lod: str, molecule: str) -> int:
    """Tamanho de marcador visual. Nao altera coordenadas."""
    if representation == "atoms":
        return 3 if visual_lod == region_nav.LOD_LOW else 4
    nucleic = molecule in {"DNA", "RNA", "COMPLEX"}
    if visual_lod == region_nav.LOD_LOW:
        return 5 if nucleic else 4
    if visual_lod == region_nav.LOD_MEDIUM:
        return 7 if nucleic else 6
    return 10 if nucleic else 8


def _visual_line_width(visual_lod: str, molecule: str) -> int:
    """Espessura de linha visual. Nao altera coordenadas."""
    nucleic = molecule in {"DNA", "RNA", "COMPLEX"}
    if visual_lod == region_nav.LOD_LOW:
        return 3 if nucleic else 2
    if visual_lod == region_nav.LOD_MEDIUM:
        return 5 if nucleic else 4
    return 8 if nucleic else 6


def _fmt_coord(value: object) -> str:
    number = _finite_float(value)
    if number is None:
        return "N/A"
    return f"{number:.3f}"


def _collect_xyz(backbone: Sequence[Mapping[str, Any]], atom_points: Optional[Mapping[str, Any]]) -> List[Tuple[float, float, float]]:
    points: List[Tuple[float, float, float]] = []
    for trace in backbone:
        for x, y, z in zip(trace.get("x") or [], trace.get("y") or [], trace.get("z") or []):
            points.append((float(x), float(y), float(z)))
    if atom_points:
        for x, y, z in zip(atom_points.get("x") or [], atom_points.get("y") or [], atom_points.get("z") or []):
            points.append((float(x), float(y), float(z)))
    return points


def _camera_from_points(
    points: Sequence[Tuple[float, float, float]],
    selected_query_index: Optional[int],
    mappings: Sequence[Mapping[str, Any]],
    chain_id: str,
    model_number: int,
) -> dict:
    xs = [item[0] for item in points]
    ys = [item[1] for item in points]
    zs = [item[2] for item in points]
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs), 1.0)
    center = {"x": 0.0, "y": 0.0, "z": 0.0}
    if selected_query_index is not None:
        for item in mappings:
            if item.get("query_index_0based") is None or int(item.get("query_index_0based")) != int(selected_query_index):
                continue
            if chain_id and str(item.get("chain_id") or "") != chain_id:
                continue
            if item.get("model") is not None and int(item.get("model") or 1) != model_number:
                continue
            point = _backbone_xyz(item)
            if point is not None:
                center = {
                    "x": (float(point["x"]) - cx) / span,
                    "y": (float(point["y"]) - cy) / span,
                    "z": (float(point["z"]) - cz) / span,
                }
                break
    return {
        "center": center,
        "eye": {"x": 1.6, "y": 1.6, "z": 1.25},
        "up": {"x": 0.0, "y": 0.0, "z": 1.0},
        "centroid": {"x": cx, "y": cy, "z": cz},
        "span": span,
        "note": "Camera is a view transform. It does not reindex scientific coordinates.",
    }


def _legend(color_mode: str, backbone: Sequence[Mapping[str, Any]], chain_ids: Sequence[str]) -> List[dict]:
    if color_mode == "plddt":
        return [
            {"label": "pLDDT < 50", "color": "#EF4444"},
            {"label": "pLDDT 50-70", "color": "#F59E0B"},
            {"label": "pLDDT 70-90", "color": "#FACC15"},
            {"label": "pLDDT > 90", "color": "#6893D0"},
        ]
    if color_mode == "hydropathy":
        return [
            {"label": "Hydrophilic (Kyte-Doolittle)", "color": "#6893D0"},
            {"label": "Hydrophobic (Kyte-Doolittle)", "color": "#F97316"},
        ]
    if color_mode == "conservation":
        return [
            {"label": "Lower MSA conservation", "color": "#1E293B"},
            {"label": "Higher MSA conservation", "color": "#34D399"},
        ]
    if color_mode == "charge":
        return [
            {"label": "Negative (D/E formal)", "color": "#F472B6"},
            {"label": "Neutral / His excluded", "color": "#94A3B8"},
            {"label": "Positive (K/R formal)", "color": "#6893D0"},
        ]
    items = []
    seen = set()
    for trace in backbone:
        chain = trace.get("chain_id")
        if chain in seen:
            continue
        seen.add(chain)
        items.append({"label": f"Chain {chain}", "color": trace.get("color")})
    if not items:
        for chain in chain_ids:
            dummy = {"chain_id": chain, "query_residue": "X"}
            items.append(
                {
                    "label": f"Chain {chain}",
                    "color": _color_for_residue(
                        dummy,
                        color_mode="chain",
                        highlights=set(),
                        selected_query_index=None,
                        overlay=None,
                        kind="",
                    ),
                }
            )
    return items


def _representation_note(representation: str) -> str:
    if representation == "backbone":
        return (
            "Backbone trace through deposited CA coordinates in polymer order. "
            "Not a DSSP/STRIDE cartoon."
        )
    if representation == "ca":
        return "CA atoms and CA-CA backbone from deposited coordinates."
    return "ATOM records as deposited. HETATM ligands are not drawn in this version."


def _finite_float(value: object) -> Optional[float]:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def _optional_int(value: object) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _hex_rgb(color: str) -> Tuple[int, int, int]:
    text = color.lstrip("#")
    if len(text) != 6:
        return (148, 163, 184)
    return int(text[0:2], 16), int(text[2:4], 16), int(text[4:6], 16)
