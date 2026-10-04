"""Plotly 3D viewer for a protein SceneStructure.

Builds figures only from a scene model. Does not fetch PDB/AlphaFold/NCBI,
does not parse mmCIF, and does not invent coordinates.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import plotly.graph_objects as go

from ui.charts import _apply_base_layout
from ui.tokens import (
    FONT_SANS,
    HS_ICE_BLUE,
    HS_MID_BLUE,
    HS_TEXT_SECONDARY,
    HS_WHITE_BLUE,
    PLOTLY_3D_BG,
    STRUCTURE_VIEWPORT_HEIGHT,
)

STRUCTURE_VIEWER_API: int = 14
"""Sentinel de contrato. Bump quando funcoes publicas de selecao/cena mudarem."""


def plotly_safe_text(text: object) -> str:
    """Remove marcacao HTML de texto enviado ao hover/titulo Plotly."""
    return str("" if text is None else text).replace("<", "").replace(">", "")


def count_finite_xyz(x: Sequence[Any], y: Sequence[Any], z: Sequence[Any]) -> int:
    """Conta pontos finitos. None e corte de linha, nao coordenada.

    Args:
        x: Coordenadas X (None permitido como quebra).
        y: Coordenadas Y.
        z: Coordenadas Z.

    Returns:
        Numero de ternos finitos.

    Raises:
        ValueError: Insufficient data se comprimentos diferirem.
    """
    if len(x) != len(y) or len(y) != len(z):
        raise ValueError("Insufficient data")
    n = 0
    for a, b, c in zip(x, y, z):
        if a is None and b is None and c is None:
            continue
        try:
            fa = float(a)
            fb = float(b)
            fc = float(c)
        except (TypeError, ValueError) as exc:
            raise ValueError("Insufficient data") from exc
        if fa != fa or fb != fb or fc != fc:
            raise ValueError("Insufficient data")
        if abs(fa) == float("inf") or abs(fb) == float("inf") or abs(fc) == float("inf"):
            raise ValueError("Insufficient data")
        n += 1
    return n


def figure_from_scene(scene: Mapping[str, Any], *, center_selection: bool = False) -> go.Figure:
    """Constroi uma figura Plotly 3D a partir do scene model.

    Args:
        scene: Saida de structure_scene.build_scene.
        center_selection: Se verdadeiro, usa o centro de camera da selecao.

    Returns:
        Figura Scatter3d. Um trace por segmento de backbone; atomos num unico
        scatter (nao uma mesh por atomo). Pares nucleicos numa unica trace.

    Raises:
        ValueError: Insufficient data se o scene nao estiver READY ou se
            arrays XYZ estiverem vazios, dessincronizados ou nao finitos.

    Nota:
        Cores, brilho e opacidade sao visuais. Hover nao dispara fetch.
        Identidades de base/residuo no hover vem do mapping, nao de geometria
        inventada.
    """
    if str(scene.get("status") or "") != "READY":
        raise ValueError("Insufficient data")
    objects = scene.get("objects") or {}
    traces = list(objects.get("backbone_traces") or [])
    atoms = objects.get("atom_points")
    if not traces and not atoms:
        raise ValueError("Insufficient data")

    finite_points = 0
    fig = go.Figure()
    for trace in traces:
        xs = list(trace.get("x") or [])
        ys = list(trace.get("y") or [])
        zs = list(trace.get("z") or [])
        finite_points += count_finite_xyz(xs, ys, zs)
        marker_colors = list(trace.get("colors") or [])
        if not marker_colors:
            marker_colors = [trace.get("color") or HS_MID_BLUE] * len(xs)
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="lines+markers",
                line=dict(
                    color=trace.get("color") or HS_MID_BLUE,
                    width=int((scene.get("materials") or {}).get("line_width") or 5),
                ),
                marker=dict(
                    size=int((scene.get("materials") or {}).get("marker_size") or 5),
                    color=marker_colors,
                    opacity=0.98,
                    line=dict(width=0.4, color="#020617"),
                ),
                customdata=list(trace.get("customdata") or []),
                text=[plotly_safe_text(item) for item in list(trace.get("hovertext") or [])],
                hovertemplate="%{text}<extra></extra>",
                name=plotly_safe_text(trace.get("name") or "backbone"),
                showlegend=True,
            )
        )
    pair_traces = list(objects.get("pair_traces") or [])
    for index, trace in enumerate(pair_traces):
        xs = list(trace.get("x") or [])
        ys = list(trace.get("y") or [])
        zs = list(trace.get("z") or [])
        finite_points += count_finite_xyz(xs, ys, zs)
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="lines",
                line=dict(color=trace.get("color") or "#94A3B8", width=2),
                text=[plotly_safe_text(item) for item in list(trace.get("hovertext") or [])],
                hovertemplate="%{text}<extra></extra>",
                name="base pairs" if index == 0 else "pair",
                showlegend=index == 0,
                legendgroup="pairs",
            )
        )
    if atoms:
        xs = list(atoms.get("x") or [])
        ys = list(atoms.get("y") or [])
        zs = list(atoms.get("z") or [])
        finite_points += count_finite_xyz(xs, ys, zs)
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="markers",
                marker=dict(
                    size=list(atoms.get("sizes") or []),
                    color=list(atoms.get("colors") or []),
                    opacity=float((scene.get("materials") or {}).get("opacity") or 0.95),
                    line=dict(width=0),
                ),
                customdata=list(atoms.get("customdata") or []),
                text=[plotly_safe_text(item) for item in list(atoms.get("hovertext") or [])],
                hovertemplate="%{text}<extra></extra>",
                name=plotly_safe_text(atoms.get("name") or "atoms"),
                showlegend=False,
            )
        )
    marker = objects.get("selection_marker")
    if marker:
        xs = list(marker.get("x") or [])
        ys = list(marker.get("y") or [])
        zs = list(marker.get("z") or [])
        if xs:
            finite_points += count_finite_xyz(xs, ys, zs)
            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    mode="markers",
                    marker=dict(size=11, color=HS_WHITE_BLUE, opacity=1.0, line=dict(width=1, color=HS_ICE_BLUE)),
                    customdata=list(marker.get("customdata") or []),
                    text=[plotly_safe_text(item) for item in list(marker.get("hovertext") or [])],
                    hovertemplate="%{text}<extra></extra>",
                    name="selection",
                    showlegend=False,
                )
            )
    for terminus in list(objects.get("terminus_markers") or []):
        xs = list(terminus.get("x") or [])
        ys = list(terminus.get("y") or [])
        zs = list(terminus.get("z") or [])
        if not xs:
            continue
        finite_points += count_finite_xyz(xs, ys, zs)
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode="markers",
                marker=dict(size=9, color="#FBBF24", opacity=1.0, symbol="diamond", line=dict(width=1, color="#020617")),
                text=[plotly_safe_text(item) for item in list(terminus.get("hovertext") or [])],
                hovertemplate="%{text}<extra></extra>",
                name=plotly_safe_text(terminus.get("name") or "termini"),
                showlegend=True,
            )
        )
    if finite_points <= 0 or len(fig.data) <= 0:
        raise ValueError("Insufficient data")
    camera = dict(scene.get("camera") or {})
    eye = camera.get("eye") or {"x": 1.75, "y": 1.55, "z": 1.35}
    center = camera.get("center") if center_selection else {"x": 0.0, "y": 0.0, "z": 0.0}
    fig = _apply_base_layout(fig, cartesian=False)
    fig.update_layout(
        scene=dict(
            aspectmode="data",
            xaxis=_hidden_axis("X (A)"),
            yaxis=_hidden_axis("Y (A)"),
            zaxis=_hidden_axis("Z (A)"),
            camera=dict(
                eye=eye,
                up=camera.get("up") or {"x": 0.0, "y": 0.0, "z": 1.0},
                center=center or {"x": 0.0, "y": 0.0, "z": 0.0},
            ),
            bgcolor=PLOTLY_3D_BG,
            hovermode="closest",
        ),
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.14,
            x=0,
            font=dict(size=11, color=HS_TEXT_SECONDARY),
            bgcolor="rgba(10, 10, 10, 0.55)",
            itemwidth=40,
        ),
        autosize=False,
        height=STRUCTURE_VIEWPORT_HEIGHT,
        paper_bgcolor=PLOTLY_3D_BG,
        plot_bgcolor=PLOTLY_3D_BG,
        margin=dict(l=8, r=8, t=48, b=72),
        uirevision="helixscope-structure-3d" if not center_selection else None,
        hoverlabel=dict(
            bgcolor="rgba(10,10,10,0.94)",
            bordercolor="rgba(255,255,255,0.20)",
            font=dict(color="#FFFFFF", size=12),
            align="left",
        ),
        title=dict(
            text=plotly_safe_text(_title(scene)),
            font=dict(size=12, color=HS_TEXT_SECONDARY),
            y=0.98,
            yanchor="top",
        ),
    )
    fig.update_scenes(dragmode="orbit")
    return fig


def figure_from_superposition(
    bundle: Mapping[str, Any],
    *,
    max_points_per_structure: int = 4000,
    view_mode: str = "overlay",
) -> go.Figure:
    """Figura de referencia + alvo transformado. Nao altera coordenadas origem.

    Args:
        bundle: Saida de structure_superposition.superposition_bundle.
        max_points_per_structure: Teto de pontos CA por estrutura (RESOURCE_LIMIT).
        view_mode: overlay (ambos), reference/A, ou target/B. Apenas visibilidade.

    Returns:
        Scatter3d com traces visiveis. Coordenadas origem nao sao modificadas.

    Raises:
        ValueError: Insufficient data se nao houver CA finitos.

    Nota:
        Overlay visual nao e o alinhamento. RMSD/TM-score vivem no objeto
        cientifico. Alinhamento flexivel multi-bloco chega com visual PARTIAL.
    """
    ref_atoms = list(bundle.get("reference_atoms") or [])
    tgt_atoms = list(bundle.get("transformed_target_atoms") or [])
    ref_xyz, ref_text, ref_custom = _ca_trace_arrays(
        ref_atoms, side="reference", cap=max_points_per_structure
    )
    tgt_xyz, tgt_text, tgt_custom = _ca_trace_arrays(
        tgt_atoms, side="target", cap=max_points_per_structure
    )
    if not ref_xyz[0] and not tgt_xyz[0]:
        raise ValueError("Insufficient data")
    finite = 0
    if ref_xyz[0]:
        finite += count_finite_xyz(*ref_xyz)
    if tgt_xyz[0]:
        finite += count_finite_xyz(*tgt_xyz)
    if finite < 2:
        raise ValueError("Insufficient data")
    mode = str(view_mode or "overlay").strip().lower()
    show_ref = mode in {"overlay", "reference", "a"}
    show_tgt = mode in {"overlay", "target", "b"}
    fig = go.Figure()
    if ref_xyz[0] and show_ref:
        fig.add_trace(
            go.Scatter3d(
                x=ref_xyz[0],
                y=ref_xyz[1],
                z=ref_xyz[2],
                mode="markers",
                name=plotly_safe_text(
                    f"reference {bundle.get('reference_entry') or ''} "
                    f"{bundle.get('reference_kind') or ''}"
                ),
                marker=dict(size=3, color=HS_MID_BLUE),
                text=ref_text,
                customdata=ref_custom,
                hovertemplate="%{text}<extra></extra>",
            )
        )
    if tgt_xyz[0] and show_tgt:
        fig.add_trace(
            go.Scatter3d(
                x=tgt_xyz[0],
                y=tgt_xyz[1],
                z=tgt_xyz[2],
                mode="markers",
                name=plotly_safe_text(
                    f"target (transformed copy) {bundle.get('target_entry') or ''} "
                    f"{bundle.get('target_kind') or ''}"
                ),
                marker=dict(size=3, color="#F472B6"),
                text=tgt_text,
                customdata=tgt_custom,
                hovertemplate="%{text}<extra></extra>",
            )
        )
    if not fig.data:
        raise ValueError("Insufficient data")
    visual = str(bundle.get("visual_status") or "")
    title = (
        f"{bundle.get('method') or 'alignment'} | "
        f"visual {visual} | "
        f"CA pairs {bundle.get('n_aligned_residue_pairs')} | "
        "source coordinates unchanged"
    )
    fig.update_layout(
        paper_bgcolor=PLOTLY_3D_BG,
        plot_bgcolor=PLOTLY_3D_BG,
        font=dict(color="#FFFFFF", family=FONT_SANS),
        margin=dict(l=8, r=8, t=36, b=8),
        height=STRUCTURE_VIEWPORT_HEIGHT,
        legend=dict(font=dict(size=11, color=HS_TEXT_SECONDARY), bgcolor="rgba(10,10,10,0.6)"),
        scene=dict(
            xaxis=_hidden_axis("X (A)"),
            yaxis=_hidden_axis("Y (A)"),
            zaxis=_hidden_axis("Z (A)"),
            aspectmode="data",
            bgcolor=PLOTLY_3D_BG,
        ),
        hoverlabel=dict(
            bgcolor="rgba(10,10,10,0.94)",
            bordercolor="rgba(255,255,255,0.20)",
            font=dict(color="#FFFFFF", size=12),
            align="left",
        ),
        title=dict(
            text=plotly_safe_text(title),
            font=dict(size=12, color=HS_TEXT_SECONDARY),
            y=0.98,
            yanchor="top",
        ),
    )
    fig.update_scenes(dragmode="orbit")
    return fig


def _ca_trace_arrays(
    atoms: Sequence[Mapping[str, Any]],
    *,
    side: str,
    cap: int,
) -> tuple[tuple[list, list, list], list[str], list[list]]:
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    texts: list[str] = []
    custom: list[list] = []
    for atom in atoms:
        if not isinstance(atom, Mapping):
            continue
        if str(atom.get("atom_name") or "").strip().upper() != "CA":
            continue
        if str(atom.get("group") or "") != "ATOM":
            continue
        try:
            x = float(atom["x"])
            y = float(atom["y"])
            z = float(atom["z"])
        except (TypeError, ValueError, KeyError):
            continue
        if not all(map(lambda v: v == v and abs(v) != float("inf"), (x, y, z))):
            continue
        chain = str(atom.get("label_asym_id") or atom.get("auth_asym_id") or "")
        seq = atom.get("label_seq_id")
        xs.append(x)
        ys.append(y)
        zs.append(z)
        texts.append(
            plotly_safe_text(
                f"{side} {chain}:{seq} {atom.get('comp_id') or ''} "
                f"{'transformed-copy' if atom.get('transformed') else 'original-copy'}"
            )
        )
        custom.append([side, chain, seq])
        if len(xs) >= int(cap):
            break
    return (xs, ys, zs), texts, custom


def parse_superposition_selection(event: object) -> Optional[dict]:
    """Le customdata [side, chain, seq] do overlay de superposicao.

    Args:
        event: PlotlyState, dict, ou None.

    Returns:
        Dict side, asym_id, label_seq_id, ou None.

    Raises:
        Nenhum.
    """
    if event is None:
        return None
    selection = getattr(event, "selection", None)
    if selection is None and isinstance(event, Mapping):
        selection = event.get("selection")
    if selection is None:
        return None
    points = getattr(selection, "points", None)
    if points is None and isinstance(selection, Mapping):
        points = selection.get("points")
    if not points:
        return None
    point = points[0]
    payload = point.get("customdata") if isinstance(point, Mapping) else getattr(point, "customdata", None)
    if isinstance(payload, (list, tuple)) and len(payload) >= 3:
        try:
            seq = int(payload[2]) if payload[2] is not None else None
        except (TypeError, ValueError):
            seq = None
        return {
            "side": str(payload[0] or ""),
            "asym_id": str(payload[1] or ""),
            "label_seq_id": seq,
        }
    return None


def parse_plotly_selection(event: object) -> Optional[dict]:
    """Extrai customdata de um evento Streamlit Plotly. Nao faz fetch.

    Args:
        event: PlotlyState, dict, ou None.

    Returns:
        Payload de selecao ou None.
    """
    if event is None:
        return None
    selection = getattr(event, "selection", None)
    if selection is None and isinstance(event, Mapping):
        selection = event.get("selection")
    if selection is None:
        return None
    points = getattr(selection, "points", None)
    if points is None and isinstance(selection, Mapping):
        points = selection.get("points")
    if not points:
        return None
    point = points[0]
    if not isinstance(point, Mapping):
        custom = getattr(point, "customdata", None)
        payload = custom
    else:
        payload = point.get("customdata")
    from helixscope_core.compat import structure_scene

    return structure_scene.parse_selection_payload(payload)


def parse_plotly_selection_indices(event: object) -> list[int]:
    """Indices de query 0-based de todos os pontos Plotly mapeados.

    Args:
        event: PlotlyState, dict, ou None.

    Returns:
        Lista de indices unicos, ordenados. Vazia se nao houver mapping.
    """
    if event is None:
        return []
    selection = getattr(event, "selection", None)
    if selection is None and isinstance(event, Mapping):
        selection = event.get("selection")
    if selection is None:
        return []
    points = getattr(selection, "points", None)
    if points is None and isinstance(selection, Mapping):
        points = selection.get("points")
    if not points:
        return []
    from helixscope_core.compat import structure_scene

    indices: list[int] = []
    seen: set[int] = set()
    for point in points:
        if not isinstance(point, Mapping):
            payload = getattr(point, "customdata", None)
        else:
            payload = point.get("customdata")
        parsed = structure_scene.parse_selection_payload(payload)
        if not parsed:
            continue
        idx = parsed.get("query_index_0based")
        if idx is None:
            continue
        index = int(idx)
        if index in seen:
            continue
        seen.add(index)
        indices.append(index)
    indices.sort()
    return indices


def _hidden_axis(title: str) -> dict:
    return dict(
        title=title,
        showbackground=False,
        showgrid=True,
        gridcolor="rgba(255,255,255,0.04)",
        zeroline=False,
        showspikes=False,
        backgroundcolor="rgba(0,0,0,0)",
        color=HS_TEXT_SECONDARY,
    )


def _title(scene: Mapping[str, Any]) -> str:
    structure_id = scene.get("structure_id") or "structure"
    kind = scene.get("kind_label") or scene.get("kind") or ""
    scope = scene.get("view_label") or ""
    return f"{structure_id} — {kind} — {scope}"
