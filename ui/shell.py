"""Phase 19 Part 1 application shell: header, disclosures, workspace layout.

Structural containers only. Visual polish belongs to Part 2.
These helpers consume scientific result objects; they do not recompute them.
"""

from __future__ import annotations

from typing import Any, Callable, Optional, Sequence

import pandas as pd
import streamlit as st

from helixscope_core.compat import (
    genome_store,
    provenance,
    resource_admission,
    tool_registry,
    usalign,
)
from ui.command_palette import render_command_controls
from ui.components import html_escape, status_badge
from ui.navigation import (
    ModuleSpec,
    get_module,
    grouped_modules,
    render_back_control,
    set_active_module,
    workspace_object_catalog,
)
from ui.workspace import empty_structure_html

SHELL_CSS: str = """
<style>
/* Phase 19 Part 1 structural shell. Replaceable in Part 2. */
section[data-testid="stSidebar"] {
    display: flex !important;
    visibility: visible !important;
    transform: none !important;
    opacity: 1 !important;
    min-width: 220px !important;
    border: none !important;
    border-right: none !important;
    outline: none !important;
    box-shadow: none !important;
}
section[data-testid="stSidebar"] > div:first-child {
    width: 228px;
    margin-left: 0 !important;
    border-top-left-radius: 0 !important;
    border-bottom-left-radius: 0 !important;
}
[data-testid="stToolbar"] {
    display: flex !important;
    visibility: visible !important;
}
[data-testid="stSidebarCollapsedControl"],
[data-testid="collapsedControl"],
[data-testid="stExpandSidebarButton"],
[data-testid="stSidebarCollapseButton"] {
    display: none !important;
    visibility: hidden !important;
    pointer-events: none !important;
}
[data-testid="stAppViewContainer"] [data-testid="stMainBlockContainer"],
[data-testid="stMainBlockContainer"],
.block-container,
.stMainBlockContainer {
    max-width: 100% !important;
    width: 100% !important;
    padding-top: 1.1rem !important;
    padding-left: 1.25rem !important;
    padding-right: 1.25rem !important;
    padding-bottom: 2.5rem !important;
}
.helix-workspace {
    width: 100%;
    min-width: 0;
}
.helix-context-header {
    width: 100%;
    padding: 4px 0 10px 0;
    margin: 0 0 8px 0;
    border-bottom: none;
    box-shadow: inset 0 -1px 0 var(--hs-glass-highlight);
}
.helix-context-title {
    font-family: var(--hs-font-display), 'Space Grotesk', sans-serif;
    font-size: 1.12rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    text-transform: none;
    color: var(--hs-text);
    margin: 0;
}
.helix-context-summary {
    color: var(--hs-text-secondary);
    font-size: 13px;
    line-height: 1.45;
    margin: 6px 0 0 0;
}
.helix-nav-brand {
    font-family: var(--hs-font-display), 'Space Grotesk', sans-serif;
    font-size: 1.05rem;
    font-weight: 700;
    letter-spacing: 0.16em;
    text-transform: uppercase;
    color: var(--hs-text);
    margin: 0 0 2px 0;
}
.helix-nav-tag {
    color: var(--hs-text-muted);
    font-size: 11px;
    margin: 0 0 12px 0;
}
.helix-command-slot {
    min-height: 0;
}
.helix-viewport-slot {
    width: 100%;
    min-width: 0;
    min-height: 480px;
}
.helix-split {
    display: grid;
    grid-template-columns: minmax(0, 1fr) minmax(220px, 300px);
    gap: 1rem;
    width: 100%;
    min-width: 0;
}
.helix-three-col {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 1rem;
    width: 100%;
    min-width: 0;
}
.helix-inspector-slot {
    min-width: 0;
}
@media (max-width: 1365px) {
    .helix-split {
        grid-template-columns: 1fr;
    }
    .helix-three-col {
        grid-template-columns: 1fr;
    }
}
@media (min-width: 1366px) and (max-width: 1599px) {
    .helix-three-col {
        grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
    }
}
[data-testid="stAppViewContainer"] {
    overflow-x: hidden;
}
.helix-seq-viewer,
[data-testid="stDataFrame"] {
    overflow-x: auto;
}
</style>
"""


def inject_shell_css() -> None:
    """Inject structural layout CSS for the workstation shell.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.html(SHELL_CSS)
    st.markdown(SHELL_CSS, unsafe_allow_html=True)


def render_topbar(module_id: str) -> None:
    """Compact application bar: context path, version, command-palette slot.

    Args:
        module_id: Active semantic module.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    spec = get_module(module_id)
    left, version, nav, search = st.columns([5.2, 1.1, 2.8, 1.5], gap="small")
    with left:
        st.caption(f"{spec.group} / {spec.label}")
    with version:
        st.caption(str(provenance.HELIXSCOPE_VERSION))
    with nav:
        render_back_control()
    with search:
        render_command_controls()


def contextual_header(
    spec: ModuleSpec,
    *,
    status: Optional[str] = None,
) -> None:
    """Compact page header plus method/scope/privacy/provenance expanders.

    Args:
        spec: Module metadata from the navigation registry.
        status: Optional scientific status token to display unchanged.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    status_html = status_badge(status) if status else ""
    st.markdown(
        '<div class="helix-context-header helix-workspace">'
        f'<p class="helix-context-title">{html_escape(spec.title)}</p>'
        f'<p class="helix-context-summary">{html_escape(spec.summary)}</p>'
        f"{status_html}"
        "</div>",
        unsafe_allow_html=True,
    )
    render_scientific_disclosures(spec)


def render_scientific_disclosures(spec: ModuleSpec) -> None:
    """Reusable Scientific Method / Scope / Privacy / Provenance hierarchy.

    Args:
        spec: Module metadata. Empty strings are omitted.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    blocks: list[tuple[str, str, bool]] = []
    if spec.method:
        blocks.append(("Scientific method", spec.method, False))
    if spec.scope:
        blocks.append(("Scope and limitations", spec.scope, False))
    if spec.privacy:
        blocks.append(
            (
                "External services / privacy",
                spec.privacy,
                bool(spec.privacy_expanded),
            )
        )
    if spec.provenance:
        blocks.append(("Provenance", spec.provenance, False))
    if not blocks:
        return
    st.caption(
        "Scientific method · Scope and limitations · External services / privacy · Provenance"
    )
    cols = st.columns(len(blocks))
    for column, (title, body, expanded) in zip(cols, blocks):
        with column:
            with st.expander(title, expanded=expanded):
                st.markdown(body)


def workspace_container() -> None:
    """Mark the main scientific workspace as full width.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.markdown('<div class="helix-workspace"></div>', unsafe_allow_html=True)


def split_workspace() -> tuple[Any, Any]:
    """Main + inspector columns for Part 2 structure pages.

    Returns:
        Streamlit column pair (main, inspector).

    Raises:
        Nenhum.
    """
    return st.columns([3, 1], gap="medium")


def metric_row(n: int = 4) -> Sequence[Any]:
    """Equal columns for summary metrics.

    Args:
        n: Column count.

    Returns:
        Streamlit columns.

    Raises:
        Nenhum.
    """
    count = max(1, int(n))
    return st.columns(count)


def three_column_science() -> Sequence[Any]:
    """Foundation for TREE | MSA | STRUCTURE.

    Returns:
        Three Streamlit columns.

    Raises:
        Nenhum.
    """
    return st.columns(3, gap="medium")


def inspector_slot(sections: Optional[dict] = None) -> None:
    """Optional inspector region. Empty until a page supplies sections.

    Args:
        sections: Mapping passed to ui.components.inspector_panel, or None.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.markdown('<div class="helix-inspector-slot"></div>', unsafe_allow_html=True)
    if sections:
        from ui.components import inspector_panel

        st.markdown(inspector_panel(sections), unsafe_allow_html=True)


def _object_present(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, dict):
        return bool(value)
    if isinstance(value, (list, tuple)):
        return len(value) > 0
    return True


def render_analysis_history() -> None:
    """Lista as analises desta sessao. Nao recalcula e nao grava disco.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    from modules.session_history import (
        clear_analyses,
        delete_analysis,
        ensure_history,
        live_fingerprints,
        restore_analysis,
        sync_analysis_history,
    )

    sync_analysis_history(st.session_state)
    history = ensure_history(st.session_state)
    st.markdown("#### History")
    st.caption(
        "Analyses already completed in this browser session. "
        "Open restores the stored result. It does not launch a new engine run."
    )
    if history.get("refused"):
        st.caption(str(history["refused"]))
    entries = list(history.get("entries") or [])
    if not entries:
        st.caption("No completed analysis is stored yet.")
    for entry in reversed(entries):
        entry_id = str(entry.get("id") or "")
        when = str(entry.get("created_at") or "")
        clock = when[11:16] if len(when) >= 16 else when
        label_col, open_col, delete_col = st.columns([4.2, 1, 1], gap="small")
        with label_col:
            st.markdown(
                f"{html_escape(str(entry.get('module') or '').upper())} · "
                f"{html_escape(str(entry.get('analysis_type') or ''))} · "
                f"{html_escape(str(entry.get('summary') or ''))} · "
                f"{html_escape(clock)} · "
                f"{html_escape(str(entry.get('status') or ''))}"
            )
        with open_col:
            if st.button("Abrir", key=f"history_open_{entry_id}"):
                restore_analysis(history, entry_id, st.session_state)
                set_active_module(str(entry.get("module") or ""))
                st.rerun()
        with delete_col:
            if st.button("Excluir", key=f"history_delete_{entry_id}"):
                delete_analysis(history, entry_id)
                st.session_state["helix_analysis_history"] = history
                st.rerun()
    if st.button("Clear history", key="history_clear", disabled=not entries):
        clear_analyses(history, live_fingerprints(st.session_state))
        st.session_state["helix_analysis_history"] = history
        st.rerun()


def render_overview() -> None:
    """Useful workspace home: active objects, engines, module jump.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    render_analysis_history()
    st.markdown("#### Active scientific objects")
    st.caption(
        "This list is session state already stored by scientific modules. "
        "Opening Overview does not rerun BLAST, VEP, IQ-TREE, Cas-OFFinder, "
        "or the RCSB Alignment API."
    )
    rows = []
    for key, label, module_id in workspace_object_catalog():
        present = key in st.session_state and _object_present(st.session_state.get(key))
        status = "present" if present else "empty"
        extra = ""
        payload = st.session_state.get(key)
        if present and isinstance(payload, dict):
            extra = str(
                payload.get("status")
                or payload.get("kind")
                or payload.get("structure_id")
                or ""
            )
        rows.append(
            {
                "Object": label,
                "Session key": key,
                "State": status,
                "Status/kind": extra or "N/A",
                "Module": get_module(module_id).label,
            }
        )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    st.markdown("#### Open a module")
    for group, items in grouped_modules():
        st.caption(group)
        cols = st.columns(min(4, max(1, len(items))))
        for index, spec in enumerate(items):
            with cols[index % len(cols)]:
                if st.button(spec.label, key=f"overview_go_{spec.id}", width="stretch"):
                    set_active_module(spec.id)
                    st.rerun()

    st.markdown("#### Scientific engines (compact)")
    _render_engine_table()


def _engine_display_status(row: dict) -> str:
    live = row.get("live") if isinstance(row.get("live"), dict) else None
    available = bool(row.get("available"))
    remote = str(row.get("engine_location") or "") == "remote"
    if isinstance(live, dict) and live.get("ok") and (available or remote):
        return str(live.get("status") or row.get("status") or "")
    return str(row.get("status") or "not_installed")


def _backend_label(row: dict) -> str:
    if str(row.get("engine_location") or "") == "remote":
        return "REMOTE"
    if row.get("available"):
        return "LOCAL"
    return "UNAVAILABLE"


def engine_row_is_product_visible(row: dict) -> bool:
    """True only for a live or remote validation the product can stand behind.

    Args:
        row: One tool-registry record.

    Returns:
        False for NOT_INSTALLED, UNAVAILABLE, DETECTED and any other
        unvalidated state. Those names stay in the registry and stay off
        this screen.

    Raises:
        Nenhum.
    """
    status = _engine_display_status(row).strip().lower().replace("-", "_").replace(" ", "_")
    return status in {"live_validated", "remote_validated"}


def _render_engine_summary() -> None:
    snap = tool_registry.collect_tool_snapshot()
    tools = snap.get("tools") or {}
    for name, row in tools.items():
        if not isinstance(row, dict) or not engine_row_is_product_visible(row):
            continue
        status = _engine_display_status(row)
        version = row.get("version") or (row.get("live") or {}).get("version") or "not reported"
        backend = _backend_label(row)
        st.markdown(
            f'<div class="hs-engine-row hs-content">'
            f'<span class="hs-engine-name">{html_escape(name)}</span> '
            f"{status_badge(status)} "
            f'<span class="hs-caption">{html_escape(backend)} · version '
            f"{html_escape(version)}</span></div>",
            unsafe_allow_html=True,
        )


def _render_engine_table() -> None:
    _render_engine_summary()


def render_settings() -> None:
    """Full scientific-engine registry, previously on every page.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.caption(
        "Engines without a live or remote validation are omitted from this list. "
        "LIVE_VALIDATED requires a real accepted run. "
        "DETECTED is not LIVE_VALIDATED. Remote APIs are REMOTE_VALIDATED, never local live. "
        "FastTree, when present, is approximate maximum-likelihood and is not IQ-TREE."
    )
    st.markdown("#### Engine summary")
    _render_engine_summary()
    snap = tool_registry.collect_tool_snapshot()
    tools = snap.get("tools") or {}
    with st.expander("Validation details", expanded=False):
        st.caption(
            "Timestamp, fixture and version for engines that actually ran in this environment."
        )
        any_record = False
        for name, row in tools.items():
            if not isinstance(row, dict) or not engine_row_is_product_visible(row):
                continue
            live = row.get("live") if isinstance(row.get("live"), dict) else None
            if not isinstance(live, dict):
                continue
            any_record = True
            details = live.get("details") if isinstance(live.get("details"), dict) else {}
            st.markdown(
                f"{html_escape(name)} {status_badge(live.get('status'))}",
                unsafe_allow_html=True,
            )
            st.caption(
                f"version {live.get('version') or 'not reported'} · "
                f"validated {live.get('validated_at_utc') or 'unknown'} · "
                f"ok {live.get('ok')}"
            )
            if details:
                st.caption(
                    " · ".join(f"{key}={details[key]}" for key in list(details)[:8])
                )
        if not any_record:
            st.caption("No live/remote validation records are stored in this process.")
    with st.expander("Advanced detector table", expanded=False):
        rows = []
        for name, row in tools.items():
            if not isinstance(row, dict) or not engine_row_is_product_visible(row):
                continue
            live = row.get("live") if isinstance(row.get("live"), dict) else {}
            rows.append(
                {
                    "Tool": name,
                    "Detector status": row.get("status"),
                    "Validation": (live or {}).get("status") or "none",
                    "Version": row.get("version") or "not reported",
                    "Available": row.get("available"),
                    "Backend": _backend_label(row),
                    "Source": row.get("source"),
                    "Path": row.get("path") or "",
                    "Reason": row.get("reason") or "",
                }
            )
        if rows:
            st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)


def render_references() -> None:
    """Surface existing genome-store catalog without inventing a download flow.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.caption(
        "Cas-OFFinder TEST REFERENCE search remains in CRISPR Design. "
        "This page does not install GRCh38 and does not claim genome-wide scores."
    )
    rows = []
    for item in genome_store.list_local_references():
        rows.append(
            {
                "Assembly": item.get("id"),
                "Kind": item.get("kind"),
                "Local status": item.get("local_status"),
                "Ready": item.get("ready"),
                "Contigs": item.get("n_contigs"),
                "SHA-256": (str(item.get("file_sha256") or "")[:12] or "N/A"),
                "Not a public assembly": item.get("not_a_public_assembly"),
            }
        )
    if rows:
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        for item in genome_store.list_local_references():
            if item.get("not_a_public_assembly"):
                st.markdown(
                    f"{html_escape(item.get('id'))} {status_badge('TEST_ONLY')} "
                    f"<span class='hs-caption'>synthetic TEST REFERENCE; "
                    f"not a public genome. Local status "
                    f"{html_escape(item.get('local_status'))}.</span>",
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    f"{html_escape(item.get('id'))} "
                    f"{status_badge(item.get('local_status') or 'UNAVAILABLE')} "
                    f"<span class='hs-caption'>public catalog entry; "
                    f"READY is checksum-verified, not a silent fallback.</span>",
                    unsafe_allow_html=True,
                )
    memory = resource_admission.system_memory()
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        engine_available=True,
    )
    st.markdown("#### GRCh38.p14 resource admission")
    st.caption("Accession GCF_000001405.40. RESOURCE_LIMIT is not hidden.")
    st.markdown(
        f"Decision: {html_escape(decision.get('decision'))}. "
        f"Reason: {html_escape(decision.get('reason'))}. "
        f"Measured available bytes: {html_escape(memory.get('available_bytes'))}.",
        unsafe_allow_html=True,
    )


def render_structure_hub() -> None:
    """Wide 3D foundation: inventory plus viewport slot for Part 2.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.markdown('<div class="helix-viewport-slot helix-workspace hs-viewport"></div>', unsafe_allow_html=True)
    na_status = usalign.nucleic_structural_comparison_availability()
    st.caption(
        f"DNA/RNA structural comparison: {na_status.get('status')}: "
        f"{na_status.get('reason')} "
        "Protein 3D remains in Protein Analysis. Superposition remains in Compare. "
        "This hub lists session structures; it does not fetch PDB files."
    )
    catalog = (
        ("dna_3d_representation", "DNA 3D", "dna"),
        ("dna_struct_result", "DNA deposited structure", "dna"),
        ("rna_3d_representation", "RNA 3D", "rna"),
        ("prot_struct_result", "Protein structure", "protein"),
        ("crispr_complex", "CRISPR complex", "crispr"),
        ("molecule_selection", "Current selection", "compare"),
    )
    rows = []
    for key, label, module_id in catalog:
        payload = st.session_state.get(key)
        present = _object_present(payload)
        kind = ""
        status = ""
        structure_id = ""
        if isinstance(payload, dict):
            kind = str(payload.get("kind") or "")
            status = str(payload.get("status") or "")
            structure_id = str(payload.get("structure_id") or "")
        rows.append(
            {
                "Object": label,
                "Present": present,
                "Kind": kind or "N/A",
                "Status": status or "N/A",
                "Structure ID": structure_id or "N/A",
                "Open in": get_module(module_id).label,
            }
        )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    main, inspector = split_workspace()
    with main:
        loaded = any(row.get("Present") for row in rows)
        if loaded:
            st.caption(
                "Open DNA, RNA, Protein or Compare to inspect the loaded "
                "coordinates. This hub does not remount those Plotly keys."
            )
        else:
            st.markdown(empty_structure_html(), unsafe_allow_html=True)
    with inspector:
        inspector_slot()
        st.caption("Residue inspector is populated inside DNA, RNA, Protein and Compare when a mapped residue is selected.")
    jump_cols = st.columns(4)
    for column, module_id in zip(jump_cols, ("dna", "rna", "protein", "compare")):
        spec = get_module(module_id)
        with column:
            if st.button(spec.label, key=f"hub_go_{module_id}", width="stretch"):
                set_active_module(module_id)
                st.rerun()


def render_module_page(module_id: str, renderer: Callable[[], None]) -> None:
    """Header + disclosures + module renderer.

    Args:
        module_id: Semantic id.
        renderer: Existing scientific UI function. Must not be a recalculation.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    workspace_container()
    contextual_header(get_module(module_id))
    renderer()
    from modules.session_history import sync_analysis_history

    sync_analysis_history(st.session_state)
