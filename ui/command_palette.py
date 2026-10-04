"""HelixScope command palette: real module navigation.

Opens from the topbar Pesquisar button or Ctrl/Cmd+K. Filters the navigation
registry. Does not search remote databases and does not call an LLM.
"""

from __future__ import annotations

import streamlit as st

from ui.navigation import grouped_modules, module_ids, set_active_module
from ui.pointer import mount_command_hotkey


def _matches(query: str) -> list[tuple[str, str, str]]:
    needle = query.strip().lower()
    rows: list[tuple[str, str, str]] = []
    for group, items in grouped_modules():
        for spec in items:
            blob = f"{spec.label} {spec.group} {spec.title} {spec.id}".lower()
            if not needle or needle in blob:
                rows.append((spec.id, spec.label, group))
    return rows


@st.dialog("Pesquisar")
def _command_dialog() -> None:
    """Modal: type to filter, click to navigate.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    st.caption("Navigate HelixScope modules. This is not scientific search.")
    query = st.text_input("Filter", key="helix_cmd_query", placeholder="DNA, Variant, Compare")
    rows = _matches(query)
    if not rows:
        st.caption("No matching module.")
        return
    for module_id, label, group in rows:
        if st.button(f"{group} / {label}", key=f"helix_cmd_go_{module_id}", width="stretch"):
            set_active_module(module_id)
            st.rerun()


def render_command_controls() -> None:
    """Topbar trigger plus optional keyboard open.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    opened = False
    if st.button("Pesquisar", key="helix_cmd_open", help="Open module search palette"):
        opened = True
    result = mount_command_hotkey()
    trigger = getattr(result, "open", None) if result is not None else None
    if trigger:
        opened = True
    if opened:
        _command_dialog()


def registered_command_targets() -> tuple[str, ...]:
    """Return navigable module ids for tests.

    Returns:
        Registry ids in sidebar order.

    Raises:
        Nenhum.
    """
    return module_ids()
