"""Reusable HelixScope workspace surfaces.

These helpers format scientific result objects for display. They do not
recompute GC, RMSD, alignments, or variant consequences.
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import streamlit as st

from ui.components import html_escape, status_badge
from ui.tokens import canonical_status


def unavailable_text(value: object) -> str:
    """Render missing scientific fields as Unavailable, never as 0.

    Args:
        value: Field from a result object.

    Returns:
        Display string.

    Raises:
        Nenhum.
    """
    if value is None or value == "":
        return "Unavailable"
    text = str(value)
    if text.strip().upper() in {"N/A", "NONE", "NULL"}:
        return "Unavailable"
    return text


def metric_display_value(value: object) -> str:
    """Format a metric cell. Empty becomes Unavailable; N/A stays N/A, never 0.

    Args:
        value: Field from a result object.

    Returns:
        Display string.

    Raises:
        Nenhum.
    """
    if value is None or value == "":
        return "Unavailable"
    text = str(value)
    if text.strip().upper() in {"NONE", "NULL"}:
        return "Unavailable"
    return text


def result_header_html(
    title: str,
    *,
    status: object = "",
    lines: Sequence[str] = (),
) -> str:
    """Compact post-analysis header. Values must already be computed.

    Args:
        title: Module or entity title.
        status: Scientific status token.
        lines: Preformatted metric lines from the scientific layer.

    Returns:
        HTML string.

    Raises:
        Nenhum.
    """
    status_html = status_badge(str(status)) if status else ""
    body = "".join(
        f'<div class="hs-result-line">{html_escape(line)}</div>' for line in lines if line
    )
    return (
        f'<div class="hs-result-header hs-content">'
        f'<div class="hs-result-title-row">'
        f'<h2 class="hs-result-title">{html_escape(title)}</h2>{status_html}'
        f"</div>{body}</div>"
    )


def metric_grid_html(
    items: Sequence[tuple[str, str, str]],
    *,
    columns: int | None = None,
) -> str:
    """Balanced scientific metric grid. Missing values stay Unavailable, never 0.

    Args:
        items: Sequence of (label, value, unit).
        columns: Explicit column count. Default: 4 when n>=4, else n (min 1).
            Five-plus items wrap; leftover cells do not stretch to full width.

    Returns:
        HTML string.

    Raises:
        Nenhum.
    """
    cells: list[str] = []
    for label, value, unit in items:
        shown = metric_display_value(value)
        unit_html = ""
        if unit and shown not in {"Unavailable", "N/A", "None"}:
            unit_html = f'<span class="hs-metric-unit">{html_escape(unit)}</span>'
        cells.append(
            '<div class="hs-metric">'
            f'<div class="hs-metric-value"><span class="hs-metric-number">'
            f"{html_escape(shown)}</span>{unit_html}</div>"
            f'<div class="hs-metric-label">{html_escape(label)}</div>'
            "</div>"
        )
    n = len(cells)
    if columns is None:
        if n <= 0:
            cols = 1
        elif n <= 3:
            cols = n
        else:
            cols = 4
    else:
        cols = max(1, min(6, int(columns)))
    inner = "".join(cells)
    return (
        f'<div class="hs-metric-grid hs-content" data-cols="{cols}" '
        f'style="--hs-metric-cols:{cols};">{inner}</div>'
    )


def render_metric_grid(
    items: Sequence[tuple[str, str, str]],
    *,
    columns: int | None = None,
) -> None:
    """Streamlit wrapper for metric_grid_html.

    Args:
        items: Sequence of (label, value, unit).
        columns: Optional column count.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    if not items:
        return
    st.markdown(metric_grid_html(items, columns=columns), unsafe_allow_html=True)


def metric_strip_html(items: Sequence[tuple[str, str, str]]) -> str:
    """Compact scientific metric strip (not finance KPIs).

    Args:
        items: Sequence of (label, value, unit). Empty values become Unavailable.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    return metric_grid_html(items)


def compare_hero_html(alignment: Mapping[str, object]) -> str:
    """Structure Compare header from an alignment record. Does not merge RMSD.

    Args:
        alignment: structure_alignment comparison object.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    ref = alignment.get("reference") if isinstance(alignment.get("reference"), Mapping) else {}
    tgt = alignment.get("target") if isinstance(alignment.get("target"), Mapping) else {}
    rmsd_g = alignment.get("rmsd_global_angstrom")
    rmsd_b = alignment.get("rmsd_block0_angstrom")
    n_pairs = alignment.get("n_aligned_residue_pairs")
    coverage = alignment.get("aln_coverage_percent") or []
    tm_scores = alignment.get("tm_scores") or []
    tm_text = "Unavailable"
    if isinstance(tm_scores, Sequence) and tm_scores:
        first = tm_scores[0] if isinstance(tm_scores[0], Mapping) else {}
        if first.get("value") is not None:
            tm_text = str(first.get("value"))
    cov_text = "Unavailable"
    if isinstance(coverage, Sequence) and len(coverage) >= 2:
        cov_text = f"{coverage[0]}% / {coverage[1]}%"
    method = str(alignment.get("method_display") or alignment.get("method") or "")
    engine = str(alignment.get("engine") or "RCSB")
    ref_id = html_escape(ref.get("entry_id"))
    tgt_id = html_escape(tgt.get("entry_id"))
    metrics = metric_strip_html(
        (
            ("Block RMSD", unavailable_text(None if rmsd_b is None else f"{rmsd_b}"), "A"),
            ("Global RMSD", unavailable_text(None if rmsd_g is None else f"{rmsd_g}"), "A"),
            ("TM-score", tm_text, ""),
            ("Coverage", cov_text, ""),
            ("Aligned residue pairs", unavailable_text(None if n_pairs is None else str(n_pairs)), ""),
        )
    )
    return (
        '<div class="hs-compare-hero hs-content">'
        '<div class="hs-compare-pair">'
        f'<div class="hs-compare-side"><div class="hs-mono">{ref_id}</div>'
        f'{status_badge(str(ref.get("kind_label") or "EXPERIMENTAL"))}</div>'
        f'<div class="hs-compare-method">{html_escape(method)}<div class="hs-caption">'
        f"{html_escape(engine)} · {html_escape(str(alignment.get('status') or 'RETRIEVED'))}"
        "</div></div>"
        f'<div class="hs-compare-side"><div class="hs-mono">{tgt_id}</div>'
        f'{status_badge(str(tgt.get("kind_label") or "EXPERIMENTAL"))}</div>'
        f"</div>{metrics}</div>"
    )


def variant_journey_html(steps: Sequence[tuple[str, str, str]]) -> str:
    """Genome-to-structure journey. Missing links stay UNMAPPED/UNAVAILABLE.

    Args:
        steps: (label, value, status). Status UNAVAILABLE/UNMAPPED breaks the chain visually.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    parts = []
    for index, (label, value, status) in enumerate(steps):
        broken = canonical_status(status) in {"UNMAPPED", "UNAVAILABLE", "ERROR"}
        klass = "hs-journey-step hs-journey-break" if broken else "hs-journey-step"
        parts.append(
            f'<div class="{klass}"><div class="hs-journey-label">{html_escape(label)}</div>'
            f'<div class="hs-mono">{html_escape(unavailable_text(value))}</div>'
            f"{status_badge(status)}</div>"
        )
        if index < len(steps) - 1:
            next_status = canonical_status(steps[index + 1][2])
            next_broken = next_status in {"UNMAPPED", "UNAVAILABLE", "ERROR"}
            connector = "hs-journey-gap" if broken or next_broken else "hs-journey-arrow"
            parts.append(f'<div class="{connector}" aria-hidden="true"></div>')
    return f'<div class="hs-journey">{"".join(parts)}</div>'


def evidence_item_html(item: Mapping[str, object]) -> str:
    """One evidence row. Does not invent confidence.

    Args:
        item: evidence_workspace item.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    field = html_escape(item.get("field"))
    value = html_escape(unavailable_text(item.get("value")))
    source = html_escape(item.get("source"))
    ident = html_escape(item.get("identifier"))
    status = str(item.get("evidence_status") or "UNAVAILABLE")
    stamp = html_escape(item.get("retrieved_at_utc") or item.get("computed_at_utc") or "")
    mapping = html_escape(item.get("mapping_status") or "")
    return (
        f'<div class="hs-evidence-item hs-content">'
        f'<div class="hs-evidence-field">{field}</div>'
        f'<div class="hs-mono">{value}</div>'
        f'<div class="hs-caption">{source} · {ident}</div>'
        f'<div class="hs-evidence-meta">{status_badge(status)}'
        f'<span class="hs-caption">{stamp} {mapping}</span></div>'
        f"</div>"
    )


CHAIN_ROLE_LABELS: dict[str, str] = {
    "cas_protein": "Cas protein",
    "guide_rna": "Guide RNA",
    "target_dna": "Target DNA",
    "nontarget_dna": "Nontarget DNA",
    "dna": "DNA",
    "rna": "RNA",
    "protein": "Protein",
    "unknown": "Unknown",
}


def chain_role_legend_html(chains: Sequence[Mapping[str, object]]) -> str:
    """Legend from already-classified chain roles. Does not invent roles.

    Args:
        chains: Complex annotation chain dicts with role and chain_id.

    Returns:
        HTML, or empty string when there are no chains.

    Raises:
        Nenhum.
    """
    if not chains:
        return ""
    cells = []
    for item in chains:
        role = str(item.get("role") or "unknown").strip().lower() or "unknown"
        label = CHAIN_ROLE_LABELS.get(role, role.replace("_", " "))
        chain_id = unavailable_text(item.get("chain_id"))
        cells.append(
            f'<div class="hs-legend-item">'
            f'<span class="hs-legend-swatch hs-legend-{html_escape(role)}"></span>'
            f'<span class="hs-mono">{html_escape(chain_id)}</span>'
            f'<span class="hs-caption">{html_escape(label)}</span>'
            f"</div>"
        )
    return f'<div class="hs-legend hs-content">{"".join(cells)}</div>'


def job_state_html(status: object, detail: str = "") -> str:
    """BLAST/MSA job chrome. Status text is primary; color is secondary.

    Args:
        status: Job status token already stored on the result.
        detail: Provenance line (RID, tool, etc.).

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    extra = (
        f'<span class="hs-caption">{html_escape(detail)}</span>' if detail else ""
    )
    return (
        f'<div class="hs-job-state hs-content">{status_badge(status)}{extra}</div>'
    )


def empty_structure_html() -> str:
    """Honest empty viewport. No fake molecule.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    return (
        '<div class="hs-empty-hero hs-content">'
        "<p>No structure loaded</p>"
        "<p class=\"hs-caption\">Load an experimental PDB or predicted AlphaFold "
        "structure to begin structural inspection.</p>"
        "</div>"
    )


def render_explanation(kind: str, values: Mapping[str, object] | None = None) -> None:
    """Two-layer scientific interpretation. Does not recompute results.

    Args:
        kind: Explanation kind (dna, compare, variant, ...).
        values: Already-computed result fields.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    from helixscope_core.explain import explain as build_explanation

    try:
        record = build_explanation(kind, values or {})
    except ValueError:
        return
    st.markdown(
        (
            f'<div class="hs-explain hs-content">'
            f'<div class="hs-explain-kicker">What this means</div>'
            f'<p class="hs-explain-body">{html_escape(record.plain_meaning)}</p>'
            f"</div>"
        ),
        unsafe_allow_html=True,
    )
    with st.expander("How to interpret this result"):
        st.caption(record.interpretation)
        st.caption(f"Why this result: {record.why_this_result}")
        st.caption(f"Important: {record.limitations}")
        st.caption(f"Method: {record.method}. Source: {record.source}. Status: {record.status}.")


def render_copy_caption(label: str, value: object) -> None:
    """Compact copyable scientific identifier.

    Args:
        label: Field name.
        value: Identifier (PDB, accession, hash).

    Returns:
        None.

    Raises:
        Nenhum.
    """
    text = unavailable_text(value)
    st.caption(f"{label}: `{text}`")
