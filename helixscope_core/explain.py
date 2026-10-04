"""Deterministic explainability. Not an LLM. Belongs in the scientific core."""

from __future__ import annotations

from modules.explain import (
    DIAGNOSIS_MARKERS,
    EXPERIMENTAL_CLAIM_MARKERS,
    Explanation,
    display_number,
    display_text,
    explain,
    explain_metrics as explain_metrics_records,
    explanation_as_dict,
)

__all__ = (
    "DIAGNOSIS_MARKERS",
    "EXPERIMENTAL_CLAIM_MARKERS",
    "Explanation",
    "display_number",
    "display_text",
    "explain",
    "explain_metrics",
    "explain_result",
    "explanation_as_dict",
)


def explain_result(kind: str, values: dict | None = None) -> dict:
    """Structured explanation dict for API/UI. Same text as Streamlit today."""
    return explanation_as_dict(explain(kind, values))


def explain_metrics(domain: str, values: dict | None = None) -> dict:
    """Per-metric two-layer explanations. Does not recompute science."""
    return explain_metrics_records(domain, values)
