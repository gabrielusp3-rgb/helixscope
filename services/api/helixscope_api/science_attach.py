"""Attach deterministic explanations. Does not recompute science."""

from __future__ import annotations

from typing import Any, Mapping

from helixscope_core.explain import explain_metrics, explain_result


def attach_explanations(
    payload: dict[str, Any],
    *,
    include: bool,
    kind: str,
    domain: str,
    values: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Add domain explanation plus per-metric explanations when requested."""
    if not include:
        return payload
    source = dict(values) if isinstance(values, Mapping) else payload
    payload["explanation"] = explain_result(kind, source)
    payload["metric_explanations"] = explain_metrics(domain, payload)
    return payload
