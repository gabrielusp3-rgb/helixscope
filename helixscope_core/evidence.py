"""Evidence items and conflicts. confidence_score remains None."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from modules.evidence_workspace import (
    evidence_item,
    export_pack,
    source_conflicts,
)

__all__ = ("evidence_conflict_pack", "evidence_item", "export_pack", "source_conflicts")


def evidence_conflict_pack(
    items: Sequence[Mapping[str, Any]],
    *,
    kind: str = "variants",
    inputs: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Export pack plus listed conflicts. Does not vote."""
    materialized = [dict(item) for item in items]
    pack = export_pack(kind=kind, inputs=dict(inputs or {}), items=materialized)
    pack["conflicts"] = source_conflicts(materialized)
    return pack
