"""Request helpers. No science."""

from __future__ import annotations

from typing import Any

from fastapi import Request

from helixscope_api.schemas import ApiEnvelope
from helixscope_api.serialization import drop_heavy_fields, to_transport


def _string_list(value: Any) -> list[str]:
    if not value:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def envelope(
    request: Request,
    result: Any,
    *,
    scientific_status: str | None = None,
    execution_status: str | None = None,
    input_hash: str | None = None,
    warnings: list[str] | None = None,
    limitations: list[str] | None = None,
    drop_heavy: bool = False,
) -> ApiEnvelope:
    """Wrap a core value after transport encoding."""
    payload = drop_heavy_fields(result) if drop_heavy else result
    transported = to_transport(payload)
    if not isinstance(transported, dict):
        transported = {"value": transported}
    status = scientific_status or transported.get("status")
    hash_value = (
        input_hash
        or transported.get("sequence_hash")
        or transported.get("input_hash")
        or transported.get("identity_hash")
        or transported.get("alignment_hash")
    )
    warning_list = warnings if warnings is not None else _string_list(transported.get("warnings"))
    limit_list = limitations if limitations is not None else _string_list(transported.get("limitations"))
    return ApiEnvelope(
        request_id=str(getattr(request.state, "request_id", "") or ""),
        scientific_status=str(status) if status is not None else None,
        execution_status=execution_status,
        input_hash=str(hash_value) if hash_value else None,
        warnings=warning_list,
        limitations=limit_list,
        result=transported,
    )
