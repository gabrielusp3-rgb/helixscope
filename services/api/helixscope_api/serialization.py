"""Strict JSON scientific transport. Does not change core values, only encoding.

Finite numbers including 0 stay JSON numbers.
NaN / Inf / -Inf become {value: null, value_state: NAN|INF|-INF}.
None stays null without a value_state (unavailable / not computed).
"""

from __future__ import annotations

import math
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from helixscope_core.serialize import dataframe_to_records, is_nonfinite_number, redact_filesystem_path

VALUE_STATE_NAN: str = "NAN"
VALUE_STATE_INF: str = "INF"
VALUE_STATE_NINF: str = "-INF"


def nonfinite_token(value: object) -> dict[str, Any] | None:
    """Wrapper for NaN/Inf. None if the value is finite or non-numeric."""
    if not is_nonfinite_number(value):
        return None
    number = float(value)  # type: ignore[arg-type]
    if math.isnan(number):
        return {"value": None, "value_state": VALUE_STATE_NAN}
    if number > 0:
        return {"value": None, "value_state": VALUE_STATE_INF}
    return {"value": None, "value_state": VALUE_STATE_NINF}


def encode_array(array: np.ndarray) -> dict[str, Any]:
    """dtype/shape plus nested Python data. Does not mutate the source array."""
    arr = np.asarray(array)
    snapshot = arr.copy()
    return {
        "dtype": str(snapshot.dtype),
        "shape": [int(dim) for dim in snapshot.shape],
        "c_contiguous": bool(arr.flags.c_contiguous),
        "data": to_transport(snapshot.tolist()),
    }


def to_transport(value: Any) -> Any:
    """Walk a core result into strict-JSON-safe Python objects."""
    token = nonfinite_token(value)
    if token is not None:
        return token
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        return value
    if isinstance(value, (datetime, date)):
        text = value.isoformat()
        if isinstance(value, datetime) and value.tzinfo is None:
            return text + "Z"
        return text
    if isinstance(value, Path):
        return redact_filesystem_path(value)
    if isinstance(value, np.ndarray):
        return encode_array(value)
    if isinstance(value, np.generic):
        return to_transport(value.item())
    if isinstance(value, pd.DataFrame):
        payload = dataframe_to_records(value)
        payload["transport_string_dtype"] = "string"
        payload["dtypes"] = {
            name: ("string" if str(dtype) in {"object", "str", "string"} else str(dtype))
            for name, dtype in payload["dtypes"].items()
        }
        payload["records"] = to_transport(payload["records"])
        return payload
    if isinstance(value, dict):
        out: dict[str, Any] = {}
        for key, item in value.items():
            low = str(key).lower()
            if any(part in low for part in ("password", "api_key", "secret", "token", "email")):
                continue
            if low == "path" or low.endswith("_path") or low.endswith("filepath"):
                out[str(key)] = redact_filesystem_path(item)
            else:
                out[str(key)] = to_transport(item)
        return out
    if isinstance(value, (list, tuple)):
        return [to_transport(item) for item in value]
    if hasattr(value, "__fspath__"):
        return redact_filesystem_path(value)
    return str(value)[:500]


HEAVY_STRUCTURE_KEYS: frozenset[str] = frozenset(
    {"raw_text", "raw_structure_text", "structure_text", "mmcif_text", "pdb_text"}
)


def drop_heavy_fields(value: Any) -> Any:
    """Remove bulky coordinate-file blobs from API payloads. Values stay otherwise intact."""
    if isinstance(value, dict):
        return {
            str(key): drop_heavy_fields(item)
            for key, item in value.items()
            if str(key) not in HEAVY_STRUCTURE_KEYS
        }
    if isinstance(value, list):
        return [drop_heavy_fields(item) for item in value]
    return value


def assert_strict_json(payload: Any) -> None:
    """Raise ValueError if payload cannot dump with allow_nan=False."""
    import json

    json.dumps(payload, allow_nan=False)
