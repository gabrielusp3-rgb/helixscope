"""Adapters for JSON-hostile internal types. Does not change stored science.

NaN is never rewritten as 0 or the string \"NaN\" inside science objects.
Transport encoding for FastAPI is Prompt 2.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np
import pandas as pd


def is_nonfinite_number(value: object) -> bool:
    """True for NaN or Inf. False for None, ints, finite floats, non-numbers."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating)):
        return False
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return not math.isfinite(number)


def dataframe_to_records(frame: pd.DataFrame) -> dict[str, Any]:
    """Canonical column/record view. Cell NaN stays float NaN in Python.

    Args:
        frame: Codon usage or similar table.

    Returns:
        columns, dtypes, records (row dicts), row_count. Index is positional
        only; a named index is recorded separately and not used as science.

    Raises:
        TypeError: If frame is not a DataFrame.
    """
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("dataframe_to_records expects a pandas DataFrame.")
    records = frame.to_dict(orient="records")
    return {
        "columns": [str(col) for col in frame.columns],
        "dtypes": {str(col): str(frame[col].dtype) for col in frame.columns},
        "records": records,
        "row_count": int(len(frame)),
        "index_name": frame.index.name,
        "ordered": True,
    }


def array_identity(array: np.ndarray) -> dict[str, Any]:
    """Shape/dtype fingerprint. Does not copy coordinates into lists."""
    arr = np.asarray(array)
    return {
        "dtype": str(arr.dtype),
        "shape": [int(dim) for dim in arr.shape],
        "c_contiguous": bool(arr.flags.c_contiguous),
        "nbytes": int(arr.nbytes),
    }


def arrays_equivalent(
    left: np.ndarray,
    right: np.ndarray,
    *,
    rtol: float = 0.0,
    atol: float = 0.0,
    equal_nan: bool = True,
) -> bool:
    """Numeric array parity helper. Default is exact including NaN positions."""
    return bool(
        np.array_equal(np.asarray(left), np.asarray(right), equal_nan=equal_nan)
        if rtol == 0.0 and atol == 0.0
        else np.allclose(
            np.asarray(left),
            np.asarray(right),
            rtol=rtol,
            atol=atol,
            equal_nan=equal_nan,
        )
    )


def redact_filesystem_path(path: object) -> str:
    """Basename only. Empty if missing. Never a home directory dump."""
    text = str(path or "").strip().replace("\\", "/")
    if not text:
        return ""
    return text.rsplit("/", 1)[-1]


def provenance_compare_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    """Fields that golden tests may compare. Timestamps omitted."""
    return {
        "status": record.get("status"),
        "source": record.get("source"),
        "algorithm": record.get("algorithm"),
        "software": record.get("software"),
        "software_version": record.get("software_version"),
        "input_hash": record.get("input_hash") or record.get("alignment_hash"),
        "database": record.get("database"),
        "accession": record.get("accession"),
    }
