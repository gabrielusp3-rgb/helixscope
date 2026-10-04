"""Serialization adapters must not mutate science or stringify NaN."""

from __future__ import annotations

import math

import pandas as pd

from helixscope_core.serialize import arrays_equivalent, dataframe_to_records, provenance_compare_fields
from helixscope_core.rna import codon_usage_table, codon_usage_records
import numpy as np


def test_dataframe_records_preserve_nan_as_float() -> None:
    frame = pd.DataFrame({"value": [1.0, float("nan")]})
    payload = dataframe_to_records(frame)
    assert payload["row_count"] == 2
    cell = payload["records"][1]["value"]
    assert isinstance(cell, float)
    assert math.isnan(cell)
    assert cell != "NaN"
    assert cell != 0


def test_codon_usage_adapter_does_not_mutate_source() -> None:
    seq = "AUGGCUGCUUAA"
    original = codon_usage_table(seq)
    before = original.copy(deep=True)
    _ = codon_usage_records(seq)
    pd.testing.assert_frame_equal(original, before)


def test_array_equivalence_and_provenance_omits_timestamp() -> None:
    left = np.array([1.0, float("nan")], dtype=float)
    right = np.array([1.0, float("nan")], dtype=float)
    assert arrays_equivalent(left, right)
    fields = provenance_compare_fields(
        {"status": "COMPUTED", "source": "user", "computed_at_utc": "2026-01-01T00:00:00Z"}
    )
    assert "computed_at_utc" not in fields
    assert fields["status"] == "COMPUTED"
