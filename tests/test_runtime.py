"""Contratos de runtime: stale modules, sem reload em massa."""

from __future__ import annotations

from modules import provenance, runtime
from ui import structure_viewer


def test_structure_viewer_api_sentinel_and_selection_parser():
    assert structure_viewer.STRUCTURE_VIEWER_API >= 14
    assert callable(structure_viewer.parse_plotly_selection)
    assert callable(structure_viewer.parse_plotly_selection_indices)
    assert callable(structure_viewer.figure_from_superposition)
    assert runtime.module_is_current("ui.structure_viewer", "figure_from_superposition")


def test_ensure_runtime_is_ok_for_current_version():
    runtime.reset_runtime_cache()
    report = runtime.ensure_runtime(expected_version=provenance.HELIXSCOPE_VERSION)
    assert report["ok"] is True
    assert report["needs_restart"] is False
    assert report["missing"] == []
    again = runtime.ensure_runtime(expected_version=provenance.HELIXSCOPE_VERSION)
    assert again["ok"] is True
    assert again["reloaded"] == []


def test_reload_if_missing_does_not_reload_when_complete():
    result = runtime.reload_if_missing(
        "ui.structure_viewer",
        ("parse_plotly_selection_indices", "figure_from_scene"),
    )
    assert result["reloaded"] is False
    assert result["missing_after"] == []


def test_require_attr_returns_existing_function():
    fn = runtime.require_attr(structure_viewer, "parse_plotly_selection_indices")
    assert fn is structure_viewer.parse_plotly_selection_indices
