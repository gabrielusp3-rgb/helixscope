"""Reusable scientific metric grid: balanced rows, missing values, XSS."""

from __future__ import annotations

from ui.components import html_escape
from ui.workspace import metric_display_value, metric_grid_html, metric_strip_html


def test_metric_grid_uses_fixed_columns_not_stretching_orphans() -> None:
    html = metric_grid_html(
        [
            ("Length", "23", "bp"),
            ("GC Content", "41.49", "%"),
            ("AT Content", "58.51", "%"),
            ("Tm", "Unavailable", ""),
            ("MW", "7.08", "kDa"),
            ("Shannon entropy", "1.95", "bits"),
            ("AT skew", "0.01", ""),
            ("GC skew", "-0.02", ""),
        ]
    )
    assert 'data-cols="4"' in html
    assert "--hs-metric-cols:4" in html
    assert "auto-fit" not in html
    assert html.count('class="hs-metric"') == 8
    assert "41.49" in html
    assert "hs-metric-unit" in html
    assert "Unavailable" in html


def test_metric_grid_preserves_na_and_never_invents_zero() -> None:
    html = metric_grid_html(
        [
            ("Tm", "N/A", "C"),
            ("ENC (Nc)", None, ""),
            ("Count", 0, ""),
        ]
    )
    assert "N/A" in html
    assert "Unavailable" in html
    assert ">0<" in html or ">0</span>" in html
    assert "0.00" not in html
    assert metric_display_value(None) == "Unavailable"
    assert metric_display_value("N/A") == "N/A"
    assert metric_display_value(0) == "0"


def test_metric_strip_delegates_to_grid() -> None:
    html = metric_strip_html((("Block RMSD", "0.8", "A"), ("Global RMSD", "0.99", "A")))
    assert "hs-metric-grid" in html
    assert "0.8" in html
    assert "0.99" in html


def test_metric_grid_escapes_untrusted_text() -> None:
    payload = '<img src=x onerror="alert(1)">'
    html = metric_grid_html([(payload, payload, payload)])
    assert "<img" not in html
    assert "&lt;" in html
    assert html_escape(payload) in html


def test_nine_metrics_do_not_use_five_plus_three_columns() -> None:
    items = [(f"M{i}", str(i), "") for i in range(9)]
    html = metric_grid_html(items)
    assert 'data-cols="4"' in html
    assert 'data-cols="5"' not in html
    assert 'data-cols="3"' not in html
