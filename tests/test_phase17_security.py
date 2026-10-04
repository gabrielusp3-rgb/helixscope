"""Red team Fase 17: SSRF Alignment API, UUID, entry_id, subprocess, XSS overlay."""

from __future__ import annotations

import inspect

import pytest

from modules import comparative, rcsb_alignment, structure_superposition, usalign
from ui import components, structure_viewer


def test_alignment_client_refuses_user_urls_and_private_targets():
    hostile = [
        "file:///etc/passwd",
        "gopher://alignment.rcsb.org/api/v1/structures/submit",
        "ftp://alignment.rcsb.org/api/v1/structures/submit",
        "http://alignment.rcsb.org/api/v1/structures/submit",
        "https://169.254.169.254/latest/meta-data/",
        "https://127.0.0.1/api/v1/structures/submit",
        "https://[::1]/api/v1/structures/submit",
        "https://alignment.rcsb.org/api/v1/../etc/passwd",
    ]
    for url in hostile:
        assert rcsb_alignment.url_is_allowed(url) is False, url


def test_submit_never_takes_a_url_argument():
    source = inspect.getsource(rcsb_alignment.submit_pairwise)
    assert "url=" not in source
    assert "SUBMIT_URL" in source


def test_transform_rejects_non_finite_matrix():
    with pytest.raises(structure_superposition.SuperpositionError) as exc:
        structure_superposition.apply_column_major_4x4(0, 0, 0, [float("nan")] * 16)
    assert exc.value.category == "INVALID_INPUT"


def test_plotly_superposition_text_strips_html():
    safe = structure_viewer.plotly_safe_text('<img src=x onerror="alert(1)">')
    assert "<" not in safe
    assert ">" not in safe
    card = components.glass_card('<script>alert(1)</script>', "body")
    assert "<script>" not in card


def test_usalign_and_comparative_have_no_shell_true():
    assert "shell=True" not in inspect.getsource(usalign)
    assert "shell=True" not in inspect.getsource(comparative)
    assert "shell=True" not in inspect.getsource(rcsb_alignment)
