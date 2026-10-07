"""Intro overlay is inserted once, muted, and does not replace the app mount."""

from __future__ import annotations

from pathlib import Path

from ui.intro_splash import INTRO_ASSET, INTRO_MARKER, patch_index_html

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = """\
<!DOCTYPE html>
<html>
  <body>
    <div id="root"></div>
  </body>
</html>
"""


def test_intro_asset_is_the_tracked_clip() -> None:
    clip = ROOT / "static" / "helixscope-intro.mp4"
    assert clip.is_file()
    assert 100_000 < clip.stat().st_size < 2_000_000


def test_static_serving_is_enabled() -> None:
    text = (ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8")
    assert "enableStaticServing = true" in text


def test_patch_inserts_muted_overlay_once() -> None:
    once = patch_index_html(SAMPLE)
    assert once.count(INTRO_MARKER) == 1
    assert once.count('<div id="root"></div>') == 1
    assert once.index(INTRO_MARKER) < once.index('<div id="root"></div>')
    assert 'src="/app/static/helixscope-intro.mp4"' in once
    assert INTRO_ASSET == "/app/static/helixscope-intro.mp4"
    assert "muted" in once
    assert "autoplay" in once
    assert "object-fit: contain" in once
    assert patch_index_html(once) == once


def test_patch_rejects_an_index_without_a_root_mount() -> None:
    try:
        patch_index_html("<html></html>")
    except ValueError as exc:
        assert "root mount" in str(exc)
    else:
        raise AssertionError("expected ValueError")
