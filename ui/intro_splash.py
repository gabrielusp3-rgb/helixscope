"""Full-viewport intro played over Streamlit's boot spinner.

The clip has to be in the document Streamlit serves before React mounts.
A widget created after ``app.py`` runs appears only once that spinner is
already on screen, so the markup is inserted into Streamlit's ``index.html``
at image build time. A full page load (first visit, refresh, or a new tab)
plays it once, muted. Navigation inside the running app does not reload
``index.html``, so the clip does not replay on every click.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

INTRO_MARKER = 'id="hs-intro"'
ROOT_MOUNT = '<div id="root"></div>'
INTRO_ASSET = "/app/static/helixscope-intro.mp4"

# The safety cap is longer than the clip so a slow first byte does not cut
# playback. If autoplay is refused or the file fails, the overlay closes
# immediately and the app stays reachable.
_SAFETY_MS = 20000

INTRO_MARKUP = f"""\
    <style>
      #hs-intro {{
        position: fixed;
        inset: 0;
        z-index: 2147483647;
        background: #000000;
        display: flex;
        align-items: center;
        justify-content: center;
      }}
      #hs-intro video {{
        width: 100%;
        height: 100%;
        object-fit: contain;
        background: #000000;
      }}
    </style>
    <div id="hs-intro">
      <video id="hs-intro-video" autoplay muted playsinline preload="auto">
        <source src="{INTRO_ASSET}" type="video/mp4" />
      </video>
    </div>
    <script>
      (function () {{
        var box = document.getElementById("hs-intro");
        var video = document.getElementById("hs-intro-video");
        if (!box || !video) return;
        var closed = false;
        function closeIntro() {{
          if (closed) return;
          closed = true;
          box.style.transition = "opacity 240ms linear";
          box.style.opacity = "0";
          window.setTimeout(function () {{
            try {{
              video.pause();
              video.removeAttribute("src");
              var source = video.querySelector("source");
              if (source) source.removeAttribute("src");
              video.load();
            }} catch (e) {{}}
            box.remove();
          }}, 260);
        }}
        video.muted = true;
        video.defaultMuted = true;
        video.volume = 0;
        video.playsInline = true;
        video.setAttribute("playsinline", "");
        video.addEventListener("ended", closeIntro);
        video.addEventListener("error", closeIntro);
        var started = video.play();
        if (started && typeof started.catch === "function") {{
          started.catch(closeIntro);
        }}
        window.setTimeout(closeIntro, {_SAFETY_MS});
      }})();
    </script>
"""


def patch_index_html(text: str) -> str:
    """Insert the intro overlay in front of Streamlit's root mount.

    Args:
        text: Contents of Streamlit ``static/index.html``.

    Returns:
        The same text when the marker is already present, otherwise the
        text with the overlay inserted once.

    Raises:
        ValueError: The root mount Streamlit uses is not in ``text``.
    """
    if INTRO_MARKER in text:
        return text
    if ROOT_MOUNT not in text:
        raise ValueError("Streamlit index.html has no root mount")
    return text.replace(ROOT_MOUNT, INTRO_MARKUP + "    " + ROOT_MOUNT, 1)


def installed_index_path() -> Path:
    """Return the ``index.html`` shipped inside the installed Streamlit package."""
    spec = importlib.util.find_spec("streamlit")
    if spec is None or not spec.origin:
        raise RuntimeError("streamlit is not installed")
    return Path(spec.origin).resolve().parent / "static" / "index.html"


def patch_installed_streamlit(index_path: Path | None = None) -> Path:
    """Write the intro overlay into Streamlit's ``index.html`` if it is absent.

    Args:
        index_path: File to patch. Defaults to the installed package copy.

    Returns:
        The path that was checked.

    Raises:
        FileNotFoundError: ``index_path`` does not exist.
        ValueError: The file has no root mount and no intro marker.
    """
    path = installed_index_path() if index_path is None else index_path
    if not path.is_file():
        raise FileNotFoundError(path)
    original = path.read_text(encoding="utf-8")
    patched = patch_index_html(original)
    if patched != original:
        path.write_text(patched, encoding="utf-8", newline="\n")
    return path
