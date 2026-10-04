"""Pointer-reactive glass highlight via Streamlit Components v2.

Trusted application script only. Does not receive NCBI, ClinVar, FASTA, or
user strings. isolate_styles is False so the listener can reach glass nodes
in the app document. Replaces the deprecated Components v1 HTML iframe.
Registration happens at mount time so Streamlit AppTest can run the script.
"""

from __future__ import annotations

from typing import Any, Optional

import streamlit as st

POINTER_JS: str = """
export default function(component) {
  const root = document;
  function onMove(event) {
    if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      return;
    }
    const cards = root.querySelectorAll(
      "[data-hs-glass], .hs-glass, .hs-glass-thin, .hs-glass-medium, .hs-glass-thick, [data-testid='stHeader'], section[data-testid='stSidebar'] > div:first-child, [data-testid='stDialog']"
    );
    cards.forEach(function (el) {
      const box = el.getBoundingClientRect();
      if (
        event.clientX < box.left ||
        event.clientX > box.right ||
        event.clientY < box.top ||
        event.clientY > box.bottom
      ) {
        return;
      }
      const mx = ((event.clientX - box.left) / Math.max(box.width, 1)) * 100;
      const my = ((event.clientY - box.top) / Math.max(box.height, 1)) * 100;
      el.style.setProperty("--mx", mx.toFixed(2) + "%");
      el.style.setProperty("--my", my.toFixed(2) + "%");
    });
  }
  root.addEventListener("pointermove", onMove, { passive: true });
  const uploadMessage = "The upload was rejected before analysis. Use a simple file name such as sample.fasta. The file was not stored or executed.";
  function softenUploadErrors(node) {
    if (!node || !node.querySelectorAll) {
      return;
    }
    const alerts = node.querySelectorAll("[data-testid='stFileUploader'] [data-testid='stAlert'], [data-testid='stFileUploader'] [role='alert']");
    alerts.forEach(function (el) {
      const text = el.textContent || "";
      if (text.indexOf("AxiosError") !== -1 || text.indexOf("status code 403") !== -1) {
        if (el.textContent !== uploadMessage) {
          el.textContent = uploadMessage;
        }
      }
    });
  }
  softenUploadErrors(root);
  const observer = new MutationObserver(function () {
    softenUploadErrors(root);
  });
  if (root.body) {
    observer.observe(root.body, { childList: true, subtree: true, characterData: true });
  }
  return function () {
    root.removeEventListener("pointermove", onMove);
    observer.disconnect();
  };
}
"""

HOTKEY_JS: str = """
export default function(component) {
  const { setTriggerValue } = component;
  function onKey(event) {
    if ((event.ctrlKey || event.metaKey) && String(event.key).toLowerCase() === "k") {
      event.preventDefault();
      setTriggerValue("open", Date.now());
    }
  }
  document.addEventListener("keydown", onKey);
  return function () {
    document.removeEventListener("keydown", onKey);
  };
}
"""


def _safe_mount(name: str, js: str, *, key: str, extra: Optional[dict] = None) -> Any:
    try:
        renderer = st.components.v2.component(name, js=js, isolate_styles=False)
        kwargs = dict(extra or {})
        return renderer(key=key, **kwargs)
    except Exception:
        return None


def mount_pointer_effect() -> None:
    """Attach the glass specular listener. No-op if Components v2 cannot mount.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    _safe_mount("helixscope_pointer_glass", POINTER_JS, key="helix_pointer_glass")


def mount_command_hotkey() -> Any:
    """Listen for Ctrl/Cmd+K. Returns None when the component cannot mount.

    Args:
        Nenhum.

    Returns:
        ComponentResult or None.

    Raises:
        Nenhum.
    """
    return _safe_mount(
        "helixscope_command_hotkey",
        HOTKEY_JS,
        key="helix_command_hotkey",
        extra={"on_open_change": lambda: None},
    )
