"""Patch the installed Streamlit index so the intro plays during boot."""

from __future__ import annotations

import sys
from pathlib import Path

_APP = Path(__file__).resolve().parents[1]
if str(_APP) not in sys.path:
    sys.path.insert(0, str(_APP))

from ui.intro_splash import patch_installed_streamlit


def main() -> None:
    """Patch Streamlit's index.html in place."""
    path = patch_installed_streamlit()
    print(f"HelixScope intro: {path}")


if __name__ == "__main__":
    main()
