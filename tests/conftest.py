"""Ensure the API package is importable from the repository root."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
API_ROOT = ROOT / "services" / "api"
for path in (str(ROOT), str(API_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)
