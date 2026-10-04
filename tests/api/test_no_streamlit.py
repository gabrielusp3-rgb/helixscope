"""Importing the API must not load Streamlit."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
API = ROOT / "services" / "api"


def test_api_import_does_not_load_streamlit() -> None:
    script = (
        "import sys\n"
        "from helixscope_api.main import create_app\n"
        "create_app()\n"
        "mods = [name for name in sys.modules if name == 'streamlit' or name.startswith('streamlit.')]\n"
        "raise SystemExit(0 if not mods else 'STREAMLIT=' + ','.join(mods))\n"
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT), str(API), env.get("PYTHONPATH", "")])
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
