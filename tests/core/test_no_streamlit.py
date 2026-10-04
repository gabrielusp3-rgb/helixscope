"""Static and runtime guards: helixscope_core never imports Streamlit."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "helixscope_core"


def test_core_source_has_zero_streamlit_imports() -> None:
    offenders: list[str] = []
    for path in sorted(CORE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "streamlit" or alias.name.startswith("streamlit."):
                        offenders.append(str(path.relative_to(ROOT)))
            if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("streamlit"):
                offenders.append(str(path.relative_to(ROOT)))
    assert offenders == []


def test_public_core_import_does_not_load_streamlit() -> None:
    script = (
        "import sys\n"
        "import helixscope_core\n"
        "from helixscope_core import analyze_dna, get_capabilities, explain_result\n"
        "from helixscope_core.dna import gc_content\n"
        "from helixscope_core.phylogeny import infer_tree\n"
        "from helixscope_core.structure import classify_mapping_status\n"
        "from helixscope_core.compare import parse_alignment_payload\n"
        "mods = [name for name in sys.modules if name == 'streamlit' or name.startswith('streamlit.')]\n"
        "raise SystemExit(0 if not mods else 'STREAMLIT=' + ','.join(mods))\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "streamlit" not in proc.stderr.lower()
