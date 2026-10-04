"""AppTest helpers for Phase 19 sidebar routing.

Not collected as tests. Existing AppTests import open_module / switch_module
so they render the scientific widgets of one module instead of the old tab strip.
"""

from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest

from ui.navigation import ACTIVE_MODULE_KEY

ROOT = Path(__file__).resolve().parents[1]


def open_module(module_id: str, *, timeout: int = 90) -> AppTest:
    """Start app.py on one semantic module.

    Args:
        module_id: Registry id such as dna, protein, msa.
        timeout: AppTest default timeout in seconds.

    Returns:
        A run AppTest instance.

    Raises:
        Nenhum. Callers assert app.exception.
    """
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=timeout)
    app.session_state[ACTIVE_MODULE_KEY] = module_id
    app.run()
    return app


def switch_module(app: AppTest, module_id: str) -> AppTest:
    """Switch the active module without clearing scientific session keys.

    Args:
        app: Running AppTest.
        module_id: Registry id.

    Returns:
        The same AppTest after rerun.

    Raises:
        Nenhum.
    """
    app.session_state[ACTIVE_MODULE_KEY] = module_id
    app.run()
    return app
