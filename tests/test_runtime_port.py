"""PORT parsing for the container entrypoint. No shell and no network."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "helix_runtime_port", ROOT / "docker" / "runtime_port.py"
)
assert SPEC is not None and SPEC.loader is not None
runtime_port = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runtime_port)


def test_absent_port_uses_local_default() -> None:
    assert runtime_port.listen_port(None) == 8501


@pytest.mark.parametrize("raw", ["8501", "9000", "10000", "1", "65535"])
def test_decimal_ports_are_accepted(raw: str) -> None:
    assert runtime_port.listen_port(raw) == int(raw)


@pytest.mark.parametrize(
    "raw",
    [
        "",
        "abc",
        "-1",
        "0",
        "65536",
        "8501;command",
        "8501 && something",
        "8501\n",
        "8501\r",
        " 8501",
        "8501 ",
        "8501\x00",
        "$(id)",
        "8501`id`",
    ],
)
def test_malformed_port_is_rejected(raw: str) -> None:
    with pytest.raises(runtime_port.PortError):
        runtime_port.listen_port(raw)


def test_launcher_sources_do_not_call_a_shell() -> None:
    for name in ("entrypoint.py", "healthcheck.py", "runtime_port.py"):
        text = (ROOT / "docker" / name).read_text(encoding="utf-8")
        assert "shell=True" not in text
        assert "os.system" not in text
        assert "subprocess" not in text
