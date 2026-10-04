"""Replace this process with Streamlit on the configured PORT.

Uses os.execv so the Streamlit process keeps PID 1's signal handling under
Docker's init. No shell is started.
"""

from __future__ import annotations

import os
import sys

from runtime_port import PortError, listen_port


def main() -> None:
    """Exec Streamlit bound to 0.0.0.0 and the resolved port.

    Raises:
        SystemExit: PORT is set and invalid. os.execv does not return.
    """
    try:
        port = listen_port(os.environ.get("PORT"))
    except PortError as exc:
        print(f"HelixScope: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    os.execv(
        sys.executable,
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            "app.py",
            "--server.address=0.0.0.0",
            f"--server.port={port}",
            "--server.headless=true",
            "--server.fileWatcherType=none",
            "--browser.gatherUsageStats=false",
        ],
    )


if __name__ == "__main__":
    main()
