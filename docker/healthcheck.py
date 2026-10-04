"""Probe Streamlit on the same PORT the entrypoint binds.

Exit 0 only for an HTTP 2xx response from /_stcore/health. Uses the standard
library. No shell.
"""

from __future__ import annotations

import os
import sys
import urllib.error
import urllib.request

from runtime_port import PortError, listen_port

HEALTH_PATH: str = "/_stcore/health"
TIMEOUT_S: float = 3.0


def main() -> int:
    """Return 0 when the local health endpoint answers 2xx.

    Returns:
        Process status: 0 healthy, 1 otherwise.
    """
    try:
        port = listen_port(os.environ.get("PORT"))
    except PortError:
        return 1
    url = f"http://127.0.0.1:{port}{HEALTH_PATH}"
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT_S) as response:
            status = getattr(response, "status", None)
            if status is None:
                status = response.getcode()
            if 200 <= int(status) < 300:
                return 0
    except (OSError, urllib.error.URLError, TimeoutError, ValueError):
        return 1
    return 1


if __name__ == "__main__":
    sys.exit(main())
