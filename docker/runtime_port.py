"""TCP port used by the HelixScope container.

The local Compose service leaves PORT unset or sets 8501. A public host may
inject another PORT. This parser accepts only a decimal integer in the TCP
range. It does not strip, evaluate, or pass the value to a shell.
"""

from __future__ import annotations

DEFAULT_LISTEN_PORT: int = 8501
MIN_LISTEN_PORT: int = 1
MAX_LISTEN_PORT: int = 65535


class PortError(ValueError):
    """PORT is present and is not a usable TCP port."""


def listen_port(raw: str | None) -> int:
    """Resolve the port Streamlit should bind.

    Args:
        raw: Value of the PORT environment variable. None means the variable
            is absent and the local default is used.

    Returns:
        An integer from 1 to 65535.

    Raises:
        PortError: The value is empty, padded, non-decimal, or outside
            1..65535.
    """
    if raw is None:
        return DEFAULT_LISTEN_PORT
    if not isinstance(raw, str) or raw == "" or any(char not in "0123456789" for char in raw):
        raise PortError("PORT must be an integer from 1 to 65535.")
    value = int(raw)
    if value < MIN_LISTEN_PORT or value > MAX_LISTEN_PORT:
        raise PortError("PORT must be an integer from 1 to 65535.")
    return value
