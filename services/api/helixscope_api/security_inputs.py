"""Reject URL/path-shaped identifiers. Destinations stay allowlisted in core."""

from __future__ import annotations

from helixscope_api.errors import ApiError


def reject_url_or_path(value: str, *, field: str) -> None:
    """User input is data, never a fetch target or executable path."""
    text = str(value or "")
    lowered = text.lower().strip()
    if (
        "://" in text
        or ".." in text
        or "\\" in text
        or text.startswith("/")
        or lowered.startswith("file:")
        or ".exe" in lowered
    ):
        raise ApiError(400, "INVALID_INPUT", f"{field} must not be a URL, path, or executable.")
