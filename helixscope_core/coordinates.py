"""Sequence/genome coordinate conversions. Internal is 0-based half-open."""

from __future__ import annotations

from modules.genome_coordinates import (
    DISPLAY_CONVENTION,
    INTERNAL_CONVENTION,
    CoordinateError,
    display_to_internal,
    internal_to_display,
    roundtrip_ok,
)

__all__ = (
    "DISPLAY_CONVENTION",
    "INTERNAL_CONVENTION",
    "CoordinateError",
    "display_to_internal",
    "internal_to_display",
    "roundtrip_ok",
)
