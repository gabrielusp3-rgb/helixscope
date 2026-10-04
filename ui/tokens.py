"""HelixScope design tokens.

Visual language only. These values do not compute GC, RMSD, Tm, or any
scientific quantity. Semantic status colors are identity, not good-vs-bad.

UI identity is monochrome OLED liquid glass: black and white only.
Scientific series hexes (DNA/RNA/protein, 3D traces) stay frozen.
Decorative chrome does not use those hexes. Status is labeled in text.
--hs-cyan is kept as a CSS variable name for compatibility; its value is
ice (scientific), not turquoise.
"""

from __future__ import annotations

from typing import Mapping

# ---------------------------------------------------------------------------
# UI identity (monochrome OLED)
# ---------------------------------------------------------------------------
HS_VOID: str = "#000000"
HS_BG_DEEP: str = "#000000"
HS_BG: str = "#000000"
HS_BG_RAISED: str = "#000000"
HS_SURFACE_1: str = "#000000"
HS_SURFACE_2: str = "#000000"
HS_SURFACE_TABLE: str = "#000000"

HS_TEXT: str = "#FFFFFF"
HS_TEXT_SECONDARY: str = "rgba(255, 255, 255, 0.78)"
HS_TEXT_BODY: str = "rgba(255, 255, 255, 0.92)"
HS_TEXT_SUBTITLE: str = "rgba(255, 255, 255, 0.86)"
HS_TEXT_MUTED: str = "rgba(255, 255, 255, 0.55)"
HS_TEXT_META: str = "rgba(255, 255, 255, 0.45)"

HS_AMBIENT: str = "rgba(255, 255, 255, 0.04)"
HS_GRAIN: str = "rgba(255, 255, 255, 0.028)"

# Layer 2 — glass / controls (thickness variants)
HS_GLASS_THIN: str = "rgba(0, 0, 0, 0.42)"
HS_GLASS_MEDIUM: str = "rgba(0, 0, 0, 0.60)"
HS_GLASS_THICK: str = "rgba(0, 0, 0, 0.72)"
HS_GLASS_FALLBACK: str = "#000000"
HS_GLASS: str = HS_GLASS_MEDIUM
HS_GLASS_STRONG: str = HS_GLASS_THICK
HS_GLASS_BORDER: str = "rgba(255, 255, 255, 0.08)"
HS_GLASS_HIGHLIGHT: str = "rgba(255, 255, 255, 0.12)"
HS_GLASS_SHADOW: str = "rgba(0, 0, 0, 0.72)"
HS_GLASS_BLEED: str = "rgba(255, 255, 255, 0.06)"
HS_GLASS_ACTIVE: str = "rgba(255, 255, 255, 0.08)"
HS_BORDER: str = "rgba(255, 255, 255, 0.08)"
HS_BORDER_HOVER: str = "rgba(255, 255, 255, 0.20)"
HS_HIGHLIGHT: str = "rgba(255, 255, 255, 0.12)"
HS_FOCUS: str = "rgba(255, 255, 255, 0.55)"
HS_SHADOW: str = "rgba(0, 0, 0, 0.72)"

HS_PRIMARY: str = "#FFFFFF"
HS_BACKGROUND_DEEP: str = HS_VOID
HS_BACKGROUND: str = HS_BG_DEEP

# ---------------------------------------------------------------------------
# Frozen scientific series (charts, 3D traces, molecule legends, sequence)
# ---------------------------------------------------------------------------
HS_DEEP_NAVY: str = "#001F59"
HS_NAVY: str = "#093986"
HS_ELECTRIC_BLUE: str = "#2355A0"
HS_MID_BLUE: str = "#6893D0"
HS_ICE_BLUE: str = "#91B4E4"
HS_FROST: str = "#AFCBF1"
HS_WHITE_BLUE: str = "#D9E4FC"
HS_CYAN: str = HS_ICE_BLUE
HS_SERIES: str = HS_MID_BLUE

# Spacing (px)
SPACE_1: int = 4
SPACE_2: int = 8
SPACE_3: int = 12
SPACE_4: int = 16
SPACE_5: int = 24
SPACE_6: int = 32

# Radius (px)
RADIUS_SM: int = 8
RADIUS_MD: int = 12
RADIUS_LG: int = 20

# Motion
MOTION_FAST_MS: int = 140
MOTION_STANDARD_MS: int = 220
MOTION_CONTEXTUAL_MS: int = 280
MOTION_EASING: str = "cubic-bezier(0.22, 1, 0.36, 1)"

# Type
FONT_SANS: str = "Inter, system-ui, sans-serif"
FONT_DISPLAY: str = "'Space Grotesk', Inter, system-ui, sans-serif"
FONT_MONO: str = "'JetBrains Mono', ui-monospace, SFMono-Regular, Consolas, monospace"

# Plotly chrome only. Scientific arrays are never restyled here.
PLOTLY_PAPER: str = "#000000"
PLOTLY_PLOT: str = "#000000"
PLOTLY_FONT: str = "#FFFFFF"
PLOTLY_GRID: str = "rgba(255, 255, 255, 0.06)"
PLOTLY_3D_BG: str = HS_VOID
STRUCTURE_VIEWPORT_HEIGHT: int = 640
STRUCTURE_VIEWPORT_HEIGHT_EXPANDED: int = 820

# Status is text plus opacity. Hue is not the signal. 3D scientific colors
# live in the frozen series above and are not reused here.
STATUS_COLORS: dict[str, dict[str, str]] = {
    "EXPERIMENTAL": {
        "fg": "rgba(255, 255, 255, 1)",
        "bg": "rgba(255, 255, 255, 0.10)",
        "border": "rgba(255, 255, 255, 0.70)",
        "mark": "[+]",
    },
    "PREDICTED": {
        "fg": "rgba(255, 255, 255, 0.78)",
        "bg": "rgba(255, 255, 255, 0.08)",
        "border": "rgba(255, 255, 255, 0.46)",
        "mark": "[~]",
    },
    "ILLUSTRATIVE": {
        "fg": "rgba(255, 255, 255, 0.62)",
        "bg": "rgba(0, 0, 0, 0.40)",
        "border": "rgba(255, 255, 255, 0.28)",
        "mark": "[.]",
    },
    "COMPUTED": {
        "fg": "rgba(255, 255, 255, 0.96)",
        "bg": "rgba(255, 255, 255, 0.08)",
        "border": "rgba(255, 255, 255, 0.34)",
        "mark": "[=]",
    },
    "RETRIEVED": {
        "fg": "rgba(255, 255, 255, 0.88)",
        "bg": "rgba(255, 255, 255, 0.07)",
        "border": "rgba(255, 255, 255, 0.32)",
        "mark": "[>]",
    },
    "HEURISTIC": {
        "fg": "rgba(255, 255, 255, 0.70)",
        "bg": "rgba(0, 0, 0, 0.35)",
        "border": "rgba(255, 255, 255, 0.30)",
        "mark": "[h]",
    },
    "PARTIAL": {
        "fg": "rgba(255, 255, 255, 0.58)",
        "bg": "rgba(0, 0, 0, 0.35)",
        "border": "rgba(255, 255, 255, 0.24)",
        "mark": "[/]",
    },
    "UNMAPPED": {
        "fg": "rgba(255, 255, 255, 0.42)",
        "bg": "rgba(0, 0, 0, 0.50)",
        "border": "rgba(255, 255, 255, 0.18)",
        "mark": "[?]",
    },
    "UNAVAILABLE": {
        "fg": "rgba(255, 255, 255, 0.40)",
        "bg": "rgba(0, 0, 0, 0.55)",
        "border": "rgba(255, 255, 255, 0.16)",
        "mark": "[-]",
    },
    "NOT_INSTALLED": {
        "fg": "rgba(255, 255, 255, 0.36)",
        "bg": "rgba(0, 0, 0, 0.55)",
        "border": "rgba(255, 255, 255, 0.14)",
        "mark": "[o]",
    },
    "RESOURCE_LIMIT": {
        "fg": "rgba(255, 255, 255, 0.66)",
        "bg": "rgba(0, 0, 0, 0.40)",
        "border": "rgba(255, 255, 255, 0.26)",
        "mark": "[!]",
    },
    "ERROR": {
        "fg": "rgba(255, 255, 255, 0.90)",
        "bg": "rgba(0, 0, 0, 0.45)",
        "border": "rgba(255, 255, 255, 0.80)",
        "mark": "[x]",
    },
    "LIVE_VALIDATED": {
        "fg": "rgba(255, 255, 255, 0.98)",
        "bg": "rgba(255, 255, 255, 0.14)",
        "border": "rgba(255, 255, 255, 0.64)",
        "mark": "[v]",
    },
    "REMOTE_VALIDATED": {
        "fg": "rgba(255, 255, 255, 0.74)",
        "bg": "rgba(255, 255, 255, 0.06)",
        "border": "rgba(255, 255, 255, 0.36)",
        "mark": "[r]",
    },
    "STALE": {
        "fg": "rgba(255, 255, 255, 0.52)",
        "bg": "rgba(0, 0, 0, 0.40)",
        "border": "rgba(255, 255, 255, 0.22)",
        "mark": "[s]",
    },
    "TEST_ONLY": {
        "fg": "rgba(255, 255, 255, 0.48)",
        "bg": "rgba(0, 0, 0, 0.45)",
        "border": "rgba(255, 255, 255, 0.20)",
        "mark": "[t]",
    },
}

STATUS_ALIASES: dict[str, str] = {
    "VALIDATED": "EXPERIMENTAL",
    "WARNING": "PARTIAL",
    "NO_HITS": "PARTIAL",
    "WAITING": "PARTIAL",
    "READY": "RETRIEVED",
    "TIMEOUT": "ERROR",
    "RATE_LIMITED": "RESOURCE_LIMIT",
    "SERVICE_UNAVAILABLE": "ERROR",
    "INVALID_INPUT": "ERROR",
    "INVALID_DATABASE": "ERROR",
    "JOB_FAILED": "ERROR",
    "PARSING_ERROR": "ERROR",
    "FAILED": "ERROR",
    "CACHED": "RETRIEVED",
    "QUEUED": "PARTIAL",
    "RUNNING": "PARTIAL",
    "COMPLETED": "COMPUTED",
    "CANCELLED": "PARTIAL",
    "NOT_FOUND": "ERROR",
    "NETWORK_ERROR": "ERROR",
    "INSUFFICIENT_DATA": "ERROR",
    "NO_STRUCTURE": "UNAVAILABLE",
    "TOOL_NOT_INSTALLED": "NOT_INSTALLED",
    "TOOL_FAILED": "ERROR",
    "SEQUENCE_UNAVAILABLE": "UNAVAILABLE",
    "EXACT": "RETRIEVED",
    "ALIGNED": "COMPUTED",
    "BEST_EFFORT": "PARTIAL",
    "UNCERTAIN": "PARTIAL",
    "TEST_REFERENCE": "TEST_ONLY",
    "PARTIAL_VIEW": "PARTIAL",
    "SCOPE_MISMATCH": "PARTIAL",
    "NOT_RUN": "UNAVAILABLE",
    "INCOMPLETE": "PARTIAL",
    "VERIFIED": "COMPUTED",
    "AVAILABLE": "RETRIEVED",
    "SKIPPED": "UNAVAILABLE",
    "NOT_COMPUTED": "UNAVAILABLE",
    "MAPPED": "COMPUTED",
    "NORMALIZED": "COMPUTED",
    "IDENTIFIER_ONLY": "PARTIAL",
    "REF_VALIDATED": "COMPUTED",
    "REF_NOT_CHECKED": "PARTIAL",
    "REFERENCE_MISMATCH": "ERROR",
    "UNSUPPORTED": "UNAVAILABLE",
}


def canonical_status(status: object) -> str:
    """Map a status string onto a tokenized semantic family.

    Args:
        status: Raw status from a scientific result object.

    Returns:
        A key present in STATUS_COLORS, or UNAVAILABLE.

    Raises:
        Nenhum.
    """
    key = str(status or "UNAVAILABLE").strip().upper().replace(" ", "_")
    if key in STATUS_COLORS:
        return key
    return STATUS_ALIASES.get(key, "UNAVAILABLE")


def status_style(status: object) -> Mapping[str, str]:
    """Return fg/bg/border for a scientific status.

    Args:
        status: Raw status token.

    Returns:
        Mapping with fg, bg, border.

    Raises:
        Nenhum.
    """
    return STATUS_COLORS[canonical_status(status)]


def css_variables() -> str:
    """Emit :root CSS custom properties from this module.

    Returns:
        A CSS string. No scientific values.

    Raises:
        Nenhum.
    """
    status_lines = []
    for name, colors in STATUS_COLORS.items():
        slug = name.lower().replace("_", "-")
        status_lines.append(f"    --hs-status-{slug}-fg: {colors['fg']};")
        status_lines.append(f"    --hs-status-{slug}-bg: {colors['bg']};")
        status_lines.append(f"    --hs-status-{slug}-border: {colors['border']};")
    joined = "\n".join(status_lines)
    return f"""
:root {{
    --hs-void: {HS_VOID};
    --hs-bg-deep: {HS_BG_DEEP};
    --hs-bg: {HS_BG};
    --hs-bg-raised: {HS_BG_RAISED};
    --hs-blue-deep: {HS_DEEP_NAVY};
    --hs-blue: {HS_ELECTRIC_BLUE};
    --hs-blue-mid: {HS_MID_BLUE};
    --hs-navy: {HS_NAVY};
    --hs-ice: {HS_ICE_BLUE};
    --hs-frost: {HS_FROST};
    --hs-white-blue: {HS_WHITE_BLUE};
    --hs-surface-1: {HS_SURFACE_1};
    --hs-surface-2: {HS_SURFACE_2};
    --hs-surface: {HS_SURFACE_1};
    --hs-surface-table: {HS_SURFACE_TABLE};
    --hs-glass: {HS_GLASS};
    --hs-glass-strong: {HS_GLASS_STRONG};
    --hs-glass-thin: {HS_GLASS_THIN};
    --hs-glass-medium: {HS_GLASS_MEDIUM};
    --hs-glass-thick: {HS_GLASS_THICK};
    --hs-glass-fallback: {HS_GLASS_FALLBACK};
    --hs-glass-border: {HS_GLASS_BORDER};
    --hs-glass-highlight: {HS_GLASS_HIGHLIGHT};
    --hs-glass-shadow: {HS_GLASS_SHADOW};
    --hs-glass-bleed: {HS_GLASS_BLEED};
    --hs-glass-active: {HS_GLASS_ACTIVE};
    --hs-primary: {HS_PRIMARY};
    --hs-cyan: {HS_CYAN};
    --hs-text: {HS_TEXT};
    --hs-text-secondary: {HS_TEXT_SECONDARY};
    --hs-text-body: {HS_TEXT_BODY};
    --hs-text-subtitle: {HS_TEXT_SUBTITLE};
    --hs-text-muted: {HS_TEXT_MUTED};
    --hs-text-meta: {HS_TEXT_META};
    --hs-border: {HS_BORDER};
    --hs-border-hover: {HS_BORDER_HOVER};
    --hs-highlight: {HS_HIGHLIGHT};
    --hs-shadow: {HS_SHADOW};
    --hs-focus: {HS_FOCUS};
    --hs-ambient: {HS_AMBIENT};
    --hs-grain: {HS_GRAIN};
    --hs-space-1: {SPACE_1}px;
    --hs-space-2: {SPACE_2}px;
    --hs-space-3: {SPACE_3}px;
    --hs-space-4: {SPACE_4}px;
    --hs-space-5: {SPACE_5}px;
    --hs-space-6: {SPACE_6}px;
    --hs-radius-sm: {RADIUS_SM}px;
    --hs-radius-md: {RADIUS_MD}px;
    --hs-radius-lg: {RADIUS_LG}px;
    --hs-motion-fast: {MOTION_FAST_MS}ms;
    --hs-motion-standard: {MOTION_STANDARD_MS}ms;
    --hs-motion-contextual: {MOTION_CONTEXTUAL_MS}ms;
    --hs-ease: {MOTION_EASING};
    --hs-font-sans: {FONT_SANS};
    --hs-font-display: {FONT_DISPLAY};
    --hs-font-mono: {FONT_MONO};
    --background-color: {HS_VOID};
    --secondary-background-color: {HS_SURFACE_1};
    --text-color: {HS_TEXT};
    --primary-color: {HS_PRIMARY};
{joined}
}}
"""
