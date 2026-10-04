"""HelixScope global CSS: tokens, chrome, liquid glass vs scientific content.

Visual language only. Pointer effect is Components v2 (ui.pointer).
The deprecated Components v1 HTML iframe is not used.

Glass is reserved for navigation and floating controls. Sequences, tables,
charts, and long scientific text stay on solid OLED surfaces.
"""

from __future__ import annotations

import streamlit as st

from ui.pointer import mount_pointer_effect
from ui.tokens import css_variables

_CHROME_CSS: str = r"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&family=JetBrains+Mono:wght@400;500;700&display=swap');

html, body, .stApp, [data-testid="stAppViewContainer"], [data-testid="stApp"],
section.main, section[data-testid="stMain"], .main,
[data-testid="stMainBlockContainer"] {
    background: transparent !important;
    background-color: transparent !important;
    color: var(--hs-text) !important;
    font-family: var(--hs-font-sans) !important;
}

html, body {
    background: var(--hs-void) !important;
}

/* Environment: OLED void with a single restrained specular bloom. No particles. */
.stApp::before {
    content: "";
    position: fixed;
    inset: 0;
    z-index: 0;
    pointer-events: none;
    background:
        radial-gradient(ellipse 70% 48% at 12% -8%, var(--hs-glass-bleed) 0%, transparent 58%),
        radial-gradient(ellipse 42% 32% at 88% 4%, var(--hs-highlight) 0%, transparent 52%),
        linear-gradient(180deg, var(--hs-void) 0%, var(--hs-bg-deep) 42%, var(--hs-bg) 100%);
}

.stApp::after {
    content: "";
    position: fixed;
    inset: 0;
    z-index: 0;
    pointer-events: none;
    opacity: 0.22;
    background-image:
        radial-gradient(circle at 18% 22%, var(--hs-grain, rgba(255,255,255,0.035)) 0, transparent 1.2px),
        radial-gradient(circle at 78% 68%, rgba(255,255,255,0.025) 0, transparent 1px);
    background-size: 160px 160px, 110px 110px;
}

.block-container {
    max-width: 100% !important;
    padding-top: 1.1rem !important;
    padding-bottom: 2.75rem !important;
    padding-left: 1.25rem !important;
    padding-right: 1.25rem !important;
    position: relative;
    z-index: 1;
}

#MainMenu, footer,
[data-testid="stDecoration"], [data-testid="stStatusWidget"] {
    visibility: hidden !important;
    display: none !important;
}

/* Topbar: thin floating black glass. No hard rule. Command stays rightmost via key. */
[data-testid="stHeader"] {
    background: var(--hs-glass-thin) !important;
    backdrop-filter: blur(12px) saturate(118%) !important;
    -webkit-backdrop-filter: blur(12px) saturate(118%) !important;
    border-bottom: none !important;
    box-shadow:
        inset 0 -1px 0 var(--hs-glass-highlight),
        0 8px 24px var(--hs-glass-shadow) !important;
}

.st-key-helix_cmd_open {
    display: flex !important;
    justify-content: flex-end !important;
    width: 100% !important;
}
.st-key-helix_cmd_open button {
    margin-left: auto !important;
    min-width: 7.25rem;
    border-radius: 999px !important;
    padding: 0.38rem 1.1rem !important;
}
.st-key-helix_back button,
.st-key-helix_forward button {
    min-width: 6.6rem !important;
    min-height: 2.35rem !important;
    width: 100% !important;
    padding: 0.35rem 0.7rem !important;
    font-size: 15px !important;
    font-weight: 650 !important;
    letter-spacing: 0.02em !important;
    border-radius: 999px !important;
    background: rgba(0, 0, 0, 0.72) !important;
    color: #FFFFFF !important;
    border: 1px solid rgba(255, 255, 255, 0.62) !important;
    box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.16) !important;
}
.st-key-helix_back button:disabled,
.st-key-helix_forward button:disabled {
    opacity: 1 !important;
    color: rgba(255, 255, 255, 0.55) !important;
    border-color: rgba(255, 255, 255, 0.28) !important;
    background: rgba(0, 0, 0, 0.45) !important;
}
.st-key-helix_back button:hover:not(:disabled),
.st-key-helix_forward button:hover:not(:disabled),
.st-key-helix_back button:focus-visible,
.st-key-helix_forward button:focus-visible {
    border-color: #FFFFFF !important;
    background: #FFFFFF !important;
    color: #000000 !important;
}

/* Sidebar: black glass rail flush to the left. No floating left silhouette. */
section[data-testid="stSidebar"] {
    display: flex !important;
    visibility: visible !important;
    min-width: 220px;
    background: transparent !important;
    position: relative;
    border: none !important;
    border-left: none !important;
    border-right: none !important;
    outline: none !important;
    box-shadow: none !important;
    padding-left: 0 !important;
    margin-left: 0 !important;
}
section[data-testid="stSidebar"] > div:first-child {
    position: relative;
    z-index: 1;
    isolation: isolate;
    margin: 0 6px 0 0 !important;
    border-radius: 0 var(--hs-radius-lg) var(--hs-radius-lg) 0 !important;
    background:
        linear-gradient(90deg, transparent 0%, transparent 88%, var(--hs-highlight)),
        var(--hs-glass-thick) !important;
    backdrop-filter: blur(16px) saturate(120%) !important;
    -webkit-backdrop-filter: blur(16px) saturate(120%) !important;
    border: none !important;
    border-left: none !important;
    border-right: none !important;
    outline: none !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        14px 0 32px var(--hs-glass-shadow) !important;
    height: 100vh !important;
    max-height: 100vh !important;
    min-height: 0 !important;
    display: flex !important;
    flex-direction: column !important;
    overflow: hidden !important;
}
section[data-testid="stSidebar"] > div:first-child::before {
    content: "";
    position: absolute;
    inset: 0;
    pointer-events: none;
    z-index: 1;
    border-radius: inherit;
    background:
        linear-gradient(180deg, var(--hs-glass-highlight) 0%, transparent 22%),
        radial-gradient(120% 70% at var(--mx, 30%) var(--my, 8%), var(--hs-highlight), transparent 44%);
    mix-blend-mode: screen;
    opacity: 0.8;
}
section[data-testid="stSidebar"] > div:first-child > * {
    position: relative;
    z-index: 2;
    min-height: 0 !important;
}
section[data-testid="stSidebar"] > div:first-child > *:last-child {
    flex: 1 1 auto !important;
    overflow-x: hidden !important;
    overflow-y: auto !important;
    overscroll-behavior: contain !important;
}
[data-testid="stSidebarCollapsedControl"],
[data-testid="collapsedControl"],
[data-testid="stExpandSidebarButton"],
[data-testid="stSidebarCollapseButton"] {
    display: none !important;
    visibility: hidden !important;
    pointer-events: none !important;
    width: 0 !important;
    height: 0 !important;
    overflow: hidden !important;
}
section[data-testid="stSidebar"] [data-testid="stSidebarContent"] {
    border: none !important;
    border-right: none !important;
    outline: none !important;
    box-shadow: none !important;
    background: transparent !important;
    flex: 1 1 auto !important;
    min-height: 0 !important;
    max-height: 100% !important;
    overflow-x: hidden !important;
    overflow-y: auto !important;
    overscroll-behavior: contain !important;
    scrollbar-width: thin;
    scrollbar-color: rgba(255, 255, 255, 0.45) rgba(0, 0, 0, 0);
}
section[data-testid="stSidebar"] [data-testid="stSidebarContent"]::-webkit-scrollbar {
    width: 8px;
}
section[data-testid="stSidebar"] [data-testid="stSidebarContent"]::-webkit-scrollbar-thumb {
    background: rgba(255, 255, 255, 0.45);
    border-radius: 8px;
}
section[data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
    border: none !important;
    outline: none !important;
    box-shadow: none !important;
    background: transparent !important;
    overflow: visible !important;
    height: auto !important;
    max-height: none !important;
}

[data-testid="stToolbar"] {
    display: flex !important;
    visibility: visible !important;
}

h1, h2, h3 {
    font-family: var(--hs-font-display) !important;
    text-transform: none !important;
    letter-spacing: 0.02em !important;
    font-weight: 700 !important;
    color: var(--hs-text) !important;
}
h1 { font-size: 1.55rem !important; }
h2 { font-size: 1.22rem !important; letter-spacing: 0.03em !important; }
h3 { font-size: 1.05rem !important; font-weight: 600 !important; color: var(--hs-text-subtitle) !important; }

/* Scientific content: solid, not glass-on-glass */
.glass-card, .helix-panel, .hs-content, .module-shell {
    position: relative;
    background: var(--hs-surface-1) !important;
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-md) !important;
    box-shadow: 0 8px 18px var(--hs-shadow) !important;
    padding: var(--hs-space-4) !important;
    margin-bottom: var(--hs-space-4) !important;
    backdrop-filter: none !important;
    -webkit-backdrop-filter: none !important;
    transform: none !important;
}

/* Shared glass stack: blur + black tint + inner highlight + specular. */
.hs-glass, [data-hs-glass], .hs-glass-medium {
    position: relative;
    isolation: isolate;
    background:
        linear-gradient(180deg, var(--hs-highlight), transparent 28%),
        radial-gradient(
            420px circle at var(--mx, 50%) var(--my, 0%),
            var(--hs-glass-highlight),
            transparent 46%
        ),
        var(--hs-glass-medium) !important;
    backdrop-filter: blur(14px) saturate(120%) !important;
    -webkit-backdrop-filter: blur(14px) saturate(120%) !important;
    border: 1px solid var(--hs-glass-border) !important;
    border-radius: var(--hs-radius-md) !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        0 12px 28px var(--hs-glass-shadow) !important;
}
.hs-glass-thin {
    position: relative;
    background:
        radial-gradient(
            280px circle at var(--mx, 50%) var(--my, 0%),
            var(--hs-highlight),
            transparent 48%
        ),
        var(--hs-glass-thin) !important;
    backdrop-filter: blur(10px) saturate(118%) !important;
    -webkit-backdrop-filter: blur(10px) saturate(118%) !important;
    border: 1px solid var(--hs-glass-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight) !important;
}
.hs-glass-thick {
    position: relative;
    background: var(--hs-glass-thick) !important;
    backdrop-filter: blur(16px) saturate(120%) !important;
    -webkit-backdrop-filter: blur(16px) saturate(120%) !important;
    border: 1px solid var(--hs-glass-border) !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        0 18px 40px var(--hs-glass-shadow) !important;
}

.glass-card--dna { border-color: rgba(255, 255, 255, 0.22) !important; }
.glass-card--rna { border-color: rgba(255, 255, 255, 0.16) !important; }
.glass-card--prot { border-color: rgba(255, 255, 255, 0.12) !important; }

.hs-result-header { margin-bottom: var(--hs-space-4); }
.hs-result-title-row {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: var(--hs-space-3);
    flex-wrap: wrap;
}
.hs-result-title {
    font-family: var(--hs-font-display);
    font-size: 1.12rem;
    margin: 0;
    color: var(--hs-text);
}
.hs-result-line {
    color: var(--hs-text-secondary);
    font-size: 13px;
    margin-top: 4px;
    font-family: var(--hs-font-mono);
}

.hs-metric-strip,
.hs-metric-grid {
    display: grid;
    grid-template-columns: repeat(var(--hs-metric-cols, 4), minmax(0, 1fr));
    gap: var(--hs-space-3);
    width: 100%;
    margin: var(--hs-space-3) 0;
    align-items: stretch;
}
@media (max-width: 1100px) {
    .hs-metric-strip,
    .hs-metric-grid {
        grid-template-columns: repeat(2, minmax(0, 1fr));
    }
}
.hs-metric, .metric-chip {
    background:
        linear-gradient(180deg, var(--hs-highlight) 0%, var(--hs-surface-1) 42%) !important;
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    padding: 10px 12px !important;
    min-height: 76px;
    display: flex !important;
    flex-direction: column;
    justify-content: flex-end;
    text-align: center;
    box-sizing: border-box;
    backdrop-filter: none !important;
    -webkit-backdrop-filter: none !important;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight) !important;
}
.hs-metric-value, .metric-chip .value {
    font-family: var(--hs-font-display) !important;
    font-size: 1.15rem !important;
    font-weight: 700 !important;
    color: var(--hs-text) !important;
    line-height: 1.2 !important;
    display: flex;
    align-items: baseline;
    justify-content: center;
    flex-wrap: wrap;
    gap: 4px;
    white-space: nowrap;
}
.hs-metric-number {
    font: inherit;
    color: var(--hs-text);
}
.hs-metric-unit, .metric-chip .unit {
    font-size: 12px !important;
    font-weight: 500 !important;
    color: var(--hs-text-secondary) !important;
    margin-left: 0;
    line-height: 1.2;
}
.hs-metric-label, .metric-chip .label {
    font-size: 10px !important;
    font-weight: 600 !important;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: var(--hs-text-muted) !important;
    margin-top: 6px !important;
    min-height: 2.2em;
    line-height: 1.15;
    overflow-wrap: anywhere;
}

.hs-mono, .helix-insp-v, .helix-seq-block {
    font-family: var(--hs-font-mono) !important;
}

.hs-caption { color: var(--hs-text-secondary); font-size: 12px; line-height: 1.45; }

.hs-explain {
    margin: var(--hs-space-3) 0 var(--hs-space-4) 0;
    padding: 0;
    background: var(--hs-surface-1);
    border: 1px solid var(--hs-border);
    border-radius: var(--hs-radius-sm);
    max-width: 100%;
    overflow: hidden;
    backdrop-filter: none;
}
.hs-explain-kicker {
    font-size: 10px;
    letter-spacing: 0.10em;
    text-transform: uppercase;
    color: var(--hs-text-muted);
    margin: 0;
    padding: 8px 12px 0 12px;
    background: var(--hs-glass-thin);
}
.hs-explain-title {
    font-family: var(--hs-font-display);
    font-size: 0.95rem;
    font-weight: 600;
    margin: 0;
    padding: 4px 12px 8px 12px;
    color: var(--hs-text-subtitle);
    background: var(--hs-glass-thin);
}
.hs-explain-body,
.hs-explain p {
    margin: 0;
    padding: 10px 12px 12px 12px;
    font-size: 13px;
    line-height: 1.55;
    color: var(--hs-text-body);
    background: var(--hs-surface-1);
    text-align: left;
}
.hs-engine-row {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px;
    padding: 6px 0;
    border-bottom: 1px solid var(--hs-border);
}
.hs-engine-name {
    min-width: 160px;
    font-weight: 600;
}

.hs-compare-hero { padding: var(--hs-space-4); }
.hs-compare-pair {
    display: grid;
    grid-template-columns: 1fr minmax(140px, 0.8fr) 1fr;
    gap: var(--hs-space-4);
    align-items: center;
    margin-bottom: var(--hs-space-3);
}
.hs-compare-side { text-align: center; }
.hs-compare-method { text-align: center; font-size: 13px; color: var(--hs-text-secondary); }

.hs-journey {
    display: flex;
    flex-wrap: wrap;
    align-items: stretch;
    gap: 0;
    margin: var(--hs-space-3) 0 var(--hs-space-4) 0;
}
.hs-journey-step {
    min-width: 92px;
    flex: 1 1 92px;
    background: var(--hs-surface-2);
    border: 1px solid var(--hs-border);
    border-radius: var(--hs-radius-sm);
    padding: 8px 10px;
}
.hs-journey-break { opacity: 0.72; border-style: dashed; }
.hs-journey-label {
    font-size: 10px;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: var(--hs-text-muted);
}
.hs-journey-arrow, .hs-journey-gap {
    width: 16px;
    align-self: center;
    height: 1px;
    background: var(--hs-border);
}
.hs-journey-gap { background: repeating-linear-gradient(90deg, var(--hs-border) 0 4px, transparent 4px 8px); }

.hs-evidence-item { margin-bottom: var(--hs-space-3); }
.hs-evidence-field {
    font-size: 11px;
    letter-spacing: 0.07em;
    text-transform: uppercase;
    color: var(--hs-text-muted);
}
.hs-empty-hero {
    min-height: 220px;
    display: flex;
    flex-direction: column;
    justify-content: center;
    text-align: center;
}

.hs-legend {
    display: flex;
    flex-wrap: wrap;
    gap: var(--hs-space-3);
    align-items: center;
    padding: var(--hs-space-3);
}
.hs-legend-item {
    display: flex;
    align-items: center;
    gap: 8px;
}
.hs-legend-swatch {
    width: 10px;
    height: 10px;
    border-radius: 2px;
    background: var(--hs-text-secondary);
    border: 1px solid var(--hs-border);
}
.hs-legend-cas_protein { background: #e2e8f0; }
.hs-legend-guide_rna { background: #c4b5fd; }
.hs-legend-target_dna { background: var(--hs-blue); }
.hs-legend-nontarget_dna { background: var(--hs-ice); }
.hs-legend-dna { background: var(--hs-blue); }
.hs-legend-rna { background: var(--hs-blue-mid); }
.hs-legend-protein { background: var(--hs-white-blue); }
.hs-legend-unknown { background: #64748b; }

.hs-job-state {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: var(--hs-space-3);
    padding: var(--hs-space-3);
    margin-bottom: var(--hs-space-3);
}

.hs-msa-viewer {
    font-family: var(--hs-font-mono);
    background: var(--hs-surface-table);
    border: 1px solid var(--hs-border);
    border-radius: var(--hs-radius-sm);
    padding: var(--hs-space-3);
    overflow-x: auto;
}

.helix-seq-block {
    background: var(--hs-surface-table) !important;
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    padding: var(--hs-space-4) !important;
}

.helix-meta-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 10px 14px;
    margin: 8px 0 14px 0;
    width: 100%;
}
.helix-meta-item {
    min-width: 0;
    padding: 8px 10px;
    border: 1px solid var(--hs-border);
    border-radius: var(--hs-radius-sm);
    background: var(--hs-surface-2);
}
.helix-meta-label {
    color: var(--hs-text-muted);
    font-size: 10px;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    margin-bottom: 4px;
}
.helix-meta-value {
    color: var(--hs-text-body);
    font-size: 12px;
    line-height: 1.4;
    overflow-wrap: anywhere;
    word-break: break-word;
    font-family: var(--hs-font-mono);
}
.helix-badge-row {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    align-items: center;
    margin: 6px 0 12px 0;
}

.helix-badge, .hs-status-badge {
    display: inline-flex;
    align-items: center;
    border-radius: 999px;
    padding: 2px 8px;
    font-size: 10px;
    letter-spacing: 0.06em;
    font-weight: 600;
    backdrop-filter: blur(8px) saturate(118%);
    -webkit-backdrop-filter: blur(8px) saturate(118%);
}

.helix-inspector, .hs-glass.helix-inspector {
    max-height: 640px;
    overflow-y: auto;
    overflow-x: hidden;
    padding: 10px 12px;
    z-index: 2;
    position: relative;
}
.helix-insp-section { margin-bottom: 12px; }
.helix-insp-title {
    color: var(--hs-text-subtitle);
    font-size: 11px;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    margin-bottom: 6px;
}
.helix-insp-row {
    display: grid;
    grid-template-columns: minmax(72px, 38%) minmax(0, 1fr);
    gap: 8px;
    padding: 4px 0;
    border-bottom: 1px solid var(--hs-border);
}
.helix-insp-k { color: var(--hs-text-muted); font-size: 11px; }
.helix-insp-v { color: var(--hs-text); font-size: 11px; overflow-wrap: anywhere; }

/* Navigation items: quiet until hover/active. Not a stack of neon rectangles. */
section[data-testid="stSidebar"] [data-testid="stCaptionContainer"],
section[data-testid="stSidebar"] .stCaption {
    font-size: 10px !important;
    letter-spacing: 0.14em !important;
    text-transform: uppercase !important;
    color: var(--hs-text-muted) !important;
    font-weight: 600 !important;
    opacity: 0.9;
    margin-top: 10px !important;
    margin-bottom: 4px !important;
}
section[data-testid="stSidebar"] .stButton > button {
    border-radius: var(--hs-radius-sm) !important;
    padding: 0.38rem 0.7rem !important;
    font-size: 12px !important;
    letter-spacing: 0.03em !important;
    box-shadow: none !important;
    min-height: 2rem !important;
    text-transform: none !important;
    font-weight: 500 !important;
    background: transparent !important;
    border-color: transparent !important;
    color: var(--hs-text-secondary) !important;
    backdrop-filter: none !important;
    -webkit-backdrop-filter: none !important;
    transition: background var(--hs-motion-fast) var(--hs-ease),
                color var(--hs-motion-fast) var(--hs-ease),
                border-color var(--hs-motion-fast) var(--hs-ease),
                box-shadow var(--hs-motion-fast) var(--hs-ease) !important;
}
section[data-testid="stSidebar"] .stButton > button:hover {
    background: var(--hs-glass-thin) !important;
    border-color: var(--hs-glass-border) !important;
    color: var(--hs-text) !important;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight) !important;
}
section[data-testid="stSidebar"] .stButton > button[kind="primary"] {
    background: var(--hs-glass-active) !important;
    border-color: var(--hs-border-hover) !important;
    color: var(--hs-text) !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        inset 3px 0 0 var(--hs-text) !important;
}

.helix-nav-brand {
    font-family: var(--hs-font-display) !important;
    font-size: 1.02rem !important;
    font-weight: 700 !important;
    letter-spacing: 0.16em !important;
    text-transform: uppercase !important;
    color: var(--hs-text) !important;
    margin: 4px 0 2px 0 !important;
}
.helix-nav-tag {
    color: var(--hs-text-muted) !important;
    font-size: 11px !important;
    letter-spacing: 0.04em !important;
    margin: 0 0 14px 0 !important;
}

.stTextArea textarea, .stTextInput input, .stNumberInput input,
[data-baseweb="input"], [data-baseweb="textarea"],
[data-baseweb="select"] > div, [data-baseweb="select"] > div > div {
    background-color: rgba(8, 8, 8, 0.88) !important;
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    color: var(--hs-text) !important;
    box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.55) !important;
    backdrop-filter: none !important;
}

.stTextArea textarea:focus, .stTextInput input:focus,
[data-baseweb="input"]:focus-within, [data-baseweb="textarea"]:focus-within,
[data-baseweb="select"] > div:focus-within,
.stButton > button:focus-visible {
    outline: 1px solid var(--hs-text) !important;
    outline-offset: 2px !important;
    border-color: var(--hs-border-hover) !important;
    box-shadow: 0 0 0 3px var(--hs-highlight) !important;
}

.stTextInput label, .stTextArea label, .stSelectbox label,
.stSlider label, .stCheckbox label, .stRadio label, .stFileUploader label {
    color: var(--hs-text-secondary) !important;
    font-size: 12px !important;
    font-weight: 600 !important;
    letter-spacing: 0.04em;
    text-transform: uppercase;
}

.stButton > button, .stDownloadButton > button {
    background: var(--hs-bg-raised) !important;
    border: 1px solid var(--hs-border-hover) !important;
    color: var(--hs-text) !important;
    border-radius: var(--hs-radius-sm) !important;
    font-family: var(--hs-font-display) !important;
    font-weight: 600 !important;
    letter-spacing: 0.03em !important;
    text-transform: none !important;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight) !important;
    backdrop-filter: none !important;
    -webkit-backdrop-filter: none !important;
    transition: border-color var(--hs-motion-fast) var(--hs-ease),
                background var(--hs-motion-fast) var(--hs-ease),
                box-shadow var(--hs-motion-fast) var(--hs-ease),
                color var(--hs-motion-fast) var(--hs-ease) !important;
}
.stButton > button:hover, .stDownloadButton > button:hover {
    border-color: var(--hs-border-hover) !important;
    background: var(--hs-glass-active) !important;
    transform: none !important;
    color: var(--hs-text) !important;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight) !important;
}
.stButton > button[kind="primary"] {
    background: var(--hs-glass-medium) !important;
    border-color: var(--hs-border-hover) !important;
    color: var(--hs-text) !important;
    backdrop-filter: blur(8px) saturate(118%) !important;
    -webkit-backdrop-filter: blur(8px) saturate(118%) !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        0 8px 18px var(--hs-shadow) !important;
}
.stButton > button:active, .stDownloadButton > button:active {
    background: var(--hs-text) !important;
    color: var(--hs-void) !important;
    border-color: var(--hs-text) !important;
}
section[data-testid="stSidebar"] .stButton > button:active,
section[data-testid="stSidebar"] .stButton > button[kind="primary"]:active {
    background: var(--hs-glass-active) !important;
    color: var(--hs-text) !important;
    border-color: var(--hs-border-hover) !important;
}

/* Dialog / command palette: premium floating glass. One blur surface. */
[data-testid="stDialog"],
div[role="dialog"] {
    animation: hs-materialize var(--hs-motion-contextual) var(--hs-ease);
}
[data-testid="stDialog"] > div,
div[role="dialog"] > div {
    background: var(--hs-glass-thick) !important;
    backdrop-filter: blur(20px) saturate(122%) !important;
    -webkit-backdrop-filter: blur(20px) saturate(122%) !important;
    border: 1px solid var(--hs-glass-border) !important;
    border-radius: var(--hs-radius-lg) !important;
    box-shadow:
        inset 0 1px 0 var(--hs-glass-highlight),
        0 24px 48px var(--hs-glass-shadow) !important;
}
@keyframes hs-materialize {
    from { opacity: 0; transform: translateY(8px) scale(0.985); }
    to { opacity: 1; transform: none; }
}

[data-testid="stExpander"] {
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    background: var(--hs-surface-1) !important;
    backdrop-filter: none !important;
    margin-bottom: 10px !important;
}
[data-testid="stExpander"] [data-testid="stIconMaterial"] { display: none !important; }
[data-testid="stExpander"] summary {
    display: flex !important;
    align-items: center !important;
    gap: 10px !important;
}
[data-testid="stExpander"] summary::before {
    content: "";
    display: block;
    width: 7px;
    height: 7px;
    border-right: 2px solid var(--hs-text-secondary);
    border-bottom: 2px solid var(--hs-text-secondary);
    transform: rotate(-45deg);
    flex-shrink: 0;
    transition: transform var(--hs-motion-fast) var(--hs-ease);
}
[data-testid="stExpander"] details[open] summary::before,
[data-testid="stExpander"] summary[aria-expanded="true"]::before {
    transform: rotate(45deg);
}
[data-testid="stExpander"] summary p {
    font-family: var(--hs-font-display) !important;
    font-weight: 600 !important;
    color: var(--hs-text) !important;
    margin: 0 !important;
}

[data-testid="stToggle"] { margin-bottom: 8px !important; }
[data-testid="stVerticalBlockBorderWrapper"] {
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-md) !important;
    background: var(--hs-surface-1) !important;
    backdrop-filter: none !important;
    margin-bottom: 12px !important;
    padding: 12px 16px !important;
}
[data-testid="stFileUploaderDropzone"] {
    background: var(--hs-surface-2) !important;
    border: 1px dashed var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
}
[data-testid="stDataFrame"], [data-testid="stTable"] {
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
    overflow: hidden;
    background: var(--hs-surface-table) !important;
}
[data-testid="stAlert"], .stAlert {
    background: var(--hs-surface-2) !important;
    border: 1px solid var(--hs-border) !important;
    border-radius: var(--hs-radius-sm) !important;
}
[data-testid="stSlider"] [role="slider"] {
    background: var(--hs-text) !important;
    box-shadow: none !important;
}
.stRadio label, .stCheckbox label { color: var(--hs-text) !important; }

[data-testid="stCaptionContainer"], .stCaption {
    overflow-wrap: anywhere !important;
    word-break: break-word !important;
    max-width: 100% !important;
    white-space: normal !important;
    line-height: 1.45 !important;
    color: var(--hs-text-secondary) !important;
}

.js-plotly-plot .plotly .hoverlayer { z-index: 6 !important; }
[data-testid="stPlotlyChart"],
.stPlotlyChart,
.js-plotly-plot,
iframe[title*="plotly" i],
iframe[title*="streamlit_plotly" i] {
    min-height: 480px !important;
    width: 100% !important;
    overflow: visible !important;
    transform: none !important;
    filter: none !important;
    backdrop-filter: none !important;
}

.hs-viewport {
    position: relative;
    border-radius: var(--hs-radius-lg);
    padding: 6px;
    background: var(--hs-void);
    box-shadow:
        0 0 64px var(--hs-glass-shadow),
        inset 0 1px 0 var(--hs-glass-highlight);
}
.hs-viewport::after {
    content: "";
    pointer-events: none;
    position: absolute;
    inset: 0;
    border-radius: inherit;
    box-shadow: inset 0 0 72px rgba(0, 0, 0, 0.35);
}
.hs-viewport [data-testid="stPlotlyChart"],
.hs-viewport .js-plotly-plot {
    min-height: 640px !important;
}
.js-plotly-plot .plotly,
.js-plotly-plot .svg-container,
.js-plotly-plot .gl-container { min-height: 480px !important; }

@media (max-width: 1365px) {
    .hs-compare-pair { grid-template-columns: 1fr; }
    .helix-inspector { max-height: 360px; }
    .hs-viewport [data-testid="stPlotlyChart"],
    .hs-viewport .js-plotly-plot { min-height: 420px !important; }
}

.hs-explain, .hs-explain-body, .hs-explain p {
    text-align: left;
}
.hs-metric-grid.hs-content {
    width: 100%;
    box-sizing: border-box;
}
.hs-hand-bar {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 10px 14px;
    margin: 0 0 8px 0;
    font-size: 12px;
    color: var(--hs-text-secondary);
}
.hs-hand-root {
    min-height: 0;
}
.hs-hand-hud {
    font-size: 11px;
    line-height: 1.4;
    color: var(--hs-text-secondary);
    border: 1px solid var(--hs-glass-border);
    background: var(--hs-glass-thin);
    backdrop-filter: blur(10px) saturate(118%);
    -webkit-backdrop-filter: blur(10px) saturate(118%);
    border-radius: var(--hs-radius-sm);
    padding: 8px 10px;
    margin: 0 0 8px 0;
    max-width: 28rem;
    text-align: left;
    box-shadow: inset 0 1px 0 var(--hs-glass-highlight);
}
.hs-hand-hud strong {
    color: var(--hs-text);
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    font-size: 10px;
}
.hs-hand-preview {
    display: none;
    width: 160px;
    height: 120px;
    object-fit: cover;
    border-radius: var(--hs-radius-sm);
    border: 1px solid var(--hs-border);
    margin-top: 6px;
    transform: scaleX(-1);
}
.hs-hand-preview.is-on {
    display: block;
}

@supports not ((backdrop-filter: blur(1px)) or (-webkit-backdrop-filter: blur(1px))) {
    [data-testid="stHeader"],
    section[data-testid="stSidebar"] > div:first-child,
    .hs-glass, [data-hs-glass],
    .hs-glass-thin, .hs-glass-medium, .hs-glass-thick,
    .stButton > button, .stDownloadButton > button,
    [data-testid="stDialog"] > div,
    div[role="dialog"] > div,
    .hs-hand-hud {
        background: var(--hs-glass-fallback) !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
    }
}

@media (prefers-reduced-transparency: reduce) {
    [data-testid="stHeader"],
    section[data-testid="stSidebar"] > div:first-child,
    .hs-glass, [data-hs-glass],
    .hs-glass-thin, .hs-glass-medium, .hs-glass-thick,
    .stButton > button, .stDownloadButton > button,
    [data-testid="stDialog"] > div,
    .hs-hand-hud {
        background: var(--hs-glass-fallback) !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
    }
}

@media (prefers-reduced-motion: reduce) {
    *, *::before, *::after {
        animation: none !important;
        transition: none !important;
    }
    [data-testid="stAppViewContainer"],
    [data-testid="stSidebar"],
    [data-testid="stHeader"],
    [data-hs-glass],
    .hs-glass,
    .hs-glass-thin,
    .hs-glass-medium,
    .hs-glass-thick,
    .stButton > button,
    .hs-hand-hud,
    [data-testid="stDialog"] > div {
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
    }
}
</style>
"""


def inject_custom_css() -> None:
    """Inject tokenized HelixScope CSS and the pointer listener.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    sheet = "<style>" + css_variables() + "</style>" + _CHROME_CSS
    st.html(sheet)
    st.markdown(sheet, unsafe_allow_html=True)
    mount_pointer_effect()
