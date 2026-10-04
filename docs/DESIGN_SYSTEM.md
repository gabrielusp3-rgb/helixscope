# HelixScope Design System (Phase 19 Part 2)

Visual language only. Tokens live in `ui/tokens.py`. Scientific numbers are never computed here.

## Layers

| Layer | Role | Material |
|---|---|---|
| 0 | Environment | Near-black navy, static radial wash, faint grain |
| 1 | Scientific content | Solid `HS_SURFACE_*` (sequences, tables, plots, results) |
| 2 | Interface | Liquid Glass on sidebar, header, inspector, command dialog |

Glass is not applied to every card. Content stays opaque for contrast.

## Palette

Defined as Python constants and CSS variables (`--hs-*`):

- `HS_BACKGROUND_DEEP` / `HS_BACKGROUND`
- `HS_SURFACE_1` / `HS_SURFACE_2`
- `HS_PRIMARY` (electric blue) / `HS_CYAN`
- `HS_TEXT` / `HS_TEXT_SECONDARY` / `HS_BORDER`

No per-page hex chaos. Semantic status colors are in `STATUS_COLORS`.

## Semantic status

Status is always text (and usually a badge). Color is secondary.

- **EXPERIMENTAL** cool silver — not “good”
- **PREDICTED** violet — not experimental
- **ILLUSTRATIVE** muted cyan — not a deposited structure
- **UNAVAILABLE** / **NOT_INSTALLED** inactive gray — not error red
- **ERROR** red
- **LIVE_VALIDATED** vs **REMOTE_VALIDATED** remain distinct families
- **RESOURCE_LIMIT** muted amber

`UNAVAILABLE` must not look like `ERROR`.

## Typography

- Display: Space Grotesk (module titles)
- Body: Inter / system-ui
- Scientific data: JetBrains Mono (sequences, accessions, hashes, HGVS, Newick)

Metric values are compact (~1.15rem), not finance-style KPIs.

## Spacing and radius

4 / 8 / 12 / 16 / 24 / 32 px. Radius 6 / 10 / 16. No one-off 7px/23px hacks in new CSS.

## Motion

FAST 120 ms, STANDARD 200 ms, CONTEXTUAL 280 ms. No continuous background animation. `prefers-reduced-motion` disables transitions.

## Focus

Visible cyan halo. `outline: none` is not used without a replacement.

## Components

- Command palette: `ui/command_palette.py` (`st.dialog` + Ctrl/Cmd+K via Components v2)
- Pointer glass: `ui/pointer.py` (Components v2, `isolate_styles=False`, trusted script only)
- Workspace HTML: `ui/workspace.py` (headers, metric strip, Compare hero, variant journey, evidence items, CRISPR chain-role legend, BLAST/MSA job state)
- Plotly theme: `ui/charts.py` `_apply_base_layout` + `ui/structure_viewer.py` (arrays unchanged)
- Chain-role legend: deposited roles only (Cas protein / guide RNA / target DNA). Unknown is not inferred.
- Job states: WAITING / RUNNING / COMPLETED / NO_HITS / ERROR as text badges. Color is secondary.
- TEST_ONLY is inactive gray, distinct from ERROR. It is never a public genome.

## Module atmosphere

Subtle ambient tint is allowed only through shared `--hs-*` tokens. Modules are not separate themes.
