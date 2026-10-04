"""Componentes visuais HTML do HelixScope.

Cada funcao retorna uma string HTML que usa as classes CSS definidas em
ui/styles.py (glass-card, metric-chip, badge) combinadas a estilos inline para
garantir aparencia consistente. Nenhuma string contem emojis.
"""

from __future__ import annotations

import html
import re
from typing import Dict, Mapping, Optional, Sequence

from ui.tokens import HS_FROST, HS_ICE_BLUE, HS_MID_BLUE, HS_SURFACE_TABLE, HS_TEXT, status_style

MAX_SEQUENCE_DISPLAY_RESIDUES: int = 12_000
"""Limite do visualizador HTML. As metricas continuam usando a sequencia completa."""


def html_escape(text: object) -> str:
    """Escapa texto para interpolacao em HTML gerado pela interface.

    Args:
        text: Valor a exibir; None vira string vazia.

    Returns:
        String segura para inserir em atributos e nos nos de texto HTML.

    Raises:
        Nenhum.
    """
    return html.escape("" if text is None else str(text), quote=True)


def safe_download_filename(stem: str, fallback: str = "download") -> str:
    """Remove caracteres inseguros de um nome de arquivo de download.

    Args:
        stem: Nome proposto (accession, rotulo do usuario, etc.).
        fallback: Nome usado quando o stem nao restar utilizavel.

    Returns:
        Identificador com letras, digitos, ponto, hifen e sublinhado, no maximo
        80 caracteres, sem travessia de diretorio.

    Raises:
        Nenhum.
    """
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(stem or ""))
    cleaned = cleaned.strip("._-")
    if len(cleaned) > 80:
        cleaned = cleaned[:80].rstrip("._-")
    return cleaned or fallback


CARD_VARIANTS: Dict[str, str] = {
    "dna": "glass-card--dna",
    "rna": "glass-card--rna",
    "prot": "glass-card--prot",
}
"""Mapeia o nome da variante de cartao para a classe modificadora correspondente."""

BADGE_STYLES: Dict[str, Dict[str, str]] = {
    "dna": {"bg": "rgba(145,180,228,0.12)", "color": HS_ICE_BLUE, "border": "rgba(145,180,228,0.35)"},
    "rna": {"bg": "rgba(104,147,208,0.12)", "color": HS_MID_BLUE, "border": "rgba(104,147,208,0.35)"},
    "prot": {"bg": "rgba(175,203,241,0.12)", "color": HS_FROST, "border": "rgba(175,203,241,0.35)"},
    "warn": {"bg": "rgba(251,191,36,0.12)", "color": "#fbbf24", "border": "rgba(251,191,36,0.35)"},
    "err": {"bg": "rgba(248,113,113,0.12)", "color": "#F87171", "border": "rgba(248,113,113,0.35)"},
    "computed": {"bg": "rgba(255,255,255,0.08)", "color": "#E0E0E0", "border": "rgba(255,255,255,0.28)"},
    "retrieved": {"bg": "rgba(255,255,255,0.07)", "color": "#DADADA", "border": "rgba(255,255,255,0.26)"},
    "predicted": {"bg": "rgba(251,191,36,0.12)", "color": "#fbbf24", "border": "rgba(251,191,36,0.35)"},
    "unavailable": {"bg": "rgba(111,111,111,0.14)", "color": "#8A8A8A", "border": "rgba(111,111,111,0.35)"},
}
"""Estilos de cor (fundo a 8%, texto solido e borda sutil) por variante de badge."""

DNA_COLORS: Dict[str, str] = {
    "A": HS_ICE_BLUE, "T": "#F472B6", "G": "#34D399", "C": "#FBBF24", "N": "#64748b",
}
"""Cores por base para exibicao de sequencias de DNA."""

RNA_COLORS: Dict[str, str] = {
    "A": HS_ICE_BLUE, "U": "#818CF8", "G": "#34D399", "C": "#FBBF24",
}
"""Cores por base para exibicao de sequencias de RNA."""

PROTEIN_GROUP_COLORS: Dict[str, str] = {
    "hydrophobic": "#f97316",
    "polar": HS_ICE_BLUE,
    "positive": "#F87171",
    "negative": "#818CF8",
}
"""Cores por grupo fisico-quimico de aminoacidos."""

PROTEIN_GROUPS: Dict[str, str] = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic", "I": "hydrophobic",
    "M": "hydrophobic", "F": "hydrophobic", "W": "hydrophobic", "P": "hydrophobic",
    "G": "hydrophobic",
    "S": "polar", "T": "polar", "C": "polar", "Y": "polar", "N": "polar", "Q": "polar",
    "H": "positive", "K": "positive", "R": "positive",
    "D": "negative", "E": "negative",
}
"""Classificacao de cada aminoacido padrao em grupo fisico-quimico."""

DEFAULT_RESIDUE_COLOR: str = HS_TEXT
"""Cor usada para residuos sem classificacao definida."""

SEQUENCE_BACKGROUND: str = HS_SURFACE_TABLE
"""Cor de fundo do bloco de exibicao de sequencias (superficie cientifica solida)."""


def glass_card(title: str, content: str, variant: str = "default") -> str:
    """Monta um cartao glassmorphism com titulo e conteudo HTML.

    Args:
        title: Titulo do cartao, exibido em caixa alta e cor suave.
        content: Conteudo do corpo, inserido como HTML interno sem escape.
        variant: Variante visual: "default", "dna", "rna" ou "prot".

    Returns:
        String HTML de um elemento div com a classe glass-card e, quando
        aplicavel, a classe modificadora da variante.

    Raises:
        Nenhum.

    Nota:
        O titulo e escapado. O corpo continua sendo HTML montado pela aplicacao;
        campos vindos do usuario ou do NCBI devem ser passados por html_escape
        antes de entrar em content.
    """
    modifier = CARD_VARIANTS.get(variant, "")
    classes = "glass-card" + (f" {modifier}" if modifier else "")
    title_style = (
        "font-family:'Space Grotesk',sans-serif;font-size:13px;text-transform:uppercase;"
        "letter-spacing:0.08em;color:var(--hs-text-muted);margin-bottom:10px;font-weight:700;"
    )
    return (
        f'<div class="{classes}">'
        f'<div style="{title_style}">{html_escape(title)}</div>'
        f'<div>{content}</div>'
        f"</div>"
    )


def metric_chip(label: str, value: str, unit: str = "") -> str:
    """Monta um chip de metrica com valor em destaque, unidade e rotulo.

    Args:
        label: Rotulo da metrica, exibido em caixa alta abaixo do valor.
        value: Valor principal, exibido em destaque.
        unit: Unidade opcional, exibida em linha apos o valor.

    Returns:
        String HTML de um elemento div com a classe metric-chip e os estilos
        inline especificados.

    Raises:
        Nenhum.

    Nota:
        Use para exibir grandezas calculadas (por exemplo, conteudo GC ou massa
        molecular) de forma compacta e legivel.
    """
    container_style = (
        "background:var(--hs-surface-2);border:1px solid var(--hs-border);"
        "border-radius:6px;padding:10px 12px;"
    )
    value_style = (
        "font-family:var(--hs-font-display);font-size:1.15rem;font-weight:700;"
        "line-height:1.2;color:var(--hs-text);display:flex;align-items:baseline;"
        "justify-content:center;gap:4px;white-space:nowrap;"
    )
    unit_style = (
        "font-family:var(--hs-font-sans);font-size:12px;font-weight:500;color:var(--hs-text-secondary);"
    )
    label_style = (
        "font-family:var(--hs-font-sans);font-size:10px;font-weight:600;color:var(--hs-text-secondary);"
        "text-transform:uppercase;letter-spacing:0.08em;margin-top:6px;min-height:2.2em;"
    )
    unit_html = (
        f'<span class="unit hs-metric-unit" style="{unit_style}">{html_escape(unit)}</span>'
        if unit
        else ""
    )
    return (
        f'<div class="metric-chip" style="{container_style}">'
        f'<div class="value hs-metric-value" style="{value_style}">'
        f'<span class="hs-metric-number">{html_escape(value)}</span>{unit_html}</div>'
        f'<div class="label" style="{label_style}">{html_escape(label)}</div>'
        f"</div>"
    )


def badge(text: str, variant: str) -> str:
    """Monta um badge em formato de pilula com cor por variante.

    Args:
        text: Texto exibido dentro do badge.
        variant: Variante de cor: "dna", "rna", "prot", "warn", "err",
            "computed", "retrieved", "predicted" ou "unavailable".

    Returns:
        String HTML de um elemento span estilizado como pilula.

    Raises:
        ValueError: Se a variante nao for uma das suportadas.

    Nota:
        As variantes "warn" e "err" comunicam alertas e erros; as demais
        identificam o tipo de molecula associada.
    """
    if variant not in BADGE_STYLES:
        suportadas = ", ".join(sorted(BADGE_STYLES))
        raise ValueError(f"Variante de badge invalida. Use uma de: {suportadas}.")
    style = BADGE_STYLES[variant]
    badge_style = (
        f"background:{style['bg']};color:{style['color']};"
        f"border:1px solid {style['border']};border-radius:999px;"
        "padding:3px 10px;font-size:11px;font-weight:700;display:inline-block;"
    )
    return f'<span class="badge" style="{badge_style}">{html_escape(text)}</span>'


STATUS_BADGE_VARIANTS: Dict[str, str] = {
    "COMPUTED": "computed",
    "RETRIEVED": "retrieved",
    "PREDICTED": "predicted",
    "VALIDATED": "prot",
    "WARNING": "warn",
    "ERROR": "err",
    "UNAVAILABLE": "unavailable",
    "NO_HITS": "warn",
    "WAITING": "predicted",
    "READY": "retrieved",
    "TIMEOUT": "err",
    "RATE_LIMITED": "warn",
    "SERVICE_UNAVAILABLE": "err",
    "INVALID_INPUT": "err",
    "INVALID_DATABASE": "err",
    "JOB_FAILED": "err",
    "PARSING_ERROR": "err",
    "RESOURCE_LIMIT": "err",
    "FAILED": "err",
    "CACHED": "retrieved",
    "QUEUED": "predicted",
    "RUNNING": "predicted",
    "COMPLETED": "computed",
    "CANCELLED": "warn",
    "EXPERIMENTAL": "retrieved",
    "NOT_FOUND": "err",
    "NETWORK_ERROR": "err",
    "INSUFFICIENT_DATA": "err",
    "NO_STRUCTURE": "unavailable",
    "TOOL_NOT_INSTALLED": "unavailable",
    "TOOL_FAILED": "err",
    "SEQUENCE_UNAVAILABLE": "unavailable",
    "ILLUSTRATIVE": "predicted",
    "EXACT": "retrieved",
    "ALIGNED": "predicted",
    "PARTIAL": "warn",
    "BEST_EFFORT": "warn",
    "UNMAPPED": "unavailable",
    "UNCERTAIN": "warn",
    "STALE": "warn",
    "REMOTE_VALIDATED": "retrieved",
    "LIVE_VALIDATED": "retrieved",
    "TEST_ONLY": "warn",
    "PARTIAL_VIEW": "warn",
    "SCOPE_MISMATCH": "warn",
    "NOT_RUN": "unavailable",
    "INCOMPLETE": "warn",
    "VERIFIED": "computed",
    "AVAILABLE": "retrieved",
    "SKIPPED": "unavailable",
    "NOT_COMPUTED": "unavailable",
    "MAPPED": "predicted",
    "NORMALIZED": "computed",
    "IDENTIFIER_ONLY": "warn",
    "REF_VALIDATED": "computed",
    "REF_NOT_CHECKED": "warn",
    "REFERENCE_MISMATCH": "err",
    "UNSUPPORTED": "unavailable",
}


def status_badge(status: str) -> str:
    """Monta um badge semantico com cor tokenizada e texto obrigatorio.

    Args:
        status: Token cientifico (COMPUTED, EXPERIMENTAL, PREDICTED, ...).

    Returns:
        HTML do badge. A cor nao substitui o texto: EXPERIMENTAL e PREDICTED
        permanecem distinguiveis pelo rotulo.

    Raises:
        Nenhum.
    """
    key = str(status or "UNAVAILABLE").strip().upper()
    style = status_style(key)
    label = key.replace("_", " ")
    badge_style = (
        f"background:{style['bg']};color:{style['fg']};"
        f"border:1px solid {style['border']};border-radius:999px;"
        "padding:3px 8px;font-size:10px;font-weight:700;display:inline-block;"
        "letter-spacing:0.06em;"
    )
    return (
        f'<span class="badge hs-status-badge" style="{badge_style}" '
        f'aria-label="{html_escape(label)}">{html_escape(label)}</span>'
    )


def _residue_color(residue: str, seq_type: str) -> str:
    """Determina a cor de um residuo conforme o tipo de sequencia (uso interno).

    Args:
        residue: Caractere do residuo em maiuscula.
        seq_type: Tipo de sequencia normalizado ("DNA", "RNA" ou "PROT").

    Returns:
        Codigo hexadecimal de cor para o residuo.

    Nota biologica:
        A coloracao por tipo destaca padroes de composicao: bases para acidos
        nucleicos e grupos fisico-quimicos para aminoacidos.
    """
    if seq_type == "DNA":
        return DNA_COLORS.get(residue, DEFAULT_RESIDUE_COLOR)
    if seq_type == "RNA":
        return RNA_COLORS.get(residue, DEFAULT_RESIDUE_COLOR)
    group = PROTEIN_GROUPS.get(residue)
    if group is None:
        return DEFAULT_RESIDUE_COLOR
    return PROTEIN_GROUP_COLORS[group]


def sequence_display(
    seq: str,
    seq_type: str,
    line_width: int = 60,
    highlight_indices: Optional[Sequence[int]] = None,
    view_start: int = 0,
    view_end: Optional[int] = None,
) -> str:
    """Monta um bloco HTML colorido para exibir uma sequencia em linhas.

    Quebra a sequencia em linhas de line_width caracteres, numera a posicao
    inicial de cada linha a esquerda e colore cada residuo conforme o tipo de
    sequencia.

    Args:
        seq: Sequencia de residuos (sera normalizada para maiusculas).
        seq_type: Tipo da sequencia: "DNA", "RNA" ou "PROT" (sem distincao de
            caixa).
        line_width: Numero de residuos por linha; deve ser positivo.
        highlight_indices: Indices 0-based a destacar. Visual apenas.
        view_start: Inicio 0-based da janela do viewer.
        view_end: Fim exclusivo; None = min(view_start+MAX, length).

    Returns:
        String HTML de um bloco monoespacado com numeracao e coloracao por
        residuo. A visualizacao e limitada a MAX_SEQUENCE_DISPLAY_RESIDUES
        por janela; isso nao altera os calculos feitos sobre a sequencia completa.

    Raises:
        ValueError: Se a sequencia for vazia, line_width nao for positivo ou
            seq_type nao for "DNA", "RNA" nem "PROT".

    Nota biologica:
        A visualizacao posicionada e colorida facilita inspecionar composicao,
        motivos e regioes de interesse ao longo da sequencia.
    """
    cleaned = "".join(seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia esta vazia.")
    if line_width <= 0:
        raise ValueError("line_width deve ser positivo.")
    normalized_type = seq_type.strip().upper()
    if normalized_type not in {"DNA", "RNA", "PROT"}:
        raise ValueError("seq_type deve ser 'DNA', 'RNA' ou 'PROT'.")

    total = len(cleaned)
    try:
        start = max(0, int(view_start or 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("view_start must be an integer.") from exc
    if start >= total:
        raise ValueError(f"view_start {start} is outside 0..{total}.")
    if view_end is None:
        stop = min(total, start + MAX_SEQUENCE_DISPLAY_RESIDUES)
    else:
        try:
            stop = int(view_end)
        except (TypeError, ValueError) as exc:
            raise ValueError("view_end must be an integer.") from exc
        if stop <= start:
            raise ValueError(f"Region [{start}, {stop}) is empty.")
        stop = min(total, stop)
    if stop - start > MAX_SEQUENCE_DISPLAY_RESIDUES:
        stop = start + MAX_SEQUENCE_DISPLAY_RESIDUES
    window = cleaned[start:stop]
    truncated_note = ""
    if start != 0 or stop != total:
        truncated_note = (
            f'<div style="color:var(--hs-text-secondary);font-size:12px;margin-bottom:8px;'
            f'overflow-wrap:anywhere;">'
            f"Showing residues {start:,}–{stop:,} of {total:,} "
            f"(0-based, end exclusive). Metrics use the full sequence.</div>"
        )
    highlighted = {int(index) for index in (highlight_indices or [])}
    number_width = len(str(total))
    lines_html = []
    for offset in range(0, len(window), line_width):
        block = window[offset : offset + line_width]
        abs_offset = start + offset
        position = abs_offset + 1
        residues_html = []
        for local_index, residue in enumerate(block):
            abs_index = abs_offset + local_index
            extra = ""
            if abs_index in highlighted:
                extra = (
                    "background:rgba(255,255,255,0.12);outline:1px solid rgba(255,255,255,0.45);"
                    "border-radius:2px;"
                )
            residues_html.append(
                f'<span style="color:{_residue_color(residue, normalized_type)};{extra}">'
                f"{html_escape(residue)}</span>"
            )
        residues_joined = "".join(residues_html)
        number_html = (
            f'<span style="color:var(--hs-text-meta);">{str(position).rjust(number_width)}</span>'
        )
        lines_html.append(f"{number_html}  {residues_joined}")

    body = "\n".join(lines_html)
    container_style = (
        f"background:{SEQUENCE_BACKGROUND};padding:16px;border-radius:12px;"
        "border:1px solid rgba(255,255,255,0.06);"
        "overflow-x:auto;font-family:'JetBrains Mono','Fira Code',monospace;"
        "font-size:13px;line-height:1.6;white-space:pre;"
    )
    return f'<div class="helix-seq-viewer" style="{container_style}">{truncated_note}{body}</div>'


def meta_grid(rows: Sequence[tuple]) -> str:
    """Grelha de metadados com wrapping. Evita uma unica linha com pipes.

    Args:
        rows: Pares (label, value).

    Returns:
        HTML de definition list.

    Raises:
        Nenhum.
    """
    cells = []
    for label, value in rows:
        text = "N/A" if value in (None, "") else str(value)
        cells.append(
            "<div class='helix-meta-item'>"
            f"<div class='helix-meta-label'>{html_escape(label)}</div>"
            f"<div class='helix-meta-value' title='{html_escape(text)}'>"
            f"{html_escape(text)}</div></div>"
        )
    return f"<div class='helix-meta-grid'>{''.join(cells)}</div>"


def badge_row(statuses: Sequence[str]) -> str:
    """Linha de badges com wrap. Espaco independente por estado.

    Args:
        statuses: Lista de status_badge keys.

    Returns:
        HTML.

    Raises:
        Nenhum.
    """
    parts = [status_badge(str(item)) for item in statuses if str(item or "").strip()]
    return f"<div class='helix-badge-row'>{''.join(parts)}</div>"


def inspector_panel(sections: Mapping[str, Sequence[tuple]]) -> str:
    """Painel de inspector seccionado. Nao interpreta ciencia extra.

    Args:
        sections: Nome da secao -> lista de (campo, valor).

    Returns:
        HTML scrollavel.

    Raises:
        Nenhum.
    """
    blocks = []
    for title, rows in sections.items():
        items = []
        for label, value in rows:
            text = "N/A" if value in (None, "") else str(value)
            items.append(
                f"<div class='helix-insp-row'><span class='helix-insp-k'>"
                f"{html_escape(label)}</span><span class='helix-insp-v' title='"
                f"{html_escape(text)}'>{html_escape(text)}</span></div>"
            )
        blocks.append(
            f"<div class='helix-insp-section'><div class='helix-insp-title'>"
            f"{html_escape(title)}</div>{''.join(items)}</div>"
        )
    return f"<div class='helix-inspector hs-glass' data-hs-glass='1'>{''.join(blocks)}</div>"

