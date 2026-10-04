"""Engine grafico do HelixScope baseado em Plotly (serie 6.x).

Todas as funcoes retornam plotly.graph_objects.Figure com um layout base
comum (fundo translucido, fonte Inter, margens e gridlines suaves). Nenhuma
string de label ou hover contem emojis.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from helixscope_core.compat import scale_profile

from ui.tokens import (
    FONT_SANS,
    HS_VOID,
    PLOTLY_FONT,
    PLOTLY_GRID,
    PLOTLY_PAPER,
    PLOTLY_PLOT,
)

GRID_COLOR: str = PLOTLY_GRID
"""Cor quase imperceptivel das gridlines em eixos cartesianos e polares."""

NUCLEOTIDE_COLORS: Dict[str, str] = {
    "A": "rgba(255,255,255,1)",
    "T": "rgba(255,255,255,0.72)",
    "U": "rgba(255,255,255,0.72)",
    "G": "rgba(255,255,255,0.50)",
    "C": "rgba(255,255,255,0.34)",
    "N": "rgba(255,255,255,0.18)",
}
"""Luminancia por base em graficos 2D. A identidade da base esta no rotulo."""

STANDARD_AMINO_ACIDS: List[str] = list("ACDEFGHIKLMNPQRSTVWY")
"""Os 20 aminoacidos padrao em ordem alfabetica."""

CODON_COLORSCALE = [[0.0, "#000000"], [1.0, "#FFFFFF"]]
"""Luminancia do heatmap de codons. O valor numerico permanece no hover."""

DINUCLEOTIDE_COLORSCALE = [[0.0, "#000000"], [1.0, "#FFFFFF"]]
"""Luminancia do heatmap de dinucleotideos. O valor numerico permanece no hover."""

AMINO_ACID_COLORS: Dict[str, str] = {
    "G": "rgba(255,255,255,0.92)", "A": "rgba(255,255,255,0.92)",
    "P": "rgba(255,255,255,0.92)", "V": "rgba(255,255,255,0.92)",
    "L": "rgba(255,255,255,0.92)", "I": "rgba(255,255,255,0.92)",
    "M": "rgba(255,255,255,0.92)",
    "F": "rgba(255,255,255,0.70)", "Y": "rgba(255,255,255,0.70)",
    "W": "rgba(255,255,255,0.70)",
    "S": "rgba(255,255,255,0.52)", "T": "rgba(255,255,255,0.52)",
    "C": "rgba(255,255,255,0.52)", "N": "rgba(255,255,255,0.52)",
    "Q": "rgba(255,255,255,0.52)",
    "K": "rgba(255,255,255,0.36)", "R": "rgba(255,255,255,0.36)",
    "H": "rgba(255,255,255,0.36)",
    "D": "rgba(255,255,255,0.24)", "E": "rgba(255,255,255,0.24)",
    "*": "rgba(255,255,255,0.12)", "X": "rgba(255,255,255,0.12)",
}
"""Luminancia por classe de cadeia lateral em graficos 2D. O rotulo carrega a identidade."""

CATEGORY_COLORS: Dict[str, str] = {
    "nonpolar_aliphatic": "rgba(255,255,255,0.92)",
    "aromatic": "rgba(255,255,255,0.70)",
    "polar_uncharged": "rgba(255,255,255,0.52)",
    "positively_charged": "rgba(255,255,255,0.36)",
    "negatively_charged": "rgba(255,255,255,0.24)",
    "charged_total": "rgba(255,255,255,0.30)",
    "hydrophobic_kyte_doolittle": "rgba(255,255,255,0.80)",
}
"""Cor por categoria de aminoacido usada no grafico de categorias."""

CATEGORY_LABELS: Dict[str, str] = {
    "nonpolar_aliphatic": "Nonpolar aliphatic",
    "aromatic": "Aromatic",
    "polar_uncharged": "Polar uncharged",
    "positively_charged": "Positively charged",
    "negatively_charged": "Negatively charged",
    "charged_total": "Charged (total)",
    "hydrophobic_kyte_doolittle": "Hydrophobic (Kyte-Doolittle)",
}
"""Rotulo legivel de cada categoria de aminoacido."""

DOTPLOT_COLORSCALE = [[0.0, HS_VOID], [1.0, "#FFFFFF"]]
"""Escala de cores binaria do dotplot de alinhamento."""


def _apply_base_layout(fig: go.Figure, cartesian: bool = True) -> go.Figure:
    """Aplica o layout base obrigatorio a uma figura (uso interno).

    Args:
        fig: Figura Plotly a estilizar.
        cartesian: Se verdadeiro, aplica gridlines aos eixos x e y cartesianos.

    Returns:
        A mesma figura, com fundo, fonte, margens e gridlines padronizados.

    Nota:
        Centraliza o tema visual exigido para manter coerencia entre todos os
        graficos do projeto.
    """
    fig.update_layout(
        paper_bgcolor=PLOTLY_PAPER,
        plot_bgcolor=PLOTLY_PLOT,
        font=dict(family=FONT_SANS, color=PLOTLY_FONT),
        margin=dict(l=40, r=20, t=40, b=40),
        hoverlabel=dict(
            bgcolor="rgba(0,0,0,0.94)",
            bordercolor="rgba(255,255,255,0.28)",
            font=dict(color="#FFFFFF", size=12),
            align="left",
        ),
    )
    if cartesian:
        fig.update_xaxes(gridcolor=GRID_COLOR, zeroline=False)
        fig.update_yaxes(gridcolor=GRID_COLOR, zeroline=False)
    return fig


def _count_and_frequency(value: object) -> tuple[float, float]:
    """Extrai contagem e frequencia de um valor de composicao (uso interno).

    Args:
        value: Valor associado a um simbolo: um dict com "count" e "frequency"
            ou um numero interpretado como contagem.

    Returns:
        Tupla (count, frequency); frequency e 0.0 quando nao disponivel.

    Nota:
        Permite aceitar tanto a saida estruturada de nucleotide_composition
        quanto um simples mapeamento de base para contagem.
    """
    if isinstance(value, dict):
        return float(value.get("count", 0)), float(value.get("frequency", 0.0))
    return float(value), 0.0


def nucleotide_bar_chart(composition: dict) -> go.Figure:
    """Cria um grafico de barras horizontais da composicao de nucleotideos.

    Args:
        composition: Mapeamento de base para contagem, ou para um dict com
            "count" e "frequency" (saida de nucleotide_composition).

    Returns:
        Figura Plotly de barras horizontais coloridas por base, com hover de
        contagem e percentual.

    Raises:
        ValueError: Se nenhuma base tiver contagem positiva (Insufficient data).

    Nota biologica:
        A composicao de bases revela vies de GC e a abundancia relativa de cada
        nucleotideo na sequencia.
    """
    bases: List[str] = []
    counts: List[float] = []
    freqs: List[float] = []
    colors: List[str] = []
    total = 0.0
    for base, value in composition.items():
        count, freq = _count_and_frequency(value)
        total += count
        bases.append(base)
        counts.append(count)
        freqs.append(freq)
    if total > 0 and all(f == 0.0 for f in freqs):
        freqs = [round((c / total) * 100.0, 2) for c in counts]

    filtered = [(b, c, f) for b, c, f in zip(bases, counts, freqs) if c > 0]
    if not filtered:
        raise ValueError("Insufficient data")
    bases = [b for b, _, _ in filtered]
    counts = [c for _, c, _ in filtered]
    freqs = [f for _, _, f in filtered]
    colors = ["rgba(255,255,255,0.92)"] * len(bases)
    patterns = [
        {"A": "", "C": "/", "G": ".", "T": "x", "U": "-", "N": "+"}.get(base, "")
        for base in bases
    ]

    fig = go.Figure(
        go.Bar(
            x=counts,
            y=bases,
            orientation="h",
            marker=dict(color=colors, pattern=dict(shape=patterns)),
            customdata=freqs,
            hovertemplate="Base %{y}: %{x:.0f} (%{customdata:.2f}%)<extra></extra>",
        )
    )
    fig.update_layout(xaxis_title="Count", yaxis_title="Base")
    return _apply_base_layout(fig)


def gc_gauge(gc_percent: float) -> go.Figure:
    """Cria um medidor (gauge) do conteudo GC.

    Args:
        gc_percent: Conteudo GC em porcentagem, no intervalo [0, 100].

    Returns:
        Figura Plotly com um go.Indicator em modo gauge+number, faixas de cor
        para GC baixo, medio e alto, e sufixo de porcentagem.

    Raises:
        ValueError: Se gc_percent for NaN (Insufficient data) ou estiver fora
            do intervalo [0, 100].

    Nota biologica:
        O conteudo GC influencia a estabilidade termica do DNA e o desenho de
        oligonucleotideos; faixas extremas dificultam amplificacao e hibridacao.
    """
    if isinstance(gc_percent, float) and math.isnan(gc_percent):
        raise ValueError("Insufficient data")
    if not 0.0 <= gc_percent <= 100.0:
        raise ValueError("gc_percent deve estar entre 0 e 100.")
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=gc_percent,
            number={"suffix": "%"},
            gauge={
                "axis": {"range": [0, 100]},
                "bar": {"color": "#FFFFFF"},
                "steps": [
                    {"range": [0, 40], "color": "rgba(255,255,255,0.16)"},
                    {"range": [40, 60], "color": "rgba(255,255,255,0.34)"},
                    {"range": [60, 100], "color": "rgba(255,255,255,0.55)"},
                ],
            },
        )
    )
    return _apply_base_layout(fig, cartesian=False)


def gc_skew_plot(
    skew_values: List[float], window: int, y_title: str = "GC Skew"
) -> go.Figure:
    """Cria um grafico de GC skew com area preenchida por sinal.

    Args:
        skew_values: Valores de GC skew por janela (podem conter NaN).
        window: Tamanho da janela usada, para converter indices em posicao
            genomica (posicao = window * indice).

    Returns:
        Figura Plotly com a linha de GC skew. A area positiva e mais clara e a
        negativa e mais escura; a legenda nomeia o sinal. A linha de referencia
        fica em y igual a zero.

    Raises:
        Nenhum.

    Nota biologica:
        O GC skew ajuda a localizar origens e terminos de replicacao, pois as
        fitas lider e tardia acumulam vieses opostos de G e C.
    """
    x = [window * i for i in range(len(skew_values))]
    clean = [None if (v is None or (isinstance(v, float) and math.isnan(v))) else v for v in skew_values]
    positive = []
    negative = []
    for value in clean:
        if value is None:
            positive.append(None)
            negative.append(None)
        else:
            positive.append(value if value > 0 else 0)
            negative.append(value if value < 0 else 0)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x, y=positive, mode="lines", line=dict(width=0),
            fill="tozeroy", fillcolor="rgba(255,255,255,0.28)",
            name="Positive skew",
            hoverinfo="skip", showlegend=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x, y=negative, mode="lines", line=dict(width=0),
            fill="tozeroy", fillcolor="rgba(255,255,255,0.08)",
            name="Negative skew",
            hoverinfo="skip", showlegend=True,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x, y=clean, mode="lines", line=dict(color="#FFFFFF", width=1.5, dash="dot"),
            hovertemplate="Position %{x} bp: %{y:.4f}<extra></extra>",
            showlegend=False,
        )
    )
    fig.add_hline(y=0, line_color="rgba(255,255,255,0.2)")
    fig.update_layout(xaxis_title="Genomic Position (bp)", yaxis_title=y_title)
    return _apply_base_layout(fig)


def codon_usage_heatmap(
    codon_df: pd.DataFrame, value_column: str = "Synonymous_pct"
) -> go.Figure:
    """Cria um heatmap de uso de codons por aminoacido.

    Args:
        codon_df: DataFrame com as colunas "Codon", "AminoAcid" e a coluna de
            valores indicada (saida de codon_usage_table).
        value_column: Nome da coluna numerica a colorir; por padrao a porcentagem
            dentro da familia sinonima.

    Returns:
        Figura Plotly com um go.Heatmap de codons (eixo x) por aminoacidos
        (eixo y), colorido pelo valor escolhido.

    Raises:
        ValueError: Se faltarem colunas obrigatorias no DataFrame.

    Nota biologica:
        O mapa evidencia o vies de uso de codons sinonimos, relevante para
        eficiencia de traducao e desenho de genes sinteticos.
    """
    required = {"Codon", "AminoAcid", value_column}
    if not required.issubset(codon_df.columns):
        faltando = ", ".join(sorted(required - set(codon_df.columns)))
        raise ValueError(f"Colunas ausentes no DataFrame: {faltando}.")

    pivot = codon_df.pivot_table(
        index="AminoAcid", columns="Codon", values=value_column, aggfunc="sum"
    )
    fig = go.Figure(
        go.Heatmap(
            z=pivot.values,
            x=list(pivot.columns),
            y=list(pivot.index),
            colorscale=CODON_COLORSCALE,
            hovertemplate=(
                "Codon: %{x}<br>Amino acid: %{y}<br>Frequency: %{z:.2f}%<extra></extra>"
            ),
        )
    )
    fig.update_layout(xaxis_title="Codon", yaxis_title="Amino acid")
    return _apply_base_layout(fig)


def codon_usage_bar_chart(
    codon_df: pd.DataFrame, value_column: str = "Absolute_pct"
) -> go.Figure:
    """Cria um grafico de barras de uso de codons agrupado por aminoacido.

    Args:
        codon_df: DataFrame com as colunas "Codon", "AminoAcid" e a coluna de
            valores indicada (saida de codon_usage_table).
        value_column: Nome da coluna numerica a plotar no eixo y.

    Returns:
        Figura Plotly de barras verticais, uma por codon, rotuladas como
        "aminoacido/codon" e coloridas pela natureza da cadeia lateral do
        aminoacido.

    Raises:
        ValueError: Se faltarem colunas obrigatorias no DataFrame.

    Nota biologica:
        Agrupar os codons pelo aminoacido que especificam permite comparar
        diretamente os sinonimos entre si, que e a comparacao biologicamente
        relevante: a escolha entre eles nao altera a proteina, apenas a velocidade
        e a fidelidade da traducao.
    """
    required = {"Codon", "AminoAcid", value_column}
    if not required.issubset(codon_df.columns):
        faltando = ", ".join(sorted(required - set(codon_df.columns)))
        raise ValueError(f"Colunas ausentes no DataFrame: {faltando}.")

    ordered = codon_df.sort_values(by=["AminoAcid", "Codon"]).reset_index(drop=True)
    labels = [
        f"{row.AminoAcid}/{row.Codon}" for row in ordered.itertuples(index=False)
    ]
    values = ordered[value_column].tolist()
    colors = [AMINO_ACID_COLORS.get(a, "#64748B") for a in ordered["AminoAcid"]]

    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color=colors,
            customdata=ordered["Count"].tolist() if "Count" in ordered else values,
            hovertemplate=(
                "%{x}<br>" + value_column + ": %{y:.2f}<br>"
                "Count: %{customdata}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        xaxis_title="Amino acid / Codon",
        yaxis_title=value_column.replace("_", " "),
        bargap=0.15,
    )
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=9))
    return _apply_base_layout(fig)


def rscu_bar_chart(rscu_df: pd.DataFrame) -> go.Figure:
    """Cria um grafico de barras de RSCU com linha de referencia em 1.0.

    Args:
        rscu_df: DataFrame com as colunas "Codon", "AminoAcid" e "RSCU" (saida de
            relative_synonymous_codon_usage).

    Returns:
        Figura Plotly de barras verticais por codon, coloridas em ciano quando o
        RSCU esta acima de 1.0 (codon preferido) e em ambar quando esta abaixo
        (codon evitado), com linha tracejada em y igual a 1.0.

    Raises:
        ValueError: Se faltarem colunas obrigatorias no DataFrame.

    Nota biologica:
        A linha em 1.0 marca a ausencia de vies, isto e, o uso de um codon na
        proporcao esperada entre seus sinonimos. Barras muito acima de 1.0
        identificam os codons otimos do organismo, que costumam corresponder aos
        tRNAs mais abundantes.
    """
    required = {"Codon", "AminoAcid", "RSCU"}
    if not required.issubset(rscu_df.columns):
        faltando = ", ".join(sorted(required - set(rscu_df.columns)))
        raise ValueError(f"Colunas ausentes no DataFrame: {faltando}.")

    ordered = rscu_df.sort_values(by=["AminoAcid", "Codon"]).reset_index(drop=True)
    labels = [
        f"{row.AminoAcid}/{row.Codon}" for row in ordered.itertuples(index=False)
    ]
    values = [float(v) for v in ordered["RSCU"]]
    colors = ["#6893D0" if v >= 1.0 else "#FBBF24" for v in values]

    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color=colors,
            hovertemplate="%{x}<br>RSCU: %{y:.3f}<extra></extra>",
        )
    )
    fig.add_hline(
        y=1.0, line_color="rgba(255,255,255,0.35)", line_dash="dash"
    )
    fig.update_layout(
        xaxis_title="Amino acid / Codon", yaxis_title="RSCU", bargap=0.15
    )
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=9))
    return _apply_base_layout(fig)


def dinucleotide_heatmap(frequencies: dict) -> go.Figure:
    """Cria um heatmap 4x4 da razao observado/esperado dos dinucleotideos.

    Args:
        frequencies: Saida de dna_analysis.dinucleotide_frequencies(), mapeando
            cada dinucleotideo para um dict com "count", "frequency" e
            "observed_expected".

    Returns:
        Figura Plotly com um go.Heatmap divergente centrado em 1.0, com a primeira
        base no eixo y e a segunda no eixo x, anotando o valor de cada celula.

    Raises:
        ValueError: Se frequencies estiver vazio.

    Nota biologica:
        Celulas proximas de 1.0 indicam bases adjacentes independentes. O desvio
        classico e a deplecao de CpG em vertebrados, com razao proxima de 0.2,
        consequencia da desaminacao de citosinas metiladas nesse contexto.
    """
    if not frequencies:
        raise ValueError("frequencies esta vazio.")

    bases = list("ACGT")
    matrix: List[List[float]] = []
    text: List[List[str]] = []
    for first in bases:
        row_values: List[float] = []
        row_text: List[str] = []
        for second in bases:
            entry = frequencies.get(first + second, {})
            odds = entry.get("observed_expected") if isinstance(entry, dict) else None
            if odds is None or (isinstance(odds, float) and math.isnan(odds)):
                row_values.append(None)
                row_text.append("N/A")
            else:
                number = float(odds)
                row_values.append(number)
                row_text.append(f"{number:.2f}")
        matrix.append(row_values)
        text.append(row_text)

    fig = go.Figure(
        go.Heatmap(
            z=matrix,
            x=bases,
            y=bases,
            text=text,
            texttemplate="%{text}",
            textfont=dict(size=12, color="#E2E8F0"),
            colorscale=DINUCLEOTIDE_COLORSCALE,
            zmid=1.0,
            hovertemplate=(
                "Dinucleotide %{y}%{x}<br>Observed/Expected: %{z:.3f}<extra></extra>"
            ),
        )
    )
    fig.update_layout(xaxis_title="Second base", yaxis_title="First base")
    return _apply_base_layout(fig)


def gc_sliding_window_plot(profile: List[dict]) -> go.Figure:
    """Cria o grafico de linha do conteudo GC em janela deslizante.

    Args:
        profile: Saida de dna_analysis.gc_sliding_window(), lista de dicts com
            "midpoint" e "gc_percent".

    Returns:
        Figura Plotly com a linha de GC por janela, area preenchida e linha
        tracejada na media do perfil.

    Raises:
        ValueError: Se profile estiver vazio.

    Nota biologica:
        Desvios locais acentuados em relacao a media sinalizam ilhas CpG,
        isocoros, regioes adquiridas por transferencia horizontal ou fronteiras
        entre exons e introns.
    """
    if not profile:
        raise ValueError("profile esta vazio.")

    x = [int(point["midpoint"]) for point in profile]
    y: List[Optional[float]] = []
    for point in profile:
        value = float(point["gc_percent"])
        y.append(None if math.isnan(value) else value)
    defined = [value for value in y if value is not None]
    if not defined:
        raise ValueError("Insufficient data")
    mean_gc = sum(defined) / len(defined)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color="#6893D0", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(56,189,248,0.12)",
            hovertemplate="Position %{x} bp<br>GC %{y:.2f}%<extra></extra>",
            showlegend=False,
        )
    )
    fig.add_hline(
        y=mean_gc,
        line_color="rgba(255,255,255,0.35)",
        line_dash="dash",
        annotation_text=f"Mean {mean_gc:.2f}%",
        annotation_position="top left",
    )
    fig.update_layout(
        xaxis_title="Genomic Position (bp)",
        yaxis_title="GC Content (%)",
        yaxis_range=[0, 100],
    )
    return _apply_base_layout(fig)


def orf_length_histogram(lengths: List[int]) -> go.Figure:
    """Cria o histograma da distribuicao de tamanho das regioes codificadoras.

    Args:
        lengths: Comprimentos em nucleotideos das ORFs ou CDS selecionados.

    Returns:
        Figura Plotly com um go.Histogram dos comprimentos e linha tracejada na
        mediana.

    Raises:
        ValueError: Se lengths estiver vazio.

    Nota biologica:
        Genes reais concentram-se em uma faixa de comprimento caracteristica do
        organismo, enquanto ORFs geradas por acaso produzem uma cauda exponencial
        dominada por quadros curtos. Uma distribuicao com muitos quadros no limite
        inferior do filtro sugere que parte das ORFs e artefato estatistico.
    """
    if not lengths:
        raise ValueError("lengths esta vazio.")

    values = [int(v) for v in lengths]
    ordered = sorted(values)
    middle = len(ordered) // 2
    median = (
        float(ordered[middle])
        if len(ordered) % 2 == 1
        else (ordered[middle - 1] + ordered[middle]) / 2.0
    )

    fig = go.Figure(
        go.Histogram(
            x=values,
            marker_color="#6893D0",
            marker_line=dict(color="rgba(255,255,255,0.15)", width=1),
            hovertemplate="Length %{x} nt<br>Count: %{y}<extra></extra>",
        )
    )
    fig.add_vline(
        x=median,
        line_color="rgba(255,255,255,0.35)",
        line_dash="dash",
        annotation_text=f"Median {median:.0f} nt",
        annotation_position="top right",
    )
    fig.update_layout(
        xaxis_title="Coding Region Length (nt)",
        yaxis_title="Number of Regions",
        bargap=0.05,
    )
    return _apply_base_layout(fig)


def amino_acid_bar_chart(composition: dict) -> go.Figure:
    """Cria um grafico de barras da composicao de aminoacidos.

    Args:
        composition: Mapeamento de aminoacido para contagem, ou para um dict com
            "count" e "frequency" (saida de amino_acid_composition).

    Returns:
        Figura Plotly de barras verticais sobre os 20 aminoacidos padrao, coloridas
        pela natureza da cadeia lateral, com hover de contagem e percentual.

    Raises:
        Nenhum.

    Nota biologica:
        A leitura em barras permite comparar magnitudes absolutas entre residuos,
        algo que o grafico de radar dificulta; desvios marcantes na composicao
        apontam vies funcional, como enriquecimento de cisteina em proteinas
        secretadas ou de prolina em regioes desordenadas.
    """
    aminos: List[str] = []
    counts: List[float] = []
    freqs: List[float] = []
    total = 0.0
    for amino in STANDARD_AMINO_ACIDS:
        count, freq = _count_and_frequency(composition.get(amino, 0))
        aminos.append(amino)
        counts.append(count)
        freqs.append(freq)
        total += count
    if total > 0 and all(f == 0.0 for f in freqs):
        freqs = [round((c / total) * 100.0, 2) for c in counts]

    colors = [AMINO_ACID_COLORS.get(a, "#64748B") for a in aminos]

    fig = go.Figure(
        go.Bar(
            x=aminos,
            y=counts,
            marker_color=colors,
            customdata=freqs,
            hovertemplate=(
                "Amino acid %{x}: %{y:.0f} (%{customdata:.2f}%)<extra></extra>"
            ),
        )
    )
    fig.update_layout(xaxis_title="Amino acid", yaxis_title="Count", bargap=0.2)
    return _apply_base_layout(fig)


def amino_acid_category_chart(categories: dict) -> go.Figure:
    """Cria um grafico de barras horizontais das categorias de aminoacidos.

    Args:
        categories: Saida de protein_analysis.amino_acid_categories(), mapeando
            nome de categoria para um dict com "count", "frequency" e "residues".

    Returns:
        Figura Plotly de barras horizontais com a frequencia percentual de cada
        categoria presente no dicionario.

    Raises:
        ValueError: Se categories estiver vazio.

    Nota biologica:
        As cinco primeiras categorias sao mutuamente exclusivas e somam 100 por
        cento; "Charged (total)" e "Hydrophobic (Kyte-Doolittle)" sao derivadas e
        se sobrepoem as anteriores, por isso nao devem ser somadas as demais.
    """
    if not categories:
        raise ValueError("categories esta vazio.")

    names: List[str] = []
    values: List[float] = []
    counts: List[float] = []
    colors: List[str] = []
    for key, entry in categories.items():
        count, freq = _count_and_frequency(entry)
        names.append(CATEGORY_LABELS.get(key, key.replace("_", " ")))
        values.append(freq)
        counts.append(count)
        colors.append(CATEGORY_COLORS.get(key, "#94A3B8"))

    fig = go.Figure(
        go.Bar(
            x=values,
            y=names,
            orientation="h",
            marker_color=colors,
            customdata=counts,
            hovertemplate="%{y}: %{x:.2f}% (%{customdata:.0f} residues)<extra></extra>",
        )
    )
    fig.update_layout(xaxis_title="Frequency (%)", yaxis_title="Category")
    return _apply_base_layout(fig)


def amino_acid_radar(composition: dict) -> go.Figure:
    """Cria um grafico de radar da composicao de aminoacidos.

    Args:
        composition: Mapeamento de aminoacido para contagem, ou para um dict com
            "count" e "frequency" (saida de amino_acid_composition).

    Returns:
        Figura Plotly com um go.Scatterpolar sobre os 20 aminoacidos padrao em
        ordem alfabetica, com preenchimento e eixo radial ajustado ao maximo.

    Raises:
        Nenhum.

    Nota biologica:
        O perfil de aminoacidos resume propriedades globais da proteina, como
        predominancia de residuos carregados, polares ou hidrofobicos.
    """
    values: List[float] = []
    for amino in STANDARD_AMINO_ACIDS:
        count, _ = _count_and_frequency(composition.get(amino, 0))
        values.append(count)
    max_value = max(values) if values else 0.0

    theta = STANDARD_AMINO_ACIDS + [STANDARD_AMINO_ACIDS[0]]
    radial = values + [values[0]]

    fig = go.Figure(
        go.Scatterpolar(
            r=radial,
            theta=theta,
            fill="toself",
            fillcolor="rgba(56,189,248,0.1)",
            line=dict(color="#6893D0"),
            hovertemplate="Amino acid %{theta}: %{r:.0f}<extra></extra>",
        )
    )
    fig.update_layout(
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(range=[0, max_value * 1.1 if max_value > 0 else 1], gridcolor=GRID_COLOR),
            angularaxis=dict(gridcolor=GRID_COLOR),
        ),
    )
    return _apply_base_layout(fig, cartesian=False)


def hydrophobicity_plot(profile: List[float], protein_seq: str) -> go.Figure:
    """Cria o grafico do perfil de hidrofobicidade de Kyte-Doolittle.

    Args:
        profile: Valores medios de hidropatia por janela.
        protein_seq: Sequencia de proteina usada para anotar o residuo no hover.

    Returns:
        Figura Plotly com linha suavizada (spline), areas positiva (laranja) e
        negativa (ciano) e linha de referencia em y igual a zero.

    Raises:
        Nenhum.

    Nota biologica:
        Picos positivos prolongados sugerem regioes transmembrana, enquanto
        vales indicam segmentos hidrofilicos expostos ao solvente.
    """
    protein = "".join(protein_seq.split()).upper()
    x = list(range(1, len(profile) + 1))
    residues = [protein[i] if i < len(protein) else "" for i in range(len(profile))]
    positive = [v if v > 0 else 0 for v in profile]
    negative = [v if v < 0 else 0 for v in profile]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x, y=positive, mode="lines", line=dict(width=0, shape="spline"),
            fill="tozeroy", fillcolor="rgba(249,115,22,0.2)",
            hoverinfo="skip", showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x, y=negative, mode="lines", line=dict(width=0, shape="spline"),
            fill="tozeroy", fillcolor="rgba(56,189,248,0.2)",
            hoverinfo="skip", showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=x, y=profile, mode="lines", line=dict(color="#E2E8F0", width=1.5, shape="spline"),
            customdata=residues,
            hovertemplate=(
                "Position %{x}<br>Residue %{customdata}<br>"
                "Kyte-Doolittle %{y:.3f}<extra></extra>"
            ),
            showlegend=False,
        )
    )
    fig.add_hline(y=0, line_color="rgba(255,255,255,0.2)")
    fig.update_layout(xaxis_title="Window Position", yaxis_title="Hydropathy")
    return _apply_base_layout(fig)


GUIDE_RISK_COLORS: Dict[str, str] = {
    "Low": "#34D399",
    "Moderate": "#FBBF24",
    "High": "#F87171",
    "Unknown": "#64748B",
}
"""Cor por nivel de risco global de off-target de um guia."""


def guide_ranking_scatter(guides: List[dict]) -> go.Figure:
    """Cria o grafico de dispersao eficiencia x especificidade dos guias.

    Args:
        guides: Lista de guias enriquecidos por crispr.evaluate_guides, com
            "doench_score", "specificity_proxy" e, opcionalmente, "risk" e
            "guide_sequence".

    Returns:
        Figura Plotly com eficiencia heuristica no eixo x, proxy local de
        especificidade no eixo y e cor por nivel de risco. Apenas guias com proxy
        de especificidade calculado sao plotados.

    Raises:
        ValueError: Se nenhum guia possuir proxy de especificidade.

    Nota biologica:
        O canto superior direito reune os guias mais desejaveis: alta eficiencia
        prevista e alta especificidade local. Como os eixos usam metricas
        heuristicas e um proxy local, a leitura e comparativa, nao absoluta.
    """
    plotted = [g for g in guides if g.get("specificity_proxy") is not None]
    if not plotted:
        raise ValueError(
            "Nenhum guia possui proxy de especificidade; rode a analise de "
            "off-target para habilitar este grafico."
        )

    x = [float(g["doench_score"]) for g in plotted]
    y = [float(g["specificity_proxy"]) for g in plotted]
    colors = [
        GUIDE_RISK_COLORS.get((g.get("risk") or {}).get("level", "Unknown"), "#64748B")
        for g in plotted
    ]
    labels = [g.get("guide_sequence", "") for g in plotted]
    ranks = [g.get("rank", 0) for g in plotted]

    fig = go.Figure(
        go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(size=12, color=colors, line=dict(width=1, color="#0B1220")),
            customdata=list(zip(labels, ranks)),
            hovertemplate=(
                "Rank %{customdata[1]}<br>%{customdata[0]}<br>"
                "Efficiency %{x:.3f}<br>Specificity proxy %{y:.3f}<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        xaxis_title="Heuristic efficiency (0 to 1)",
        yaxis_title="Local specificity proxy (0 to 1)",
        xaxis_range=[0, 1],
        yaxis_range=[0, 1.02],
    )
    return _apply_base_layout(fig)


def off_target_mismatch_histogram(summary: Dict[int, int]) -> go.Figure:
    """Cria o grafico de barras da distribuicao de off-targets por mismatches.

    Args:
        summary: Dicionario de crispr.off_target_mismatch_summary, mapeando numero
            de mismatches para contagem de sitios.

    Returns:
        Figura Plotly de barras verticais, uma por numero de mismatches, com a
        contagem de sitios locais.

    Raises:
        ValueError: Se summary estiver vazio.

    Nota biologica:
        Sitios com menos mismatches, sobretudo na regiao proxima ao PAM, sao os de
        maior risco de clivagem; a barra em zero e um mismatch merece a maior
        atencao.
    """
    if not summary:
        raise ValueError("summary esta vazio.")

    keys = sorted(summary.keys())
    values = [summary[k] for k in keys]
    colors = ["#F87171" if k <= 1 else "#FBBF24" if k == 2 else "#6893D0" for k in keys]

    fig = go.Figure(
        go.Bar(
            x=[str(k) for k in keys],
            y=values,
            marker_color=colors,
            hovertemplate="%{x} mismatch(es): %{y} sites<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis_title="Number of mismatches",
        yaxis_title="Verified sites under selected method",
        bargap=0.25,
    )
    return _apply_base_layout(fig)


def off_target_contig_histogram(hits: List[dict]) -> go.Figure:
    """Barras de hits verificados por contig/cromossomo da referencia fornecida.

    Args:
        hits: Lista de off-targets com chromosome_or_contig.

    Returns:
        Figura Plotly. Cada barra e uma contagem real de hits da busca.

    Raises:
        ValueError: Se hits estiver vazio.

    Nota biologica:
        A distribuicao descreve os contigs do FASTA carregado, nao um cariotipo
        humano ilustrativo. Sem hits reais o grafico nao e desenhado.
    """
    if not hits:
        raise ValueError("No verified off-target hits to plot.")
    counts: Dict[str, int] = {}
    for hit in hits:
        name = str(hit.get("chromosome_or_contig") or "").strip() or "unnamed"
        counts[name] = counts.get(name, 0) + 1
    labels = sorted(counts.keys())
    values = [counts[name] for name in labels]
    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color="#6893D0",
            hovertemplate="%{x}: %{y} verified hits<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis_title="Contig in provided reference",
        yaxis_title="Verified hits",
        bargap=0.25,
    )
    return _apply_base_layout(fig)


def off_target_pam_histogram(counts: Dict[str, int]) -> go.Figure:
    """Barras de PAMs realmente observados nos hits verificados.

    Args:
        counts: Mapa classe/sequencia de PAM -> contagem.

    Returns:
        Figura Plotly. Sem dados reais nao e chamada.

    Raises:
        ValueError: Se counts estiver vazio.
    """
    if not counts:
        raise ValueError("No verified PAM counts to plot.")
    labels = sorted(counts.keys())
    values = [counts[name] for name in labels]
    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color="#34D399",
            hovertemplate="%{x}: %{y} verified hits<extra></extra>",
        )
    )
    fig.update_layout(
        xaxis_title="PAM class in verified hits",
        yaxis_title="Verified hits",
        bargap=0.25,
    )
    return _apply_base_layout(fig)


def off_target_score_histogram(
    scores: List[float],
    *,
    title: str,
    x_title: str,
) -> go.Figure:
    """Histograma de scores realmente calculados nos hits verificados.

    Args:
        scores: Lista de scores numericos (Hsu 0-100 ou CFD 0-1).
        title: Titulo que declara o metodo.
        x_title: Eixo X com a escala.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Se scores estiver vazio.
    """
    if not scores:
        raise ValueError("No verified scores to plot.")
    fig = go.Figure(
        go.Histogram(
            x=scores,
            marker_color="#6893D0",
            hovertemplate="count=%{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title="Verified hits",
        bargap=0.05,
    )
    return _apply_base_layout(fig)


def cut_site_map(sequence_length: int, guide: dict, cas_system: str) -> go.Figure:
    """Desenha a posicao do protoespacador, do PAM e do sitio de corte.

    Args:
        sequence_length: Comprimento total da sequencia alvo, em pares de base.
        guide: Dicionario de guia com "position", "strand", "guide_sequence",
            "pam_sequence" e, opcionalmente, "cut_site".
        cas_system: Nome do sistema Cas, usado apenas no titulo do hover.

    Returns:
        Figura Plotly com um eixo horizontal representando a sequencia, uma barra
        para o protoespacador, um marcador para o PAM e uma linha vertical no
        sitio de corte quando aplicavel.

    Raises:
        ValueError: Se sequence_length nao for positivo.

    Nota biologica:
        Visualizar onde a nuclease corta em relacao a regiao alvo ajuda a prever
        qual parte do gene sera afetada e a posicionar primers de validacao.
    """
    if sequence_length <= 0:
        raise ValueError("sequence_length deve ser positivo.")

    start = int(guide.get("position", 1))
    glen = len(str(guide.get("guide_sequence", "")))
    end = start + glen
    strand = guide.get("strand", "+")
    cut = guide.get("cut_site")

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[1, sequence_length],
            y=[0, 0],
            mode="lines",
            line=dict(color="rgba(148,163,184,0.5)", width=2),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[start, end],
            y=[0, 0],
            mode="lines",
            line=dict(color="#6893D0", width=12),
            hovertemplate=(
                f"Protospacer ({cas_system})<br>Strand {strand}<br>"
                f"{start}-{end}<extra></extra>"
            ),
            showlegend=False,
        )
    )
    pam_seq = str(guide.get("pam_sequence", ""))
    if pam_seq:
        fig.add_trace(
            go.Scatter(
                x=[end],
                y=[0],
                mode="markers",
                marker=dict(size=12, color="#FBBF24", symbol="square"),
                hovertemplate=f"PAM {pam_seq}<extra></extra>",
                showlegend=False,
            )
        )
    if cut is not None:
        fig.add_vline(
            x=int(cut),
            line_color="#F87171",
            line_dash="dash",
            annotation_text=f"Cut {int(cut)}",
            annotation_position="top",
        )
    fig.update_layout(
        xaxis_title="Position in target sequence (bp)",
        yaxis=dict(visible=False, range=[-1, 1]),
        xaxis_range=[0, sequence_length + 1],
        height=220,
    )
    return _apply_base_layout(fig)


def guide_comparison_chart(guides: List[dict]) -> go.Figure:
    """Compara os melhores guias em eficiencia, especificidade e GC.

    Args:
        guides: Lista de guias enriquecidos por crispr.evaluate_guides, ja
            ordenados; apenas os primeiros sao comparados.

    Returns:
        Figura Plotly de barras agrupadas com series normalizadas a [0, 1]:
        eficiencia heuristica, proxy de especificidade (somente quando calculado)
        e conteudo GC dividido por 100. Especificidade ausente nao e plotada como
        zero.

    Raises:
        ValueError: Se a lista de guias estiver vazia.

    Nota biologica:
        Comparar lado a lado evita escolher um guia apenas pela eficiencia quando
        outro, quase tao eficiente, e bem mais especifico.
    """
    if not guides:
        raise ValueError("A lista de guias esta vazia.")

    labels = [f"#{g.get('rank', i + 1)} {g.get('guide_sequence', '')[:8]}" for i, g in enumerate(guides)]
    efficiency = [float(g.get("doench_score", 0.0)) for g in guides]
    specificity = [
        float(g["specificity_proxy"]) if g.get("specificity_proxy") is not None else None
        for g in guides
    ]
    gc_scaled = [float(g.get("gc_content", 0.0)) / 100.0 for g in guides]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            name="Heuristic efficiency",
            x=labels,
            y=efficiency,
            marker_color="#6893D0",
        )
    )
    if any(value is not None for value in specificity):
        fig.add_trace(
            go.Bar(
                name="Local specificity proxy",
                x=labels,
                y=specificity,
                marker_color="#34D399",
            )
        )
    fig.add_trace(
        go.Bar(name="GC / 100", x=labels, y=gc_scaled, marker_color="#A78BFA")
    )
    fig.update_layout(
        barmode="group",
        xaxis_title="Guide",
        yaxis_title="Normalised value (0 to 1)",
        yaxis_range=[0, 1.02],
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return _apply_base_layout(fig)


def alignment_dotplot(matrix: np.ndarray, label1: str, label2: str) -> go.Figure:
    """Cria um dotplot de alinhamento a partir de uma matriz binaria.

    Args:
        matrix: Matriz binaria (0/1) de identidade entre janelas das sequencias.
        label1: Rotulo da primeira sequencia (linhas, eixo y).
        label2: Rotulo da segunda sequencia (colunas, eixo x).

    Returns:
        Figura Plotly com um go.Heatmap binario sem barra de cores.

    Raises:
        Nenhum.

    Nota biologica:
        O dotplot evidencia regioes de similaridade, repeticoes e inversoes
        entre duas sequencias.
    """
    fig = go.Figure(
        go.Heatmap(
            z=np.asarray(matrix),
            colorscale=DOTPLOT_COLORSCALE,
            showscale=False,
            zmin=0,
            zmax=1,
            hovertemplate=f"{label2}: %{{x}}<br>{label1}: %{{y}}<extra></extra>",
        )
    )
    fig.update_layout(xaxis_title=label2, yaxis_title=label1)
    return _apply_base_layout(fig)


def entropy_profile_plot(profile: List[dict]) -> go.Figure:
    """Grafico de entropia de Shannon em janela deslizante.

    Args:
        profile: Saida de dna_analysis.entropy_sliding_window(), com
            "midpoint" e "entropy".

    Returns:
        Figura Plotly da entropia ao longo da sequencia.

    Raises:
        ValueError: Se profile estiver vazio.
    """
    if not profile:
        raise ValueError("Insufficient data")
    x = [int(point["midpoint"]) for point in profile]
    y = [float(point["entropy"]) for point in profile]
    fig = go.Figure(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color="#818CF8", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(129,140,248,0.12)",
            hovertemplate="Position %{x} bp<br>Entropy %{y:.4f} bits<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(
        xaxis_title="Position (bp)",
        yaxis_title="Shannon entropy (bits)",
        yaxis_range=[0, 2.05],
    )
    fig.update_layout(title="Entropy profile (computed from loaded sequence)")
    return _apply_base_layout(fig)


def kmer_bar_chart(counts: dict, top_n: int = 30) -> go.Figure:
    """Grafico dos k-mers mais abundantes.

    Args:
        counts: Saida de kmer_counts (k-mer -> count/frequency).
        top_n: Maximo de barras; deve ser positivo.

    Returns:
        Figura Plotly de barras horizontais.

    Raises:
        ValueError: Se counts estiver vazio ou top_n nao for positivo.
    """
    if not counts:
        raise ValueError("Insufficient data")
    if top_n <= 0:
        raise ValueError("top_n deve ser positivo.")
    ranked = sorted(
        counts.items(),
        key=lambda item: int(item[1]["count"]),
        reverse=True,
    )[:top_n]
    labels = [name for name, _ in ranked][::-1]
    values = [int(item["count"]) for _, item in ranked][::-1]
    freqs = [float(item["frequency"]) for _, item in ranked][::-1]
    fig = go.Figure(
        go.Bar(
            x=values,
            y=labels,
            orientation="h",
            marker_color="#6893D0",
            customdata=freqs,
            hovertemplate="k-mer %{y}: %{x:.0f} (%{customdata:.4f}%)<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Most abundant k-mers (top {len(labels)})",
        xaxis_title="Count",
        yaxis_title="k-mer",
    )
    return _apply_base_layout(fig)


def charge_profile_plot(profile: List[float], protein_seq: str) -> go.Figure:
    """Grafico do perfil de carga formal de cadeias laterais.

    Args:
        profile: Saida de protein_analysis.charge_count_profile().
        protein_seq: Sequencia usada no hover.

    Returns:
        Figura Plotly do perfil de carga formal.

    Raises:
        ValueError: Se profile estiver vazio.
    """
    if not profile:
        raise ValueError("Insufficient data")
    protein = "".join(protein_seq.split()).upper()
    x = list(range(1, len(profile) + 1))
    residues = [protein[i] if i < len(protein) else "" for i in range(len(profile))]
    fig = go.Figure(
        go.Scatter(
            x=x,
            y=profile,
            mode="lines",
            line=dict(color="#34D399", width=1.5),
            customdata=residues,
            hovertemplate=(
                "Window %{x}<br>Residue %{customdata}<br>"
                "Formal charge %{y:.3f}<extra></extra>"
            ),
            showlegend=False,
        )
    )
    fig.add_hline(y=0, line_color="rgba(255,255,255,0.2)")
    fig.update_layout(
        title="Side-chain formal charge profile (K,R minus D,E)",
        xaxis_title="Window position",
        yaxis_title="Mean formal charge",
    )
    return _apply_base_layout(fig)


def sequence_feature_map(
    sequence_length: int, tracks: List[dict]
) -> go.Figure:
    """Mapa horizontal de features com coordenadas reais.

    Args:
        sequence_length: Comprimento da sequencia de referencia.
        tracks: Lista de trilhas {"name": str, "features": [{"start", "end",
            "label"}]}, coordenadas 0-based semiabertas na fita de entrada.

    Returns:
        Figura Plotly com segmentos por feature.

    Raises:
        ValueError: Se o comprimento for invalido ou nao houver features.
    """
    if sequence_length <= 0:
        raise ValueError("sequence_length deve ser positivo.")
    features_found = False
    fig = go.Figure()
    colors = ["#6893D0", "#818CF8", "#34D399", "#FBBF24"]
    for index, track in enumerate(tracks):
        name = str(track.get("name") or f"track {index + 1}")
        color = colors[index % len(colors)]
        for feature in track.get("features") or []:
            start = int(feature["start"])
            end = int(feature["end"])
            if end <= start:
                continue
            features_found = True
            label = str(feature.get("label") or name)
            fig.add_trace(
                go.Scatter(
                    x=[start, end],
                    y=[name, name],
                    mode="lines",
                    line=dict(color=color, width=14),
                    hovertemplate=(
                        f"{label}<br>start {start}<br>end {end}"
                        "<extra></extra>"
                    ),
                    showlegend=False,
                )
            )
    if not features_found:
        raise ValueError("Insufficient data")
    fig.update_layout(
        title="Sequence map (coordinates from computed features)",
        xaxis_title="Position on input strand (0-based)",
        yaxis_title="Track",
        xaxis_range=[0, sequence_length],
    )
    return _apply_base_layout(fig)


def alignment_column_map(classes: List[str]) -> go.Figure:
    """Mapa 1D das colunas do alinhamento: match, mismatch ou gap.

    Args:
        classes: Lista gerada por alignment.classify_alignment_columns.

    Returns:
        Figura Plotly com uma barra por coluna.

    Raises:
        ValueError: Se classes estiver vazia ou contiver categoria desconhecida.
    """
    if not classes:
        raise ValueError("Insufficient data")
    allowed = {"match", "mismatch", "gap"}
    for item in classes:
        if item not in allowed:
            raise ValueError("Insufficient data")
    color_map = {
        "match": "#34D399",
        "mismatch": "#F87171",
        "gap": "#64748B",
    }
    fig = go.Figure()
    for category in ("match", "mismatch", "gap"):
        xs = [i for i, value in enumerate(classes) if value == category]
        if not xs:
            continue
        fig.add_trace(
            go.Bar(
                x=xs,
                y=[1] * len(xs),
                name=category,
                marker_color=color_map[category],
                hovertemplate=f"column %{{x}}<br>{category}<extra></extra>",
            )
        )
    if not fig.data:
        raise ValueError("Insufficient data")
    fig.update_layout(
        title="Alignment column map (match / mismatch / gap)",
        xaxis_title="Alignment column (0-based)",
        yaxis_title="",
        yaxis_visible=False,
        barmode="overlay",
        bargap=0,
        showlegend=True,
    )
    return _apply_base_layout(fig)


def _blast_hit_labels(hits: List[dict]) -> List[str]:
    """Rotulos curtos de accession NCBI para eixos de grafico BLAST."""
    labels: List[str] = []
    for index, hit in enumerate(hits, start=1):
        acc = str(hit.get("accession") or hit.get("hit_id") or f"hit_{index}")
        labels.append(acc[:40])
    return labels


def blast_hit_ranking_chart(hits: List[dict]) -> go.Figure:
    """Barras horizontais de bit score NCBI, na ordem dos hits.

    Args:
        hits: Hits ja validados do XML BlastOutput.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Se hits estiver vazio (Insufficient data).
    """
    if not hits:
        raise ValueError("Insufficient data")
    labels = _blast_hit_labels(hits)
    scores = [float(hit["bit_score"]) for hit in hits]
    fig = go.Figure(
        go.Bar(
            x=scores[::-1],
            y=labels[::-1],
            orientation="h",
            marker_color="#6893D0",
            hovertemplate="Hit %{y}<br>Bit score %{x}<extra></extra>",
        )
    )
    fig.update_layout(
        title="BLAST hit ranking (NCBI bit score)",
        xaxis_title="Bit score (NCBI)",
        yaxis_title="Subject accession",
    )
    return _apply_base_layout(fig)


def blast_evalue_chart(hits: List[dict]) -> go.Figure:
    """E-values NCBI por hit. Eixo log apenas quando todos os valores sao > 0.

    Args:
        hits: Hits ja validados.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data.
    """
    if not hits:
        raise ValueError("Insufficient data")
    labels = _blast_hit_labels(hits)
    values = [float(hit["evalue"]) for hit in hits]
    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color="#818CF8",
            hovertemplate="Hit %{x}<br>E-value %{y}<extra></extra>",
        )
    )
    if all(value > 0 for value in values):
        y_title = "E-value (NCBI, log scale)"
        fig.update_yaxes(type="log")
    else:
        y_title = (
            "E-value (NCBI; 0 is underflow in the XML, not a computed substitute)"
        )
    fig.update_layout(
        title="BLAST E-value (NCBI)",
        xaxis_title="Subject accession",
        yaxis_title=y_title,
    )
    return _apply_base_layout(fig)


def blast_identity_chart(hits: List[dict]) -> go.Figure:
    """Percentual de identidade derivado das contagens NCBI (identities/align_len).

    Args:
        hits: Hits ja validados.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data.
    """
    if not hits:
        raise ValueError("Insufficient data")
    labels = _blast_hit_labels(hits)
    values = [float(hit["identity_pct"]) for hit in hits]
    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color="#34D399",
            hovertemplate="Hit %{x}<br>Identity %{y:.2f}%<extra></extra>",
        )
    )
    fig.update_layout(
        title="BLAST identity percent (NCBI identity counts)",
        xaxis_title="Subject accession",
        yaxis_title="Identity %",
        yaxis_range=[0, 100],
    )
    return _apply_base_layout(fig)


def blast_coverage_chart(hits: List[dict]) -> go.Figure:
    """Cobertura da query a partir das coordenadas NCBI do HSP.

    Args:
        hits: Hits ja validados.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data.
    """
    if not hits:
        raise ValueError("Insufficient data")
    labels = _blast_hit_labels(hits)
    values = [float(hit["query_coverage_pct"]) for hit in hits]
    fig = go.Figure(
        go.Bar(
            x=labels,
            y=values,
            marker_color="#FBBF24",
            hovertemplate="Hit %{x}<br>Query coverage %{y:.2f}%<extra></extra>",
        )
    )
    fig.update_layout(
        title="BLAST query coverage (NCBI HSP coordinates)",
        xaxis_title="Subject accession",
        yaxis_title="Query coverage %",
        yaxis_range=[0, 100],
    )
    return _apply_base_layout(fig)


def blast_alignment_overview(hits: List[dict], query_length: int) -> go.Figure:
    """Intervalos query_from-query_to reportados pelo NCBI.

    Args:
        hits: Hits ja validados.
        query_length: Comprimento da query no XML.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data ou comprimento invalido.
    """
    if not hits or query_length <= 0:
        raise ValueError("Insufficient data")
    fig = go.Figure()
    for index, hit in enumerate(hits):
        start = int(hit["query_from"])
        end = int(hit["query_to"])
        low = min(start, end)
        high = max(start, end)
        label = str(hit.get("accession") or hit.get("hit_id") or f"hit_{index + 1}")
        fig.add_trace(
            go.Scatter(
                x=[low, high],
                y=[index, index],
                mode="lines",
                line=dict(color="#6893D0", width=8),
                name=label[:40],
                hovertemplate=f"{label}<br>query {start}-{end}<extra></extra>",
                showlegend=False,
            )
        )
    fig.update_layout(
        title="BLAST alignment overview (query coordinates from NCBI)",
        xaxis_title="Query position (1-based, NCBI)",
        yaxis_title="Hit index",
        xaxis_range=[1, query_length],
    )
    return _apply_base_layout(fig)


def msa_conservation_chart(
    columns: List[int], scores: List[float], consensus: List[str]
) -> go.Figure:
    """Perfil de conservacao Shannon por coluna do MSA validado.

    Args:
        columns: Indices 0-based da janela visivel.
        scores: Scores [0, 1] ou NaN (coluna so de gaps).
        consensus: Simbolo de consenso por coluna, mesmo comprimento.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data se vazio ou comprimentos divergentes.
    """
    if not columns or not scores or len(columns) != len(scores):
        raise ValueError("Insufficient data")
    if consensus and len(consensus) != len(columns):
        raise ValueError("Insufficient data")
    xs: List[int] = []
    ys: List[float] = []
    labels: List[str] = []
    for index, (column, score) in enumerate(zip(columns, scores)):
        if isinstance(score, float) and math.isnan(score):
            continue
        xs.append(int(column))
        ys.append(float(score))
        symbol = consensus[index] if consensus else ""
        labels.append(symbol)
    if not xs:
        raise ValueError("Insufficient data")
    fig = go.Figure(
        go.Scatter(
            x=xs,
            y=ys,
            mode="lines+markers",
            line=dict(color="#6893D0", width=1.5),
            marker=dict(size=5),
            customdata=labels,
            hovertemplate=(
                "Column %{x}<br>Conservation %{y:.3f}<br>Consensus %{customdata}"
                "<extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title="MSA conservation (1 - Shannon entropy / log2(alphabet))",
        xaxis_title="Alignment column (0-based)",
        yaxis_title="Conservation [0, 1]",
        yaxis_range=[0, 1],
    )
    return _apply_base_layout(fig)


def rna_arc_diagram(
    sequence: str,
    paths: List[dict],
    selected: Optional[int] = None,
) -> go.Figure:
    """Diagrama de arcos derivado dos pares MFE. Sem pares ficticios.

    Args:
        sequence: RNA dobrada.
        paths: Saida de rna_folding.arc_paths.
        selected: Posicao 0-based a destacar, ou None.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data se a sequencia for vazia.
    """
    residues = str(sequence or "")
    if not residues:
        raise ValueError("Insufficient data")
    fig = go.Figure()
    xs = list(range(len(residues)))
    ys = [0.0] * len(residues)
    colors = [NUCLEOTIDE_COLORS.get(base, "#94A3B8") for base in residues]
    sizes = [12 if selected is not None and index == selected else 8 for index in xs]
    hover = [
        f"pos {index} (0-based)<br>{base}"
        for index, base in enumerate(residues)
    ]
    for path in paths:
        fig.add_trace(
            go.Scatter(
                x=list(path.get("x") or []),
                y=list(path.get("y") or []),
                mode="lines",
                line=dict(color="rgba(56,189,248,0.7)", width=1.5),
                hovertemplate=(
                    f"{path.get('base_i')}{path.get('i')} - "
                    f"{path.get('base_j')}{path.get('j')}"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers+text",
            marker=dict(size=sizes, color=colors),
            text=list(residues),
            textposition="bottom center",
            hovertext=hover,
            hoverinfo="text",
            showlegend=False,
        )
    )
    fig.update_layout(
        title="Predicted MFE base pairs (arc diagram, 0-based). Not a 3D structure.",
        xaxis_title="Sequence position (0-based)",
        yaxis_title="Pair span",
        yaxis=dict(visible=False),
        dragmode="zoom",
        hovermode="closest",
    )
    return _apply_base_layout(fig)


def rna_circular_pairs(
    layout: dict,
    selected: Optional[int] = None,
) -> go.Figure:
    """Diagrama circular das cordas de pareamento preditas. Nao e RNAPlot.

    Args:
        layout: Saida de rna_folding.circular_layout.
        selected: Posicao 0-based a destacar, ou None.

    Returns:
        Figura Plotly.

    Raises:
        ValueError: Insufficient data se o layout estiver vazio.
    """
    xs = list(layout.get("x") or [])
    ys = list(layout.get("y") or [])
    bases = list(layout.get("bases") or [])
    if not xs or not bases or len(xs) != len(bases):
        raise ValueError("Insufficient data")
    fig = go.Figure()
    for chord in list(layout.get("chords") or []):
        fig.add_trace(
            go.Scatter(
                x=list(chord.get("x") or []),
                y=list(chord.get("y") or []),
                mode="lines",
                line=dict(color="rgba(52,211,153,0.8)", width=1.5),
                hovertemplate=(
                    f"{chord.get('base_i')}{chord.get('i')} - "
                    f"{chord.get('base_j')}{chord.get('j')}"
                    "<extra></extra>"
                ),
                showlegend=False,
            )
        )
    colors = [NUCLEOTIDE_COLORS.get(base, "#94A3B8") for base in bases]
    sizes = [14 if selected is not None and index == selected else 9 for index in range(len(bases))]
    hover = [f"pos {index} (0-based)<br>{base}" for index, base in enumerate(bases)]
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="markers+text",
            marker=dict(size=sizes, color=colors),
            text=bases,
            textposition="top center",
            hovertext=hover,
            hoverinfo="text",
            showlegend=False,
        )
    )
    fig.update_layout(
        title="Predicted pairing topology (circular). Not RNAPlot and not 3D.",
        xaxis=dict(visible=False, scaleanchor="y"),
        yaxis=dict(visible=False),
        dragmode="zoom",
        hovermode="closest",
    )
    return _apply_base_layout(fig, cartesian=False)


def phylogenetic_tree_figure(
    layout: dict,
    *,
    selected_tree_id: str = "",
    color_groups: Optional[Dict[str, str]] = None,
    show_branch_lengths: bool = True,
    show_support: bool = True,
    rooted_label: str = "",
) -> go.Figure:
    """Desenha a arvore a partir dos edges reais do layout. Sem dendrograma.

    Args:
        layout: Saida de phylogeny.layout_phylogram (nodes, edges, scale_bar).
        selected_tree_id: Folha destacada, se houver.
        color_groups: tree_id da folha -> categoria taxonomica real.
        show_branch_lengths: Mostra escala so se scale_bar.present.
        show_support: Rotula suporte apenas quando o valor existe.
        rooted_label: Texto de rooting (unrooted / midpoint / ...).

    Returns:
        Figura Plotly. Cada linha e um edge da inferencia.

    Raises:
        ValueError: Insufficient data se faltar node/edge.

    Nota biologica:
        Layout visual nao altera a topologia. Cor taxonomica nao e distancia
        filogenetica. Escala ausente se nao houver branch lengths.
    """
    nodes = list((layout or {}).get("nodes") or [])
    edges = list((layout or {}).get("edges") or [])
    if not nodes or not edges:
        raise ValueError("Insufficient data")
    by_id = {str(item.get("id")): item for item in nodes}
    fig = go.Figure()
    for edge in edges:
        parent = by_id.get(str(edge.get("parent")))
        child = by_id.get(str(edge.get("child")))
        if parent is None or child is None:
            raise ValueError("Insufficient data")
        x0 = float(edge.get("x0"))
        y0 = float(edge.get("y0"))
        x1 = float(edge.get("x1"))
        y1 = float(edge.get("y1"))
        length = edge.get("length")
        hover = f"length {length}" if length is not None else "branch length N/A"
        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(color="rgba(148,163,184,0.85)", width=1.6),
                hovertemplate=hover + "<extra></extra>",
                showlegend=False,
            )
        )
    groups = color_groups or {}
    palette = [
        "#6893D0",
        "#34D399",
        "#FBBF24",
        "#F472B6",
        "#A78BFA",
        "#FB923C",
        "#2DD4BF",
        "#F87171",
    ]
    unique_groups: List[str] = []
    for item in nodes:
        if not item.get("is_leaf"):
            continue
        label = str(item.get("label") or "")
        group = groups.get(label, "")
        if group and group not in unique_groups:
            unique_groups.append(group)
    group_color = {
        name: palette[index % len(palette)] for index, name in enumerate(unique_groups)
    }
    leaf_x: List[float] = []
    leaf_y: List[float] = []
    leaf_text: List[str] = []
    leaf_color: List[str] = []
    leaf_size: List[int] = []
    leaf_hover: List[str] = []
    support_x: List[float] = []
    support_y: List[float] = []
    support_text: List[str] = []
    for item in nodes:
        if item.get("is_leaf"):
            label = str(item.get("label") or "")
            group = groups.get(label, "")
            color = group_color.get(group, "#F8FAFC")
            size = 14 if selected_tree_id and label == selected_tree_id else 9
            leaf_x.append(float(item["x"]))
            leaf_y.append(float(item["y"]))
            leaf_text.append(label)
            leaf_color.append(color)
            leaf_size.append(size)
            hover = label
            if group:
                hover += f"<br>{group}"
            leaf_hover.append(hover)
            continue
        support = item.get("support")
        if show_support and support is not None:
            support_x.append(float(item["x"]))
            support_y.append(float(item["y"]))
            support_text.append(f"{float(support):.0f}")
    n_leaves = len(leaf_x)
    show_leaf_text = n_leaves <= scale_profile.MAX_TREE_LEAF_LABELS
    fig.add_trace(
        go.Scatter(
            x=leaf_x,
            y=leaf_y,
            mode="markers+text" if show_leaf_text else "markers",
            marker=dict(size=leaf_size, color=leaf_color, line=dict(width=0)),
            text=leaf_text if show_leaf_text else None,
            textposition="middle right",
            hovertext=leaf_hover,
            hoverinfo="text",
            name="Leaves",
            showlegend=False,
        )
    )
    if support_x:
        fig.add_trace(
            go.Scatter(
                x=support_x,
                y=support_y,
                mode="text",
                text=support_text,
                textfont=dict(size=10, color="#FBBF24"),
                hoverinfo="skip",
                showlegend=False,
                name="Bootstrap support",
            )
        )
    scale = (layout or {}).get("scale_bar") or {}
    if show_branch_lengths and scale.get("present") and scale.get("length"):
        length = float(scale["length"])
        fig.add_trace(
            go.Scatter(
                x=[0, length],
                y=[-1.2, -1.2],
                mode="lines+text",
                line=dict(color="#F8FAFC", width=2),
                text=["", str(scale.get("label") or f"{length:g}")],
                textposition="bottom right",
                hoverinfo="skip",
                showlegend=False,
            )
        )
    if unique_groups:
        for name in unique_groups:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(size=8, color=group_color[name]),
                    name=name,
                )
            )
    title = "Phylogenetic inference"
    if rooted_label:
        title = f"Phylogenetic inference ({rooted_label})"
    if not show_leaf_text:
        title = (
            f"{title} — {n_leaves} leaves; labels in hover only "
            f"(cap {scale_profile.MAX_TREE_LEAF_LABELS} on-canvas labels)"
        )
    fig.update_xaxes(zeroline=False)
    fig = _apply_base_layout(fig)
    fig.update_layout(
        title=title,
        xaxis_title="Branch length (algorithm units)" if scale.get("present") else "",
        yaxis=dict(visible=False, zeroline=False),
        dragmode="pan",
        hovermode="closest",
        legend_title_text="Taxonomic group (visual only)" if unique_groups else "",
        height=min(900, max(360, 240 + n_leaves * 18)),
        margin=dict(l=40, r=40, t=56, b=48),
    )
    return fig


