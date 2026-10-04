"""Perfis de escala para analise grande. Nao altera formulas.

50_000 nucleotidos e um alvo de ANALISE, nao de render atomico completo.
Limites abaixo sao tetos tecnicos explicitos, nao regras biologicas.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Dict

TARGET_DNA_NT: int = 50_000
FULL_SCALE_TARGET_NT: int = 150_000
MAX_WINDOW_PROFILE_ROWS: int = 2_500
MAX_PLOT_POINTS: int = 2_500
MAX_SEQUENCE_VIEWER_RESIDUES: int = 12_000
MAX_TREE_LEAF_LABELS: int = 16
MAX_HELIX_NT: int = 400


def recommended_window_step(length: int, window: int, *, max_rows: int = MAX_WINDOW_PROFILE_ROWS) -> int:
    """Passo minimo para nao exceder max_rows janelas.

    Args:
        length: Comprimento da sequencia.
        window: Tamanho da janela.
        max_rows: Teto de linhas (tecnico).

    Returns:
        Step >= 1. Nao muda a definicao de GC; muda a resolucao da grelha.

    Raises:
        ValueError: Se length/window/max_rows nao forem positivos.
    """
    if length <= 0 or window <= 0 or max_rows <= 0:
        raise ValueError("length, window and max_rows must be positive.")
    if length < window:
        return 1
    n_if_step_one = length - window + 1
    if n_if_step_one <= max_rows:
        return 1
    step = (n_if_step_one + max_rows - 1) // max_rows
    return max(1, int(step))


def window_count(length: int, window: int, step: int) -> int:
    """Numero de janelas na grelha inclusiva.

    Args:
        length: Comprimento.
        window: Janela.
        step: Passo.

    Returns:
        Inteiro >= 0.

    Raises:
        ValueError: Se window ou step nao forem positivos.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive.")
    if length < window:
        return 0
    return 1 + (length - window) // step


def analysis_plan(length: int, *, molecule: str = "DNA") -> Dict[str, object]:
    """Recomendacoes de UI para uma sequencia deste comprimento.

    Args:
        length: Residuos.
        molecule: DNA, RNA ou PROTEIN.

    Returns:
        Dict com tetos e passos sugeridos. Nao executa analise.

    Raises:
        Nenhum.
    """
    n = max(0, int(length))
    mol = str(molecule or "DNA").upper()
    window = 100 if n >= 100 else max(1, n)
    if n >= TARGET_DNA_NT:
        window = min(2000, max(200, n // 50))
    elif n >= 10_000:
        window = min(1000, max(100, n // 50))
    step = recommended_window_step(n, window)
    return {
        "length": n,
        "molecule": mol,
        "target_analysis_nt": TARGET_DNA_NT,
        "suggested_window": window,
        "suggested_step": step,
        "max_window_profile_rows": MAX_WINDOW_PROFILE_ROWS,
        "max_plot_points": MAX_PLOT_POINTS,
        "max_sequence_viewer_residues": MAX_SEQUENCE_VIEWER_RESIDUES,
        "max_illustrative_helix_nt": MAX_HELIX_NT,
        "full_scale_target_nt": FULL_SCALE_TARGET_NT,
        "notes": (
            "50k nt sizes the default window grid. Linear DNA statistics run on "
            "the full input up to the technical residue cap, including the "
            "150k target. Sliding-window plots use an explicit step so the grid "
            "stays within the display cap. Illustrative 3D uses a consecutive "
            "window of at most 400 nt and does not shorten the analysis. "
            "Protein residue counts are not treated as equivalent to nucleotide counts."
        ),
    }


def enforce_window_budget(
    length: int,
    window: int,
    step: int,
    *,
    max_rows: int = MAX_WINDOW_PROFILE_ROWS,
) -> int:
    """Recusa grelhas que excedam o teto tecnico. Nao subamostra em silencio.

    Args:
        length: Comprimento da sequencia.
        window: Tamanho da janela.
        step: Passo.
        max_rows: Teto de linhas.

    Returns:
        Numero de janelas se estiver dentro do teto.

    Raises:
        ValueError: RESOURCE_LIMIT textual se a grelha for grande demais.
    """
    n = window_count(length, window, step)
    if n > max_rows:
        suggested = recommended_window_step(length, window, max_rows=max_rows)
        raise ValueError(
            f"Sliding-window grid would produce {n:,} rows "
            f"(window={window}, step={step}). Technical display/memory limit is "
            f"{max_rows:,}. Increase the step (suggested >= {suggested}). "
            "The GC/entropy formulas are unchanged. This is not a biological rule."
        )
    return n