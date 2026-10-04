"""Navegacao de regioes e LOD visual. Nao altera dados cientificos.

Uma sequencia de 50 kb pode ser analisada por completo. O viewer 3D ilustrativo
mostra apenas uma janela contigua limitada. Subamostrar bases nao-adjacentes
numa helice e proibido: isso inventaria geometria.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

LOD_LOW: str = "low"
LOD_MEDIUM: str = "medium"
LOD_HIGH: str = "high"
LOD_LEVELS: tuple[str, ...] = (LOD_LOW, LOD_MEDIUM, LOD_HIGH)

# Limites do fragmento ilustrativo (nucleotidos consecutivos reais).
MAX_HELIX_NT_HIGH: int = 400
MAX_HELIX_NT_MEDIUM: int = 200
MAX_HELIX_NT_LOW: int = 80

# Passo visual do backbone 3D (residuos desenhados). Coordenadas restantes
# continuam as depositadas/ilustrativas; pontos omitidos nao sao interpolados.
# Nome deliberadamente NAO e STRIDE: isto nao e o programa Frishman/Argos.
LOD_BACKBONE_STEP: Dict[str, int] = {
    LOD_LOW: 4,
    LOD_MEDIUM: 2,
    LOD_HIGH: 1,
}


class RegionError(Exception):
    """Intervalo invalido. Nunca e corrigido em silencio.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def normalize_lod(lod: object) -> str:
    """Normaliza o nivel de detalhe visual.

    Args:
        lod: low, medium, high (caixa indiferente).

    Returns:
        Um de LOD_LEVELS.

    Raises:
        RegionError: INVALID_INPUT.
    """
    key = str(lod or LOD_HIGH).strip().lower()
    if key not in LOD_LEVELS:
        raise RegionError("LOD must be low, medium or high.", "INVALID_INPUT")
    return key


def helix_limit_for_lod(lod: object) -> int:
    """Maximo de nucleotidos consecutivos no modelo ilustrativo.

    Args:
        lod: Nivel visual.

    Returns:
        Inteiro positivo.

    Raises:
        RegionError.
    """
    level = normalize_lod(lod)
    if level == LOD_LOW:
        return MAX_HELIX_NT_LOW
    if level == LOD_MEDIUM:
        return MAX_HELIX_NT_MEDIUM
    return MAX_HELIX_NT_HIGH


def clamp_region(length: int, start: int, end: Optional[int] = None) -> Tuple[int, int]:
    """Valida [start, end) 0-based sobre uma sequencia.

    Args:
        length: Comprimento da sequencia.
        start: Inicio inclusivo.
        end: Fim exclusivo; None = ate o fim.

    Returns:
        (start, end) validos.

    Raises:
        RegionError: INVALID_INPUT se vazio ou fora.
    """
    if length <= 0:
        raise RegionError("Sequence length must be positive.", "INVALID_INPUT")
    try:
        start_i = int(start)
        end_i = int(length if end is None else end)
    except (TypeError, ValueError) as exc:
        raise RegionError("Region bounds must be integers.", "INVALID_INPUT") from exc
    if start_i < 0 or end_i > length or end_i <= start_i:
        raise RegionError(
            f"Region [{start_i}, {end_i}) is outside 0..{length} or empty.",
            "INVALID_INPUT",
        )
    return start_i, end_i


def contiguous_view(
    length: int,
    *,
    start: int = 0,
    end: Optional[int] = None,
    lod: object = LOD_HIGH,
    center: Optional[int] = None,
) -> dict:
    """Janela contigua para visualizacao. Nao interpola bases.

    Args:
        length: Comprimento total da sequencia analisada.
        start: Pedido do utilizador (inclusivo).
        end: Pedido do utilizador (exclusivo).
        lod: Limite de comprimento do fragmento 3D.
        center: Se definido, centra a janela neste indice 0-based.

    Returns:
        Dict start, end, n, max_nt, lod, partial_view, reason.

    Raises:
        RegionError.

    Nota biologica:
        A analise da sequencia completa nao depende desta janela. O 3D
        ilustrativo de 50 kb nao desenha 50 kb de atomos.
    """
    level = normalize_lod(lod)
    cap = helix_limit_for_lod(level)
    start_i, end_i = clamp_region(length, start, end)
    if center is not None:
        try:
            mid = int(center)
        except (TypeError, ValueError) as exc:
            raise RegionError("Center index must be an integer.", "INVALID_INPUT") from exc
        if mid < 0 or mid >= length:
            raise RegionError("Center index is outside the sequence.", "INVALID_INPUT")
        half = cap // 2
        start_i = max(0, mid - half)
        end_i = min(length, start_i + cap)
        start_i = max(0, end_i - cap)
    span = end_i - start_i
    if span > cap:
        end_i = start_i + cap
        span = cap
    partial = start_i != 0 or end_i != length
    return {
        "start": start_i,
        "end": end_i,
        "n": span,
        "max_nt": cap,
        "lod": level,
        "partial_view": partial or span < length,
        "reason": (
            "Visual window of consecutive residues. Omitted flanks are not "
            "drawn and are not interpolated. Analysis metrics use the full "
            "sequence."
            if (partial or span < length)
            else "Visual window covers the full sequence."
        ),
        "full_length": length,
    }


def backbone_stride(lod: object) -> int:
    """Passo visual do backbone. 1 = todos os residuos com coordenadas.

    Args:
        lod: Nivel.

    Returns:
        Inteiro >= 1.

    Raises:
        RegionError.
    """
    return int(LOD_BACKBONE_STEP[normalize_lod(lod)])


def downsample_indices(count: int, stride: int) -> List[int]:
    """Indices 0-based mantendo o primeiro e o ultimo.

    Args:
        count: Numero de pontos.
        stride: Passo visual.

    Returns:
        Lista de indices. Vazia se count < 1.

    Raises:
        RegionError: INVALID_INPUT se stride < 1.
    """
    if stride < 1:
        raise RegionError("Visual stride must be >= 1.", "INVALID_INPUT")
    if count <= 0:
        return []
    if stride == 1 or count <= 2:
        return list(range(count))
    kept = list(range(0, count, stride))
    if kept[-1] != count - 1:
        kept.append(count - 1)
    return kept
