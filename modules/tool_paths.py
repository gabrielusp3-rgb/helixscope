"""Candidatos de executaveis no diretorio local de ferramentas.

Nao percorre C:\\. Apenas HELIXSCOPE_TOOLS_DIR (se definido) e a pasta
tools/ do repositorio, profundidade maxima 3. Basename continua allowlisted
em tool_detection.resolve_allowlisted_executable.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import os
from typing import Sequence

from . import tool_detection

TOOLS_ENV: str = "HELIXSCOPE_TOOLS_DIR"
MAX_SCAN_DEPTH: int = 3
MAX_SCAN_FILES: int = 80


def repository_tools_dir() -> str:
    """Pasta tools/ na raiz do repositorio (relativa a este modulo).

    Args:
        Nenhum.

    Returns:
        Caminho absoluto. O diretorio pode nao existir.

    Raises:
        Nenhum.
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "tools"))


def tools_roots() -> tuple[str, ...]:
    """Raizes permitidas: env + tools/ do projeto.

    Args:
        Nenhum.

    Returns:
        Tupla de caminhos absolutos, sem duplicados.

    Raises:
        Nenhum.
    """
    roots: list[str] = []
    env_path = (os.environ.get(TOOLS_ENV) or "").strip().strip('"')
    if env_path:
        roots.append(os.path.abspath(env_path))
    repo = repository_tools_dir()
    if repo not in roots:
        roots.append(repo)
    return tuple(roots)


def _scan_root(root: str, allowed: set[str]) -> list[str]:
    hits: list[str] = []
    if not os.path.isdir(root):
        return hits
    root_abs = os.path.abspath(root)
    for dirpath, dirnames, filenames in os.walk(root_abs):
        rel = os.path.relpath(dirpath, root_abs)
        depth = 0 if rel == os.curdir else rel.count(os.sep) + 1
        if depth > MAX_SCAN_DEPTH:
            dirnames[:] = []
            continue
        dirnames[:] = [name for name in dirnames if name not in {".git", "__pycache__", ".venv"}]
        for filename in filenames:
            key = filename.lower()
            stem = tool_detection.normalized_tool_basename(filename)
            if key not in allowed and stem not in allowed:
                continue
            hits.append(os.path.abspath(os.path.join(dirpath, filename)))
            if len(hits) >= MAX_SCAN_FILES:
                return hits
    return hits


def candidates_for(names: Sequence[str]) -> tuple[str, ...]:
    """Ficheiros allowlisted sob tools/ e HELIXSCOPE_TOOLS_DIR.

    Args:
        names: Basenames permitidos (iqtree3.exe, mkdssp, ...).

    Returns:
        Caminhos absolutos existentes. Lista vazia se nada estiver instalado.

    Raises:
        Nenhum.
    """
    allowed = tool_detection.allowlisted_basenames(names)
    found: list[str] = []
    seen: set[str] = set()
    for root in tools_roots():
        for path in _scan_root(root, allowed):
            if path in seen:
                continue
            seen.add(path)
            found.append(path)
    return tuple(found)
