"""Deteccao de executaveis cientificos opcionais sem varrer o disco.

Usa shutil.which, basename allowlisted, caminhos de ficheiro explicitos e
variaveis de ambiente validadas. Nao percorre C:\\ de forma recursiva.
Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import os
import platform
import shutil
from typing import Callable, Optional, Sequence

WhichFn = Callable[[str], Optional[str]]

VIENNARNA_WINDOWS_CANDIDATES: tuple[str, ...] = (
    r"C:\Program Files\ViennaRNA\RNAfold.exe",
    r"C:\Program Files (x86)\ViennaRNA\RNAfold.exe",
)
"""Caminhos de ficheiro conhecidos no Windows. Nao sao diretorios a varrer."""

TOOL_STATUS_INSTALLED: str = "installed"
TOOL_STATUS_NOT_INSTALLED: str = "not_installed"
TOOL_STATUS_INVALID: str = "invalid_executable"
TOOL_STATUS_VERSION_UNAVAILABLE: str = "version_unavailable"


def normalized_tool_basename(path: str) -> str:
    """Basename em minusculas, sem .exe/.bat/.cmd.

    Args:
        path: Caminho ou nome de ficheiro.

    Returns:
        Basename normalizado.

    Raises:
        Nenhum.
    """
    text = str(path or "").replace("\\", "/")
    base = os.path.basename(text).lower()
    for suffix in (".exe", ".bat", ".cmd"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def allowlisted_basenames(names: Sequence[str]) -> set[str]:
    """Conjunto de basenames permitidos, com e sem extensao Windows.

    Args:
        names: Nomes allowlisted (RNAfold, RNAfold.exe, ...).

    Returns:
        Conjunto em minusculas.

    Raises:
        Nenhum.
    """
    allowed: set[str] = set()
    for item in names:
        low = str(item).lower()
        allowed.add(low)
        allowed.add(normalized_tool_basename(low))
    return allowed


def sanitize_tool_path(path: str) -> str:
    """Devolve um caminho reportavel sem expandir o perfil do utilizador.

    Args:
        path: Caminho absoluto detectado.

    Returns:
        Basename, ou o caminho se estiver sob Program Files.

    Raises:
        Nenhum.
    """
    text = str(path or "").replace("\\", "/")
    if "\x00" in text or ".." in text.split("/"):
        return ""
    resolved = os.path.abspath(text)
    base = os.path.basename(text)
    if base in {"", ".", ".."}:
        return ""
    lowered = text.lower()
    program_files = os.environ.get("ProgramFiles", r"C:\Program Files").replace("\\", "/").lower()
    program_files_x86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)").replace("\\", "/").lower()
    if lowered.startswith(program_files) or lowered.startswith(program_files_x86):
        return resolved
    return base


def classify_tool_record(
    *,
    available: bool,
    path: str,
    version: str,
    invalid: bool = False,
) -> str:
    """Classifica o estado de uma ferramenta opcional.

    Args:
        available: True se o executavel passou a allowlist e existe.
        path: Caminho resolvido ou vazio.
        version: Versao detectada ou vazia.
        invalid: True se um candidato existia mas foi recusado.

    Returns:
        installed, not_installed, invalid_executable ou version_unavailable.

    Raises:
        Nenhum.
    """
    if invalid:
        return TOOL_STATUS_INVALID
    if not available or not path:
        return TOOL_STATUS_NOT_INSTALLED
    if not str(version or "").strip():
        return TOOL_STATUS_VERSION_UNAVAILABLE
    return TOOL_STATUS_INSTALLED


def resolve_allowlisted_executable(
    names: Sequence[str],
    *,
    extra_file_candidates: Sequence[str] = (),
    env_var: str = "",
    which_fn: Optional[WhichFn] = None,
) -> Optional[str]:
    """Resolve um executavel allowlisted sem varrer o disco.

    Ordem: variavel de ambiente (ficheiro explicito), shutil.which, depois
    candidatos de ficheiro conhecidos. Cada caminho e aceite so se o basename
    estiver na allowlist e for um ficheiro regular.

    Args:
        names: Basenames permitidos.
        extra_file_candidates: Caminhos de ficheiro explicitos, nao pastas.
        env_var: Nome de variavel de ambiente com um caminho de ficheiro.
        which_fn: shutil.which injetavel.

    Returns:
        Caminho absoluto ou None.

    Raises:
        Nenhum.
    """
    allowed = allowlisted_basenames(names)
    finder = which_fn or shutil.which
    candidates: list[str] = []
    if env_var:
        env_path = (os.environ.get(env_var) or "").strip().strip('"')
        if env_path:
            candidates.append(env_path)
    for name in names:
        base = os.path.basename(str(name))
        if base.lower() not in allowed:
            continue
        found = finder(base)
        if found:
            candidates.append(found)
    for extra in extra_file_candidates:
        if extra:
            candidates.append(str(extra))
    seen: set[str] = set()
    for raw in candidates:
        resolved = os.path.abspath(raw)
        if resolved in seen:
            continue
        seen.add(resolved)
        found_name = os.path.basename(resolved).lower()
        if found_name not in allowed and normalized_tool_basename(found_name) not in allowed:
            continue
        if not os.path.isfile(resolved):
            continue
        return resolved
    return None


def environment_report(*, tools: Optional[dict] = None) -> dict:
    """Diagnostico de plataforma sem dados pessoais.

    Args:
        tools: Mapa opcional nome -> dict de disponibilidade.

    Returns:
        platform, python_version, architecture e ferramentas.

    Raises:
        Nenhum.
    """
    records = {}
    for name, info in dict(tools or {}).items():
        path = str((info or {}).get("path") or "")
        records[str(name)] = {
            "available": bool((info or {}).get("available")),
            "version": str((info or {}).get("version") or ""),
            "version_status": str((info or {}).get("version_status") or ""),
            "status": classify_tool_record(
                available=bool((info or {}).get("available")),
                path=path,
                version=str((info or {}).get("version") or ""),
                invalid=bool((info or {}).get("invalid")),
            ),
            "path": sanitize_tool_path(path) if path else "",
        }
    return {
        "platform": platform.system(),
        "platform_release": platform.release(),
        "python_version": platform.python_version(),
        "architecture": platform.machine(),
        "tools": records,
    }
