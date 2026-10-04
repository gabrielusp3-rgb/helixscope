"""Estado live-validated de motores externos neste processo.

DETECTED (binario encontrado) nao e LIVE_VALIDATED. A validacao live so e
registada apos execucao real com output aceite pelo parser. Mocks com run_fn
nao devem chamar record_live_validation.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, Mapping, Optional

from . import provenance, tool_detection

STATUS_NOT_INSTALLED: str = "not_installed"
STATUS_DETECTED: str = "detected"
STATUS_LIVE_VALIDATED: str = "live_validated"
STATUS_REMOTE_VALIDATED: str = "remote_validated"
STATUS_BROKEN: str = "broken"
STATUS_UNAVAILABLE: str = "unavailable"

_LIVE: Dict[str, Dict[str, Any]] = {}
LIVE_STORE_NAME: str = "engine_live.json"


def _live_store_path() -> str:
    from . import genome_store

    return os.path.join(genome_store.jobs_root(), LIVE_STORE_NAME)


def _persist_record(record: Mapping[str, Any]) -> None:
    path = _live_store_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload: Dict[str, Any] = {}
    if os.path.isfile(path):
        try:
            with open(path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
            if isinstance(loaded, dict):
                payload = loaded
        except (OSError, json.JSONDecodeError):
            payload = {}
    payload[str(record.get("tool") or "")] = dict(record)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    os.replace(tmp, path)


def _load_store() -> Dict[str, Dict[str, Any]]:
    path = _live_store_path()
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(loaded, dict):
        return {}
    rows: Dict[str, Dict[str, Any]] = {}
    for key, value in loaded.items():
        if isinstance(value, dict) and str(value.get("tool") or key):
            rows[str(value.get("tool") or key)] = dict(value)
    return rows


def record_live_validation(
    tool: str,
    *,
    ok: bool,
    version: str = "",
    details: Optional[Dict[str, Any]] = None,
) -> dict:
    """Regista o resultado de uma execucao live neste processo.

    Args:
        tool: Nome canonico (IQ-TREE, DSSP, FastTree, BLAST+).
        ok: True se output passou validacao cientifica.
        version: Versao reportada pelo binario.
        details: Metadados (modelo, n_leaves, n_residues, ...).

    Returns:
        Registro armazenado.

    Raises:
        Nenhum.
    """
    name = str(tool or "").strip()
    record = {
        "tool": name,
        "ok": bool(ok),
        "version": str(version or ""),
        "status": STATUS_LIVE_VALIDATED if ok else STATUS_BROKEN,
        "validated_at_utc": provenance.utc_now(),
        "details": dict(details or {}),
    }
    _LIVE[name] = record
    try:
        _persist_record(record)
    except OSError:
        pass
    return dict(record)


def record_remote_validation(
    tool: str,
    *,
    ok: bool,
    version: str = "",
    details: Optional[Dict[str, Any]] = None,
) -> dict:
    """Regista uma validacao remota. Nao e LIVE_VALIDATED local.

    Args:
        tool: Nome canonico (ex. RCSB Alignment API).
        ok: True se o JSON remoto passou o parser cientifico.
        version: Versao da API ou do servico.
        details: Metadados (ticket, RMSD, pares).

    Returns:
        Registro armazenado com status remote_validated ou broken.

    Raises:
        Nenhum.

    Nota:
        Um servico HTTPS nao e um motor local. remote_validated significa que
        este processo obteve um resultado COMPLETE aceite; nao implica que o
        binario rode nesta maquina.
    """
    name = str(tool or "").strip()
    payload = dict(details or {})
    payload.setdefault("engine_location", "remote")
    record = {
        "tool": name,
        "ok": bool(ok),
        "version": str(version or ""),
        "status": STATUS_REMOTE_VALIDATED if ok else STATUS_BROKEN,
        "engine_location": "remote",
        "validated_at_utc": provenance.utc_now(),
        "details": payload,
    }
    _LIVE[name] = record
    try:
        _persist_record(record)
    except OSError:
        pass
    return dict(record)


def live_record(tool: str) -> Optional[dict]:
    """Copia o registro live, ou None.

    Args:
        tool: Nome canonico.

    Returns:
        Dict ou None.

    Raises:
        Nenhum.
    """
    name = str(tool or "").strip()
    row = _LIVE.get(name)
    if isinstance(row, dict):
        return dict(row)
    stored = _load_store().get(name)
    if isinstance(stored, dict):
        _LIVE[name] = dict(stored)
        return dict(stored)
    return None


def merge_status(*, detected: bool, version: str, tool: str, broken: bool = False) -> str:
    """Combina deteccao com validacao live.

    Args:
        detected: Executavel allowlisted encontrado.
        version: Versao lida do binario (pode ser vazia).
        tool: Nome para lookup live.
        broken: True se a deteccao ja falhou (version command, etc.).

    Returns:
        not_installed, detected, live_validated, broken ou unavailable.

    Raises:
        Nenhum.
    """
    live = live_record(tool)
    if live:
        return str(live.get("status") or STATUS_BROKEN)
    if broken:
        return STATUS_BROKEN
    if not detected:
        return STATUS_NOT_INSTALLED
    if not str(version or "").strip():
        return tool_detection.TOOL_STATUS_VERSION_UNAVAILABLE
    return STATUS_DETECTED


def reset_live_validations() -> None:
    """Limpa registros live (testes).

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    _LIVE.clear()
