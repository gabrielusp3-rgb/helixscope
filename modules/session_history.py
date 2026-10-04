"""Historico de navegacao e de analises da sessao atual.

Nao grava disco, nao compartilha dados entre usuarios e nao calcula biologia.
O texto cientifico continua nos modulos que produziram o resultado.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from typing import Any, Mapping, MutableMapping, Optional, Sequence

MAX_MODULE_HISTORY: int = 24
MAX_ANALYSIS_ENTRIES: int = 100
MAX_INPUT_CHARS: int = 2_000_000
MAX_PAYLOAD_CHARS: int = 60_000

ANALYSIS_HISTORY_KEY: str = "helix_analysis_history"

_OBSERVERS: tuple[dict[str, str], ...] = (
    {"module": "dna", "analysis_type": "DNA composition", "key": "dna_input", "status": "COMPUTED"},
    {"module": "rna", "analysis_type": "RNA composition", "key": "rna_input", "status": "COMPUTED"},
    {"module": "protein", "analysis_type": "Protein profile", "key": "protein_input", "status": "COMPUTED"},
    {"module": "alignment", "analysis_type": "Pairwise alignment", "key": "align_result", "status": "COMPUTED"},
    {"module": "msa", "analysis_type": "Multiple sequence alignment", "key": "msa_result", "status": "COMPUTED"},
    {"module": "phylogeny", "analysis_type": "Phylogenetic tree", "key": "phylo_result", "status": "COMPUTED"},
    {"module": "crispr", "analysis_type": "CRISPR guides", "key": "crispr_result", "status": "COMPUTED"},
)


def blank_history() -> dict[str, Any]:
    """Estado vazio do historico de analises.

    Returns:
        Dict com entries, inputs, suppressed, next_id, refused e input_chars.

    Raises:
        Nenhum.

    Nota biologica:
        Este objeto nao e um resultado cientifico. Ele apenas aponta para
        analises ja calculadas nesta sessao.
    """
    return {
        "entries": [],
        "inputs": {},
        "suppressed": [],
        "next_id": 1,
        "refused": "",
        "input_chars": 0,
    }


def visit_module(
    back: Sequence[str],
    forward: Sequence[str],
    current: str,
    target: str,
    *,
    known: set[str],
) -> tuple[list[str], list[str], str]:
    """Abre um modulo e descarta o caminho de avancar.

    Args:
        back: Modulos anteriores, do mais antigo ao mais recente.
        forward: Modulos abandonados por Voltar.
        current: Modulo ativo.
        target: Modulo escolhido na barra ou na paleta.
        known: Identificadores registrados.

    Returns:
        Novo back, forward vazio e o modulo ativo.

    Raises:
        Nenhum. Um alvo desconhecido nao muda o estado.

    Nota biologica:
        Trocar de modulo nao recalcula sequencias, arvores ou estruturas.
    """
    if target not in known or target == current:
        return list(back), list(forward), current
    previous = list(back)
    if current in known:
        previous.append(current)
    return previous[-MAX_MODULE_HISTORY:], [], target


def go_back(
    back: Sequence[str],
    forward: Sequence[str],
    current: str,
    *,
    known: set[str],
) -> tuple[list[str], list[str], str]:
    """Volta um modulo na ordem real da visita.

    Args:
        back: Pilha de modulos anteriores.
        forward: Pilha de modulos a frente.
        current: Modulo ativo.
        known: Identificadores registrados.

    Returns:
        Pilhas atualizadas e o modulo restaurado. Sem historico, nada muda.

    Raises:
        Nenhum.

    Nota biologica:
        Voltar nao apaga a analise do modulo que ficou para tras.
    """
    previous = list(back)
    ahead = list(forward)
    if not previous:
        return previous, ahead, current
    target = previous.pop()
    if current in known:
        ahead.append(current)
    if target not in known:
        return previous, ahead, current
    return previous[-MAX_MODULE_HISTORY:], ahead[-MAX_MODULE_HISTORY:], target


def go_forward(
    back: Sequence[str],
    forward: Sequence[str],
    current: str,
    *,
    known: set[str],
) -> tuple[list[str], list[str], str]:
    """Avanca para o modulo abandonado pelo ultimo Voltar.

    Args:
        back: Pilha de modulos anteriores.
        forward: Pilha de modulos a frente.
        current: Modulo ativo.
        known: Identificadores registrados.

    Returns:
        Pilhas atualizadas e o modulo restaurado. Sem destino, nada muda.

    Raises:
        Nenhum.

    Nota biologica:
        Avancar nao dispara uma analise nova.
    """
    previous = list(back)
    ahead = list(forward)
    if not ahead:
        return previous, ahead, current
    target = ahead.pop()
    if current in known:
        previous.append(current)
    if target not in known:
        return previous, ahead, current
    return previous[-MAX_MODULE_HISTORY:], ahead[-MAX_MODULE_HISTORY:], target


def _digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest()


def _fingerprint(module: str, analysis_type: str, input_hash: str) -> str:
    return f"{module}|{analysis_type}|{input_hash}"


def _summary_for(value: object) -> tuple[int, str, str, Optional[dict]]:
    if isinstance(value, str):
        cleaned = "".join(value.split())
        return len(cleaned), _digest(cleaned), f"{len(cleaned):,} residues", None
    if isinstance(value, dict):
        raw = json.dumps(value, default=str, sort_keys=True)
        length = int(value.get("n_leaves") or value.get("length") or value.get("n_sequences") or 0)
        identity = str(value.get("alignment_hash") or value.get("sequence_hash") or _digest(raw[:4000]))
        status = str(value.get("status") or "")
        label = status or "stored result"
        payload = json.loads(raw) if len(raw) <= MAX_PAYLOAD_CHARS else None
        return length, identity, label, payload
    return 0, "", "", None


def ensure_history(state: MutableMapping[str, Any]) -> dict[str, Any]:
    """Garante o dict de historico dentro do estado da sessao.

    Args:
        state: Mapa mutavel da sessao.

    Returns:
        O historico existente ou um historico vazio gravado em state.

    Raises:
        Nenhum.

    Nota biologica:
        Nao interpreta sequencias.
    """
    current = state.get(ANALYSIS_HISTORY_KEY)
    if isinstance(current, dict) and isinstance(current.get("entries"), list):
        current.setdefault("inputs", {})
        current.setdefault("suppressed", [])
        current.setdefault("next_id", 1)
        current.setdefault("refused", "")
        current.setdefault("input_chars", 0)
        return current
    fresh = blank_history()
    state[ANALYSIS_HISTORY_KEY] = fresh
    return fresh


def sync_analysis_history(state: MutableMapping[str, Any]) -> dict[str, Any]:
    """Registra analises ja presentes na sessao, sem duplicar nem recriar exclusoes.

    Args:
        state: Sessao com as chaves cientificas e o historico.

    Returns:
        O historico atualizado.

    Raises:
        Nenhum.

    Nota biologica:
        Uma entrada COMPUTED ou PREDICTED aponta para um resultado ja produzido.
        Este registro nao recalcula GC, arvore, MFE ou alinhamento.
    """
    history = ensure_history(state)
    suppressed = set(str(item) for item in history.get("suppressed") or [])
    known = {
        _fingerprint(str(item.get("module")), str(item.get("analysis_type")), str(item.get("input_hash")))
        for item in history["entries"]
    }
    for spec in _OBSERVERS:
        value = state.get(spec["key"])
        if value is None or value == "" or value == {} or value == []:
            continue
        length, identity, label, payload = _summary_for(value)
        if not identity:
            continue
        finger = _fingerprint(spec["module"], spec["analysis_type"], identity)
        if finger in suppressed or finger in known:
            continue
        if len(history["entries"]) >= MAX_ANALYSIS_ENTRIES:
            history["refused"] = (
                f"History is limited to {MAX_ANALYSIS_ENTRIES} entries. "
                "Delete or clear an entry before another analysis is recorded. "
                "Older entries were not removed."
            )
            break
        inputs = history["inputs"]
        if isinstance(value, str):
            cleaned = "".join(value.split())
            if identity not in inputs:
                if int(history["input_chars"]) + len(cleaned) > MAX_INPUT_CHARS:
                    label = f"{label}; full input not retained"
                else:
                    inputs[identity] = cleaned
                    history["input_chars"] = int(history["input_chars"]) + len(cleaned)
        entry = {
            "id": f"h{int(history['next_id']):04d}",
            "module": spec["module"],
            "analysis_type": spec["analysis_type"],
            "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "input_length": length,
            "input_hash": identity,
            "parameters": {},
            "summary": label,
            "provenance": "session result already computed",
            "status": str(value.get("status") or spec["status"]) if isinstance(value, dict) else spec["status"],
            "warnings": [],
            "restore_key": spec["key"],
            "payload": payload,
        }
        history["next_id"] = int(history["next_id"]) + 1
        history["entries"].append(entry)
        known.add(finger)
    state[ANALYSIS_HISTORY_KEY] = history
    return history


def delete_analysis(history: MutableMapping[str, Any], entry_id: str) -> dict[str, Any]:
    """Remove uma entrada e impede que o mesmo resultado volte no proximo rerun.

    Args:
        history: Historico mutavel.
        entry_id: Identificador hNNNN.

    Returns:
        O mesmo historico.

    Raises:
        Nenhum. Um id ausente nao altera as entradas.

    Nota biologica:
        Apagar o registro nao apaga a sequencia ainda aberta no modulo.
    """
    kept = []
    removed = None
    for entry in list(history.get("entries") or []):
        if str(entry.get("id")) == entry_id and removed is None:
            removed = entry
            continue
        kept.append(entry)
    history["entries"] = kept
    if removed is None:
        return dict(history)
    finger = _fingerprint(
        str(removed.get("module")),
        str(removed.get("analysis_type")),
        str(removed.get("input_hash")),
    )
    suppressed = list(history.get("suppressed") or [])
    if finger not in suppressed:
        suppressed.append(finger)
    history["suppressed"] = suppressed
    _drop_unused_input(history, str(removed.get("input_hash") or ""))
    return dict(history)


def clear_analyses(history: MutableMapping[str, Any], live_fingerprints: Sequence[str] = ()) -> dict[str, Any]:
    """Esvazia o historico e suprime o que ainda esta na sessao.

    Args:
        history: Historico mutavel.
        live_fingerprints: Impressoes dos resultados atuais, para nao voltarem no rerun.

    Returns:
        Historico vazio com a lista de supressao.

    Raises:
        Nenhum.

    Nota biologica:
        Limpar o historico nao recalcula e nao apaga o modulo aberto.
    """
    suppressed = list(history.get("suppressed") or [])
    for finger in list(live_fingerprints):
        if finger not in suppressed:
            suppressed.append(str(finger))
    for entry in list(history.get("entries") or []):
        finger = _fingerprint(
            str(entry.get("module")),
            str(entry.get("analysis_type")),
            str(entry.get("input_hash")),
        )
        if finger not in suppressed:
            suppressed.append(finger)
    history["entries"] = []
    history["inputs"] = {}
    history["input_chars"] = 0
    history["suppressed"] = suppressed
    history["refused"] = ""
    return dict(history)


def live_fingerprints(state: Mapping[str, Any]) -> list[str]:
    """Impressoes dos resultados atualmente na sessao.

    Args:
        state: Sessao.

    Returns:
        Lista de impressoes module|tipo|hash.

    Raises:
        Nenhum.

    Nota biologica:
        Serve apenas para a supressao apos limpar o historico.
    """
    found = []
    for spec in _OBSERVERS:
        value = state.get(spec["key"])
        if value is None or value == "" or value == {} or value == []:
            continue
        _length, identity, _label, _payload = _summary_for(value)
        if identity:
            found.append(_fingerprint(spec["module"], spec["analysis_type"], identity))
    return found


def restore_analysis(history: Mapping[str, Any], entry_id: str, state: MutableMapping[str, Any]) -> Optional[dict]:
    """Devolve o input ou o payload guardado para a chave original.

    Args:
        history: Historico.
        entry_id: Identificador.
        state: Sessao que recebera o valor, se ele foi retido.

    Returns:
        A entrada, ou None se o id nao existe.

    Raises:
        Nenhum.

    Nota biologica:
        Restaurar nao chama motor, API nem folding. Se o texto nao foi retido,
        a entrada continua visivel e o modulo abre sem inventar sequencia.
    """
    entry = None
    for item in history.get("entries") or []:
        if str(item.get("id")) == entry_id:
            entry = item
            break
    if entry is None:
        return None
    key = str(entry.get("restore_key") or "")
    payload = entry.get("payload")
    stored = (history.get("inputs") or {}).get(str(entry.get("input_hash") or ""))
    if key and isinstance(stored, str):
        state[key] = stored
    elif key and isinstance(payload, dict):
        state[key] = payload
    return dict(entry)


def _drop_unused_input(history: MutableMapping[str, Any], input_hash: str) -> None:
    if not input_hash:
        return
    still_used = any(str(item.get("input_hash")) == input_hash for item in history.get("entries") or [])
    inputs = history.get("inputs") or {}
    if still_used or input_hash not in inputs:
        return
    text = inputs.pop(input_hash)
    history["input_chars"] = max(0, int(history.get("input_chars") or 0) - len(text))


def history_memory_chars(history: Mapping[str, Any]) -> int:
    """Conta caracteres retidos no historico.

    Args:
        history: Historico.

    Returns:
        Soma dos inputs e do JSON das entradas.

    Raises:
        Nenhum.

    Nota biologica:
        Nao mede RAM do processo. Mede o texto cientifico duplicado no historico.
    """
    entries = json.dumps(list(history.get("entries") or []), default=str)
    return int(history.get("input_chars") or 0) + len(entries)
