"""Rastreabilidade cientifica do HelixScope.

Define os estados de evidencia e monta registros de proveniencia para que um
resultado possa ser auditado: de onde veio, como foi calculado e com quais
parametros. Nenhuma funcao inventa valores cientificos.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import hashlib
import math
import os
from datetime import datetime, timezone
from typing import Any, Dict, Mapping, Optional

HELIXSCOPE_VERSION: str = "0.24.3-19"
"""Versao do software incluida em registros de proveniencia."""

EVIDENCE_STATUSES: tuple[str, ...] = (
    "RETRIEVED",
    "COMPUTED",
    "PREDICTED",
    "HEURISTIC",
    "EXPERIMENTAL",
    "ILLUSTRATIVE",
    "UNAVAILABLE",
    "ERROR",
    "PARTIAL",
    "UNMAPPED",
    "STALE",
    "RESOURCE_LIMIT",
)
"""Estados de origem permitidos para um resultado cientifico.

PARTIAL, UNMAPPED, STALE e RESOURCE_LIMIT sao estados de primeira classe
(nao apenas notas de mapping/cena). Um valor nesta tupla que nenhum caminho
de codigo atribua e um estado morto e deve ser reportado na auditoria.
"""

STRUCTURE_KINDS: tuple[str, ...] = (
    "experimental",
    "retrieved",
    "predicted",
    "illustrative",
    "unavailable",
)
"""Classificacao de uma estrutura molecular para o futuro renderer 3D.

O renderer nunca deve promover uma sequencia a estrutura experimental.
"""

SECRET_ENV_KEYS: frozenset[str] = frozenset(
    {
        "NCBI_API_KEY",
        "ENTREZ_API_KEY",
        "API_KEY",
        "TOKEN",
        "SECRET",
        "PASSWORD",
    }
)
"""Nomes de variaveis de ambiente que nunca entram em logs de diagnostico."""


def utc_now() -> str:
    """Devolve o instante atual em UTC no formato ISO 8601.

    Args:
        Nenhum.

    Returns:
        Timestamp com sufixo Z, por exemplo 2026-08-26T12:00:00Z.

    Raises:
        Nenhum.
    """
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_provenance(
    *,
    status: str,
    algorithm: str,
    parameters: Optional[Mapping[str, Any]] = None,
    source: str = "HelixScope local computation",
    database: str = "",
    accession: str = "",
    version: str = "",
    organism: str = "",
    input_identifier: str = "",
) -> Dict[str, Any]:
    """Monta um registro de proveniencia para um resultado cientifico.

    Args:
        status: Um de EVIDENCE_STATUSES.
        algorithm: Nome do metodo realmente usado.
        parameters: Parametros efetivos do calculo; None vira dict vazio.
        source: Origem (computacao local, NCBI Entrez, etc.).
        database: Banco externo quando houver.
        accession: Accession ou identificador recuperado.
        version: Versao do registro externo, se conhecida.
        organism: Organismo quando a fonte o informar.
        input_identifier: Identificador da entrada do usuario, se houver.

    Returns:
        Dicionario JSON-serializavel com origem, metodo, parametros, software e
        timestamp. status desconhecido e gravado como ERROR sem inventar outro.

    Raises:
        Nenhum.

    Nota biologica:
        Proveniencia nao e um score de confianca. Ela apenas documenta como o
        numero ou a anotacao foi produzido.
    """
    normalized = str(status or "ERROR").strip().upper()
    if normalized not in EVIDENCE_STATUSES:
        normalized = "ERROR"
    return {
        "status": normalized,
        "source": source,
        "database": database or "",
        "accession": accession or "",
        "version": version or "",
        "organism": organism or "",
        "retrieval_timestamp": utc_now() if normalized in {"RETRIEVED", "EXPERIMENTAL"} else "",
        "computed_timestamp": utc_now()
        if normalized
        in {
            "COMPUTED",
            "PREDICTED",
            "HEURISTIC",
            "ILLUSTRATIVE",
            "PARTIAL",
            "UNMAPPED",
            "STALE",
            "RESOURCE_LIMIT",
        }
        else "",
        "algorithm": algorithm,
        "algorithm_version": "",
        "parameters": dict(parameters or {}),
        "input_identifier": input_identifier or "",
        "input_hash": "",
        "software": "HelixScope",
        "software_version": HELIXSCOPE_VERSION,
    }


def hashes_match(stored: object, current: object) -> bool:
    """Compara hashes de sequencia sem tratar vazio como coincidencia.

    Args:
        stored: Digest gravado no resultado.
        current: Digest da entrada atual.

    Returns:
        True somente se ambos forem strings nao vazias e iguais.

    Raises:
        Nenhum.
    """
    left = str(stored or "")
    right = str(current or "")
    return bool(left) and left == right


def sequence_digest(sequence: str) -> str:
    """Hash SHA-256 da sequencia ja normalizada (sem espacos, maiusculas).

    Args:
        sequence: Sequencia que o chamador ja normalizou. Esta funcao nao
            altera alfabeto nem remove ambiguidades.

    Returns:
        Digest hexadecimal de 64 caracteres. String vazia produz o hash de "".

    Raises:
        Nenhum.

    Nota:
        O hash identifica o texto exatamente como foi analisado. Nao use esta
        funcao sobre a entrada bruta se a analise operou sobre a forma limpa.
    """
    payload = (sequence or "").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def csv_cell(value: Any) -> Any:
    """Converte valores ausentes para texto CSV inequívoco.

    Args:
        value: Celula de exportacao.

    Returns:
        O proprio valor quando e um numero finito, string ou inteiro; "N/A"
        para None e NaN; "ERROR" para infinitos. Zero permanece 0.

    Raises:
        Nenhum.
    """
    if value is None:
        return "N/A"
    if isinstance(value, float):
        if math.isnan(value):
            return "N/A"
        if math.isinf(value):
            return "ERROR"
    return value


def workspace_snapshot(
    *,
    sequence: str,
    molecule: str,
    source: str,
    identifier: str = "",
    selected_start: int | None = None,
    selected_end: int | None = None,
    selected_strand: str = "+",
    accession: str = "",
    version: str = "",
    organism: str = "",
) -> Dict[str, Any]:
    """Registra a sequencia ativa do workspace sem duplicar copias extras.

    Args:
        sequence: Sequencia ja normalizada que alimenta as analises.
        molecule: DNA, RNA ou PROTEIN.
        source: Origem (user paste, FASTA, NCBI).
        identifier: Identificador FASTA ou accession, se houver.
        selected_start: Inicio 0-based da selecao, ou None para a sequencia toda.
        selected_end: Fim exclusivo da selecao.
        selected_strand: "+" ou "-".
        accession: Accession NCBI, se a sequencia veio de Entrez.
        version: Accession versionado NCBI, se conhecido.
        organism: Organismo do registro, se conhecido.

    Returns:
        Dict com molecule, length, hash, source, identifier e selecao. A
        sequencia em si nao e copiada neste registro; o chamador guarda uma
        unica copia no session_state.

    Raises:
        Nenhum.
    """
    cleaned = sequence or ""
    start = 0 if selected_start is None else selected_start
    end = len(cleaned) if selected_end is None else selected_end
    return {
        "molecule": str(molecule or "").strip().upper(),
        "length": len(cleaned),
        "input_hash": sequence_digest(cleaned),
        "source": source,
        "identifier": identifier or "",
        "selected_start": start,
        "selected_end": end,
        "selected_strand": selected_strand,
        "timestamp": utc_now(),
        "software_version": HELIXSCOPE_VERSION,
        "accession": accession or "",
        "version": version or "",
        "organism": organism or "",
    }


def json_safe(value: Any) -> Any:
    """Converte NaN/Inf em None para JSON sem fingir que o valor e zero.

    Args:
        value: Numero, None, string ou estrutura aninhada.

    Returns:
        Copia JSON-serializavel; NaN e infinitos viram None.

    Raises:
        Nenhum.
    """
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    return value


def analysis_envelope(
    *,
    module: str,
    payload: Mapping[str, Any],
    status: str,
    algorithm: str,
    parameters: Optional[Mapping[str, Any]] = None,
    sequence: str = "",
    source: str = "HelixScope local computation",
    accession: str = "",
    organism: str = "",
    input_identifier: str = "",
) -> Dict[str, Any]:
    """Empacota um resultado analitico com proveniencia e hash da sequencia.

    Args:
        module: Nome do modulo (DNA, RNA, PROTEIN, ALIGNMENT, MOTIF).
        payload: Campos cientificos ja calculados.
        status: Um de EVIDENCE_STATUSES.
        algorithm: Metodo realmente usado.
        parameters: Parametros efetivos.
        sequence: Sequencia normalizada da qual o payload foi derivado.
        source: Origem do calculo ou recuperacao.
        accession: Accession quando houver.
        organism: Organismo quando houver.
        input_identifier: Identificador da entrada.

    Returns:
        Dict JSON-safe com proveniencia, hash, modulo e payload.

    Raises:
        Nenhum.
    """
    record = build_provenance(
        status=status,
        algorithm=algorithm,
        parameters=parameters,
        source=source,
        accession=accession,
        organism=organism,
        input_identifier=input_identifier,
    )
    if sequence:
        record["input_hash"] = sequence_digest(sequence)
        record["input_length"] = len(sequence)
    return json_safe(
        {
            **record,
            "module": module,
            **dict(payload),
        }
    )


def structure_disclaimer(kind: str) -> str:
    """Texto obrigatorio para qualquer visualizacao estrutural futura.

    Args:
        kind: Um de STRUCTURE_KINDS.

    Returns:
        Frase em ingles descrevendo o que a visualizacao e e o que ela nao e.

    Raises:
        ValueError: Se kind nao for reconhecido.

    Nota biologica:
        Uma fita de DNA nao define coordenadas atomicas. Tratar uma sequencia
        como estrutura experimental seria fabricacao.
    """
    key = str(kind or "").strip().lower()
    if key not in STRUCTURE_KINDS:
        raise ValueError(
            "Structure kind must be experimental, retrieved, predicted, "
            "illustrative or unavailable."
        )
    messages = {
        "experimental": (
            "Experimental structure. Coordinates come from a resolved experiment "
            "(for example X-ray, NMR or cryo-EM) identified by the source record."
        ),
        "retrieved": (
            "Retrieved structure. Coordinates were downloaded from a structural "
            "database; they were not calculated by HelixScope."
        ),
        "predicted": (
            "Predicted structure. Coordinates come from a computational model. "
            "They are not an experimental measurement."
        ),
        "illustrative": (
            "Illustrative representation. This image is a visual aid and must not "
            "be interpreted as atomic coordinates."
        ),
        "unavailable": (
            "Structure unavailable. HelixScope has no experimental, retrieved or "
            "predicted coordinates for this molecule."
        ),
    }
    return messages[key]


def redact_secrets(mapping: Mapping[str, Any]) -> Dict[str, Any]:
    """Copia um mapeamento omitindo chaves que parecem segredos.

    Args:
        mapping: Dict de diagnostico.

    Returns:
        Novo dict sem chaves cujo nome indica credencial, e sem valores iguais a
        NCBI_API_KEY do ambiente.

    Raises:
        Nenhum.
    """
    secret_values = {
        (os.environ.get(name) or "").strip()
        for name in SECRET_ENV_KEYS
        if (os.environ.get(name) or "").strip()
    }
    redacted: Dict[str, Any] = {}
    for key, value in mapping.items():
        key_upper = str(key).upper()
        if any(token in key_upper for token in SECRET_ENV_KEYS):
            redacted[str(key)] = "[redacted]"
            continue
        text = "" if value is None else str(value)
        if text and text in secret_values:
            redacted[str(key)] = "[redacted]"
            continue
        redacted[str(key)] = value
    return redacted


def diagnostic_event(
    *,
    module: str,
    operation: str,
    status: str,
    exception: str = "",
) -> Dict[str, Any]:
    """Monta um evento de diagnostico sem sequencias nem credenciais.

    Args:
        module: Nome do modulo (dna_analysis, ncbi_fetch, etc.).
        operation: Operacao (fetch, align, motif_search, etc.).
        status: Estado curto (ok, error, unavailable).
        exception: Classe ou mensagem curta; nao incluir e-mail nem API key.

    Returns:
        Registro interno para observabilidade.

    Raises:
        Nenhum.
    """
    return redact_secrets(
        {
            "module": module,
            "operation": operation,
            "status": status,
            "exception": exception,
            "timestamp": utc_now(),
            "software_version": HELIXSCOPE_VERSION,
        }
    )
