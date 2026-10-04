"""Alinhamento multiplo real e analise comparativa derivada do MSA.

Submete conjuntos de sequencias ao Clustal Omega hospedado pela EMBL-EBI
(Job Dispatcher REST) ou a executaveis locais allowlisted (clustalo, mafft,
muscle) quando presentes. Nao calcula um MSA proprio nem rebatiza pairwise
como alinhamento multiplo.

Nenhuma funcao aqui importa Streamlit. Destinos HTTP sao fixos (sem URL do
usuario). Subprocesso usa lista de argumentos, nunca shell=True.
"""

from __future__ import annotations

import hashlib
import importlib
import io
import math
import os
import re
import shutil
import socket
import subprocess
import tempfile
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from . import dna_analysis, provenance, scientific_checks

EBI_CLUSTALO_BASE: str = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
"""Unico prefixo HTTP permitido para Clustal Omega remoto."""

EBI_TOOL: str = "HelixScope"
"""Valor de User-Agent enviado ao Job Dispatcher."""

HTTP_TIMEOUT_S: float = 45.0
"""Timeout de cada requisicao HTTP a EMBL-EBI."""

MIN_POLL_INTERVAL_S: float = 5.0
"""Intervalo minimo entre consultas de status do mesmo job EBI."""

MAX_JOB_WAIT_S: float = 300.0
"""Tempo maximo desde o submit ate TIMEOUT, sem converter em UNAVAILABLE."""

MAX_SEQUENCES: int = 30
"""Teto de sequencias no conjunto MSA (protecao de processo, nao regra biologica)."""

MAX_RESIDUES_PER_SEQUENCE: int = 8_000
"""Teto por sequencia, alinhado ao pairwise."""

MAX_RESIDUES_TOTAL: int = 80_000
"""Soma maxima de residuos nao alinhados no conjunto."""

MAX_OUTPUT_BYTES: int = 2_000_000
"""Teto do FASTA alinhado devolvido pela ferramenta."""

MAX_ALIGNMENT_COLUMNS: int = 20_000
"""Teto de colunas apos o parse."""

MAX_DISPLAY_ROWS: int = 20
"""Linhas visiveis por pagina no viewer."""

MAX_DISPLAY_COLUMNS: int = 80
"""Colunas visiveis por pagina no viewer."""

MAX_BLAST_HITS_TO_FETCH: int = 10
"""Maximo de subjects BLAST buscados no NCBI por acao."""

LOCAL_TIMEOUT_S: float = 30.0
"""Timeout do executavel local."""

JOB_ID_PATTERN: re.Pattern[str] = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{4,80}$")
"""Job ID EBI. Recusa URL, path ou comando."""

GAP_CHARS: frozenset[str] = frozenset("-.")
"""Simbolos de gap aceitos no FASTA alinhado; ponto e normalizado para hifen."""

LOCAL_TOOL_NAMES: Dict[str, Tuple[str, ...]] = {
    "clustalo": ("clustalo", "clustalo.exe"),
    "mafft": ("mafft", "mafft.bat", "mafft.exe"),
    "muscle": ("muscle", "muscle.exe", "muscle5", "muscle5.exe"),
}
"""Basenames allowlisted. Nenhum path arbitrario do usuario."""

DNA_IUPAC_FROM_SET: Dict[frozenset[str], str] = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}
"""Codigo IUPAC de DNA para um empate de bases, sem inventar outros simbolos."""


class MsaError(RuntimeError):
    """Falha classificada de MSA. Nunca deve ser reescrita como UNAVAILABLE.

    Attributes:
        category: TOOL_NOT_INSTALLED, TOOL_FAILED, TIMEOUT, RESOURCE_LIMIT,
            PARSING_ERROR, INVALID_INPUT, RATE_LIMITED, SERVICE_UNAVAILABLE
            ou JOB_FAILED.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "TOOL_FAILED")


def tool_availability() -> dict:
    """Declara backends de MSA reais: EBI Clustal Omega e binarios locais.

    Args:
        Nenhum.

    Returns:
        Dict com any_available, backends (lista) e summary_reason. Pairwise
        nao e listado como MSA.

    Raises:
        Nenhum.

    Nota biologica:
        Clustal Omega, MAFFT e MUSCLE sao algoritmos publicados. HelixScope
        nao os reproduz: executa o servico EBI ou o binario detectado.
    """
    local = detect_local_tools()
    backends = [
        {
            "id": "ebi_clustalo",
            "tool": "Clustal Omega",
            "method": "EMBL-EBI Job Dispatcher REST",
            "available": True,
            "version": "",
            "version_status": "reported_after_job_if_present",
            "reason": (
                "Remote Clustal Omega at EMBL-EBI. The algorithm runs on EBI "
                "infrastructure, not inside HelixScope."
            ),
            "endpoint": EBI_CLUSTALO_BASE,
        }
    ]
    for tool_id, info in local.items():
        backends.append(info)
    any_local = any(item["available"] for item in local.values())
    return {
        "any_available": True,
        "ebi_clustalo": True,
        "local_any": any_local,
        "backends": backends,
        "summary_reason": (
            "MSA uses EMBL-EBI Clustal Omega and, when installed, local "
            "clustalo/MAFFT/MUSCLE. Pairwise Needleman-Wunsch is not MSA."
        ),
        "pairwise_is_not_msa": True,
        "blast_is_not_msa": True,
        "phylogeny": False,
        "phylogeny_module": "modules.phylogeny",
        "structure": False,
    }


def detect_local_tools() -> Dict[str, dict]:
    """Procura apenas basenames allowlisted via shutil.which.

    Args:
        Nenhum.

    Returns:
        Mapa tool_id -> dict available, version (detectada ou vazia), path.

    Raises:
        Nenhum.
    """
    found: Dict[str, dict] = {}
    for tool_id, names in LOCAL_TOOL_NAMES.items():
        executable = _resolve_allowlisted_executable(names)
        if executable is None:
            found[tool_id] = {
                "id": f"local_{tool_id}",
                "tool": tool_id,
                "method": "local executable",
                "available": False,
                "version": "",
                "reason": f"{tool_id} executable not installed.",
                "path": "",
            }
            continue
        version = _probe_local_version(executable, tool_id)
        found[tool_id] = {
            "id": f"local_{tool_id}",
            "tool": tool_id,
            "method": "local executable",
            "available": True,
            "version": version,
            "reason": f"Found allowlisted executable {os.path.basename(executable)}.",
            "path": executable,
        }
    return found


def available_backend_ids() -> List[str]:
    """IDs de backend que a UI pode oferecer.

    Args:
        Nenhum.

    Returns:
        Lista com ebi_clustalo e locais realmente instalados.

    Raises:
        Nenhum.
    """
    ids = ["ebi_clustalo"]
    local = detect_local_tools()
    for tool_id, info in local.items():
        if info["available"]:
            ids.append(f"local_{tool_id}")
    return ids


def parse_unaligned_fasta(text: str) -> List[dict]:
    """Le FASTA multiplo para o conjunto MSA sem descartar registros.

    Args:
        text: FASTA ou sequencia crua.

    Returns:
        Lista de dicts identifier/sequence ainda nao validados como molecula.

    Raises:
        ValueError: Parse FASTA invalido ou excesso de registros.

    Nota:
        Usa dna_analysis.parse_fasta_records quando o modulo em memoria a
        expoe. Se o processo Streamlit ainda tiver um dna_analysis antigo,
        recarrega o modulo ou parseia com SeqIO aqui. Nao escolhe o primeiro
        registro em silencio.
    """
    parser = getattr(dna_analysis, "parse_fasta_records", None)
    if parser is None:
        importlib.reload(dna_analysis)
        parser = getattr(dna_analysis, "parse_fasta_records", None)
    if parser is not None:
        parsed = parser(text, max_records=MAX_SEQUENCES)
        return list(parsed.get("records") or [])
    parsed = _parse_fasta_records_local(text, max_records=MAX_SEQUENCES)
    return list(parsed.get("records") or [])


def _parse_fasta_records_local(text: str, *, max_records: int) -> dict:
    """Parser FASTA local (SeqIO) se dna_analysis estiver stale (interno)."""
    if not isinstance(text, str):
        raise TypeError("A sequencia deve ser uma string.")
    stripped = (text or "").strip()
    if not stripped:
        return {"format": "empty", "record_count": 0, "records": []}
    if not stripped.startswith(">"):
        return {
            "format": "raw",
            "record_count": 1,
            "records": [
                {"identifier": "", "sequence": stripped, "length": len(stripped)}
            ],
        }
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not records:
        raise ValueError("FASTA header found but no sequence records could be parsed.")
    if len(records) > int(max_records):
        raise ValueError(
            f"FASTA contains {len(records)} records; the limit is {max_records}. "
            "The first record was not selected automatically."
        )
    parsed = [
        {
            "identifier": str(record.id or ""),
            "sequence": str(record.seq or ""),
            "length": len(str(record.seq or "")),
        }
        for record in records
    ]
    return {"format": "fasta", "record_count": len(parsed), "records": parsed}


def build_collection_member(
    *,
    sequence: str,
    identifier: str = "",
    source: str = "user input",
    molecule: str = "",
    accession: str = "",
    version: str = "",
    organism: str = "",
    database: str = "",
    retrieved_at: str = "",
    hash_value: str = "",
    blast_rid: str = "",
    blast_program: str = "",
    blast_database: str = "",
    query_hash: str = "",
) -> dict:
    """Normaliza um membro do conjunto comparativo.

    Args:
        sequence: Sequencia nao alinhada, inalterada apos validacao.
        identifier: ID FASTA ou rotulo.
        source: user input, workspace, NCBI Entrez ou BLAST-derived selection.
        molecule: DNA, RNA ou PROTEIN; vazio detecta.
        accession, version, organism, database, retrieved_at: proveniencia NCBI.
        hash_value: SHA-256 se ja calculado; senao e derivado da sequencia.
        blast_rid, blast_program, blast_database, query_hash: proveniencia BLAST.

    Returns:
        Membro com sequence, hash, molecule e metadados. Sem recodificacao.

    Raises:
        MsaError: INVALID_INPUT se a sequencia for invalida.
    """
    info = dna_analysis.validate_sequence(sequence)
    if not info["is_valid"]:
        raise MsaError(
            info.get("rejection_reason") or "Sequence is invalid.",
            "INVALID_INPUT",
        )
    detected = str(info["type"])
    requested = str(molecule or "").strip().upper()
    if requested and requested != detected:
        raise MsaError(
            f"Declared molecule {requested} does not match detected {detected}. "
            "No silent conversion was applied.",
            "INVALID_INPUT",
        )
    residues = str(info["sequence"])
    digest = hash_value or provenance.sequence_digest(residues)
    if hash_value and hash_value != provenance.sequence_digest(residues):
        raise MsaError(
            "Provided sequence hash does not match the sequence. "
            "The sequence was not altered to force a match.",
            "INVALID_INPUT",
        )
    label = str(identifier or accession or digest[:12]).strip() or digest[:12]
    return {
        "identifier": label,
        "sequence": residues,
        "length": len(residues),
        "molecule": detected,
        "hash": digest,
        "source": source,
        "accession": accession or "",
        "version": version or "",
        "organism": organism or "",
        "database": database or "",
        "retrieved_at": retrieved_at or "",
        "blast_rid": blast_rid or "",
        "blast_program": blast_program or "",
        "blast_database": blast_database or "",
        "query_hash": query_hash or "",
        "status": "READY",
    }


def member_from_ncbi_record(record: Mapping[str, Any]) -> dict:
    """Constroi um membro a partir de um registro Entrez ja recuperado.

    Args:
        record: Payload de ncbi_fetch.fetch_by_accession.

    Returns:
        Membro com sequencia verbatim e metadados NCBI.

    Raises:
        MsaError: INVALID_INPUT se a sequencia nao estiver disponivel.
    """
    from . import ncbi_fetch

    try:
        payload = ncbi_fetch.sequence_for_workspace(dict(record))
    except ValueError as exc:
        raise MsaError(str(exc), "INVALID_INPUT") from exc
    return build_collection_member(
        sequence=str(payload["sequence"]),
        identifier=str(payload.get("version") or payload.get("accession") or ""),
        source="NCBI Entrez",
        molecule=str(payload.get("molecule") or ""),
        accession=str(payload.get("accession") or ""),
        version=str(payload.get("version") or ""),
        organism=str(payload.get("organism") or ""),
        database=str(record.get("database") or ""),
        retrieved_at=str(payload.get("timestamp") or ""),
    )


def pending_member_from_blast_hit(
    hit: Mapping[str, Any],
    blast_result: Mapping[str, Any],
) -> dict:
    """Metadados de um hit BLAST sem inventar a sequencia do subject.

    Args:
        hit: Hit validado do XML BLAST.
        blast_result: Envelope BLAST (status, rid, program, database).

    Returns:
        Membro com status SEQUENCE_UNAVAILABLE e sequence vazia.

    Raises:
        MsaError: INVALID_INPUT se o BLAST nao estiver READY ou faltar accession.
    """
    status = str(blast_result.get("status") or "")
    if status in {"TIMEOUT", "RATE_LIMITED", "SERVICE_UNAVAILABLE", "JOB_FAILED"}:
        raise MsaError(
            f"BLAST result status is {status}. Subject sequences were not invented.",
            "INVALID_INPUT",
        )
    if status not in {"READY", ""}:
        raise MsaError(
            f"BLAST hits can be selected only from a READY result, not {status}.",
            "INVALID_INPUT",
        )
    accession = str(hit.get("accession") or "").strip()
    if not accession:
        raise MsaError(
            "BLAST hit has no accession. The HSP alignment is not a full sequence.",
            "INVALID_INPUT",
        )
    return {
        "identifier": accession,
        "sequence": "",
        "length": 0,
        "molecule": "",
        "hash": "",
        "source": "BLAST-derived selection",
        "accession": accession,
        "version": "",
        "organism": str(hit.get("organism") or ""),
        "database": str(blast_result.get("database") or ""),
        "retrieved_at": "",
        "blast_rid": str(blast_result.get("rid") or ""),
        "blast_program": str(blast_result.get("program") or ""),
        "blast_database": str(blast_result.get("database") or ""),
        "query_hash": str(blast_result.get("query_hash") or ""),
        "status": "SEQUENCE_UNAVAILABLE",
        "unavailable_reason": (
            "BLAST HSP (qseq/hseq) is not the full subject sequence. "
            "Retrieve the NCBI record by accession before adding to MSA."
        ),
    }


def complete_blast_member_with_ncbi(
    pending: Mapping[str, Any], record: Mapping[str, Any]
) -> dict:
    """Anexa a sequencia NCBI a um hit BLAST selecionado, sem usar o HSP.

    Args:
        pending: Saida de pending_member_from_blast_hit.
        record: Registro Entrez recuperado.

    Returns:
        Membro READY com sequencia NCBI e proveniencia BLAST+NCBI.

    Raises:
        MsaError: INVALID_INPUT se accession divergir ou a sequencia faltar.
    """
    member = member_from_ncbi_record(record)
    expected = str(pending.get("accession") or "").strip()
    got = str(member.get("accession") or "")
    version = str(member.get("version") or "")
    if expected and expected not in {got, version, member.get("identifier")}:
        if not (
            got.startswith(expected.rstrip("."))
            or version.startswith(expected)
            or expected.startswith(got)
        ):
            raise MsaError(
                f"NCBI record {version or got} does not match BLAST accession "
                f"{expected}. The sequence was not added.",
                "INVALID_INPUT",
            )
    member["source"] = "BLAST-derived selection"
    member["blast_rid"] = str(pending.get("blast_rid") or "")
    member["blast_program"] = str(pending.get("blast_program") or "")
    member["blast_database"] = str(pending.get("blast_database") or "")
    member["query_hash"] = str(pending.get("query_hash") or "")
    if pending.get("organism") and not member.get("organism"):
        member["organism"] = str(pending.get("organism") or "")
    return member


def validate_collection(members: Sequence[Mapping[str, Any]]) -> dict:
    """Valida o conjunto antes de executar o MSA.

    Args:
        members: Lista de membros READY.

    Returns:
        Dict com molecule, n_sequences, total_residues, identical_groups,
        identifiers e hashes na ordem de entrada.

    Raises:
        MsaError: INVALID_INPUT ou RESOURCE_LIMIT.
    """
    if not isinstance(members, (list, tuple)):
        raise MsaError("MSA collection must be a list.", "INVALID_INPUT")
    ready: List[dict] = []
    for item in members:
        if not isinstance(item, Mapping):
            raise MsaError("Each MSA member must be a mapping.", "INVALID_INPUT")
        if str(item.get("status") or "READY") != "READY":
            raise MsaError(
                str(item.get("unavailable_reason") or "A selected sequence is unavailable."),
                "INVALID_INPUT",
            )
        if not str(item.get("sequence") or ""):
            raise MsaError(
                "A selected member has no sequence. HSP fragments were not used.",
                "INVALID_INPUT",
            )
        ready.append(dict(item))
    if len(ready) < 2:
        raise MsaError(
            "MSA requires at least two sequences. Pairwise alignment is a different analysis.",
            "INVALID_INPUT",
        )
    if len(ready) > MAX_SEQUENCES:
        raise MsaError(
            f"MSA is limited to {MAX_SEQUENCES} sequences.",
            "RESOURCE_LIMIT",
        )
    molecules = {str(item.get("molecule") or "") for item in ready}
    if len(molecules) != 1 or "" in molecules:
        raise MsaError(
            "All sequences must share one molecule type (DNA, RNA or PROTEIN). "
            "No silent recoding between alphabets.",
            "INVALID_INPUT",
        )
    molecule = molecules.pop()
    total = 0
    for item in ready:
        length = len(str(item.get("sequence") or ""))
        if length > MAX_RESIDUES_PER_SEQUENCE:
            raise MsaError(
                f"Each MSA sequence is limited to {MAX_RESIDUES_PER_SEQUENCE:,} residues.",
                "RESOURCE_LIMIT",
            )
        total += length
    if total > MAX_RESIDUES_TOTAL:
        raise MsaError(
            f"MSA input is limited to {MAX_RESIDUES_TOTAL:,} total residues.",
            "RESOURCE_LIMIT",
        )
    return {
        "molecule": molecule,
        "n_sequences": len(ready),
        "total_residues": total,
        "members": ready,
        "identifiers": [str(item.get("identifier") or "") for item in ready],
        "hashes": [str(item.get("hash") or "") for item in ready],
        "identical_groups": identical_sequence_groups(ready),
        "order": list(range(len(ready))),
    }


def identical_sequence_groups(members: Sequence[Mapping[str, Any]]) -> List[dict]:
    """Agrupa identificadores que compartilham o mesmo hash. Nao deduplica.

    Args:
        members: Membros da colecao.

    Returns:
        Lista de grupos com hash e identifiers quando n > 1.

    Raises:
        Nenhum.
    """
    buckets: Dict[str, List[str]] = {}
    for item in members:
        digest = str(item.get("hash") or "")
        if not digest:
            continue
        buckets.setdefault(digest, []).append(str(item.get("identifier") or ""))
    return [
        {"hash": digest, "identifiers": labels, "n": len(labels)}
        for digest, labels in buckets.items()
        if len(labels) > 1
    ]


def cache_key(
    *,
    hashes: Sequence[str],
    backend: str,
    tool_version: str,
    parameters: Mapping[str, Any],
) -> str:
    """Chave de cache: ordem + hashes + ferramenta + versao + parametros.

    Args:
        hashes: Hashes na ordem de entrada.
        backend: ebi_clustalo ou local_*.
        tool_version: Versao detectada; vazia entra na chave como unknown.
        parameters: Parametros efetivos.

    Returns:
        Digest SHA-256.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [
            "v1",
            backend,
            tool_version or "unknown",
            ",".join(hashes),
            repr(sorted((str(k), str(v)) for k, v in dict(parameters).items())),
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def seconds_until_next_poll(job: Mapping[str, Any], now: Optional[float] = None) -> float:
    """Segundos restantes antes do proximo Get de status EBI.

    Args:
        job: Job com last_poll_monotonic e submitted_monotonic.
        now: Relogio monotonic opcional.

    Returns:
        0.0 se o poll for permitido.

    Raises:
        Nenhum.
    """
    import time as time_module

    current = time_module.monotonic() if now is None else float(now)
    last_poll = job.get("last_poll_monotonic")
    submitted_raw = job.get("submitted_monotonic")
    submitted = current if submitted_raw is None else float(submitted_raw)
    if last_poll is None:
        remaining = (submitted + MIN_POLL_INTERVAL_S) - current
        return max(0.0, remaining)
    remaining = (float(last_poll) + MIN_POLL_INTERVAL_S) - current
    return max(0.0, remaining)


def submit_msa(
    members: Sequence[Mapping[str, Any]],
    email: str,
    *,
    backend: str = "ebi_clustalo",
    urlopen_fn=None,
) -> Dict[str, Any]:
    """Submete o MSA ao backend escolhido.

    Args:
        members: Colecao validada.
        email: E-mail exigido pela EMBL-EBI (e ignorado no local, mas validado
            se o backend for EBI).
        backend: ebi_clustalo ou local_clustalo/local_mafft/local_muscle.
        urlopen_fn: urlopen injetavel para testes.

    Returns:
        Job WAITING/RUNNING/COMPLETED conforme o backend.

    Raises:
        MsaError, ValueError.
    """
    import time as time_module

    validated = validate_collection(members)
    backend_id = str(backend or "ebi_clustalo").strip()
    allowed = available_backend_ids()
    if backend_id not in allowed:
        if backend_id.startswith("local_"):
            raise MsaError(
                f"{backend_id} is not installed. This is TOOL_NOT_INSTALLED, not a computed MSA.",
                "TOOL_NOT_INSTALLED",
            )
        raise MsaError(f"MSA backend '{backend_id}' is not offered.", "INVALID_INPUT")
    parameters = {
        "backend": backend_id,
        "outfmt": "fa",
        "order": "input",
        "stype": _stype(validated["molecule"]),
    }
    if backend_id == "ebi_clustalo":
        _require_email(email)
        fasta = _collection_fasta(validated["members"])
        job_id = _ebi_submit(fasta, email, parameters, urlopen_fn=urlopen_fn)
        return {
            "status": "QUEUED",
            "backend": backend_id,
            "tool": "Clustal Omega",
            "tool_version": "",
            "job_id": job_id,
            "parameters": parameters,
            "molecule": validated["molecule"],
            "n_sequences": validated["n_sequences"],
            "input_hashes": validated["hashes"],
            "identifiers": validated["identifiers"],
            "identical_groups": validated["identical_groups"],
            "members": validated["members"],
            "submitted_at_utc": provenance.utc_now(),
            "submitted_monotonic": time_module.monotonic(),
            "last_poll_monotonic": None,
            "source": "EMBL-EBI Job Dispatcher",
            "endpoint": EBI_CLUSTALO_BASE,
            "cache_status": "live",
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
    return _run_local_backend(validated, backend_id, parameters)


def poll_msa(
    job: Mapping[str, Any],
    email: str,
    *,
    urlopen_fn=None,
    now: Optional[float] = None,
) -> Dict[str, Any]:
    """Consulta o status de um job EBI, respeitando o intervalo minimo.

    Args:
        job: Job de submit_msa.
        email: E-mail de contato.
        urlopen_fn: urlopen injetavel.
        now: Relogio monotonic opcional.

    Returns:
        Copia do job com status QUEUED, RUNNING, COMPLETED ou FAILED.

    Raises:
        MsaError.
    """
    import time as time_module

    updated = dict(job)
    if str(updated.get("backend") or "") != "ebi_clustalo":
        return updated
    _require_email(email)
    current = time_module.monotonic() if now is None else float(now)
    job_id = str(updated.get("job_id") or "")
    if not JOB_ID_PATTERN.fullmatch(job_id):
        raise MsaError("MSA job id failed validation.", "INVALID_INPUT")
    submitted_raw = updated.get("submitted_monotonic")
    submitted_mono = current if submitted_raw is None else float(submitted_raw)
    if current - submitted_mono > MAX_JOB_WAIT_S:
        raise MsaError(
            f"MSA job {job_id} exceeded {int(MAX_JOB_WAIT_S)}s. This is TIMEOUT, not UNAVAILABLE.",
            "TIMEOUT",
        )
    wait = seconds_until_next_poll(updated, now=current)
    if wait > 0:
        raise MsaError(
            f"Wait {wait:.0f}s before polling MSA job {job_id}.",
            "RATE_LIMITED",
        )
    status = _ebi_status(job_id, urlopen_fn=urlopen_fn)
    updated["last_poll_monotonic"] = current
    updated["status_checked_at_utc"] = provenance.utc_now()
    if status in {"QUEUED", "WAITING", "PENDING"}:
        updated["status"] = "QUEUED"
    elif status in {"RUNNING"}:
        updated["status"] = "RUNNING"
    elif status in {"FINISHED", "COMPLETED"}:
        updated["status"] = "COMPLETED"
    elif status in {"ERROR", "FAILURE", "FAILED", "NOT_FOUND"}:
        raise MsaError(
            f"EMBL-EBI reported MSA job {job_id} as {status}. This is TOOL_FAILED, not UNAVAILABLE.",
            "TOOL_FAILED",
        )
    else:
        raise MsaError(
            f"Unrecognized EMBL-EBI job status '{status}'.",
            "JOB_FAILED",
        )
    return updated


def retrieve_msa(
    job: Mapping[str, Any],
    email: str,
    *,
    urlopen_fn=None,
) -> Dict[str, Any]:
    """Baixa o FASTA alinhado, valida e deriva consenso/conservacao/mapping.

    Args:
        job: Job COMPLETED (EBI) ou ja local.
        email: E-mail EBI.
        urlopen_fn: urlopen injetavel.

    Returns:
        Envelope MSA tipado.

    Raises:
        MsaError.
    """
    if str(job.get("backend") or "") == "ebi_clustalo":
        _require_email(email)
        if str(job.get("status") or "") != "COMPLETED":
            raise MsaError(
                "MSA results can be retrieved only when status is COMPLETED.",
                "JOB_FAILED",
            )
        job_id = str(job.get("job_id") or "")
        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise MsaError("MSA job id failed validation.", "INVALID_INPUT")
        raw = _ebi_result_aligned(job_id, urlopen_fn=urlopen_fn)
        tool_version = str(job.get("tool_version") or "")
        return build_msa_result(
            raw_alignment=raw,
            members=list(job.get("members") or []),
            tool="Clustal Omega",
            tool_version=tool_version,
            method="EMBL-EBI Job Dispatcher REST",
            parameters=dict(job.get("parameters") or {}),
            job_id=job_id,
            source="EMBL-EBI Job Dispatcher",
            identical_groups=list(job.get("identical_groups") or []),
        )
    if str(job.get("status") or "") == "COMPLETED" and job.get("result"):
        return dict(job["result"])
    raise MsaError("No MSA result is available to retrieve.", "JOB_FAILED")


def parse_aligned_fasta(body: str) -> List[dict]:
    """Interpreta FASTA alinhado. Nao corrige linhas de comprimentos diferentes.

    Args:
        body: Texto FASTA com gaps.

    Returns:
        Lista de identifier, aligned, ungapped.

    Raises:
        MsaError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    if body is None:
        raise MsaError("Aligned FASTA body is empty.", "PARSING_ERROR")
    raw = body if isinstance(body, str) else body.decode("utf-8", errors="replace")
    encoded = raw.encode("utf-8")
    if len(encoded) > MAX_OUTPUT_BYTES:
        raise MsaError(
            f"MSA output exceeds {MAX_OUTPUT_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    stripped = raw.strip()
    if not stripped:
        raise MsaError("Aligned FASTA body is empty.", "PARSING_ERROR")
    try:
        from Bio import SeqIO
        import io as io_module
    except ImportError as exc:
        raise MsaError("Biopython is required to parse aligned FASTA.", "PARSING_ERROR") from exc
    records = list(SeqIO.parse(io_module.StringIO(stripped), "fasta"))
    if len(records) < 2:
        raise MsaError(
            "Aligned FASTA does not contain at least two sequences.",
            "PARSING_ERROR",
        )
    if len(records) > MAX_SEQUENCES:
        raise MsaError(
            f"Aligned FASTA has {len(records)} rows; limit is {MAX_SEQUENCES}.",
            "RESOURCE_LIMIT",
        )
    parsed: List[dict] = []
    for record in records:
        aligned = _normalize_gaps(str(record.seq or ""))
        parsed.append(
            {
                "identifier": str(record.id or ""),
                "aligned": aligned,
                "ungapped": ungap(aligned),
            }
        )
    return parsed


def ungap(aligned: str) -> str:
    """Remove gaps de uma linha alinhada.

    Args:
        aligned: Sequencia com '-' ou '.'.

    Returns:
        Residuos sem gaps.

    Raises:
        Nenhum.
    """
    return "".join(char for char in str(aligned or "") if char not in GAP_CHARS)


def import_prealigned_fasta(text: str) -> dict:
    """Aceita FASTA ja alinhado como MSA COMPLETED. Nao corre Clustal/MAFFT.

    Args:
        text: FASTA com colunas de gap ja presentes.

    Returns:
        Envelope de MSA (status COMPLETED) com method declarando prealigned.

    Raises:
        MsaError: PARSING_ERROR ou INVALID_INPUT.

    Nota biologica:
        O utilizador declara que as sequencias ja estao alinhadas. HelixScope
        nao infere colunas homólogas. Isto nao e um MSA calculado e nao e uma
        arvore. Sequencias identicas com IDs distintos permanecem folhas distintas.
    """
    rows = parse_aligned_fasta(text)
    members = []
    for row in rows:
        members.append(
            build_collection_member(
                sequence=str(row.get("ungapped") or ""),
                identifier=str(row.get("identifier") or ""),
                source="user-supplied prealigned FASTA",
            )
        )
    molecules = {str(item.get("molecule") or "") for item in members}
    if len(molecules) != 1:
        raise MsaError(
            "Prealigned FASTA mixes molecule types. No silent conversion was applied.",
            "INVALID_INPUT",
        )
    raw_lines = []
    for row in rows:
        raw_lines.append(f">{row.get('identifier')}")
        raw_lines.append(str(row.get("aligned") or ""))
    return build_msa_result(
        raw_alignment="\n".join(raw_lines) + "\n",
        members=members,
        tool="user-supplied alignment",
        tool_version="",
        method="prealigned FASTA (user-declared aligned columns; not Clustal/MAFFT/MUSCLE)",
        parameters={"engine": "none", "user_declared_aligned": True},
        source="user input",
    )


def validate_aligned_rows(
    rows: Sequence[Mapping[str, Any]],
    members: Sequence[Mapping[str, Any]],
) -> dict:
    """Confere comprimento igual, alfabeto e correspondencia com as entradas.

    Args:
        rows: Saida de parse_aligned_fasta.
        members: Membros originais na ordem submetida.

    Returns:
        Dict com alignment_length e rows ordenadas como members.

    Raises:
        MsaError: PARSING_ERROR se houver inconsistencia. Nao corrige.
    """
    if len(rows) != len(members):
        raise MsaError(
            f"MSA returned {len(rows)} rows for {len(members)} input sequences.",
            "PARSING_ERROR",
        )
    lengths = {len(str(row.get("aligned") or "")) for row in rows}
    if len(lengths) != 1:
        raise MsaError(
            "Aligned rows do not share the same alignment length.",
            "PARSING_ERROR",
        )
    alignment_length = lengths.pop()
    if alignment_length <= 0:
        raise MsaError("Alignment length is not positive.", "PARSING_ERROR")
    if alignment_length > MAX_ALIGNMENT_COLUMNS:
        raise MsaError(
            f"Alignment has {alignment_length} columns; limit is {MAX_ALIGNMENT_COLUMNS:,}.",
            "RESOURCE_LIMIT",
        )
    molecule = str(members[0].get("molecule") or "")
    alphabet = _residue_alphabet(molecule)
    ordered = _match_rows_to_members(rows, members)
    for row, member in zip(ordered, members):
        aligned = str(row["aligned"])
        for char in aligned:
            if char in GAP_CHARS:
                continue
            if char not in alphabet:
                raise MsaError(
                    f"Aligned residue '{char}' is not valid for {molecule}.",
                    "PARSING_ERROR",
                )
        ungapped = ungap(aligned)
        original = str(member.get("sequence") or "")
        if ungapped.upper() != original.upper():
            raise MsaError(
                f"Ungapped MSA row for {member.get('identifier')} does not match "
                "the input sequence. The alignment was not silently repaired.",
                "PARSING_ERROR",
            )
        if not scientific_checks.msa_row_length_matches(aligned, alignment_length):
            raise MsaError("An aligned row length disagrees with alignment_length.", "PARSING_ERROR")
    return {"alignment_length": alignment_length, "rows": ordered}


def column_coordinate_maps(
    rows: Sequence[Mapping[str, Any]],
) -> List[List[Optional[int]]]:
    """Mapeia cada coluna alinhada para a posicao 0-based original, ou None se gap.

    Args:
        rows: Linhas alinhadas validadas.

    Returns:
        maps[seq_index][column] = posicao original ou None.

    Raises:
        Nenhum.
    """
    maps: List[List[Optional[int]]] = []
    for row in rows:
        aligned = str(row.get("aligned") or "")
        cursor = 0
        mapping: List[Optional[int]] = []
        for char in aligned:
            if char in GAP_CHARS:
                mapping.append(None)
            else:
                mapping.append(cursor)
                cursor += 1
        maps.append(mapping)
    return maps


def consensus_from_columns(
    aligned_rows: Sequence[str],
    molecule: str,
    *,
    majority_threshold: float = 0.5,
) -> dict:
    """Consenso por maioria nas colunas do MSA. Nao e sequencia biologica.

    Args:
        aligned_rows: Linhas alinhadas do mesmo comprimento.
        molecule: DNA, RNA ou PROTEIN.
        majority_threshold: Fracao minima entre residuos nao-gap para emitir
            o residuo majoritario; abaixo disso, N/X. Empates usam IUPAC (DNA)
            ou X (proteina). Coluna so de gaps vira '-'.

    Returns:
        Dict com sequence, method, threshold, status COMPUTED.

    Raises:
        MsaError: INVALID_INPUT se nao houver linhas.
    """
    if not aligned_rows:
        raise MsaError("Cannot compute consensus from an empty alignment.", "INVALID_INPUT")
    length = len(aligned_rows[0])
    if any(len(row) != length for row in aligned_rows):
        raise MsaError("Consensus requires equal-length aligned rows.", "PARSING_ERROR")
    try:
        threshold = float(majority_threshold)
    except (TypeError, ValueError) as exc:
        raise MsaError("majority_threshold must be a number.", "INVALID_INPUT") from exc
    if not (0.0 < threshold <= 1.0):
        raise MsaError("majority_threshold must be in (0, 1].", "INVALID_INPUT")
    symbols: List[str] = []
    for col in range(length):
        residues = [
            row[col] for row in aligned_rows if row[col] not in GAP_CHARS
        ]
        if not residues:
            symbols.append("-")
            continue
        counts: Dict[str, int] = {}
        for residue in residues:
            counts[residue] = counts.get(residue, 0) + 1
        top = max(counts.values())
        winners = sorted(base for base, count in counts.items() if count == top)
        if (top / len(residues)) < threshold:
            symbols.append("N" if molecule in {"DNA", "RNA"} else "X")
            continue
        if len(winners) > 1:
            symbols.append(_tie_symbol(winners, molecule))
            continue
        symbols.append(winners[0])
    return {
        "sequence": "".join(symbols),
        "method": (
            "Per-column majority of non-gap residues. If the top residue "
            "fraction is below majority_threshold: N (nucleic acid) or X "
            "(protein). Ties that still meet the threshold: DNA IUPAC or N; "
            "RNA N; protein X. All-gap: '-'."
        ),
        "majority_threshold": threshold,
        "gap_handling": "Gaps are excluded from the majority denominator.",
        "status": "COMPUTED",
        "disclaimer": (
            "Consensus is a computational summary of this alignment, not an "
            "experimental or biological sequence."
        ),
    }


def conservation_shannon(
    aligned_rows: Sequence[str],
    molecule: str,
) -> dict:
    """Conservacao por coluna: 1 - H / log2(A), H = entropia de Shannon.

    Args:
        aligned_rows: Linhas alinhadas.
        molecule: Define A=4 (DNA/RNA) ou A=20 (proteina).

    Returns:
        Dict com scores (0-1 ou NaN se coluna so de gaps), method, alphabet_size.

    Raises:
        MsaError: INVALID_INPUT se vazio.

    Nota biologica:
        Score 1.0 significa uma unica letra nao-gap. Nao e evidencia de sitio
        ativo nem de residuo evolutivamente essencial.
    """
    if not aligned_rows:
        raise MsaError("Cannot compute conservation from an empty alignment.", "INVALID_INPUT")
    length = len(aligned_rows[0])
    alphabet_size = 4 if molecule in {"DNA", "RNA"} else 20
    hmax = math.log2(alphabet_size)
    scores: List[float] = []
    gap_fractions: List[float] = []
    n_seq = len(aligned_rows)
    for col in range(length):
        chars = [row[col] for row in aligned_rows]
        n_gap = sum(1 for char in chars if char in GAP_CHARS)
        residues = [char for char in chars if char not in GAP_CHARS]
        gap_fractions.append(round(n_gap / n_seq, 6) if n_seq else float("nan"))
        if not residues:
            scores.append(float("nan"))
            continue
        counts: Dict[str, int] = {}
        for residue in residues:
            counts[residue] = counts.get(residue, 0) + 1
        entropy = 0.0
        total = float(len(residues))
        for count in counts.values():
            prob = count / total
            entropy -= prob * math.log2(prob)
        scores.append(max(0.0, min(1.0, 1.0 - (entropy / hmax))))
    return {
        "scores": scores,
        "gap_fractions": gap_fractions,
        "method": (
            "conservation = 1 - H / log2(A); H is Shannon entropy of non-gap "
            "residues in the column; A=4 for DNA/RNA and A=20 for protein. "
            "All-gap columns are NaN, not zero. Range is [0, 1], not 0-100."
        ),
        "alphabet_size": alphabet_size,
        "gap_treatment": "Gaps excluded from entropy; gap_fraction stored separately.",
        "ambiguity_treatment": "Ambiguous symbols count as their own letters.",
        "status": "COMPUTED",
        "disclaimer": (
            "Conservation is among the analyzed sequences only. It is not "
            "functional annotation and not a phylogenetic claim."
        ),
    }


def classify_column_variation(
    aligned_rows: Sequence[str],
) -> List[str]:
    """Classifica cada coluna: conserved, variable, gap, mixed.

    Args:
        aligned_rows: Linhas alinhadas.

    Returns:
        Lista de rotulos por coluna.

    Raises:
        MsaError: INVALID_INPUT se vazio.
    """
    if not aligned_rows:
        raise MsaError("Cannot classify columns of an empty alignment.", "INVALID_INPUT")
    length = len(aligned_rows[0])
    labels: List[str] = []
    n_seq = len(aligned_rows)
    for col in range(length):
        chars = [row[col] for row in aligned_rows]
        n_gap = sum(1 for char in chars if char in GAP_CHARS)
        residues = [char for char in chars if char not in GAP_CHARS]
        if n_gap == n_seq:
            labels.append("gap")
        elif not residues:
            labels.append("gap")
        elif n_gap:
            labels.append("mixed")
        elif len(set(residues)) == 1:
            labels.append("conserved")
        else:
            labels.append("variable")
    return labels


def variation_vs_reference(
    aligned_rows: Sequence[str],
    reference_index: int = 0,
) -> List[List[str]]:
    """Compara cada sequencia a uma referencia nas colunas do MSA.

    Args:
        aligned_rows: Linhas alinhadas.
        reference_index: Indice da referencia na ordem do MSA.

    Returns:
        classes[seq][col] em {conserved, substitution, insertion, deletion, gap}.
        insertion: referencia gap e query residuo; deletion: o inverso.
        gap: ambos gap.

    Raises:
        MsaError: INVALID_INPUT se o indice for invalido.

    Nota biologica:
        Insercao/delecao sao relativas a sequencia de referencia escolhida,
        nao a um genoma de referencia externo.
    """
    if not aligned_rows:
        raise MsaError("Cannot compare sequences of an empty alignment.", "INVALID_INPUT")
    if reference_index < 0 or reference_index >= len(aligned_rows):
        raise MsaError("Reference index is outside the MSA.", "INVALID_INPUT")
    ref = aligned_rows[reference_index]
    length = len(ref)
    table: List[List[str]] = []
    for row in aligned_rows:
        row_labels: List[str] = []
        for col in range(length):
            rchar = ref[col]
            qchar = row[col]
            r_gap = rchar in GAP_CHARS
            q_gap = qchar in GAP_CHARS
            if r_gap and q_gap:
                row_labels.append("gap")
            elif r_gap and not q_gap:
                row_labels.append("insertion")
            elif (not r_gap) and q_gap:
                row_labels.append("deletion")
            elif rchar == qchar:
                row_labels.append("conserved")
            else:
                row_labels.append("substitution")
        table.append(row_labels)
    return table


def identity_matrix_from_msa(aligned_rows: Sequence[str]) -> dict:
    """Identidade pairwise a partir das colunas do MSA, nao Needleman-Wunsch.

    Args:
        aligned_rows: Linhas alinhadas.

    Returns:
        Dict com matrix (fracao 0-1 ou NaN), method e denominator.

    Raises:
        MsaError: INVALID_INPUT se vazio.

    Nota biologica:
        Denominador = colunas em que ambas as sequencias tem residuo (nao gap).
        Gaps nao contam como mismatch. Ambiguos so coincidem se identicos.
    """
    n = len(aligned_rows)
    if n < 1:
        raise MsaError("Cannot compute an identity matrix from an empty alignment.", "INVALID_INPUT")
    matrix: List[List[float]] = []
    for i in range(n):
        row: List[float] = []
        for j in range(n):
            matches = 0
            denom = 0
            for a, b in zip(aligned_rows[i], aligned_rows[j]):
                if a in GAP_CHARS or b in GAP_CHARS:
                    continue
                denom += 1
                if a == b:
                    matches += 1
            if denom == 0:
                row.append(float("nan"))
            else:
                row.append(matches / denom)
        matrix.append(row)
    return {
        "matrix": matrix,
        "method": "Computed from MSA columns",
        "denominator": "Columns where neither sequence has a gap",
        "gap_treatment": "Gap columns ignored, not scored as mismatches",
        "ambiguous_residues": "Count as identity only when the symbols are identical",
        "not_needleman_wunsch": True,
        "identity_scale": "fraction_0_1",
        "identity_definition": (
            "Ungapped pairwise fraction from MSA columns where neither sequence has a gap"
        ),
        "status": "COMPUTED",
    }


def column_detail(
    result: Mapping[str, Any],
    column: int,
) -> dict:
    """Detalhe de uma coluna 0-based do MSA tipado.

    Args:
        result: Envelope de build_msa_result.
        column: Indice 0-based.

    Returns:
        Dict com residuos, consenso, conservacao, gaps, coordenadas originais.

    Raises:
        MsaError: INVALID_INPUT se a coluna for invalida.
    """
    aligned = [str(row.get("aligned") or "") for row in list(result.get("rows") or [])]
    if not aligned:
        raise MsaError("MSA result has no aligned rows.", "INVALID_INPUT")
    length = int(result.get("alignment_length") or 0)
    if column < 0 or column >= length:
        raise MsaError("Column index is outside the alignment.", "INVALID_INPUT")
    residues = [row[column] for row in aligned]
    n_gap = sum(1 for char in residues if char in GAP_CHARS)
    consensus = str((result.get("consensus") or {}).get("sequence") or "")
    scores = list((result.get("conservation") or {}).get("scores") or [])
    labels = list(result.get("column_classes") or [])
    maps = list(result.get("coordinate_maps") or [])
    original_positions = []
    for mapping in maps:
        if column < len(mapping):
            original_positions.append(mapping[column])
        else:
            original_positions.append(None)
    return {
        "column": column,
        "residues": residues,
        "identifiers": [str(row.get("identifier") or "") for row in result.get("rows") or []],
        "consensus": consensus[column] if column < len(consensus) else "",
        "conservation": scores[column] if column < len(scores) else float("nan"),
        "variation_class": labels[column] if column < len(labels) else "",
        "n_gaps": n_gap,
        "gap_fraction": n_gap / len(residues) if residues else float("nan"),
        "original_positions_0based": original_positions,
        "status": "COMPUTED",
    }


def viewer_window(
    result: Mapping[str, Any],
    *,
    row_start: int = 0,
    col_start: int = 0,
    n_rows: int = MAX_DISPLAY_ROWS,
    n_cols: int = MAX_DISPLAY_COLUMNS,
) -> dict:
    """Recorte paginado do MSA para o viewer. Nao reinterpreta strings cruas.

    Args:
        result: Envelope MSA.
        row_start: Primeira linha 0-based.
        col_start: Primeira coluna 0-based.
        n_rows, n_cols: Tamanho da janela, limitado aos tetos de display.

    Returns:
        Dict com slice de linhas, consenso e scores da janela.

    Raises:
        MsaError: INVALID_INPUT se o resultado nao for um MSA validado.
    """
    rows = list(result.get("rows") or [])
    length = int(result.get("alignment_length") or 0)
    if not rows or length <= 0:
        raise MsaError("MSA viewer requires a validated alignment.", "INVALID_INPUT")
    rows_n = max(1, min(int(n_rows), MAX_DISPLAY_ROWS, len(rows)))
    cols_n = max(1, min(int(n_cols), MAX_DISPLAY_COLUMNS, length))
    r0 = max(0, min(int(row_start), max(0, len(rows) - 1)))
    c0 = max(0, min(int(col_start), max(0, length - 1)))
    r1 = min(len(rows), r0 + rows_n)
    c1 = min(length, c0 + cols_n)
    consensus = str((result.get("consensus") or {}).get("sequence") or "")
    scores = list((result.get("conservation") or {}).get("scores") or [])
    slice_rows = []
    for row in rows[r0:r1]:
        aligned = str(row.get("aligned") or "")
        slice_rows.append(
            {
                "identifier": row.get("identifier"),
                "source": row.get("source"),
                "aligned_slice": aligned[c0:c1],
            }
        )
    return {
        "row_start": r0,
        "row_end": r1,
        "col_start": c0,
        "col_end": c1,
        "n_rows": len(rows),
        "alignment_length": length,
        "rows": slice_rows,
        "consensus_slice": consensus[c0:c1],
        "conservation_slice": scores[c0:c1],
    }


def spans_for_column(
    result: Mapping[str, Any],
    column: int,
    sequence_index: int,
) -> Optional[dict]:
    """Prepara um intervalo 0-based de um residuo para CRISPR/3D futuros.

    Args:
        result: Envelope MSA.
        column: Coluna alinhada.
        sequence_index: Indice da sequencia.

    Returns:
        scientific_checks.make_span de comprimento 1, ou None se a coluna for gap.

    Raises:
        MsaError: INVALID_INPUT se os indices forem invalidos.
    """
    maps = list(result.get("coordinate_maps") or [])
    members = list(result.get("input_members") or [])
    if sequence_index < 0 or sequence_index >= len(maps) or sequence_index >= len(members):
        raise MsaError("Sequence index is outside the MSA.", "INVALID_INPUT")
    detail = column_detail(result, column)
    pos = detail["original_positions_0based"][sequence_index]
    if pos is None:
        return None
    sequence = str(members[sequence_index].get("sequence") or "")
    return scientific_checks.make_span(
        start=int(pos),
        end=int(pos) + 1,
        sequence=sequence,
        strand="+",
        source="MSA column mapping",
        label=f"msa_col_{column}",
        kind="msa_residue",
        status="COMPUTED",
    )


def build_msa_result(
    *,
    raw_alignment: str,
    members: Sequence[Mapping[str, Any]],
    tool: str,
    tool_version: str,
    method: str,
    parameters: Mapping[str, Any],
    job_id: str = "",
    source: str = "",
    identical_groups: Optional[Sequence[Mapping[str, Any]]] = None,
) -> dict:
    """Parse + validacao + consenso/conservacao/variacao/mapping.

    Args:
        raw_alignment: FASTA alinhado da ferramenta.
        members: Entradas na ordem submetida.
        tool, tool_version, method, parameters, job_id, source: proveniencia.
        identical_groups: Duplicatas de hash ja detectadas.

    Returns:
        Envelope tipado. tool_version vazio permanece vazio (nao e inventado).

    Raises:
        MsaError: PARSING_ERROR ou INVALID_INPUT.
    """
    rows = parse_aligned_fasta(raw_alignment)
    validated = validate_aligned_rows(rows, members)
    ordered_rows = validated["rows"]
    aligned_strings = [str(row["aligned"]) for row in ordered_rows]
    molecule = str(members[0].get("molecule") or "")
    maps = column_coordinate_maps(ordered_rows)
    consensus = consensus_from_columns(aligned_strings, molecule)
    conservation = conservation_shannon(aligned_strings, molecule)
    column_classes = classify_column_variation(aligned_strings)
    identity = identity_matrix_from_msa(aligned_strings)
    variation = variation_vs_reference(aligned_strings, 0)
    decorated_rows = []
    for row, member, mapping in zip(ordered_rows, members, maps):
        decorated_rows.append(
            {
                "identifier": member.get("identifier"),
                "aligned": row["aligned"],
                "ungapped": row["ungapped"],
                "hash": member.get("hash"),
                "source": member.get("source"),
                "accession": member.get("accession"),
                "version": member.get("version"),
                "organism": member.get("organism"),
                "coordinate_map": mapping,
            }
        )
    envelope = provenance.analysis_envelope(
        module="MSA",
        payload={
            "tool": tool,
            "tool_version": tool_version or "",
            "method": method,
            "job_id": job_id,
            "alignment_length": validated["alignment_length"],
            "n_sequences": len(members),
            "molecule": molecule,
            "identical_groups": list(identical_groups or []),
        },
        status="COMPUTED",
        algorithm=f"{tool} {method}".strip(),
        parameters=dict(parameters),
        source=source or tool,
        input_identifier=",".join(str(item.get("hash") or "") for item in members),
    )
    envelope["status"] = "COMPLETED"
    envelope["tool"] = tool
    envelope["tool_version"] = tool_version or ""
    envelope["method"] = method
    envelope["parameters"] = dict(parameters)
    envelope["job_id"] = job_id
    envelope["molecule"] = molecule
    envelope["alignment_length"] = validated["alignment_length"]
    envelope["rows"] = decorated_rows
    envelope["input_members"] = [dict(item) for item in members]
    envelope["input_hashes"] = [str(item.get("hash") or "") for item in members]
    envelope["input_order"] = [str(item.get("identifier") or "") for item in members]
    envelope["coordinate_maps"] = maps
    envelope["consensus"] = consensus
    envelope["conservation"] = conservation
    envelope["column_classes"] = column_classes
    envelope["identity_matrix"] = identity
    envelope["variation_vs_reference"] = variation
    envelope["reference_index"] = 0
    envelope["identical_groups"] = list(identical_groups or [])
    envelope["raw_alignment"] = raw_alignment
    envelope["cache_status"] = "live"
    envelope["retrieved_at_utc"] = provenance.utc_now()
    envelope["disclaimer"] = (
        "Multiple sequence alignment is computational. Consensus is not a "
        "biological sequence. Conservation is not function. This MSA is not "
        "itself a phylogenetic tree; inference is a separate step. BLAST hits "
        "are not an MSA and not a 3D structure."
    )
    identity = hashlib.sha256()
    identity.update(str(molecule).encode("utf-8"))
    identity.update(b"\n")
    identity.update(str(validated["alignment_length"]).encode("utf-8"))
    identity.update(b"\n")
    for digest in envelope["input_hashes"]:
        identity.update(str(digest).encode("utf-8"))
        identity.update(b"\n")
    for row in decorated_rows:
        identity.update(str(row.get("identifier") or "").encode("utf-8"))
        identity.update(b"\n")
        identity.update(str(row.get("aligned") or "").encode("utf-8"))
        identity.update(b"\n")
        identity.update(str(row.get("hash") or "").encode("utf-8"))
        identity.update(b"\n")
    envelope["alignment_hash"] = identity.hexdigest()
    return envelope


def export_aligned_fasta(result: Mapping[str, Any]) -> str:
    """FASTA alinhado com comentario de proveniencia.

    Args:
        result: Envelope MSA.

    Returns:
        Texto FASTA.

    Raises:
        MsaError: INVALID_INPUT se nao houver linhas.
    """
    rows = list(result.get("rows") or [])
    if not rows:
        raise MsaError("No aligned rows to export.", "INVALID_INPUT")
    header = (
        f"; HelixScope MSA export\n"
        f"; tool={result.get('tool')} version={result.get('tool_version') or 'not reported'}\n"
        f"; method={result.get('method')}\n"
        f"; parameters={result.get('parameters')}\n"
        f"; timestamp={result.get('retrieved_at_utc')}\n"
        f"; input_hashes={','.join(result.get('input_hashes') or [])}\n"
    )
    chunks = [header]
    for row in rows:
        chunks.append(f">{row.get('identifier')}\n")
        aligned = str(row.get("aligned") or "")
        for i in range(0, len(aligned), 60):
            chunks.append(aligned[i : i + 60] + "\n")
    return "".join(chunks)


def export_clustal_like(result: Mapping[str, Any]) -> str:
    """Texto estilo CLUSTAL a partir do MSA validado, nao de um parser inverso.

    Args:
        result: Envelope MSA.

    Returns:
        Texto CLUSTAL com cabecalho de proveniencia.

    Raises:
        MsaError: INVALID_INPUT.
    """
    rows = list(result.get("rows") or [])
    if not rows:
        raise MsaError("No aligned rows to export.", "INVALID_INPUT")
    tool = str(result.get("tool") or "MSA")
    version = str(result.get("tool_version") or "version not reported")
    lines = [
        f"CLUSTAL {tool} ({version}) multiple sequence alignment",
        f"! source={result.get('method')}",
        f"! parameters={result.get('parameters')}",
        f"! timestamp={result.get('retrieved_at_utc')}",
        f"! input_hashes={','.join(result.get('input_hashes') or [])}",
        "",
    ]
    width = 60
    length = int(result.get("alignment_length") or 0)
    ids = [str(row.get("identifier") or "")[:15].ljust(16) for row in rows]
    aligned = [str(row.get("aligned") or "") for row in rows]
    for start in range(0, length, width):
        end = min(length, start + width)
        for label, seq in zip(ids, aligned):
            lines.append(f"{label}{seq[start:end]}")
        lines.append("")
    return "\n".join(lines) + "\n"


def export_column_csv(result: Mapping[str, Any]) -> str:
    """CSV por coluna: consenso, conservacao, gaps, classe.

    Args:
        result: Envelope MSA.

    Returns:
        CSV UTF-8.

    Raises:
        Nenhum.
    """
    import csv
    import io as io_module

    consensus = str((result.get("consensus") or {}).get("sequence") or "")
    scores = list((result.get("conservation") or {}).get("scores") or [])
    gaps = list((result.get("conservation") or {}).get("gap_fractions") or [])
    classes = list(result.get("column_classes") or [])
    length = int(result.get("alignment_length") or 0)
    buffer = io_module.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=[
            "column_0based",
            "consensus",
            "conservation_shannon",
            "gap_fraction",
            "variation_class",
            "method",
        ],
    )
    writer.writeheader()
    method = str((result.get("conservation") or {}).get("method") or "")
    for col in range(length):
        score = scores[col] if col < len(scores) else ""
        if isinstance(score, float) and math.isnan(score):
            score = "N/A"
        writer.writerow(
            {
                "column_0based": col,
                "consensus": consensus[col] if col < len(consensus) else "",
                "conservation_shannon": score,
                "gap_fraction": gaps[col] if col < len(gaps) else "",
                "variation_class": classes[col] if col < len(classes) else "",
                "method": method,
            }
        )
    return buffer.getvalue()


def export_result_bundle(result: Mapping[str, Any]) -> dict:
    """Copia JSON-safe do envelope, sem NaN/Inf.

    Args:
        result: Envelope MSA.

    Returns:
        Dict serializavel.

    Raises:
        Nenhum.
    """
    return provenance.json_safe(dict(result))


def mark_cached_result(result: Mapping[str, Any]) -> dict:
    """Marca resultado de cache sem fingir nova execucao.

    Args:
        result: Envelope MSA.

    Returns:
        Copia com cache_status=cached.

    Raises:
        Nenhum.
    """
    copied = dict(result)
    copied["cache_status"] = "cached"
    copied["served_from_cache_at_utc"] = provenance.utc_now()
    return copied


def _collection_fasta(members: Sequence[Mapping[str, Any]]) -> str:
    lines: List[str] = []
    for index, item in enumerate(members):
        identifier = str(item.get("identifier") or f"seq{index + 1}").replace(" ", "_")
        identifier = re.sub(r"[^A-Za-z0-9._-]+", "_", identifier)[:50] or f"seq{index + 1}"
        lines.append(f">{identifier}")
        seq = str(item.get("sequence") or "")
        for i in range(0, len(seq), 60):
            lines.append(seq[i : i + 60])
    return "\n".join(lines) + "\n"


def _stype(molecule: str) -> str:
    mapping = {"DNA": "dna", "RNA": "rna", "PROTEIN": "protein"}
    key = mapping.get(str(molecule or ""))
    if key is None:
        raise MsaError("Molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")
    return key


def _residue_alphabet(molecule: str) -> frozenset[str]:
    if molecule == "DNA":
        return dna_analysis.CANONICAL_DNA_ALPHABET
    if molecule == "RNA":
        return dna_analysis.CANONICAL_RNA_ALPHABET
    if molecule == "PROTEIN":
        return dna_analysis.CANONICAL_PROTEIN_ALPHABET
    raise MsaError("Molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")


def _normalize_gaps(aligned: str) -> str:
    return "".join("-" if char in GAP_CHARS else char for char in aligned.upper())


def _tie_symbol(winners: Sequence[str], molecule: str) -> str:
    if molecule == "DNA":
        return DNA_IUPAC_FROM_SET.get(frozenset(winners), "N")
    if molecule == "RNA":
        return "N"
    return "X"


def _match_rows_to_members(
    rows: Sequence[Mapping[str, Any]],
    members: Sequence[Mapping[str, Any]],
) -> List[dict]:
    """Associa linhas do FASTA alinhado as entradas por ungapped, depois por ordem."""
    remaining = [dict(row) for row in rows]
    ordered: List[dict] = []
    for member in members:
        original = str(member.get("sequence") or "").upper()
        match_index = next(
            (
                index
                for index, row in enumerate(remaining)
                if str(row.get("ungapped") or "").upper() == original
            ),
            None,
        )
        if match_index is None:
            raise MsaError(
                f"No aligned row matches input sequence {member.get('identifier')}.",
                "PARSING_ERROR",
            )
        ordered.append(remaining.pop(match_index))
    if remaining:
        raise MsaError("Aligned FASTA has extra rows not present in the input.", "PARSING_ERROR")
    return ordered


def _require_email(email: str) -> None:
    from . import ncbi_fetch

    if not ncbi_fetch._valid_entrez_email(email or ""):
        raise ValueError(
            "Informe um e-mail de contato valido para a EMBL-EBI "
            "(formato local@dominio)."
        )


def _ebi_submit(
    fasta: str,
    email: str,
    parameters: Mapping[str, Any],
    *,
    urlopen_fn=None,
) -> str:
    payload = {
        "email": email.strip(),
        "sequence": fasta,
        "stype": str(parameters.get("stype") or "dna"),
        "outfmt": "fa",
        "order": "input",
        "title": EBI_TOOL,
    }
    body = _ebi_request(
        f"{EBI_CLUSTALO_BASE}/run",
        method="POST",
        data=payload,
        urlopen_fn=urlopen_fn,
    )
    job_id = body.strip()
    if not JOB_ID_PATTERN.fullmatch(job_id):
        raise MsaError(
            "EMBL-EBI did not return a valid Clustal Omega job id.",
            "JOB_FAILED",
        )
    return job_id


def _ebi_status(job_id: str, *, urlopen_fn=None) -> str:
    body = _ebi_request(
        f"{EBI_CLUSTALO_BASE}/status/{job_id}",
        method="GET",
        urlopen_fn=urlopen_fn,
    )
    return body.strip().upper()


def _ebi_result_aligned(job_id: str, *, urlopen_fn=None) -> str:
    """Busca FASTA alinhado nos tipos documentados do Job Dispatcher.

    Tenta aln-fasta, depois fa e out. Nao interpreta HSP nem texto de log
    como alinhamento.
    """
    last_error: Optional[MsaError] = None
    for kind in ("aln-fasta", "fa", "out"):
        try:
            body = _ebi_result(job_id, kind, urlopen_fn=urlopen_fn)
        except MsaError as exc:
            last_error = exc
            if exc.category == "INVALID_INPUT":
                continue
            raise
        if body.strip().startswith(">"):
            return body
        last_error = MsaError(
            f"EMBL-EBI result type {kind} was not aligned FASTA.",
            "PARSING_ERROR",
        )
    if last_error is not None:
        raise MsaError(
            "EMBL-EBI did not return aligned FASTA (aln-fasta/fa/out).",
            "PARSING_ERROR",
        ) from last_error
    raise MsaError("EMBL-EBI did not return aligned FASTA.", "PARSING_ERROR")


def _ebi_result(job_id: str, result_type: str, *, urlopen_fn=None) -> str:
    allowed = {"fa", "aln-fasta", "out"}
    kind = str(result_type or "fa")
    if kind not in allowed:
        raise MsaError("Unsupported EBI result type.", "INVALID_INPUT")
    return _ebi_request(
        f"{EBI_CLUSTALO_BASE}/result/{job_id}/{kind}",
        method="GET",
        urlopen_fn=urlopen_fn,
    )


def _ebi_request(
    url: str,
    *,
    method: str,
    data: Optional[Mapping[str, str]] = None,
    urlopen_fn=None,
) -> str:
    if not url.startswith(EBI_CLUSTALO_BASE + "/") and url != EBI_CLUSTALO_BASE:
        raise MsaError("Refusing a non-allowlisted MSA URL.", "INVALID_INPUT")
    opener = urlopen_fn or urlopen
    encoded = None
    headers = {"User-Agent": f"{EBI_TOOL}/msa"}
    if data is not None:
        encoded = urlencode(dict(data)).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = Request(url, data=encoded, method=method, headers=headers)
    try:
        with opener(request, timeout=HTTP_TIMEOUT_S) as handle:
            raw = handle.read()
    except HTTPError as exc:
        raise _http_error(exc) from exc
    except (URLError, TimeoutError, socket.timeout, OSError) as exc:
        message = str(exc).lower()
        if "timed out" in message or isinstance(exc, (TimeoutError, socket.timeout)):
            raise MsaError(
                "EMBL-EBI MSA request timed out. This is TIMEOUT, not UNAVAILABLE.",
                "TIMEOUT",
            ) from exc
        raise MsaError(
            f"EMBL-EBI MSA unavailable: {exc}. This is not UNAVAILABLE-as-success.",
            "SERVICE_UNAVAILABLE",
        ) from exc
    if isinstance(raw, bytes):
        if len(raw) > MAX_OUTPUT_BYTES:
            raise MsaError(
                f"MSA response exceeds {MAX_OUTPUT_BYTES:,} bytes.",
                "RESOURCE_LIMIT",
            )
        return raw.decode("utf-8", errors="replace")
    text = str(raw)
    if len(text.encode("utf-8")) > MAX_OUTPUT_BYTES:
        raise MsaError(
            f"MSA response exceeds {MAX_OUTPUT_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    return text


def _http_error(exc: HTTPError) -> MsaError:
    code = int(getattr(exc, "code", 0) or 0)
    if code == 429:
        return MsaError(
            "EMBL-EBI rate limited the MSA request (HTTP 429).",
            "RATE_LIMITED",
        )
    if code >= 500:
        return MsaError(
            f"EMBL-EBI MSA unavailable (HTTP {code}).",
            "SERVICE_UNAVAILABLE",
        )
    if code in {400, 404, 415}:
        return MsaError(f"EMBL-EBI rejected the MSA request (HTTP {code}).", "INVALID_INPUT")
    return MsaError(f"EMBL-EBI MSA HTTP {code}: {exc}.", "SERVICE_UNAVAILABLE")


def _normalized_tool_basename(path: str) -> str:
    text = str(path or "").replace("\\", "/")
    base = os.path.basename(text).lower()
    for suffix in (".exe", ".bat", ".cmd"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def _allowlisted_basenames(names: Sequence[str]) -> set[str]:
    allowed: set[str] = set()
    for item in names:
        low = str(item).lower()
        allowed.add(low)
        allowed.add(_normalized_tool_basename(low))
    return allowed


def _resolve_allowlisted_executable(names: Sequence[str]) -> Optional[str]:
    from . import tool_detection

    return tool_detection.resolve_allowlisted_executable(
        names,
        which_fn=shutil.which,
    )


def _probe_local_version(executable: str, tool_id: str) -> str:
    args = [executable]
    if tool_id == "muscle":
        args.append("-version")
    else:
        args.append("--version")
    try:
        completed = subprocess.run(
            args,
            capture_output=True,
            timeout=5,
            check=False,
            shell=False,
            text=True,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = (completed.stdout or "") + (completed.stderr or "")
    line = text.strip().splitlines()[0] if text.strip() else ""
    return line[:120]


def _run_local_backend(
    validated: Mapping[str, Any],
    backend_id: str,
    parameters: Mapping[str, Any],
) -> dict:
    import time as time_module

    tool_id = backend_id.replace("local_", "", 1)
    local = detect_local_tools().get(tool_id)
    if not local or not local.get("available"):
        raise MsaError(
            f"{tool_id} executable not installed.",
            "TOOL_NOT_INSTALLED",
        )
    executable = str(local.get("path") or "")
    fasta = _collection_fasta(list(validated["members"]))
    raw, version = _execute_local(executable, tool_id, fasta, str(validated["molecule"]))
    result = build_msa_result(
        raw_alignment=raw,
        members=list(validated["members"]),
        tool=tool_id,
        tool_version=version or str(local.get("version") or ""),
        method="local executable",
        parameters=dict(parameters),
        job_id="",
        source=f"local {tool_id}",
        identical_groups=list(validated.get("identical_groups") or []),
    )
    return {
        "status": "COMPLETED",
        "backend": backend_id,
        "tool": tool_id,
        "tool_version": result.get("tool_version") or "",
        "job_id": "",
        "parameters": dict(parameters),
        "molecule": validated["molecule"],
        "n_sequences": validated["n_sequences"],
        "input_hashes": validated["hashes"],
        "identifiers": validated["identifiers"],
        "identical_groups": validated["identical_groups"],
        "members": validated["members"],
        "submitted_at_utc": provenance.utc_now(),
        "submitted_monotonic": time_module.monotonic(),
        "last_poll_monotonic": time_module.monotonic(),
        "source": f"local {tool_id}",
        "cache_status": "live",
        "software_version": provenance.HELIXSCOPE_VERSION,
        "result": result,
    }


def _execute_local(
    executable: str, tool_id: str, fasta: str, molecule: str
) -> Tuple[str, str]:
    version = _probe_local_version(executable, tool_id)
    tmpdir = tempfile.mkdtemp(prefix="helixscope_msa_")
    try:
        infile = os.path.join(tmpdir, "input.fa")
        outfile = os.path.join(tmpdir, "output.fa")
        with open(infile, "w", encoding="utf-8") as handle:
            handle.write(fasta)
        args = _local_args(executable, tool_id, infile, outfile, molecule, version=version)
        try:
            completed = subprocess.run(
                args,
                capture_output=True,
                timeout=LOCAL_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmpdir,
                text=True,
            )
        except subprocess.TimeoutExpired as exc:
            raise MsaError(
                f"Local {tool_id} exceeded {int(LOCAL_TIMEOUT_S)}s. This is TIMEOUT.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise MsaError(f"Local {tool_id} failed to start: {exc}.", "TOOL_FAILED") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()[:300]
            raise MsaError(
                f"Local {tool_id} exited {completed.returncode}. {detail}",
                "TOOL_FAILED",
            )
        if os.path.isfile(outfile):
            with open(outfile, encoding="utf-8") as handle:
                raw = handle.read()
        else:
            raw = completed.stdout or ""
        if not raw.strip():
            raise MsaError(f"Local {tool_id} produced empty alignment output.", "PARSING_ERROR")
        if len(raw.encode("utf-8")) > MAX_OUTPUT_BYTES:
            raise MsaError(
                f"Local {tool_id} output exceeds {MAX_OUTPUT_BYTES:,} bytes.",
                "RESOURCE_LIMIT",
            )
        return raw, version
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _local_args(
    executable: str,
    tool_id: str,
    infile: str,
    outfile: str,
    molecule: str,
    version: str = "",
) -> List[str]:
    seqtype = {"DNA": "DNA", "RNA": "RNA", "PROTEIN": "Protein"}[molecule]
    if tool_id == "clustalo":
        return [
            executable,
            f"--infile={infile}",
            f"--outfile={outfile}",
            "--outfmt=fa",
            "--force",
            "--output-order=input-order",
            f"--seqtype={seqtype}",
        ]
    if tool_id == "mafft":
        return [executable, "--auto", "--inputorder", infile]
    if tool_id == "muscle":
        text = str(version or "").lower()
        if re.search(r"\b3\.\d", text) or "muscle 3" in text:
            return [executable, "-in", infile, "-out", outfile]
        return [executable, "-align", infile, "-output", outfile]
    raise MsaError(f"Local tool '{tool_id}' is not allowlisted.", "INVALID_INPUT")
