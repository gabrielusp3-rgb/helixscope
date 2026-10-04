"""Cliente do NCBI BLAST Common URL API (QBLAST / Blast.cgi).

Submete buscas reais ao servico BLAST do NCBI, acompanha o RID e interpreta
BlastOutput XML. BLAST+ local corre apenas contra o prefixo declarado em
HELIXSCOPE_BLAST_DB. Timeout remoto nao troca para local.

Nenhuma funcao aqui importa Streamlit. O destino HTTP e fixo: nao ha URL
arbitraria do usuario (protecao SSRF).
"""

from __future__ import annotations

import math
import os
import re
import socket
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from typing import Any, Dict, List, Mapping, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from . import dna_analysis, engine_validation, provenance, scientific_checks, tool_detection, tool_paths

BLAST_ENDPOINT: str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
"""Unico destino HTTP permitido para BLAST. Nao aceitar URL do usuario."""

BLAST_TOOL: str = "HelixScope"
"""Parametro tool enviado ao NCBI, conforme a politica de desenvolvedores."""

BLAST_HTTP_TIMEOUT_S: float = 45.0
"""Timeout de cada requisicao HTTP ao Blast.cgi."""

MIN_CONTACT_INTERVAL_S: float = 10.0
"""Intervalo minimo entre qualquer contato com o servidor BLAST (politica NCBI)."""

RID_POLL_INTERVAL_S: float = 60.0
"""Intervalo minimo entre polls do mesmo RID apos o primeiro Get."""

MAX_JOB_WAIT_S: float = 600.0
"""Tempo maximo desde o Put ate declarar TIMEOUT, sem converter em NO_HITS."""

MAX_QUERY_LENGTH: int = 8_000
"""Teto de residuos na query. Nao e um limite do NCBI; protege o processo."""

MAX_HITLIST_SIZE: int = 50
"""HITLIST_SIZE maximo pedido ao NCBI."""

DEFAULT_HITLIST_SIZE: int = 20
"""HITLIST_SIZE padrao."""

MAX_XML_BYTES: int = 2_000_000
"""Teto do XML de resposta; respostas maiores sao RESOURCE_LIMIT."""

MAX_DISPLAY_HITS: int = 10
"""Hits por pagina na interface; o XML pode conter ate HITLIST_SIZE."""

MAX_SESSION_CACHE: int = 8
"""Numero maximo de resultados BLAST guardados na sessao do usuario."""

MAX_CONCURRENT_JOBS: int = 1
"""Um job BLAST por sessao. Sem fila distribuida nesta fase."""

RID_PATTERN: re.Pattern[str] = re.compile(r"^[A-Za-z0-9]{8,32}$")
"""RID devolvido pelo NCBI. Recusa qualquer coisa que pareca URL ou comando."""

BLAST_PROGRAMS: Dict[str, Dict[str, Tuple[str, ...]]] = {
    "blastn": {
        "query_molecules": ("DNA", "RNA"),
        "databases": ("core_nt", "nt"),
    },
    "blastp": {
        "query_molecules": ("PROTEIN",),
        "databases": ("swissprot", "nr", "refseq_protein"),
    },
    "blastx": {
        "query_molecules": ("DNA", "RNA"),
        "databases": ("swissprot", "nr", "refseq_protein"),
    },
    "tblastn": {
        "query_molecules": ("PROTEIN",),
        "databases": ("core_nt", "nt"),
    },
    "tblastx": {
        "query_molecules": ("DNA", "RNA"),
        "databases": ("core_nt", "nt"),
    },
}
"""Programas e bancos que este cliente realmente submete. Sem opcoes ornamentais."""

DEFAULT_DATABASE: Dict[str, str] = {
    "blastn": "core_nt",
    "blastp": "swissprot",
    "blastx": "swissprot",
    "tblastn": "core_nt",
    "tblastx": "core_nt",
}
"""Banco padrao por programa: core_nt e swissprot sao nomes oficiais da URL API."""

BLASTPLUS_PROGRAMS: tuple[str, ...] = (
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
)
BLASTPLUS_ENV: Dict[str, str] = {
    "blastn": "HELIXSCOPE_BLASTN",
    "blastp": "HELIXSCOPE_BLASTP",
    "blastx": "HELIXSCOPE_BLASTX",
    "tblastn": "HELIXSCOPE_TBLASTN",
    "tblastx": "HELIXSCOPE_TBLASTX",
}
BLAST_DB_ENV: str = "HELIXSCOPE_BLAST_DB"
MAKEBLASTDB_ENV: str = "HELIXSCOPE_MAKEBLASTDB"
DECLARED_NUCL_DB_NAME: str = "helixscope_tiny_nucl"
DECLARED_PROT_DB_NAME: str = "helixscope_tiny_prot"
LOCAL_BLAST_TIMEOUT_S: float = 30.0
LOCAL_BLAST_MAX_OUTPUT_BYTES: int = 2_000_000
BLASTPLUS_VERSION_RE: re.Pattern[str] = re.compile(
    r"blast(?:n|p|x|tblastn|tblastx)?[:\s]+([0-9][0-9.\w+]*)",
    re.I,
)


class BlastError(RuntimeError):
    """Falha classificada de BLAST. Nunca deve ser exibida como NO_HITS.

    Attributes:
        category: TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE, INVALID_INPUT,
            INVALID_DATABASE, JOB_FAILED, PARSING_ERROR ou RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "JOB_FAILED")


def _blastplus_names(program: str) -> tuple[str, ...]:
    return (program, f"{program}.exe")


def _detect_blastplus_program(program: str) -> dict:
    executable = tool_detection.resolve_allowlisted_executable(
        _blastplus_names(program),
        extra_file_candidates=tool_paths.candidates_for(_blastplus_names(program)),
        env_var=BLASTPLUS_ENV[program],
    )
    if executable is None:
        return {
            "available": False,
            "program": program,
            "path": "",
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
        }
    version = ""
    try:
        completed = subprocess.run(
            [executable, "-version"],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
        blob = f"{completed.stdout or ''}\n{completed.stderr or ''}"
        match = BLASTPLUS_VERSION_RE.search(blob)
        if match:
            version = match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        version = ""
    return {
        "available": True,
        "program": program,
        "path": tool_detection.sanitize_tool_path(executable),
        "version": version,
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
    }


def declared_blast_prefix(program: str = "blastn") -> str:
    """Prefixo do banco local: env HELIXSCOPE_BLAST_DB ou fixture do projeto.

    Args:
        program: blastn/blastp/... para escolher nucl vs prot quando o env
            nao estiver definido.

    Returns:
        Caminho absoluto do prefixo, ou string vazia.

    Raises:
        Nenhum.
    """
    env = (os.environ.get(BLAST_DB_ENV) or "").strip().strip('"')
    if env:
        return os.path.abspath(env)
    tools = tool_paths.repository_tools_dir()
    name = str(program or "blastn").strip().lower()
    folder = DECLARED_PROT_DB_NAME if name in {"blastp", "blastx"} else DECLARED_NUCL_DB_NAME
    candidate = os.path.join(tools, "blast_db", folder)
    if _local_blast_database_record_for_prefix(candidate).get("available"):
        return os.path.abspath(candidate)
    return ""


def _local_blast_database_record() -> dict:
    prefix = (os.environ.get(BLAST_DB_ENV) or "").strip().strip('"')
    if prefix:
        resolved = os.path.abspath(prefix)
        suffixes = (".nin", ".nal", ".pin", ".pal", ".nsq", ".psq")
        if any(os.path.isfile(resolved + suffix) for suffix in suffixes):
            nsq = resolved + ".nsq" if os.path.isfile(resolved + ".nsq") else resolved + ".psq"
            identity = os.path.basename(resolved)
            try:
                size = os.path.getsize(nsq) if os.path.isfile(nsq) else 0
                identity = f"{os.path.basename(resolved)}:{size}"
            except OSError:
                identity = os.path.basename(resolved)
            return {
                "available": True,
                "prefix": os.path.basename(resolved),
                "identity": identity,
                "reason": "Declared BLAST+ database prefix exists as a local index file.",
            }
        return {
            "available": False,
            "prefix": os.path.basename(resolved),
            "reason": (
                "HELIXSCOPE_BLAST_DB is set but no BLAST index file was found next "
                "to that prefix (.nin/.nal/.pin/.pal/.nsq/.psq). The path is not executed."
            ),
        }
    nucl = declared_blast_prefix("blastn")
    prot = declared_blast_prefix("blastp")
    if nucl or prot:
        return {
            "available": True,
            "prefix": "helixscope_tiny",
            "identity": "helixscope_tiny_declared_fixture",
            "reason": (
                "Project-local declared BLAST fixture under tools/blast_db/. "
                "Not NCBI nt/nr."
            ),
            "not_nt": True,
            "not_nr": True,
        }
    return {
        "available": False,
        "prefix": "",
        "reason": (
            f"{BLAST_DB_ENV} is unset and no tools/blast_db fixture is present. "
            "BLAST+ binaries without a declared database prefix cannot run a local "
            "search. Databases are not downloaded automatically."
        ),
    }


def local_blast_availability() -> dict:
    """Detecta BLAST+ NCBI oficial no PATH/env. Nao descarrega bancos.

    Args:
        Nenhum.

    Returns:
        Dict executables_detected, programs, database, available (True so se
        binario E banco existirem), version, reason.

    Raises:
        Nenhum.

    Nota biologica:
        Binario sem banco nao e uma busca BLAST. NCBI remoto permanece o
        backend de homologia desta aplicacao ate o utilizador escolher
        explicitamente um banco local.
    """
    programs = {name: _detect_blastplus_program(name) for name in BLASTPLUS_PROGRAMS}
    detected = [row for row in programs.values() if row["available"]]
    database = _local_blast_database_record()
    version = ""
    for name in BLASTPLUS_PROGRAMS:
        version = str(programs[name].get("version") or "")
        if version:
            break
    executables_detected = bool(detected)
    can_search = executables_detected and bool(database.get("available"))
    if not executables_detected:
        reason = (
            "NCBI BLAST+ executables were not found on PATH or HELIXSCOPE_BLASTN. "
            "No recursive disk scan is performed. Remote NCBI BLAST remains available."
        )
        status = tool_detection.TOOL_STATUS_NOT_INSTALLED
    elif not database.get("available"):
        reason = (
            "BLAST+ executable detected; local database is UNAVAILABLE. "
            + str(database.get("reason") or "")
        )
        status = tool_detection.classify_tool_record(
            available=True,
            path=str(detected[0].get("path") or "blastn"),
            version=version,
        )
    else:
        reason = "BLAST+ executable and declared local database prefix are present."
        status = tool_detection.classify_tool_record(
            available=True,
            path=str(detected[0].get("path") or "blastn"),
            version=version,
        )
    return {
        "available": can_search,
        "executables_detected": executables_detected,
        "version": version,
        "status": status,
        "programs": programs,
        "database": database,
        "reason": reason,
        "backend": "NCBI BLAST+",
        "remote_unchanged": True,
    }


def detect_makeblastdb() -> dict:
    """Detecta makeblastdb allowlisted. Nao cria bancos automaticamente.

    Args:
        Nenhum.

    Returns:
        Dict available, version, path sanitizado.

    Raises:
        Nenhum.
    """
    names = ("makeblastdb", "makeblastdb.exe")
    executable = tool_detection.resolve_allowlisted_executable(
        names,
        extra_file_candidates=tool_paths.candidates_for(names),
        env_var=MAKEBLASTDB_ENV,
    )
    if executable is None:
        return {
            "available": False,
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "reason": "makeblastdb is not installed. HelixScope does not download nt/nr.",
        }
    version = ""
    try:
        completed = subprocess.run(
            [executable, "-version"],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
        blob = f"{completed.stdout or ''}\n{completed.stderr or ''}"
        match = BLASTPLUS_VERSION_RE.search(blob)
        if match:
            version = match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        version = ""
    return {
        "available": True,
        "version": version,
        "path": tool_detection.sanitize_tool_path(executable),
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
        "reason": "makeblastdb detected. Not run at import.",
    }


def _subprocess_exited_nonzero(completed: object) -> bool:
    """True when a subprocess did not exit 0.

    Args:
        completed: Result of subprocess.run.

    Returns:
        True if returncode is missing or not 0. Exit 0 is success: `0 or 1`
        must not be used, because bool(0) is False.

    Raises:
        Nenhum.
    """
    code = getattr(completed, "returncode", None)
    if code is None:
        return True
    try:
        return int(code) != 0
    except (TypeError, ValueError):
        return True


def build_declared_blast_database(
    *,
    fasta_text: str,
    dbtype: str,
    prefix: str,
    title: str = "helixscope_declared_db",
) -> dict:
    """Cria um banco BLAST+ pequeno a partir de FASTA declarado (testes/fixtures).

    Args:
        fasta_text: FASTA de sequencias conhecidas, nao nt/nr.
        dbtype: nucl ou prot.
        prefix: Prefixo absoluto permitido (diretorio ja existente + nome).
        title: Titulo gravado no banco.

    Returns:
        Dict available, identity, version.

    Raises:
        BlastError: UNAVAILABLE, INVALID_INPUT, TOOL_ERROR.
    """
    kind = str(dbtype or "").strip().lower()
    if kind not in {"nucl", "prot"}:
        raise BlastError("dbtype must be nucl or prot.", "INVALID_INPUT")
    text = str(fasta_text or "").strip()
    if not text.startswith(">") or len(text.encode("utf-8")) > 200_000:
        raise BlastError("Declared BLAST FASTA is empty or exceeds 200 KiB.", "INVALID_INPUT")
    info = detect_makeblastdb()
    if not info.get("available"):
        raise BlastError(str(info.get("reason") or "makeblastdb UNAVAILABLE"), "UNAVAILABLE")
    names = ("makeblastdb", "makeblastdb.exe")
    executable = tool_detection.resolve_allowlisted_executable(
        names,
        extra_file_candidates=tool_paths.candidates_for(names),
        env_var=MAKEBLASTDB_ENV,
    )
    if not executable:
        raise BlastError("makeblastdb UNAVAILABLE.", "UNAVAILABLE")
    resolved = os.path.abspath(str(prefix or ""))
    parent = os.path.dirname(resolved)
    if not parent or not os.path.isdir(parent):
        raise BlastError("BLAST database prefix directory does not exist.", "INVALID_INPUT")
    fasta_path = resolved + ".fa"
    with open(fasta_path, "w", encoding="utf-8") as handle:
        handle.write(text if text.endswith("\n") else text + "\n")
    try:
        completed = subprocess.run(
            [
                executable,
                "-in",
                fasta_path,
                "-dbtype",
                kind,
                "-out",
                resolved,
                "-title",
                str(title or "helixscope_declared_db"),
                "-parse_seqids",
            ],
            capture_output=True,
            text=True,
            timeout=LOCAL_BLAST_TIMEOUT_S,
            check=False,
            shell=False,
            cwd=parent,
        )
    except subprocess.TimeoutExpired as exc:
        raise BlastError("makeblastdb exceeded the timeout.", "TIMEOUT") from exc
    except OSError as exc:
        raise BlastError(f"makeblastdb could not run: {exc}", "TOOL_ERROR") from exc
    if _subprocess_exited_nonzero(completed):
        err = str(getattr(completed, "stderr", "") or "")[:300]
        raise BlastError(f"makeblastdb failed: {err or 'nonzero exit'}", "TOOL_ERROR")
    record = _local_blast_database_record_for_prefix(resolved)
    if not record.get("available"):
        raise BlastError("makeblastdb finished but no BLAST index file was found.", "PARSING_ERROR")
    return {
        "available": True,
        "prefix": os.path.basename(resolved),
        "identity": record.get("identity"),
        "dbtype": kind,
        "title": title,
        "tool": "makeblastdb",
        "version": str(info.get("version") or ""),
        "not_nt": True,
        "not_nr": True,
        "reason": "Declared fixture/test BLAST database. Not NCBI nt/nr/refseq.",
    }


def _local_blast_database_record_for_prefix(resolved: str) -> dict:
    suffixes = (".nin", ".nal", ".pin", ".pal", ".nsq", ".psq")
    if any(os.path.isfile(resolved + suffix) for suffix in suffixes):
        nsq = resolved + ".nsq" if os.path.isfile(resolved + ".nsq") else resolved + ".psq"
        identity = os.path.basename(resolved)
        try:
            size = os.path.getsize(nsq) if os.path.isfile(nsq) else 0
            identity = f"{os.path.basename(resolved)}:{size}"
        except OSError:
            identity = os.path.basename(resolved)
        return {"available": True, "prefix": os.path.basename(resolved), "identity": identity}
    return {"available": False, "prefix": os.path.basename(resolved), "identity": ""}


def run_local_blast(
    *,
    program: str,
    query: str,
    molecule: str,
    timeout_s: float = LOCAL_BLAST_TIMEOUT_S,
) -> dict:
    """Executa BLAST+ local contra o prefixo HELIXSCOPE_BLAST_DB.

    Args:
        program: blastn/blastp/blastx/tblastn/tblastx.
        query: Sequencia da query.
        molecule: DNA, RNA ou PROTEIN.
        timeout_s: Timeout do subprocesso.

    Returns:
        Envelope parseado (mesmo contrato XML que o BLAST remoto) com backend local.

    Raises:
        BlastError: UNAVAILABLE, INVALID_INPUT, TIMEOUT, PARSING_ERROR, RESOURCE_LIMIT.

    Nota biologica:
        Hits, bit scores e E-values vem do BLAST+ oficial contra o banco
        declarado. O banco de teste nao e nt/nr. Timeout remoto NCBI nao chama
        esta funcao. blastn usa o filtro DUST padrao do BLAST+; queries de
        baixa complexidade podem devolver NO_HITS mesmo com identidade
        literal no banco. Isso nao e falha do parser.
    """
    name = str(program or "").strip().lower()
    if name not in BLASTPLUS_PROGRAMS:
        raise BlastError(f"BLAST+ program '{program}' is not offered.", "INVALID_INPUT")
    allowed = programs_for_molecule(molecule)
    if name not in allowed:
        raise BlastError(
            f"{name} is not valid for molecule {molecule}.",
            "INVALID_INPUT",
        )
    validated_query = dna_analysis.validate_for_molecule(query, molecule)
    if not validated_query.get("is_valid"):
        raise BlastError(
            str(validated_query.get("rejection_reason") or "BLAST query is invalid."),
            "INVALID_INPUT",
        )
    query_seq = str(validated_query.get("sequence") or "")
    if len(query_seq) > MAX_QUERY_LENGTH:
        raise BlastError(
            f"Query exceeds {MAX_QUERY_LENGTH} residues.",
            "RESOURCE_LIMIT",
        )
    avail = local_blast_availability()
    if not avail.get("executables_detected"):
        raise BlastError(str(avail.get("reason") or "BLAST+ UNAVAILABLE"), "UNAVAILABLE")
    database = avail.get("database") or {}
    if not database.get("available"):
        raise BlastError(str(database.get("reason") or "Local BLAST database UNAVAILABLE"), "UNAVAILABLE")
    prefix = declared_blast_prefix(name)
    if not prefix:
        raise BlastError(f"{BLAST_DB_ENV} is unset and no declared local BLAST fixture was found.", "UNAVAILABLE")
    resolved = os.path.abspath(prefix)
    if not _local_blast_database_record_for_prefix(resolved).get("available"):
        raise BlastError("Declared BLAST+ database prefix is invalid.", "INVALID_DATABASE")
    program_info = (avail.get("programs") or {}).get(name) or {}
    executable = tool_detection.resolve_allowlisted_executable(
        _blastplus_names(name),
        extra_file_candidates=tool_paths.candidates_for(_blastplus_names(name)),
        env_var=BLASTPLUS_ENV[name],
    )
    if not executable:
        raise BlastError(f"{name} executable UNAVAILABLE.", "UNAVAILABLE")
    try:
        timeout = float(timeout_s)
    except (TypeError, ValueError) as exc:
        raise BlastError("Local BLAST timeout must be numeric.", "INVALID_INPUT") from exc
    if not math.isfinite(timeout) or timeout <= 0:
        raise BlastError("Local BLAST timeout must be positive.", "INVALID_INPUT")
    with tempfile.TemporaryDirectory(prefix="helixscope_blast_") as tmp:
        query_path = os.path.join(tmp, "query.fa")
        out_path = os.path.join(tmp, "out.xml")
        with open(query_path, "w", encoding="utf-8") as handle:
            handle.write(f">query\n{query_seq}\n")
        argv = [
            executable,
            "-query",
            query_path,
            "-db",
            resolved,
            "-outfmt",
            "5",
            "-out",
            out_path,
            "-max_target_seqs",
            str(MAX_HITLIST_SIZE),
        ]
        try:
            completed = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
                shell=False,
                cwd=tmp,
            )
        except subprocess.TimeoutExpired as exc:
            raise BlastError(
                f"Local BLAST+ exceeded {timeout:.0f}s and was stopped.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise BlastError(f"BLAST+ could not be executed: {exc}", "TOOL_ERROR") from exc
        if _subprocess_exited_nonzero(completed):
            err = str(getattr(completed, "stderr", "") or "")[:400]
            raise BlastError(f"BLAST+ failed: {err or 'nonzero exit'}", "TOOL_ERROR")
        if not os.path.isfile(out_path):
            raise BlastError("BLAST+ did not write XML output.", "PARSING_ERROR")
        size = os.path.getsize(out_path)
        if size > LOCAL_BLAST_MAX_OUTPUT_BYTES:
            raise BlastError("BLAST+ XML exceeds the output cap.", "RESOURCE_LIMIT")
        with open(out_path, "r", encoding="utf-8", errors="replace") as handle:
            xml_body = handle.read()
    parsed = parse_blast_xml(xml_body)
    hits = list(parsed.get("hits") or [])
    skipped = int(parsed.get("skipped_hsps") or 0)
    if skipped and not hits:
        raise BlastError(
            "BLAST+ XML was parsed but every HSP failed validation. "
            "This is PARSING_ERROR, not NO_HITS.",
            "PARSING_ERROR",
        )
    status = "NO_HITS" if not hits else "READY"
    envelope = provenance.analysis_envelope(
        module="BLAST",
        payload={
            "program": parsed.get("program") or name,
            "database": parsed.get("database") or os.path.basename(resolved),
            "n_hits": len(hits),
            "backend": "NCBI BLAST+ local",
        },
        status="RETRIEVED" if hits else "UNAVAILABLE",
        algorithm="NCBI BLAST+",
        parameters={
            "program": name,
            "database_identity": str(database.get("identity") or ""),
            "format": "XML outfmt 5",
            "tool_version": str(program_info.get("version") or avail.get("version") or ""),
            "complexity_filter": "BLAST+ default (DUST for blastn, SEG for blastp)",
        },
        sequence="",
        source="NCBI BLAST+ local",
        accession="",
        input_identifier="",
    )
    envelope["hits"] = hits
    envelope["status"] = status
    envelope["n_hits"] = len(hits)
    envelope["blast_version"] = parsed.get("blast_version") or ""
    envelope["program"] = parsed.get("program") or name
    envelope["database"] = parsed.get("database") or os.path.basename(resolved)
    envelope["query_length"] = parsed.get("query_length")
    envelope["parse_warnings"] = parsed.get("warnings") or []
    envelope["backend"] = "NCBI BLAST+ local"
    envelope["database_identity"] = str(database.get("identity") or database.get("prefix") or "")
    envelope["tool_version"] = str(program_info.get("version") or avail.get("version") or "")
    envelope["not_remote_ncbi"] = True
    envelope["rid"] = ""
    engine_validation.record_live_validation(
        "BLAST+",
        ok=True,
        version=str(envelope.get("tool_version") or ""),
        details={
            "program": name,
            "database_identity": envelope["database_identity"],
            "n_hits": len(hits),
        },
    )
    return envelope


def ncbi_blast_availability() -> dict:
    """Declara que a homologia usa o NCBI BLAST Common URL API.

    Args:
        Nenhum.

    Returns:
        available True para o servico NCBI remoto; local BLAST+ aninhado.

    Raises:
        Nenhum.

    Nota biologica:
        Os hits, bit scores e E-values vem do BLAST executado nos servidores
        NCBI contra o banco pedido. HelixScope nao reproduz o algoritmo BLAST.
        Timeout remoto nao troca silenciosamente para BLAST+ local.
    """
    local = local_blast_availability()
    return {
        "available": True,
        "backend": "NCBI BLAST Common URL API",
        "endpoint": BLAST_ENDPOINT,
        "reason": (
            "Homology search submits real jobs to NCBI BLAST (Blast.cgi Put/Get). "
            "HelixScope does not compute BLAST scores locally."
        ),
        "local_blast_plus": local,
        "programs": tuple(BLAST_PROGRAMS.keys()),
        "recommended_tools": [
            "NCBI BLAST web form",
            "BLAST+ against a local database",
        ],
    }


def programs_for_molecule(molecule: str) -> List[str]:
    """Lista programas BLAST compativeis com o tipo molecular da query.

    Args:
        molecule: DNA, RNA ou PROTEIN.

    Returns:
        Nomes de programa em ordem estavel.

    Raises:
        ValueError: Se molecule nao for reconhecido.
    """
    key = str(molecule or "").strip().upper()
    if key not in {"DNA", "RNA", "PROTEIN"}:
        raise ValueError("molecule must be DNA, RNA or PROTEIN.")
    return [
        name
        for name, spec in BLAST_PROGRAMS.items()
        if key in spec["query_molecules"]
    ]


def databases_for_program(program: str) -> Tuple[str, ...]:
    """Bancos permitidos para um programa.

    Args:
        program: blastn, blastp, blastx, tblastn ou tblastx.

    Returns:
        Tupla de nomes de banco NCBI.

    Raises:
        BlastError: INVALID_INPUT se o programa nao for suportado.
    """
    spec = BLAST_PROGRAMS.get(str(program or "").strip().lower())
    if spec is None:
        raise BlastError(
            f"BLAST program '{program}' is not offered. Supported: "
            + ", ".join(BLAST_PROGRAMS),
            "INVALID_INPUT",
        )
    return spec["databases"]


def prepare_query(sequence: str, program: str) -> Dict[str, Any]:
    """Valida a query contra o programa e devolve o texto enviado ao NCBI.

    Args:
        sequence: Sequencia do usuario (FASTA ou cru).
        program: Programa BLAST.

    Returns:
        Dict com sequence (normalizada), molecule, fasta, query_hash,
        rna_to_dna (bool) e length.

    Raises:
        BlastError: INVALID_INPUT, RESOURCE_LIMIT.
    """
    program_key = str(program or "").strip().lower()
    if program_key not in BLAST_PROGRAMS:
        raise BlastError(
            f"BLAST program '{program}' is not offered.",
            "INVALID_INPUT",
        )
    parsed = dna_analysis.parse_sequence_payload(sequence or "")
    raw = str(parsed.get("sequence") or "")
    info = dna_analysis.validate_sequence(raw)
    if not info["is_valid"]:
        raise BlastError(
            info.get("rejection_reason") or "Query sequence is invalid.",
            "INVALID_INPUT",
        )
    molecule = str(info["type"])
    allowed = BLAST_PROGRAMS[program_key]["query_molecules"]
    if molecule not in allowed:
        raise BlastError(
            f"Cannot run {program_key} on a {molecule} query. "
            f"Use one of: {', '.join(allowed)}.",
            "INVALID_INPUT",
        )
    residues = str(info["sequence"])
    rna_to_dna = False
    if molecule == "RNA" and program_key in {"blastn", "blastx", "tblastx"}:
        residues = residues.replace("U", "T")
        rna_to_dna = True
    if len(residues) > MAX_QUERY_LENGTH:
        raise BlastError(
            f"BLAST query is limited to {MAX_QUERY_LENGTH:,} residues. "
            "Paste a shorter region.",
            "RESOURCE_LIMIT",
        )
    if not residues:
        raise BlastError("Query contains no sequence data.", "INVALID_INPUT")
    fasta = f">helixscope_query\n{residues}"
    return {
        "sequence": residues,
        "molecule": molecule,
        "fasta": fasta,
        "query_hash": provenance.sequence_digest(residues),
        "rna_to_dna": rna_to_dna,
        "length": len(residues),
        "identifier": str(parsed.get("identifier") or ""),
    }


def cache_key(
    *,
    query_hash: str,
    program: str,
    database: str,
    expect: float,
    hitlist_size: int,
    backend: str = "ncbi_remote",
    database_identity: str = "",
    tool_version: str = "",
) -> str:
    """Chave de cache: hash da query + programa + banco + parametros.

    Args:
        query_hash: SHA-256 da sequencia enviada.
        program: Programa BLAST.
        database: Banco NCBI ou identidade do prefixo local.
        expect: limiar E-value.
        hitlist_size: HITLIST_SIZE efetivo.
        backend: ncbi_remote ou blastplus_local.
        database_identity: Identidade do banco local (vazio no remoto).
        tool_version: Versao BLAST+ local; vazio no remoto NCBI.

    Returns:
        Digest hexadecimal da chave composta.

    Raises:
        Nenhum.
    """
    import hashlib

    payload = (
        f"{query_hash}|{program}|{database}|{expect}|{hitlist_size}|"
        f"{backend}|{database_identity}|{tool_version}|{provenance.HELIXSCOPE_VERSION}".encode("utf-8")
    )
    return hashlib.sha256(payload).hexdigest()


def seconds_until_next_poll(job: Mapping[str, Any], now: Optional[float] = None) -> float:
    """Segundos restantes antes do proximo Get, segundo a politica NCBI.

    Args:
        job: Dict de job com submitted_monotonic e last_poll_monotonic.
        now: Relogio monotonic; None usa time.monotonic.

    Returns:
        0.0 se o poll for permitido agora.

    Raises:
        Nenhum.
    """
    import time as time_module

    current = time_module.monotonic() if now is None else float(now)
    last_poll = job.get("last_poll_monotonic")
    submitted_raw = job.get("submitted_monotonic")
    submitted = current if submitted_raw is None else float(submitted_raw)
    rtoe = float(job.get("rtoe_s") or MIN_CONTACT_INTERVAL_S)
    first_wait = max(MIN_CONTACT_INTERVAL_S, min(rtoe, 60.0))
    if last_poll is None:
        remaining = (submitted + first_wait) - current
        return max(0.0, remaining)
    remaining = (float(last_poll) + RID_POLL_INTERVAL_S) - current
    return max(0.0, remaining)


def submit_search(
    sequence: str,
    email: str,
    *,
    program: str = "blastn",
    database: str = "",
    expect: float = 10.0,
    hitlist_size: int = DEFAULT_HITLIST_SIZE,
    urlopen_fn=None,
) -> Dict[str, Any]:
    """Submete CMD=Put ao NCBI BLAST e devolve o RID.

    Args:
        sequence: Query (FASTA ou cru).
        email: E-mail de contato exigido pelo NCBI.
        program: Um de BLAST_PROGRAMS.
        database: Banco da allowlist do programa; vazio usa o padrao.
        expect: E-value maximo pedido ao NCBI (nao calculado aqui).
        hitlist_size: Numero de sujeitos a manter.
        urlopen_fn: urlopen injetavel para testes; None usa a rede.

    Returns:
        Job com rid, rtoe_s, status WAITING, parametros e proveniencia.

    Raises:
        BlastError: VALIDACAO ou falha HTTP/NCBI.
        ValueError: E-mail invalido.
    """
    import time as time_module

    _require_email(email)
    prepared = prepare_query(sequence, program)
    program_key = str(program).strip().lower()
    db = str(database or DEFAULT_DATABASE[program_key]).strip()
    allowed_db = databases_for_program(program_key)
    if db not in allowed_db:
        raise BlastError(
            f"Database '{db}' is not enabled for {program_key}. "
            f"Use one of: {', '.join(allowed_db)}.",
            "INVALID_DATABASE",
        )
    try:
        expect_value = float(expect)
    except (TypeError, ValueError) as exc:
        raise BlastError("expect must be a positive number.", "INVALID_INPUT") from exc
    if not math.isfinite(expect_value) or expect_value <= 0:
        raise BlastError("expect must be a positive number.", "INVALID_INPUT")
    try:
        hits = int(hitlist_size)
    except (TypeError, ValueError) as exc:
        raise BlastError("hitlist_size must be an integer.", "INVALID_INPUT") from exc
    hits = max(1, min(hits, MAX_HITLIST_SIZE))

    parameters = {
        "CMD": "Put",
        "PROGRAM": program_key,
        "DATABASE": db,
        "QUERY": prepared["fasta"],
        "EXPECT": str(expect_value),
        "HITLIST_SIZE": str(hits),
        "FORMAT_TYPE": "XML",
        "TOOL": BLAST_TOOL,
        "EMAIL": email.strip(),
    }
    body = _blast_request(parameters, urlopen_fn=urlopen_fn)
    rid, rtoe = parse_put_response(body)
    submitted = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {
        "rid": rid,
        "rtoe_s": rtoe,
        "status": "WAITING",
        "program": program_key,
        "database": db,
        "expect": expect_value,
        "hitlist_size": hits,
        "query_hash": prepared["query_hash"],
        "query_length": prepared["length"],
        "molecule": prepared["molecule"],
        "rna_to_dna": prepared["rna_to_dna"],
        "submitted_at_utc": submitted,
        "submitted_monotonic": time_module.monotonic(),
        "last_poll_monotonic": None,
        "source": "NCBI BLAST Common URL API",
        "endpoint": BLAST_ENDPOINT,
        "cache_status": "live",
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def poll_search(
    job: Mapping[str, Any],
    email: str,
    *,
    urlopen_fn=None,
    now: Optional[float] = None,
) -> Dict[str, Any]:
    """Consulta CMD=Get o status do RID, respeitando o intervalo NCBI.

    Args:
        job: Job devolvido por submit_search.
        email: E-mail de contato.
        urlopen_fn: urlopen injetavel.
        now: Relogio monotonic opcional.

    Returns:
        Copia do job com status WAITING, READY ou FAILED.

    Raises:
        BlastError: TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE, JOB_FAILED,
            INVALID_INPUT.
    """
    import time as time_module

    _require_email(email)
    current = time_module.monotonic() if now is None else float(now)
    updated = dict(job)
    rid = str(updated.get("rid") or "")
    if not RID_PATTERN.fullmatch(rid):
        raise BlastError("BLAST job RID is invalid.", "INVALID_INPUT")
    submitted_raw = updated.get("submitted_monotonic")
    submitted_mono = current if submitted_raw is None else float(submitted_raw)
    if current - submitted_mono > MAX_JOB_WAIT_S:
        raise BlastError(
            f"BLAST job {rid} exceeded {int(MAX_JOB_WAIT_S)}s without READY. "
            "This is TIMEOUT, not NO_HITS.",
            "TIMEOUT",
        )
    wait = seconds_until_next_poll(updated, now=current)
    if wait > 0:
        raise BlastError(
            f"Wait {wait:.0f}s before polling RID {rid} (NCBI poll policy).",
            "RATE_LIMITED",
        )
    parameters = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo",
        "TOOL": BLAST_TOOL,
        "EMAIL": email.strip(),
    }
    body = _blast_request(parameters, urlopen_fn=urlopen_fn)
    status = parse_status_response(body)
    updated["last_poll_monotonic"] = current
    updated["status"] = status
    updated["status_checked_at_utc"] = datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    if status == "FAILED":
        raise BlastError(
            f"NCBI reported BLAST job {rid} as FAILED. This is not NO_HITS.",
            "JOB_FAILED",
        )
    if status == "UNKNOWN":
        raise BlastError(
            f"NCBI status for RID {rid} is UNKNOWN. The job may have expired.",
            "JOB_FAILED",
        )
    return updated


def retrieve_search(
    job: Mapping[str, Any],
    email: str,
    *,
    urlopen_fn=None,
) -> Dict[str, Any]:
    """Baixa FORMAT_TYPE=XML quando o job esta READY e interpreta BlastOutput.

    Args:
        job: Job com rid e parametros.
        email: E-mail de contato.
        urlopen_fn: urlopen injetavel.

    Returns:
        Resultado validado com hits, proveniencia e status READY ou NO_HITS.

    Raises:
        BlastError: PARSING_ERROR, RESOURCE_LIMIT, falhas HTTP, job nao READY.
    """
    _require_email(email)
    rid = str(job.get("rid") or "")
    if not RID_PATTERN.fullmatch(rid):
        raise BlastError("BLAST job RID is invalid.", "INVALID_INPUT")
    if str(job.get("status") or "") != "READY":
        raise BlastError(
            "BLAST results can be retrieved only when status is READY.",
            "JOB_FAILED",
        )
    parameters = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": str(int(job.get("hitlist_size") or DEFAULT_HITLIST_SIZE)),
        "ALIGNMENTS": str(int(job.get("hitlist_size") or DEFAULT_HITLIST_SIZE)),
        "DESCRIPTIONS": str(int(job.get("hitlist_size") or DEFAULT_HITLIST_SIZE)),
        "TOOL": BLAST_TOOL,
        "EMAIL": email.strip(),
    }
    body = _blast_request(parameters, urlopen_fn=urlopen_fn)
    parsed = parse_blast_xml(body)
    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    hits = list(parsed["hits"])
    skipped = int(parsed.get("skipped_hsps") or 0)
    if skipped and not hits:
        raise BlastError(
            "BLAST XML was parsed but every HSP failed validation. "
            "This is PARSING_ERROR, not NO_HITS.",
            "PARSING_ERROR",
        )
    status = "NO_HITS" if not hits else "READY"
    envelope = provenance.analysis_envelope(
        module="BLAST",
        payload={
            "rid": rid,
            "program": parsed.get("program") or job.get("program"),
            "blast_version": parsed.get("blast_version") or "",
            "database": parsed.get("database") or job.get("database"),
            "query_id": parsed.get("query_id") or "",
            "query_def": parsed.get("query_def") or "",
            "query_length": parsed.get("query_length") or job.get("query_length"),
            "hits": hits,
            "n_hits": len(hits),
            "skipped_hsps": parsed.get("skipped_hsps") or 0,
            "parse_warnings": parsed.get("warnings") or [],
            "expect": job.get("expect"),
            "hitlist_size": job.get("hitlist_size"),
            "status": status,
            "cache_status": "live",
            "endpoint": BLAST_ENDPOINT,
        },
        status="RETRIEVED",
        algorithm="NCBI BLAST Common URL API",
        parameters={
            "program": job.get("program"),
            "database": job.get("database"),
            "expect": job.get("expect"),
            "hitlist_size": job.get("hitlist_size"),
            "format_type": "XML",
        },
        sequence="",
        source="NCBI BLAST",
        accession="",
        input_identifier=str(job.get("query_hash") or ""),
    )
    envelope["query_hash"] = job.get("query_hash")
    envelope["rid"] = rid
    envelope["retrieved_at_utc"] = retrieved
    envelope["submitted_at_utc"] = job.get("submitted_at_utc") or ""
    envelope["cache_status"] = "live"
    envelope["hits"] = hits
    envelope["status"] = status
    envelope["n_hits"] = len(hits)
    envelope["blast_version"] = parsed.get("blast_version") or ""
    envelope["program"] = parsed.get("program") or job.get("program")
    envelope["database"] = parsed.get("database") or job.get("database")
    envelope["query_length"] = parsed.get("query_length") or job.get("query_length")
    envelope["parse_warnings"] = parsed.get("warnings") or []
    return envelope


def parse_put_response(body: str) -> Tuple[str, float]:
    """Extrai RID e RTOE da pagina de confirmacao do NCBI.

    Args:
        body: Texto HTML/plain devolvido por CMD=Put.

    Returns:
        Tupla (rid, rtoe_segundos). rtoe ausente vira MIN_CONTACT_INTERVAL_S.

    Raises:
        BlastError: JOB_FAILED ou INVALID_INPUT se o NCBI devolver erro.
    """
    text = body or ""
    _raise_if_ncbi_error_page(text)
    rid_match = re.search(r"RID\s*=\s*([A-Za-z0-9]+)", text)
    rtoe_match = re.search(r"RTOE\s*=\s*([0-9]+)", text)
    if not rid_match:
        raise BlastError(
            "NCBI Put response did not include a RID. The search was not queued.",
            "JOB_FAILED",
        )
    rid = rid_match.group(1)
    if not RID_PATTERN.fullmatch(rid):
        raise BlastError("NCBI returned a RID that failed validation.", "PARSING_ERROR")
    rtoe = float(rtoe_match.group(1)) if rtoe_match else MIN_CONTACT_INTERVAL_S
    return rid, rtoe


def parse_status_response(body: str) -> str:
    """Le Status= da pagina SearchInfo.

    Args:
        body: Texto devolvido por CMD=Get FORMAT_OBJECT=SearchInfo.

    Returns:
        WAITING, READY, FAILED ou UNKNOWN (maiúsculas).

    Raises:
        BlastError: se a pagina for um erro NCBI explicito.
    """
    text = body or ""
    _raise_if_ncbi_error_page(text)
    if "<BlastOutput" in text:
        return "READY"
    match = re.search(r"Status\s*=\s*([A-Za-z]+)", text)
    if not match:
        if "Status=" not in text and text.strip() and "<html" not in text.lower():
            return "READY"
        return "WAITING"
    status = match.group(1).strip().upper()
    if status in {"WAITING", "READY", "FAILED", "UNKNOWN"}:
        return status
    return "UNKNOWN"


def parse_blast_xml(body: str) -> Dict[str, Any]:
    """Interpreta BlastOutput XML do NCBI e valida cada HSP.

    Args:
        body: Documento XML (FORMAT_TYPE=XML).

    Returns:
        Dict com program, blast_version, database, query_* , hits e warnings.
        Hits vazios sao NO_HITS, nao erro.

    Raises:
        BlastError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    if body is None:
        raise BlastError("BLAST XML body is empty.", "PARSING_ERROR")
    raw = body if isinstance(body, str) else body.decode("utf-8", errors="replace")
    encoded = raw.encode("utf-8")
    if len(encoded) > MAX_XML_BYTES:
        raise BlastError(
            f"BLAST XML exceeds {MAX_XML_BYTES:,} bytes. Resource limit reached.",
            "RESOURCE_LIMIT",
        )
    stripped = raw.strip()
    if not stripped:
        raise BlastError("BLAST XML body is empty.", "PARSING_ERROR")
    _raise_if_ncbi_error_page(stripped)
    if "<BlastOutput" not in stripped:
        status = parse_status_response(stripped)
        if status == "WAITING":
            raise BlastError(
                "BLAST XML was requested but NCBI still reports WAITING.",
                "JOB_FAILED",
            )
        raise BlastError(
            "Response is not NCBI BlastOutput XML.",
            "PARSING_ERROR",
        )
    try:
        root = ET.fromstring(stripped)
    except ET.ParseError as exc:
        raise BlastError(f"BLAST XML could not be parsed: {exc}.", "PARSING_ERROR") from exc
    if _local_tag(root.tag) != "BlastOutput":
        raise BlastError("Root element is not BlastOutput.", "PARSING_ERROR")

    program = _child_text(root, "BlastOutput_program")
    version = _child_text(root, "BlastOutput_version")
    database = _child_text(root, "BlastOutput_db")
    query_id = _child_text(root, "BlastOutput_query-ID")
    query_def = _child_text(root, "BlastOutput_query-def")
    query_len = _child_int(root, "BlastOutput_query-len")
    warnings: List[str] = []
    skipped = 0
    hits: List[dict] = []

    iterations = _find_child(root, "BlastOutput_iterations")
    if iterations is None:
        raise BlastError(
            "BlastOutput XML has no BlastOutput_iterations element.",
            "PARSING_ERROR",
        )
    for iteration in _findall_local(iterations, "Iteration"):
        message = _child_text(iteration, "Iteration_message").lower()
        if "no hits found" in message:
            continue
        hit_block = _find_child(iteration, "Iteration_hits")
        if hit_block is None:
            continue
        for hit_el in _findall_local(hit_block, "Hit"):
            parsed_hit, skipped_here, hit_warnings = _parse_hit(hit_el, query_len)
            skipped += skipped_here
            warnings.extend(hit_warnings)
            if parsed_hit is not None:
                hits.append(parsed_hit)
    return {
        "program": program,
        "blast_version": version,
        "database": database,
        "query_id": query_id,
        "query_def": query_def,
        "query_length": query_len,
        "hits": hits,
        "warnings": warnings,
        "skipped_hsps": skipped,
    }


def paginate_hits(hits: List[dict], page: int, page_size: int = MAX_DISPLAY_HITS) -> dict:
    """Recorta a lista de hits sem carregar milhares de linhas na UI.

    Args:
        hits: Lista ja validada.
        page: Pagina 1-based.
        page_size: Tamanho da pagina; limitado a MAX_DISPLAY_HITS.

    Returns:
        Dict com items, page, page_size, n_pages, n_hits.

    Raises:
        BlastError: INVALID_INPUT se page_size for invalido.
    """
    size = max(1, min(int(page_size), MAX_DISPLAY_HITS))
    total = len(hits)
    n_pages = max(1, math.ceil(total / size) if total else 1)
    current = max(1, min(int(page), n_pages))
    start = (current - 1) * size
    return {
        "items": hits[start : start + size],
        "page": current,
        "page_size": size,
        "n_pages": n_pages,
        "n_hits": total,
    }


def _parse_hit(hit_el: ET.Element, query_len: int) -> Tuple[Optional[dict], int, List[str]]:
    """Converte um elemento Hit em dict; descarta HSPs inconsistentes."""
    warnings: List[str] = []
    skipped = 0
    accession = _child_text(hit_el, "Hit_accession")
    hit_id = _child_text(hit_el, "Hit_id")
    description = _child_text(hit_el, "Hit_def")
    hit_len = _child_int(hit_el, "Hit_len")
    hsps: List[dict] = []
    hsps_el = _find_child(hit_el, "Hit_hsps")
    if hsps_el is None:
        return None, 0, ["Hit without HSPs was omitted."]
    for hsp_el in _findall_local(hsps_el, "Hsp"):
        hsp, reason = _parse_hsp(hsp_el, query_len, hit_len)
        if hsp is None:
            skipped += 1
            warnings.append(reason or "Inconsistent HSP omitted.")
            continue
        hsps.append(hsp)
    if not hsps:
        return None, skipped, warnings
    best = hsps[0]
    organism = _organism_from_def(description)
    return (
        {
            "hit_id": hit_id,
            "accession": accession or hit_id,
            "description": description,
            "organism": organism,
            "hit_length": hit_len if hit_len > 0 else None,
            "bit_score": best["bit_score"],
            "score": best["score"],
            "evalue": best["evalue"],
            "identities": best["identities"],
            "positives": best["positives"],
            "gaps": best["gaps"],
            "alignment_length": best["alignment_length"],
            "identity_pct": best["identity_pct"],
            "query_coverage_pct": best["query_coverage_pct"],
            "query_from": best["query_from"],
            "query_to": best["query_to"],
            "hit_from": best["hit_from"],
            "hit_to": best["hit_to"],
            "qseq": best["qseq"],
            "hseq": best["hseq"],
            "midline": best["midline"],
            "hsps": hsps,
        },
        skipped,
        warnings,
    )


def _parse_hsp(
    hsp_el: ET.Element, query_len: int, hit_len: int
) -> Tuple[Optional[dict], str]:
    """Valida um HSP NCBI. Nao corrige numeros invalidos."""
    align_len = _child_int(hsp_el, "Hsp_align-len")
    identities = _child_int(hsp_el, "Hsp_identity")
    positives = _child_int(hsp_el, "Hsp_positive")
    gaps = _child_int(hsp_el, "Hsp_gaps")
    qseq = _child_text(hsp_el, "Hsp_qseq")
    hseq = _child_text(hsp_el, "Hsp_hseq")
    midline = _child_text(hsp_el, "Hsp_midline")
    query_from = _child_int(hsp_el, "Hsp_query-from")
    query_to = _child_int(hsp_el, "Hsp_query-to")
    hit_from = _child_int(hsp_el, "Hsp_hit-from")
    hit_to = _child_int(hsp_el, "Hsp_hit-to")
    if align_len <= 0:
        return None, "HSP omitted: alignment_length missing or not positive."
    if identities < 0 or identities > align_len:
        return None, "HSP omitted: identities exceed alignment_length."
    if positives < 0 or (positives and positives > align_len):
        return None, "HSP omitted: positives exceed alignment_length."
    if gaps < 0 or gaps > align_len:
        return None, "HSP omitted: gaps exceed alignment_length."
    if qseq and hseq and len(qseq) != len(hseq):
        return None, "HSP omitted: query and subject aligned strings differ in length."
    if qseq and align_len and len(qseq) != align_len:
        return None, "HSP omitted: aligned query length disagrees with alignment_length."
    if midline and qseq and len(midline) != len(qseq):
        return None, "HSP omitted: midline length disagrees with aligned query."
    if query_from <= 0 or query_to <= 0:
        return None, "HSP omitted: query coordinates are not valid 1-based positions."
    if hit_from <= 0 or hit_to <= 0:
        return None, "HSP omitted: subject coordinates are not valid 1-based positions."
    if query_len > 0 and max(query_from, query_to) > query_len:
        return None, "HSP omitted: query coordinates fall outside query length."
    if hit_len > 0 and max(hit_from, hit_to) > hit_len:
        return None, "HSP omitted: subject coordinates fall outside subject length."
    identity_pct = round((identities / align_len) * 100.0, 2)
    query_span = abs(query_to - query_from) + 1
    coverage = (
        round((query_span / query_len) * 100.0, 2) if query_len > 0 else float("nan")
    )
    if not (0.0 <= identity_pct <= 100.0):
        return None, "HSP omitted: identity percent is outside 0-100."
    if isinstance(coverage, float) and not math.isnan(coverage) and not (
        0.0 <= coverage <= 100.0
    ):
        return None, "HSP omitted: query coverage percent is outside 0-100."
    bit_score = _child_float(hsp_el, "Hsp_bit-score")
    score = _child_float(hsp_el, "Hsp_score")
    evalue = _child_float(hsp_el, "Hsp_evalue")
    if bit_score is None or evalue is None:
        return None, "HSP omitted: bit score or E-value missing from NCBI XML."
    parsed = {
        "bit_score": bit_score,
        "score": score,
        "evalue": evalue,
        "identities": identities,
        "positives": positives if positives else identities,
        "gaps": gaps,
        "alignment_length": align_len,
        "identity_pct": identity_pct,
        "query_coverage_pct": coverage,
        "identity_definition": "NCBI Hsp_identity / Hsp_align-len, percent",
        "query_coverage_definition": (
            "(|query_to - query_from| + 1) / query_length, percent, NCBI 1-based HSP coordinates"
        ),
        "query_from": query_from,
        "query_to": query_to,
        "hit_from": hit_from,
        "hit_to": hit_to,
        "qseq": qseq,
        "hseq": hseq,
        "midline": midline,
        "query_frame": _child_text(hsp_el, "Hsp_query-frame"),
        "hit_frame": _child_text(hsp_el, "Hsp_hit-frame"),
    }
    if not scientific_checks.blast_hsp_is_valid(
        parsed, query_len=query_len, hit_len=hit_len
    ):
        return None, "HSP omitted: failed HelixScope BLAST invariants."
    return parsed, ""


def _organism_from_def(definition: str) -> str:
    """Tenta ler o organismo entre colchetes no Hit_def NCBI, se existir."""
    match = re.search(r"\[([^\[\]]+)\]\s*$", definition or "")
    if not match:
        return ""
    return match.group(1).strip()


def _require_email(email: str) -> None:
    from . import ncbi_fetch

    if not ncbi_fetch._valid_entrez_email(email or ""):
        raise ValueError(
            "Informe um e-mail de contato valido para o NCBI BLAST "
            "(formato local@dominio)."
        )


def _blast_request(parameters: Mapping[str, str], *, urlopen_fn=None) -> str:
    """POST application/x-www-form-urlencoded somente ao BLAST_ENDPOINT."""
    opener = urlopen_fn or urlopen
    payload = urlencode(dict(parameters)).encode("utf-8")
    request = Request(
        BLAST_ENDPOINT,
        data=payload,
        method="POST",
        headers={
            "User-Agent": f"{BLAST_TOOL}/homology",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    try:
        with opener(request, timeout=BLAST_HTTP_TIMEOUT_S) as handle:
            raw = handle.read()
    except HTTPError as exc:
        raise _http_error(exc) from exc
    except (URLError, TimeoutError, socket.timeout, OSError) as exc:
        message = str(exc).lower()
        if "timed out" in message or isinstance(exc, (TimeoutError, socket.timeout)):
            raise BlastError(
                "NCBI BLAST request timed out. This is TIMEOUT, not NO_HITS.",
                "TIMEOUT",
            ) from exc
        raise BlastError(
            f"NCBI BLAST unavailable: {exc}. This is not NO_HITS.",
            "SERVICE_UNAVAILABLE",
        ) from exc
    if isinstance(raw, bytes):
        if len(raw) > MAX_XML_BYTES:
            raise BlastError(
                f"BLAST response exceeds {MAX_XML_BYTES:,} bytes.",
                "RESOURCE_LIMIT",
            )
        return raw.decode("utf-8", errors="replace")
    text = str(raw)
    if len(text.encode("utf-8")) > MAX_XML_BYTES:
        raise BlastError(
            f"BLAST response exceeds {MAX_XML_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    return text


def _http_error(exc: HTTPError) -> BlastError:
    code = int(getattr(exc, "code", 0) or 0)
    if code == 429:
        return BlastError(
            "NCBI BLAST rate limited the request (HTTP 429). This is not NO_HITS.",
            "RATE_LIMITED",
        )
    if code >= 500:
        return BlastError(
            f"NCBI BLAST unavailable (HTTP {code}). This is not NO_HITS.",
            "SERVICE_UNAVAILABLE",
        )
    if code in {400, 404}:
        return BlastError(
            f"NCBI BLAST rejected the request (HTTP {code}).",
            "INVALID_INPUT",
        )
    return BlastError(
        f"NCBI BLAST HTTP {code}: {exc}.",
        "SERVICE_UNAVAILABLE",
    )


def _raise_if_ncbi_error_page(text: str) -> None:
    lowered = text.lower()
    if "message id#" in lowered or 'class="error' in lowered:
        match = re.search(r"Message ID#\d+\s*Error:\s*([^\n<]+)", text, re.IGNORECASE)
        detail = match.group(1).strip() if match else "NCBI BLAST returned an error page."
        category = "INVALID_INPUT"
        if "no data" in detail.lower() or "no sequence" in detail.lower():
            category = "INVALID_INPUT"
        raise BlastError(f"NCBI BLAST: {detail}", category)
    if "error: failed to read the blast query" in lowered:
        raise BlastError(
            "NCBI BLAST failed to read the query (molecule/program mismatch or empty query).",
            "INVALID_INPUT",
        )


def _local_tag(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", 1)[-1]
    return tag


def _find_child(parent: ET.Element, name: str) -> Optional[ET.Element]:
    direct = parent.find(name)
    if direct is not None:
        return direct
    for child in parent:
        if _local_tag(child.tag) == name:
            return child
    return None


def _child_text(parent: ET.Element, name: str) -> str:
    node = _find_child(parent, name)
    if node is None or node.text is None:
        return ""
    return str(node.text).strip()


def _child_int(parent: ET.Element, name: str) -> int:
    text = _child_text(parent, name)
    if not text:
        return 0
    try:
        return int(float(text))
    except ValueError:
        return 0


def _child_float(parent: ET.Element, name: str) -> Optional[float]:
    text = _child_text(parent, name)
    if not text:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def mark_cached_result(result: Mapping[str, Any]) -> Dict[str, Any]:
    """Marca um resultado servido do cache sem fingir nova recuperacao.

    Args:
        result: Envelope devolvido por retrieve_search.

    Returns:
        Copia com cache_status=cached e served_from_cache_at_utc. Mantem
        retrieved_at_utc original.

    Raises:
        Nenhum.
    """
    copied = dict(result)
    copied["cache_status"] = "cached"
    copied["served_from_cache_at_utc"] = datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    return copied


def export_hits_table(result: Mapping[str, Any]) -> str:
    """CSV dos hits NCBI. Sem hits, so cabecalho e linhas de proveniencia.

    Args:
        result: Envelope BLAST.

    Returns:
        Texto CSV UTF-8.

    Raises:
        Nenhum.
    """
    import csv
    import io

    buffer = io.StringIO()
    fieldnames = [
        "accession",
        "hit_id",
        "description",
        "organism",
        "bit_score",
        "score",
        "evalue",
        "identities",
        "positives",
        "gaps",
        "alignment_length",
        "identity_pct",
        "query_coverage_pct",
        "query_from",
        "query_to",
        "hit_from",
        "hit_to",
    ]
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for hit in list(result.get("hits") or []):
        if not isinstance(hit, dict):
            continue
        writer.writerow({key: hit.get(key, "") for key in fieldnames})
    return buffer.getvalue()


def export_result_bundle(result: Mapping[str, Any]) -> dict:
    """Copia JSON-serializavel do resultado, com proveniencia.

    Args:
        result: Envelope BLAST.

    Returns:
        Dict sem NaN/Inf.

    Raises:
        Nenhum.
    """
    return provenance.json_safe(dict(result))


def deployment_readiness() -> dict:
    """O que esta preparado vs o que ainda nao existe para um deploy publico.

    Args:
        Nenhum.

    Returns:
        Flags honestas. authentication e quotas continuam False.

    Raises:
        Nenhum.
    """
    return {
        "authentication": False,
        "per_user_quotas": False,
        "job_quotas": False,
        "ncbi_rate_isolation": False,
        "distributed_workers": False,
        "session_job_limit": MAX_CONCURRENT_JOBS,
        "session_result_cache_limit": MAX_SESSION_CACHE,
        "endpoint_allowlist": [BLAST_ENDPOINT],
        "local_blast_plus": bool(local_blast_availability().get("available")),
    }


def _findall_local(parent: ET.Element, name: str) -> List[ET.Element]:
    """Filhos com o nome local, com ou sem namespace XML."""
    found = list(parent.findall(name))
    if found:
        return found
    return [child for child in parent if _local_tag(child.tag) == name]

