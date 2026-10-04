"""Fila local de jobs CRISPR em escala genomica.

Worker em subprocesso separado (nao congela a UI). CPU limitada a 1 thread
por omissao. Job IDs sao identificadores de colisao, nao segredos. Auth,
quotas e workers publicos permanecem UNAVAILABLE.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import json
import os
import secrets
import signal
import subprocess
import sys
from typing import Any, Mapping, Optional

from . import genome_store, provenance

JOB_QUEUED: str = "QUEUED"
JOB_RUNNING: str = "RUNNING"
JOB_COMPLETED_FULL_REFERENCE: str = "COMPLETED_FULL_REFERENCE"
JOB_COMPLETED_TEST_REFERENCE: str = "COMPLETED_TEST_REFERENCE"
JOB_CANCELLED: str = "CANCELLED"
JOB_TIMEOUT: str = "TIMEOUT"
JOB_RESOURCE_LIMIT: str = "RESOURCE_LIMIT"
JOB_ERROR: str = "ERROR"
JOB_INTERRUPTED: str = "INTERRUPTED"

STAGE_PREPARING: str = "Preparing"
STAGE_RUNNING: str = "Running"
STAGE_VALIDATING: str = "Validating"
STAGE_SCORING: str = "Scoring"
STAGE_COMPLETED: str = "Completed"

DEFAULT_TIMEOUT_S: float = 3600.0
DEFAULT_MAX_HITS: int = 5_000
DEFAULT_THREADS: int = 1
MAX_THREADS: int = 4


class GenomeJobError(ValueError):
    """Falha da fila de jobs.

    Attributes:
        category: INVALID_INPUT, UNAVAILABLE, ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def new_job_id() -> str:
    """Identificador de 32 hex. Nao e um segredo de sessao.

    Args:
        Nenhum.

    Returns:
        str.

    Raises:
        Nenhum.
    """
    return secrets.token_hex(16)


def job_dir(job_id: str) -> str:
    """Pasta de um job.

    Args:
        job_id: id.

    Returns:
        Caminho absoluto.

    Raises:
        GenomeJobError: INVALID_INPUT.
    """
    ident = str(job_id or "").strip()
    if len(ident) != 32 or any(ch not in "0123456789abcdef" for ch in ident):
        raise GenomeJobError("Job id is not a HelixScope hex identifier.", "INVALID_INPUT")
    return os.path.join(genome_store.jobs_root(), ident)


def _atomic_write(path: str, payload: Mapping[str, Any]) -> None:
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(dict(payload), handle, indent=2)
    os.replace(tmp, path)


def read_json(path: str) -> Optional[dict]:
    """Le JSON ou None.

    Args:
        path: Ficheiro.

    Returns:
        dict ou None.

    Raises:
        Nenhum.
    """
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def write_status(job_id: str, status: Mapping[str, Any]) -> dict:
    """Grava status.json.

    Args:
        job_id: id.
        status: Dict.

    Returns:
        Dict gravado.

    Raises:
        GenomeJobError: INVALID_INPUT.
    """
    payload = dict(status)
    payload["job_id"] = str(job_id)
    payload["updated_at"] = provenance.utc_now()
    _atomic_write(os.path.join(job_dir(job_id), "status.json"), payload)
    return payload


def load_job(job_id: str) -> Optional[dict]:
    """Le job.json.

    Args:
        job_id: id.

    Returns:
        dict ou None.

    Raises:
        Nenhum.
    """
    try:
        return read_json(os.path.join(job_dir(job_id), "job.json"))
    except GenomeJobError:
        return None


def load_status(job_id: str) -> Optional[dict]:
    """Le status.json e promove RUNNING morto a ERROR.

    Args:
        job_id: id.

    Returns:
        dict ou None.

    Raises:
        Nenhum.
    """
    try:
        status = read_json(os.path.join(job_dir(job_id), "status.json"))
    except GenomeJobError:
        return None
    if not status:
        return None
    if str(status.get("status") or "") == JOB_RUNNING:
        pid = status.get("pid")
        if pid and not _pid_alive(int(pid)):
            status = write_status(
                job_id,
                {
                    **status,
                    "status": JOB_INTERRUPTED,
                    "reason": (
                        "INTERRUPTED: worker process is no longer running and "
                        "did not write a completed result. Not converted to "
                        "COMPLETED_FULL_REFERENCE or zero hits."
                    ),
                },
            )
    return status


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        if os.name == "nt":
            import ctypes

            handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, int(pid))
            if handle:
                ctypes.windll.kernel32.CloseHandle(handle)
                return True
            return False
        os.kill(int(pid), 0)
        return True
    except OSError:
        return False


def cache_identity(
    *,
    guide_hash: str,
    assembly_id: str,
    reference_sha256: str,
    engine: str,
    engine_version: str,
    max_mismatches: int,
    dna_bulge: int,
    rna_bulge: int,
    include_nag: bool,
) -> str:
    """Identidade de cache. Versao do engine e hash da referencia entram.

    Args:
        guide_hash: Hash do spacer.
        assembly_id: id.
        reference_sha256: Hash do FASTA READY.
        engine: cas_offinder.
        engine_version: versao do binario.
        max_mismatches: teto.
        dna_bulge: bulge DNA.
        rna_bulge: bulge RNA.
        include_nag: se NAG foi pesquisado.

    Returns:
        SHA-256 hex.

    Raises:
        Nenhum.
    """
    import hashlib

    payload = "|".join(
        [
            str(guide_hash),
            str(assembly_id),
            str(reference_sha256),
            str(engine),
            str(engine_version),
            str(int(max_mismatches)),
            str(int(dna_bulge)),
            str(int(rna_bulge)),
            "nag" if include_nag else "ngg_only",
            provenance.HELIXSCOPE_VERSION,
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def full_genome_search_preflight(assembly_id: str) -> dict:
    """FULL GENOME so com reference READY e Cas-OFFinder LIVE_VALIDATED.

    TEST REFERENCE nao e FULL GENOME; o preflight public recusa o FASTA de teste.

    Args:
        assembly_id: id de catalogo.

    Returns:
        Dict allowed, reason, engine_status, reference_status.

    Raises:
        Nenhum.
    """
    from . import engine_validation, crispr_casoffinder

    ready = genome_store.ready_record(assembly_id)
    detected = crispr_casoffinder.detect_cas_offinder()
    live = engine_validation.live_record("Cas-OFFinder")
    live_ok = bool(live and live.get("ok") and str(live.get("status") or "") == engine_validation.STATUS_LIVE_VALIDATED)
    if ready is None:
        return {
            "allowed": False,
            "reason": "Download/verify a reference assembly first. READY is required.",
            "engine_status": str((live or {}).get("status") or detected.get("status") or "not_installed"),
            "reference_status": "NOT_READY",
            "is_test_reference": False,
        }
    is_test = bool(ready.get("not_a_public_assembly"))
    if is_test:
        return {
            "allowed": False,
            "reason": (
                "This READY FASTA is the HelixScope TEST REFERENCE. It can be "
                "searched as a test, but it is not FULL GENOME / genome-wide."
            ),
            "engine_status": str((live or {}).get("status") or detected.get("status") or ""),
            "reference_status": "READY_TEST",
            "is_test_reference": True,
            "test_search_allowed": bool(detected.get("available")),
        }
    if not detected.get("available"):
        return {
            "allowed": False,
            "reason": (
                "Cas-OFFinder is not installed. HelixScope does not simulate hits."
            ),
            "engine_status": "not_installed",
            "reference_status": "READY",
            "is_test_reference": False,
        }
    if not live_ok:
        return {
            "allowed": False,
            "reason": (
                "Cas-OFFinder requires a working OpenCL runtime. The binary is "
                "DETECTED but not LIVE_VALIDATED. Run the TEST REFERENCE search "
                "after OpenCL is available. FULL GENOME stays locked."
            ),
            "engine_status": str((live or {}).get("status") or "detected"),
            "reference_status": "READY",
            "is_test_reference": False,
        }
    admission = _admit_ready_search(assembly_id, ready)
    from . import resource_admission

    allowed = admission.get("decision") == resource_admission.DECISION_PROCEED
    return {
        "allowed": allowed,
        "reason": (
            "READY public assembly and Cas-OFFinder LIVE_VALIDATED. "
            + str(admission.get("reason") or "")
            if allowed
            else str(admission.get("reason") or "RESOURCE_LIMIT: search cannot start.")
        ),
        "engine_status": engine_validation.STATUS_LIVE_VALIDATED,
        "reference_status": "READY",
        "is_test_reference": False,
        "resource_admission": admission,
        "resource_status": admission.get("decision"),
    }


def _admit_ready_search(assembly_id: str, ready: dict) -> dict:
    """Mede RAM vs tamanho da referencia READY. Nao inicia a busca.

    Args:
        assembly_id: id de catalogo.
        ready: ready_record.

    Returns:
        Dict de resource_admission.admit_reference_search.
    """
    from . import genome_fasta, resource_admission

    bases = 0
    fasta_bytes = 0
    fai_path = str(ready.get("fai_path") or "")
    if fai_path and os.path.isfile(fai_path):
        try:
            parsed = genome_fasta.load_fai(fai_path)
            bases = sum(
                int(row.get("length") or 0)
                for row in (parsed.get("contigs") or {}).values()
            )
        except genome_fasta.GenomeFastaError:
            bases = 0
    fasta_path = str(ready.get("fasta_path") or "")
    if fasta_path and os.path.isfile(fasta_path):
        fasta_bytes = int(os.path.getsize(fasta_path))
    try:
        if bases > 0:
            return resource_admission.admit_reference_search(
                assembly_id=assembly_id,
                bases=bases,
                engine_available=True,
            )
        return resource_admission.admit_reference_search(
            assembly_id=assembly_id,
            fasta_bytes=fasta_bytes,
            engine_available=True,
        )
    except ValueError as exc:
        return {
            "decision": resource_admission.DECISION_UNKNOWN,
            "assembly_id": str(assembly_id),
            "reason": str(exc),
            "blocker": "estimate",
        }


def sanitize_job_log(text: str) -> str:
    """Remove home path e tokens obvios de um log limitado.

    Args:
        text: Log cru.

    Returns:
        Texto sanitizado (teto 8000 caracteres).

    Raises:
        Nenhum.
    """
    payload = str(text or "")[:8000]
    home = os.path.expanduser("~")
    if home:
        payload = payload.replace(home, "<home>")
    for key in provenance.SECRET_ENV_KEYS:
        token = str(os.environ.get(key) or "")
        if token and token in payload:
            payload = payload.replace(token, "<redacted>")
    return payload


def submit_casoffinder_job(
    *,
    guide: Mapping[str, Any],
    assembly_id: str,
    max_mismatches: int = 3,
    dna_bulge: int = 0,
    rna_bulge: int = 0,
    include_nag: bool = True,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    max_hits: int = DEFAULT_MAX_HITS,
    threads: int = DEFAULT_THREADS,
    spawn: bool = True,
) -> dict:
    """Cria um job e, se spawn=True, arranca o worker.

    Args:
        guide: Guia SpCas9 validado.
        assembly_id: referencia READY.
        max_mismatches: 0-4.
        dna_bulge: 0 (default) ou >0 com semantica Cas-OFFinder.
        rna_bulge: 0 (default).
        include_nag: pesquisar NAG alem de NGG.
        timeout_s: timeout do worker.
        max_hits: teto honesto.
        threads: 1 por omissao; max 4.
        spawn: se False, so grava o job (testes).

    Returns:
        Status inicial.

    Raises:
        GenomeJobError: INVALID_INPUT / UNAVAILABLE.
    """
    ready = genome_store.ready_record(assembly_id)
    if ready is None:
        raise GenomeJobError(
            "Reference is not READY (verified). Genome-wide search refused.",
            "UNAVAILABLE",
        )
    spacer = str(guide.get("guide_sequence") or "")
    if len(spacer) != 20:
        raise GenomeJobError("Guide spacer must be 20 nt SpCas9.", "INVALID_INPUT")
    try:
        mm = int(max_mismatches)
        dna_b = int(dna_bulge)
        rna_b = int(rna_bulge)
        hits = int(max_hits)
        n_threads = int(threads)
        timeout = float(timeout_s)
    except (TypeError, ValueError) as exc:
        raise GenomeJobError("Numeric job parameters are invalid.", "INVALID_INPUT") from exc
    if mm < 0 or mm > 4 or dna_b < 0 or rna_b < 0 or hits < 1:
        raise GenomeJobError("Job parameter out of range.", "INVALID_INPUT")
    if dna_b or rna_b:
        raise GenomeJobError(
            "DNA/RNA bulges are UNAVAILABLE on Cas-OFFinder 2.4.1. "
            "Cas-OFFinder 3 is not used.",
            "UNAVAILABLE",
        )
    if n_threads < 1 or n_threads > MAX_THREADS:
        raise GenomeJobError(
            f"threads must be 1-{MAX_THREADS}. HelixScope does not default to all cores.",
            "INVALID_INPUT",
        )
    if timeout <= 0:
        raise GenomeJobError("timeout_s must be positive.", "INVALID_INPUT")
    gate = full_genome_search_preflight(assembly_id)
    if not ready.get("not_a_public_assembly") and not gate.get("allowed"):
        raise GenomeJobError(str(gate.get("reason") or "FULL GENOME search is blocked."), "UNAVAILABLE")
    job_id = new_job_id()
    engine_version = ""
    from . import crispr_casoffinder

    detected = crispr_casoffinder.detect_cas_offinder()
    engine_version = str(detected.get("version") or "")
    identity = cache_identity(
        guide_hash=provenance.sequence_digest(spacer),
        assembly_id=str(assembly_id),
        reference_sha256=str(ready["manifest"].get("file_sha256") or ""),
        engine="cas_offinder",
        engine_version=engine_version,
        max_mismatches=mm,
        dna_bulge=dna_b,
        rna_bulge=rna_b,
        include_nag=bool(include_nag),
    )
    cached = find_cached_job(identity)
    if cached is not None:
        status = dict(cached.get("status") or {})
        status["cache_hit"] = True
        status["job_id"] = cached.get("job_id")
        return status
    job = {
        "job_id": job_id,
        "status": JOB_QUEUED,
        "assembly_id": assembly_id,
        "fasta_path": ready["fasta_path"],
        "fai_path": ready["fai_path"],
        "reference_sha256": ready["manifest"].get("file_sha256"),
        "kind": ready["kind"],
        "not_a_public_assembly": ready["not_a_public_assembly"],
        "guide": dict(guide),
        "max_mismatches": mm,
        "dna_bulge": dna_b,
        "rna_bulge": rna_b,
        "include_nag": bool(include_nag),
        "timeout_s": timeout,
        "max_hits": hits,
        "threads": n_threads,
        "engine": "cas_offinder",
        "engine_version": engine_version,
        "cache_identity": identity,
        "created_at": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "progress": "indeterminate",
        "progress_note": (
            "Cas-OFFinder does not report percent complete. Status is RUNNING "
            "until the process exits."
        ),
        "stage": STAGE_PREPARING,
    }
    directory = job_dir(job_id)
    os.makedirs(directory, exist_ok=True)
    _atomic_write(os.path.join(directory, "job.json"), job)
    write_status(job_id, {"status": JOB_QUEUED, "reason": "Queued.", "pid": None})
    if spawn:
        _spawn_worker(job_id, n_threads)
    return load_status(job_id) or {"job_id": job_id, "status": JOB_QUEUED}


def _spawn_worker(job_id: str, threads: int) -> None:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(int(threads))
    env["HELIXSCOPE_JOB_ID"] = job_id
    log_path = os.path.join(job_dir(job_id), "worker.log")
    log_handle = open(log_path, "w", encoding="utf-8")
    popen_kwargs: dict[str, Any] = {
        "args": [sys.executable, "-m", "modules.genome_worker", "--job-id", job_id],
        "cwd": genome_store.repository_root(),
        "env": env,
        "stdout": log_handle,
        "stderr": subprocess.STDOUT,
        "shell": False,
    }
    if os.name == "nt":
        popen_kwargs["creationflags"] = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
    else:
        popen_kwargs["start_new_session"] = True
    proc = subprocess.Popen(**popen_kwargs)
    write_status(
        job_id,
        {
            "status": JOB_RUNNING,
            "reason": "Worker started.",
            "pid": proc.pid,
            "progress": "indeterminate",
        },
    )
    with open(os.path.join(job_dir(job_id), "pid"), "w", encoding="utf-8") as handle:
        handle.write(str(proc.pid))


def request_cancel(job_id: str) -> dict:
    """Pede cancelamento. Mata o PID se ainda estiver vivo.

    Args:
        job_id: id.

    Returns:
        Status.

    Raises:
        GenomeJobError: INVALID_INPUT.
    """
    directory = job_dir(job_id)
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "cancel"), "w", encoding="utf-8") as handle:
        handle.write("1")
    status = load_status(job_id) or {}
    pid = status.get("pid")
    if pid:
        _kill_pid(int(pid))
    return write_status(
        job_id,
        {
            **status,
            "status": JOB_CANCELLED,
            "reason": "Cancelled. Incomplete; not genome-wide; scores N/A.",
        },
    )


def cancel_requested(job_id: str) -> bool:
    """True se o ficheiro cancel existir.

    Args:
        job_id: id.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    try:
        return os.path.isfile(os.path.join(job_dir(job_id), "cancel"))
    except GenomeJobError:
        return False


def _kill_pid(pid: int) -> None:
    if pid <= 0:
        return
    try:
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/PID", str(pid), "/T", "/F"],
                capture_output=True,
                check=False,
                shell=False,
            )
        else:
            os.kill(pid, signal.SIGTERM)
    except OSError:
        return


def find_cached_job(identity: str) -> Optional[dict]:
    """Reutiliza um job COMPLETED_* nao truncado com a mesma identidade.

    Args:
        identity: SHA-256 de cache_identity().

    Returns:
        Dict job_id, status, result, ou None.

    Raises:
        Nenhum.
    """
    key = str(identity or "").strip()
    if not key:
        return None
    root = genome_store.jobs_root()
    if not os.path.isdir(root):
        return None
    names = sorted(os.listdir(root))[:400]
    reusable = {
        JOB_COMPLETED_FULL_REFERENCE,
        JOB_COMPLETED_TEST_REFERENCE,
    }
    for name in names:
        if len(name) != 32:
            continue
        job = load_job(name)
        if not job or str(job.get("cache_identity") or "") != key:
            continue
        status = load_status(name)
        result = load_result(name)
        if not status or not result:
            continue
        if str(status.get("status") or "") not in reusable:
            continue
        if result.get("truncated"):
            continue
        return {
            "job_id": name,
            "status": status,
            "result": result,
            "cache_hit": True,
        }
    return None


def load_result(job_id: str) -> Optional[dict]:
    """Le result.json.

    Args:
        job_id: id.

    Returns:
        dict ou None.

    Raises:
        Nenhum.
    """
    try:
        return read_json(os.path.join(job_dir(job_id), "result.json"))
    except GenomeJobError:
        return None
