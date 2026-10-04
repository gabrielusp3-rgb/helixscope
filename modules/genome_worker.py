"""Worker local de jobs genomicos (subprocesso separado da UI Streamlit).

Invocado como: python -m modules.genome_worker --job-id HEX
ou: python -m modules.genome_worker --download-assembly GRCh38.p14

Nao usa shell=True. Nao importa Streamlit. Crash do worker nao derruba a app;
o job fica ERROR se o PID morrer sem result.json.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback
from typing import Any, Mapping, Optional

from . import genome_jobs, genome_store, provenance


def job_int(job: Mapping[str, Any], key: str, default: int) -> int:
    """Le um inteiro de job.json. 0 e um valor valido, nao um default.

    Args:
        job: Payload persistido.
        key: Campo numerico.
        default: Usado so quando o campo esta ausente ou vazio.

    Returns:
        Inteiro. `max_mismatches` 0 permanece 0 (nao vira 3).

    Raises:
        ValueError: valor nao numerico.
    """
    value = job.get(key)
    if value is None or value == "":
        return int(default)
    return int(value)


def run_casoffinder_job(job_id: str) -> dict:
    """Executa Cas-OFFinder para um job gravado. Escreve result.json.

    Args:
        job_id: Identificador hex de 32 caracteres.

    Returns:
        Envelope de busca ou dict de erro.

    Raises:
        Nenhum. Falhas viram status ERROR no disco.
    """
    from . import crispr_casoffinder

    job = genome_jobs.load_job(job_id)
    if job is None:
        payload = {
            "status": genome_jobs.JOB_ERROR,
            "reason": "job.json is missing or unreadable.",
            "job_id": job_id,
        }
        try:
            genome_jobs.write_status(job_id, payload)
        except genome_jobs.GenomeJobError:
            pass
        return payload
    if genome_jobs.cancel_requested(job_id):
        status = genome_jobs.write_status(
            job_id,
            {
                "status": genome_jobs.JOB_CANCELLED,
                "reason": "Cancelled before Cas-OFFinder started. Not genome-wide.",
                "pid": os.getpid(),
            },
        )
        return status
    genome_jobs.write_status(
        job_id,
        {
            "status": genome_jobs.JOB_RUNNING,
            "reason": "Preparing Cas-OFFinder input. Progress is indeterminate.",
            "pid": os.getpid(),
            "progress": "indeterminate",
            "stage": genome_jobs.STAGE_PREPARING,
            "cache_identity": job.get("cache_identity"),
        },
    )
    genome_jobs.write_status(
        job_id,
        {
            "status": genome_jobs.JOB_RUNNING,
            "reason": "Cas-OFFinder running. Progress is indeterminate.",
            "pid": os.getpid(),
            "progress": "indeterminate",
            "stage": genome_jobs.STAGE_RUNNING,
            "cache_identity": job.get("cache_identity"),
        },
    )
    try:
        envelope = crispr_casoffinder.search_indexed_fasta(
            job["guide"],
            assembly_id=str(job.get("assembly_id") or ""),
            fasta_path=str(job.get("fasta_path") or ""),
            fai_path=str(job.get("fai_path") or ""),
            reference_sha256=str(job.get("reference_sha256") or ""),
            max_mismatches=job_int(job, "max_mismatches", 3),
            include_perfect=True,
            timeout_s=float(
                job["timeout_s"]
                if job.get("timeout_s") not in (None, "")
                else genome_jobs.DEFAULT_TIMEOUT_S
            ),
            max_hits=job_int(job, "max_hits", genome_jobs.DEFAULT_MAX_HITS),
            dna_bulge=job_int(job, "dna_bulge", 0),
            rna_bulge=job_int(job, "rna_bulge", 0),
            include_nag=bool(job.get("include_nag", True)),
            cancel_check=lambda: genome_jobs.cancel_requested(job_id),
        )
    except Exception as exc:
        payload = {
            "status": genome_jobs.JOB_ERROR,
            "reason": f"Worker exception: {exc}. Not converted to zero hits.",
            "pid": os.getpid(),
            "traceback": traceback.format_exc()[-4000:],
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
        genome_jobs._atomic_write(
            os.path.join(genome_jobs.job_dir(job_id), "result.json"),
            payload,
        )
        return genome_jobs.write_status(job_id, payload)
    genome_jobs.write_status(
        job_id,
        {
            "status": genome_jobs.JOB_RUNNING,
            "reason": "Validating and scoring Cas-OFFinder hits.",
            "pid": os.getpid(),
            "progress": "indeterminate",
            "stage": genome_jobs.STAGE_VALIDATING,
            "cache_identity": job.get("cache_identity"),
        },
    )
    genome_jobs._atomic_write(
        os.path.join(genome_jobs.job_dir(job_id), "result.json"),
        envelope if isinstance(envelope, dict) else {"status": genome_jobs.JOB_ERROR},
    )
    status_name = str(envelope.get("status") or genome_jobs.JOB_ERROR)
    stage = (
        genome_jobs.STAGE_COMPLETED
        if status_name
        in {
            genome_jobs.JOB_COMPLETED_FULL_REFERENCE,
            genome_jobs.JOB_COMPLETED_TEST_REFERENCE,
        }
        else str(envelope.get("status") or genome_jobs.STAGE_SCORING)
    )
    return genome_jobs.write_status(
        job_id,
        {
            "status": status_name,
            "reason": str(envelope.get("reason") or ""),
            "pid": os.getpid(),
            "genome_wide": bool(envelope.get("genome_wide")),
            "truncated": bool(envelope.get("truncated")),
            "verified_hit_count": envelope.get("verified_hit_count"),
            "cache_identity": job.get("cache_identity"),
            "progress": "complete" if stage == genome_jobs.STAGE_COMPLETED else "indeterminate",
            "stage": stage,
            "software_version": provenance.HELIXSCOPE_VERSION,
        },
    )


def run_download_assembly(assembly_id: str) -> dict:
    """Download explicito de uma assembly publica. Nao e genome-wide search.

    Args:
        assembly_id: id de catalogo (GRCh38.p14, T2T-CHM13v2.0, GRCm39).

    Returns:
        Manifesto ou dict de erro.

    Raises:
        Nenhum.
    """
    from . import genome_download

    try:
        return genome_download.download_and_install_assembly(assembly_id)
    except (genome_download.GenomeDownloadError, genome_store.GenomeStoreError) as exc:
        return {
            "status": getattr(exc, "category", "ERROR"),
            "reason": str(exc),
            "assembly_id": assembly_id,
        }


def main(argv: Optional[list[str]] = None) -> int:
    """Ponto de entrada CLI.

    Args:
        argv: Argumentos; default sys.argv[1:].

    Returns:
        Codigo de saida 0 ou 1.

    Raises:
        Nenhum.
    """
    parser = argparse.ArgumentParser(prog="modules.genome_worker")
    parser.add_argument("--job-id", default="", help="Cas-OFFinder job hex id")
    parser.add_argument(
        "--download-assembly",
        default="",
        help="Public catalog id to download (explicit; not startup)",
    )
    args = parser.parse_args(argv)
    job_id = str(args.job_id or os.environ.get("HELIXSCOPE_JOB_ID") or "").strip()
    assembly = str(args.download_assembly or "").strip()
    if job_id and assembly:
        print("Provide only one of --job-id or --download-assembly.", file=sys.stderr)
        return 1
    if job_id:
        result = run_casoffinder_job(job_id)
        print(str(result.get("status") or "ERROR"))
        return 0 if str(result.get("status") or "") not in {
            genome_jobs.JOB_ERROR,
        } else 1
    if assembly:
        result = run_download_assembly(assembly)
        print(str(result.get("status") or "ERROR"))
        return 0 if str(result.get("status") or "") == genome_store.STATUS_READY else 1
    print("Missing --job-id or --download-assembly.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
