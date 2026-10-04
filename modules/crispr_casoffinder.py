"""Execucao real de Cas-OFFinder, se o binario allowlisted estiver instalado.

Nao reimplementa Cas-OFFinder em Python. Se o executavel nao for encontrado,
o status e UNAVAILABLE. Nunca simula hits. Nunca aceita caminho de genoma do
usuario nem URL. A referencia e o FASTA ja validado, escrito num diretorio
temporario seguro.

CLI (Bae et al.; https://github.com/snugel/cas-offinder):
cas-offinder {input} {C|G|A} {output}

Este modulo usa apenas o dispositivo C (CPU). Nao pede GPU. Sem bulges (nao
passa parametros de DNA/RNA bulge). PAM SpCas9: duas corridas, padrao NGG e
padrao NAG. Hits sao reextraidos do contig HelixScope e recusados se o PAM ou
a sequencia nao baterem.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import time
from typing import Any, Callable, Mapping, Optional, Sequence

from . import (
    crispr,
    crispr_cfd,
    crispr_offtarget,
    crispr_reference,
    genome_fasta,
    genome_store,
    opencl_runtime,
    provenance,
    tool_detection,
    tool_paths,
)

CAS_OFFINDER_NAMES: tuple[str, ...] = ("cas-offinder", "cas-offinder.exe")
CAS_OFFINDER_ENV: str = "HELIXSCOPE_CAS_OFFINDER"
ALGORITHM_CAS_OFFINDER: str = "cas_offinder_spcas9_ngg_nag_no_bulge"
DEVICE_CPU: str = "C"
DEFAULT_TIMEOUT_S: float = 20.0
DEFAULT_MAX_HITS: int = 200
SPACER_NT: int = 20
CONTEXT_NT: int = 23
CAS_OFFINDER_STABLE_VERSION: str = "2.4.1"
CAS_OFFINDER_V3_PRODUCTION_READY: bool = False
CAS_OFFINDER_V3_WARNING: str = (
    "Cas-OFFinder 3 is not production ready yet; results can differ from "
    "Cas-OFFinder 2. For production use please use the latest Cas-OFFinder 2 "
    "instead (upstream README, confirmed 2026-08-28)."
)

RunFn = Callable[..., subprocess.CompletedProcess]


def detect_cas_offinder() -> dict:
    """Detecta o binario Cas-OFFinder allowlisted e a versao reportada.

    Args:
        Nenhum.

    Returns:
        Dict available, path (sanitizado para relatorio), version, status.

    Raises:
        Nenhum.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        CAS_OFFINDER_NAMES,
        extra_file_candidates=tool_paths.candidates_for(CAS_OFFINDER_NAMES),
        env_var=CAS_OFFINDER_ENV,
    )
    if executable is None:
        return {
            "id": "cas_offinder",
            "tool": "Cas-OFFinder",
            "available": False,
            "version": "",
            "version_status": "unavailable",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "path": "",
            "reason": (
                "Cas-OFFinder executable not installed. HelixScope does not "
                "reimplement it. Built-in PAM-index search remains available."
            ),
            "algorithm": ALGORITHM_CAS_OFFINDER,
            "bulges_supported": False,
        }
    version = _probe_version(executable)
    return {
        "id": "cas_offinder",
        "tool": "Cas-OFFinder",
        "available": True,
        "version": version,
        "version_status": "detected" if version else "unavailable",
        "status": tool_detection.classify_tool_record(
            available=True,
            path=executable,
            version=version,
        ),
        "path": executable,
        "report_path": tool_detection.sanitize_tool_path(executable),
        "reason": (
            f"Found allowlisted executable {os.path.basename(executable)}. "
            "Device C (CPU) is preferred when a CPU OpenCL ICD exists; device "
            "G is used if only GPU OpenCL is available. This Windows build of "
            "Cas-OFFinder still needs a working OpenCL ICD; missing OpenCL is "
            "ERROR, not simulated hits. GPU is not required. Cas-OFFinder 3 is "
            "not used (not production-ready per upstream)."
        ),
        "algorithm": ALGORITHM_CAS_OFFINDER,
        "bulges_supported": False,
        "stable_version_policy": cas_offinder_version_policy(),
    }


def cas_offinder_version_policy() -> dict:
    """Politica de versao estavel. Nao promove Cas-OFFinder 3.

    Args:
        Nenhum.

    Returns:
        Dict production_version, v3_production_ready, bulges.

    Raises:
        Nenhum.
    """
    return {
        "production_version": CAS_OFFINDER_STABLE_VERSION,
        "v3_production_ready": CAS_OFFINDER_V3_PRODUCTION_READY,
        "v3_warning": CAS_OFFINDER_V3_WARNING,
        "native_bulges": False,
        "bulge_status": "UNAVAILABLE",
        "bulge_reason": (
            "Stable Cas-OFFinder 2.4.1 does not provide production-native "
            "DNA/RNA bulge search. Native bulges are a Cas-OFFinder 3 beta "
            "feature. HelixScope does not switch engines to obtain bulges."
        ),
        "source": "https://github.com/snugel/cas-offinder",
        "confirmed_utc_date": "2026-08-28",
    }


def resolve_search_device(*, real_binary: bool) -> dict:
    """C se CPU OpenCL existir; G se so GPU. Mocks nao enumeram OpenCL.

    Args:
        real_binary: True quando subprocesso real (run_fn is None).

    Returns:
        Dict device, reason.

    Raises:
        Nenhum.
    """
    if not real_binary:
        return {
            "device": DEVICE_CPU,
            "reason": "Injected run_fn; OpenCL device is not used.",
            "cpu_available": False,
            "gpu_available": False,
        }
    return opencl_runtime.select_cas_offinder_device()


def full_reference_definition(assembly_id: str, fai: Mapping[str, Any]) -> dict:
    """Conjunto exacto de sequencias que constitui a referencia pesquisada.

    Args:
        assembly_id: id de catalogo.
        fai: indice carregado.

    Returns:
        Dict scope, contig_ids, n_contigs, total_nt.

    Raises:
        Nenhum.
    """
    from . import crispr_assemblies

    catalog = crispr_assemblies.catalog_by_id(assembly_id) or {}
    contigs = []
    total = 0
    for name, row in (fai.get("contigs") or {}).items():
        length = int(row.get("length") or 0)
        total += length
        contigs.append({"identifier": str(name), "length": length})
    is_test = bool(catalog.get("not_a_public_assembly"))
    return {
        "assembly_id": assembly_id,
        "assembly": catalog.get("assembly") or assembly_id,
        "accession": catalog.get("accession") or "",
        "is_public_assembly": not is_test and str(catalog.get("kind") or "") == "public_assembly",
        "scope_label": str(catalog.get("reference_scope") or catalog.get("notes") or ""),
        "mixes_primary_alt_patch_silently": False,
        "includes_all_sequences_in_downloaded_fasta": True,
        "n_contigs": len(contigs),
        "total_nt": total,
        "contigs": contigs,
        "note": (
            "The search input is the READY FASTA as indexed. HelixScope does "
            "not drop alternate loci, patches or unplaced sequences silently "
            "if they are present in that FASTA, and does not add sequences "
            "that are not in that FASTA."
        ),
    }


def build_engine_input_manifest(
    *,
    fasta_path: str,
    fai: Mapping[str, Any],
    reference_sha256: str,
    observed_sha256: str,
    device: str,
) -> dict:
    """Identidade dos ficheiros/contigs entregues ao engine.

    Args:
        fasta_path: FASTA READY.
        fai: indice.
        reference_sha256: hash declarado no job.
        observed_sha256: SHA-256 calculado agora.
        device: C/G/A.

    Returns:
        Dict manifest.

    Raises:
        Nenhum.
    """
    contig_ids = [str(name) for name in (fai.get("contigs") or {}).keys()]
    total = sum(int(row.get("length") or 0) for row in (fai.get("contigs") or {}).values())
    match = observed_sha256.lower() == str(reference_sha256 or "").lower()
    return {
        "fasta_basename": os.path.basename(os.path.abspath(fasta_path)),
        "n_input_files": 1,
        "n_contigs": len(contig_ids),
        "contig_ids": contig_ids,
        "total_nt": total,
        "declared_sha256": str(reference_sha256 or ""),
        "observed_sha256": observed_sha256,
        "sha256_match": match,
        "device": device,
        "cas_offinder_first_line": "READY FASTA path (single file, all indexed contigs)",
    }


def search_completeness_record(
    *,
    manifest: Mapping[str, Any],
    truncated: bool,
    timed_out: bool,
    cancelled: bool,
    tool_error: str,
    exit_ok: bool,
) -> dict:
    """Completion nao e so o exit code.

    Args:
        manifest: engine input manifest.
        truncated: hit cap.
        timed_out: timeout.
        cancelled: cancel.
        tool_error: stderr/reason.
        exit_ok: return code 0.

    Returns:
        Dict complete, reason, completion_status candidate.

    Raises:
        Nenhum.
    """
    sha_ok = bool(manifest.get("sha256_match"))
    n_contigs = int(manifest.get("n_contigs") or 0)
    if cancelled:
        return {
            "complete": False,
            "reason": "Cancelled. Incomplete; not genome-wide.",
            "proof": "cancel_flag",
        }
    if timed_out:
        return {
            "complete": False,
            "reason": "Timeout. Incomplete; not genome-wide.",
            "proof": "timeout",
        }
    if truncated:
        return {
            "complete": False,
            "reason": "Hit cap reached. RESOURCE_LIMIT; not COMPLETED_FULL_REFERENCE.",
            "proof": "hit_cap",
        }
    if tool_error or not exit_ok:
        return {
            "complete": False,
            "reason": str(tool_error or "Cas-OFFinder did not exit 0."),
            "proof": "engine_error",
        }
    if not sha_ok:
        return {
            "complete": False,
            "reason": "READY FASTA SHA-256 does not match the job declaration.",
            "proof": "checksum_mismatch",
        }
    if n_contigs < 1:
        return {
            "complete": False,
            "reason": "FAI contains no contigs.",
            "proof": "empty_reference",
        }
    return {
        "complete": True,
        "reason": (
            "Input identity verified (FASTA SHA-256 + FAI contig list) and "
            "Cas-OFFinder exited 0 without truncation. Cas-OFFinder 2.4.1 does "
            "not emit a per-contig completion log; HelixScope does not infer "
            "completion from hit count."
        ),
        "proof": "input_identity_and_clean_exit",
        "n_contigs_delivered": n_contigs,
    }


def search_with_cas_offinder(
    guide: Mapping[str, Any],
    reference: Mapping[str, Any],
    *,
    max_mismatches: int = 3,
    include_perfect: bool = True,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    max_hits: int = DEFAULT_MAX_HITS,
    target_hash: str = "",
    run_fn: Optional[RunFn] = None,
    which_info: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Corre Cas-OFFinder sobre o FASTA validado e revalida cada hit.

    Args:
        guide: Guia SpCas9 20 nt + NGG.
        reference: Referencia ja validada (nao um path do usuario).
        max_mismatches: 0 a 4.
        include_perfect: Inclui 0-mismatch (padrao True, igual a busca PAM-index).
        timeout_s: Timeout do subprocesso.
        max_hits: Teto; excesso e RESOURCE_LIMIT.
        target_hash: Hash do DNA de desenho para a chave de cache.
        run_fn: subprocess.run injetavel nos testes.
        which_info: Deteccao precomputada.

    Returns:
        Envelope no mesmo contrato da busca PAM-index. genome_wide e False.

    Raises:
        Nenhum. Falhas viram status no envelope.
    """
    started = time.perf_counter()
    info = dict(which_info) if which_info is not None else detect_cas_offinder()
    if not info.get("available"):
        envelope = crispr_offtarget._empty_envelope(
            status="UNAVAILABLE",
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason=str(info.get("reason") or "Cas-OFFinder is not installed."),
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        envelope["tool"] = info
        envelope["kind"] = "UNAVAILABLE"
        return envelope

    validation = crispr.validate_spcas9_guide(guide)
    if not validation["valid"]:
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="; ".join(validation["errors"]),
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        return envelope

    try:
        max_mm = int(max_mismatches)
        time_cap = float(timeout_s)
        hit_cap = int(max_hits)
    except (TypeError, ValueError):
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="max_mismatches, timeout_s and max_hits must be numeric.",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        return envelope
    if max_mm < 0 or max_mm > 4:
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mm,
            include_perfect=include_perfect,
            reason="max_mismatches must be between 0 and 4 for this method.",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        return envelope
    if hit_cap < 1:
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_RESOURCE_LIMIT,
            guide=guide,
            reference=reference,
            max_mismatches=max_mm,
            include_perfect=include_perfect,
            reason="max_hits must be at least 1.",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        return envelope

    spacer = validation["guide_sequence"]
    executable = str(info.get("path") or "")
    device_info = resolve_search_device(real_binary=run_fn is None)
    tmpdir = tempfile.mkdtemp(prefix="helixscope_casoffinder_")
    try:
        fasta_path = os.path.join(tmpdir, "reference.fa")
        fasta_text = crispr_reference.to_fasta_text(reference)
        with open(fasta_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(fasta_text)
        output_text, tool_error, timed_out = _run_queries(
            executable=executable,
            fasta_path=fasta_path,
            tmpdir=tmpdir,
            spacer=spacer,
            max_mm=max_mm,
            timeout_s=time_cap,
            run_fn=run_fn or subprocess.run,
            device=str(device_info.get("device") or DEVICE_CPU),
        )
        if timed_out:
            envelope = crispr_offtarget._empty_envelope(
                status=crispr_offtarget.JOB_TIMEOUT,
                guide=guide,
                reference=reference,
                max_mismatches=max_mm,
                include_perfect=include_perfect,
                reason=(
                    f"Cas-OFFinder exceeded {time_cap:.1f}s. Incomplete; not "
                    "converted to zero hits."
                ),
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            envelope["truncated"] = True
            return envelope
        if tool_error:
            envelope = crispr_offtarget._empty_envelope(
                status=crispr_offtarget.JOB_FAILED,
                guide=guide,
                reference=reference,
                max_mismatches=max_mm,
                include_perfect=include_perfect,
                reason=tool_error,
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            return envelope
        parsed = parse_cas_offinder_output(output_text)
        if parsed["status"] != "OK":
            envelope = crispr_offtarget._empty_envelope(
                status="PARSING_ERROR",
                guide=guide,
                reference=reference,
                max_mismatches=max_mm,
                include_perfect=include_perfect,
                reason=parsed["reason"],
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            return envelope
        hits, unverified, truncated = _verify_rows(
            rows=parsed["rows"],
            guide=guide,
            spacer=spacer,
            reference=reference,
            max_mm=max_mm,
            include_perfect=include_perfect,
            hit_cap=hit_cap,
        )
        if unverified:
            envelope = crispr_offtarget._empty_envelope(
                status="PARSING_ERROR",
                guide=guide,
                reference=reference,
                max_mismatches=max_mm,
                include_perfect=include_perfect,
                reason=(
                    f"{unverified} Cas-OFFinder line(s) could not be re-extracted "
                    "from the HelixScope FASTA. The catalogue is not reported "
                    "as complete. Not converted to zero hits."
                ),
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            envelope["unverified_lines"] = unverified
            return envelope
        status = crispr_offtarget.JOB_COMPLETED
        reason = ""
        if truncated:
            status = crispr_offtarget.JOB_RESOURCE_LIMIT
            reason = (
                f"Hit cap of {hit_cap} reached. Result is truncated and is not "
                "a complete off-target catalogue."
            )
        elif not hits:
            reason = (
                "0 verified hits under selected method/reference. This is not a "
                "safety claim and not 'no risk'."
            )
        return crispr_offtarget.build_search_envelope(
            status=status,
            reason=reason,
            guide=guide,
            reference=reference,
            hits=hits,
            max_mm=max_mm,
            include_perfect=include_perfect,
            truncated=truncated,
            hit_cap=hit_cap,
            time_cap=time_cap,
            target_hash=target_hash,
            elapsed_ms=_elapsed_ms(started),
            search_ms=_elapsed_ms(started),
            index_hash="",
            index_n_sites=None,
            index_build_ms=None,
            algorithm=ALGORITHM_CAS_OFFINDER,
            tool=_public_tool(info, device=str(device_info.get("device") or DEVICE_CPU)),
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def parse_cas_offinder_output(text: str) -> dict:
    """Interpreta o TSV/espacos do Cas-OFFinder.

    Campos esperados (documentacao do projeto):
    query chromosome position sequence strand mismatches

    Args:
        text: Stdout ou ficheiro de saida.

    Returns:
        Dict status OK|PARSING_ERROR, rows, reason.

    Raises:
        Nenhum.
    """
    rows: list[dict] = []
    for line_no, raw in enumerate(str(text or "").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6:
            return {
                "status": "PARSING_ERROR",
                "rows": [],
                "reason": (
                    f"Cas-OFFinder line {line_no} does not have 6 fields. "
                    "Output was not converted to hits."
                ),
            }
        dna_bulge = 0
        rna_bulge = 0
        bulge_tail = False
        if (
            len(parts) >= 8
            and parts[-4] in {"+", "-"}
            and parts[-5].replace("a", "A").replace("c", "C").replace("g", "G").replace("t", "T").replace("u", "U")
        ):
            try:
                int(parts[-2])
                int(parts[-1])
                int(parts[-6])
                bulge_tail = True
            except ValueError:
                bulge_tail = False
        if bulge_tail:
            try:
                dna_bulge = int(parts[-2])
                rna_bulge = int(parts[-1])
            except ValueError:
                bulge_tail = False
        if bulge_tail:
            mm_text = parts[-3]
            strand = parts[-4]
            sequence = parts[-5]
            pos_text = parts[-6]
            chrom_tokens = parts[1:-6]
        else:
            mm_text = parts[-1]
            strand = parts[-2]
            sequence = parts[-3]
            pos_text = parts[-4]
            chrom_tokens = parts[1:-4]
        if not chrom_tokens:
            return {
                "status": "PARSING_ERROR",
                "rows": [],
                "reason": f"Cas-OFFinder line {line_no} is missing a contig name.",
            }
        chrom = chrom_tokens[0]
        query = parts[0]
        try:
            position = int(pos_text)
            mismatches = int(mm_text)
        except ValueError:
            return {
                "status": "PARSING_ERROR",
                "rows": [],
                "reason": (
                    f"Cas-OFFinder line {line_no} has a non-integer position "
                    "or mismatch count."
                ),
            }
        if strand not in {"+", "-"}:
            return {
                "status": "PARSING_ERROR",
                "rows": [],
                "reason": f"Cas-OFFinder line {line_no} has strand '{strand}'.",
            }
        seq = sequence.upper().replace("U", "T")
        if len(seq) != CONTEXT_NT or set(seq) - set("ACGT"):
            return {
                "status": "PARSING_ERROR",
                "rows": [],
                "reason": (
                    f"Cas-OFFinder line {line_no} sequence is not a 23-mer DNA."
                ),
            }
        rows.append(
            {
                "query": query.upper().replace("U", "T"),
                "chromosome": chrom,
                "position_reported": position,
                "sequence": seq,
                "strand": strand,
                "mismatches_reported": mismatches,
                "dna_bulge": dna_bulge,
                "rna_bulge": rna_bulge,
                "line_no": line_no,
            }
        )
    return {"status": "OK", "rows": rows, "reason": ""}


def casoffinder_position_to_internal(
    position: int,
    *,
    origin: str,
) -> int:
    """Converte a coordenada reportada pelo Cas-OFFinder para 0-based.

    Args:
        position: Inteiro reportado.
        origin: '0-based' ou '1-based'. Unico ponto de conversao.

    Returns:
        Start 0-based do 23-mer no contig sense.

    Raises:
        ValueError: Se origin for desconhecido.
    """
    if origin == "0-based":
        return int(position)
    if origin == "1-based":
        return int(position) - 1
    raise ValueError(f"Unknown coordinate origin {origin!r}.")


def _run_queries(
    *,
    executable: str,
    fasta_path: str,
    tmpdir: str,
    spacer: str,
    max_mm: int,
    timeout_s: float,
    run_fn: RunFn,
    dna_bulge: int = 0,
    rna_bulge: int = 0,
    include_nag: bool = True,
    device: str = DEVICE_CPU,
) -> tuple[str, str, bool]:
    chunks: list[str] = []
    remaining = max(0.1, float(timeout_s))
    mode = str(device or DEVICE_CPU).strip().upper()[:1]
    if mode not in {"C", "G", "A"}:
        mode = DEVICE_CPU
    pams = [("NGG", "NGG")]
    if include_nag:
        pams.append(("NAG", "NAG"))
    for pam_pattern, query_pam in pams:
        started = time.monotonic()
        input_path = os.path.join(tmpdir, f"input_{pam_pattern}.txt")
        output_path = os.path.join(tmpdir, f"output_{pam_pattern}.txt")
        pattern = ("N" * SPACER_NT) + pam_pattern
        query = spacer + query_pam
        with open(input_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(f"{fasta_path}\n{pattern}\n")
            if int(dna_bulge) or int(rna_bulge):
                handle.write(f"{query} {max_mm} {int(dna_bulge)} {int(rna_bulge)}\n")
            else:
                handle.write(f"{query} {max_mm}\n")
        try:
            completed = run_fn(
                [executable, input_path, mode, output_path],
                capture_output=True,
                timeout=remaining,
                check=False,
                shell=False,
                cwd=tmpdir,
                text=True,
                env=opencl_runtime.subprocess_environ(),
            )
        except subprocess.TimeoutExpired:
            return "", "", True
        except OSError as exc:
            return "", f"Cas-OFFinder failed to start: {exc}.", False
        if int(getattr(completed, "returncode", 1) or 0) != 0:
            detail = (
                (getattr(completed, "stderr", None) or getattr(completed, "stdout", None) or "")
            )
            detail_text = str(detail).strip()[:300]
            return "", f"Cas-OFFinder exited {completed.returncode}. {detail_text}", False
        if os.path.isfile(output_path):
            with open(output_path, "r", encoding="utf-8") as handle:
                chunks.append(handle.read())
        remaining -= time.monotonic() - started
        if remaining <= 0:
            return "", "", True
    return "\n".join(chunks), "", False


def _verify_rows(
    *,
    rows: Sequence[Mapping[str, Any]],
    guide: Mapping[str, Any],
    spacer: str,
    reference: Mapping[str, Any],
    max_mm: int,
    include_perfect: bool,
    hit_cap: int,
) -> tuple[list[dict], int, bool]:
    hits: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    unverified = 0
    truncated = False
    for row in rows:
        site = locate_casoffinder_site(row, reference)
        if site is None:
            unverified += 1
            continue
        mismatches, positions = crispr_offtarget._mismatch_positions(
            spacer, str(site["spacer"])
        )
        if mismatches > max_mm:
            continue
        if mismatches == 0 and not include_perfect:
            continue
        try:
            reported = int(row.get("mismatches_reported"))
        except (TypeError, ValueError):
            unverified += 1
            continue
        if reported != mismatches:
            unverified += 1
            continue
        key = (str(site["contig_id"]), int(site["start_0based"]), str(site["strand"]))
        if key in seen:
            continue
        seen.add(key)
        try:
            hit = crispr_offtarget._verified_hit(
                guide=guide,
                spacer=spacer,
                site=site,
                reference=reference,
                mismatches=mismatches,
                mismatch_positions=positions,
            )
        except (crispr.CrisprError, crispr_reference.ReferenceError):
            unverified += 1
            continue
        hits.append(hit)
        if len(hits) >= hit_cap:
            truncated = True
            break
    return hits, unverified, truncated


def locate_casoffinder_site(
    row: Mapping[str, Any], reference: Mapping[str, Any]
) -> Optional[dict]:
    """Reconstroi o sitio no contig HelixScope a partir da linha do tool.

    Tenta origem 0-based e, se falhar, 1-based via helper explicito. Recusa
    intervalos que nao reproduzem a sequencia reportada.

    Args:
        row: Linha parseada.
        reference: Referencia validada.

    Returns:
        Dict no formato do indice PAM ou None.

    Raises:
        Nenhum.
    """
    contig_id = str(row.get("chromosome") or "")
    try:
        contig = crispr_reference.contig_by_id(reference, contig_id)
    except crispr_reference.ReferenceError:
        return None
    seq = str(contig.get("sequence") or "")
    reported = str(row.get("sequence") or "").upper().replace("U", "T")
    strand = str(row.get("strand") or "")
    if len(reported) != CONTEXT_NT:
        return None
    pos = int(row.get("position_reported"))
    start = None
    for origin in ("0-based", "1-based"):
        candidate = casoffinder_position_to_internal(pos, origin=origin)
        if candidate < 0 or candidate + CONTEXT_NT > len(seq):
            continue
        window = seq[candidate : candidate + CONTEXT_NT]
        if strand == "+" and window == reported:
            start = candidate
            break
        if strand == "-" and crispr._reverse_complement(window) == reported:
            start = candidate
            break
    if start is None:
        return None
    if strand == "+":
        spacer = reported[:SPACER_NT]
        pam = reported[SPACER_NT:]
        spacer_start = start
        spacer_end = start + SPACER_NT
        pam_start = start + SPACER_NT
        pam_end = start + CONTEXT_NT
        position = spacer_start + 1
    else:
        spacer = reported[:SPACER_NT]
        pam = reported[SPACER_NT:]
        spacer_start = start + 3
        spacer_end = start + CONTEXT_NT
        pam_start = start
        pam_end = start + 3
        position = spacer_start + 1
    if crispr._iupac_match(pam, "NGG"):
        pam_class = "NGG"
    elif crispr._iupac_match(pam, "NAG"):
        pam_class = "NAG_reduced_activity"
    else:
        return None
    return {
        "contig_id": contig_id,
        "strand": strand,
        "position_1based": position,
        "start_0based": spacer_start,
        "end_0based": spacer_end,
        "pam_start_0based": pam_start,
        "pam_end_0based": pam_end,
        "spacer": spacer,
        "pam_sequence": pam,
        "pam_class": pam_class,
    }


def _probe_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable],
            capture_output=True,
            timeout=5,
            check=False,
            shell=False,
            text=True,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = ((completed.stdout or "") + (completed.stderr or "")).strip()
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped[:80]
    return ""


def _public_tool(info: Mapping[str, Any], *, device: str = DEVICE_CPU) -> dict:
    return {
        "tool": "Cas-OFFinder",
        "available": bool(info.get("available")),
        "version": str(info.get("version") or ""),
        "version_status": str(info.get("version_status") or ""),
        "report_path": str(info.get("report_path") or os.path.basename(str(info.get("path") or ""))),
        "bulges_supported": False,
        "device": str(device or DEVICE_CPU),
    }


def search_indexed_fasta(
    guide: Mapping[str, Any],
    *,
    assembly_id: str,
    fasta_path: str,
    fai_path: str,
    reference_sha256: str,
    max_mismatches: int = 3,
    include_perfect: bool = True,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    max_hits: int = DEFAULT_MAX_HITS,
    dna_bulge: int = 0,
    rna_bulge: int = 0,
    include_nag: bool = True,
    run_fn: Optional[RunFn] = None,
    which_info: Optional[Mapping[str, Any]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> dict:
    """Cas-OFFinder sobre FASTA READY indexado. Reextrai cada hit via FAI.

    Args:
        guide: Guia SpCas9 20 nt + NGG.
        assembly_id: id do catalogo.
        fasta_path: Caminho no store.
        fai_path: Indice.
        reference_sha256: Hash do FASTA READY.
        max_mismatches: 0-4.
        include_perfect: inclui 0-mismatch.
        timeout_s: timeout.
        max_hits: teto; excesso RESOURCE_LIMIT.
        dna_bulge: 0 default. >0 usa a sintaxe oficial de bulge DNA.
        rna_bulge: 0 default.
        include_nag: segunda corrida NAG.
        run_fn: subprocess.run injetavel.
        which_info: deteccao precomputada.
        cancel_check: se True, CANCELLED.

    Returns:
        Envelope. genome_wide so se assembly publica READY e status
        COMPLETED_FULL_REFERENCE.

    Raises:
        Nenhum.

    Nota biologica:
        Bulges nao recebem CFD classico (Doench 2016 e Hamming 20+PAM).
        Pesquisa com bulge nao e Hamming-only.
    """
    started = time.perf_counter()
    info = dict(which_info) if which_info is not None else detect_cas_offinder()
    from . import crispr_assemblies, engine_validation

    catalog = crispr_assemblies.catalog_by_id(assembly_id) or {}
    is_test = bool(catalog.get("not_a_public_assembly")) or str(
        catalog.get("kind") or ""
    ) == "synthetic_test_reference"
    reference = {
        "identity_hash": reference_sha256,
        "scope": "test_reference" if is_test else "full_genome",
        "scope_label": (
            "HELIXSCOPE TEST REFERENCE (not a public assembly)"
            if is_test
            else f"READY assembly {catalog.get('assembly') or assembly_id}"
        ),
        "organism_declared": catalog.get("organism") or "",
        "assembly_declared": catalog.get("assembly") or assembly_id,
        "verified_assembly": "" if is_test else str(catalog.get("assembly") or ""),
        "assembly_identity_status": (
            "synthetic_test_reference"
            if is_test
            else "checksum_verified_local_file"
        ),
        "version_declared": catalog.get("release") or "",
        "accession": catalog.get("accession") or "",
        "contigs": [],
        "genome_wide": False,
    }
    if not info.get("available"):
        envelope = crispr_offtarget._empty_envelope(
            status="UNAVAILABLE",
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason=str(info.get("reason") or "Cas-OFFinder is not installed."),
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["algorithm"] = ALGORITHM_CAS_OFFINDER
        envelope["tool"] = info
        envelope["kind"] = "UNAVAILABLE"
        envelope["genome_wide"] = False
        return envelope
    if not genome_store.path_is_in_store(fasta_path):
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="FASTA path is not inside the HelixScope reference store.",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["genome_wide"] = False
        return envelope
    if cancel_check and cancel_check():
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_CANCELLED,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="Cancelled before Cas-OFFinder started.",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["genome_wide"] = False
        return envelope
    validation = crispr.validate_spcas9_guide(guide)
    if not validation["valid"]:
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="; ".join(validation["errors"]),
            elapsed_ms=_elapsed_ms(started),
        )
        return envelope
    if int(dna_bulge) or int(rna_bulge):
        policy = cas_offinder_version_policy()
        envelope = crispr_offtarget._empty_envelope(
            status="UNAVAILABLE",
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason=str(policy.get("bulge_reason") or "DNA/RNA bulges are UNAVAILABLE."),
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["genome_wide"] = False
        envelope["dna_bulge"] = int(dna_bulge)
        envelope["rna_bulge"] = int(rna_bulge)
        envelope["tool"] = _public_tool(info)
        envelope["kind"] = "UNAVAILABLE"
        return envelope
    spacer = validation["guide_sequence"]
    fai = genome_fasta.load_fai(fai_path)
    reference["contigs"] = [
        {"identifier": name, "length": row["length"]}
        for name, row in fai["contigs"].items()
    ]
    reference["full_reference"] = full_reference_definition(assembly_id, fai)
    try:
        observed_sha = genome_fasta.file_sha256(fasta_path)
    except genome_fasta.GenomeFastaError as exc:
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_FAILED,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason=f"Could not hash READY FASTA: {exc}",
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["genome_wide"] = False
        return envelope
    device_info = resolve_search_device(real_binary=run_fn is None)
    input_manifest = build_engine_input_manifest(
        fasta_path=fasta_path,
        fai=fai,
        reference_sha256=reference_sha256,
        observed_sha256=observed_sha,
        device=str(device_info.get("device") or DEVICE_CPU),
    )
    if not input_manifest.get("sha256_match"):
        envelope = crispr_offtarget._empty_envelope(
            status=crispr_offtarget.JOB_FAILED,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason=(
                "READY FASTA SHA-256 does not match the job declaration. "
                "Search refused."
            ),
            elapsed_ms=_elapsed_ms(started),
        )
        envelope["genome_wide"] = False
        envelope["engine_input_manifest"] = input_manifest
        return envelope
    executable = str(info.get("path") or "")
    tmpdir = tempfile.mkdtemp(prefix="helixscope_casoffinder_gw_")
    try:
        output_text, tool_error, timed_out = _run_queries(
            executable=executable,
            fasta_path=os.path.abspath(fasta_path),
            tmpdir=tmpdir,
            spacer=spacer,
            max_mm=int(max_mismatches),
            timeout_s=float(timeout_s),
            run_fn=run_fn or subprocess.run,
            dna_bulge=0,
            rna_bulge=0,
            include_nag=bool(include_nag),
            device=str(device_info.get("device") or DEVICE_CPU),
        )
        if timed_out:
            envelope = crispr_offtarget._empty_envelope(
                status=crispr_offtarget.JOB_TIMEOUT,
                guide=guide,
                reference=reference,
                max_mismatches=max_mismatches,
                include_perfect=include_perfect,
                reason="Cas-OFFinder exceeded timeout. Incomplete; not genome-wide.",
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["genome_wide"] = False
            envelope["truncated"] = True
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            envelope["engine_input_manifest"] = input_manifest
            envelope["search_completeness"] = search_completeness_record(
                manifest=input_manifest,
                truncated=False,
                timed_out=True,
                cancelled=False,
                tool_error="",
                exit_ok=False,
            )
            return envelope
        if tool_error:
            envelope = crispr_offtarget._empty_envelope(
                status=crispr_offtarget.JOB_FAILED,
                guide=guide,
                reference=reference,
                max_mismatches=max_mismatches,
                include_perfect=include_perfect,
                reason=tool_error,
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["genome_wide"] = False
            envelope["tool"] = _public_tool(info, device=str(device_info.get("device") or DEVICE_CPU))
            envelope["engine_input_manifest"] = input_manifest
            envelope["search_completeness"] = search_completeness_record(
                manifest=input_manifest,
                truncated=False,
                timed_out=False,
                cancelled=False,
                tool_error=tool_error,
                exit_ok=False,
            )
            envelope["opencl"] = {
                "preferred_device": str(device_info.get("device") or DEVICE_CPU),
                "selection": device_info,
                "gpu_required": False,
                "note": "Cas-OFFinder requires a working OpenCL runtime.",
            }
            return envelope
        parsed = parse_cas_offinder_output(output_text)
        if parsed["status"] != "OK":
            envelope = crispr_offtarget._empty_envelope(
                status="PARSING_ERROR",
                guide=guide,
                reference=reference,
                max_mismatches=max_mismatches,
                include_perfect=include_perfect,
                reason=parsed["reason"],
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["genome_wide"] = False
            envelope["engine_input_manifest"] = input_manifest
            return envelope
        hits, unverified, truncated = _verify_indexed_rows(
            rows=parsed["rows"],
            guide=guide,
            spacer=spacer,
            fasta_path=fasta_path,
            fai=fai,
            reference=reference,
            max_mm=int(max_mismatches),
            include_perfect=include_perfect,
            hit_cap=int(max_hits),
        )
        if unverified:
            envelope = crispr_offtarget._empty_envelope(
                status="PARSING_ERROR",
                guide=guide,
                reference=reference,
                max_mismatches=max_mismatches,
                include_perfect=include_perfect,
                reason=(
                    f"{unverified} Cas-OFFinder line(s) failed FAI re-extraction. "
                    "Not reported as a complete catalogue."
                ),
                elapsed_ms=_elapsed_ms(started),
            )
            envelope["genome_wide"] = False
            envelope["engine_input_manifest"] = input_manifest
            return envelope
        completeness = search_completeness_record(
            manifest=input_manifest,
            truncated=truncated,
            timed_out=False,
            cancelled=False,
            tool_error="",
            exit_ok=True,
        )
        if truncated:
            status = crispr_offtarget.JOB_RESOURCE_LIMIT
            reason = (
                f"Hit cap of {max_hits} reached. Truncated; not genome-wide; "
                "aggregate specificity is N/A."
            )
            genome_wide = False
        elif is_test:
            status = crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE
            reason = (
                "Completed search of the HelixScope TEST REFERENCE. "
                "Not GRCh38/T2T/GRCm39 genome-wide."
            )
            genome_wide = False
        elif completeness.get("complete"):
            status = crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE
            reason = (
                f"Completed Cas-OFFinder search of READY assembly "
                f"{catalog.get('assembly')} ({catalog.get('accession')}). "
                f"Scope: {catalog.get('reference_scope') or 'downloaded genomic FASTA'}."
            )
            genome_wide = True
        else:
            status = crispr_offtarget.JOB_FAILED
            reason = str(completeness.get("reason") or "Search was not proven complete.")
            genome_wide = False
        algorithm = ALGORITHM_CAS_OFFINDER
        envelope = crispr_offtarget.build_search_envelope(
            status=status,
            reason=reason,
            guide=guide,
            reference=reference,
            hits=hits,
            max_mm=int(max_mismatches),
            include_perfect=include_perfect,
            truncated=truncated,
            hit_cap=int(max_hits),
            time_cap=float(timeout_s),
            target_hash="",
            elapsed_ms=_elapsed_ms(started),
            search_ms=_elapsed_ms(started),
            index_hash=str(reference_sha256),
            index_n_sites=None,
            index_build_ms=None,
            algorithm=algorithm,
            tool=_public_tool(info, device=str(device_info.get("device") or DEVICE_CPU)),
            genome_wide=genome_wide,
        )
        envelope["dna_bulge"] = 0
        envelope["rna_bulge"] = 0
        envelope["mismatch_model"] = (
            "Hamming distance on 20 nt spacer (no bulges requested)"
        )
        envelope["engine_input_manifest"] = input_manifest
        envelope["search_completeness"] = completeness
        envelope["full_reference"] = reference.get("full_reference")
        envelope["opencl"] = {
            "preferred_device": str(device_info.get("device") or DEVICE_CPU),
            "selection": device_info,
            "gpu_required": False,
        }
        if genome_wide and (
            truncated
            or not completeness.get("complete")
            or status != crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE
        ):
            envelope["genome_wide"] = False
            genome_wide = False
        if run_fn is None and status in crispr_offtarget.COMPLETE_SEARCH_STATUSES:
            engine_validation.record_live_validation(
                "Cas-OFFinder",
                ok=True,
                version=str(info.get("version") or ""),
                details={
                    "n_hits": int(envelope.get("verified_hit_count") or 0),
                    "assembly_id": assembly_id,
                    "status": status,
                },
            )
        return envelope
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _verify_indexed_rows(
    *,
    rows: Sequence[Mapping[str, Any]],
    guide: Mapping[str, Any],
    spacer: str,
    fasta_path: str,
    fai: Mapping[str, Any],
    reference: Mapping[str, Any],
    max_mm: int,
    include_perfect: bool,
    hit_cap: int,
) -> tuple[list[dict], int, bool]:
    hits: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    unverified = 0
    truncated = False
    for row in rows:
        if int(row.get("dna_bulge") or 0) or int(row.get("rna_bulge") or 0):
            if len(str(row.get("sequence") or "")) != CONTEXT_NT:
                unverified += 1
                continue
        site = locate_casoffinder_site_indexed(row, fasta_path, fai)
        if site is None:
            unverified += 1
            continue
        mismatches, positions = crispr_offtarget._mismatch_positions(
            spacer, str(site["spacer"])
        )
        if mismatches > max_mm:
            continue
        if mismatches == 0 and not include_perfect:
            continue
        try:
            reported = int(row.get("mismatches_reported"))
        except (TypeError, ValueError):
            unverified += 1
            continue
        if reported != mismatches:
            unverified += 1
            continue
        key = (str(site["contig_id"]), int(site["start_0based"]), str(site["strand"]))
        if key in seen:
            continue
        seen.add(key)
        hit = _verified_hit_indexed(
            guide=guide,
            spacer=spacer,
            site=site,
            reference=reference,
            mismatches=mismatches,
            mismatch_positions=positions,
            fasta_path=fasta_path,
            fai=fai,
            dna_bulge=int(row.get("dna_bulge") or 0),
            rna_bulge=int(row.get("rna_bulge") or 0),
        )
        if hit is None:
            unverified += 1
            continue
        hits.append(hit)
        if len(hits) >= hit_cap:
            truncated = True
            break
    return hits, unverified, truncated


def locate_casoffinder_site_indexed(
    row: Mapping[str, Any],
    fasta_path: str,
    fai: Mapping[str, Any],
) -> Optional[dict]:
    """Reextrai o 23-mer via FAI. Nao carrega o contig inteiro.

    Args:
        row: Linha Cas-OFFinder.
        fasta_path: FASTA READY.
        fai: Indice parseado.

    Returns:
        Site dict ou None.

    Raises:
        Nenhum.
    """
    contig_id = str(row.get("chromosome") or "")
    fai_row = (fai.get("contigs") or {}).get(contig_id)
    if not fai_row:
        return None
    reported = str(row.get("sequence") or "").upper().replace("U", "T")
    strand = str(row.get("strand") or "")
    if len(reported) != CONTEXT_NT:
        return None
    pos = int(row.get("position_reported"))
    start = None
    for origin in ("0-based", "1-based"):
        candidate = casoffinder_position_to_internal(pos, origin=origin)
        if candidate < 0 or candidate + CONTEXT_NT > int(fai_row["length"]):
            continue
        try:
            window = genome_fasta.fetch_sequence(
                fasta_path, fai_row, candidate, candidate + CONTEXT_NT, strand="+"
            )
        except genome_fasta.GenomeFastaError:
            continue
        if strand == "+" and window == reported:
            start = candidate
            break
        if strand == "-" and crispr._reverse_complement(window) == reported:
            start = candidate
            break
    if start is None:
        return None
    if strand == "+":
        spacer = reported[:SPACER_NT]
        pam = reported[SPACER_NT:]
        spacer_start = start
        spacer_end = start + SPACER_NT
        pam_start = start + SPACER_NT
        pam_end = start + CONTEXT_NT
        position = spacer_start + 1
    else:
        spacer = reported[:SPACER_NT]
        pam = reported[SPACER_NT:]
        spacer_start = start + 3
        spacer_end = start + CONTEXT_NT
        pam_start = start
        pam_end = start + 3
        position = spacer_start + 1
    if crispr._iupac_match(pam, "NGG"):
        pam_class = "NGG"
    elif crispr._iupac_match(pam, "NAG"):
        pam_class = "NAG_reduced_activity"
    else:
        return None
    return {
        "contig_id": contig_id,
        "strand": strand,
        "position_1based": position,
        "start_0based": spacer_start,
        "end_0based": spacer_end,
        "pam_start_0based": pam_start,
        "pam_end_0based": pam_end,
        "spacer": spacer,
        "pam_sequence": pam,
        "pam_class": pam_class,
    }


def _verified_hit_indexed(
    *,
    guide: Mapping[str, Any],
    spacer: str,
    site: Mapping[str, Any],
    reference: Mapping[str, Any],
    mismatches: int,
    mismatch_positions: list[int],
    fasta_path: str,
    fai: Mapping[str, Any],
    dna_bulge: int,
    rna_bulge: int,
) -> Optional[dict]:
    fai_row = (fai.get("contigs") or {}).get(str(site["contig_id"]))
    if not fai_row:
        return None
    start = int(site["start_0based"])
    end = int(site["end_0based"])
    pam_start = int(site["pam_start_0based"])
    pam_end = int(site["pam_end_0based"])
    strand = str(site["strand"])
    try:
        sense_spacer = genome_fasta.fetch_sequence(
            fasta_path, fai_row, start, end, strand="+"
        )
        sense_pam = genome_fasta.fetch_sequence(
            fasta_path, fai_row, pam_start, pam_end, strand="+"
        )
    except genome_fasta.GenomeFastaError:
        return None
    if strand == "+":
        target_seq = sense_spacer
        pam_seq = sense_pam
    else:
        target_seq = crispr._reverse_complement(sense_spacer)
        pam_seq = crispr._reverse_complement(sense_pam)
    if target_seq != str(site["spacer"]) or pam_seq != str(site["pam_sequence"]):
        return None
    hsu = crispr.hsu_single_hit_score(spacer, target_seq)
    display = genome_fasta.display_interval_1based_inclusive(start, end)
    hit = {
        "guide_sequence": spacer,
        "guide_hash": provenance.sequence_digest(spacer),
        "reference_identity_hash": reference.get("identity_hash"),
        "organism_declared": reference.get("organism_declared") or "",
        "assembly_declared": reference.get("assembly_declared") or "",
        "version_declared": reference.get("version_declared") or "",
        "chromosome_or_contig": str(site["contig_id"]),
        "position": int(site["position_1based"]),
        "start_0based": start,
        "end_0based": end,
        "start_1based": display.get("start_1based"),
        "end_1based": display.get("end_1based"),
        "coordinate_display": display.get("convention"),
        "pam_start_0based": pam_start,
        "pam_end_0based": pam_end,
        "coordinate_system": "contig_sense_internal",
        "strand": strand,
        "target_sequence": target_seq,
        "PAM": pam_seq,
        "pam_sequence": pam_seq,
        "pam_class": site["pam_class"],
        "mismatches": mismatches,
        "mismatch_positions": mismatch_positions,
        "score": hsu["score"],
        "score_method": hsu["method"],
        "score_status": hsu["status"],
        "score_model": hsu["model"],
        "score_version": hsu["version"],
        "status": "VERIFIED",
        "dna_bulge": int(dna_bulge),
        "rna_bulge": int(rna_bulge),
        "bulges_supported": bool(dna_bulge or rna_bulge),
        "site_class": "perfect_site" if mismatches == 0 else "mismatch_site",
        "feature_annotation": "UNKNOWN",
        "annotation_note": (
            "No GFF/GTF loaded for this assembly. Coordinates are not a gene claim."
        ),
    }
    if dna_bulge or rna_bulge:
        hit["cfd_score"] = None
        hit["cfd_status"] = "UNAVAILABLE"
        hit["cfd_reason"] = (
            "Doench 2016 CFD is not defined for bulge alignments. Score is N/A."
        )
        return hit
    return crispr_cfd.attach_cfd_to_hit(hit, str(guide.get("pam_sequence") or ""))


def _elapsed_ms(started: float) -> float:
    return round((time.perf_counter() - started) * 1000.0, 3)
