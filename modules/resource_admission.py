"""Admissao de recursos antes de tarefas grandes (RAM, disco, engine).

O problema concreto: Cas-OFFinder 2.4.1 le a referencia FASTA completa para
memoria antes de procurar. Numa maquina com pouca RAM livre, iniciar uma busca
genome-wide nao produz um resultado tardio: produz paginacao severa, um sistema
instavel e, no limite, um processo morto pelo sistema operativo sem resultado.

Este modulo mede a memoria e o disco reais, estima a necessidade da tarefa a
partir de uma base declarada e devolve PROCEED ou RESOURCE_LIMIT. RESOURCE_LIMIT
e um resultado correto e final, nao uma falha a esconder: e a resposta honesta
"esta maquina nao suporta esta busca", que e preferivel a um numero global
produzido por um processo que nao terminou.

Nenhuma funcao aqui importa Streamlit. Nenhuma funcao aqui inicia a tarefa.
"""

from __future__ import annotations

import ctypes
import os
import platform
import shutil
from typing import Any, Mapping, Optional

from . import provenance

DECISION_PROCEED: str = "PROCEED"
DECISION_RESOURCE_LIMIT: str = "RESOURCE_LIMIT"
DECISION_UNKNOWN: str = "UNKNOWN"

MEMORY_SOURCE_WINDOWS: str = "GlobalMemoryStatusEx (Windows API)"
MEMORY_SOURCE_PROC: str = "/proc/meminfo (Linux kernel)"
MEMORY_SOURCE_SYSCONF: str = "os.sysconf SC_PAGE_SIZE * SC_AV_PHYS_PAGES"
MEMORY_SOURCE_UNAVAILABLE: str = "unavailable"

CAS_OFFINDER_BYTES_PER_BASE: float = 1.0
"""Bytes de memoria por base da referencia carregada por Cas-OFFinder 2.4.1.

Cas-OFFinder 2.2 e posteriores leem o conjunto de cromossomas para um buffer de
caracteres antes de comparar, um byte por base. Esta e a base declarada da
estimativa; nao e um valor medido em GRCh38 nesta maquina.
"""

CAS_OFFINDER_OVERHEAD_FACTOR: float = 1.35
"""Fator conservador para os buffers de comparacao e o runtime OpenCL.

Cobre a copia de trabalho do padrao, as listas de posicoes candidatas e o
contexto do device OpenCL. E uma margem de engenharia declarada, nao uma medicao
publicada.
"""

SAFETY_HEADROOM_BYTES: int = 512 * 1024 * 1024
"""Margem de RAM que deve continuar livre para o sistema e para a propria UI."""

FASTA_HEADER_AND_NEWLINE_OVERHEAD: float = 0.02
"""Fracao do ficheiro FASTA que nao e base (cabecalhos e fins de linha).

Um FASTA com 60 bases por linha gasta cerca de 1.7% dos bytes em fins de linha.
Usar 2% mantem a estimativa do lado conservador.
"""


class _MemoryStatusEx(ctypes.Structure):
    _fields_ = [
        ("dwLength", ctypes.c_ulong),
        ("dwMemoryLoad", ctypes.c_ulong),
        ("ullTotalPhys", ctypes.c_ulonglong),
        ("ullAvailPhys", ctypes.c_ulonglong),
        ("ullTotalPageFile", ctypes.c_ulonglong),
        ("ullAvailPageFile", ctypes.c_ulonglong),
        ("ullTotalVirtual", ctypes.c_ulonglong),
        ("ullAvailVirtual", ctypes.c_ulonglong),
        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
    ]


def system_memory() -> dict:
    """Mede a memoria fisica total e disponivel desta maquina.

    Args:
        Nenhum.

    Returns:
        Dict com total_bytes, available_bytes, load_percent, page_file_total_bytes,
        page_file_available_bytes, source e measured (bool). Valores None quando o
        sistema operativo nao expoe a informacao.

    Raises:
        Nenhum. Falha de medicao devolve measured=False, nunca um numero inventado.

    Nota:
        A memoria disponivel varia entre chamadas. Uma decisao de admissao usa a
        leitura do momento e deve ser refeita antes de iniciar a tarefa.
    """
    if platform.system() == "Windows":
        try:
            status = _MemoryStatusEx()
            status.dwLength = ctypes.sizeof(_MemoryStatusEx)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return {
                    "measured": True,
                    "source": MEMORY_SOURCE_WINDOWS,
                    "total_bytes": int(status.ullTotalPhys),
                    "available_bytes": int(status.ullAvailPhys),
                    "load_percent": int(status.dwMemoryLoad),
                    "page_file_total_bytes": int(status.ullTotalPageFile),
                    "page_file_available_bytes": int(status.ullAvailPageFile),
                }
        except (AttributeError, OSError):
            pass
    meminfo = "/proc/meminfo"
    if os.path.isfile(meminfo):
        try:
            values: dict[str, int] = {}
            with open(meminfo, "r", encoding="utf-8") as handle:
                for line in handle:
                    parts = line.split(":", 1)
                    if len(parts) != 2:
                        continue
                    number = parts[1].strip().split()
                    if number and number[0].isdigit():
                        values[parts[0].strip()] = int(number[0]) * 1024
            total = values.get("MemTotal")
            available = values.get("MemAvailable", values.get("MemFree"))
            if total and available is not None:
                return {
                    "measured": True,
                    "source": MEMORY_SOURCE_PROC,
                    "total_bytes": total,
                    "available_bytes": available,
                    "load_percent": (
                        int(round(100 * (total - available) / total)) if total else None
                    ),
                    "page_file_total_bytes": values.get("SwapTotal"),
                    "page_file_available_bytes": values.get("SwapFree"),
                }
        except OSError:
            pass
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        total_pages = os.sysconf("SC_PHYS_PAGES")
        avail_pages = os.sysconf("SC_AV_PHYS_PAGES")
        return {
            "measured": True,
            "source": MEMORY_SOURCE_SYSCONF,
            "total_bytes": int(page_size * total_pages),
            "available_bytes": int(page_size * avail_pages),
            "load_percent": None,
            "page_file_total_bytes": None,
            "page_file_available_bytes": None,
        }
    except (AttributeError, ValueError, OSError):
        return {
            "measured": False,
            "source": MEMORY_SOURCE_UNAVAILABLE,
            "total_bytes": None,
            "available_bytes": None,
            "load_percent": None,
            "page_file_total_bytes": None,
            "page_file_available_bytes": None,
            "reason": (
                "This platform does not expose physical memory through any of the "
                "supported interfaces. Memory admission is UNKNOWN and the caller "
                "must not claim the task is safe."
            ),
        }


def disk_free(path: str) -> dict:
    """Mede o espaco livre no volume de um caminho.

    Args:
        path: Caminho existente ou o seu diretorio pai mais proximo.

    Returns:
        Dict measured, path, total_bytes, free_bytes.

    Raises:
        Nenhum.
    """
    target = os.path.abspath(str(path or "."))
    while target and not os.path.isdir(target):
        parent = os.path.dirname(target)
        if parent == target:
            break
        target = parent
    try:
        usage = shutil.disk_usage(target)
    except OSError as exc:
        return {
            "measured": False,
            "path": target,
            "total_bytes": None,
            "free_bytes": None,
            "reason": f"Disk usage could not be read: {exc}",
        }
    return {
        "measured": True,
        "path": target,
        "total_bytes": int(usage.total),
        "free_bytes": int(usage.free),
    }


def estimate_bases_from_fasta_bytes(fasta_bytes: int) -> int:
    """Estima o numero de bases de um FASTA a partir do tamanho em bytes.

    Args:
        fasta_bytes: Tamanho do FASTA descomprimido.

    Returns:
        Numero estimado de bases (inteiro, arredondado para baixo).

    Raises:
        ValueError: Se fasta_bytes for negativo.

    Nota:
        Desconta uma fracao declarada para cabecalhos e fins de linha. Nao le o
        ficheiro; e uma conversao de tamanho, nao uma contagem.
    """
    size = int(fasta_bytes)
    if size < 0:
        raise ValueError("fasta_bytes cannot be negative.")
    return int(size * (1.0 - FASTA_HEADER_AND_NEWLINE_OVERHEAD))


def estimate_cas_offinder_memory(*, fasta_bytes: int = 0, bases: Optional[int] = None) -> dict:
    """Estima a RAM necessaria para uma busca Cas-OFFinder 2.4.1.

    Args:
        fasta_bytes: Tamanho do FASTA descomprimido; ignorado se bases for dado.
        bases: Numero de bases da referencia, quando conhecido pelo indice FAI.

    Returns:
        Dict com bases, bytes_per_base, overhead_factor, estimated_bytes,
        estimated_gib, basis e is_estimate=True.

    Raises:
        ValueError: Se nem fasta_bytes nem bases forem utilizaveis.

    Nota:
        E uma ESTIMATIVA de engenharia com base declarada, nao uma medicao do
        consumo real desta versao nesta maquina. A soma do FAI e a melhor entrada
        disponivel porque conta bases, nao bytes de ficheiro.
    """
    if bases is None:
        if int(fasta_bytes or 0) <= 0:
            raise ValueError("Provide fasta_bytes or bases to estimate memory.")
        base_count = estimate_bases_from_fasta_bytes(int(fasta_bytes))
        basis_input = f"{int(fasta_bytes):,} FASTA bytes minus declared header and newline overhead"
    else:
        base_count = int(bases)
        if base_count < 0:
            raise ValueError("bases cannot be negative.")
        basis_input = f"{base_count:,} bases counted from the FAI index"
    estimated = int(base_count * CAS_OFFINDER_BYTES_PER_BASE * CAS_OFFINDER_OVERHEAD_FACTOR)
    return {
        "is_estimate": True,
        "bases": base_count,
        "bytes_per_base": CAS_OFFINDER_BYTES_PER_BASE,
        "overhead_factor": CAS_OFFINDER_OVERHEAD_FACTOR,
        "estimated_bytes": estimated,
        "estimated_gib": round(estimated / (1024**3), 3),
        "basis": (
            "Cas-OFFinder 2.2 and later read the chromosome set into a character "
            f"buffer at one byte per base. Input: {basis_input}. A declared "
            f"overhead factor of {CAS_OFFINDER_OVERHEAD_FACTOR} covers the "
            "comparison buffers and the OpenCL device context. This is an "
            "engineering estimate with a stated basis, not a measurement of this "
            "build on this machine."
        ),
        "engine": "Cas-OFFinder 2.4.1",
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def admit_reference_search(
    *,
    assembly_id: str,
    fasta_bytes: int = 0,
    bases: Optional[int] = None,
    engine_available: bool = True,
    engine_reason: str = "",
    memory: Optional[Mapping[str, Any]] = None,
    allow_page_file: bool = False,
    headroom_bytes: int = SAFETY_HEADROOM_BYTES,
) -> dict:
    """Decide se uma busca na referencia pode comecar nesta maquina.

    Args:
        assembly_id: Assembly alvo, apenas para o registo.
        fasta_bytes: Tamanho do FASTA descomprimido.
        bases: Bases da referencia quando conhecidas pelo FAI.
        engine_available: False quando a engine nao esta instalada ou utilizavel.
        engine_reason: Motivo textual quando engine_available e False.
        memory: Leitura de memoria injetavel; por omissao mede o sistema.
        allow_page_file: True permite contar o ficheiro de paginacao disponivel.
            Por omissao False: paginar uma busca de 3 GiB deixa a maquina
            inutilizavel e a decisao passa a depender de swap, nao de RAM.
        headroom_bytes: RAM que deve continuar livre para o sistema e para a UI.

    Returns:
        Dict com decision (PROCEED, RESOURCE_LIMIT ou UNKNOWN), reason,
        estimate, memory, required_bytes, budget_bytes e assembly_id.

    Raises:
        ValueError: Se a estimativa nao puder ser calculada.

    Nota:
        UNKNOWN nunca deve ser tratado como PROCEED. Sem medicao de memoria, a
        decisao correta e pedir confirmacao explicita, nao assumir capacidade.
    """
    if not engine_available:
        return {
            "decision": DECISION_RESOURCE_LIMIT,
            "assembly_id": str(assembly_id),
            "reason": (
                engine_reason
                or "The search engine is not available, so no reference search can start."
            ),
            "blocker": "engine",
            "estimate": None,
            "memory": dict(memory) if memory else system_memory(),
            "required_bytes": None,
            "budget_bytes": None,
            "checked_at_utc": provenance.utc_now(),
        }
    estimate = estimate_cas_offinder_memory(fasta_bytes=fasta_bytes, bases=bases)
    reading = dict(memory) if memory is not None else system_memory()
    required = int(estimate["estimated_bytes"]) + int(headroom_bytes)
    if not reading.get("measured") or reading.get("available_bytes") is None:
        return {
            "decision": DECISION_UNKNOWN,
            "assembly_id": str(assembly_id),
            "reason": (
                "Physical memory could not be measured, so HelixScope cannot say "
                "whether this search fits. UNKNOWN is not PROCEED: starting anyway "
                "risks an out-of-memory kill with no result."
            ),
            "blocker": "memory_unmeasurable",
            "estimate": estimate,
            "memory": reading,
            "required_bytes": required,
            "budget_bytes": None,
            "checked_at_utc": provenance.utc_now(),
        }
    budget = int(reading["available_bytes"])
    page_file_used = False
    if allow_page_file and reading.get("page_file_available_bytes"):
        budget += int(reading["page_file_available_bytes"])
        page_file_used = True
    fits = budget >= required
    gib = 1024**3
    detail = (
        f"Estimated need {estimate['estimated_gib']} GiB plus "
        f"{headroom_bytes / gib:.2f} GiB of headroom = {required / gib:.2f} GiB. "
        f"Available budget {budget / gib:.2f} GiB"
        + (" including the page file" if page_file_used else " of physical RAM")
        + f" out of {int(reading['total_bytes'] or 0) / gib:.2f} GiB total."
    )
    return {
        "decision": DECISION_PROCEED if fits else DECISION_RESOURCE_LIMIT,
        "assembly_id": str(assembly_id),
        "reason": (
            f"{detail} The search can start."
            if fits
            else (
                f"{detail} This machine does not have enough free memory for a "
                "complete search of this reference. RESOURCE_LIMIT is the correct "
                "result: no genome-wide specificity score is produced, and none is "
                "guessed."
            )
        ),
        "blocker": "" if fits else "memory",
        "estimate": estimate,
        "memory": reading,
        "required_bytes": required,
        "budget_bytes": budget,
        "page_file_counted": page_file_used,
        "headroom_bytes": int(headroom_bytes),
        "checked_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def admit_download(
    *,
    assembly_id: str,
    compressed_bytes: int,
    store_path: str,
    keep_archive: bool = True,
) -> dict:
    """Decide se ha disco para descarregar e descomprimir uma assembly.

    Args:
        assembly_id: Assembly alvo.
        compressed_bytes: Tamanho do pacote comprimido anunciado pela fonte.
        store_path: Diretorio de destino, usado para medir o volume correto.
        keep_archive: True mantem o arquivo comprimido depois de descomprimir.

    Returns:
        Dict decision, reason, required_bytes, disk, assembly_id.

    Raises:
        ValueError: Se compressed_bytes nao for positivo.

    Nota:
        A necessidade e o comprimido mais o descomprimido, porque o ficheiro
        original existe em disco enquanto o gunzip corre. Assumir apenas o
        tamanho final enche o volume no meio da operacao.
    """
    size = int(compressed_bytes)
    if size <= 0:
        raise ValueError("compressed_bytes must be positive to plan a download.")
    decompressed = size * 3
    index_and_sidecar = 250 * 1024 * 1024
    required = decompressed + index_and_sidecar + (size if keep_archive else 0)
    disk = disk_free(store_path)
    if not disk.get("measured"):
        return {
            "decision": DECISION_UNKNOWN,
            "assembly_id": str(assembly_id),
            "reason": (
                "Free disk space could not be measured, so the download was not "
                "declared safe."
            ),
            "required_bytes": required,
            "disk": disk,
            "checked_at_utc": provenance.utc_now(),
        }
    free = int(disk["free_bytes"])
    fits = free >= required
    gib = 1024**3
    return {
        "decision": DECISION_PROCEED if fits else DECISION_RESOURCE_LIMIT,
        "assembly_id": str(assembly_id),
        "reason": (
            f"Need about {required / gib:.2f} GiB (compressed {size / gib:.2f} GiB, "
            f"decompressed estimated {decompressed / gib:.2f} GiB at a declared 3x "
            f"gzip ratio, plus index and sidecar). Free: {free / gib:.2f} GiB."
            + ("" if fits else " Not enough disk space; the download does not start.")
        ),
        "required_bytes": required,
        "decompressed_estimate_bytes": decompressed,
        "gzip_ratio_assumed": 3.0,
        "keep_archive": bool(keep_archive),
        "disk": disk,
        "checked_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
