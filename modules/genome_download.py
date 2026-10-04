"""Download explicito de assemblies NCBI para o store local.

HTTPS allowlisted. Resume via Range quando o servidor permite. Ficheiro
.partial nunca e READY. Checksum NCBI MD5 e obrigatorio antes de descomprimir.
Nenhuma funcao importa Streamlit. Nao arranca no import.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import zipfile
from typing import Callable, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener, urlopen

from . import crispr_assemblies, genome_fasta, genome_store, provenance

ALLOWED_DOWNLOAD_URLS: tuple[tuple[str, str], ...] = (
    ("ftp.ncbi.nlm.nih.gov", "/genomes/all/"),
    ("github.com", "/snugel/cas-offinder/releases/"),
    ("github.com", "/iqtree/iqtree3/releases/"),
    ("objects.githubusercontent.com", "/"),
    ("release-assets.githubusercontent.com", "/"),
)

USER_AGENT: str = "HelixScope/genome-reference"
HTTP_TIMEOUT_S: float = 60.0
MAX_UNCOMPRESSED_GENOME_BYTES: int = 20 * 1024 * 1024 * 1024
CHUNK: int = 1024 * 1024


class GenomeDownloadError(ValueError):
    """Falha de download ou verificacao.

    Attributes:
        category: INVALID_INPUT, NETWORK_ERROR, RESOURCE_LIMIT, CORRUPT, ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def url_is_allowed(url: str) -> bool:
    """HTTPS + host/prefixo allowlisted.

    Args:
        url: URL absoluta.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    parsed = urlparse(str(url or ""))
    if parsed.scheme != "https":
        return False
    host = (parsed.hostname or "").lower()
    path = parsed.path or ""
    for allowed_host, prefix in ALLOWED_DOWNLOAD_URLS:
        if host == allowed_host and path.startswith(prefix):
            return True
    return False


class _AllowlistRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not url_is_allowed(newurl):
            raise GenomeDownloadError(
                "Refusing HTTP redirect off the genome download allowlist.",
                "INVALID_INPUT",
            )
        return HTTPRedirectHandler.redirect_request(
            self, req, fp, code, msg, headers, newurl
        )


def _opener():
    return build_opener(_AllowlistRedirect)


def download_to_partial(
    url: str,
    partial_path: str,
    *,
    expected_md5: str = "",
    max_bytes: int = 0,
    urlopen_fn=None,
    progress_fn: Optional[Callable[[int], None]] = None,
) -> dict:
    """Descarrega para .partial com resume Range. Nao promove a READY.

    Args:
        url: HTTPS allowlisted.
        partial_path: Destino .partial.
        expected_md5: Se o ficheiro ja estiver completo, confere MD5.
        max_bytes: Teto; 0 usa 8x estimated default 8 GB.
        urlopen_fn: urlopen injetavel.
        progress_fn: Callback com bytes totais no destino.

    Returns:
        Dict path, bytes, resumed, md5.

    Raises:
        GenomeDownloadError: categorias classificadas.
    """
    if not url_is_allowed(url):
        raise GenomeDownloadError("Download URL is not on the allowlist.", "INVALID_INPUT")
    dest = os.path.abspath(partial_path)
    if not dest.endswith(genome_store.PARTIAL_SUFFIX):
        raise GenomeDownloadError(
            "Downloads must target a .partial file until checksum passes.",
            "INVALID_INPUT",
        )
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cap = int(max_bytes) if int(max_bytes) > 0 else 8 * 1024 * 1024 * 1024
    existing = os.path.getsize(dest) if os.path.isfile(dest) else 0
    headers = {"User-Agent": USER_AGENT}
    if existing:
        headers["Range"] = f"bytes={existing}-"
    request = Request(url, headers=headers, method="GET")
    opener = urlopen_fn or _opener().open
    resumed = False
    try:
        response = opener(request, timeout=HTTP_TIMEOUT_S)
    except HTTPError as exc:
        if existing and getattr(exc, "code", None) == 416:
            observed = genome_fasta.file_md5(dest)
            return {
                "path": dest,
                "bytes": existing,
                "resumed": True,
                "md5": observed,
                "complete": True,
            }
        if existing and getattr(exc, "code", None) == 200:
            existing = 0
            request = Request(url, headers={"User-Agent": USER_AGENT}, method="GET")
            try:
                response = opener(request, timeout=HTTP_TIMEOUT_S)
            except (HTTPError, URLError, TimeoutError, OSError) as retry_exc:
                raise GenomeDownloadError(
                    f"Download failed: {retry_exc}",
                    "NETWORK_ERROR",
                ) from retry_exc
        else:
            raise GenomeDownloadError(f"Download HTTP error: {exc}", "NETWORK_ERROR") from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise GenomeDownloadError(f"Download failed: {exc}", "NETWORK_ERROR") from exc
    mode = "ab" if existing else "wb"
    resumed = bool(existing)
    written = existing
    try:
        with open(dest, mode) as handle:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                written += len(chunk)
                if written > cap:
                    raise GenomeDownloadError(
                        "Download exceeded the byte cap. File remains PARTIAL.",
                        "RESOURCE_LIMIT",
                    )
                handle.write(chunk)
                if progress_fn is not None:
                    progress_fn(written)
    finally:
        try:
            response.close()
        except Exception:
            pass
    observed = genome_fasta.file_md5(dest)
    if expected_md5 and observed.lower() != str(expected_md5).lower():
        raise GenomeDownloadError(
            "Downloaded file MD5 does not match NCBI md5checksums.txt. "
            "Kept as PARTIAL/CORRUPT, not READY.",
            "CORRUPT",
        )
    return {
        "path": dest,
        "bytes": written,
        "resumed": resumed,
        "md5": observed,
        "complete": True,
    }


def decompress_gzip_to_fasta(gz_path: str, fasta_path: str) -> dict:
    """Descomprime gzip em streaming para FASTA. Teto 20 GiB.

    Args:
        gz_path: .gz completo.
        fasta_path: Destino .fa.

    Returns:
        Dict fasta_path, uncompressed_bytes.

    Raises:
        GenomeDownloadError: RESOURCE_LIMIT / PARSING_ERROR.
    """
    dest = os.path.abspath(fasta_path)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    total = 0
    try:
        with gzip.open(gz_path, "rb") as src, open(dest, "wb") as out:
            while True:
                chunk = src.read(CHUNK)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_UNCOMPRESSED_GENOME_BYTES:
                    out.close()
                    try:
                        os.remove(dest)
                    except OSError:
                        pass
                    raise GenomeDownloadError(
                        "Decompressed FASTA exceeded 20 GiB. Treated as a bomb.",
                        "RESOURCE_LIMIT",
                    )
                out.write(chunk)
    except GenomeDownloadError:
        raise
    except OSError as exc:
        raise GenomeDownloadError(
            f"gzip could not be decompressed: {exc}",
            "PARSING_ERROR",
        ) from exc
    return {"fasta_path": dest, "uncompressed_bytes": total}


def install_public_assembly_from_gz(
    assembly_id: str,
    gz_path: str,
    *,
    observed_md5: str,
) -> dict:
    """Apos MD5 do .gz, descomprime, indexa e marca READY.

    Args:
        assembly_id: id de catalogo publico.
        gz_path: genomic.fna.gz verificado.
        observed_md5: MD5 ja calculado.

    Returns:
        Manifesto READY.

    Raises:
        GenomeDownloadError / genome_store.GenomeStoreError.
    """
    catalog = crispr_assemblies.catalog_by_id(assembly_id)
    if catalog is None or catalog.get("kind") == "synthetic_test_reference":
        raise GenomeDownloadError("Unknown public assembly id.", "INVALID_INPUT")
    expected = str(catalog.get("official_compressed_md5") or "")
    if not expected:
        raise GenomeDownloadError(
            "Catalog is missing the official NCBI MD5. Refusing READY.",
            "CORRUPT",
        )
    if observed_md5.lower() != expected.lower():
        genome_store.mark_corrupt(
            assembly_id,
            reason="Compressed MD5 does not match NCBI md5checksums.txt.",
        )
        raise GenomeDownloadError("Official MD5 mismatch.", "CORRUPT")
    fasta_path = os.path.join(genome_store.assembly_dir(assembly_id), genome_store.FASTA_NAME)
    decompress_gzip_to_fasta(gz_path, fasta_path)
    return genome_store.finalize_ready_from_fasta(
        assembly_id,
        fasta_path,
        kind="public_assembly",
        not_a_public_assembly=False,
        official_compressed_md5=expected,
        compressed_md5_observed=observed_md5,
    )


def download_and_install_assembly(
    assembly_id: str,
    *,
    urlopen_fn=None,
    progress_fn=None,
) -> dict:
    """Download explicito: .partial → MD5 → decompress → FAI → READY.

    Args:
        assembly_id: GRCh38.p14 / T2T-CHM13v2.0 / GRCm39.
        urlopen_fn: urlopen injetavel.
        progress_fn: callback de bytes.

    Returns:
        Manifesto READY.

    Raises:
        GenomeDownloadError: disco, rede, checksum.
    """
    catalog = crispr_assemblies.catalog_by_id(assembly_id)
    if catalog is None or catalog.get("kind") == "synthetic_test_reference":
        raise GenomeDownloadError("Assembly is not a downloadable public catalog entry.", "INVALID_INPUT")
    need = int(catalog.get("estimated_compressed_bytes") or 0) * 3
    disk = genome_store.disk_status(need_bytes=need)
    if not disk.get("ok"):
        raise GenomeDownloadError(str(disk.get("reason") or "Low disk."), "RESOURCE_LIMIT")
    url = str(catalog.get("download_url") or "")
    expected = str(catalog.get("official_compressed_md5") or "")
    if not url or not expected:
        raise GenomeDownloadError("Catalog download URL or official MD5 is missing.", "INVALID_INPUT")
    directory = genome_store.assembly_dir(assembly_id)
    os.makedirs(directory, exist_ok=True)
    partial = os.path.join(directory, catalog["official_fasta_name"] + genome_store.PARTIAL_SUFFIX)
    genome_store.write_manifest(
        assembly_id,
        {
            **genome_store.empty_manifest(assembly_id),
            "status": genome_store.STATUS_DOWNLOADING,
            "downloaded_at": provenance.utc_now(),
        },
    )
    try:
        result = download_to_partial(
            url,
            partial,
            expected_md5=expected,
            max_bytes=int(catalog.get("estimated_compressed_bytes") or 0) * 3
            or 8 * 1024 * 1024 * 1024,
            urlopen_fn=urlopen_fn,
            progress_fn=progress_fn,
        )
    except GenomeDownloadError as exc:
        if exc.category == "CORRUPT":
            genome_store.mark_corrupt(assembly_id, reason=str(exc))
        else:
            genome_store.mark_partial(assembly_id, reason=str(exc))
        raise
    genome_store.write_manifest(
        assembly_id,
        {
            **genome_store.empty_manifest(assembly_id),
            "status": genome_store.STATUS_VERIFYING,
            "compressed_md5_observed": result["md5"],
        },
    )
    gz_final = os.path.join(directory, catalog["official_fasta_name"])
    if os.path.isfile(gz_final):
        os.remove(gz_final)
    shutil.move(result["path"], gz_final)
    return install_public_assembly_from_gz(
        assembly_id, gz_final, observed_md5=str(result["md5"])
    )


CAS_OFFINDER_RELEASE: str = "2.4.1"
CAS_OFFINDER_LICENSE: str = "BSD-3-Clause (upstream cas-offinder)"
CAS_OFFINDER_SOURCE: str = "https://github.com/snugel/cas-offinder"
CAS_OFFINDER_ZIP_URL: str = (
    "https://github.com/snugel/cas-offinder/releases/download/2.4.1/"
    "cas-offinder_windows_x86-64.zip"
)
CAS_OFFINDER_ZIP_SHA256: str = (
    "786a7fc68c474b98f4072836b469edf5074bf032f7cd059c9b19cbf1c05afb87"
)
"""SHA-256 of cas-offinder_windows_x86-64.zip from GitHub release 2.4.1, hashed locally 2026-08-28."""


def extract_zip_safe(zip_path: str, dest_dir: str, *, max_files: int = 40) -> list[str]:
    """Extrai um zip recusando path traversal e teto de ficheiros.

    Args:
        zip_path: Arquivo zip.
        dest_dir: Destino.
        max_files: Teto.

    Returns:
        Lista de caminhos extraidos.

    Raises:
        GenomeDownloadError: INVALID_INPUT / RESOURCE_LIMIT.
    """
    dest = os.path.abspath(dest_dir)
    os.makedirs(dest, exist_ok=True)
    extracted: list[str] = []
    try:
        archive = zipfile.ZipFile(zip_path)
    except zipfile.BadZipFile as exc:
        raise GenomeDownloadError("Zip is not a valid archive.", "PARSING_ERROR") from exc
    try:
        infos = archive.infolist()
        if len(infos) > int(max_files):
            raise GenomeDownloadError(
                "Zip has too many members. Refusing extraction.",
                "RESOURCE_LIMIT",
            )
        for info in infos:
            name = str(info.filename or "").replace("\\", "/")
            if not name or name.endswith("/"):
                continue
            parts = [p for p in name.split("/") if p not in {"", "."}]
            if any(part == ".." for part in parts) or name.startswith("/"):
                raise GenomeDownloadError(
                    "Zip member path is not allowed (traversal).",
                    "INVALID_INPUT",
                )
            target = os.path.abspath(os.path.join(dest, *parts))
            try:
                common = os.path.commonpath([dest, target])
            except ValueError as exc:
                raise GenomeDownloadError(
                    "Zip member resolves outside the destination.",
                    "INVALID_INPUT",
                ) from exc
            if common != dest:
                raise GenomeDownloadError(
                    "Zip member resolves outside the destination.",
                    "INVALID_INPUT",
                )
            if info.file_size > 50 * 1024 * 1024:
                raise GenomeDownloadError(
                    "Zip member exceeds 50 MiB. Refusing extraction.",
                    "RESOURCE_LIMIT",
                )
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with archive.open(info) as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)
            extracted.append(target)
    finally:
        archive.close()
    return extracted


def download_cas_offinder_windows(
    *,
    urlopen_fn=None,
    dest_dir: str = "",
) -> dict:
    """Descarrega o zip oficial Cas-OFFinder 2.4.1 e extrai sob tools/.

    Args:
        urlopen_fn: urlopen injetavel.
        dest_dir: Destino; default tools/cas-offinder.

    Returns:
        Dict path, version, sha256.

    Raises:
        GenomeDownloadError: checksum / rede / zip.
    """
    expected = str(CAS_OFFINDER_ZIP_SHA256 or "").strip().lower()
    if not expected:
        raise GenomeDownloadError(
            "Cas-OFFinder zip SHA-256 is not recorded. Refusing install.",
            "CORRUPT",
        )
    dest = dest_dir or os.path.join(genome_store.repository_root(), "tools", "cas-offinder")
    os.makedirs(dest, exist_ok=True)
    partial = os.path.join(dest, "cas-offinder_windows_x86-64.zip.partial")
    result = download_to_partial(
        CAS_OFFINDER_ZIP_URL,
        partial,
        expected_md5="",
        max_bytes=5 * 1024 * 1024,
        urlopen_fn=urlopen_fn,
    )
    zip_final = os.path.join(dest, "cas-offinder_windows_x86-64.zip")
    if os.path.isfile(zip_final):
        os.remove(zip_final)
    shutil.move(result["path"], zip_final)
    observed = genome_fasta.file_sha256(zip_final)
    if observed.lower() != expected:
        try:
            os.remove(zip_final)
        except OSError:
            pass
        raise GenomeDownloadError(
            "Cas-OFFinder zip SHA-256 does not match the recorded official asset.",
            "CORRUPT",
        )
    extracted = extract_zip_safe(zip_final, dest)
    exe = ""
    for path in extracted:
        if os.path.basename(path).lower() in {"cas-offinder.exe", "cas-offinder"}:
            exe = path
            break
    return {
        "status": "INSTALLED" if exe else "EXTRACTED_NO_EXECUTABLE",
        "zip_sha256": observed,
        "zip_path": zip_final,
        "executable": exe,
        "version_expected": CAS_OFFINDER_RELEASE,
        "source": CAS_OFFINDER_SOURCE,
        "license": CAS_OFFINDER_LICENSE,
        "extracted": [os.path.basename(p) for p in extracted],
    }


def spawn_download_assembly(assembly_id: str) -> dict:
    """Arranca download em subprocesso. Nao bloqueia a UI.

    Args:
        assembly_id: id publico.

    Returns:
        Dict pid, assembly_id, status DOWNLOADING.

    Raises:
        GenomeDownloadError: INVALID_INPUT / RESOURCE_LIMIT.
    """
    catalog = crispr_assemblies.catalog_by_id(assembly_id)
    if catalog is None or catalog.get("kind") == "synthetic_test_reference":
        raise GenomeDownloadError("Assembly is not a downloadable public catalog entry.", "INVALID_INPUT")
    need = int(catalog.get("estimated_compressed_bytes") or 0) * 3
    disk = genome_store.disk_status(need_bytes=need)
    if not disk.get("ok"):
        raise GenomeDownloadError(str(disk.get("reason") or "Low disk."), "RESOURCE_LIMIT")
    ident = str(catalog["id"])
    genome_store.write_manifest(
        ident,
        {
            **genome_store.empty_manifest(ident),
            "status": genome_store.STATUS_DOWNLOADING,
            "downloaded_at": provenance.utc_now(),
            "reason": "Download worker started. Opening CRISPR does not start this.",
        },
    )
    log_path = os.path.join(genome_store.assembly_dir(ident), "download.log")
    log_handle = open(log_path, "w", encoding="utf-8")
    popen_kwargs: dict = {
        "args": [
            sys.executable,
            "-m",
            "modules.genome_worker",
            "--download-assembly",
            ident,
        ],
        "cwd": genome_store.repository_root(),
        "stdout": log_handle,
        "stderr": subprocess.STDOUT,
        "shell": False,
    }
    if os.name == "nt":
        popen_kwargs["creationflags"] = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
    else:
        popen_kwargs["start_new_session"] = True
    proc = subprocess.Popen(**popen_kwargs)
    pid_path = os.path.join(genome_store.assembly_dir(ident), "download.pid")
    with open(pid_path, "w", encoding="utf-8") as handle:
        handle.write(str(proc.pid))
    return {
        "assembly_id": ident,
        "pid": proc.pid,
        "status": genome_store.STATUS_DOWNLOADING,
        "estimated_size_label": catalog.get("estimated_size_label"),
        "source": catalog.get("download_url"),
        "destination": genome_store.assembly_dir(ident),
    }
