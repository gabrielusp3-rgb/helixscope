"""Armazenamento local de referencias genomicas verificadas.

Estados: NOT_INSTALLED, DOWNLOADING, VERIFYING, READY, CORRUPT, PARTIAL, ERROR.
READY exige checksum/content validation. Ficheiros parciais nunca sao READY.
Nenhuma funcao importa Streamlit. Genomas gigantes nao sao commitados.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from typing import Any, Optional

from . import crispr_assemblies, genome_fasta, provenance

STATUS_NOT_INSTALLED: str = "NOT_INSTALLED"
STATUS_DOWNLOADING: str = "DOWNLOADING"
STATUS_VERIFYING: str = "VERIFYING"
STATUS_READY: str = "READY"
STATUS_CORRUPT: str = "CORRUPT"
STATUS_PARTIAL: str = "PARTIAL"
STATUS_ERROR: str = "ERROR"

STORE_ENV: str = "HELIXSCOPE_REFERENCE_DIR"
JOBS_ENV: str = "HELIXSCOPE_JOBS_DIR"
MANIFEST_NAME: str = "manifest.json"
FASTA_NAME: str = "genome.fa"
FAI_NAME: str = "genome.fa.fai"
PARTIAL_SUFFIX: str = ".partial"
_ASSEMBLY_ID_RE: re.Pattern[str] = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}\Z")

MIN_FREE_BYTES: int = 200 * 1024 * 1024


class GenomeStoreError(ValueError):
    """Falha do repositorio local.

    Attributes:
        category: INVALID_INPUT, RESOURCE_LIMIT, CORRUPT, ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def repository_root() -> str:
    """Raiz do repositorio HelixScope.

    Args:
        Nenhum.

    Returns:
        Caminho absoluto.

    Raises:
        Nenhum.
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def store_root() -> str:
    """Diretorio de referencias: HELIXSCOPE_REFERENCE_DIR ou data/references.

    Args:
        Nenhum.

    Returns:
        Caminho absoluto (pode ainda nao existir).

    Raises:
        Nenhum.
    """
    env = (os.environ.get(STORE_ENV) or "").strip().strip('"')
    if env:
        return os.path.abspath(env)
    return os.path.join(repository_root(), "data", "references")


def jobs_root() -> str:
    """Diretorio de jobs locais: HELIXSCOPE_JOBS_DIR ou data/jobs.

    Args:
        Nenhum.

    Returns:
        Caminho absoluto.

    Raises:
        Nenhum.
    """
    env = (os.environ.get(JOBS_ENV) or "").strip().strip('"')
    if env:
        return os.path.abspath(env)
    return os.path.join(repository_root(), "data", "jobs")


def assembly_dir(assembly_id: str) -> str:
    """Pasta de uma assembly no store.

    Args:
        assembly_id: id de catalogo.

    Returns:
        Caminho absoluto.

    Raises:
        GenomeStoreError: INVALID_INPUT.
    """
    ident = _safe_id(assembly_id)
    root = os.path.realpath(store_root())
    path = os.path.realpath(os.path.join(root, ident))
    try:
        shared = os.path.commonpath([root, path])
    except ValueError as exc:
        raise GenomeStoreError("Assembly path leaves the reference store.", "INVALID_INPUT") from exc
    if os.path.normcase(shared) != os.path.normcase(root):
        raise GenomeStoreError("Assembly path leaves the reference store.", "INVALID_INPUT")
    return path


def _safe_id(assembly_id: str) -> str:
    """Map a request to one catalog id. The returned text is the catalog constant.

    Args:
        assembly_id: Caller-supplied label.

    Returns:
        The catalog id, not the raw caller string.

    Raises:
        GenomeStoreError: INVALID_INPUT when the label is not in the catalog.
    """
    ident = str(assembly_id or "")
    if ident == "" or ident != ident.strip() or ".." in ident or "\x00" in ident:
        raise GenomeStoreError("Assembly id is not a safe directory name.", "INVALID_INPUT")
    if _ASSEMBLY_ID_RE.fullmatch(ident) is None:
        raise GenomeStoreError("Assembly id contains illegal path characters.", "INVALID_INPUT")
    if ident == "GRCh38.p14":
        return "GRCh38.p14"
    if ident == "T2T-CHM13v2.0":
        return "T2T-CHM13v2.0"
    if ident == "GRCm39":
        return "GRCm39"
    if ident == crispr_assemblies.TEST_REFERENCE_ID:
        return crispr_assemblies.TEST_REFERENCE_ID
    raise GenomeStoreError("Assembly id is not in the public catalog.", "INVALID_INPUT")


def path_is_in_store(path: str) -> bool:
    """True se o caminho absoluto estiver sob o store ou fixtures de teste.

    Args:
        path: Caminho de ficheiro.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    abs_path = os.path.abspath(path)
    roots = [os.path.abspath(store_root()), os.path.join(repository_root(), "tests", "fixtures")]
    for root in roots:
        try:
            common = os.path.commonpath([abs_path, root])
        except ValueError:
            continue
        if common == root:
            return True
    return False


def disk_status(*, need_bytes: int = 0) -> dict:
    """Espaco livre no volume do store.

    Args:
        need_bytes: Bytes adicionais pretendidos.

    Returns:
        Dict free_bytes, total_bytes, ok.

    Raises:
        Nenhum.
    """
    root = store_root()
    parent = root if os.path.isdir(root) else os.path.dirname(root) or root
    if not os.path.isdir(parent):
        parent = repository_root()
    usage = shutil.disk_usage(parent)
    required = max(int(need_bytes), 0) + MIN_FREE_BYTES
    return {
        "free_bytes": int(usage.free),
        "total_bytes": int(usage.total),
        "need_bytes": int(need_bytes),
        "ok": int(usage.free) >= required,
        "reason": (
            ""
            if int(usage.free) >= required
            else "Not enough free disk space for this reference operation."
        ),
    }


def empty_manifest(assembly_id: str) -> dict:
    """Manifesto NOT_INSTALLED.

    Args:
        assembly_id: id.

    Returns:
        Dict manifesto.

    Raises:
        Nenhum.
    """
    catalog = crispr_assemblies.catalog_by_id(assembly_id) or {}
    return {
        "assembly_id": str(assembly_id or ""),
        "accession": str(catalog.get("accession") or ""),
        "assembly": str(catalog.get("assembly") or assembly_id),
        "organism": str(catalog.get("organism") or ""),
        "source": str(catalog.get("source") or ""),
        "version": str(catalog.get("release") or ""),
        "kind": str(catalog.get("kind") or "public_assembly"),
        "file": "",
        "fai": "",
        "checksum": "",
        "file_sha256": "",
        "official_compressed_md5": str(catalog.get("official_compressed_md5") or ""),
        "compressed_md5_observed": "",
        "downloaded_at": "",
        "verified_at": "",
        "status": STATUS_NOT_INSTALLED,
        "n_contigs": 0,
        "not_a_public_assembly": bool(catalog.get("not_a_public_assembly")),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "catalog_snapshot": crispr_assemblies.CATALOG_SNAPSHOT,
    }


def read_manifest(assembly_id: str) -> dict:
    """Le o manifesto local ou devolve NOT_INSTALLED.

    Args:
        assembly_id: id.

    Returns:
        Dict manifesto.

    Raises:
        Nenhum.
    """
    path = os.path.join(assembly_dir(assembly_id), MANIFEST_NAME)
    if not os.path.isfile(path):
        return empty_manifest(assembly_id)
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        row = empty_manifest(assembly_id)
        row["status"] = STATUS_ERROR
        row["reason"] = "Manifest JSON is unreadable."
        return row
    if not isinstance(data, dict):
        row = empty_manifest(assembly_id)
        row["status"] = STATUS_ERROR
        row["reason"] = "Manifest JSON is not an object."
        return row
    return data


def write_manifest(assembly_id: str, manifest: dict) -> dict:
    """Grava manifesto.

    Args:
        assembly_id: id.
        manifest: Dict.

    Returns:
        Copia gravada.

    Raises:
        GenomeStoreError: ERROR.
    """
    ident = _safe_id(assembly_id)
    directory = assembly_dir(ident)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, MANIFEST_NAME)
    payload = dict(manifest)
    payload["assembly_id"] = ident
    payload["software_version"] = provenance.HELIXSCOPE_VERSION
    try:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
    except OSError as exc:
        raise GenomeStoreError(f"Could not write manifest: {exc}", "ERROR") from exc
    return payload


def list_local_references() -> list[dict]:
    """Catalogo publico + estado local, mais o FASTA de teste se instalado.

    Args:
        Nenhum.

    Returns:
        Lista de dicts para a UI.

    Raises:
        Nenhum.
    """
    rows = []
    for item in crispr_assemblies.list_catalog(include_test=True):
        ident = str(item.get("id") or "")
        local = read_manifest(ident)
        row = dict(item)
        row["local_status"] = str(local.get("status") or STATUS_NOT_INSTALLED)
        row["ready"] = str(local.get("status") or "") == STATUS_READY
        row["file_sha256"] = str(local.get("file_sha256") or "")
        row["n_contigs"] = int(local.get("n_contigs") or 0)
        row["verified_at"] = str(local.get("verified_at") or "")
        row["path_public"] = ident if row["ready"] else ""
        rows.append(row)
    return rows


def ready_record(assembly_id: str) -> Optional[dict]:
    """Registo READY com caminhos, ou None.

    Args:
        assembly_id: id.

    Returns:
        Dict fasta_path, fai_path, manifest, ou None.

    Raises:
        Nenhum.
    """
    manifest = read_manifest(assembly_id)
    if str(manifest.get("status") or "") != STATUS_READY:
        return None
    fasta = str(manifest.get("file") or "")
    fai = str(manifest.get("fai") or "")
    if not fasta or not os.path.isfile(fasta) or not os.path.isfile(fai):
        return None
    if not path_is_in_store(fasta):
        return None
    expected = str(manifest.get("file_sha256") or "")
    if expected and not genome_fasta.fai_is_current(fasta, fai, expected):
        return None
    return {
        "assembly_id": assembly_id,
        "fasta_path": fasta,
        "fai_path": fai,
        "manifest": manifest,
        "kind": str(manifest.get("kind") or ""),
        "not_a_public_assembly": bool(manifest.get("not_a_public_assembly")),
    }


def any_public_assembly_ready() -> bool:
    """True se alguma assembly publica do catalogo estiver READY e verificada.

    Args:
        Nenhum.

    Returns:
        bool. O FASTA de teste nunca conta.

    Raises:
        Nenhum.
    """
    for item in crispr_assemblies.ASSEMBLY_CATALOG:
        record = ready_record(str(item.get("id") or ""))
        if record is not None and not record.get("not_a_public_assembly"):
            return True
    return False


def mark_partial(assembly_id: str, *, reason: str) -> dict:
    """Marca PARTIAL. Nunca READY.

    Args:
        assembly_id: id.
        reason: Motivo.

    Returns:
        Manifesto.

    Raises:
        Nenhum.
    """
    manifest = read_manifest(assembly_id)
    manifest["status"] = STATUS_PARTIAL
    manifest["reason"] = str(reason or "Partial file.")
    manifest["verified_at"] = ""
    return write_manifest(assembly_id, manifest)


def mark_corrupt(assembly_id: str, *, reason: str) -> dict:
    """Marca CORRUPT apos checksum falhar.

    Args:
        assembly_id: id.
        reason: Motivo.

    Returns:
        Manifesto.

    Raises:
        Nenhum.
    """
    manifest = read_manifest(assembly_id)
    manifest["status"] = STATUS_CORRUPT
    manifest["reason"] = str(reason or "Checksum mismatch.")
    manifest["verified_at"] = ""
    return write_manifest(assembly_id, manifest)


def install_test_reference(*, fixture_path: str = "") -> dict:
    """Copia o FASTA de teste para o store, indexa e marca READY.

    Args:
        fixture_path: FASTA de teste; default tests/fixtures.

    Returns:
        Manifesto READY.

    Raises:
        GenomeStoreError: INVALID_INPUT.
    """
    ident = crispr_assemblies.TEST_REFERENCE_ID
    source = fixture_path or os.path.join(
        repository_root(), "tests", "fixtures", "helixscope_test_reference.fa"
    )
    if not os.path.isfile(source):
        raise GenomeStoreError("Test reference FASTA fixture is missing.", "INVALID_INPUT")
    directory = assembly_dir(ident)
    os.makedirs(directory, exist_ok=True)
    dest = os.path.join(directory, FASTA_NAME)
    shutil.copyfile(source, dest)
    return finalize_ready_from_fasta(
        ident,
        dest,
        kind="synthetic_test_reference",
        not_a_public_assembly=True,
        official_compressed_md5="",
        compressed_md5_observed="",
    )


def finalize_ready_from_fasta(
    assembly_id: str,
    fasta_path: str,
    *,
    kind: str,
    not_a_public_assembly: bool,
    official_compressed_md5: str,
    compressed_md5_observed: str,
) -> dict:
    """Indexa FASTA, grava hashes e so entao marca READY.

    Args:
        assembly_id: id.
        fasta_path: FASTA ja completo no store.
        kind: public_assembly ou synthetic_test_reference.
        not_a_public_assembly: True para o FASTA de teste.
        official_compressed_md5: MD5 NCBI do .gz, ou vazio.
        compressed_md5_observed: MD5 observado do .gz, ou vazio.

    Returns:
        Manifesto READY.

    Raises:
        GenomeStoreError: CORRUPT / ERROR.
    """
    ident = _safe_id(assembly_id)
    if not path_is_in_store(fasta_path):
        raise GenomeStoreError("FASTA is outside the reference store.", "INVALID_INPUT")
    if kind == "public_assembly" and (
        not str(official_compressed_md5 or "").strip()
        or not str(compressed_md5_observed or "").strip()
    ):
        raise GenomeStoreError(
            "Public assembly READY requires the NCBI compressed MD5 and the "
            "observed MD5. Empty checksums are refused.",
            "CORRUPT",
        )
    catalog = crispr_assemblies.catalog_by_id(ident) or {}
    write_manifest(
        ident,
        {
            **empty_manifest(ident),
            "status": STATUS_VERIFYING,
            "file": os.path.abspath(fasta_path),
        },
    )
    try:
        indexed = genome_fasta.build_fai_for_fasta(fasta_path)
    except genome_fasta.GenomeFastaError as exc:
        mark_corrupt(ident, reason=str(exc))
        raise GenomeStoreError(str(exc), "CORRUPT") from exc
    if official_compressed_md5 and compressed_md5_observed:
        if compressed_md5_observed.lower() != official_compressed_md5.lower():
            mark_corrupt(
                ident,
                reason=(
                    "Observed MD5 of genomic.fna.gz does not match NCBI "
                    "md5checksums.txt. File is not READY."
                ),
            )
            raise GenomeStoreError("Official compressed MD5 mismatch.", "CORRUPT")
    manifest = {
        **empty_manifest(ident),
        "status": STATUS_READY,
        "file": os.path.abspath(fasta_path),
        "fai": indexed["fai_path"],
        "file_sha256": indexed["file_sha256"],
        "checksum": indexed["file_sha256"],
        "official_compressed_md5": official_compressed_md5
        or str(catalog.get("official_compressed_md5") or ""),
        "compressed_md5_observed": compressed_md5_observed,
        "n_contigs": indexed["n_contigs"],
        "contig_manifest": [
            {
                "identifier": str(row.get("name") or ""),
                "length": int(row.get("length") or 0),
            }
            for row in list(indexed.get("rows") or [])
        ],
        "kind": kind,
        "not_a_public_assembly": bool(not_a_public_assembly),
        "downloaded_at": provenance.utc_now(),
        "verified_at": provenance.utc_now(),
        "reason": "",
    }
    return write_manifest(ident, manifest)
