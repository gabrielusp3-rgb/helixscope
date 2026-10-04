"""Local assembly/reference catalog. Does not expose filesystem paths."""

from __future__ import annotations

from typing import Any

from modules import genome_store, resource_admission
from modules.crispr_assemblies import catalog_by_id
from modules.genome_download import GenomeDownloadError, spawn_download_assembly
from modules.resource_admission import admit_download, admit_reference_search

__all__ = (
    "admit_catalog_download",
    "admit_reference_search",
    "list_references",
    "start_catalog_download",
)


def list_references() -> list[dict[str, Any]]:
    """GRCh38/T2T/GRCm39/TEST REFERENCE metadata with local status.

    Returns:
        Rows with assembly_id, status, ready, kind, label, checksum if known.
        Absolute paths are omitted. TEST REFERENCE is flagged not_a_public_assembly.

    Raises:
        None. Catalog failures yield an empty list.
    """
    rows: list[dict[str, Any]] = []
    try:
        catalog = genome_store.list_local_references()
    except Exception:
        return rows
    for row in catalog:
        if not isinstance(row, dict):
            continue
        rows.append(
            {
                "assembly_id": row.get("id") or row.get("assembly_id"),
                "status": row.get("local_status") or row.get("status"),
                "ready": bool(row.get("ready")),
                "kind": row.get("kind"),
                "label": row.get("label") or row.get("name"),
                "not_a_public_assembly": bool(row.get("not_a_public_assembly")),
                "file_sha256": row.get("file_sha256") or "",
                "n_contigs": int(row.get("n_contigs") or 0),
            }
        )
    return rows


def admit_catalog_download(assembly_id: str) -> dict[str, Any]:
    """Resource admission for a catalogued public assembly. Does not download.

    Args:
        assembly_id: Catalog id such as GRCh38.p14.

    Returns:
        Combined disk download admission plus Cas-OFFinder search admission.
        TEST REFERENCE is not a downloadable public assembly.

    Raises:
        GenomeDownloadError: Unknown id.
    """
    ident = str(assembly_id or "").strip()
    catalog = catalog_by_id(ident)
    if catalog is None:
        raise GenomeDownloadError("Assembly is not in the HelixScope catalog.", "INVALID_INPUT")
    if catalog.get("kind") == "synthetic_test_reference" or catalog.get("not_a_public_assembly"):
        return {
            "status": "TEST_ONLY",
            "assembly_id": ident,
            "decision": "INVALID_INPUT",
            "download": {
                "decision": "INVALID_INPUT",
                "reason": "TEST REFERENCE is not a public genome and is not downloaded.",
            },
            "search": None,
            "not_a_public_assembly": True,
        }
    compressed = int(catalog.get("estimated_compressed_bytes") or 0)
    if compressed <= 0:
        compressed = 1
    download = admit_download(
        assembly_id=ident,
        compressed_bytes=compressed,
        store_path=genome_store.store_root(),
        keep_archive=True,
    )
    search = admit_reference_search(
        assembly_id=ident,
        fasta_bytes=max(compressed * 3, 1),
    )
    decision = str(download.get("decision") or "")
    if str(search.get("decision") or "") == resource_admission.DECISION_RESOURCE_LIMIT:
        decision = resource_admission.DECISION_RESOURCE_LIMIT
    return {
        "status": decision,
        "assembly_id": ident,
        "decision": decision,
        "download": download,
        "search": search,
        "label": catalog.get("label") or catalog.get("name"),
        "kind": catalog.get("kind"),
        "not_a_public_assembly": False,
    }


def start_catalog_download(assembly_id: str) -> dict[str, Any]:
    """Start an explicit catalog download only after admission PROCEED.

    Args:
        assembly_id: Catalog id.

    Returns:
        Job-like dict with admission and spawn metadata. Does not bypass RESOURCE_LIMIT.

    Raises:
        GenomeDownloadError: Invalid id or RESOURCE_LIMIT.
    """
    admission = admit_catalog_download(assembly_id)
    if admission.get("not_a_public_assembly"):
        raise GenomeDownloadError(str((admission.get("download") or {}).get("reason")), "INVALID_INPUT")
    if admission.get("decision") != resource_admission.DECISION_PROCEED:
        raise GenomeDownloadError(
            str((admission.get("download") or {}).get("reason") or "Download is not admitted."),
            "RESOURCE_LIMIT",
        )
    spawned = spawn_download_assembly(str(admission.get("assembly_id") or assembly_id))
    return {
        "status": spawned.get("status") or "DOWNLOADING",
        "assembly_id": admission.get("assembly_id"),
        "admission": admission,
        "spawn": {
            "pid": spawned.get("pid"),
            "status": spawned.get("status"),
        },
        "note": "Download started after resource_admission PROCEED. READY only after checksum and index.",
    }
