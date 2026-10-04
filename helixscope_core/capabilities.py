"""Sanitized engine capability snapshot for API/CLI. No secret values."""

from __future__ import annotations

from typing import Any

from helixscope_core.serialize import redact_filesystem_path
from modules import genome_store, opencl_runtime, tool_registry


def get_capabilities(*, refresh: bool = False) -> dict[str, Any]:
    """Runtime engine/reference capability data with paths redacted.

    Args:
        refresh: Re-run binary detection.

    Returns:
        tools (name -> row without absolute paths), opencl (no WMI dump),
        references summary statuses,         scan note.

    Raises:
        None. Detector failures are recorded in the payload.
    """
    snap = tool_registry.collect_tool_snapshot(refresh=refresh)
    tools_out: dict[str, Any] = {}
    for name, row in (snap.get("tools") or {}).items():
        if not isinstance(row, dict):
            continue
        live = row.get("live")
        live_out = None
        if isinstance(live, dict):
            live_out = {
                "ok": live.get("ok"),
                "status": live.get("status"),
                "version": live.get("version"),
                "engine_location": live.get("engine_location"),
                "validated_at_utc": live.get("validated_at_utc"),
            }
        tools_out[str(name)] = {
            "name": row.get("name"),
            "version": row.get("version"),
            "status": row.get("status"),
            "available": row.get("available"),
            "source": row.get("source"),
            "capabilities": list(row.get("capabilities") or []),
            "reason": row.get("reason"),
            "engine_location": row.get("engine_location"),
            "binary": redact_filesystem_path(row.get("path")),
            "live": live_out,
        }
    try:
        opencl = opencl_runtime.diagnose(refresh=refresh)
    except Exception as exc:
        opencl = {"ready": False, "blocker": type(exc).__name__, "reason": type(exc).__name__}
    if not isinstance(opencl, dict):
        opencl = {"ready": False, "blocker": "invalid diagnose payload"}
    refs = []
    try:
        for row in genome_store.list_local_references():
            if not isinstance(row, dict):
                continue
            refs.append(
                {
                    "assembly_id": row.get("id") or row.get("assembly_id"),
                    "status": row.get("local_status") or row.get("status"),
                    "kind": row.get("kind"),
                    "label": row.get("label") or row.get("name"),
                    "ready": bool(row.get("ready")),
                    "not_a_public_assembly": bool(row.get("not_a_public_assembly")),
                }
            )
    except Exception:
        refs = []
    return {
        "tools": tools_out,
        "opencl": {
            "ready": bool(opencl.get("ready") or opencl.get("opencl_ready_for_cas_offinder")),
            "blocker": str(opencl.get("blocker") or opencl.get("reason") or ""),
            "loader_present": bool(opencl.get("loader_present")),
            "loader_basename": str(opencl.get("loader_basename") or ""),
        },
        "references": refs,
        "platform": snap.get("platform"),
        "python_version": snap.get("python_version"),
        "scan": snap.get("scan"),
        "note": snap.get("note"),
        "remote_providers": remote_provider_inventory(),
    }


def remote_provider_inventory() -> list[dict[str, Any]]:
    """Sanitized remote scientific providers. Does not claim live success.

    Returns:
        Rows with status CONFIGURED (client present) or AUTH_REQUIRED notes.
        Never REMOTE_VALIDATED without a controlled live interaction.
        Secrets are not included.
    """
    from modules.ncbi_fetch import api_key_from_environment

    key_present = bool(api_key_from_environment())
    return [
        {
            "provider": "NCBI E-utilities",
            "service": "Entrez (nucleotide, protein, gene, pubmed, ClinVar)",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "api_key_configured": key_present,
            "limitation": (
                "Contact email is required per request. Official limit is 3 "
                "requests/s without a key and 10 requests/s with a configured "
                "server-side NCBI API key. REMOTE_VALIDATED only after a "
                "controlled live call."
            ),
            "privacy": "Accession, search term, or identifier is sent to NCBI. Sequences are sent only when the user fetches or BLASTs.",
        },
        {
            "provider": "NCBI BLAST URL API",
            "service": "Blast.cgi",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "limitation": "Do not poll a RID more than once per minute. Do not contact more than once per 10 seconds.",
            "privacy": "The BLAST query sequence is sent to NCBI when remote BLAST is chosen.",
        },
        {
            "provider": "UniProt REST",
            "service": "rest.uniprot.org/uniprotkb",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "limitation": "Requires an explicit UniProt accession. Anonymous sequence is not mapped automatically.",
            "privacy": "Only the requested accession is sent, not a pasted workstation sequence.",
        },
        {
            "provider": "InterPro API",
            "service": "www.ebi.ac.uk/interpro/api",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "limitation": "Domains are retrieved annotations, not predicted from sequence appearance.",
            "privacy": "Only the UniProt accession is sent.",
        },
        {
            "provider": "Ensembl VEP REST",
            "service": "rest.ensembl.org",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "limitation": "GRCh38 on rest.ensembl.org; GRCh37 uses the archive host. Not a diagnosis.",
            "privacy": "Variant identity (region/allele) is sent when VEP is enabled.",
        },
        {
            "provider": "EMBL-EBI Job Dispatcher",
            "service": "Clustal Omega REST",
            "local_or_remote": "remote",
            "status": "CONFIGURED",
            "limitation": "Fair-use remote MSA. Engine UNAVAILABLE if the service or contact rules fail.",
            "privacy": "Unaligned sequences are sent to EMBL-EBI when that backend is chosen.",
        },
    ]
