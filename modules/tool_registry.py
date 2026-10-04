"""Registo central de ferramentas cientificas opcionais.

Agrega detectores ja existentes. Nao percorre o disco. Nao executa um caminho
enviado pela UI. Snapshot em memoria do processo para evitar dezenas de which().

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, Optional

from . import (
    blast_search,
    crispr_casoffinder,
    engine_validation,
    phylogeny,
    protein_structure,
    rcsb_alignment,
    rna_folding,
    tool_detection,
    usalign,
)

_SNAPSHOT: Optional[Dict[str, Any]] = None


def _overlay_live(snapshot: Dict[str, Any]) -> Dict[str, Any]:
    """Refresh live/remote validation onto a detector snapshot copy.

    Cached detector rows must not freeze a stale LIVE_VALIDATED /
    REMOTE_VALIDATED badge. Detector availability itself stays as scanned.
    """
    tools = snapshot.get("tools")
    if not isinstance(tools, dict):
        return snapshot
    for name, row in tools.items():
        if not isinstance(row, dict):
            continue
        row["live"] = engine_validation.live_record(str(row.get("name") or name))
    return snapshot


def _tool_row(
    *,
    name: str,
    record: Dict[str, Any],
    source: str,
    capabilities: tuple[str, ...],
) -> dict:
    available = bool(record.get("available"))
    version = str(record.get("version") or "")
    path = str(record.get("path") or "")
    broken = str(record.get("status") or "") == tool_detection.TOOL_STATUS_INVALID
    status = engine_validation.merge_status(
        detected=available,
        version=version,
        tool=name,
        broken=broken,
    )
    if not available:
        status = engine_validation.STATUS_NOT_INSTALLED
    return {
        "name": name,
        "path": path,
        "version": version,
        "status": status,
        "available": available,
        "source": source,
        "capabilities": list(capabilities),
        "reason": str(record.get("reason") or ""),
        "live": engine_validation.live_record(name),
    }


def collect_tool_snapshot(*, refresh: bool = False) -> dict:
    """Detecta ferramentas allowlisted uma vez por processo.

    Args:
        refresh: Se True, ignora o snapshot em memoria e volta a detectar.

    Returns:
        Dict tools (nome -> linha), platform, note. Caminhos ja sanitizados
        pelos detectores.

    Raises:
        Nenhum.
    """
    global _SNAPSHOT
    if _SNAPSHOT is not None and not refresh:
        return _overlay_live(copy.deepcopy(_SNAPSHOT))

    fasttree = phylogeny.detect_fasttree()
    iqtree = phylogeny.detect_iqtree()
    dssp = protein_structure.dssp_availability()
    stride = protein_structure.stride_availability()
    surface = protein_structure.surface_availability()
    vienna = rna_folding.tool_availability()
    casoff = crispr_casoffinder.detect_cas_offinder()
    blast_local = blast_search.local_blast_availability()
    rnafold = vienna.get("rnafold") or rna_folding.detect_rnafold_executable()
    python_rna = vienna.get("python_rna") or {}

    tools = {
        "FastTree": _tool_row(
            name="FastTree",
            record=fasttree,
            source="PATH / HELIXSCOPE_FASTTREE",
            capabilities=("approximate_ml_phylogeny",),
        ),
        "IQ-TREE": _tool_row(
            name="IQ-TREE",
            record=iqtree,
            source="PATH / HELIXSCOPE_IQTREE",
            capabilities=("ml_phylogeny", "modelfinder", "ufboot", "sh_alrt"),
        ),
        "DSSP": _tool_row(
            name="DSSP",
            record=dssp,
            source="PATH / HELIXSCOPE_DSSP",
            capabilities=("secondary_structure_assignment",),
        ),
        "DSSP PDB-REDO API": {
            "name": "DSSP PDB-REDO API",
            "path": "",
            "version": str((protein_structure.dssp_remote_availability() or {}).get("version") or ""),
            "status": str(
                (engine_validation.live_record(protein_structure.DSSP_REMOTE_TOOL) or {}).get("status")
                or "detected"
            ),
            "available": True,
            "source": protein_structure.PDB_REDO_DSSP_DO,
            "capabilities": ["secondary_structure_assignment"],
            "reason": str(protein_structure.dssp_remote_availability().get("reason") or ""),
            "engine_location": "remote",
            "live": engine_validation.live_record(protein_structure.DSSP_REMOTE_TOOL),
        },
        "STRIDE": _tool_row(
            name="STRIDE",
            record=stride,
            source="PATH / HELIXSCOPE_STRIDE",
            capabilities=("secondary_structure_assignment",),
        ),
        "MSMS": _tool_row(
            name="MSMS",
            record=protein_structure.msms_availability(),
            source="PATH / HELIXSCOPE_MSMS / tools/msms/",
            capabilities=("molecular_surface",),
        ),
        "EDTSurf": _tool_row(
            name="EDTSurf",
            record=(
                surface
                if "edtsurf" in str(surface.get("tool") or "").lower()
                else {
                    "available": False,
                    "version": "",
                    "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
                    "reason": str(surface.get("reason") or "EDTSurf is not installed."),
                }
            ),
            source="PATH / HELIXSCOPE_SURFACE / tools/edtsurf/",
            capabilities=("molecular_surface",),
        ),
        "RNAfold": _tool_row(
            name="RNAfold",
            record=rnafold if isinstance(rnafold, dict) else {"available": False},
            source="PATH / HELIXSCOPE_RNAFOLD / tools/viennarna / Program Files ViennaRNA",
            capabilities=("rna_mfe_folding",),
        ),
        "ViennaRNA Python": _tool_row(
            name="ViennaRNA Python",
            record=python_rna if isinstance(python_rna, dict) else {"available": False},
            source="import RNA",
            capabilities=("rna_mfe_folding",),
        ),
        "Cas-OFFinder": _tool_row(
            name="Cas-OFFinder",
            record=casoff,
            source="PATH / HELIXSCOPE_CAS_OFFINDER",
            capabilities=("off_target_search",),
        ),
        "BLAST+": {
            "name": "BLAST+",
            "path": str((blast_local.get("programs") or {}).get("blastn", {}).get("path") or ""),
            "version": str(blast_local.get("version") or ""),
            "status": engine_validation.merge_status(
                detected=bool(blast_local.get("executables_detected")),
                version=str(blast_local.get("version") or ""),
                tool="BLAST+",
            ),
            "available": bool(blast_local.get("executables_detected")),
            "source": "PATH / HELIXSCOPE_BLASTN / tools/",
            "capabilities": ["local_blast_search"],
            "reason": str(blast_local.get("reason") or ""),
            "database_available": bool((blast_local.get("database") or {}).get("available")),
            "live": engine_validation.live_record("BLAST+"),
        },
        "RCSB Alignment API": {
            "name": "RCSB Alignment API",
            "path": "",
            "version": "api/v1",
            "status": (
                str(
                    (engine_validation.live_record("RCSB Alignment API") or {}).get("status")
                    or ""
                )
                or "detected"
            ),
            "available": True,
            "source": rcsb_alignment.DOCS_URL,
            "capabilities": ["protein_structure_alignment"],
            "reason": (
                "Remote RCSB PDB Structure Alignment API. Not a local binary. "
                "Success is remote_validated, never live_validated."
            ),
            "engine_location": "remote",
            "methods": [row["api_name"] for row in rcsb_alignment.supported_methods()],
            "live": engine_validation.live_record("RCSB Alignment API"),
        },
        "US-align": _tool_row(
            name="US-align",
            record=usalign.detect_usalign(),
            source="PATH / HELIXSCOPE_USALIGN / https://zhanggroup.org/US-align/",
            capabilities=("protein_na_complex_structure_alignment",),
        ),
    }
    platform = tool_detection.environment_report(tools={})
    snapshot = {
        "tools": tools,
        "platform": platform.get("platform"),
        "python_version": platform.get("python_version"),
        "architecture": platform.get("architecture"),
        "scan": "PATH, explicit env files, known Program Files files. No recursive disk walk.",
        "note": (
            "Detected is not live_validated. Live_validated requires a real run "
            "with accepted output. Optional binaries never run at import."
        ),
    }
    _SNAPSHOT = snapshot
    return _overlay_live(copy.deepcopy(_SNAPSHOT))


def tool_row(name: str, *, refresh: bool = False) -> dict:
    """Uma linha do registo, ou UNAVAILABLE se o nome nao existir.

    Args:
        name: Chave em collect_tool_snapshot()['tools'].
        refresh: Re-detectar.

    Returns:
        Linha do tool ou available False.

    Raises:
        Nenhum.
    """
    snap = collect_tool_snapshot(refresh=refresh)
    row = (snap.get("tools") or {}).get(name)
    if isinstance(row, dict):
        return dict(row)
    return {
        "name": name,
        "path": "",
        "version": "",
        "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
        "available": False,
        "source": "",
        "capabilities": [],
        "reason": "Unknown tool name in the HelixScope registry.",
    }
