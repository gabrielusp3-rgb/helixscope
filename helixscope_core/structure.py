"""Structure science. Coordinates stay NumPy/lists from modules; no Plotly camera."""

from __future__ import annotations

import json
from typing import Any

from helixscope_core.errors import StructureError
from modules.crispr_structure_catalog import (
    list_experimental_complexes,
    mapping_contract,
)
from modules.protein_structure import (
    BUNDLED_MMCIF_DIR,
    classify_mapping_status,
    dssp_availability,
    dssp_remote_availability,
    fetch_alphafold_prediction,
    generate_molecular_surface,
    load_deposited_macromolecule,
    load_protein_structure,
    parse_rcsb_entry_metadata,
    read_bundled_mmcif,
    search_rcsb_by_sequence,
    search_structure_sources,
    shrake_rupley_sasa,
    stride_availability,
    structure_3d_input,
    surface_availability,
    assign_secondary_structure_dssp,
    assign_secondary_structure_dssp_remote,
    assign_secondary_structure_stride,
)

BUNDLED_STRUCTURE_IDS: frozenset[str] = frozenset({"1CRN", "1BNA", "1RNA"})

__all__ = (
    "BUNDLED_MMCIF_DIR",
    "BUNDLED_STRUCTURE_IDS",
    "analyze_structure_coordinates",
    "bundled_coordinate_file",
    "classify_mapping_status",
    "fetch_alphafold_prediction",
    "list_experimental_complexes",
    "load_bundled_deposited",
    "load_bundled_mapped",
    "load_catalog_complex",
    "load_deposited_macromolecule",
    "load_protein_structure",
    "mapping_contract",
    "read_bundled_mmcif",
    "search_rcsb_by_sequence",
    "search_structure_sources",
    "structure_3d_input",
    "structure_engine_capabilities",
    "viewer_payload",
)


def _bundled_metadata(pdb_id: str) -> dict[str, Any]:
    """Local RCSB entry metadata for a bundled fixture. Never fetches."""
    fallback = {"kind": "experimental", "source": "RCSB PDB", "structure_id": pdb_id}
    meta_path = BUNDLED_MMCIF_DIR / f"rcsb_entry_{pdb_id}.json"
    if not meta_path.is_file():
        return fallback
    try:
        raw = json.loads(meta_path.read_text(encoding="utf-8"))
        parsed = parse_rcsb_entry_metadata(raw)
    except Exception:
        return fallback
    if not isinstance(parsed, dict) or not parsed:
        return fallback
    parsed.setdefault("kind", "experimental")
    parsed.setdefault("source", "RCSB PDB")
    return parsed


def load_bundled_deposited(structure_id: str) -> dict[str, Any]:
    """Load a bundled RCSB mmCIF copy (1CRN/1BNA/1RNA). No network fetch.

    Args:
        structure_id: PDB ID of a bundled fixture.

    Returns:
        Deposited macromolecule dict from ``load_deposited_macromolecule``.

    Raises:
        StructureError: Unknown id or missing fixture.
    """
    pdb_id = str(structure_id or "").strip().upper()
    if pdb_id not in BUNDLED_STRUCTURE_IDS:
        raise StructureError(
            "Bundled structure endpoint only serves 1CRN, 1BNA, or 1RNA.",
            "INVALID_INPUT",
        )
    text = read_bundled_mmcif(f"{pdb_id}.cif")
    return load_deposited_macromolecule(
        source="RCSB PDB",
        structure_id=pdb_id,
        structure_text=text,
        metadata=_bundled_metadata(pdb_id),
    )


def load_bundled_mapped(
    structure_id: str,
    query_sequence: str,
    *,
    molecule: str = "PROTEIN",
) -> dict[str, Any]:
    """Map a query onto a bundled deposited structure. No network fetch.

    Args:
        structure_id: PDB ID of a bundled fixture.
        query_sequence: Polymer sequence to map.
        molecule: PROTEIN, DNA, or RNA.

    Returns:
        Mapped structure dict from ``load_protein_structure``.

    Raises:
        StructureError: Invalid id, molecule, or mapping input.
    """
    pdb_id = str(structure_id or "").strip().upper()
    if pdb_id not in BUNDLED_STRUCTURE_IDS:
        raise StructureError(
            "Bundled structure endpoint only serves 1CRN, 1BNA, or 1RNA.",
            "INVALID_INPUT",
        )
    text = read_bundled_mmcif(f"{pdb_id}.cif")
    return load_protein_structure(
        query_sequence,
        source="RCSB PDB",
        structure_id=pdb_id,
        structure_text=text,
        metadata=_bundled_metadata(pdb_id),
        molecule=molecule,
    )


def load_catalog_complex(structure_id: str) -> dict[str, Any]:
    """Load a catalogued Cas-guide-DNA complex from RCSB Files API.

    Args:
        structure_id: PDB identifier that must already appear in
            ``list_experimental_complexes`` (currently 4UN3 or 4OO8).

    Returns:
        Experimental deposited macromolecule plus a viewer payload. Mapping of
        a user guide onto this structure remains UNCERTAIN unless a later
        classify call proves an exact polymer match. Coordinates are RCSB
        experimental, never a synthetic Cas complex.

    Raises:
        StructureError: Unknown catalog id or fetch/parse failure.
    """
    pdb_id = str(structure_id or "").strip().upper()
    catalog = {
        str(item.get("pdb_id") or "").upper(): item
        for item in list_experimental_complexes()
    }
    if pdb_id not in catalog:
        raise StructureError(
            "Only catalogued Cas-guide-DNA PDB identifiers can be fetched. "
            "Unknown ids are refused so this endpoint cannot become an open SSRF.",
            "INVALID_INPUT",
        )
    parsed = load_deposited_macromolecule(source="RCSB PDB", structure_id=pdb_id)
    scene = viewer_payload(parsed)
    for banned in ("camera", "eye", "zoom", "mediapipe", "gesture", "rotation"):
        scene.pop(banned, None)
    return {
        "status": "EXPERIMENTAL",
        "structure_id": pdb_id,
        "kind": parsed.get("kind") or "experimental",
        "source": parsed.get("source") or f"RCSB PDB {pdb_id}",
        "mapping_status": "UNCERTAIN",
        "mapping_note": (
            "Guide and PAM highlights require an exact subsequence match on the "
            "deposited polymer. A protospacer span on pasted DNA is a sequence "
            "interval, not a 3D correspondence."
        ),
        "catalog": catalog[pdb_id],
        "structure_hash": parsed.get("structure_hash") or parsed.get("content_hash"),
        "sequence_hash": parsed.get("sequence_hash"),
        "n_atoms": scene.get("n_atoms") or len(scene.get("coordinates") or scene.get("atoms") or []),
        "viewer": scene,
        "disclaimer": catalog[pdb_id].get("notes"),
    }


def bundled_coordinate_file(structure_id: str) -> dict[str, Any]:
    """Return allowlisted bundled mmCIF text for the molecular viewer.

    Args:
        structure_id: Bundled PDB id (1CRN, 1BNA, or 1RNA).

    Returns:
        Experimental mmCIF transport. Not a network fetch.

    Raises:
        StructureError: Unknown or non-bundled identifier.
    """
    pdb_id = str(structure_id or "").strip().upper()
    if pdb_id not in BUNDLED_STRUCTURE_IDS:
        raise StructureError(
            "Coordinate-file endpoint only serves bundled 1CRN, 1BNA, or 1RNA. "
            "Arbitrary PDB ids and URLs are refused.",
            "INVALID_INPUT",
        )
    text = read_bundled_mmcif(f"{pdb_id}.cif")
    return {
        "status": "EXPERIMENTAL",
        "structure_id": pdb_id,
        "kind": "experimental",
        "source": "bundled_fixture",
        "format": "mmcif",
        "mmcif": text,
        "n_chars": len(text),
        "limitation": "Bundled RCSB mmCIF copy. Not a user URL. Not AlphaFold.",
    }


def viewer_payload(result: dict[str, Any]) -> dict[str, Any]:
    """3D transport envelope. No camera or gesture fields.

    Args:
        result: Structure dict from load helpers.

    Returns:
        ``structure_3d_input`` scene with atoms copied from chains if needed.

    Raises:
        None.
    """
    adapted = dict(result)
    atoms = list(adapted.get("atoms") or [])
    if not atoms:
        for chain in list(adapted.get("chains") or []):
            if isinstance(chain, dict):
                atoms.extend(list(chain.get("atoms") or []))
        adapted["atoms"] = atoms
    if not adapted.get("selected_chain"):
        chains = list(adapted.get("chains") or [])
        if chains and isinstance(chains[0], dict):
            adapted["selected_chain"] = chains[0]
    scene = structure_3d_input(adapted)
    for banned in ("camera", "eye", "zoom", "mediapipe", "gesture", "rotation"):
        scene.pop(banned, None)
    return scene


def structure_engine_capabilities() -> dict[str, Any]:
    """Detector snapshot for protein-structure engines. No invented availability.

    Returns:
        Named engine rows with available/status/version/limitation.

    Raises:
        None. Detector failures stay on the row.
    """
    from modules import engine_validation

    def _row(name: str, record: dict[str, Any], *, limitation: str) -> dict[str, Any]:
        live = engine_validation.live_record(name) or {}
        available = bool(record.get("available"))
        return {
            "name": name,
            "available": available,
            "status": record.get("status") or ("NOT_INSTALLED" if not available else "detected"),
            "version": record.get("version") or live.get("version") or "",
            "engine_location": record.get("engine_location") or live.get("engine_location") or "local",
            "reason": record.get("reason") or "",
            "limitation": limitation,
            "live": {
                "ok": live.get("ok"),
                "status": live.get("status"),
                "validated_at_utc": live.get("validated_at_utc"),
            },
        }

    dssp = dssp_availability() if callable(dssp_availability) else {}
    dssp_remote = dssp_remote_availability()
    stride = stride_availability()
    surface = surface_availability()
    return {
        "status": "COMPUTED",
        "dssp_local": _row(
            "DSSP",
            dssp if isinstance(dssp, dict) else {},
            limitation="Geometric assignment of deposited coordinates, not sequence SS prediction.",
        ),
        "dssp_pdb_redo": {
            "name": "DSSP PDB-REDO API",
            "available": True,
            "status": "CONFIGURED",
            "version": str((dssp_remote or {}).get("version") or ""),
            "engine_location": "remote",
            "reason": str((dssp_remote or {}).get("reason") or ""),
            "limitation": "Sends coordinates to pdb-redo.eu. Not a local mkdssp run.",
            "live": engine_validation.live_record("DSSP PDB-REDO API"),
        },
        "stride": _row(
            "STRIDE",
            stride if isinstance(stride, dict) else {},
            limitation="STRIDE is separate from DSSP. Results are not mixed.",
        ),
        "edtsurf": _row(
            "EDTSurf",
            surface if isinstance(surface, dict) else {},
            limitation="Triangulated molecular surface. Not Shrake-Rupley SASA. Not MSMS.",
        ),
        "sasa": {
            "name": "Shrake-Rupley SASA",
            "available": True,
            "status": "COMPUTED",
            "version": "",
            "engine_location": "local",
            "reason": "",
            "limitation": "Geometric SASA on ATOM records. Not a molecular surface mesh.",
            "live": None,
        },
        "alphafold_retrieval": {
            "name": "AlphaFold DB",
            "available": True,
            "status": "CONFIGURED",
            "engine_location": "remote",
            "limitation": "Retrieves deposited AFDB predictions. Does not run AlphaFold locally.",
        },
        "rcsb_search": {
            "name": "RCSB Search API",
            "available": True,
            "status": "CONFIGURED",
            "engine_location": "remote",
            "limitation": "Official sequence search (MMseqs2). Not an invented structural score.",
        },
    }


def analyze_structure_coordinates(
    parsed: dict[str, Any],
    *,
    run_sasa: bool = False,
    run_dssp: bool = False,
    run_dssp_remote: bool = False,
    run_stride: bool = False,
    run_edtsurf: bool = False,
    structure_text: str = "",
) -> dict[str, Any]:
    """Run requested coordinate analyses. Missing engines stay UNAVAILABLE.

    Args:
        parsed: Structure dict with atoms.
        run_sasa: Shrake-Rupley.
        run_dssp: Local mkdssp when installed.
        run_dssp_remote: PDB-REDO DSSP HTTP.
        run_stride: Local STRIDE when installed.
        run_edtsurf: Local EDTSurf when installed.
        structure_text: mmCIF/PDB text required by DSSP/STRIDE/EDTSurf.

    Returns:
        Named analysis envelopes. No sequence SS invention.

    Raises:
        StructureError: Invalid parsed input when an analysis is requested.
    """
    payload: dict[str, Any] = {
        "status": "COMPUTED",
        "kind": parsed.get("kind"),
        "source": parsed.get("source"),
        "structure_id": parsed.get("structure_id") or parsed.get("id"),
        "sasa": {"status": "NOT_COMPUTED"},
        "dssp": {"status": "NOT_COMPUTED"},
        "dssp_remote": {"status": "NOT_COMPUTED"},
        "stride": {"status": "NOT_COMPUTED"},
        "edtsurf": {"status": "NOT_COMPUTED"},
    }
    if run_sasa:
        payload["sasa"] = shrake_rupley_sasa(parsed)
    if run_dssp:
        payload["dssp"] = assign_secondary_structure_dssp(structure_text)
    if run_dssp_remote:
        payload["dssp_remote"] = assign_secondary_structure_dssp_remote(structure_text)
    if run_stride:
        payload["stride"] = assign_secondary_structure_stride(structure_text)
    if run_edtsurf:
        payload["edtsurf"] = generate_molecular_surface(
            kind=str(parsed.get("kind") or "experimental"),
            parsed=parsed,
            structure_hash=str(parsed.get("structure_hash") or parsed.get("content_hash") or ""),
        )
    return payload
