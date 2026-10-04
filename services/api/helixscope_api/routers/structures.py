"""Structure scientific data. No Plotly figures or camera state."""

from __future__ import annotations

from fastapi import APIRouter, Request

from helixscope_core.structure import (
    bundled_coordinate_file,
    classify_mapping_status,
    list_experimental_complexes,
    load_bundled_deposited,
    load_bundled_mapped,
    load_catalog_complex,
    mapping_contract,
    viewer_payload,
)

from helixscope_api.dependencies import envelope
from helixscope_api.schemas import (
    ApiEnvelope,
    StructureComplexRequest,
    StructureCoordinateFileRequest,
    StructureFixtureRequest,
    StructureMappedRequest,
    StructureMappingRequest,
)

router = APIRouter(prefix="/api/v1/structures", tags=["structures"])


@router.get(
    "/catalog/experimental-complexes",
    operation_id="structureExperimentalCatalog",
    summary="Published Cas-guide-DNA PDB pointers without fabricated coordinates",
    response_model=ApiEnvelope,
)
def experimental_catalog(request: Request) -> ApiEnvelope:
    """4UN3 remains a pointer. Mapping is not promoted."""
    rows = list_experimental_complexes()
    contract = mapping_contract()
    return envelope(
        request,
        {
            "complexes": rows,
            "mapping_contract": contract,
            "status": contract.get("status") or "UNAVAILABLE",
        },
    )


@router.post(
    "/fixture",
    operation_id="structureLoadFixture",
    summary="Load bundled 1CRN/1BNA/1RNA mmCIF (no network)",
    response_model=ApiEnvelope,
)
def structure_fixture(request: Request, body: StructureFixtureRequest) -> ApiEnvelope:
    """Coordinates from the deposited copy. No camera fields."""
    parsed = load_bundled_deposited(body.structure_id)
    scene = viewer_payload(parsed)
    kind = str(parsed.get("kind") or scene.get("kind") or "experimental")
    status = {"experimental": "EXPERIMENTAL", "predicted": "PREDICTED", "illustrative": "ILLUSTRATIVE"}.get(
        kind.lower(),
        "COMPUTED",
    )
    payload = {
        "status": status,
        "structure_id": body.structure_id,
        "kind": kind,
        "source": parsed.get("source"),
        "structure_hash": parsed.get("structure_hash") or parsed.get("content_hash"),
        "sequence_hash": parsed.get("sequence_hash"),
        "n_atoms": scene.get("n_atoms") or len(scene.get("coordinates") or []),
        "viewer": scene,
    }
    return envelope(request, payload, drop_heavy=True)


@router.post(
    "/fixture/mapped",
    operation_id="structureLoadFixtureMapped",
    summary="Map a query onto a bundled deposited structure",
    response_model=ApiEnvelope,
)
def structure_fixture_mapped(request: Request, body: StructureMappedRequest) -> ApiEnvelope:
    """Preserves mapping_status including UNCERTAIN. No promotion."""
    mapped = load_bundled_mapped(body.structure_id, body.query_sequence, molecule=body.molecule)
    scene = viewer_payload(mapped)
    payload = {
        "status": mapped.get("status") or mapped.get("mapping_status") or "COMPUTED",
        "structure_id": body.structure_id,
        "mapping_status": mapped.get("mapping_status"),
        "mapping_status_note": mapped.get("mapping_status_note"),
        "kind": mapped.get("kind"),
        "viewer": scene,
        "sequence_hash": mapped.get("sequence_hash"),
        "structure_hash": mapped.get("content_hash") or mapped.get("structure_hash"),
    }
    return envelope(request, payload, drop_heavy=True)


@router.post(
    "/mapping/classify",
    operation_id="structureClassifyMapping",
    summary="Classify mapping status without promoting best-effort mapping",
    response_model=ApiEnvelope,
)
def structure_classify_mapping(request: Request, body: StructureMappingRequest) -> ApiEnvelope:
    """polymer_length 0 stays UNCERTAIN. Used for 4UN3-class uncertainty."""
    payload = classify_mapping_status(
        sequences_identical=body.sequences_identical,
        polymer_index_mode=body.polymer_index_mode,
        observed_coordinate_residues=body.observed_coordinate_residues,
        polymer_length=body.polymer_length,
    )
    return envelope(request, payload, scientific_status=str(payload.get("mapping_status") or ""))


@router.post(
    "/experimental-complex",
    operation_id="structureLoadExperimentalComplex",
    summary="Fetch a catalogued Cas-guide-DNA complex from RCSB",
    response_model=ApiEnvelope,
)
def structure_load_experimental_complex(request: Request, body: StructureComplexRequest) -> ApiEnvelope:
    """4UN3/4OO8 experimental coordinates. Guide mapping stays UNCERTAIN."""
    payload = load_catalog_complex(body.structure_id)
    status = str(payload.get("status") or "EXPERIMENTAL")
    return envelope(request, payload, scientific_status=status, drop_heavy=True)


@router.post(
    "/coordinate-file",
    operation_id="structureCoordinateFile",
    summary="Allowlisted bundled mmCIF text for the molecular viewer",
    response_model=ApiEnvelope,
)
def structure_coordinate_file(request: Request, body: StructureCoordinateFileRequest) -> ApiEnvelope:
    """Serves only 1CRN/1BNA/1RNA bundled copies. No user URL."""
    payload = bundled_coordinate_file(body.structure_id)
    return envelope(request, payload, scientific_status="EXPERIMENTAL", drop_heavy=False)
