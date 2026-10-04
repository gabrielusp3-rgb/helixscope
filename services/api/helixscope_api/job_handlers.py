"""Job handlers. Call helixscope_core only."""

from __future__ import annotations

from typing import Any

from helixscope_core.blast import run_local_blast, submit_search
from helixscope_core.compare import (
    compare_structures_remote,
    compare_structures_usalign,
    load_parsed_structure_for_entry,
)
from helixscope_core.genome import detect_cas_offinder, submit_casoffinder_job
from helixscope_core.msa import (
    build_collection_member,
    import_prealigned_fasta,
    parse_unaligned_fasta,
    submit_msa,
    validate_collection,
)
from helixscope_core.phylogeny import infer_tree
from helixscope_core.structure import read_bundled_mmcif
from helixscope_core.variants import explore_variant, identify_variant

from helixscope_api.jobs import JobRecord


def run_phylogeny_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    msa = import_prealigned_fasta(str(params.get("fasta") or ""))
    extra: dict[str, Any] = {}
    if params.get("distance_model"):
        extra["distance_model"] = params["distance_model"]
    payload = infer_tree(msa, method=str(params.get("method") or ""), **extra)
    payload["source_msa_hash"] = params.get("alignment_hash") or msa.get("alignment_hash")
    payload["source_msa_n_sequences"] = msa.get("n_sequences")
    return payload


def run_blast_local_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    return run_local_blast(
        program=str(params.get("program") or "blastn"),
        query=str(params.get("query") or ""),
        molecule=str(params.get("molecule") or "DNA"),
    )


def run_blast_remote_submit(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    return submit_search(
        str(params.get("query") or ""),
        str(params.get("email") or ""),
        program=str(params.get("program") or "blastn"),
    )


def run_msa_engine_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    rows = parse_unaligned_fasta(str(params.get("fasta") or ""))
    members = [
        build_collection_member(
            sequence=str(row.get("sequence") or row.get("ungapped") or ""),
            identifier=str(row.get("identifier") or f"seq_{index}"),
            source="api upload",
        )
        for index, row in enumerate(rows, start=1)
    ]
    validate_collection(members)
    return submit_msa(
        members,
        str(params.get("email") or ""),
        backend=str(params.get("backend") or "ebi_clustalo"),
    )


def run_casoffinder_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    from helixscope_core.errors import OpenCLError

    detected = detect_cas_offinder()
    if not detected.get("available"):
        raise OpenCLError(str(detected.get("reason") or "Cas-OFFinder is not available."), "UNAVAILABLE")
    guide = {
        "guide_sequence": str(params.get("guide_sequence") or ""),
        "pam_sequence": str(params.get("pam") or "NGG"),
        "cas_system": str(params.get("cas_system") or "SpCas9"),
    }
    return submit_casoffinder_job(guide=guide, assembly_id=str(params.get("assembly_id") or ""))


def run_compare_remote_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    fetch_coordinates = bool(params.get("fetch_coordinates", True))
    reference_entry = str(params.get("reference_entry") or "")
    target_entry = str(params.get("target_entry") or "")
    reference_parsed = None
    target_parsed = None
    reference_meta: dict[str, Any] = {"requested": fetch_coordinates, "status": "NOT_REQUESTED"}
    target_meta: dict[str, Any] = {"requested": fetch_coordinates, "status": "NOT_REQUESTED"}
    if fetch_coordinates:
        reference_load = load_parsed_structure_for_entry(reference_entry)
        target_load = load_parsed_structure_for_entry(target_entry)
        reference_meta = {key: value for key, value in reference_load.items() if key != "parsed"}
        target_meta = {key: value for key, value in target_load.items() if key != "parsed"}
        if reference_load.get("status") == "AVAILABLE":
            reference_parsed = reference_load.get("parsed")
        if target_load.get("status") == "AVAILABLE":
            target_parsed = target_load.get("parsed")
    payload = compare_structures_remote(
        reference_entry=reference_entry,
        target_entry=target_entry,
        reference_chain=str(params.get("reference_chain") or "A"),
        target_chain=str(params.get("target_chain") or "A"),
        method=str(params.get("method") or "tm-align"),
        reference_parsed=reference_parsed,
        target_parsed=target_parsed,
    )
    overlay_available = payload.get("superposition") is not None
    payload["coordinate_fetch"] = {
        "requested": fetch_coordinates,
        "reference": reference_meta,
        "target": target_meta,
        "overlay_status": "AVAILABLE" if overlay_available else "UNAVAILABLE",
        "note": (
            "Cartesian overlay uses helixscope_core superposition_bundle on "
            "allowlisted RCSB/AlphaFold coordinate files. Residue-pair indexes "
            "are not converted into invented XYZ."
        ),
    }
    return payload


def run_usalign_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    reference_entry = str(params.get("reference_entry") or "").upper()
    target_entry = str(params.get("target_entry") or "").upper()
    reference_parsed = {"raw_structure_text": read_bundled_mmcif(f"{reference_entry}.cif")}
    target_parsed = {"raw_structure_text": read_bundled_mmcif(f"{target_entry}.cif")}
    payload = compare_structures_usalign(
        reference_parsed=reference_parsed,
        target_parsed=target_parsed,
        reference_entry=reference_entry,
        target_entry=target_entry,
        reference_chain=str(params.get("reference_chain") or "A"),
        target_chain=str(params.get("target_chain") or "A"),
    )
    payload["coordinate_fetch"] = {
        "requested": False,
        "overlay_status": "UNAVAILABLE",
        "note": (
            "US-align rotates structure 1 onto 2. This path does not provide "
            "RCSB target-transform overlay atoms."
        ),
    }
    return payload


def run_variant_explore_job(params: dict[str, Any], record: JobRecord) -> dict[str, Any]:
    variant = identify_variant(str(params.get("text") or ""), assembly=str(params.get("assembly") or ""))
    return explore_variant(
        variant=variant,
        email=str(params.get("email") or ""),
        enable_vep=bool(params.get("enable_vep")),
        enable_clinvar=bool(params.get("enable_clinvar")),
        enable_domains=bool(params.get("enable_domains")),
    )
