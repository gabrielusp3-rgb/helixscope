"""Structure/variant/guide comparison. Does not hard-code RMSD or TM-score."""

from __future__ import annotations

from modules.comparative import (
    compare_guides,
    compare_proteins,
    compare_structures_remote,
    compare_structures_usalign,
    compare_variants,
    group_associated_positions,
    load_parsed_structure_for_entry,
    msa_column_to_structure,
    structure_residue_to_msa_column,
)
from modules.msa import column_detail, import_prealigned_fasta
from modules.rcsb_alignment import supported_methods
from modules.structure_alignment import lookup_correspondence, parse_alignment_payload
from modules.variant_core import build_variant

__all__ = (
    "compare_evolution",
    "compare_evidence_from_variants",
    "compare_guides",
    "compare_modes",
    "compare_proteins",
    "compare_structures_remote",
    "compare_structures_usalign",
    "compare_variants",
    "compare_variants_from_text",
    "group_associated_positions",
    "load_parsed_structure_for_entry",
    "lookup_correspondence",
    "msa_column_to_structure",
    "parse_alignment_payload",
    "structure_alignment_methods",
    "structure_residue_to_msa_column",
)


def compare_modes() -> dict:
    """Declared comparison modes. Not a generic JSON blob.

    Returns:
        Mode descriptors with dedicated endpoints. Evidence is the pack from
        identified variants, not a decorative empty panel.
    """
    return {
        "status": "COMPUTED",
        "modes": [
            {
                "id": "structures",
                "label": "Structures",
                "endpoints": [
                    "/api/v1/compare/parse",
                    "/api/v1/compare/jobs/remote",
                    "/api/v1/compare/jobs/usalign",
                ],
                "note": "RCSB Alignment API and local US-align. Source XYZ is never mutated.",
            },
            {
                "id": "variants",
                "label": "Variants",
                "endpoints": ["/api/v1/compare/variants"],
                "note": "Side-by-side identity and evidence. No pathogenicity ranking.",
            },
            {
                "id": "proteins",
                "label": "Proteins",
                "endpoints": ["/api/v1/compare/proteins"],
                "note": "Sequence identity is not RMSD and not InterPro.",
            },
            {
                "id": "guides",
                "label": "Guides",
                "endpoints": ["/api/v1/compare/guides"],
                "note": "Compares submitted CRISPR guide properties. No invented universal score.",
            },
            {
                "id": "evolution",
                "label": "Evolution",
                "endpoints": ["/api/v1/compare/evolution"],
                "note": "MSA column/group inspection. Not dN/dS.",
            },
            {
                "id": "evidence",
                "label": "Evidence",
                "endpoints": ["/api/v1/compare/evidence"],
                "note": "Evidence pack from two identified variants. Confidence stays null.",
            },
        ],
    }


def compare_evidence_from_variants(
    text_a: str,
    assembly_a: str,
    text_b: str,
    assembly_b: str,
) -> dict:
    """Evidence-mode view of two identified variants. Does not vote.

    Args:
        text_a: Variant A text.
        assembly_a: Mandatory assembly A.
        text_b: Variant B text.
        assembly_b: Mandatory assembly B.

    Returns:
        mode evidence plus the export pack from compare_variants.

    Raises:
        VariantInputError: Invalid identity.
        ComparativeError: Missing variant objects.
    """
    compared = compare_variants_from_text(text_a, assembly_a, text_b, assembly_b)
    export = compared.get("export") if isinstance(compared.get("export"), dict) else {}
    return {
        "mode": "evidence",
        "status": export.get("status") or "RETRIEVED",
        "ranking": None,
        "confidence_score": None,
        "export": export,
        "side_by_side": compared.get("side_by_side"),
        "conflicts": export.get("conflicts") or compared.get("conflicts") or [],
        "note": (
            "Evidence items are retrieved or computed records already attached to "
            "the two variants. HelixScope does not rank pathogenicity or resolve conflicts."
        ),
    }


def structure_alignment_methods() -> list[dict]:
    """RCSB Alignment API methods that HelixScope may submit."""
    return supported_methods()


def compare_variants_from_text(
    text_a: str,
    assembly_a: str,
    text_b: str,
    assembly_b: str,
) -> dict:
    """Identify two variants and compare explorer-shaped envelopes without remotes.

    Args:
        text_a: Variant A text.
        assembly_a: Mandatory assembly A.
        text_b: Variant B text.
        assembly_b: Mandatory assembly B.

    Returns:
        compare_variants payload. No pathogenicity ranking.

    Raises:
        VariantInputError: Invalid identity.
        ComparativeError: Missing variant objects.
    """
    va = build_variant(text=text_a, assembly=assembly_a)
    vb = build_variant(text=text_b, assembly=assembly_b)
    return compare_variants(
        {"variant": va, "layers": {}},
        {"variant": vb, "layers": {}},
    )


def compare_evolution(
    fasta: str,
    *,
    group_ids: list[str] | None = None,
    column: int = 0,
    member_id: str = "",
) -> dict:
    """Column/group inspection of a user-declared MSA. Not dN/dS.

    Args:
        fasta: Prealigned FASTA.
        group_ids: Optional member identifiers for group-associated positions.
        column: 0-based MSA column.
        member_id: Optional member for residue-to-column mapping.

    Returns:
        mode evolution, column detail, optional associated positions.

    Raises:
        MsaError: Invalid alignment.
    """
    msa_result = import_prealigned_fasta(fasta)
    detail = column_detail(msa_result, int(column))
    associated = None
    ids = [str(item).strip() for item in (group_ids or []) if str(item).strip()]
    if ids:
        associated = group_associated_positions(msa_result=msa_result, group_ids=ids)
    mapped = None
    if str(member_id or "").strip():
        mapped = structure_residue_to_msa_column(
            msa_result=msa_result,
            member_id=str(member_id).strip(),
            ungapped_index=0,
        )
    return {
        "mode": "evolution",
        "status": "COMPUTED",
        "alignment_hash": msa_result.get("alignment_hash"),
        "n_sequences": msa_result.get("n_sequences"),
        "column": detail,
        "group_associated": associated,
        "member_mapping": mapped,
        "note": (
            "Conservation is Shannon among the analyzed sequences, not functional "
            "essentiality. Group-associated residue patterns are not dN/dS."
        ),
        "topology_unchanged": True,
    }
