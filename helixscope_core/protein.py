"""Protein public surface. Sequence SS predictors remain UNAVAILABLE."""

from __future__ import annotations

from typing import Any

from helixscope_core.sequence_input import resolve_pasted_sequence
from modules import provenance
from modules.protein_analysis import (
    STANDARD_AMINO_ACIDS,
    aliphatic_index,
    amino_acid_categories,
    amino_acid_composition,
    aromaticity,
    charge_count_profile,
    gravy_index,
    hydrophobicity_profile,
    hydrophobicity_profile_records,
    instability_index,
    isoelectric_point,
    molecular_weight_protein,
    net_charge,
    physicochemical_report,
    secondary_structure_availability,
    translate,
)

__all__ = (
    "aliphatic_index",
    "amino_acid_categories",
    "amino_acid_composition",
    "analyze_protein",
    "aromaticity",
    "charge_count_profile",
    "gravy_index",
    "hydrophobicity_profile",
    "hydrophobicity_profile_records",
    "instability_index",
    "isoelectric_point",
    "molecular_weight_protein",
    "net_charge",
    "physicochemical_report",
    "retrieve_uniprot_entry",
    "secondary_structure_availability",
    "translate",
)


def analyze_protein(sequence: str, *, ph: float = 7.0) -> dict[str, Any]:
    """Standard-protein physicochemical package plus provenance.

    Args:
        sequence: Amino acid sequence or single-record FASTA. Headers are not
            residues.
        ph: Net-charge pH.

    Returns:
        physicochemical_report fields plus composition and provenance.
        Raises from protein_analysis on invalid input.

    Raises:
        TypeError: Non-string input.
        ValueError: Multi-FASTA, empty residues after parse, or symbols outside
            the 20 canonical amino acids.
    """
    parsed = resolve_pasted_sequence(sequence)
    residues = "".join(str(parsed.get("sequence") or "").split()).upper()
    if not residues:
        raise ValueError(
            "Expected a protein sequence containing the 20 canonical amino-acid "
            "symbols. The input was empty after FASTA and whitespace parsing."
        )
    invalid = sorted({residue for residue in residues if residue not in STANDARD_AMINO_ACIDS})
    if invalid:
        raise ValueError(
            "Expected a protein sequence containing the 20 canonical amino-acid "
            f"symbols. Found unsupported residues: {', '.join(invalid)}."
        )
    report = physicochemical_report(residues, ph=ph)
    composition = amino_acid_composition(residues)
    payload: dict[str, Any] = {
        **report,
        "sequence": residues,
        "composition": composition,
        "ss_prediction": secondary_structure_availability(),
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
        "sequence_hash": provenance.sequence_digest(residues),
        "provenance": provenance.build_provenance(
            status=str(report.get("status") or "COMPUTED"),
            algorithm="protein_analysis.physicochemical_report",
            parameters={"ph": ph, "input_format": parsed.get("format")},
            input_identifier=str(parsed.get("identifier") or ""),
        ),
    }
    try:
        payload["hydropathy_profile"] = hydrophobicity_profile_records(residues)
        payload["hydropathy_status"] = "COMPUTED"
    except ValueError as exc:
        payload["hydropathy_profile"] = []
        payload["hydropathy_status"] = "UNAVAILABLE"
        payload["hydropathy_note"] = str(exc)
    try:
        scores = charge_count_profile(residues)
        window = 9
        payload["charge_profile"] = [
            {
                "start": index,
                "end": index + window,
                "midpoint": index + window // 2,
                "score": score,
                "window": window,
                "method": "formal side-chain count (K+R-D-E)/window; His excluded",
            }
            for index, score in enumerate(scores)
        ]
        payload["charge_profile_status"] = "COMPUTED"
    except ValueError as exc:
        payload["charge_profile"] = []
        payload["charge_profile_status"] = "UNAVAILABLE"
        payload["charge_profile_note"] = str(exc)
    try:
        payload["amino_acid_categories"] = amino_acid_categories(residues)
        payload["amino_acid_categories_status"] = "COMPUTED"
    except ValueError as exc:
        payload["amino_acid_categories"] = {}
        payload["amino_acid_categories_status"] = "UNAVAILABLE"
        payload["amino_acid_categories_note"] = str(exc)
    payload["identity_scope"] = "sequence_derived"
    payload["identity_note"] = (
        "Physicochemical values were computed from the submitted residue string. "
        "UniProt names, organism and domains appear only after an explicit "
        "accession retrieval, not by guessing identity from sequence appearance."
    )
    return payload


def retrieve_uniprot_entry(accession: str, *, include_domains: bool = False) -> dict[str, Any]:
    """Retrieve UniProtKB annotation for an explicit accession. Makes network.

    Args:
        accession: UniProtKB accession, not an anonymous sequence.
        include_domains: If True, also query InterPro for that accession.

    Returns:
        UniProt fields with evidence_status RETRIEVED. Optional InterPro block.

    Raises:
        DomainError: Invalid accession or remote failure categories.
    """
    from modules import protein_domains

    entry = protein_domains.fetch_uniprot_protein(accession=accession)
    payload: dict[str, Any] = dict(entry)
    payload["status"] = "RETRIEVED"
    if include_domains:
        payload["interpro"] = protein_domains.fetch_protein_domains(accession=accession)
    return payload
