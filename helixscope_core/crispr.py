"""CRISPR facade. TEST REFERENCE is not GRCh38. Unavailable scores stay None."""

from __future__ import annotations

from typing import Any, List

from helixscope_core.sequence_input import resolve_pasted_sequence
from modules.crispr import (
    CAS_SYSTEMS,
    GUIDE_LENGTH_BY_CAS,
    PAM_ORIENTATION,
    PAM_PATTERNS,
    PFS_FORBIDDEN_CAS13,
    design_primer_pair,
    evaluate_guides,
    find_guides,
    find_off_targets,
    model_availability,
    rank_guides,
    system_info,
)

__all__ = (
    "analyze_crispr",
    "design_guides",
    "design_primer_pair",
    "evaluate_guides",
    "find_guides",
    "find_off_targets",
    "list_cas_systems",
    "model_availability",
    "rank_guides",
    "require_cas_system",
)


def require_cas_system(cas_system: str) -> str:
    """Return the canonical CAS_SYSTEMS key. Never alias informal tokens.

    Args:
        cas_system: Exact Core identifier.

    Returns:
        The same key if registered.

    Raises:
        ValueError: Unknown token. Informal names such as Cas12a are not mapped.
    """
    key = str(cas_system or "").strip()
    if key not in CAS_SYSTEMS:
        supported = ", ".join(CAS_SYSTEMS)
        raise ValueError(
            f"Sistema Cas invalido. Use um de: {supported}. "
            "Informal tokens such as Cas12a, CBE, ABE, or Prime Editor are not "
            "accepted; the UI must send the canonical_key from GET /api/v1/crispr/systems."
        )
    return key


def list_cas_systems() -> dict[str, Any]:
    """Canonical Cas-system registry for API/UI. Not a React-side list.

    Returns:
        systems rows plus shared model_availability. Status COMPUTED for the
        registry itself; published on-target models remain UNAVAILABLE.

    Raises:
        None.
    """
    models = model_availability()
    unavailable_published = [
        {
            "metric": "on_target_doench_ruleset2",
            "status": "UNAVAILABLE",
            "reason": (models.get("on_target_doench_ruleset2") or {}).get("reason"),
        },
        {
            "metric": "on_target_azimuth",
            "status": "UNAVAILABLE",
            "reason": (models.get("on_target_azimuth") or {}).get("reason"),
        },
        {
            "metric": "on_target_deephf",
            "status": "UNAVAILABLE",
            "reason": (models.get("on_target_deephf") or {}).get("reason"),
        },
        {
            "metric": "genome_wide_cfd_mit",
            "status": "UNAVAILABLE",
            "reason": "Genome-wide CFD/MIT requires COMPLETED_FULL_REFERENCE.",
        },
    ]
    systems: list[dict[str, Any]] = []
    for key in CAS_SYSTEMS:
        info = system_info(key)
        backbone = info.get("pam_system")
        target = str(info.get("target_molecule") or "DNA")
        editor = str(info.get("editor") or "nuclease")
        if target == "RNA":
            guide_length = int(GUIDE_LENGTH_BY_CAS.get("Cas13") or 28)
            pam = None
            pam_orientation = None
            pfs = {
                "kind": "PFS",
                "forbidden_3prime": sorted(PFS_FORBIDDEN_CAS13),
                "note": str(info.get("pam_display") or "PFS 3' non-G (LwaCas13a)"),
            }
        else:
            guide_length = int(GUIDE_LENGTH_BY_CAS.get(str(backbone or key)) or 20)
            pam = PAM_PATTERNS.get(str(backbone)) if backbone else None
            pam_orientation = PAM_ORIENTATION.get(str(backbone)) if backbone else None
            pfs = None
        supported = [
            "guide_sequence",
            "gc_content",
            "poly_t",
            "self_complementarity",
            "heuristic_efficiency",
        ]
        unsupported = list(unavailable_published)
        mode = "nuclease"
        limitations = [str(info.get("notes") or "")]
        if editor == "nuclease_rna":
            mode = "rna_nuclease"
            supported.append("pfs_filter")
            unsupported.append(
                {
                    "metric": "dna_cut_site",
                    "status": "METHOD_NOT_APPLICABLE",
                    "reason": "Cas13 targets RNA and does not produce a site-specific DNA DSB.",
                }
            )
            limitations.append("DNA target material is rejected. Provide RNA (A/U/C/G).")
        elif editor == "base_editor":
            mode = "base_editor"
            supported.extend(["edit_from", "edit_to", "edit_window"])
            unsupported.append(
                {
                    "metric": "double_strand_break",
                    "status": "METHOD_NOT_APPLICABLE",
                    "reason": "Base editors nick; they are not nucleases.",
                }
            )
        elif editor == "prime_editor":
            mode = "prime_editor"
            supported.append("prime_editing_notes")
            unsupported.append(
                {
                    "metric": "pegrna_design",
                    "status": "PARTIAL",
                    "reason": "Core records nick/PBS/RTT notes; it does not invent a pegRNA.",
                }
            )
        else:
            supported.append("cut_site")
        systems.append(
            {
                "canonical_key": key,
                "display_name": str(info.get("label") or key),
                "target_molecule": target,
                "guide_length": guide_length,
                "pam": pam,
                "pam_display": info.get("pam_display"),
                "pam_orientation": pam_orientation,
                "pfs": pfs,
                "mode": mode,
                "editor": editor,
                "cut_type": info.get("cut_type"),
                "edit_from": info.get("edit_from"),
                "edit_to": info.get("edit_to"),
                "edit_window": list(info["edit_window"]) if info.get("edit_window") else None,
                "capabilities": supported,
                "supported_metrics": supported,
                "unsupported_metrics": unsupported,
                "scientific_limitations": [item for item in limitations if item],
                "heuristic_efficiency_status": "HEURISTIC",
            }
        )
    return {
        "status": "COMPUTED",
        "systems": systems,
        "model_availability": models,
        "note": (
            "canonical_key must be sent unchanged. Informal UI tokens are invalid. "
            "No system silently falls back to SpCas9."
        ),
    }


def design_guides(
    sequence: str,
    cas_system: str = "SpCas9",
    *,
    max_mismatches: int = 3,
    run_off_target: bool = False,
) -> List[dict[str, Any]]:
    """Find and evaluate guides on the provided sequence (not genome-wide).

    Args:
        sequence: Target DNA (or Cas13 RNA) or single-record FASTA. A 20-nt
            spacer alone is not a substitute for a PAM-containing locus.
        cas_system: Canonical CAS_SYSTEMS key.
        max_mismatches: Off-target mismatch cap when scanning the same sequence.
        run_off_target: Sequence-local scan only. Not genome-wide.

    Returns:
        Evaluated guide list from Core.

    Raises:
        TypeError: Non-string input.
        ValueError: FASTA parse failure, invalid target alphabet, or unknown Cas.
    """
    key = require_cas_system(cas_system)
    parsed = resolve_pasted_sequence(sequence)
    target = str(parsed.get("sequence") or "")
    info = system_info(key)
    backbone = str(info.get("pam_system") or key)
    glen = int(GUIDE_LENGTH_BY_CAS.get(backbone) or GUIDE_LENGTH_BY_CAS.get(key) or 20)
    found = find_guides(target, key, guide_length=glen)
    return evaluate_guides(
        found,
        target,
        key,
        max_mismatches=max_mismatches,
        run_off_target=run_off_target,
    )


def analyze_crispr(
    sequence: str,
    cas_system: str = "SpCas9",
    *,
    max_mismatches: int = 3,
    run_off_target: bool = False,
) -> dict[str, Any]:
    """Complete CRISPR investigation on a provided PAM-containing locus.

    Args:
        sequence: Target locus DNA (or Cas13 RNA) or single-record FASTA.
        cas_system: Canonical CAS_SYSTEMS key.
        max_mismatches: Off-target mismatch cap for a sequence-local scan.
        run_off_target: Sequence-local scan only. Not genome-wide.

    Returns:
        Domain result with guides, scope flags, model availability, and
        optional validation primers. Status HEURISTIC.

    Raises:
        TypeError: Non-string input.
        ValueError: FASTA parse failure, invalid target alphabet, or unknown Cas.
    """
    key = require_cas_system(cas_system)
    info = system_info(key)
    parsed = resolve_pasted_sequence(sequence)
    target = "".join(str(parsed.get("sequence") or "").split())
    guides = design_guides(
        sequence,
        key,
        max_mismatches=max_mismatches,
        run_off_target=run_off_target,
    )
    target_mol = str(info.get("target_molecule") or "DNA")
    if target_mol == "RNA":
        input_note = (
            "Cas13 requires RNA (A/U/C/G) including PFS context. DNA (T) is not "
            "transcribed automatically. This is not a DNA DSB designer."
        )
    else:
        pam = str(info.get("pam_display") or "PAM")
        input_note = (
            f"{key} needs PAM-containing target DNA ({pam}), not a spacer pasted "
            "alone. A spacer without PAM context is not a locus."
        )
    payload: dict[str, Any] = {
        "status": "HEURISTIC",
        "cas_system": key,
        "cas_display_name": str(info.get("label") or key),
        "target_molecule": target_mol,
        "editor": info.get("editor"),
        "guides": guides,
        "n_guides": len(guides),
        "target_length": len(target),
        "target_sequence": target,
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
        "off_target_scope": "provided_sequence" if run_off_target else "not_scanned",
        "not_genome_wide": True,
        "model_availability": model_availability(),
        "spacer_is_not_target": True,
        "input_note": input_note,
        "system": {
            row["canonical_key"]: row
            for row in list_cas_systems()["systems"]
            if row["canonical_key"] == key
        }.get(key)
        or {},
    }
    if guides and target_mol != "RNA":
        try:
            payload["validation_primers"] = design_primer_pair(
                str(guides[0].get("guide_sequence") or ""),
                target,
            )
        except (TypeError, ValueError) as exc:
            payload["validation_primers"] = {
                "status": "UNAVAILABLE",
                "reason": str(exc),
            }
    elif target_mol == "RNA":
        payload["validation_primers"] = {
            "status": "METHOD_NOT_APPLICABLE",
            "reason": "DNA PCR validation primers are not designed for Cas13 RNA targeting.",
        }
    return payload
