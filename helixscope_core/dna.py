"""DNA public surface. Formulas live in ``modules.dna_analysis``."""

from __future__ import annotations

from typing import Any

from helixscope_core.sequence_input import resolve_pasted_sequence
from modules import provenance, thermodynamics
from modules.dna_analysis import (
    at_content,
    at_skew,
    composition_is_consistent,
    cpg_islands,
    dinucleotide_frequencies,
    entropy_sliding_window,
    find_orfs,
    find_restriction_sites,
    gc_content,
    gc_skew,
    gc_skew_global,
    kmer_counts,
    kmer_summary,
    melting_temperature,
    melting_temperature_report,
    molecular_weight,
    nucleotide_composition,
    orf_coordinates_on_input,
    parse_sequence_payload,
    reverse_complement,
    select_non_overlapping_orfs,
    shannon_entropy,
    validate_for_molecule,
    validate_sequence,
    windowed_profiles,
)

__all__ = (
    "analyze_dna",
    "at_content",
    "at_skew",
    "composition_is_consistent",
    "cpg_islands",
    "dinucleotide_frequencies",
    "dna_3d_viewer",
    "entropy_sliding_window",
    "find_orfs",
    "find_restriction_sites",
    "gc_content",
    "gc_skew",
    "gc_skew_global",
    "kmer_counts",
    "kmer_summary",
    "melting_temperature",
    "melting_temperature_report",
    "molecular_weight",
    "nucleotide_composition",
    "orf_coordinates_on_input",
    "parse_sequence_payload",
    "reverse_complement",
    "select_non_overlapping_orfs",
    "shannon_entropy",
    "santalucia_tm",
    "validate_for_molecule",
    "validate_sequence",
    "windowed_profiles",
)


def santalucia_tm(
    sequence: str,
    *,
    na_m: float = 0.05,
    oligo_nm: float = 250.0,
    mg_m: float = 0.0,
) -> dict:
    """SantaLucia nearest-neighbor Tm. Same callable as thermodynamics."""
    return thermodynamics.santalucia_nearest_neighbor_tm(
        sequence, na_m=na_m, oligo_nm=oligo_nm, mg_m=mg_m
    )


def analyze_dna(
    sequence: str,
    *,
    min_orf_length: int = 150,
    include_orfs: bool = False,
    include_restriction: bool = False,
    drop_overlapping_orfs: bool = True,
    include_profiles: bool = False,
    include_cpg: bool = False,
    include_kmers: bool = False,
    include_santalucia: bool = False,
    profile_window: int = 100,
    profile_step: int = 50,
    kmer_k: int = 3,
    kmer_top_n: int = 30,
    santalucia_na_m: float = 0.05,
    santalucia_oligo_nm: float = 250.0,
) -> dict[str, Any]:
    """Validate DNA and compute the domain investigation package.

    Args:
        sequence: Raw DNA or single-record FASTA (Streamlit pasted-sequence contract).
        min_orf_length: ORF length floor when include_orfs is True.
        include_orfs: If True, run find_orfs and optionally non-overlapping selection.
        include_restriction: If True, run restriction scan.
        drop_overlapping_orfs: If True, keep one ORF per locus (legacy default).
        include_profiles: If True, compute windowed GC/AT/skew/entropy arrays.
        include_cpg: If True, run the CpG-island detector.
        include_kmers: If True, rank k-mers.
        include_santalucia: If True, run SantaLucia nearest-neighbor Tm when applicable.
        profile_window: Sliding window length in bp.
        profile_step: Window step in bp.
        kmer_k: k-mer length.
        kmer_top_n: How many abundant k-mers to keep.
        santalucia_na_m: Na+ molarity for SantaLucia.
        santalucia_oligo_nm: Oligo concentration in nM for SantaLucia.

    Returns:
        Dict with validation, metrics, optional feature arrays, sequence_hash,
        provenance. GC may be NaN. status ERROR when the sequence is not valid
        DNA after FASTA parse. Anonymous sequence analysis does not invent
        organism, gene, or chromosome identity.

    Raises:
        TypeError: Propagated from validators.
        ValueError: Multi-FASTA or FASTA parse failure (same as Streamlit).
    """
    parsed = resolve_pasted_sequence(sequence)
    residues = str(parsed.get("sequence") or "")
    info = validate_for_molecule(residues, "DNA")
    clean = str(info.get("sequence") or "")
    if not info.get("is_valid"):
        return {
            "status": "ERROR",
            "validation": info,
            "sequence": clean,
            "input_format": parsed.get("format"),
            "input_identifier": parsed.get("identifier") or "",
            "sequence_hash": provenance.sequence_digest(clean) if clean else "",
            "provenance": provenance.build_provenance(
                status="ERROR",
                algorithm="dna_analysis.validate_for_molecule",
                parameters={"molecule": "DNA", "input_format": parsed.get("format")},
                input_identifier=str(parsed.get("identifier") or ""),
            ),
        }
    payload: dict[str, Any] = {
        "status": "COMPUTED",
        "validation": info,
        "sequence": clean,
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
        "length": int(info.get("length") or len(clean)),
        "sequence_hash": provenance.sequence_digest(clean),
        "gc_content": gc_content(clean),
        "at_content": at_content(clean),
        "gc_skew": gc_skew_global(clean),
        "at_skew": at_skew(clean),
        "composition": nucleotide_composition(clean),
        "dinucleotides": dinucleotide_frequencies(clean),
        "reverse_complement": reverse_complement(clean),
        "shannon_entropy": shannon_entropy(clean),
        "molecular_weight_dna": molecular_weight(clean, "DNA"),
        "tm_report": melting_temperature_report(clean),
        "composition_consistent": composition_is_consistent(clean),
        "provenance": provenance.build_provenance(
            status="COMPUTED",
            algorithm="dna_analysis core metrics",
            parameters={"molecule": "DNA", "input_format": parsed.get("format")},
            input_identifier=str(parsed.get("identifier") or ""),
        ),
    }
    payload["provenance"]["input_hash"] = payload["sequence_hash"]
    if include_orfs:
        raw_orfs = find_orfs(clean, min_orf_length)
        payload["orfs_raw_count"] = len(raw_orfs)
        selected = (
            select_non_overlapping_orfs(raw_orfs, min_length=min_orf_length)
            if drop_overlapping_orfs
            else raw_orfs
        )
        payload["orfs"] = selected
        payload["orfs_drop_overlapping"] = bool(drop_overlapping_orfs)
        features: list[dict[str, Any]] = []
        length = int(payload["length"])
        for orf in selected:
            if not isinstance(orf, dict):
                continue
            try:
                start, end = orf_coordinates_on_input(orf, length)
            except (ValueError, TypeError, KeyError):
                continue
            features.append(
                {
                    "start": start,
                    "end": end,
                    "label": f"Predicted ORF frame {orf.get('frame')}",
                    "kind": "orf",
                    "strand": orf.get("strand"),
                    "frame": orf.get("frame"),
                    "status": orf.get("status"),
                }
            )
        payload["sequence_features"] = features
    if include_restriction:
        sites = find_restriction_sites(clean)
        payload["restriction_sites"] = sites
        enzyme_rows: list[dict[str, Any]] = []
        hits: list[dict[str, Any]] = []
        if isinstance(sites, dict):
            for enzyme, data in sites.items():
                spec = data if isinstance(data, dict) else {}
                enzyme_rows.append(
                    {
                        "enzyme": enzyme,
                        "pattern": spec.get("pattern"),
                        "count": spec.get("count"),
                        "sites": spec.get("sites"),
                    }
                )
                raw_hits = spec.get("hits") or []
                if isinstance(raw_hits, list):
                    hits.extend(item for item in raw_hits if isinstance(item, dict))
        payload["restriction_enzyme_rows"] = enzyme_rows
        payload["restriction_hits"] = hits
    if include_profiles:
        try:
            payload["windowed_profiles"] = windowed_profiles(
                clean,
                window=profile_window,
                step=profile_step,
            )
            payload["profile_window"] = profile_window
            payload["profile_step"] = profile_step
        except ValueError as exc:
            payload["windowed_profiles"] = []
            payload["profiles_status"] = "UNAVAILABLE"
            payload["profiles_note"] = str(exc)
    if include_cpg:
        payload["cpg_islands"] = cpg_islands(clean)
    if include_kmers:
        try:
            payload["kmer_summary"] = kmer_summary(clean, k=kmer_k, top_n=kmer_top_n)
        except ValueError as exc:
            payload["kmer_summary"] = {"status": "UNAVAILABLE", "top": [], "note": str(exc)}
    if include_santalucia:
        payload["santalucia_tm"] = santalucia_tm(
            clean,
            na_m=santalucia_na_m,
            oligo_nm=santalucia_oligo_nm,
        )
    tracks: list[dict[str, Any]] = []
    orf_features = payload.get("sequence_features") if isinstance(payload.get("sequence_features"), list) else []
    if orf_features:
        tracks.append({"name": "Predicted ORFs", "features": orf_features})
    if include_cpg:
        cpg_features: list[dict[str, Any]] = []
        for island in payload.get("cpg_islands") or []:
            if not isinstance(island, dict):
                continue
            cpg_features.append(
                {
                    "start": island.get("start"),
                    "end": island.get("end"),
                    "label": "CpG island",
                    "kind": "cpg",
                }
            )
        if cpg_features:
            tracks.append({"name": "CpG islands", "features": cpg_features})
    if tracks:
        payload["sequence_feature_tracks"] = tracks
    payload["identity_scope"] = "sequence_derived"
    payload["identity_note"] = (
        "This analysis used only the submitted DNA string. Organism, gene, "
        "chromosome and function are unknown unless a separate NCBI or BLAST "
        "retrieval provides them."
    )
    return payload


def dna_3d_viewer(
    sequence: str,
    *,
    prefer: str = "illustrative",
    region_start: int = 0,
    region_end: int | None = None,
    lod: str = "high",
    inspect_position: int = 0,
    center_on_selected: bool = False,
) -> dict[str, Any]:
    """DNA 3D coordinates from Core. ILLUSTRATIVE helix is not experimental.

    Args:
        sequence: Raw DNA or single-record FASTA.
        prefer: ``illustrative``, ``experimental``, ``predicted``, or ``auto``.
        region_start: Inclusive 0-based start of the illustrative helix window.
        region_end: Exclusive end. ``None`` or 0 lets Core use the LOD cap.
        lod: Visual LOD (low/medium/high). Caps consecutive nt; does not
            interpolate omitted bases.
        inspect_position: 0-based index on the full analysis sequence.
        center_on_selected: If True, Core recenters the illustrative window.

    Returns:
        Status/kind plus a viewer payload of deposited or illustrative XYZ,
        plus position mapping for ``inspect_position``. PREDICTED remains
        UNAVAILABLE without a deposited envelope.

    Raises:
        TypeError: Non-string input.
        ValueError: Multi-FASTA / FASTA parse failure.

    Nota biologica:
        Implementacao heuristica simplificada inspirada em parametros
        helicoidais B-DNA (rise 3.38 A, 10 bp/volta, twist 36 deg; fibre
        B-DNA classico, nao 10.5 bp/volta em solucao);
        nao reproduz o modelo/algoritmo original publicado de um
        refinador cristalografico. PREDICTED sem envelope depositado
        permanece UNAVAILABLE.
    """
    from helixscope_core.structure import viewer_payload
    from modules import dna_3d, region_nav

    parsed = resolve_pasted_sequence(sequence)
    seq = str(parsed.get("sequence") or "")
    mode = str(prefer or "illustrative").strip().lower()
    prefer_arg: str | None
    if mode in {"", "auto"}:
        prefer_arg = None
    else:
        prefer_arg = mode
    end = None if region_end in (None, 0) else int(region_end)
    visual_lod = region_nav.normalize_lod(lod)
    representation = dna_3d.representation_for_sequence(
        seq,
        prefer=prefer_arg,
        region_start=int(region_start),
        region_end=end,
        lod=visual_lod,
        center=int(inspect_position) if center_on_selected else None,
    )
    envelope = representation.get("envelope")
    viewer: dict[str, Any] = {}
    if isinstance(envelope, dict):
        viewer = viewer_payload(envelope)
        for banned in ("camera", "eye", "zoom", "mediapipe", "gesture", "rotation"):
            viewer.pop(banned, None)
    coords = list(viewer.get("coordinates") or viewer.get("atoms") or [])
    try:
        mapping = dna_3d.position_to_3d(representation, int(inspect_position))
    except dna_3d.Dna3dError:
        mapping = {
            "status": "ERROR",
            "query_index_0based": inspect_position,
            "mapping": None,
        }
    base = seq[inspect_position] if 0 <= int(inspect_position) < len(seq) else None
    return {
        "status": representation.get("status"),
        "kind": representation.get("kind"),
        "kind_label": representation.get("kind_label"),
        "disclaimer": representation.get("disclaimer"),
        "reason": representation.get("reason"),
        "sequence_hash": representation.get("sequence_hash"),
        "n_atoms": viewer.get("n_atoms") or len(coords),
        "partial_view": representation.get("partial_view"),
        "visual_lod": visual_lod,
        "helix_limit_nt": region_nav.helix_limit_for_lod(visual_lod),
        "region_start": representation.get("region_start", region_start),
        "region_end": representation.get("region_end", end),
        "inspect_position": int(inspect_position),
        "inspect_base": base,
        "position_mapping": mapping,
        "viewer": viewer,
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
    }
