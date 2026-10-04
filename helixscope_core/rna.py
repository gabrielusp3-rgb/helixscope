"""RNA public surface. Folding stays ViennaRNA; no invented structure."""

from __future__ import annotations

from typing import Any

from helixscope_core.sequence_input import resolve_pasted_sequence
from helixscope_core.serialize import dataframe_to_records
from modules import provenance
from modules.dna_analysis import (
    at_content,
    dinucleotide_frequencies,
    gc_content,
    kmer_summary,
    nucleotide_composition,
    windowed_profiles,
)
from modules.rna_analysis import (
    codon_adaptation_index,
    codon_adaptation_index_report,
    codon_usage_table,
    effective_number_of_codons,
    effective_number_of_codons_report,
    kmer_counts,
    relative_synonymous_codon_usage,
    resolve_coding_codons,
    shannon_entropy,
    transcribe,
)
from modules.rna_folding import FoldingError, circular_layout, fold_rna

__all__ = (
    "analyze_rna",
    "codon_adaptation_index",
    "codon_adaptation_index_report",
    "codon_usage_table",
    "codon_usage_records",
    "effective_number_of_codons",
    "effective_number_of_codons_report",
    "fold_rna",
    "kmer_counts",
    "relative_synonymous_codon_usage",
    "resolve_coding_codons",
    "rna_3d_viewer",
    "rscu_records",
    "shannon_entropy",
    "transcribe",
    "transcribe_dna",
    "translate_coding_composition",
)


def transcribe_dna(dna_seq: str) -> str:
    """DNA T->U transcription. Same as ``rna_analysis.transcribe``.

    Args:
        dna_seq: Coding-strand DNA or single-record FASTA.

    Returns:
        RNA with U in place of T.

    Raises:
        TypeError: Non-string input.
        ValueError: FASTA parse failure or invalid DNA for transcription.
    """
    parsed = resolve_pasted_sequence(dna_seq)
    return transcribe(str(parsed.get("sequence") or ""))


def codon_usage_records(mrna_seq: str) -> dict[str, Any]:
    """Codon-usage table as records. DataFrame remains available internally."""
    return dataframe_to_records(codon_usage_table(mrna_seq))


def rscu_records(mrna_seq: str) -> dict[str, Any]:
    """RSCU table as records."""
    return dataframe_to_records(relative_synonymous_codon_usage(mrna_seq))


def analyze_rna(
    sequence: str,
    *,
    fold: bool = False,
    backend: str = "auto",
    include_codon_metrics: bool = False,
    include_profiles: bool = False,
    include_kmers: bool = False,
    profile_window: int = 100,
    profile_step: int = 50,
    kmer_k: int = 3,
) -> dict[str, Any]:
    """Validate RNA and compute the domain investigation package.

    Args:
        sequence: RNA input (A, U, C, G).
        fold: If True, call ViennaRNA MFE (may be UNAVAILABLE).
        backend: fold_rna backend token.
        include_codon_metrics: If True, compute codon usage, RSCU, ENC and CAI
            from this RNA string as a codon table. That does not make the
            string a validated CDS.
        include_profiles: If True, compute windowed composition arrays.
        include_kmers: If True, rank k-mers after U to T for the DNA counter.
        profile_window: Sliding window length.
        profile_step: Window step.
        kmer_k: k-mer length.

    Returns:
        status COMPUTED or ERROR. Fold sub-object keeps PREDICTED/UNAVAILABLE.
        Codon metrics, when requested, stay labeled as codon-table statistics.

    Raises:
        TypeError: From validators.
        FoldingError: If fold=True and the engine fails (not missing).
    """
    from modules.dna_analysis import validate_for_molecule

    parsed = resolve_pasted_sequence(sequence)
    info = validate_for_molecule(str(parsed.get("sequence") or ""), "RNA")
    clean = str(info.get("sequence") or "")
    if not info.get("is_valid"):
        return {
            "status": "ERROR",
            "validation": info,
            "sequence": clean,
            "input_format": parsed.get("format"),
            "input_identifier": parsed.get("identifier") or "",
            "provenance": provenance.build_provenance(
                status="ERROR",
                algorithm="validate_for_molecule",
                parameters={"molecule": "RNA", "input_format": parsed.get("format")},
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
        "shannon_entropy": shannon_entropy(clean),
        "composition": nucleotide_composition(clean),
        "dinucleotides": dinucleotide_frequencies(clean),
        "gc_content": gc_content(clean),
        "au_content": at_content(clean),
        "provenance": provenance.build_provenance(
            status="COMPUTED",
            algorithm="rna_analysis core metrics",
            parameters={"molecule": "RNA", "input_format": parsed.get("format")},
            input_identifier=str(parsed.get("identifier") or ""),
        ),
    }
    payload["provenance"]["input_hash"] = payload["sequence_hash"]
    if fold:
        folded = fold_rna(clean, identifier="helixscope_core.analyze_rna", backend=backend)
        if isinstance(folded, dict) and folded.get("sequence"):
            try:
                folded["circular_layout"] = circular_layout(
                    str(folded.get("sequence") or ""),
                    list(folded.get("base_pairs") or []),
                )
            except FoldingError:
                folded["circular_layout"] = None
        payload["fold"] = folded
    if include_codon_metrics:
        payload["codon_metrics_scope"] = (
            "Codon usage, RSCU, ENC and CAI are computed from this RNA string "
            "as a codon table. That is not independent evidence that the string "
            "is a coding sequence of a known gene."
        )
        try:
            payload["codon_usage"] = codon_usage_records(clean)
            payload["rscu"] = rscu_records(clean)
            payload["cai"] = codon_adaptation_index_report(clean)
            payload["enc"] = effective_number_of_codons_report(clean)
        except (ValueError, TypeError) as exc:
            payload["codon_metrics_status"] = "UNAVAILABLE"
            payload["codon_metrics_note"] = str(exc)
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
    if include_kmers:
        try:
            payload["kmer_summary"] = kmer_summary(clean.replace("U", "T"), k=kmer_k)
            payload["kmer_alphabet_note"] = (
                "RNA k-mers use the same Core counter as DNA after U to T."
            )
        except ValueError as exc:
            payload["kmer_summary"] = {"status": "UNAVAILABLE", "top": [], "note": str(exc)}
    payload["identity_scope"] = "sequence_derived"
    payload["identity_note"] = (
        "This analysis used only the submitted RNA string. Organism, gene and "
        "transcript identity are unknown unless a separate database retrieval "
        "provides them."
    )
    return payload


def rna_3d_viewer(
    sequence: str,
    *,
    prefer: str = "illustrative",
    region_start: int = 0,
    region_end: int | None = None,
    lod: str = "high",
    inspect_position: int = 0,
    center_on_selected: bool = False,
) -> dict[str, Any]:
    """RNA 3D coordinates from Core. ILLUSTRATIVE A-RNA is not experimental.

    Args:
        sequence: Raw RNA or single-record FASTA.
        prefer: ``illustrative``, ``experimental``, ``predicted``, or ``auto``.
        region_start: Inclusive 0-based start of the illustrative helix window.
        region_end: Exclusive end. ``None`` or 0 lets Core use the LOD cap.
        lod: Visual LOD (low/medium/high).
        inspect_position: 0-based index on the full analysis sequence.
        center_on_selected: If True, Core recenters the illustrative window.

    Returns:
        Status/kind plus a viewer payload of deposited or illustrative XYZ.
        ViennaRNA MFE is never converted into coordinates. PREDICTED 3D remains
        UNAVAILABLE without a deposited envelope.

    Raises:
        TypeError: Non-string input.
        ValueError: Multi-FASTA / FASTA parse failure.

    Nota biologica:
        Implementação heurística simplificada inspirada em parâmetros
        helicoidais A-RNA; não reproduz o modelo/algoritmo original publicado
        de um refinador cristalográfico. MFE/dot-bracket não é estrutura 3D.
    """
    from helixscope_core.structure import viewer_payload
    from modules import dna_3d, region_nav, rna_3d

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
    representation = rna_3d.representation_for_sequence(
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
        "mfe_does_not_yield_3d": True,
        "viewer": viewer,
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
    }


def translate_coding_composition(sequence: str) -> dict[str, Any]:
    """Translate resolved coding RNA/CDS and report amino-acid composition.

    Args:
        sequence: RNA (preferred) or DNA that can be viewed as coding-strand.

    Returns:
        translation, composition, coding-context method, limitations. Status
        COMPUTED, UNAVAILABLE, or ERROR. Does not imply arbitrary RNA is a gene.

    Raises:
        TypeError: Non-string input.
    """
    from modules.protein_analysis import amino_acid_composition, translate as translate_rna

    parsed = resolve_pasted_sequence(sequence)
    raw = str(parsed.get("sequence") or "")
    try:
        resolution = resolve_coding_codons(raw, source="auto", allow_whole_sequence=False)
    except ValueError as exc:
        return {
            "status": "UNAVAILABLE",
            "reason": str(exc),
            "coding_context": "none",
            "limitations": [
                "Arbitrary RNA is not assumed to be protein-coding. "
                "HelixScope translates concatenated resolved CDS/ORFs only."
            ],
            "translation": "",
            "composition": {},
        }
    codons = str(resolution.get("codons") or "")
    rna_codons = "".join(codons.split()).upper().replace("T", "U")
    if not rna_codons or "U" not in rna_codons:
        return {
            "status": "UNAVAILABLE",
            "reason": "Resolved coding nucleotides did not yield RNA codons.",
            "coding_context": str(resolution.get("method") or ""),
            "translation": "",
            "composition": {},
            "resolution": resolution,
        }
    try:
        protein = translate_rna(rna_codons)
        composition = amino_acid_composition(protein)
    except ValueError as exc:
        return {
            "status": "ERROR",
            "reason": str(exc),
            "coding_context": str(resolution.get("method") or ""),
            "translation": "",
            "composition": {},
            "resolution": resolution,
        }
    return {
        "status": "COMPUTED",
        "coding_context": str(resolution.get("method") or "prediction"),
        "method": str(resolution.get("method_label") or resolution.get("method") or ""),
        "frame": 1,
        "stop_semantics": "to_first_stop_of_concatenated_coding_codons",
        "n_cds": resolution.get("n_cds"),
        "n_codons": resolution.get("n_codons"),
        "translation": protein,
        "composition": composition,
        "resolution": {
            "method": resolution.get("method"),
            "method_label": resolution.get("method_label"),
            "n_cds": resolution.get("n_cds"),
            "n_codons": resolution.get("n_codons"),
            "warnings": resolution.get("warnings") or [],
        },
        "limitations": [
            "This summarises amino-acid usage of concatenated coding regions, "
            "not the sequence of any one experimental protein.",
            "ORF prediction is not gene annotation.",
        ],
        "input_format": parsed.get("format"),
        "input_identifier": parsed.get("identifier") or "",
    }
