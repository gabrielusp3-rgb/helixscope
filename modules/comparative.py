"""Comparacoes cientificas: estruturas, variantes, guias, evolucao.

Reutiliza alinhamento de sequencias (alignment/msa) e Variant Explorer.
Nao inventa patogenicidade, dN/dS, genome-wide MIT, nem estruturas mutantes.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import copy
import tempfile
from typing import Any, Mapping, Optional, Sequence

from . import (
    alignment,
    evidence_workspace,
    msa as msa_module,
    provenance,
    protein_structure,
    rcsb_alignment,
    structure_alignment,
    structure_superposition,
    usalign,
)

COMPARISON_MODES: tuple[str, ...] = (
    "structures",
    "variants",
    "proteins",
    "guides",
    "evolution",
)


class ComparativeError(RuntimeError):
    """Comparacao recusada ou incompleta.

    Attributes:
        category: INVALID_INPUT, UNAVAILABLE, SCOPE_MISMATCH ou ERROR.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper()


def compare_structures_remote(
    *,
    reference_entry: str,
    target_entry: str,
    reference_chain: str = "A",
    target_chain: str = "A",
    method: str = rcsb_alignment.METHOD_TMALIGN,
    reference_parsed: Optional[Mapping[str, Any]] = None,
    target_parsed: Optional[Mapping[str, Any]] = None,
    urlopen_fn=None,
    sleep_fn=None,
) -> dict:
    """Alinha via RCSB API e, se houver coordenadas locais, valida a matriz.

    Args:
        reference_entry: PDB/CSM referencia.
        target_entry: PDB/CSM alvo.
        reference_chain: Cadeia referencia.
        target_chain: Cadeia alvo.
        method: Nome de API RCSB.
        reference_parsed: Estrutura parseada opcional (imutavel).
        target_parsed: Estrutura parseada opcional.
        urlopen_fn: HTTP injetavel.
        sleep_fn: sleep injetavel.

    Returns:
        Dict alignment, superposition (ou None), evidence, export.

    Raises:
        rcsb_alignment.AlignmentApiError, structure_alignment.StructureAlignmentError.

    Nota biologica:
        PDB experimental vs AlphaFold previsto: os kind_label permanecem
        distintos. Coverage desigual e esperada quando o modelo AF e full-length
        e o PDB e um fragmento.
    """
    raw_wrap = rcsb_alignment.align_pairwise(
        reference_entry=reference_entry,
        target_entry=target_entry,
        reference_chain=reference_chain,
        target_chain=target_chain,
        method=method,
        urlopen_fn=urlopen_fn,
        sleep_fn=sleep_fn,
    )
    alignment = structure_alignment.parse_alignment_payload(
        raw_wrap["raw"],
        method=raw_wrap["method"],
        method_kind=raw_wrap["method_kind"],
        query=raw_wrap.get("query"),
    )
    superposition = None
    rmsd_check = None
    if reference_parsed is not None and target_parsed is not None:
        superposition = structure_superposition.superposition_bundle(
            alignment=alignment,
            reference_parsed=reference_parsed,
            target_parsed=target_parsed,
            block_index=0,
        )
        rmsd_check = structure_superposition.rmsd_agrees(
            superposition["independent_rmsd"],
            alignment.get("rmsd_global_angstrom"),
        )
    items = evidence_workspace.collect_from_alignment(alignment)
    export = evidence_workspace.export_pack(
        kind="structures",
        inputs={
            "reference_entry": reference_entry,
            "target_entry": target_entry,
            "reference_chain": reference_chain,
            "target_chain": target_chain,
            "method": method,
            "backend": "RCSB Alignment API",
        },
        items=items,
        extra={
            "identity_hash": alignment.get("identity_hash"),
            "ticket": alignment.get("ticket"),
            "rmsd_agreement": rmsd_check,
        },
    )
    return {
        "mode": "structures",
        "alignment": alignment,
        "superposition": superposition,
        "rmsd_agreement": rmsd_check,
        "evidence": items,
        "export": export,
        "network": raw_wrap.get("network"),
        "engine_location": "remote",
    }


def compare_structures_usalign(
    *,
    reference_parsed: Mapping[str, Any],
    target_parsed: Mapping[str, Any],
    reference_entry: str = "",
    target_entry: str = "",
    reference_chain: str = "A",
    target_chain: str = "A",
    mol: str = "prot",
    record_validation: bool = False,
) -> dict:
    """Alinha duas estruturas ja carregadas com US-align local.

    Args:
        reference_parsed: Envelope ou parse com raw_structure_text.
        target_parsed: Idem para o alvo.
        reference_entry: Identificador de exibicao (PDB).
        target_entry: Identificador de exibicao.
        reference_chain: -chain1.
        target_chain: -chain2.
        mol: auto, prot ou RNA.
        record_validation: Grava LIVE_VALIDATED se o parse passar.

    Returns:
        Dict alignment, superposition None, evidence, engine_location local.

    Raises:
        usalign.USAlignError, ComparativeError.

    Nota biologica:
        Este caminho nao chama a RCSB Alignment API. Overlay 3D RCSB nao e
        reutilizado: US-align roda a estrutura 1 sobre a 2.
    """
    text1 = str(
        reference_parsed.get("raw_structure_text")
        or reference_parsed.get("raw_text")
        or ""
    )
    text2 = str(
        target_parsed.get("raw_structure_text")
        or target_parsed.get("raw_text")
        or ""
    )
    if not text1.strip() or not text2.strip():
        raise ComparativeError(
            "US-align needs loaded coordinate text for both structures. "
            "Download coordinates first. This is not an RCSB API alignment.",
            "UNAVAILABLE",
        )
    with tempfile.TemporaryDirectory(prefix="helixscope_usalign_cmp_") as tmp:
        file1 = usalign.write_structure_text(
            text1,
            tmp,
            f"{str(reference_entry or 'ref').replace('/', '_')}.cif",
        )
        file2 = usalign.write_structure_text(
            text2,
            tmp,
            f"{str(target_entry or 'tgt').replace('/', '_')}.cif",
        )
        alignment = usalign.align_structure_files(
            file1,
            file2,
            mol=mol,
            chain1=reference_chain,
            chain2=target_chain,
            record_validation=record_validation,
        )
    if reference_entry:
        alignment["reference"]["entry_id"] = reference_entry
        alignment["reference"]["kind_label"] = str(
            reference_parsed.get("kind") or reference_parsed.get("kind_label") or "EXPERIMENTAL"
        ).upper()
        if str(alignment["reference"]["kind_label"]).lower() == "predicted":
            alignment["reference"]["kind_label"] = "PREDICTED"
        if str(alignment["reference"]["kind_label"]).lower() == "illustrative":
            alignment["reference"]["kind_label"] = "ILLUSTRATIVE"
    if target_entry:
        alignment["target"]["entry_id"] = target_entry
        alignment["target"]["kind_label"] = str(
            target_parsed.get("kind") or target_parsed.get("kind_label") or "EXPERIMENTAL"
        ).upper()
        if str(alignment["target"]["kind_label"]).lower() == "predicted":
            alignment["target"]["kind_label"] = "PREDICTED"
        if str(alignment["target"]["kind_label"]).lower() == "illustrative":
            alignment["target"]["kind_label"] = "ILLUSTRATIVE"
    items = evidence_workspace.collect_from_alignment(alignment)
    export = evidence_workspace.export_pack(
        kind="structures",
        inputs={
            "reference_entry": reference_entry,
            "target_entry": target_entry,
            "reference_chain": reference_chain,
            "target_chain": target_chain,
            "method": "US-align",
            "backend": "US-align",
        },
        items=items,
        extra={"identity_hash": alignment.get("identity_hash")},
    )
    return {
        "mode": "structures",
        "alignment": alignment,
        "superposition": None,
        "rmsd_agreement": None,
        "evidence": items,
        "export": export,
        "network": None,
        "engine_location": "local",
        "backend": "US-align",
    }


def compare_variants(
    result_a: Mapping[str, Any],
    result_b: Mapping[str, Any],
) -> dict:
    """Compara duas saidas do Variant Explorer. Sem ranking de patogenicidade.

    Args:
        result_a: explore_variant A.
        result_b: explore_variant B.

    Returns:
        Dict side_by_side, property_delta (se ambos tiverem aa), evidence.

    Raises:
        ComparativeError: INVALID_INPUT se faltar variant.

    Nota biologica:
        Duas variantes no mesmo gene nao sao ordenadas como "pior". ClinVar
        divergente e mostrado. Sem coordenadas mutantes.
    """
    va = result_a.get("variant")
    vb = result_b.get("variant")
    if not isinstance(va, Mapping) or not isinstance(vb, Mapping):
        raise ComparativeError("Both sides must be Variant Explorer results.", "INVALID_INPUT")
    layers_a = result_a.get("layers") or {}
    layers_b = result_b.get("layers") or {}
    delta = None
    res_a = (layers_a.get("residue_properties") or {}).get("data")
    res_b = (layers_b.get("residue_properties") or {}).get("data")
    if isinstance(res_a, Mapping) and isinstance(res_b, Mapping):
        delta = {
            "status": "COMPUTED",
            "note": (
                "Physico-chemical deltas already computed per variant against the "
                "reference amino acid. HelixScope does not emit ddG or a "
                "pathogenicity rank."
            ),
            "variant_a": {
                "reference_aa": res_a.get("reference_amino_acid"),
                "variant_aa": res_a.get("variant_amino_acid"),
                "hydropathy_difference": res_a.get("hydropathy_difference"),
                "charge_class_changed": res_a.get("charge_class_changed"),
                "volume_difference_a3": res_a.get("volume_difference_a3"),
                "stability_claim": res_a.get("stability_claim"),
            },
            "variant_b": {
                "reference_aa": res_b.get("reference_amino_acid"),
                "variant_aa": res_b.get("variant_amino_acid"),
                "hydropathy_difference": res_b.get("hydropathy_difference"),
                "charge_class_changed": res_b.get("charge_class_changed"),
                "volume_difference_a3": res_b.get("volume_difference_a3"),
                "stability_claim": res_b.get("stability_claim"),
            },
        }
    items = evidence_workspace.collect_from_variant_result(result_a, label="Variant A")
    items.extend(evidence_workspace.collect_from_variant_result(result_b, label="Variant B"))
    export = evidence_workspace.export_pack(
        kind="variants",
        inputs={
            "a": va.get("identity_hash"),
            "b": vb.get("identity_hash"),
            "assembly_a": va.get("assembly"),
            "assembly_b": vb.get("assembly"),
        },
        items=items,
    )
    return {
        "mode": "variants",
        "ranking": None,
        "ranking_note": "HelixScope does not rank variants by pathogenicity.",
        "side_by_side": {
            "a": _variant_card(result_a),
            "b": _variant_card(result_b),
        },
        "property_delta": delta,
        "mutant_structure": {
            "status": "UNAVAILABLE",
            "reason": (
                "No mutant atomic model is generated. An experimental or predicted "
                "structure explicitly corresponding to the substitution may be "
                "compared in Structures mode if the user supplies that PDB/AF id."
            ),
        },
        "evidence": items,
        "conflicts": export["conflicts"],
        "export": export,
    }


def compare_guides(
    guide_a: Mapping[str, Any],
    guide_b: Mapping[str, Any],
    *,
    search_a: Optional[Mapping[str, Any]] = None,
    search_b: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Compara dois guias. MIT/CFD so se o ambito for o mesmo.

    Args:
        guide_a: Dict de guia (guide_sequence, pam, ...).
        guide_b: Dict de guia.
        search_a: Envelope de busca off-target opcional.
        search_b: Envelope de busca off-target opcional.

    Returns:
        Dict comparable_specificity, side_by_side, warnings.

    Raises:
        ComparativeError: INVALID_INPUT.
    """
    seq_a = str(guide_a.get("guide_sequence") or guide_a.get("sequence") or "").strip().upper()
    seq_b = str(guide_b.get("guide_sequence") or guide_b.get("sequence") or "").strip().upper()
    if not seq_a or not seq_b:
        raise ComparativeError("Both guides need a guide_sequence.", "INVALID_INPUT")
    spec_ok, spec_reason = _specificity_comparable(search_a, search_b)
    return {
        "mode": "guides",
        "side_by_side": {
            "a": _guide_card(guide_a, search_a),
            "b": _guide_card(guide_b, search_b),
        },
        "specificity_comparable": spec_ok,
        "specificity_reason": spec_reason,
        "genome_wide_comparison": False,
        "genome_wide_note": (
            "Genome-wide MIT/CFD comparison requires COMPLETED_FULL_REFERENCE on "
            "the same public assembly. TEST REFERENCE comparisons are TEST-ONLY."
        ),
        "export": evidence_workspace.export_pack(
            kind="guides",
            inputs={"a": seq_a, "b": seq_b},
            items=[
                evidence_workspace.evidence_item(
                    field="guide_a",
                    value=seq_a,
                    source="crispr",
                    evidence_status="COMPUTED",
                ),
                evidence_workspace.evidence_item(
                    field="guide_b",
                    value=seq_b,
                    source="crispr",
                    evidence_status="COMPUTED",
                ),
            ],
        ),
    }


def msa_column_to_structure(
    *,
    msa_result: Mapping[str, Any],
    column: int,
    member_id: str,
    parsed_structure: Mapping[str, Any],
    chain: str,
) -> dict:
    """Coluna de MSA -> residuo da estrutura se a sequencia sem gaps casar.

    Args:
        msa_result: MSA COMPLETED.
        column: Indice 0-based da coluna.
        member_id: Identificador da linha.
        parsed_structure: Estrutura parseada.
        chain: label_asym_id.

    Returns:
        Dict status AVAILABLE/UNMAPPED, query residue, CA se existir.

    Raises:
        ComparativeError: INVALID_INPUT.

    Nota biologica:
        Conservacao da coluna nao e essencialidade funcional.
    """
    try:
        info = msa_module.column_detail(msa_result, int(column))
    except Exception as exc:  # noqa: BLE001 - MSA errors become UNMAPPED/INVALID
        raise ComparativeError(str(exc), "INVALID_INPUT") from exc
    scores = list((msa_result.get("conservation") or {}).get("scores") or [])
    conservation = scores[int(column)] if 0 <= int(column) < len(scores) else None
    residue_in_row = None
    ungapped_index = 0
    aligned_rows = msa_result.get("aligned_rows") or msa_result.get("sequences") or []
    # fallback: members + raw alignment is complex; use inspect if it has per-id
    for member in msa_result.get("rows") or msa_result.get("members") or []:
        if str(member.get("identifier") or member.get("id") or "") != str(member_id):
            continue
        aligned = str(member.get("aligned") or member.get("sequence_aligned") or "")
        if aligned and 0 <= int(column) < len(aligned):
            symbol = aligned[int(column)]
            if symbol not in {"-", "."}:
                residue_in_row = symbol.upper()
                ungapped_index = sum(1 for ch in aligned[: int(column)] if ch not in {"-", "."})
        break
    if residue_in_row is None:
        return {
            "status": "UNMAPPED",
            "column": int(column),
            "conservation_shannon": conservation,
            "conservation_note": (
                "MSA conservation is 1 - H/log2(A) among analyzed sequences. "
                "Highly conserved is not functionally essential."
            ),
            "reason": "This MSA row has a gap in the column or the member id was not found.",
            "column_info": info,
        }
    mapped = _structure_residue_at_ungapped(
        parsed_structure, chain=chain, ungapped_index=ungapped_index
    )
    if mapped is None:
        return {
            "status": "UNMAPPED",
            "column": int(column),
            "msa_residue": residue_in_row,
            "ungapped_index": ungapped_index,
            "conservation_shannon": conservation,
            "conservation_note": (
                "MSA conservation is 1 - H/log2(A) among analyzed sequences. "
                "Highly conserved is not functionally essential."
            ),
            "reason": "No C-alpha at this ungapped index on the requested chain.",
        }
    return {
        "status": "AVAILABLE",
        "column": int(column),
        "msa_residue": residue_in_row,
        "ungapped_index": ungapped_index,
        "conservation_shannon": conservation,
        "conservation_note": (
            "MSA conservation is 1 - H/log2(A) among analyzed sequences. "
            "Highly conserved is not functionally essential."
        ),
        "structure_residue": mapped,
        "legend": {
            "scale": "Shannon conservation 0-1 (gaps excluded from the alphabet count)",
            "not_functional_annotation": True,
        },
    }


def structure_residue_to_msa_column(
    *,
    msa_result: Mapping[str, Any],
    member_id: str,
    ungapped_index: int,
) -> dict:
    """Residuo (indice sem gaps) -> coluna do MSA.

    Args:
        msa_result: MSA COMPLETED.
        member_id: Linha.
        ungapped_index: 0-based na sequencia sem gaps.

    Returns:
        Dict status, column.

    Raises:
        ComparativeError: INVALID_INPUT.
    """
    for member in msa_result.get("rows") or msa_result.get("members") or []:
        if str(member.get("identifier") or member.get("id") or "") != str(member_id):
            continue
        aligned = str(member.get("aligned") or member.get("sequence_aligned") or "")
        seen = -1
        for col, symbol in enumerate(aligned):
            if symbol in {"-", "."}:
                continue
            seen += 1
            if seen == int(ungapped_index):
                return {
                    "status": "AVAILABLE",
                    "column": col,
                    "ungapped_index": int(ungapped_index),
                    "residue": symbol.upper(),
                    "member_id": member_id,
                }
        return {
            "status": "UNMAPPED",
            "reason": "ungapped_index is past the end of this MSA row.",
            "member_id": member_id,
        }
    raise ComparativeError(f"MSA member {member_id!r} was not found.", "INVALID_INPUT")


def group_associated_positions(
    *,
    msa_result: Mapping[str, Any],
    group_ids: Sequence[str],
    min_group_identity: float = 1.0,
) -> dict:
    """Posicoes onde o grupo escolhido e uniforme e o resto difere.

    Args:
        msa_result: MSA COMPLETED.
        group_ids: Identificadores das folhas/linhas do grupo.
        min_group_identity: Fracao minima de identidade dentro do grupo (1.0 =
            todos iguais, sem gap).

    Returns:
        Dict positions, statement (group-associated residue pattern).

    Raises:
        ComparativeError: INVALID_INPUT.

    Nota biologica:
        Nao e mutacao adaptativa, nao e dN/dS, nao e selecao positiva.
    """
    wanted = {str(item) for item in group_ids if str(item).strip()}
    if len(wanted) < 1:
        raise ComparativeError("Select at least one MSA member for the group.", "INVALID_INPUT")
    members = list(msa_result.get("rows") or msa_result.get("members") or [])
    rows = []
    others = []
    for member in members:
        ident = str(member.get("identifier") or member.get("id") or "")
        aligned = str(member.get("aligned") or member.get("sequence_aligned") or "")
        if not aligned:
            continue
        if ident in wanted:
            rows.append((ident, aligned))
        else:
            others.append((ident, aligned))
    if not rows:
        raise ComparativeError("None of the group identifiers matched MSA members.", "INVALID_INPUT")
    length = len(rows[0][1])
    if any(len(aligned) != length for _ident, aligned in rows + others):
        raise ComparativeError("MSA rows have unequal aligned length.", "INVALID_INPUT")
    hits: list[dict] = []
    for col in range(length):
        group_chars = [aligned[col].upper() for _ident, aligned in rows]
        if any(ch in {"-", "."} for ch in group_chars):
            continue
        majority = group_chars[0]
        frac = sum(1 for ch in group_chars if ch == majority) / len(group_chars)
        if frac < float(min_group_identity):
            continue
        other_chars = [
            aligned[col].upper()
            for _ident, aligned in others
            if aligned[col] not in {"-", "."}
        ]
        if not other_chars:
            continue
        if all(ch == majority for ch in other_chars):
            continue
        hits.append(
            {
                "column": col,
                "group_residue": majority,
                "group_fraction": frac,
                "rest_residues": sorted({ch for ch in other_chars}),
                "label": "group-associated residue pattern",
            }
        )
    return {
        "status": "COMPUTED",
        "n_positions": len(hits),
        "positions": hits,
        "group_ids": sorted(wanted),
        "statement": (
            "These columns show a residue shared inside the selected group and "
            "a different residue set in the remaining rows. This is a "
            "group-associated residue pattern, not an adaptive mutation and not "
            "a test of positive selection. dN/dS is not computed."
        ),
        "not_positive_selection": True,
        "not_dn_ds": True,
    }


def _variant_card(result: Mapping[str, Any]) -> dict:
    variant = result.get("variant") or {}
    layers = result.get("layers") or {}
    cons = layers.get("consequences") or {}
    protein = layers.get("protein_mapping") or {}
    clinical = layers.get("clinical_evidence") or {}
    domains = layers.get("domains") or {}
    data_cons = cons.get("data") or {}
    data_prot = protein.get("data") or {}
    transcripts = data_prot.get("transcripts") or []
    primary = transcripts[0] if transcripts else {}
    return {
        "assembly": variant.get("assembly"),
        "contig": variant.get("contig"),
        "position_1based": variant.get("position_1based"),
        "ref": variant.get("ref"),
        "alt": variant.get("alt"),
        "identity_hash": variant.get("identity_hash"),
        "consequences_status": cons.get("status"),
        "most_severe_consequence": data_cons.get("most_severe_consequence"),
        "transcript_count": data_cons.get("transcript_count"),
        "protein_status": protein.get("status"),
        "hgvsp": primary.get("hgvsp"),
        "amino_acids": primary.get("amino_acids"),
        "protein_position": primary.get("protein_start"),
        "clinvar_status": (clinical.get("data") or {}).get("status") or clinical.get("status"),
        "domains_status": domains.get("status"),
        "helixscope_effect": None,
    }


def _guide_card(guide: Mapping[str, Any], search: Optional[Mapping[str, Any]]) -> dict:
    spec = (search or {}).get("guide_specificity") or {}
    mit = spec.get("mit_sguide") or {}
    cfd = spec.get("cfd_sguide") or {}
    return {
        "guide_sequence": guide.get("guide_sequence") or guide.get("sequence"),
        "pam": guide.get("pam") or guide.get("pam_sequence"),
        "system": guide.get("cas_system") or guide.get("system") or "SpCas9",
        "search_status": None if search is None else search.get("status"),
        "genome_wide": bool((search or {}).get("genome_wide")),
        "search_scope": (search or {}).get("search_completeness") or mit.get("search_scope"),
        "n_hits": len(search.get("hits") or []) if search else None,
        "mit_sguide": mit.get("score") if search else None,
        "mit_genome_wide": bool(mit.get("genome_wide")) if search else False,
        "cfd_sguide": cfd.get("score") if search else None,
        "cfd_genome_wide": bool(cfd.get("genome_wide")) if search else False,
        "test_only": str((search or {}).get("status") or "").endswith("TEST_REFERENCE")
        or str(search.get("assembly_id") if search else "") == "HELIXSCOPE_TEST_REF",
    }


def _specificity_comparable(
    search_a: Optional[Mapping[str, Any]],
    search_b: Optional[Mapping[str, Any]],
) -> tuple[bool, str]:
    if search_a is None or search_b is None:
        return False, "Specificity scores are only compared when both guides have a completed search envelope."
    if bool(search_a.get("genome_wide")) != bool(search_b.get("genome_wide")):
        return False, "One search is genome-wide and the other is not. MIT/CFD are not comparable."
    if bool(search_a.get("genome_wide")):
        return False, (
            "No COMPLETED_FULL_REFERENCE comparison is available on this machine. "
            "Genome-wide guide comparison is refused."
        )
    status_a = str(search_a.get("status") or "")
    status_b = str(search_b.get("status") or "")
    if status_a != status_b:
        return False, f"Search statuses differ ({status_a} vs {status_b})."
    hash_a = str((search_a.get("engine_input_manifest") or {}).get("declared_sha256") or "")
    hash_b = str((search_b.get("engine_input_manifest") or {}).get("declared_sha256") or "")
    if hash_a and hash_b and hash_a != hash_b:
        return False, "Reference FASTA SHA-256 differs; MIT/CFD are not on the same scope."
    return True, "Same search status and reference identity; scoped MIT/CFD may be compared."


def _structure_residue_at_ungapped(
    parsed: Mapping[str, Any],
    *,
    chain: str,
    ungapped_index: int,
) -> Optional[dict]:
    cas = []
    for atom in parsed.get("atoms") or []:
        if not isinstance(atom, Mapping):
            continue
        if str(atom.get("atom_name") or "").upper() != "CA":
            continue
        if str(atom.get("group") or "") != "ATOM":
            continue
        atom_chain = str(atom.get("label_asym_id") or atom.get("auth_asym_id") or "")
        if atom_chain != str(chain):
            continue
        seq = atom.get("label_seq_id")
        if seq is None:
            continue
        cas.append(atom)
    cas.sort(key=lambda item: int(item.get("label_seq_id") or 0))
    if ungapped_index < 0 or ungapped_index >= len(cas):
        return None
    atom = cas[ungapped_index]
    return {
        "chain": str(atom.get("label_asym_id") or ""),
        "label_seq_id": atom.get("label_seq_id"),
        "auth_seq_id": atom.get("auth_seq_id"),
        "comp_id": atom.get("comp_id"),
        "x": atom.get("x"),
        "y": atom.get("y"),
        "z": atom.get("z"),
    }


def load_parsed_structure_for_entry(
    entry_id: str,
    *,
    urlopen_fn=None,
    cache: Optional[dict] = None,
) -> dict:
    """Carrega coordenadas para superposicao sem mutar o objeto original.

    Args:
        entry_id: PDB de 4 caracteres, AF-P01308-F1, ou CSM AF_AF...Fn.
        urlopen_fn: HTTP injetavel.
        cache: Cache injetavel do protein_structure.

    Returns:
        Dict status AVAILABLE/UNAVAILABLE, parsed, kind_label.

    Raises:
        Nenhum. Falhas de rede/parser viram UNAVAILABLE.

    Nota biologica:
        PDB experimental e AlphaFold previsto permanecem etiquetados. Gene
        symbol sozinho nao e identidade de proteina.
    """
    try:
        canonical = rcsb_alignment.validate_entry_id(entry_id)
    except rcsb_alignment.AlignmentApiError as exc:
        return {
            "status": "INVALID_INPUT",
            "parsed": None,
            "kind_label": None,
            "reason": str(exc),
        }
    uniprot_frag = rcsb_alignment.parse_rcsb_csm_uniprot(canonical)
    try:
        if len(canonical) == 4:
            parsed = protein_structure.load_deposited_macromolecule(
                source="RCSB PDB",
                structure_id=canonical,
                urlopen_fn=urlopen_fn,
                cache=cache,
            )
            kind_label = "EXPERIMENTAL"
            source = "RCSB PDB"
        elif uniprot_frag:
            accession, _fragment = uniprot_frag
            meta = protein_structure.fetch_alphafold_prediction(
                accession, urlopen_fn=urlopen_fn
            )
            cif_url = str(meta.get("cif_url") or "")
            if not cif_url:
                return {
                    "status": "UNAVAILABLE",
                    "parsed": None,
                    "kind_label": "PREDICTED",
                    "entry_id": canonical,
                    "reason": (
                        "AlphaFold metadata has no cif_url. Alignment metrics "
                        "from the RCSB API may still be AVAILABLE."
                    ),
                }
            parsed = protein_structure.load_deposited_macromolecule(
                source="AlphaFold DB",
                structure_id=canonical,
                metadata={**meta, "kind": "predicted"},
                cif_url=cif_url,
                urlopen_fn=urlopen_fn,
                cache=cache,
            )
            parsed["kind"] = "predicted"
            kind_label = "PREDICTED"
            source = "AlphaFold DB"
        else:
            return {
                "status": "UNAVAILABLE",
                "parsed": None,
                "kind_label": (
                    "PREDICTED"
                    if canonical.upper().startswith(("AF_", "MA_"))
                    else "EXPERIMENTAL"
                ),
                "entry_id": canonical,
                "reason": (
                    "No allowlisted coordinate download is configured for this "
                    "identifier. Alignment table remains usable; superposition "
                    "is UNAVAILABLE."
                ),
            }
    except Exception as exc:  # noqa: BLE001 - coordinate fetch is optional
        return {
            "status": "UNAVAILABLE",
            "parsed": None,
            "kind_label": None,
            "entry_id": canonical,
            "reason": str(exc),
        }
    return {
        "status": "AVAILABLE",
        "parsed": parsed,
        "kind_label": kind_label,
        "source": source,
        "entry_id": canonical,
        "n_atoms": len(parsed.get("atoms") or []),
        "original_coordinates_preserved": True,
    }


def compare_proteins(
    *,
    sequence_a: str,
    sequence_b: str,
    identifier_a: str = "protein_a",
    identifier_b: str = "protein_b",
    domains_a: Optional[Mapping[str, Any]] = None,
    domains_b: Optional[Mapping[str, Any]] = None,
    structure_comparison: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Compara duas proteinas: identidade de sequencia + InterPro + estrutura.

    Args:
        sequence_a: Sequencia aminoacidica A.
        sequence_b: Sequencia aminoacidica B.
        identifier_a: Rotulo A.
        identifier_b: Rotulo B.
        domains_a: Resultado InterPro real ou None.
        domains_b: Resultado InterPro real ou None.
        structure_comparison: Saida de compare_structures_remote ou None.

    Returns:
        Dict identity, domains, structure, disclaimer.

    Raises:
        ComparativeError: INVALID_INPUT.

    Nota biologica:
        Dominios nao sao inferidos pelo nome da proteina. InterPro ausente e
        UNAVAILABLE, nao um dominio inventado. Identidade de sequencia nao e
        RMSD.
    """
    seq_a = str(sequence_a or "").strip().upper()
    seq_b = str(sequence_b or "").strip().upper()
    if not seq_a or not seq_b:
        raise ComparativeError("Both protein sequences are required.", "INVALID_INPUT")
    aligned = alignment.pairwise_global(seq_a, seq_b)

    def _domain_card(payload: Optional[Mapping[str, Any]], sequence: str) -> dict:
        if not isinstance(payload, Mapping):
            return {
                "status": "UNAVAILABLE",
                "reason": (
                    "No InterPro retrieval was supplied. Domains are not guessed "
                    "from the protein name."
                ),
                "entries": [],
            }
        entries = list(payload.get("entries") or [])
        return {
            "status": payload.get("status") or "UNAVAILABLE",
            "source": payload.get("source") or "InterPro (EMBL-EBI)",
            "uniprot_accession": payload.get("uniprot_accession"),
            "count": payload.get("count"),
            "sequence_hash": payload.get("sequence_hash")
            or provenance.sequence_digest(sequence),
            "entries": [
                {
                    "accession": row.get("accession"),
                    "name": row.get("name"),
                    "source_database": row.get("source_database") or row.get("database"),
                    "type": row.get("type"),
                    "ranges": row.get("ranges"),
                    "release": row.get("release") or row.get("database_version"),
                }
                for row in entries[:40]
                if isinstance(row, Mapping)
            ],
        }

    return {
        "mode": "proteins",
        "identity": {
            "status": "COMPUTED",
            "method": aligned.get("method"),
            "identity_pct": aligned.get("identity_pct"),
            "ungapped_identity_pct": aligned.get("ungapped_identity_pct"),
            "coverage_seq1_pct": aligned.get("coverage_seq1_pct"),
            "coverage_seq2_pct": aligned.get("coverage_seq2_pct"),
            "identity_definition": aligned.get("identity_definition"),
            "hash_a": provenance.sequence_digest(seq_a),
            "hash_b": provenance.sequence_digest(seq_b),
            "identifier_a": identifier_a,
            "identifier_b": identifier_b,
        },
        "domains": {
            "a": _domain_card(domains_a, seq_a),
            "b": _domain_card(domains_b, seq_b),
            "note": (
                "Domain lists are InterPro retrievals bound to a UniProt "
                "accession/sequence hash. A shared domain name is not proof "
                "that the same residues are homologous."
            ),
        },
        "structure": structure_comparison,
        "sasa": {
            "status": "UNAVAILABLE",
            "reason": (
                "SASA comparison requires loaded coordinate envelopes on both "
                "sides. It is not estimated from sequence."
            ),
        },
        "disclaimer": (
            "Sequence identity, InterPro domains and structural RMSD answer "
            "different questions. None of them is function or pathogenicity."
        ),
        "export": evidence_workspace.export_pack(
            kind="proteins",
            inputs={
                "hash_a": provenance.sequence_digest(seq_a),
                "hash_b": provenance.sequence_digest(seq_b),
            },
            items=[
                evidence_workspace.evidence_item(
                    field="ungapped_identity_pct",
                    value=aligned.get("ungapped_identity_pct"),
                    source="alignment.pairwise_global",
                    evidence_status="COMPUTED",
                )
            ],
        ),
    }


_ = copy
