"""Modelo imutavel de comparacao estrutural e parser do JSON RCSB.

Nao chama a rede. Nao aplica transformacoes 3D (ver structure_superposition).
Nao inventa RMSD nem TM-score: so extrai os scores tipados do JSON e declara
atomos (C-alpha), numero de pares e se o valor e do bloco ou do summary global.

Numeracao RCSB: beg_seq_id e 1-based sequential (label_seq_id mmCIF), nao
auth_seq_id. O mapeamento A<->B so usa essas correspondencias; nunca o residuo
espacialmente mais proximo.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Mapping, Optional, Sequence

from . import provenance, rcsb_alignment

ATOMS_USED: str = "CA"
"""Atomos de ajustamento declarados pela Alignment API (C-alpha de proteina)."""

SCORE_RMSD: str = "RMSD"
SCORE_TM: str = "TM-score"
SCORE_IDENTITY: str = "sequence-identity"
SCORE_SIMILARITY: str = "sequence-similarity"
SCORE_SIMILARITY_FATCAT: str = "similarity-score"


class StructureAlignmentError(RuntimeError):
    """JSON de alinhamento invalido ou incompleto.

    Attributes:
        category: PARSING_ERROR, INVALID_INPUT, UNAVAILABLE ou PARTIAL.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "PARSING_ERROR").strip().upper()


def scores_by_type(scores: Any) -> list[dict]:
    """Preserva cada score com o seu tipo. Nao funde TM-scores.

    Args:
        scores: Lista {value, type} da API.

    Returns:
        Lista de dicts value, type, unit.

    Raises:
        Nenhum. Entradas malformadas sao ignoradas.
    """
    rows: list[dict] = []
    if not isinstance(scores, Sequence) or isinstance(scores, (str, bytes)):
        return rows
    for item in scores:
        if not isinstance(item, Mapping):
            continue
        score_type = str(item.get("type") or "").strip()
        if not score_type:
            continue
        try:
            value = float(item.get("value"))
        except (TypeError, ValueError):
            continue
        unit = "angstrom" if score_type.upper() == SCORE_RMSD.upper() or score_type == SCORE_RMSD else ""
        if "TM" in score_type.upper():
            unit = "dimensionless_0_1"
        if score_type in {SCORE_IDENTITY, SCORE_SIMILARITY}:
            unit = "fraction_0_1"
        rows.append(
            {
                "type": score_type,
                "value": value,
                "unit": unit,
                "source": "RCSB Alignment API summary.scores",
            }
        )
    return rows


def score_of_type(scores: Sequence[Mapping[str, Any]], score_type: str) -> Optional[float]:
    """Primeiro valor do tipo pedido, ou None (N/A, nunca 0).

    Args:
        scores: Saida de scores_by_type.
        score_type: Tipo exacto (RMSD, TM-score, ...).

    Returns:
        float ou None.

    Raises:
        Nenhum.
    """
    wanted = str(score_type or "").strip()
    matches = [float(row["value"]) for row in scores if str(row.get("type")) == wanted]
    if not matches:
        return None
    return matches[0]


def all_scores_of_type(scores: Sequence[Mapping[str, Any]], prefix: str) -> list[dict]:
    """Todos os scores cujo tipo comeca por prefixo (ex. TM-score).

    Args:
        scores: Saida de scores_by_type.
        prefix: Prefixo case-sensitive ou "TM".

    Returns:
        Lista filtrada, nao agregada.

    Raises:
        Nenhum.
    """
    needle = str(prefix or "").strip().upper()
    return [dict(row) for row in scores if needle in str(row.get("type") or "").upper()]


def expand_region_pairs(block: Mapping[str, Any]) -> list[dict]:
    """Pares de residuos estruturalmente equivalentes de um bloco.

    Args:
        block: Objeto structure_alignment[i].

    Returns:
        Lista de dicts com label_seq_id_ref, label_seq_id_tgt, chain ids,
        block_index nao incluido (o caller adiciona).

    Raises:
        StructureAlignmentError: PARSING_ERROR se as regioes nao tiverem o
            mesmo comprimento alinhado.

    Nota biologica:
        Cada posicao i do bloco alinha o residuo sequencial (1-based) da
        referencia com o da estrutura alvo. Gaps de sequencia nao entram
        nestas regioes estruturais; usam sequence_alignment.gaps.
    """
    regions = block.get("regions")
    if not isinstance(regions, Sequence) or len(regions) != 2:
        raise StructureAlignmentError(
            "structure_alignment block must have regions for exactly two structures.",
            "PARSING_ERROR",
        )
    ref_ranges = _range_list(regions[0])
    tgt_ranges = _range_list(regions[1])
    ref_ids = _expand_ranges(ref_ranges)
    tgt_ids = _expand_ranges(tgt_ranges)
    if len(ref_ids) != len(tgt_ids):
        raise StructureAlignmentError(
            "Aligned region lengths differ between reference and target. "
            "HelixScope refuses to invent missing residue pairs.",
            "PARSING_ERROR",
        )
    pairs: list[dict] = []
    for (ref_chain, ref_seq), (tgt_chain, tgt_seq) in zip(ref_ids, tgt_ids):
        pairs.append(
            {
                "reference_asym_id": ref_chain,
                "target_asym_id": tgt_chain,
                "reference_label_seq_id": ref_seq,
                "target_label_seq_id": tgt_seq,
                "numbering": "mmCIF label_seq_id (1-based sequential), not auth_seq_id",
                "atoms": ATOMS_USED,
            }
        )
    return pairs


def parse_alignment_payload(
    payload: Mapping[str, Any],
    *,
    method: str = "",
    method_kind: str = "",
    query: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Converte o JSON COMPLETE da RCSB num objeto de comparacao.

    Args:
        payload: Corpo /results com status COMPLETE.
        method: Nome de API pedido (se o meta nao o trouxer).
        method_kind: rigid ou flexible.
        query: Query original, para proveniencia.

    Returns:
        Dict imutavel (copia profunda) com estruturas, blocos, pares, scores
        globais e de bloco, coverage, status.

    Raises:
        StructureAlignmentError: PARSING_ERROR / INVALID_INPUT.

    Nota biologica:
        RMSD do summary global e RMSD de cada bloco sao valores distintos.
        TM-score, quando presente, e o valor da engine; HelixScope nao o
        renormaliza. Coverage e a fracao de residuos modelados alinhados
        (aln_coverage), nao identidade de sequencia.
    """
    if not isinstance(payload, Mapping):
        raise StructureAlignmentError("Alignment payload is not a JSON object.", "PARSING_ERROR")
    info = payload.get("info") if isinstance(payload.get("info"), Mapping) else {}
    status = str(info.get("status") or "").upper()
    if status and status != "COMPLETE":
        raise StructureAlignmentError(
            f"Alignment payload status is {status}, not COMPLETE.",
            "INVALID_INPUT",
        )
    meta = payload.get("meta") if isinstance(payload.get("meta"), Mapping) else {}
    api_method = str(meta.get("alignment_method") or method or "").strip()
    mode = str(meta.get("alignment_mode") or rcsb_alignment.MODE_PAIRWISE)
    results = payload.get("results")
    if not isinstance(results, Sequence) or not results:
        raise StructureAlignmentError("Alignment payload has no results array.", "PARSING_ERROR")
    first = results[0]
    if not isinstance(first, Mapping):
        raise StructureAlignmentError("Alignment result is not an object.", "PARSING_ERROR")
    structures = first.get("structures")
    if not isinstance(structures, Sequence) or len(structures) != 2:
        raise StructureAlignmentError(
            "Pairwise alignment must list exactly two structures (reference, target).",
            "PARSING_ERROR",
        )
    reference = _structure_ref(structures[0], role="reference")
    target = _structure_ref(structures[1], role="target")
    blocks_raw = first.get("structure_alignment")
    if not isinstance(blocks_raw, Sequence) or not blocks_raw:
        raise StructureAlignmentError("structure_alignment is missing.", "PARSING_ERROR")
    kind = str(method_kind or "").strip() or (
        "flexible"
        if api_method in rcsb_alignment.FLEXIBLE_METHODS
        else "rigid"
        if api_method in rcsb_alignment.RIGID_METHODS
        else "unknown"
    )
    if kind == "rigid" and len(blocks_raw) != 1:
        # still parse, but mark PARTIAL rather than pretending one transform is the whole
        kind_note = (
            "Rigid methods are documented to return one block. This payload has "
            f"{len(blocks_raw)} blocks; HelixScope does not collapse them."
        )
    else:
        kind_note = ""
    blocks: list[dict] = []
    all_pairs: list[dict] = []
    for index, raw_block in enumerate(blocks_raw):
        if not isinstance(raw_block, Mapping):
            raise StructureAlignmentError("structure_alignment block is not an object.", "PARSING_ERROR")
        pairs = expand_region_pairs(raw_block)
        for pair in pairs:
            pair["block_index"] = index
        transforms = raw_block.get("transformations")
        if not isinstance(transforms, Sequence) or len(transforms) != 2:
            raise StructureAlignmentError(
                "Each block must contain two 4x4 transformations (reference, target).",
                "PARSING_ERROR",
            )
        block_summary = raw_block.get("summary") if isinstance(raw_block.get("summary"), Mapping) else {}
        block_scores = scores_by_type(block_summary.get("scores"))
        block = {
            "block_index": index,
            "n_pairs": len(pairs),
            "pairs": pairs,
            "reference_transform_column_major_4x4": _as_float_list(transforms[0]),
            "target_transform_column_major_4x4": _as_float_list(transforms[1]),
            "transform_convention": (
                "4x4 column-major, index j*4+i, as documented by RCSB Alignment API"
            ),
            "scores": block_scores,
            "rmsd_angstrom": score_of_type(block_scores, SCORE_RMSD),
            "rmsd_scope": "aligned_CA_pairs_of_this_block",
            "n_aln_residue_pairs": _optional_int(block_summary.get("n_aln_residue_pairs")),
        }
        blocks.append(block)
        all_pairs.extend(pairs)
    global_summary = first.get("summary") if isinstance(first.get("summary"), Mapping) else {}
    global_scores = scores_by_type(global_summary.get("scores"))
    tm_scores = all_scores_of_type(global_scores, "TM")
    coverage = global_summary.get("aln_coverage")
    coverage_list = (
        [float(item) for item in coverage]
        if isinstance(coverage, Sequence) and not isinstance(coverage, (str, bytes))
        else []
    )
    modeled = global_summary.get("n_modeled_residues")
    modeled_list = (
        [int(item) for item in modeled]
        if isinstance(modeled, Sequence) and not isinstance(modeled, (str, bytes))
        else []
    )
    n_pairs = int(global_summary.get("n_aln_residue_pairs") or len(all_pairs))
    visual_status = "AVAILABLE"
    if kind == "flexible" and len(blocks) > 1:
        visual_status = "PARTIAL"
    mapping_ab = {
        (row["reference_asym_id"], row["reference_label_seq_id"]): (
            row["target_asym_id"],
            row["target_label_seq_id"],
        )
        for row in all_pairs
    }
    mapping_ba = {
        (row["target_asym_id"], row["target_label_seq_id"]): (
            row["reference_asym_id"],
            row["reference_label_seq_id"],
        )
        for row in all_pairs
    }
    identity_hash = hashlib.sha256(
        json.dumps(
            {
                "reference": reference,
                "target": target,
                "method": api_method,
                "n_pairs": n_pairs,
                "n_blocks": len(blocks),
            },
            sort_keys=True,
            default=str,
        ).encode("utf-8")
    ).hexdigest()
    record = {
        "kind": "structural_comparison",
        "status": "RETRIEVED",
        "evidence_status": "RETRIEVED",
        "source": "RCSB PDB Structure Alignment API",
        "docs": rcsb_alignment.DOCS_URL,
        "engine": "RCSB PDB Structure Alignment API",
        "engine_location": "remote",
        "engine_version": "api/v1",
        "method": api_method,
        "method_display": (
            rcsb_alignment.display_name_for(api_method)
            if api_method in {row["api_name"] for row in rcsb_alignment.METHOD_CATALOG}
            else api_method
        ),
        "method_kind": kind,
        "method_kind_note": kind_note,
        "alignment_mode": mode,
        "atoms_fitted": ATOMS_USED,
        "atoms_fitted_note": (
            "RCSB Alignment API fits backbone C-alpha atoms of the first model, "
            "first altloc, first microheterogeneous residue. Minimum 10 CA residues "
            "per structure (API documentation, 2026-08-28)."
        ),
        "reference": reference,
        "target": target,
        "n_blocks": len(blocks),
        "blocks": blocks,
        "residue_pairs": all_pairs,
        "n_aligned_residue_pairs": n_pairs,
        "mapping_reference_to_target": {
            f"{chain}:{seq}": f"{other[0]}:{other[1]}"
            for (chain, seq), other in mapping_ab.items()
        },
        "mapping_target_to_reference": {
            f"{chain}:{seq}": f"{other[0]}:{other[1]}"
            for (chain, seq), other in mapping_ba.items()
        },
        "unmapped_policy": (
            "A residue absent from mapping_* is UNMAPPED. HelixScope does not "
            "assign the spatially nearest residue."
        ),
        "global_scores": global_scores,
        "rmsd_global_angstrom": score_of_type(global_scores, SCORE_RMSD),
        "rmsd_block0_angstrom": blocks[0]["rmsd_angstrom"] if blocks else None,
        "rmsd_note": (
            "RMSD is over structurally equivalent C-alpha pairs, not all atoms "
            "and not unaligned residues. Block RMSD and global summary RMSD are "
            "kept separate and are not averaged."
        ),
        "tm_scores": tm_scores,
        "tm_score_note": (
            "TM-score values are those returned by the engine. Distinct "
            "normalizations are not fused. Absence is N/A, not 0."
        ),
        "sequence_identity": score_of_type(global_scores, SCORE_IDENTITY),
        "sequence_similarity": score_of_type(global_scores, SCORE_SIMILARITY),
        "aln_coverage_percent": coverage_list,
        "coverage_note": (
            "aln_coverage is percent of modeled residues included in the "
            "alignment, one value per structure [reference, target]. An AlphaFold "
            "full-length model vs a PDB fragment will show unequal coverage; RMSD "
            "without this field is misleading."
        ),
        "n_modeled_residues": modeled_list,
        "seq_aln_len": _optional_int(global_summary.get("seq_aln_len")),
        "sequence_alignment": first.get("sequence_alignment"),
        "superposition_visual": visual_status,
        "superposition_visual_note": (
            "Flexible alignments with multiple blocks cannot be shown as a single "
            "rigid transform. The table of blocks remains AVAILABLE; the 3D overlay "
            "is PARTIAL (first block only) or UNAVAILABLE."
            if visual_status == "PARTIAL"
            else "Rigid alignments have one transform applied to a coordinate copy."
        ),
        "ticket": str(info.get("uuid") or ""),
        "query": copy.deepcopy(query) if query else None,
        "identity_hash": identity_hash,
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "original_coordinates_preserved": True,
        "disclaimer": (
            "Structural similarity is not homology, not function, and not a "
            "pathogenicity claim. Experimental PDB and predicted AlphaFold models "
            "remain labelled as such."
        ),
    }
    return copy.deepcopy(record)


def lookup_correspondence(
    record: Mapping[str, Any],
    *,
    side: str,
    asym_id: str,
    label_seq_id: int,
) -> dict:
    """A -> B ou B -> A. Sem vizinho espacial.

    Args:
        record: Saida de parse_alignment_payload.
        side: reference ou target (o residuo clicado).
        asym_id: Cadeia.
        label_seq_id: label_seq_id 1-based.

    Returns:
        Dict status AVAILABLE ou UNMAPPED, e o par se existir.

    Raises:
        StructureAlignmentError: INVALID_INPUT.
    """
    role = str(side or "").strip().lower()
    if role not in {"reference", "target"}:
        raise StructureAlignmentError("side must be reference or target.", "INVALID_INPUT")
    key = f"{asym_id}:{int(label_seq_id)}"
    if role == "reference":
        mapped = (record.get("mapping_reference_to_target") or {}).get(key)
    else:
        mapped = (record.get("mapping_target_to_reference") or {}).get(key)
    if not mapped:
        return {
            "status": "UNMAPPED",
            "clicked_side": role,
            "clicked": key,
            "corresponding": None,
            "reason": (
                "This residue is not in the structural correspondence returned by "
                "the alignment engine. It is UNMAPPED, not assigned to the nearest "
                "C-alpha."
            ),
        }
    other_side = "target" if role == "reference" else "reference"
    chain, _, seq = str(mapped).partition(":")
    return {
        "status": "AVAILABLE",
        "clicked_side": role,
        "clicked": key,
        "corresponding_side": other_side,
        "corresponding": mapped,
        "corresponding_asym_id": chain,
        "corresponding_label_seq_id": int(seq) if seq else None,
    }


def comparison_summary_rows(record: Mapping[str, Any]) -> list[tuple[str, str]]:
    """Linhas para a grelha da UI.

    Args:
        record: Objeto de comparacao.

    Returns:
        Lista (label, value).

    Raises:
        Nenhum.
    """
    ref = record.get("reference") or {}
    tgt = record.get("target") or {}
    tm = record.get("tm_scores") or []
    tm_text = (
        ", ".join(f"{item.get('type')}={item.get('value')}" for item in tm)
        if tm
        else "N/A"
    )
    coverage = record.get("aln_coverage_percent") or []
    coverage_text = (
        f"reference {coverage[0]}%; target {coverage[1]}%"
        if len(coverage) >= 2
        else "N/A"
    )
    rmsd_g = record.get("rmsd_global_angstrom")
    rmsd_b = record.get("rmsd_block0_angstrom")
    identity = record.get("sequence_identity")
    return [
        ("Reference", f"{ref.get('entry_id')} chain {ref.get('asym_id')}"),
        ("Target", f"{tgt.get('entry_id')} chain {tgt.get('asym_id')}"),
        ("Method", str(record.get("method_display") or record.get("method"))),
        ("Method kind", str(record.get("method_kind"))),
        ("Atoms fitted", str(record.get("atoms_fitted"))),
        ("Aligned residue pairs", str(record.get("n_aligned_residue_pairs"))),
        ("Blocks", str(record.get("n_blocks"))),
        ("Coverage", coverage_text),
        ("RMSD global (CA, aligned)", "N/A" if rmsd_g is None else f"{rmsd_g} A"),
        ("RMSD block 0 (CA, aligned)", "N/A" if rmsd_b is None else f"{rmsd_b} A"),
        ("TM-score(s)", tm_text),
        ("Sequence identity", "N/A" if identity is None else str(identity)),
        ("Superposition visual", str(record.get("superposition_visual"))),
        ("Status", str(record.get("status"))),
        ("Engine", str(record.get("engine"))),
    ]


def _structure_ref(item: Mapping[str, Any], *, role: str) -> dict:
    selection = item.get("selection") if isinstance(item.get("selection"), Mapping) else {}
    entry = str(item.get("entry_id") or "")
    kind = "predicted" if entry.upper().startswith(("AF_", "MA_")) else "experimental"
    return {
        "role": role,
        "entry_id": entry,
        "asym_id": str(selection.get("asym_id") or ""),
        "kind": kind,
        "kind_label": "PREDICTED" if kind == "predicted" else "EXPERIMENTAL",
        "source": "RCSB PDB / RCSB CSM" if kind == "predicted" else "RCSB PDB",
    }


def _range_list(side: Any) -> list[dict]:
    if isinstance(side, Mapping):
        return [side]
    if isinstance(side, Sequence) and not isinstance(side, (str, bytes)):
        return [item for item in side if isinstance(item, Mapping)]
    raise StructureAlignmentError("Region list is malformed.", "PARSING_ERROR")


def _expand_ranges(ranges: Sequence[Mapping[str, Any]]) -> list[tuple[str, int]]:
    residues: list[tuple[str, int]] = []
    for item in ranges:
        chain = str(item.get("asym_id") or "")
        start = item.get("beg_seq_id")
        length = item.get("length")
        if start is None or length is None:
            raise StructureAlignmentError("Region lacks beg_seq_id or length.", "PARSING_ERROR")
        start_i = int(start)
        length_i = int(length)
        if length_i < 0:
            raise StructureAlignmentError("Region length cannot be negative.", "PARSING_ERROR")
        for offset in range(length_i):
            residues.append((chain, start_i + offset))
    return residues


def _as_float_list(values: Any) -> list[float]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise StructureAlignmentError("Transformation is not a numeric list.", "PARSING_ERROR")
    if len(values) != 16:
        raise StructureAlignmentError(
            "Transformation must be a 4x4 matrix flattened to 16 numbers.",
            "PARSING_ERROR",
        )
    try:
        return [float(item) for item in values]
    except (TypeError, ValueError) as exc:
        raise StructureAlignmentError("Transformation contains a non-numeric entry.", "PARSING_ERROR") from exc


def _optional_int(value: object) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
