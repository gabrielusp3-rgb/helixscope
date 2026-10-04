"""Selecao molecular compartilhada entre sequencia, MSA, NCBI, motivo e 3D.

Eventos carregam referencias (hash, indices, fonte). Nao recalculam ciencia e
nao inventam residuos. Um residuo de proteina nao e tratado como nucleotideo.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, List, Mapping, Optional, Sequence

from . import protein_structure, provenance

MOLECULE_KINDS: tuple[str, ...] = ("DNA", "RNA", "PROTEIN")
SELECTION_SOURCES: tuple[str, ...] = (
    "sequence",
    "structure",
    "msa",
    "motif",
    "ncbi",
    "protein_analysis",
    "blast",
    "tree",
)
GAP_MESSAGE: str = (
    "No structural residue mapped to this sequence at this alignment column."
)


class SelectionError(Exception):
    """Falha classificada de selecao. Nunca deve virar um residuo inventado.

    Attributes:
        category: INVALID_INPUT, UNMAPPED, STALE, UNAVAILABLE.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def molecules_are_compatible(left: object, right: object) -> bool:
    """Proteina nao e DNA/RNA. DNA e RNA tambem nao se substituem.

    Args:
        left: Tipo da selecao.
        right: Tipo do alvo.

    Returns:
        True somente se ambos forem o mesmo de MOLECULE_KINDS.
    """
    a = str(left or "").strip().upper()
    b = str(right or "").strip().upper()
    return a in MOLECULE_KINDS and a == b


def make_selection(
    *,
    molecule: str,
    sequence_hash: str,
    query_indices: Sequence[int],
    source: str,
    label: str = "",
    structure_id: object = None,
    structure_hash: object = None,
    chain_id: object = None,
    model: object = None,
    msa_column: object = None,
    feature_label: object = None,
    atom_name: object = None,
    message: str = "",
    status: str = "READY",
    extra: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Constroi uma MoleculeSelection. Indices vazios nao sao preenchidos.

    Args:
        molecule: DNA, RNA ou PROTEIN.
        sequence_hash: Hash da sequencia de analise.
        query_indices: Posicoes 0-based na sequencia de analise.
        source: Um de SELECTION_SOURCES.
        label: Rotulo visual.
        structure_id: Identificador estrutural, se a selecao for estrutural.
        structure_hash: Hash do ficheiro estrutural, se conhecido.
        chain_id: Chain, se conhecida.
        model: Modelo NMR, se conhecido.
        msa_column: Coluna de MSA, se a fonte for MSA.
        feature_label: Rotulo NCBI/motivo.
        atom_name: Atomo, se a selecao for atomica.
        message: Texto honesto (gap, unavailable).
        status: READY, UNMAPPED, STALE, UNAVAILABLE, INVALID_INPUT.
        extra: Metadados ja calculados noutro modulo (conservacao, metodo).

    Returns:
        Dict da selecao. query_indices e uma lista de int unicos ordenados.
        selection_hash e um SHA-256 dos campos de identidade (nao e um score).

    Raises:
        SelectionError: INVALID_INPUT se molecule ou source forem invalidos.
    """
    kind = str(molecule or "").strip().upper()
    if kind not in MOLECULE_KINDS:
        raise SelectionError(f"Unknown molecule kind: {molecule}.", "INVALID_INPUT")
    origin = str(source or "").strip().lower()
    if origin not in SELECTION_SOURCES:
        raise SelectionError(f"Unknown selection source: {source}.", "INVALID_INPUT")
    indices: List[int] = []
    seen = set()
    for item in query_indices:
        value = int(item)
        if value < 0 or value in seen:
            continue
        seen.add(value)
        indices.append(value)
    indices.sort()
    payload = {
        "status": str(status or "READY"),
        "molecule": kind,
        "sequence_hash": str(sequence_hash or ""),
        "query_indices": indices,
        "source": origin,
        "label": str(label or ""),
        "structure_id": structure_id,
        "structure_hash": structure_hash,
        "chain_id": chain_id,
        "model": model,
        "msa_column": msa_column,
        "feature_label": feature_label,
        "atom_name": atom_name,
        "message": str(message or ""),
        "event_kind": _event_kind(origin, indices),
    }
    if extra:
        for key, value in extra.items():
            if key not in payload:
                payload[str(key)] = value
    payload["selection_hash"] = selection_identity_hash(payload)
    return payload


def selection_identity_hash(selection: Mapping[str, Any]) -> str:
    """Hash estavel da identidade da selecao. Nao e um score cientifico.

    Args:
        selection: MoleculeSelection ou campos equivalentes.

    Returns:
        SHA-256 hex. Campos cientificos duplicados em extra nao entram.

    Raises:
        Nenhum.
    """
    payload = {
        "molecule": selection.get("molecule"),
        "sequence_hash": selection.get("sequence_hash"),
        "source": selection.get("source"),
        "status": selection.get("status"),
        "query_indices": list(selection.get("query_indices") or []),
        "msa_column": selection.get("msa_column"),
        "structure_id": selection.get("structure_id"),
        "structure_hash": selection.get("structure_hash"),
        "feature_label": selection.get("feature_label"),
        "atom_name": selection.get("atom_name"),
        "chain_id": selection.get("chain_id"),
        "model": selection.get("model"),
    }
    encoded = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def query_indices_of(selection: Optional[Mapping[str, Any]]) -> List[int]:
    """Devolve os indices 0-based da selecao, ou lista vazia."""
    if not isinstance(selection, Mapping):
        return []
    values: List[int] = []
    for item in list(selection.get("query_indices") or []):
        try:
            values.append(int(item))
        except (TypeError, ValueError):
            continue
    return values


def primary_index(selection: Optional[Mapping[str, Any]]) -> Optional[int]:
    """Primeiro indice da selecao, ou None se vazia/unmapped."""
    indices = query_indices_of(selection)
    if not indices:
        return None
    return indices[0]


def selection_applies_to(
    selection: Optional[Mapping[str, Any]],
    *,
    sequence_hash: str,
    molecule: str,
    structure_id: object = None,
    structure_hash: object = None,
) -> bool:
    """True somente se a selecao for da mesma sequencia, molecula e estrutura.

    Hashes vazios nunca coincidem. Uma selecao de estrutura A nao se aplica a B.
    """
    if not isinstance(selection, Mapping):
        return False
    if not provenance.hashes_match(selection.get("sequence_hash"), sequence_hash):
        return False
    if not molecules_are_compatible(selection.get("molecule"), molecule):
        return False
    stored_structure = selection.get("structure_id")
    if stored_structure not in (None, "") and structure_id not in (None, ""):
        if str(stored_structure) != str(structure_id):
            return False
    stored_hash = selection.get("structure_hash")
    if stored_hash not in (None, "") and structure_hash not in (None, ""):
        if not provenance.hashes_match(stored_hash, structure_hash):
            return False
    return True


def selection_from_sequence_position(
    *,
    molecule: str,
    sequence: str,
    position: int,
    source: str = "sequence",
    label: str = "",
) -> dict:
    """Selecao de uma posicao 0-based da sequencia de analise."""
    cleaned = str(sequence or "")
    digest = provenance.sequence_digest(cleaned)
    if position < 0 or position >= len(cleaned):
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source=source,
            label=label or f"sequence position {position}",
            status="INVALID_INPUT",
            message="Position is outside the analysis sequence.",
        )
    return make_selection(
        molecule=molecule,
        sequence_hash=digest,
        query_indices=[position],
        source=source,
        label=label or f"sequence position {position}",
    )


def selection_from_structure_pick(
    *,
    sequence_hash: str,
    structure_id: object,
    structure_hash: object = None,
    chain_id: object = None,
    model: object = None,
    query_index: object = None,
    atom_name: object = None,
    label_seq_id: object = None,
    auth_seq_id: object = None,
) -> dict:
    """Selecao a partir de um ponto 3D ja mapeado. Nao adivinha o indice."""
    if query_index is None:
        return make_selection(
            molecule="PROTEIN",
            sequence_hash=sequence_hash,
            query_indices=[],
            source="structure",
            structure_id=structure_id,
            structure_hash=structure_hash,
            chain_id=chain_id,
            model=model,
            atom_name=atom_name,
            status="UNMAPPED",
            message="No mapped structural residue",
            extra={"label_seq_id": label_seq_id, "auth_seq_id": auth_seq_id},
        )
    return make_selection(
        molecule="PROTEIN",
        sequence_hash=sequence_hash,
        query_indices=[int(query_index)],
        source="structure",
        structure_id=structure_id,
        structure_hash=structure_hash,
        chain_id=chain_id,
        model=model,
        atom_name=atom_name,
        extra={"label_seq_id": label_seq_id, "auth_seq_id": auth_seq_id},
    )


def selection_from_msa_column(
    msa_result: Mapping[str, Any],
    column: int,
    *,
    sequence: str,
    molecule: str = "PROTEIN",
    structure_result: Optional[Mapping[str, Any]] = None,
) -> dict:
    """MSA column -> posicao original. Gap nao cria residuo 3D.

    Args:
        msa_result: Envelope de MSA.
        column: Coluna 0-based.
        sequence: Sequencia de analise (proteina/DNA/RNA) cujo hash deve estar no MSA.
        molecule: Tipo da sequencia de analise.
        structure_result: Envelope estrutural opcional para anexar mapping.

    Returns:
        MoleculeSelection. status UNMAPPED quando a coluna e um gap nessa sequencia.
    """
    digest = provenance.sequence_digest(str(sequence or ""))
    hashes = [str(item) for item in list(msa_result.get("input_hashes") or [])]
    if digest not in hashes:
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source="msa",
            msa_column=column,
            status="STALE",
            message="This MSA does not include the current analysis sequence.",
        )
    row_index = _msa_row_index(msa_result, digest)
    maps = list(msa_result.get("coordinate_maps") or [])
    if row_index is None or row_index >= len(maps):
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source="msa",
            msa_column=column,
            status="UNAVAILABLE",
            message="MSA coordinate map is unavailable for this sequence.",
        )
    mapping_row = maps[row_index]
    if column < 0 or column >= len(mapping_row):
        raise SelectionError("MSA column is outside the coordinate map.", "INVALID_INPUT")
    position = mapping_row[column]
    conservation = msa_result.get("conservation") or {}
    scores = list(conservation.get("scores") or [])
    score = scores[column] if column < len(scores) else None
    extra = {
        "conservation": score,
        "conservation_method": conservation.get("method") or "",
        "conservation_disclaimer": (
            "Conservation among analyzed sequences. Not a functional residue annotation."
        ),
        "msa_tool": msa_result.get("tool"),
    }
    if position is None:
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source="msa",
            label=f"MSA column {column}",
            msa_column=column,
            status="UNMAPPED",
            message=GAP_MESSAGE,
            extra=extra,
        )
    mapped = None
    structure_id = None
    structure_hash = None
    has_coordinates = None
    if isinstance(structure_result, Mapping):
        if not provenance.hashes_match(structure_result.get("sequence_hash"), digest):
            return make_selection(
                molecule=molecule,
                sequence_hash=digest,
                query_indices=[int(position)],
                source="msa",
                label=f"MSA column {column}",
                msa_column=column,
                status="STALE",
                message="Loaded structure does not match the current analysis sequence.",
                extra=extra,
            )
        mapped = protein_structure.map_msa_column_to_structure_residue(
            mapping_row, column, structure_result
        )
        structure_id = structure_result.get("structure_id")
        structure_hash = structure_result.get("content_hash")
        if mapped is not None:
            has_coordinates = bool(mapped.get("has_coordinates"))
            extra["label_seq_id"] = mapped.get("label_seq_id")
            extra["auth_seq_id"] = mapped.get("auth_seq_id")
            extra["chain_id"] = mapped.get("chain_id")
            extra["mapping_row_status"] = mapped.get("mapping_row_status")
    if mapped is None and isinstance(structure_result, Mapping):
        extra["structure_message"] = "No mapped structural residue"
    extra["has_coordinates"] = has_coordinates
    return make_selection(
        molecule=molecule,
        sequence_hash=digest,
        query_indices=[int(position)],
        source="msa",
        label=f"MSA column {column}",
        msa_column=column,
        structure_id=structure_id,
        structure_hash=structure_hash,
        chain_id=None if mapped is None else mapped.get("chain_id"),
        extra=extra,
        message=(
            ""
            if has_coordinates is not False
            else "Coordinates unavailable for this residue"
        ),
    )


def selection_from_motif_hit(
    hit: Mapping[str, Any],
    *,
    molecule: str,
    sequence: str,
) -> dict:
    """Motif hit [start, end) -> indices da sequencia. Nao adivinha mapping 3D."""
    digest = provenance.sequence_digest(str(sequence or ""))
    try:
        start = int(hit.get("start"))
        end = int(hit.get("end"))
    except (TypeError, ValueError) as exc:
        raise SelectionError("Motif hit has no valid span.", "INVALID_INPUT") from exc
    if start < 0 or end > len(sequence) or start >= end:
        raise SelectionError("Motif hit falls outside the analysis sequence.", "INVALID_INPUT")
    indices = list(range(start, end))
    return make_selection(
        molecule=molecule,
        sequence_hash=digest,
        query_indices=indices,
        source="motif",
        label=str(hit.get("pattern") or hit.get("match") or "motif"),
        feature_label=str(hit.get("match") or ""),
        extra={"start": start, "end": end, "strand": hit.get("strand") or "+"},
    )


def selection_from_ncbi_spans(
    spans: Sequence[Mapping[str, Any]],
    *,
    sequence: str,
    molecule: str,
    accession: object = None,
) -> dict:
    """Spans NCBI 0-based -> indices. Recusa molecula incompatível no chamador.

    Args:
        spans: Saida de feature_spans_for_analysis / make_span.
        sequence: Sequencia de analise (mesmo alfabeto do registro).
        molecule: DNA, RNA ou PROTEIN.
        accession: Accession NCBI, se conhecido.

    Returns:
        Selecao com a uniao dos spans que cabem na sequencia.
    """
    digest = provenance.sequence_digest(str(sequence or ""))
    if not spans:
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source="ncbi",
            status="UNAVAILABLE",
            message="No NCBI feature spans are stored.",
            extra={"accession": accession},
        )
    indices: List[int] = []
    labels = []
    for span in spans:
        try:
            start = int(span.get("start"))
            end = int(span.get("end"))
        except (TypeError, ValueError):
            continue
        if start < 0 or end > len(sequence) or start >= end:
            continue
        indices.extend(range(start, end))
        labels.append(str(span.get("label") or span.get("kind") or "NCBI feature"))
    if not indices:
        return make_selection(
            molecule=molecule,
            sequence_hash=digest,
            query_indices=[],
            source="ncbi",
            status="UNMAPPED",
            message="Stored NCBI feature spans do not fit the current analysis sequence.",
            extra={"accession": accession},
        )
    return make_selection(
        molecule=molecule,
        sequence_hash=digest,
        query_indices=indices,
        source="ncbi",
        label=labels[0] if labels else "NCBI feature",
        feature_label=labels[0] if labels else "",
        extra={
            "accession": accession,
            "feature_source": "NCBI",
            "n_spans": len(spans),
        },
    )


def selection_from_tree_leaf(
    tree_result: Mapping[str, Any],
    tree_id: str,
    *,
    sequence: str = "",
    molecule: str = "DNA",
) -> dict:
    """Selecao ao nivel da sequencia a partir de uma folha da arvore.

    Args:
        tree_result: Envelope de inferencia filogenetica.
        tree_id: Identificador Newick da folha.
        sequence: Sequencia de analise opcional; se o hash nao coincidir, STALE.
        molecule: DNA, RNA ou PROTEIN.

    Returns:
        MoleculeSelection. Indices vazios: a selecao e a sequencia inteira
        (folha), nao um residuo inventado.

    Raises:
        SelectionError: INVALID_INPUT se a folha nao existir.

    Nota biologica:
        Destacar uma folha nao cria relacao CRISPR nem estrutura 3D.
    """
    alignment_hash = str(tree_result.get("alignment_hash") or "")
    tree_hash = str(tree_result.get("tree_hash") or "")
    target = str(tree_id or "").strip()
    leaf = None
    for item in list(tree_result.get("leaves") or []):
        if str(item.get("tree_id") or "") == target:
            leaf = item
            break
    if leaf is None:
        raise SelectionError(
            f"Tree leaf '{tree_id}' is not in this phylogenetic result.",
            "INVALID_INPUT",
        )
    digest = str(leaf.get("sequence_hash") or "")
    if sequence:
        live = provenance.sequence_digest(str(sequence))
        if live != digest:
            return make_selection(
                molecule=molecule,
                sequence_hash=live,
                query_indices=[],
                source="tree",
                label=str(leaf.get("original_id") or target),
                status="STALE",
                message="This tree leaf does not match the current analysis sequence.",
                extra={
                    "tree_id": target,
                    "tree_hash": tree_hash,
                    "alignment_hash": alignment_hash,
                    "leaf_sequence_hash": digest,
                },
            )
    return make_selection(
        molecule=str(tree_result.get("molecule") or molecule),
        sequence_hash=digest,
        query_indices=[],
        source="tree",
        label=str(leaf.get("original_id") or target),
        extra={
            "tree_id": target,
            "original_id": leaf.get("original_id"),
            "tree_hash": tree_hash,
            "alignment_hash": alignment_hash,
            "accession": leaf.get("accession") or "",
            "organism": leaf.get("organism_display") or "",
        },
    )


def _msa_row_index(msa_result: Mapping[str, Any], digest: str) -> Optional[int]:
    rows = list(msa_result.get("rows") or [])
    for index, row in enumerate(rows):
        if str(row.get("hash") or "") == digest:
            return index
    hashes = [str(item) for item in list(msa_result.get("input_hashes") or [])]
    try:
        return hashes.index(digest)
    except ValueError:
        return None


def ungapped_sequence_for_hash(msa_result: Mapping[str, Any], digest: str) -> str:
    """Sequencia sem gaps do MSA para um hash. Vazio se o mapping nao existir.

    Args:
        msa_result: Envelope MSA.
        digest: sequence_hash da folha.

    Returns:
        Sequencia ungapped ou string vazia.

    Raises:
        Nenhum.
    """
    target = str(digest or "")
    if not target:
        return ""
    for row in list(msa_result.get("rows") or []):
        if str(row.get("hash") or "") == target:
            return str(row.get("ungapped") or "")
    for member in list(msa_result.get("input_members") or []):
        if str(member.get("hash") or "") == target:
            return str(member.get("sequence") or "")
    return ""


def structure_availability_for_hash(
    sequence_hash: str,
    *,
    protein_result: Optional[Mapping[str, Any]] = None,
    nucleic_result: Optional[Mapping[str, Any]] = None,
    complex_result: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Declara se uma estrutura ja validada na sessao corresponde a este hash.

    Args:
        sequence_hash: Hash da folha/sequencia.
        protein_result: Envelope PDB/AlphaFold da aba Protein.
        nucleic_result: Envelope DNA/RNA depositado ou ilustrativo.
        complex_result: Complexo CRISPR depositado (nao implica o hash da folha).

    Returns:
        Dict status AVAILABLE ou UNAVAILABLE. Nao faz fetch. Nao inventa PDB.

    Raises:
        Nenhum.

    Nota biologica:
        Uma folha de arvore sem estrutura mapeada permanece
        "No validated structure available". O complexo CRISPR so conta se o
        hash da query coincidir com um polimero depositado.
    """
    digest = str(sequence_hash or "")
    candidates = [
        ("protein", protein_result),
        ("nucleic", nucleic_result),
        ("complex", complex_result),
    ]
    for source_name, envelope in candidates:
        if not isinstance(envelope, Mapping):
            continue
        stored = str(envelope.get("sequence_hash") or "")
        if digest and stored and stored == digest:
            kind = str(envelope.get("kind") or "")
            return {
                "status": "AVAILABLE",
                "source_layer": source_name,
                "kind": kind,
                "kind_label": protein_structure.structure_kind_label(kind),
                "structure_id": envelope.get("structure_id") or "",
                "source": envelope.get("source") or "",
                "sequence_hash": stored,
            }
    return {
        "status": "UNAVAILABLE",
        "source_layer": "",
        "kind": "unavailable",
        "kind_label": "Coordinates Unavailable",
        "structure_id": "",
        "source": "",
        "sequence_hash": digest,
        "reason": "No validated structure available.",
    }


def _event_kind(source: str, indices: Sequence[int]) -> str:
    if source == "msa":
        return "structureFeatureSelected"
    if source == "motif":
        return "structureFeatureSelected"
    if source == "ncbi":
        return "structureFeatureSelected"
    if source == "tree":
        return "structureFeatureSelected"
    if len(indices) > 1:
        return "structureFeatureSelected"
    return "sequencePositionSelected"
