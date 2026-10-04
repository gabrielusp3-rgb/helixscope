"""Caminho estavel DNA sequence -> representacao 3D.

Camadas distintas (nao o mesmo codigo tentando as tres coisas):

* ANALYSIS: fora deste modulo (dna_analysis).
* NAVIGATION: janela contigua, salto, overview com teto de objetos.
* ATOMIC_DETAIL: helice ILLUSTRATIVE ou envelope EXPERIMENTAL/PREDICTED depositado.

PREDICTED DNA 3D sem modelo externo permanece UNAVAILABLE. Nao ha motor
interno de fold de DNA. ILLUSTRATIVE nunca e relabelado como EXPERIMENTAL.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from . import (
    dna_analysis,
    motif_search,
    na_structure_catalog,
    nucleic_geometry,
    protein_structure,
    provenance,
    region_nav,
    structure_scene,
)

_FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures"

LAYER_ANALYSIS: str = "analysis"
LAYER_NAVIGATION: str = "navigation"
LAYER_ATOMIC_DETAIL: str = "atomic_detail"

KIND_EXPERIMENTAL: str = "experimental"
KIND_PREDICTED: str = "predicted"
KIND_ILLUSTRATIVE: str = "illustrative"
KIND_UNAVAILABLE: str = "unavailable"
DNA_3D_KINDS: tuple[str, ...] = (
    KIND_EXPERIMENTAL,
    KIND_PREDICTED,
    KIND_ILLUSTRATIVE,
    KIND_UNAVAILABLE,
)

MAX_BASE_LABELS: int = 24
"""Teto de rotulos de base simultaneos. Independente do comprimento da sequencia."""

PREFER_EXPERIMENTAL: str = "experimental"
PREFER_PREDICTED: str = "predicted"
PREFER_ILLUSTRATIVE: str = "illustrative"
PREFER_MODES: tuple[str, ...] = (
    PREFER_EXPERIMENTAL,
    PREFER_PREDICTED,
    PREFER_ILLUSTRATIVE,
)


class Dna3dError(Exception):
    """Falha classificada do caminho DNA 3D.

    Attributes:
        category: INVALID_INPUT, UNAVAILABLE, STALE, ERROR.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def normalize_kind(kind: object) -> str:
    """Normaliza o kind estrutural DNA 3D.

    Args:
        kind: experimental, predicted, illustrative ou unavailable.

    Returns:
        Um de DNA_3D_KINDS.

    Raises:
        Dna3dError: INVALID_INPUT.
    """
    key = str(kind or "").strip().lower()
    if key not in DNA_3D_KINDS:
        raise Dna3dError("DNA 3D kind must be experimental, predicted, illustrative or unavailable.")
    return key


def predicted_availability() -> dict:
    """Declara que nao ha preditor 3D de DNA integrado.

    Args:
        Nenhum.

    Returns:
        Dict status=UNAVAILABLE, kind=unavailable, ferramenta ausente.

    Raises:
        Nenhum.

    Nota biologica:
        AlphaFold-Multimer / NA e motores tipo 3DNA nao estao neste app.
        Nao se inventa um fold de DNA para preencher PREDICTED.
    """
    return provenance.analysis_envelope(
        module="DNA_3D",
        payload={
            "kind": KIND_UNAVAILABLE,
            "tool": None,
            "model": None,
            "reason": (
                "No DNA 3D prediction engine is integrated (no AlphaFold-NA, "
                "no 3DNA/fiber model beyond the ILLUSTRATIVE canonical helix). "
                "PREDICTED remains UNAVAILABLE unless a deposited predicted "
                "envelope is supplied by the caller."
            ),
        },
        status="UNAVAILABLE",
        algorithm="none",
        parameters={},
        source="HelixScope DNA 3D",
    ) | {
        "kind": KIND_UNAVAILABLE,
        "status": "UNAVAILABLE",
        "reason": (
            "No DNA 3D prediction engine is integrated (no AlphaFold-NA, "
            "no 3DNA/fiber model beyond the ILLUSTRATIVE canonical helix). "
            "PREDICTED remains UNAVAILABLE unless a deposited predicted "
            "envelope is supplied by the caller."
        ),
        "tool": None,
        "model": None,
    }


def bundled_experimental_for_sequence(sequence: str, molecule: str) -> Optional[dict]:
    """Envelope EXPERIMENTAL a partir da copia mmCIF no repositorio.

    Args:
        sequence: Polimero de analise ja validado.
        molecule: DNA ou RNA.

    Returns:
        structure_3d_input se a sequencia for identica ao polimero do
        catalogo e o ficheiro tests/fixtures/{pdb}.cif existir. None se
        nao houver match. None nao e uma estrutura.

    Raises:
        Nenhum. Falhas de parse voltam None.

    Nota biologica:
        1BNA e 1RNA sao copias publicas RCSB, nao geometria procedural.
        Nao e fetch live. O source declara bundled deposited mmCIF copy.
        Implementacao heuristica simplificada inspirada em carregamento
        local de mmCIF; nao reproduz o modelo/algoritmo original publicado
        de um refinador cristalografico.
    """
    mol = str(molecule or "").strip().upper()
    seq = str(sequence or "").strip().upper()
    if mol == "DNA":
        catalog = na_structure_catalog.list_experimental_dna_structures()
    elif mol == "RNA":
        catalog = na_structure_catalog.list_experimental_rna_structures()
    else:
        return None
    match = next(
        (
            item
            for item in catalog
            if str(item.get("polymer_sequence") or "").upper() == seq
        ),
        None,
    )
    if match is None:
        return None
    pdb_id = str(match.get("pdb_id") or "").upper()
    cif_path = _FIXTURES / f"{pdb_id}.cif"
    meta_path = _FIXTURES / f"rcsb_entry_{pdb_id}.json"
    if not cif_path.is_file() or not meta_path.is_file():
        return None

    def _no_network(*_args: Any, **_kwargs: Any) -> Any:
        raise protein_structure.StructureError(
            "Bundled mmCIF is local; network fetch was not used.",
            "NETWORK_ERROR",
        )

    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        result = protein_structure.load_protein_structure(
            seq,
            source="RCSB PDB",
            structure_id=pdb_id,
            structure_text=cif_path.read_text(encoding="utf-8"),
            metadata=meta,
            urlopen_fn=_no_network,
            molecule=mol,
        )
        envelope = protein_structure.structure_3d_input(result)
    except (protein_structure.StructureError, TypeError, ValueError, KeyError, OSError):
        return None
    envelope["bundled_fixture"] = cif_path.name
    envelope["source"] = f"RCSB PDB {pdb_id} (bundled deposited mmCIF copy)"
    envelope["molecule_type"] = mol
    return envelope


def overview_object_budget(length: int, lod: object = region_nav.LOD_HIGH) -> dict:
    """Teto de objetos de cena para overview. Nao escala com 50k nt.

    Args:
        length: Comprimento da sequencia de analise.
        lod: Nivel visual da janela de detalhe.

    Returns:
        Dict layer=navigation, max_nt, max_atoms, max_pairs, max_labels.

    Raises:
        Dna3dError: INVALID_INPUT se length < 0.

    Nota:
        50k e escala de analise/navegacao. Atomic detail usa no maximo
        helix_limit_for_lod nucleotidos consecutivos.
    """
    if int(length) < 0:
        raise Dna3dError("Sequence length must be >= 0.")
    cap = region_nav.helix_limit_for_lod(lod)
    window = min(int(length), cap) if length else 0
    return {
        "layer": LAYER_NAVIGATION,
        "full_length": int(length),
        "max_nt_atomic_detail": cap,
        "window_nt": window,
        "max_atoms": window * 2,
        "max_pairs": window,
        "max_strands": 2,
        "max_base_labels": MAX_BASE_LABELS,
        "scales_linearly_with_full_length": False,
        "note": (
            "Overview/detail object count is capped by the contiguous helix "
            "window, not by full sequence length."
        ),
    }


def navigate(
    length: int,
    *,
    action: str,
    start: int = 0,
    end: Optional[int] = None,
    position: int = 0,
    lod: object = region_nav.LOD_HIGH,
) -> dict:
    """Navegacao de janela contigua. Nao altera a sequencia.

    Args:
        length: Comprimento total.
        action: overview, jump, region, zoom, reset.
        start: Pedido de inicio (region/zoom).
        end: Pedido de fim exclusivo.
        position: Alvo do jump (0-based).
        lod: LOD da janela de detalhe.

    Returns:
        Dict layer=navigation com start/end/n/partial_view.

    Raises:
        Dna3dError: INVALID_INPUT.

    Nota:
        start>end, vazio e fora dos limites sao RegionError convertidos.
        Nao ha clamp silencioso.
    """
    key = str(action or "").strip().lower()
    if key not in {"overview", "jump", "region", "zoom", "reset"}:
        raise Dna3dError("Navigation action must be overview, jump, region, zoom or reset.")
    try:
        if key in {"overview", "reset"}:
            view = region_nav.contiguous_view(length, start=0, end=length, lod=lod)
        elif key == "jump":
            view = region_nav.contiguous_view(
                length, start=0, end=length, lod=lod, center=int(position)
            )
        else:
            view = region_nav.contiguous_view(length, start=start, end=end, lod=lod)
    except region_nav.RegionError as exc:
        raise Dna3dError(str(exc), exc.category) from exc
    view = dict(view)
    view["layer"] = LAYER_NAVIGATION
    view["action"] = key
    return view


def representation_for_sequence(
    sequence: str,
    *,
    deposited: Optional[Mapping[str, Any]] = None,
    prefer: Optional[str] = None,
    region_start: int = 0,
    region_end: Optional[int] = None,
    lod: str = region_nav.LOD_HIGH,
    center: Optional[int] = None,
    highlight_indices: Optional[Sequence[int]] = None,
) -> dict:
    """Caminho estavel: DNA validado -> representacao 3D.

    Args:
        sequence: DNA 5'->3'.
        deposited: Envelope estrutural ja validado (RCSB/parser). Nao faz fetch.
        prefer: experimental, predicted ou illustrative. None = experimental
            se depositado compativel, senao illustrative.
        region_start: Janela visual (so ILLUSTRATIVE).
        region_end: Fim exclusivo da janela.
        lod: LOD da helice ilustrativa.
        center: Centro opcional.
        highlight_indices: Posicoes 0-based da sequencia completa.

    Returns:
        Dict status, kind, layer, envelope, model, partial_view, counts.
        Sequencia valida com prefer default nunca devolve kind=unavailable.

    Raises:
        Nenhum. Falhas voltam como status ERROR/UNAVAILABLE no dict.

    Nota biologica:
        Implementacao heuristica simplificada inspirada em parametros
        helicoidais B-DNA (rise 3.38 A, 10 bp/volta, twist 36 deg; fibre
        B-DNA classico, nao 10.5 bp/volta em solucao);
        nao reproduz o modelo/algoritmo original publicado de um
        refinador cristalografico. PREDICTED sem envelope depositado
        permanece UNAVAILABLE.
    """
    info = dna_analysis.validate_for_molecule(sequence, "DNA")
    if not info["is_valid"]:
        return _unavailable(
            "ERROR",
            KIND_UNAVAILABLE,
            str(info.get("rejection_reason") or "Sequence is not DNA."),
            sequence="",
        )
    seq = str(info["sequence"])
    digest = provenance.sequence_digest(seq)
    mode = str(prefer or "").strip().lower() or None
    if deposited is None and mode in {None, PREFER_EXPERIMENTAL}:
        deposited = bundled_experimental_for_sequence(seq, "DNA")
    if mode is not None and mode not in PREFER_MODES:
        return _unavailable("ERROR", KIND_UNAVAILABLE, "prefer must be experimental, predicted or illustrative.", seq)

    if mode == PREFER_PREDICTED:
        deposited_ok = _deposited_usable(deposited, digest, KIND_PREDICTED)
        if not deposited_ok:
            payload = predicted_availability()
            payload["sequence_hash"] = digest
            payload["layer"] = LAYER_ATOMIC_DETAIL
            payload["envelope"] = None
            payload["model"] = None
            payload["n_scene_objects"] = 0
            return payload

    if deposited is not None:
        dep_hash = str(deposited.get("sequence_hash") or "")
        dep_kind = str(deposited.get("kind") or "").strip().lower()
        if dep_hash and dep_hash != digest:
            if mode == PREFER_ILLUSTRATIVE:
                pass
            else:
                return _unavailable(
                    "UNMAPPED",
                    KIND_UNAVAILABLE,
                    "Deposited structure sequence hash does not match the analysis sequence. "
                    "HelixScope does not force a mapping or substitute ILLUSTRATIVE automatically.",
                    seq,
                )
        elif dep_kind in {KIND_EXPERIMENTAL, KIND_PREDICTED} and (
            mode in {None, PREFER_EXPERIMENTAL, PREFER_PREDICTED}
        ):
            if mode == PREFER_PREDICTED and dep_kind != KIND_PREDICTED:
                payload = predicted_availability()
                payload["sequence_hash"] = digest
                payload["layer"] = LAYER_ATOMIC_DETAIL
                payload["envelope"] = None
                payload["model"] = None
                payload["n_scene_objects"] = 0
                return payload
            if mode == PREFER_EXPERIMENTAL and dep_kind != KIND_EXPERIMENTAL:
                return _unavailable(
                    "UNAVAILABLE",
                    KIND_UNAVAILABLE,
                    "Attached structure is not experimental. "
                    "ILLUSTRATIVE helix is available via prefer=illustrative or prefer=None.",
                    seq,
                )
            try:
                if deposited.get("mappings") is None:
                    envelope = protein_structure.structure_3d_input(deposited)
                else:
                    envelope = dict(deposited)
            except (protein_structure.StructureError, TypeError, ValueError, KeyError) as exc:
                return _unavailable("ERROR", KIND_UNAVAILABLE, str(exc), seq)
            status = "EXPERIMENTAL" if dep_kind == KIND_EXPERIMENTAL else "PREDICTED"
            n_atoms = len(list(envelope.get("atoms") or []))
            return provenance.analysis_envelope(
                module="DNA_3D",
                payload={"kind": dep_kind, "n_atoms": n_atoms},
                status=status,
                algorithm=str((envelope.get("metadata") or {}).get("method") or "deposited structure"),
                parameters={"prefer": mode or "auto"},
                source=str(envelope.get("source") or "deposited"),
                sequence=seq,
            ) | {
                "status": status,
                "kind": dep_kind,
                "kind_label": envelope.get("kind_label") or dep_kind,
                "layer": LAYER_ATOMIC_DETAIL,
                "sequence_hash": digest,
                "envelope": envelope,
                "model": None,
                "partial_view": bool(envelope.get("partial_view")),
                "n_strands": len(list(envelope.get("chains") or [])),
                "n_pairs": len(list(envelope.get("pairs") or [])),
                "n_atoms": n_atoms,
                "n_scene_objects": n_atoms,
                "max_base_labels": MAX_BASE_LABELS,
                "disclaimer": envelope.get("disclaimer") or "",
            }

    if mode == PREFER_EXPERIMENTAL:
        return _unavailable(
            "UNAVAILABLE",
            KIND_UNAVAILABLE,
            "No deposited experimental DNA structure is attached. "
            "ILLUSTRATIVE helix is available via prefer=illustrative or prefer=None.",
            seq,
        )

    return _illustrative_payload(
        seq,
        region_start=region_start,
        region_end=region_end,
        lod=lod,
        center=center,
        highlight_indices=highlight_indices,
    )


def scene_from_representation(
    representation: Mapping[str, Any],
    *,
    selected_query_index: Optional[int] = None,
    highlight_query_indices: Optional[Sequence[int]] = None,
    lod: str = region_nav.LOD_HIGH,
) -> dict:
    """Constroi o scene model a partir da representacao. Sem ciencia nova.

    Args:
        representation: Saida de representation_for_sequence.
        selected_query_index: Base destacada.
        highlight_query_indices: Outros indices.
        lod: LOD visual do backbone.

    Returns:
        Scene model ou dict status UNAVAILABLE.

    Raises:
        Nenhum.
    """
    envelope = representation.get("envelope")
    if not isinstance(envelope, Mapping) or str(representation.get("kind") or "") == KIND_UNAVAILABLE:
        return {
            "status": "UNAVAILABLE",
            "kind": KIND_UNAVAILABLE,
            "message": str(representation.get("reason") or "No DNA 3D envelope."),
        }
    try:
        scene = structure_scene.build_scene(
            envelope,
            expected_sequence_hash=str(representation.get("sequence_hash") or ""),
            representation="backbone",
            color_mode="chain",
            selected_query_index=selected_query_index,
            highlight_query_indices=highlight_query_indices,
            lod=lod,
        )
    except (structure_scene.SceneError, TypeError, ValueError, KeyError) as exc:
        return {
            "status": "ERROR",
            "kind": KIND_UNAVAILABLE,
            "message": str(exc),
        }
    return scene


def highlight_from_motif(sequence: str, pattern: str) -> dict:
    """Motivo real na sequencia -> indices 3D. Sem match = lista vazia.

    Args:
        sequence: DNA.
        pattern: Padrao IUPAC (ex. GAATTC).

    Returns:
        Dict status, indices, spans. Nao inventa hits.

    Raises:
        ValueError: Padrao invalido (repassado de motif_search).
    """
    info = dna_analysis.validate_for_molecule(sequence, "DNA")
    if not info["is_valid"]:
        raise Dna3dError(str(info.get("rejection_reason") or "not DNA."))
    seq = str(info["sequence"])
    hits = motif_search.find_motif(seq, pattern)
    indices: List[int] = []
    spans = []
    for hit in hits:
        start = int(hit["start"])
        end = int(hit["end"])
        indices.extend(range(start, end))
        spans.append({"start": start, "end": end, "match": hit.get("match")})
    return {
        "status": "COMPUTED",
        "kind": "sequence_span",
        "indices": sorted(set(indices)),
        "spans": spans,
        "n_hits": len(hits),
        "pattern": pattern,
        "sequence_hash": provenance.sequence_digest(seq),
    }


def highlight_from_guide(sequence: str, guide: Mapping[str, Any]) -> dict:
    """Protospacer+PAM so se forem substring contigua real.

    Args:
        sequence: DNA alvo.
        guide: Dict com guide_sequence, start_0based, pam_sequence, pam_start_0based, strand.

    Returns:
        Dict status EXACT ou UNMAPPED. UNMAPPED nao devolve highlight forçado.

    Raises:
        Dna3dError: INVALID_INPUT se a sequencia nao for DNA.
    """
    info = dna_analysis.validate_for_molecule(sequence, "DNA")
    if not info["is_valid"]:
        raise Dna3dError(str(info.get("rejection_reason") or "not DNA."))
    seq = str(info["sequence"])
    spacer = str(guide.get("guide_sequence") or "").upper()
    pam = str(guide.get("pam_sequence") or "").upper()
    try:
        start = int(guide.get("start_0based"))
        pam_start = int(guide.get("pam_start_0based")) if guide.get("pam_start_0based") is not None else None
    except (TypeError, ValueError) as exc:
        raise Dna3dError("Guide coordinates must be integers.") from exc
    strand = str(guide.get("strand") or "+")
    if not spacer or start < 0 or start + len(spacer) > len(seq):
        return _unmapped_guide("Guide coordinates fall outside the sequence.", seq)
    observed = seq[start : start + len(spacer)]
    if strand == "+" and observed != spacer:
        return _unmapped_guide("Protospacer is not an exact contiguous substring.", seq)
    if strand == "-" and observed != spacer:
        return _unmapped_guide("Protospacer is not an exact contiguous substring on this strand encoding.", seq)
    indices = list(range(start, start + len(spacer)))
    if pam and pam_start is not None:
        if pam_start < 0 or pam_start + len(pam) > len(seq):
            return _unmapped_guide("PAM coordinates fall outside the sequence.", seq)
        window = seq[pam_start : pam_start + len(pam)]
        if not _pam_matches(window, pam):
            return _unmapped_guide("PAM is not present at the declared adjacent position.", seq)
        if strand == "+" and pam_start != start + len(spacer):
            return _unmapped_guide("PAM is not adjacent 3' of the SpCas9 protospacer.", seq)
        indices.extend(range(pam_start, pam_start + len(pam)))
    return {
        "status": "EXACT",
        "kind": KIND_ILLUSTRATIVE,
        "indices": sorted(set(indices)),
        "protospacer_start": start,
        "pam_start": pam_start,
        "sequence_hash": provenance.sequence_digest(seq),
    }


def highlight_from_ncbi_spans(
    sequence: str,
    spans: Sequence[Mapping[str, Any]],
    *,
    accession: object = None,
) -> dict:
    """Spans NCBI ja parseados -> indices. Nao faz fetch.

    Args:
        sequence: Sequencia de analise.
        spans: Saida de ncbi_fetch.feature_spans_for_analysis.
        accession: Accession se conhecido.

    Returns:
        Dict status, indices. Spans que nao cabem sao ignorados, nao inventados.
    """
    from . import molecule_selection

    selection = molecule_selection.selection_from_ncbi_spans(
        spans, sequence=sequence, molecule="DNA", accession=accession
    )
    return {
        "status": selection.get("status"),
        "kind": KIND_ILLUSTRATIVE,
        "indices": list(selection.get("query_indices") or []),
        "selection": selection,
        "sequence_hash": selection.get("sequence_hash"),
    }


def highlight_from_msa_column(
    msa_result: Mapping[str, Any],
    column: int,
    sequence: str,
) -> dict:
    """Coluna MSA -> posicao original. Gap = UNMAPPED.

    Args:
        msa_result: Envelope MSA.
        column: Coluna 0-based.
        sequence: DNA cujo hash esta no MSA.

    Returns:
        Dict status, indices.
    """
    from . import molecule_selection

    selection = molecule_selection.selection_from_msa_column(
        msa_result, int(column), sequence=sequence, molecule="DNA"
    )
    return {
        "status": selection.get("status"),
        "kind": KIND_ILLUSTRATIVE,
        "indices": list(selection.get("query_indices") or []),
        "msa_column": column,
        "selection": selection,
        "sequence_hash": selection.get("sequence_hash"),
    }


def highlight_from_window(length: int, start: int, end: int) -> dict:
    """Janela de grafico [start, end) -> indices. Limites inclusive/exclusivo.

    Args:
        length: Comprimento da sequencia.
        start: Inclusivo.
        end: Exclusivo.

    Returns:
        Dict indices. Erro se vazio ou fora.

    Raises:
        Dna3dError.
    """
    try:
        start_i, end_i = region_nav.clamp_region(length, start, end)
    except region_nav.RegionError as exc:
        raise Dna3dError(str(exc), exc.category) from exc
    return {
        "status": "COMPUTED",
        "kind": KIND_ILLUSTRATIVE,
        "indices": list(range(start_i, end_i)),
        "start": start_i,
        "end": end_i,
        "layer": LAYER_NAVIGATION,
    }


def export_representation(
    representation: Mapping[str, Any],
    *,
    expected_sequence_hash: Optional[str] = None,
) -> dict:
    """Export honesto. UNAVAILABLE/STALE nao vira sucesso.

    Args:
        representation: Saida de representation_for_sequence.
        expected_sequence_hash: Hash da sequencia atual. Divergencia = STALE.

    Returns:
        Dict JSON-safe com provenance, status, hashes, limitations.

    Raises:
        Nenhum.
    """
    status = str(representation.get("status") or "UNAVAILABLE")
    kind = str(representation.get("kind") or KIND_UNAVAILABLE)
    stored_hash = str(representation.get("sequence_hash") or "")
    expected = str(expected_sequence_hash or "")
    stale = bool(expected and stored_hash and expected != stored_hash)
    if stale:
        status = "STALE"
        kind = KIND_UNAVAILABLE
    payload = {
        "status": status,
        "kind": kind,
        "source": representation.get("source"),
        "method": representation.get("algorithm") or representation.get("method"),
        "model": representation.get("model_name") or representation.get("algorithm"),
        "version": provenance.HELIXSCOPE_VERSION,
        "sequence_hash": stored_hash or None,
        "structure_hash": (representation.get("envelope") or {}).get("structure_hash")
        if isinstance(representation.get("envelope"), Mapping)
        else None,
        "hashes": {
            "sequence_hash": stored_hash or None,
            "structure_hash": (representation.get("envelope") or {}).get("structure_hash")
            if isinstance(representation.get("envelope"), Mapping)
            else None,
        },
        "limitations": (
            "stale/incompatible structure"
            if stale
            else representation.get("disclaimer") or representation.get("reason") or ""
        ),
        "successful_result": (not stale)
        and status not in {"UNAVAILABLE", "ERROR", "STALE"},
    }
    return provenance.json_safe(payload)


def position_to_3d(representation: Mapping[str, Any], position: int) -> dict:
    """Indice da sequencia -> residuo 3D se estiver na janela mapeada.

    Args:
        representation: Saida de representation_for_sequence.
        position: Indice 0-based na sequencia completa.

    Returns:
        Dict status READY ou UNMAPPED. Nao inventa coordenadas.

    Raises:
        Dna3dError: INVALID_INPUT se position nao for um inteiro.
    """
    try:
        index = int(position)
    except (TypeError, ValueError) as exc:
        raise Dna3dError("Position must be an integer.") from exc
    envelope = representation.get("envelope")
    if not isinstance(envelope, Mapping):
        return {"status": "UNMAPPED", "query_index_0based": index, "mapping": None}
    for row in list(envelope.get("mappings") or []):
        if row.get("query_index_0based") is None:
            continue
        if int(row["query_index_0based"]) == index:
            return {
                "status": "READY",
                "query_index_0based": index,
                "mapping": dict(row),
                "has_coordinates": bool(row.get("has_coordinates")),
            }
    return {"status": "UNMAPPED", "query_index_0based": index, "mapping": None}


def pick_to_position(payload: object) -> dict:
    """Selecao 3D (customdata Plotly) -> indice da sequencia.

    Args:
        payload: customdata de um ponto da cena.

    Returns:
        Dict status READY ou UNMAPPED.

    Raises:
        Nenhum.
    """
    parsed = structure_scene.parse_selection_payload(payload)
    if parsed is None or parsed.get("query_index_0based") is None:
        return {"status": "UNMAPPED", "query_index_0based": None, "pick": parsed}
    return {
        "status": "READY",
        "query_index_0based": int(parsed["query_index_0based"]),
        "pick": parsed,
    }


def base_labels_for_view(representation: Mapping[str, Any]) -> list:
    """Rotulos de base da janela, no maximo MAX_BASE_LABELS. Visual.

    Args:
        representation: Saida ILLUSTRATIVE.

    Returns:
        Lista de dicts index_0based, base. Nao altera a sequencia.

    Raises:
        Nenhum.
    """
    residues = list((representation.get("model") or {}).get("residues") or [])
    if len(residues) <= MAX_BASE_LABELS:
        return [
            {"index_0based": int(item["index_0based"]), "base": str(item.get("base") or "")}
            for item in residues
        ]
    if not residues:
        return []
    step = max(1, len(residues) // MAX_BASE_LABELS)
    picked = residues[::step][: MAX_BASE_LABELS - 1]
    last = residues[-1]
    if picked and int(picked[-1]["index_0based"]) != int(last["index_0based"]):
        picked.append(last)
    return [
        {"index_0based": int(item["index_0based"]), "base": str(item.get("base") or "")}
        for item in picked[:MAX_BASE_LABELS]
    ]


def window_integrity(representation: Mapping[str, Any], sequence: str) -> dict:
    """Confere janela contigua e bases identicas a substring da sequencia.

    Args:
        representation: Saida ILLUSTRATIVE.
        sequence: Sequencia de analise.

    Returns:
        Dict ok, start, end, n, fragment.

    Raises:
        Dna3dError: INVALID_INPUT se a representacao nao tiver janela.
    """
    start = representation.get("view_start")
    end = representation.get("view_end")
    fragment = str(representation.get("fragment") or "")
    if start is None or end is None:
        raise Dna3dError("Representation has no contiguous window.")
    start_i = int(start)
    end_i = int(end)
    expected = sequence[start_i:end_i]
    residues = list((representation.get("model") or {}).get("residues") or [])
    bases = "".join(str(item.get("base") or "") for item in residues)
    indices = [int(item["index_0based"]) for item in residues]
    contiguous = indices == list(range(start_i, end_i)) if residues else fragment == expected
    return {
        "ok": expected == fragment == bases and contiguous and (end_i - start_i) == len(fragment),
        "start": start_i,
        "end": end_i,
        "n": len(fragment),
        "fragment": fragment,
        "expected": expected,
        "residue_bases": bases,
        "contiguous": contiguous,
    }


def _deposited_usable(deposited: Optional[Mapping[str, Any]], digest: str, kind: str) -> bool:
    if not isinstance(deposited, Mapping):
        return False
    if str(deposited.get("sequence_hash") or "") != digest:
        return False
    return str(deposited.get("kind") or "").strip().lower() == kind


def _illustrative_payload(
    seq: str,
    *,
    region_start: int,
    region_end: Optional[int],
    lod: str,
    center: Optional[int],
    highlight_indices: Optional[Sequence[int]],
) -> dict:
    try:
        model = nucleic_geometry.illustrative_bdna_model(
            seq,
            highlight_indices=list(highlight_indices or []),
            region_start=region_start,
            region_end=region_end,
            lod=lod,
            center=center,
        )
        envelope = nucleic_geometry.structure_3d_input_from_illustrative(model)
    except nucleic_geometry.NucleicGeometryError as exc:
        return _unavailable("ERROR", KIND_UNAVAILABLE, str(exc), seq)
    pairs = list(model.get("pairs") or [])
    residues = list(model.get("residues") or [])
    bases = "".join(str(item.get("base") or "") for item in residues)
    fragment = str(model.get("fragment") or "")
    budget = overview_object_budget(len(seq), lod)
    n_atoms = int(model.get("n_atoms") or 0)
    meta = provenance.analysis_envelope(
        module="DNA_3D",
        payload={
            "kind": KIND_ILLUSTRATIVE,
            "n_atoms": n_atoms,
            "n_pairs": len(pairs),
            "view_start": model.get("view_start"),
            "view_end": model.get("view_end"),
        },
        status="ILLUSTRATIVE",
        algorithm=str(model.get("algorithm") or "canonical B-DNA helix illustration"),
        parameters=dict(model.get("parameters") or {}),
        source="HelixScope canonical helix illustration",
        sequence=seq,
    )
    return meta | {
        "status": "ILLUSTRATIVE",
        "kind": KIND_ILLUSTRATIVE,
        "kind_label": "Illustrative helical geometry (not a physical structure)",
        "layer": LAYER_ATOMIC_DETAIL,
        "sequence_hash": model.get("sequence_hash"),
        "envelope": envelope,
        "model": model,
        "partial_view": bool(model.get("partial_view")),
        "view_start": model.get("view_start"),
        "view_end": model.get("view_end"),
        "n_strands": 2,
        "n_pairs": len(pairs),
        "n_residues_drawn": len(residues),
        "n_atoms": n_atoms,
        "n_scene_objects": n_atoms + len(pairs),
        "max_base_labels": MAX_BASE_LABELS,
        "object_budget": budget,
        "fragment_bases": bases,
        "fragment": fragment,
        "disclaimer": model.get("disclaimer"),
        "algorithm": model.get("algorithm") or meta.get("algorithm"),
        "parameters": model.get("parameters") or model.get("provenance", {}).get("parameters"),
        "source": "HelixScope canonical helix illustration",
        "provenance": model.get("provenance") or meta.get("provenance"),
    }


def _unavailable(status: str, kind: str, reason: str, sequence: str) -> dict:
    digest = provenance.sequence_digest(sequence) if sequence else ""
    return {
        "status": status,
        "kind": kind,
        "kind_label": "Coordinates Unavailable",
        "layer": LAYER_ATOMIC_DETAIL,
        "reason": reason,
        "sequence_hash": digest,
        "envelope": None,
        "model": None,
        "partial_view": False,
        "n_strands": 0,
        "n_pairs": 0,
        "n_atoms": 0,
        "n_scene_objects": 0,
        "disclaimer": reason,
    }


def _unmapped_guide(reason: str, seq: str) -> dict:
    return {
        "status": "UNMAPPED",
        "kind": KIND_UNAVAILABLE,
        "indices": [],
        "reason": reason,
        "sequence_hash": provenance.sequence_digest(seq),
    }


def _pam_matches(window: str, pam: str) -> bool:
    if len(window) != len(pam):
        return False
    for observed, expected in zip(window, pam):
        if expected == "N":
            continue
        if expected == "G" and observed != "G":
            return False
        if expected == "A" and observed != "A":
            return False
        if expected == "C" and observed != "C":
            return False
        if expected == "T" and observed != "T":
            return False
        if expected not in "NACGT":
            return False
        if expected in "ACGT" and observed != expected:
            return False
    return True
