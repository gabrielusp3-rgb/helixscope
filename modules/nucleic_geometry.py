"""Geometria 3D ilustrativa de acidos nucleicos. Nao e estrutura experimental.

B-DNA e A-RNA usam parametros helicoidais canonicos publicados. A sequencia
real controla a identidade das bases e a ordem. As coordenadas XYZ nao sao
a conformacao fisica da molecula no celular.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional

from . import dna_analysis, provenance, region_nav

KIND_ILLUSTRATIVE: str = "illustrative"
MAX_ILLUSTRATIVE_NT: int = 400
B_DNA_RISE_A: float = 3.38
B_DNA_TWIST_DEG: float = 36.0
B_DNA_RADIUS_A: float = 9.0
A_RNA_RISE_A: float = 2.8
A_RNA_TWIST_DEG: float = 32.7
A_RNA_RADIUS_A: float = 8.7
COMPLEMENT_DNA: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "U": "A",
}
COMPLEMENT_RNA: Dict[str, str] = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "T": "A",
}


class NucleicGeometryError(Exception):
    """Falha classificada da geometria ilustrativa.

    Attributes:
        category: INVALID_INPUT, RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def illustrative_bdna_model(
    sequence: str,
    *,
    highlight_indices: Optional[List[int]] = None,
    region_start: int = 0,
    region_end: Optional[int] = None,
    lod: str = region_nav.LOD_HIGH,
    center: Optional[int] = None,
) -> dict:
    """Helice B-DNA canonica controlada pela sequencia. ILLUSTRATIVE.

    Args:
        sequence: DNA na ordem 5'->3' da fita plus.
        highlight_indices: Posicoes 0-based da sequencia completa a marcar.
        region_start: Inicio da janela visual (0-based).
        region_end: Fim exclusivo; None = pedido ate ao fim, depois limitado pelo LOD.
        lod: low/medium/high (comprimento maximo da janela, nao geometria inventada).
        center: Se definido, centra a janela neste indice.

    Returns:
        Envelope com atoms (P de cada fita), residues, kind=illustrative,
        parameters e disclaimer. Complemento Watson-Crick da fita minus.

    Raises:
        NucleicGeometryError: INVALID_INPUT ou RESOURCE_LIMIT.

    Nota biologica:
        Implementacao heuristica simplificada inspirada em parametros
        helicoidais B-DNA (rise 3.38 A, 10 bp/volta, raio ~9 A);
        nao reproduz o modelo/algoritmo original publicado de um
        refinador cristalografico nem a conformacao celular. Nao e
        EXPERIMENTAL. Nao e estrutura fisica real desta molecula.
        Sequencias longas usam uma janela contigua; bases omitidas nao
        sao interpoladas na helice.
    """
    info = dna_analysis.validate_for_molecule(sequence, "DNA")
    if not info["is_valid"]:
        raise NucleicGeometryError(
            str(info.get("rejection_reason") or "Sequence is not DNA."),
            "INVALID_INPUT",
        )
    seq = str(info["sequence"])
    return _helix_model(
        seq,
        molecule="DNA",
        rise=B_DNA_RISE_A,
        twist_deg=B_DNA_TWIST_DEG,
        radius=B_DNA_RADIUS_A,
        complement_map=COMPLEMENT_DNA,
        form_name="canonical B-DNA helix",
        highlight_indices=highlight_indices or [],
        region_start=region_start,
        region_end=region_end,
        lod=lod,
        center=center,
    )


def illustrative_arna_model(
    sequence: str,
    *,
    highlight_indices: Optional[List[int]] = None,
    region_start: int = 0,
    region_end: Optional[int] = None,
    lod: str = region_nav.LOD_HIGH,
    center: Optional[int] = None,
) -> dict:
    """Helice A-RNA linear canonica. ILLUSTRATIVE, nao o fold MFE.

    Args:
        sequence: RNA 5'->3'.
        highlight_indices: Posicoes 0-based a marcar.

    Returns:
        Envelope ilustrativo. Nao usa ViennaRNA e nao converte dot-bracket.

    Raises:
        NucleicGeometryError.

    Nota biologica:
        Implementacao heuristica simplificada inspirada em parametros
        A-RNA (rise ~2.8 A, ~11 nt/volta); nao reproduz o modelo/algoritmo
        original publicado de previsao 3D (FARFAR, RNAComposer, SimRNA).
        MFE/dot-bracket nao geram estas coordenadas.
    """
    info = dna_analysis.validate_for_molecule(sequence, "RNA")
    if not info["is_valid"]:
        raise NucleicGeometryError(
            str(info.get("rejection_reason") or "Sequence is not RNA."),
            "INVALID_INPUT",
        )
    seq = str(info["sequence"])
    return _helix_model(
        seq,
        molecule="RNA",
        rise=A_RNA_RISE_A,
        twist_deg=A_RNA_TWIST_DEG,
        radius=A_RNA_RADIUS_A,
        complement_map=COMPLEMENT_RNA,
        form_name="canonical linear A-RNA helix",
        highlight_indices=highlight_indices or [],
        region_start=region_start,
        region_end=region_end,
        lod=lod,
        center=center,
    )


def availability() -> dict:
    """Declara o que esta geometria e e o que nao e.

    Args:
        Nenhum.

    Returns:
        Dict illustrative_available, experimental_requires_pdb.

    Raises:
        Nenhum.
    """
    return {
        "illustrative_bdna": True,
        "illustrative_arna_linear": True,
        "experimental": False,
        "predicted_3d_engine": False,
        "mfe_does_not_yield_3d": True,
        "reason": (
            "Helical illustrations use published canonical parameters and the "
            "pasted sequence. Experimental DNA/RNA 3D requires a deposited "
            "structure. ViennaRNA MFE is secondary structure only."
        ),
    }


def _helix_model(
    seq: str,
    *,
    molecule: str,
    rise: float,
    twist_deg: float,
    radius: float,
    complement_map: Mapping[str, str],
    form_name: str,
    highlight_indices: List[int],
    region_start: int = 0,
    region_end: Optional[int] = None,
    lod: str = region_nav.LOD_HIGH,
    center: Optional[int] = None,
) -> dict:
    if len(seq) < 2:
        raise NucleicGeometryError(
            "An illustrative helix needs at least two nucleotides.",
            "INVALID_INPUT",
        )
    try:
        view = region_nav.contiguous_view(
            len(seq),
            start=region_start,
            end=region_end,
            lod=lod,
            center=center,
        )
    except region_nav.RegionError as exc:
        raise NucleicGeometryError(str(exc), exc.category) from exc
    start = int(view["start"])
    end = int(view["end"])
    fragment = seq[start:end]
    if len(fragment) < 2:
        raise NucleicGeometryError(
            "An illustrative helix needs at least two nucleotides in the window.",
            "INVALID_INPUT",
        )
    highlighted = {
        int(item)
        for item in highlight_indices
        if start <= int(item) < end
    }
    twist = math.radians(twist_deg)
    atoms: List[dict] = []
    residues: List[dict] = []
    pairs: List[dict] = []
    atom_id = 1
    for local, base in enumerate(fragment):
        index = start + local
        angle = local * twist
        plus = {
            "x": radius * math.cos(angle),
            "y": radius * math.sin(angle),
            "z": local * rise,
        }
        minus = {
            "x": radius * math.cos(angle + math.pi),
            "y": radius * math.sin(angle + math.pi),
            "z": local * rise,
        }
        partner = complement_map.get(base, "N")
        atoms.append(
            _atom(
                atom_id,
                "P",
                "P",
                base,
                "A",
                index + 1,
                plus,
                molecule,
            )
        )
        atom_id += 1
        atoms.append(
            _atom(
                atom_id,
                "P",
                "P",
                partner,
                "B",
                index + 1,
                minus,
                molecule,
            )
        )
        atom_id += 1
        residues.append(
            {
                "index_0based": index,
                "local_index_0based": local,
                "base": base,
                "complement": partner,
                "chain_id": "A",
                "auth_seq_id": index + 1,
                "highlighted": index in highlighted,
                "p": plus,
            }
        )
        pairs.append(
            {
                "plus_index_0based": index,
                "plus_base": base,
                "minus_base": partner,
                "plus": dict(plus),
                "minus": dict(minus),
                "kind": KIND_ILLUSTRATIVE,
            }
        )
    digest = provenance.sequence_digest(seq)
    partial = bool(view["partial_view"])
    return provenance.analysis_envelope(
        module="NUCLEIC_GEOMETRY",
        payload={
            "kind": KIND_ILLUSTRATIVE,
            "molecule": molecule,
            "n_residues": len(fragment),
            "n_atoms": len(atoms),
            "n_strands": 2,
            "n_pairs": len(pairs),
            "view_start": start,
            "view_end": end,
            "lod": view["lod"],
        },
        status="ILLUSTRATIVE",
        algorithm=form_name,
        parameters={
            "rise_angstrom": rise,
            "twist_degrees": twist_deg,
            "radius_angstrom": radius,
            "max_nt": int(view["max_nt"]),
            "lod": view["lod"],
            "view_start": start,
            "view_end": end,
            "full_length": len(seq),
        },
        source="HelixScope canonical helix illustration",
        input_identifier=digest,
        sequence=seq,
    ) | {
        "status": "ILLUSTRATIVE",
        "kind": KIND_ILLUSTRATIVE,
        "kind_label": "Illustrative helical geometry (not a physical structure)",
        "molecule": molecule,
        "sequence": seq,
        "fragment": fragment,
        "sequence_hash": digest,
        "n_residues": len(fragment),
        "n_atoms": len(atoms),
        "n_strands": 2,
        "n_pairs": len(pairs),
        "atoms": atoms,
        "residues": residues,
        "pairs": pairs,
        "view_start": start,
        "view_end": end,
        "lod": view["lod"],
        "partial_view": partial,
        "view_scope": "Partial structure view" if partial else "Full sequence window",
        "chains": [
            {
                "chain_id": "A",
                "role": "plus_strand",
                "sequence": fragment,
                "molecule": molecule,
            },
            {
                "chain_id": "B",
                "role": "minus_strand",
                "sequence": "".join(complement_map.get(base, "N") for base in fragment),
                "molecule": molecule,
            },
        ],
        "highlight_indices": sorted(highlighted),
        "branch_lengths_present": False,
        "disclaimer": (
            "ILLUSTRATIVE canonical helix. The sequence order and base identity "
            "in this window are real. The XYZ coordinates are not an experimental "
            "structure, not a predicted fold, and not the cellular conformation. "
            "ViennaRNA MFE does not produce these coordinates. "
            + (
                f"PARTIAL VIEW: residues {start}–{end} of {len(seq)} "
                f"(LOD {view['lod']}, cap {view['max_nt']} nt). "
                "Omitted sequence is not drawn and not interpolated."
                if partial
                else f"Window covers the full sequence (LOD {view['lod']})."
            )
        ),
    }


def _atom(
    atom_id: int,
    name: str,
    element: str,
    comp: str,
    chain: str,
    seq_id: int,
    xyz: Mapping[str, float],
    molecule: str,
) -> dict:
    return {
        "atom_id": atom_id,
        "atom_name": name,
        "element": element,
        "comp_id": comp,
        "auth_asym_id": chain,
        "label_asym_id": chain,
        "label_seq_id": seq_id,
        "auth_seq_id": seq_id,
        "x": float(xyz["x"]),
        "y": float(xyz["y"]),
        "z": float(xyz["z"]),
        "occupancy": 1.0,
        "b_iso": None,
        "model": 1,
        "group": "ATOM",
        "molecule": molecule,
        "kind": KIND_ILLUSTRATIVE,
    }


def structure_3d_input_from_illustrative(model: Mapping[str, Any]) -> dict:
    """Envelope 3D a partir da helice ilustrativa. Nao muda o kind.

    Args:
        model: Saida de illustrative_bdna_model / illustrative_arna_model.

    Returns:
        Contrato compativel com structure_scene.build_scene. kind=illustrative.

    Raises:
        NucleicGeometryError: INVALID_INPUT se o modelo nao for ilustrativo.
    """
    if str(model.get("kind") or "") != KIND_ILLUSTRATIVE:
        raise NucleicGeometryError(
            "Only an ILLUSTRATIVE nucleic model can enter this envelope.",
            "INVALID_INPUT",
        )
    atoms = [dict(item) for item in list(model.get("atoms") or [])]
    mappings: List[dict] = []
    for residue in list(model.get("residues") or []):
        point = dict(residue.get("p") or {})
        has = all(point.get(axis) is not None for axis in ("x", "y", "z"))
        mappings.append(
            {
                "query_index_0based": residue.get("index_0based"),
                "query_residue": residue.get("base"),
                "chain_index_0based": residue.get("index_0based"),
                "chain_residue": residue.get("base"),
                "chain_id": "A",
                "label_seq_id": residue.get("auth_seq_id"),
                "auth_seq_id": residue.get("auth_seq_id"),
                "insertion_code": "",
                "comp_id": residue.get("base"),
                "has_coordinates": has,
                "coordinate_status": "available" if has else "unavailable",
                "ca": None,
                "p": point if has else None,
                "c1": None,
                "model": 1,
                "molecule_type": model.get("molecule"),
            }
        )
    minus_seq = ""
    for chain in list(model.get("chains") or []):
        if str(chain.get("chain_id") or "") == "B":
            minus_seq = str(chain.get("sequence") or "")
    for index, base in enumerate(minus_seq):
        atom = next(
            (
                item
                for item in atoms
                if str(item.get("auth_asym_id") or "") == "B"
                and int(item.get("auth_seq_id") or 0) == index + 1
            ),
            None,
        )
        point = None
        if atom is not None:
            point = {"x": atom.get("x"), "y": atom.get("y"), "z": atom.get("z")}
        has = bool(point) and all(point.get(axis) is not None for axis in ("x", "y", "z"))
        mappings.append(
            {
                "query_index_0based": None,
                "query_residue": None,
                "chain_index_0based": index,
                "chain_residue": base,
                "chain_id": "B",
                "label_seq_id": index + 1,
                "auth_seq_id": index + 1,
                "insertion_code": "",
                "comp_id": base,
                "has_coordinates": has,
                "coordinate_status": "available" if has else "unavailable",
                "ca": None,
                "p": point if has else None,
                "c1": None,
                "model": 1,
                "molecule_type": model.get("molecule"),
            }
        )
    n_with = sum(1 for item in mappings if item.get("has_coordinates") and item.get("query_index_0based") is not None)
    n_query = len(str(model.get("sequence") or ""))
    return {
        "source": model.get("source") or "HelixScope canonical helix illustration",
        "kind": KIND_ILLUSTRATIVE,
        "kind_label": str(model.get("kind_label") or "Illustrative Structure"),
        "structure_id": "",
        "sequence_hash": model.get("sequence_hash"),
        "structure_hash": model.get("sequence_hash"),
        "chain_id": "A",
        "model_number": 1,
        "models": [1],
        "chains": list(model.get("chains") or []),
        "n_atoms": len(atoms),
        "n_residues": n_query,
        "atoms": atoms,
        "residues": list(model.get("residues") or []),
        "coordinates": [
            {
                "atom": item.get("atom_name"),
                "element": item.get("element"),
                "residue": item.get("comp_id"),
                "chain": item.get("auth_asym_id"),
                "label_seq_id": item.get("label_seq_id"),
                "auth_seq_id": item.get("auth_seq_id"),
                "insertion_code": "",
                "alt_id": "",
                "x": item.get("x"),
                "y": item.get("y"),
                "z": item.get("z"),
                "occupancy": item.get("occupancy"),
                "b_iso": item.get("b_iso"),
                "model": item.get("model"),
                "group": item.get("group"),
            }
            for item in atoms
        ],
        "mappings": mappings,
        "mapping_status": "mapped",
        "mapping_status_note": model.get("disclaimer"),
        "coverage_query_with_coordinates": n_with / n_query if n_query else 0.0,
        "n_with_coordinates": n_with,
        "n_without_coordinates": max(0, n_query - n_with),
        "residues_without_coordinates": [],
        "confidence_data": {},
        "molecule_type": model.get("molecule"),
        "pairs": [dict(item) for item in list(model.get("pairs") or [])],
        "n_strands": int(model.get("n_strands") or 2),
        "n_pairs": len(list(model.get("pairs") or [])),
        "partial_view": bool(model.get("partial_view")),
        "view_start": model.get("view_start"),
        "view_end": model.get("view_end"),
        "metadata": {
            "method": model.get("algorithm") or "canonical helix illustration",
            "resolution_angstrom": None,
            "models": [1],
            "disclaimer": model.get("disclaimer"),
        },
        "provenance": dict(model.get("provenance") or {}),
        "renderer": "data-only",
        "status": "ILLUSTRATIVE",
        "disclaimer": model.get("disclaimer"),
    }
