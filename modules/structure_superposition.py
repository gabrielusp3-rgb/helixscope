"""Superposicao geometrica sobre copias de coordenadas.

O renderer e as estruturas parseadas originais nao sao mutadas. A transformacao
4x4 (column-major, convencao RCSB) e aplicada a uma copia dos atomos CA do alvo.

RMSD independente: distancia Euclidiana media dos pares alinhados APOS a
transformacao, calculada aqui com math.sqrt, sem reutilizar o valor da API.
Um teste deve comparar este numero com o RMSD reportado pela RCSB, nao com
uma funcao que apenas ecoa o JSON.

Alinhamentos flexiveis com mais de um bloco: a vista 3D e PARTIAL (bloco 0)
ou UNAVAILABLE; nunca se aplica uma unica matriz e se afirma que representa
o alinhamento flexivel completo.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import copy
import math
from typing import Any, Mapping, Optional, Sequence

from . import provenance, structure_alignment

ATOM_CA: str = "CA"


class SuperpositionError(RuntimeError):
    """Falha ao aplicar ou validar uma transformacao.

    Attributes:
        category: INVALID_INPUT, UNMAPPED, PARTIAL ou ERROR.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper()


def apply_column_major_4x4(
    x: float,
    y: float,
    z: float,
    matrix: Sequence[float],
) -> tuple[float, float, float]:
    """Aplica M * [x, y, z, 1] com M 4x4 column-major (indice j*4+i).

    Args:
        x: Coordenada X original.
        y: Coordenada Y original.
        z: Coordenada Z original.
        matrix: 16 floats.

    Returns:
        (x', y', z') na mesma unidade (angstrom).

    Raises:
        SuperpositionError: INVALID_INPUT se a matriz nao tiver 16 finitos.

    Nota:
        Implementacao independente da API RCSB: so a convencao de armazenamento
        e a documentada. Nao chama o parser de alinhamento.
    """
    if len(matrix) != 16:
        raise SuperpositionError("Transformation must have 16 entries.", "INVALID_INPUT")
    try:
        m = [float(item) for item in matrix]
    except (TypeError, ValueError) as exc:
        raise SuperpositionError("Transformation is not numeric.", "INVALID_INPUT") from exc
    if any(not math.isfinite(item) for item in m):
        raise SuperpositionError("Transformation contains a non-finite number.", "INVALID_INPUT")
    xp = m[0] * x + m[4] * y + m[8] * z + m[12]
    yp = m[1] * x + m[5] * y + m[9] * z + m[13]
    zp = m[2] * x + m[6] * y + m[10] * z + m[14]
    return (xp, yp, zp)


def ca_index(parsed: Mapping[str, Any], *, prefer_label: bool = True) -> dict[tuple[str, int], dict]:
    """Indice (chain, label_seq_id) -> atomo CA (primeira altloc).

    Args:
        parsed: Estrutura parseada (atoms).
        prefer_label: Se True, a chave usa label_asym_id; senao auth_asym_id.

    Returns:
        Dict. Residuos sem CA ficam de fora (UNMAPPED no RMSD, nao fabricados).

    Raises:
        SuperpositionError: INVALID_INPUT se atoms faltar.
    """
    atoms = parsed.get("atoms")
    if not isinstance(atoms, Sequence):
        raise SuperpositionError("Parsed structure has no atoms list.", "INVALID_INPUT")
    index: dict[tuple[str, int], dict] = {}
    for atom in atoms:
        if not isinstance(atom, Mapping):
            continue
        if str(atom.get("group") or "") != "ATOM":
            continue
        if str(atom.get("atom_name") or "").strip().upper() != ATOM_CA:
            continue
        seq = atom.get("label_seq_id")
        if seq is None:
            continue
        chain = str(
            atom.get("label_asym_id") if prefer_label else atom.get("auth_asym_id")
            or atom.get("label_asym_id")
            or atom.get("auth_asym_id")
            or ""
        )
        key = (chain, int(seq))
        if key in index:
            continue
        if atom.get("x") is None:
            continue
        index[key] = dict(atom)
    return index


def transform_atom_copies(
    atoms: Sequence[Mapping[str, Any]],
    matrix: Sequence[float],
) -> list[dict]:
    """Copia ATOM/HETATM e aplica a matriz. Originais nao sao tocados.

    Args:
        atoms: Lista de atomos.
        matrix: 4x4 column-major.

    Returns:
        Nova lista com x,y,z transformados e flag transformed=True.

    Raises:
        SuperpositionError: INVALID_INPUT.
    """
    copies: list[dict] = []
    for atom in atoms:
        if not isinstance(atom, Mapping):
            continue
        row = copy.deepcopy(dict(atom))
        try:
            x = float(row["x"])
            y = float(row["y"])
            z = float(row["z"])
        except (TypeError, ValueError, KeyError):
            copies.append(row)
            continue
        xp, yp, zp = apply_column_major_4x4(x, y, z, matrix)
        row["x"] = xp
        row["y"] = yp
        row["z"] = zp
        row["transformed"] = True
        row["original_x"] = x
        row["original_y"] = y
        row["original_z"] = z
        copies.append(row)
    return copies


def independent_aligned_ca_rmsd(
    *,
    reference_parsed: Mapping[str, Any],
    target_parsed: Mapping[str, Any],
    pairs: Sequence[Mapping[str, Any]],
    target_matrix: Sequence[float],
) -> dict:
    """RMSD CA dos pares alinhados apos transformar o alvo. COMPUTED.

    Args:
        reference_parsed: Estrutura referencia (coordenadas originais).
        target_parsed: Estrutura alvo (coordenadas originais).
        pairs: residue_pairs do objeto de alinhamento.
        target_matrix: Transformacao 4x4 do alvo (bloco).

    Returns:
        Dict status, rmsd_angstrom, n_pairs_used, n_pairs_missing, method.

    Raises:
        SuperpositionError: INVALID_INPUT se nao houver pares.

    Nota biologica:
        Este RMSD nao e o valor da API. E um calculo local sobre as coordenadas
        depositadas desta maquina, para validar a matriz. Pares sem CA nas duas
        estruturas sao contados como missing, nao como distancia 0.
    """
    if not pairs:
        raise SuperpositionError("No aligned pairs were supplied for RMSD.", "INVALID_INPUT")
    ref_index = ca_index(reference_parsed)
    tgt_index = ca_index(target_parsed)
    squares: list[float] = []
    missing = 0
    used_pairs: list[dict] = []
    for pair in pairs:
        if not isinstance(pair, Mapping):
            missing += 1
            continue
        ref_key = (str(pair.get("reference_asym_id") or ""), int(pair["reference_label_seq_id"]))
        tgt_key = (str(pair.get("target_asym_id") or ""), int(pair["target_label_seq_id"]))
        ref_atom = ref_index.get(ref_key)
        tgt_atom = tgt_index.get(tgt_key)
        if ref_atom is None or tgt_atom is None:
            missing += 1
            continue
        xt, yt, zt = apply_column_major_4x4(
            float(tgt_atom["x"]),
            float(tgt_atom["y"]),
            float(tgt_atom["z"]),
            target_matrix,
        )
        dx = float(ref_atom["x"]) - xt
        dy = float(ref_atom["y"]) - yt
        dz = float(ref_atom["z"]) - zt
        squares.append(dx * dx + dy * dy + dz * dz)
        used_pairs.append(
            {
                "reference": f"{ref_key[0]}:{ref_key[1]}",
                "target": f"{tgt_key[0]}:{tgt_key[1]}",
                "distance_angstrom": math.sqrt(dx * dx + dy * dy + dz * dz),
            }
        )
    if not squares:
        return {
            "status": "UNMAPPED",
            "rmsd_angstrom": None,
            "n_pairs_used": 0,
            "n_pairs_missing": missing,
            "n_pairs_declared": len(pairs),
            "atoms": ATOM_CA,
            "reason": (
                "None of the aligned label_seq_id pairs had C-alpha coordinates "
                "in both parsed structures. RMSD is N/A, not 0."
            ),
            "method": "independent_euclidean_rmsd_after_target_transform",
        }
    rmsd = math.sqrt(sum(squares) / len(squares))
    return {
        "status": "COMPUTED",
        "rmsd_angstrom": rmsd,
        "n_pairs_used": len(squares),
        "n_pairs_missing": missing,
        "n_pairs_declared": len(pairs),
        "atoms": ATOM_CA,
        "method": (
            "sqrt(mean(squared Euclidean CA-CA distance)) after applying the "
            "RCSB 4x4 column-major transform to a copy of the target CA atoms. "
            "Not the API summary value; compared to it for validation."
        ),
        "evidence_status": "COMPUTED",
        "software_version": provenance.HELIXSCOPE_VERSION,
        "pair_distances_angstrom_head": used_pairs[:12],
    }


def superposition_bundle(
    *,
    alignment: Mapping[str, Any],
    reference_parsed: Mapping[str, Any],
    target_parsed: Mapping[str, Any],
    block_index: int = 0,
) -> dict:
    """Copia transformada do alvo + RMSD independente + estado visual.

    Args:
        alignment: Objeto parse_alignment_payload.
        reference_parsed: Estrutura referencia original.
        target_parsed: Estrutura alvo original.
        block_index: Bloco a visualizar.

    Returns:
        Dict visual_status, transformed_target_atoms, independent_rmsd,
        reference_atoms (copia rasa nao transformada).

    Raises:
        SuperpositionError: INVALID_INPUT se o bloco nao existir.
    """
    blocks = alignment.get("blocks") or []
    if not isinstance(blocks, Sequence) or block_index < 0 or block_index >= len(blocks):
        raise SuperpositionError("Requested superposition block does not exist.", "INVALID_INPUT")
    block = blocks[block_index]
    kind = str(alignment.get("method_kind") or "")
    visual = str(alignment.get("superposition_visual") or "AVAILABLE")
    if kind == "flexible" and len(blocks) > 1:
        visual = "PARTIAL"
    matrix = block.get("target_transform_column_major_4x4") or []
    target_atoms = list(target_parsed.get("atoms") or [])
    transformed = transform_atom_copies(target_atoms, matrix)
    ref_atoms = copy.deepcopy(list(reference_parsed.get("atoms") or []))
    rmsd = independent_aligned_ca_rmsd(
        reference_parsed=reference_parsed,
        target_parsed=target_parsed,
        pairs=block.get("pairs") or alignment.get("residue_pairs") or [],
        target_matrix=matrix,
    )
    return {
        "status": "COMPUTED",
        "visual_status": visual,
        "block_index": block_index,
        "n_blocks": len(blocks),
        "method": alignment.get("method_display") or alignment.get("method"),
        "method_kind": kind,
        "atoms_fitted": structure_alignment.ATOMS_USED,
        "reference_entry": (alignment.get("reference") or {}).get("entry_id"),
        "target_entry": (alignment.get("target") or {}).get("entry_id"),
        "reference_kind": (alignment.get("reference") or {}).get("kind_label"),
        "target_kind": (alignment.get("target") or {}).get("kind_label"),
        "reference_atoms": ref_atoms,
        "transformed_target_atoms": transformed,
        "original_target_untouched": True,
        "independent_rmsd": rmsd,
        "api_rmsd_block_angstrom": block.get("rmsd_angstrom"),
        "api_rmsd_global_angstrom": alignment.get("rmsd_global_angstrom"),
        "n_aligned_residue_pairs": alignment.get("n_aligned_residue_pairs"),
        "aln_coverage_percent": alignment.get("aln_coverage_percent"),
        "tm_scores": alignment.get("tm_scores"),
        "note": (
            "Reference coordinates are a copy of the parsed structure. Target "
            "coordinates in this bundle are a transformed copy. The parsed "
            "inputs are not modified."
            + (
                " Flexible multi-block overlay is PARTIAL: only this block's "
                "transform is drawn."
                if visual == "PARTIAL"
                else ""
            )
        ),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def rmsd_agrees(
    independent: Mapping[str, Any],
    api_rmsd: Optional[float],
    *,
    abs_tolerance_angstrom: float = 0.15,
    rel_tolerance: float = 0.1,
) -> dict:
    """Compara RMSD local com o da API. Tolerancia documentada.

    Args:
        independent: Saida de independent_aligned_ca_rmsd.
        api_rmsd: RMSD do bloco ou summary.
        abs_tolerance_angstrom: Max(|a-b|) aceite em A.
        rel_tolerance: Max |a-b|/max(a,b) quando ambos > 0.

    Returns:
        Dict agree, delta, tolerance.

    Raises:
        Nenhum.
    """
    local = independent.get("rmsd_angstrom")
    if local is None or api_rmsd is None:
        return {
            "comparable": False,
            "agree": None,
            "reason": "One or both RMSD values are N/A; agreement is not claimed.",
            "abs_tolerance_angstrom": abs_tolerance_angstrom,
            "rel_tolerance": rel_tolerance,
        }
    delta = abs(float(local) - float(api_rmsd))
    scale = max(abs(float(local)), abs(float(api_rmsd)), 1e-9)
    agree = delta <= abs_tolerance_angstrom or (delta / scale) <= rel_tolerance
    return {
        "comparable": True,
        "agree": agree,
        "local_rmsd_angstrom": float(local),
        "api_rmsd_angstrom": float(api_rmsd),
        "delta_angstrom": delta,
        "abs_tolerance_angstrom": abs_tolerance_angstrom,
        "rel_tolerance": rel_tolerance,
        "n_pairs_used": independent.get("n_pairs_used"),
        "note": (
            "Tolerances cover coordinate rounding, altloc choice and the API "
            "using the first model. A disagreement larger than both tolerances "
            "is a real mismatch, not silently ignored."
        ),
    }
