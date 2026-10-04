"""Catalogo de estruturas experimentais Cas-guia-DNA e contrato 3D futuro.

Nao descarrega coordenadas automaticamente. Os PDB listados sao
identificadores publicados. O ficheiro mmCIF e obtido da RCSB Files API
somente apos acao explicita; HelixScope nao embute XYZ nem monta Cas9
sintetico.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional

from . import provenance

EXPERIMENTAL_CAS9_COMPLEXES: tuple[dict, ...] = (
    {
        "pdb_id": "4UN3",
        "title": "SpCas9-sgRNA-DNA (Anders et al. 2014)",
        "kind": "experimental",
        "organism": "Streptococcus pyogenes Cas9",
        "source": "RCSB PDB 4UN3",
        "notes": (
            "Deposited crystal structure of a Cas9-sgRNA-DNA complex. "
            "Coordinates are not embedded. Fetch from RCSB Files API is an "
            "explicit user action."
        ),
        "helixscope_3d": "FETCH_ON_DEMAND",
    },
    {
        "pdb_id": "4OO8",
        "title": "SpCas9-sgRNA-DNA (Nishimasu et al. 2014)",
        "kind": "experimental",
        "organism": "Streptococcus pyogenes Cas9",
        "source": "RCSB PDB 4OO8",
        "notes": (
            "Deposited crystal structure. Coordinates are not embedded. "
            "Fetch from RCSB Files API is an explicit user action."
        ),
        "helixscope_3d": "FETCH_ON_DEMAND",
    },
)
"""Identificadores PDB reais. Nenhuma coordenada e embutida ou inventada."""


def list_experimental_complexes() -> list[dict]:
    """Copia do catalogo de ponteiros PDB.

    Args:
        Nenhum.

    Returns:
        Lista de dicts sem coordenadas.

    Raises:
        Nenhum.
    """
    return [dict(item) for item in EXPERIMENTAL_CAS9_COMPLEXES]


def mapping_contract(
    *,
    guide: Optional[Mapping[str, Any]] = None,
    hit: Optional[Mapping[str, Any]] = None,
    reference: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Contrato de mapeamento CRISPR -> estrutura futura.

    Args:
        guide: Guia SpCas9, se houver.
        hit: Off-target ou on-target verificado, se houver.
        reference: Referencia da busca, se houver.

    Returns:
        Dict com sequencias e coordenadas quando existirem, e
        structure_mapping UNAVAILABLE. Nunca inclui xyz fabricados.

    Raises:
        Nenhum.
    """
    guide = dict(guide or {})
    hit = dict(hit or {})
    reference = dict(reference or {})
    return provenance.json_safe(
        {
            "status": "UNAVAILABLE",
            "kind": "unavailable",
            "helixscope_3d": "FETCH_ON_DEMAND",
            "reason": (
                "HelixScope does not fabricate a Cas9-guide-DNA complex. "
                "Deposited mmCIF files (catalog PDB IDs) can be fetched from "
                "RCSB on an explicit action. Guide/target/PAM highlights require "
                "an exact subsequence match on the deposited polymer. Off-target "
                "genomic coordinates are not a 3D structure."
            ),
            "guide_sequence": guide.get("guide_sequence") or hit.get("guide_sequence") or "",
            "target_sequence": hit.get("target_sequence") or guide.get("guide_sequence") or "",
            "PAM": hit.get("PAM") or guide.get("pam_sequence") or "",
            "strand": hit.get("strand") or guide.get("strand") or "",
            "chromosome_or_contig": hit.get("chromosome_or_contig") or "",
            "start_0based": hit.get("start_0based"),
            "end_0based": hit.get("end_0based"),
            "coordinate_system": hit.get("coordinate_system") or "",
            "reference_identity_hash": (
                hit.get("reference_identity_hash")
                or reference.get("identity_hash")
                or ""
            ),
            "assembly_declared": (
                hit.get("assembly_declared") or reference.get("assembly_declared") or ""
            ),
            "verified_assembly": reference.get("verified_assembly") or "",
            "experimental_complex_pointers": [
                item["pdb_id"] for item in EXPERIMENTAL_CAS9_COMPLEXES
            ],
            "structure_mapping": None,
        }
    )
