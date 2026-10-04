"""Catalogo de estruturas experimentais de DNA/RNA depositadas.

Ponteiros RCSB PDB reais. Coordenadas nao sao inventadas: o mmCIF e parseado
do ficheiro (fixture publica ou fetch allowlisted). Nenhuma funcao importa
Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional

from . import provenance

EXPERIMENTAL_DNA_STRUCTURES: tuple[dict, ...] = (
    {
        "pdb_id": "1BNA",
        "title": "Dickerson-Drew B-DNA dodecamer (Drew et al. 1981)",
        "kind": "experimental",
        "molecule": "DNA",
        "source": "RCSB PDB 1BNA",
        "polymer_sequence": "CGCGAATTCGCG",
        "notes": (
            "Public X-ray structure of the self-complementary dodecamer "
            "CGCGAATTCGCG. Coordinates come from the mmCIF, not from a "
            "procedural helix."
        ),
        "helixscope_3d": "FIXTURE_OR_FETCH",
    },
)
"""Identificadores PDB reais de DNA. Sem XYZ embutidos neste modulo."""

EXPERIMENTAL_RNA_STRUCTURES: tuple[dict, ...] = (
    {
        "pdb_id": "1RNA",
        "title": "[U(UA)6A]2 RNA duplex (Dock-Bregeon et al. 1988/1991)",
        "kind": "experimental",
        "molecule": "RNA",
        "source": "RCSB PDB 1RNA",
        "polymer_sequence": "UUAUAUAUAUAUAA",
        "notes": (
            "Public X-ray RNA duplex. ViennaRNA MFE is not these coordinates."
        ),
        "helixscope_3d": "FIXTURE_OR_FETCH",
    },
)
"""Identificadores PDB reais de RNA. Sem XYZ embutidos neste modulo."""


def list_experimental_dna_structures() -> list[dict]:
    """Copia do catalogo DNA sem coordenadas.

    Args:
        Nenhum.

    Returns:
        Lista de dicts ponteiro.

    Raises:
        Nenhum.
    """
    return [dict(item) for item in EXPERIMENTAL_DNA_STRUCTURES]


def list_experimental_rna_structures() -> list[dict]:
    """Copia do catalogo RNA sem coordenadas.

    Args:
        Nenhum.

    Returns:
        Lista de dicts ponteiro.

    Raises:
        Nenhum.
    """
    return [dict(item) for item in EXPERIMENTAL_RNA_STRUCTURES]


def catalog_entry(pdb_id: str) -> Optional[dict]:
    """Localiza um ponteiro DNA ou RNA pelo PDB ID.

    Args:
        pdb_id: Identificador de 4 caracteres.

    Returns:
        Copia do entry ou None. None nao e uma estrutura.

    Raises:
        Nenhum.
    """
    key = str(pdb_id or "").strip().upper()
    for item in (*EXPERIMENTAL_DNA_STRUCTURES, *EXPERIMENTAL_RNA_STRUCTURES):
        if str(item.get("pdb_id") or "").upper() == key:
            return dict(item)
    return None


def rnacentral_is_appropriate_for_3d() -> dict:
    """Declara por que RNAcentral nao alimenta coordenadas 3D neste app.

    Args:
        Nenhum.

    Returns:
        Dict appropriate=False, reason.

    Raises:
        Nenhum.

    Nota biologica:
        RNAcentral agrega sequencias de RNA e IDs. Nao e um banco de
        coordenadas cartesianas. Usar um URS como se fosse um mmCIF seria
        fabricacao estrutural.
    """
    return provenance.analysis_envelope(
        module="NA_STRUCTURE_CATALOG",
        payload={
            "appropriate_for_3d_coordinates": False,
            "may_validate_sequence_identity": True,
            "reason": (
                "RNAcentral provides sequence records and cross-references, "
                "not deposited Cartesian coordinates. HelixScope uses RCSB "
                "PDB / PDBe mmCIF for experimental RNA 3D."
            ),
        },
        status="UNAVAILABLE",
        algorithm="source-scope declaration",
        parameters={},
        source="HelixScope NA structure catalog",
    )


def rna_3d_predictor_research() -> list[dict]:
    """Resultado documentado da investigacao de preditores 3D de RNA.

    Args:
        Nenhum.

    Returns:
        Uma linha por ferramenta considerada. Nenhuma e integrada.

    Raises:
        Nenhum.
    """
    return [
        {
            "tool": "RhoFold+",
            "source": "https://github.com/ml4bio/RhoFold",
            "status": "active (GitHub ml4bio/RhoFold)",
            "maintenance": "documented Linux conda env; README states MacOS is not supported",
            "model": "RNA-FM language model + 3D module (Shen et al., Nature Methods)",
            "version": "not installed in this environment (binary/checkpoint absent)",
            "license": "Apache-2.0 (source); training data on Hugging Face described as non-commercial",
            "reproducibility": "requires pretrained checkpoint and optional MSA construction",
            "windows_compatibility": "not supported natively (Linux env file only)",
            "decision": "not integrated",
            "why": (
                "Windows is the HelixScope runtime; RhoFold documents Linux-only "
                "local install. Integrating would add a giant DL stack and a "
                "checkpoint HelixScope does not vendor."
            ),
        },
        {
            "tool": "trRosettaRNA / trRosettaRNA2",
            "source": "https://github.com/YangLab-SDU/trRosettaRNA2 ; https://yanglab.qd.sdu.edu.cn/trRosettaRNA/",
            "status": "active (YangLab; Nature Protocols 2026 server paper)",
            "maintenance": "GitHub YangLab-SDU/trRosettaRNA2 latest release noted as v2.0.4 (2025-08-08)",
            "model": "transformer 2D geometries plus optional PyRosetta minimization",
            "version": "not installed in this environment",
            "license": "Apache-2.0 on the trRosettaRNA2 repository",
            "reproducibility": "weights plus optional Rosetta energy minimization",
            "windows_compatibility": "not verificable as a first-class Windows binary in this session",
            "decision": "not integrated",
            "why": (
                "Local use pulls a deep-learning protocol and optional PyRosetta. "
                "HelixScope does not bundle those engines (no dependency bloat)."
            ),
        },
        {
            "tool": "FARFAR2 (Rosetta RNA)",
            "source": "https://rosie.rosettacommons.org/farfar2",
            "status": "academic Rosetta/ROSIE server",
            "maintenance": "Rosetta Commons; not a small standalone CLI",
            "model": "fragment assembly of RNA with Rosetta energy function",
            "version": "not installed (Rosetta license/binary absent)",
            "license": "Rosetta academic license (not a drop-in open CLI)",
            "reproducibility": "requires Rosetta installation and database",
            "windows_compatibility": "not verificable; typical deployments are Linux",
            "decision": "not integrated",
            "why": (
                "FARFAR2 is not a lightweight library. Shipping Rosetta would "
                "violate the no-giant-engine rule and is not licensed for bundling."
            ),
        },
        {
            "tool": "RNAComposer",
            "source": "https://rnacomposer.cs.put.poznan.pl/",
            "status": "web server",
            "maintenance": "Poznan server (not a local library)",
            "model": "element-based 3D RNA modeling from secondary structure",
            "version": "server-side; not a HelixScope dependency",
            "license": "web service terms; not vendored",
            "reproducibility": "depends on remote server; not offline",
            "windows_compatibility": "N/A (HTTP service, not a Windows binary)",
            "decision": "not integrated",
            "why": (
                "RNAComposer is a remote server. HelixScope 3D must remain "
                "offline-capable and must not POST user sequences to unallowlisted "
                "hosts."
            ),
        },
    ]
