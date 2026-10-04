"""Estruturas de proteina recuperadas de fontes reais (RCSB PDB, AlphaFold DB, UniProt).

Nao inventa coordenadas, resolucao, pLDDT nem cobertura. Uma sequencia sozinha
nao implica estrutura. source e kind sao campos separados: um registo RETRIEVED
pode ser experimental (PDB) ou predicted (AlphaFold DB).

HTTP so para hosts allowlisted. Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import io
import json
import math
import os
import re
import secrets
import shutil
import socket
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
from types import MappingProxyType
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener, urlopen

from . import alignment, dna_analysis, engine_validation, provenance, scientific_checks, tool_detection, tool_paths

HTTP_TIMEOUT_S: float = 30.0
STRUCTURE_FILE_TIMEOUT_S: float = 90.0
DSSP_REMOTE_TIMEOUT_S: float = 60.0
MAX_METADATA_BYTES: int = 500_000
MAX_STRUCTURE_BYTES: int = 8_000_000
MAX_ATOMS: int = 50_000
MAX_RESIDUES: int = 5_000
MAX_SEARCH_HITS: int = 15
MAX_SEARCH_SEQUENCE: int = 2_000
MIN_SEARCH_SEQUENCE: int = 10
MAX_UNIPROT_XREFS: int = 40
DEFAULT_IDENTITY_CUTOFF: float = 0.9
DEFAULT_EVALUE_CUTOFF: float = 0.1
USER_AGENT: str = "HelixScope/protein-structure"

RCSB_DATA: str = "https://data.rcsb.org"
RCSB_SEARCH: str = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILES: str = "https://files.rcsb.org"
ALPHAFOLD_API: str = "https://alphafold.ebi.ac.uk"
UNIPROT_REST: str = "https://rest.uniprot.org"
PDB_REDO_DSSP_DO: str = "https://pdb-redo.eu/dssp/do"
PDB_REDO_DSSP_DB: str = "https://pdb-redo.eu/dssp/db/"
STRUCTURE_CACHE_ENV: str = "HELIXSCOPE_STRUCTURE_CACHE"
DSSP_REMOTE_TOOL: str = "DSSP PDB-REDO API"

ALLOWED_URLS: tuple[tuple[str, str], ...] = (
    ("data.rcsb.org", "/rest/v1/core/entry/"),
    ("data.rcsb.org", "/rest/v1/core/polymer_entity/"),
    ("search.rcsb.org", "/rcsbsearch/v2/query"),
    ("files.rcsb.org", "/download/"),
    ("alphafold.ebi.ac.uk", "/api/prediction/"),
    ("alphafold.ebi.ac.uk", "/files/"),
    ("rest.uniprot.org", "/uniprotkb/"),
    ("pdb-redo.eu", "/dssp/do"),
    ("pdb-redo.eu", "/dssp/db/"),
)
BUNDLED_MMCIF_DIR: Path = Path(__file__).resolve().parents[1] / "tests" / "fixtures"

PDB_ID_PATTERN: re.Pattern[str] = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
UNIPROT_PATTERN: re.Pattern[str] = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?$",
    re.IGNORECASE,
)
ALPHAFOLD_ID_PATTERN: re.Pattern[str] = re.compile(
    r"^AF-([A-Z0-9]+)-F[0-9]+$",
    re.IGNORECASE,
)

AA3_TO_1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}
NA3_TO_1: Dict[str, str] = {
    "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U",
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
    "ADE": "A", "CYT": "C", "GUA": "G", "THY": "T", "URA": "U",
}
DSSP_NAMES: tuple[str, ...] = ("mkdssp", "mkdssp.exe", "dssp", "dssp.exe")
DSSP_ENV: str = "HELIXSCOPE_DSSP"
STRIDE_NAMES: tuple[str, ...] = ("stride", "stride.exe")
STRIDE_ENV: str = "HELIXSCOPE_STRIDE"
EDTSURF_NAMES: tuple[str, ...] = ("edtsurf", "edtsurf.exe")
MSMS_NAMES: tuple[str, ...] = ("msms", "msms.exe")
SURFACE_NAMES: tuple[str, ...] = EDTSURF_NAMES + MSMS_NAMES
SURFACE_ENV: str = "HELIXSCOPE_SURFACE"
MSMS_ENV: str = "HELIXSCOPE_MSMS"
DSSP_TIMEOUT_S: float = 30.0
SURFACE_TIMEOUT_S: float = 45.0
MAX_SURFACE_VERTICES: int = 80_000
MAX_SURFACE_FACES: int = 160_000
EDTSURF_SURFACE_CODES: dict[int, str] = {
    1: "VDW",
    2: "SAS",
    3: "MS",
    4: "SES",
}
MAPPING_STATUSES: tuple[str, ...] = (
    "EXACT",
    "ALIGNED",
    "PARTIAL",
    "BEST_EFFORT",
    "UNMAPPED",
)

_CACHE: Dict[str, dict] = {}
_SASA_CACHE: Dict[str, dict] = {}
_DSSP_CACHE: Dict[str, dict] = {}
MAX_CACHE: int = 8
SASA_ALGO_VERSION: str = "shrake_rupley_1973_neighbor_list_v1"


class StructureError(Exception):
    """Falha classificada de estrutura. Nunca deve virar coordenadas falsas.

    Attributes:
        category: TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE, NETWORK_ERROR,
            INVALID_INPUT, PARSING_ERROR, RESOURCE_LIMIT, NOT_FOUND, NO_STRUCTURE,
            TOOL_FAILED.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "TOOL_FAILED")


def residue_one_letter(comp_id: str) -> str:
    """Traduz residuo PDB de 3 letras para 1 letra. Nucleotidos reais, nao X.

    Args:
        comp_id: Nome do composto (ALA, DA, G, ...).

    Returns:
        Uma letra, ou X se desconhecido.

    Raises:
        Nenhum.
    """
    key = str(comp_id or "").strip().upper()
    if key in AA3_TO_1:
        return AA3_TO_1[key]
    if key in NA3_TO_1:
        return NA3_TO_1[key]
    return "X"


def dssp_availability() -> dict:
    """Detecta DSSP/mkdssp allowlisted. Ausente = UNAVAILABLE, nao heuristica.

    Args:
        Nenhum.

    Returns:
        Dict available, version, reason.

    Raises:
        Nenhum.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        DSSP_NAMES,
        extra_file_candidates=tool_paths.candidates_for(DSSP_NAMES),
        env_var=DSSP_ENV,
    )
    if executable is None:
        return {
            "available": False,
            "tool": "DSSP",
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "reason": (
                "DSSP/mkdssp is not installed. HelixScope does not invent "
                "secondary-structure assignments. Deposited HELIX/SHEET records "
                "are a separate annotation when present in the file."
            ),
        }
    version = ""
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
        blob = f"{completed.stdout or ''}\n{completed.stderr or ''}"
        match = re.search(r"(\d+\.\d+(?:\.\d+)?)", blob)
        if match:
            version = match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        version = ""
    return {
        "available": True,
        "tool": "DSSP",
        "version": version,
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
        "path": tool_detection.sanitize_tool_path(executable),
        "reason": "DSSP executable detected. Not run automatically.",
    }


def dssp_remote_availability() -> dict:
    """Declares the PDB-REDO DSSP HTTP API. Does not ping the network.

    Args:
        Nenhum.

    Returns:
        Dict describing the remote mkdssp backend. Never LIVE_VALIDATED local.

    Raises:
        Nenhum.

    Nota biologica:
        DSSP assigns secondary structure from atomic coordinates. This remote
        path is the official PDB-REDO mkdssp service, not a HelixScope
        reimplementation and not sequence prediction.
    """
    live = engine_validation.live_record(DSSP_REMOTE_TOOL)
    status = str((live or {}).get("status") or "detected")
    return {
        "available": True,
        "tool": DSSP_REMOTE_TOOL,
        "version": str((live or {}).get("version") or "mkdssp via PDB-REDO"),
        "status": status,
        "engine_location": "remote",
        "endpoint": PDB_REDO_DSSP_DO,
        "reason": (
            "Official PDB-REDO DSSP HTTP API (https://pdb-redo.eu/dssp/do). "
            "Not a local mkdssp binary. Success is REMOTE_VALIDATED, never "
            "LIVE_VALIDATED local. DSSP is geometric assignment, not prediction."
        ),
        "live": dict(live) if isinstance(live, dict) else None,
    }


def stride_availability() -> dict:
    """Detecta o binario STRIDE (Frishman/Argos). Ausente = UNAVAILABLE.

    Args:
        Nenhum.

    Returns:
        Dict available, tool, version, reason. Nunca devolve atribuicao
        geometrica inventada.

    Raises:
        Nenhum.

    Nota biologica:
        STRIDE nao e predicao a partir da sequencia. HelixScope nao implementa
        uma imitacao. BACKBONE visual stride (LOD) e outro conceito.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        STRIDE_NAMES,
        extra_file_candidates=tool_paths.candidates_for(STRIDE_NAMES),
        env_var=STRIDE_ENV,
    )
    if executable is None:
        return {
            "available": False,
            "tool": "STRIDE",
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "reason": (
                "STRIDE (Frishman and Argos 1995) is not installed. HelixScope "
                "does not ship a house STRIDE. Deposited mmCIF _struct_conf / "
                "_struct_sheet_range remain a separate depositor annotation."
            ),
        }
    version = ""
    try:
        completed = subprocess.run(
            [executable],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
        blob = f"{completed.stdout or ''}\n{completed.stderr or ''}"
        match = re.search(r"(\d+\.\d+(?:\.\d+)?)", blob)
        if match:
            version = match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        version = ""
    return {
        "available": True,
        "tool": "STRIDE",
        "version": version,
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
        "path": tool_detection.sanitize_tool_path(executable),
        "reason": "STRIDE executable detected. Not run automatically.",
    }


def surface_tool_research() -> list[dict]:
    """Pesquisa documentada de motores de superficie molecular.

    Args:
        Nenhum.

    Returns:
        Uma linha por ferramenta. Nenhuma e integrada nesta fase.

    Raises:
        Nenhum.
    """
    return [
        {
            "tool": "MSMS",
            "algorithm": "Solvent-excluded surface (Sanner)",
            "version": "bioconda scripps-msms 2.6.1 listed; not installed here",
            "license": "Free for academic use; commercial requires contacting Sanner",
            "performance": "widely used; not measured in this environment (binary absent)",
            "output": "vertex/face SES mesh",
            "windows_compatibility": "Official Windows zip on ccsb.scripps.edu/msms/downloads/",
            "decision": "integrated when tools/msms/msms.exe is present",
            "why": (
                "MSMS is free for academic/teaching use; commercial research requires "
                "contacting Dr Sanner. HelixScope does not commit the binary. The "
                "official Windows zip from ccsb.scripps.edu may be installed locally "
                "under gitignored tools/msms/."
            ),
        },
        {
            "tool": "EDTSurf",
            "algorithm": "Euclidean distance transform macromolecular surfaces (Xu/Zhang)",
            "version": "official Windows executable from zhanggroup.org/EDTSurf/",
            "license": "permissive Zhang-lab disclaimer (retain copyright and citation)",
            "performance": "not measured (binary absent)",
            "output": "PLY triangulated VWS/SAS/MS",
            "windows_compatibility": "Official Windows executable advertised by Zhang lab",
            "decision": "integrated when tools/edtsurf/edtsurf.exe is present",
            "why": (
                "Official HTTPS download from zhanggroup.org/EDTSurf/EDTSurf.exe. "
                "Authors did not publish a checksum; HelixScope records a local SHA-256 "
                "after TLS retrieval and PE header check. Distinct from Shrake-Rupley SASA."
            ),
        },
        {
            "tool": "Mol*",
            "algorithm": "Gaussian/molecular-surface representation in the Mol* viewer",
            "version": "Mol* is a TypeScript web library, not a HelixScope Python engine",
            "license": "Molstar (typically Apache-2.0 / MIT on GitHub molstar/molstar)",
            "performance": "GPU viewer; HelixScope renderer is Plotly data-only without Mol*",
            "output": "GPU mesh inside the Mol* plugin, not a portable vertex array here",
            "windows_compatibility": "browser/JS, not a local scientific CLI",
            "decision": "not integrated",
            "why": (
                "HelixScope 3D is an isolated Plotly scene. Bundling Mol* would be "
                "a renderer rewrite (Phase 10B+), not a scientific surface engine."
            ),
        },
    ]


def edtsurf_availability() -> dict:
    """Detecta so o EDTSurf oficial. Nao trata MSMS como EDTSurf.

    Args:
        Nenhum.

    Returns:
        Dict available, tool, path, reason.

    Raises:
        Nenhum.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        EDTSURF_NAMES,
        extra_file_candidates=tool_paths.candidates_for(EDTSURF_NAMES),
        env_var=SURFACE_ENV,
    )
    if executable is None:
        return {
            "available": False,
            "tool": "EDTSurf",
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "reason": "EDTSurf is not installed.",
        }
    return {
        "available": True,
        "tool": "edtsurf",
        "version": "",
        "status": tool_detection.TOOL_STATUS_INSTALLED,
        "path": tool_detection.sanitize_tool_path(executable),
        "reason": "EDTSurf executable detected. Not run automatically.",
    }


def msms_availability() -> dict:
    """Detecta o MSMS oficial (Sanner). Ausente = NOT_INSTALLED, nao malha inventada.

    Args:
        Nenhum.

    Returns:
        Dict available, tool, path, reason. Licenca academica/ensino.

    Raises:
        Nenhum.

    Nota biologica:
        MSMS calcula superficie excluida do solvente (SES). Nao e Shrake-Rupley
        SASA nem EDTSurf. Nao e DSSP ACC.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        MSMS_NAMES,
        extra_file_candidates=tool_paths.candidates_for(MSMS_NAMES),
        env_var=MSMS_ENV,
    )
    if executable is None:
        return {
            "available": False,
            "tool": "MSMS",
            "version": "",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "reason": (
                "MSMS is not installed. Academic/teaching use of the official "
                "Scripps CCSB Windows zip is allowed; commercial research needs "
                "a separate agreement. HelixScope does not invent a SES mesh."
            ),
        }
    return {
        "available": True,
        "tool": "MSMS",
        "version": "",
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=""
        ),
        "path": tool_detection.sanitize_tool_path(executable),
        "reason": (
            "MSMS executable detected. Academic/teaching license. "
            "Version is recorded after a live SES run, not guessed from the "
            "filename. Not run automatically. Not EDTSurf."
        ),
    }


def surface_availability() -> dict:
    """Detecta um motor de superficie. Prefere EDTSurf se ambos existirem.

    Args:
        Nenhum.

    Returns:
        Dict available=False se nenhum binario allowlisted.

    Raises:
        Nenhum.
    """
    research = surface_tool_research()
    edt = edtsurf_availability()
    if edt.get("available"):
        return {**edt, "research": research}
    msms = msms_availability()
    if msms.get("available"):
        return {**msms, "research": research}
    return {
        "available": False,
        "tool": None,
        "version": "",
        "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
        "research": research,
        "reason": (
            "No molecular-surface engine (MSMS, EDTSurf) is installed. "
            "UI glass/transparency is not a molecular surface."
        ),
    }


def atoms_to_pdb_text(atoms: Sequence[Mapping[str, Any]]) -> str:
    """Write a minimal PDB from parsed ATOM coordinates. No invented XYZ.

    Args:
        atoms: Parsed atom dicts with finite x, y, z.

    Returns:
        PDB text with ATOM/HETATM rows.

    Raises:
        StructureError: INVALID_INPUT se nao houver atomos com coordenadas.
    """
    lines: List[str] = ["HEADER    HELIXSCOPE SURFACE INPUT"]
    serial = 0
    for atom in atoms:
        if not scientific_checks.atom_coordinate_is_finite(atom.get("x"), atom.get("y"), atom.get("z")):
            continue
        serial += 1
        if serial > MAX_ATOMS:
            raise StructureError(
                f"Surface input exceeds {MAX_ATOMS:,} atoms.",
                "RESOURCE_LIMIT",
            )
        group = str(atom.get("group") or "ATOM").ljust(6)[:6]
        name = str(atom.get("atom_name") or " C").strip()[:4]
        name_field = name.rjust(4) if len(name) < 4 else name[:4]
        res = str(atom.get("comp_id") or "UNK")[:3].rjust(3)
        chain = str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "A")[:1] or "A"
        seq = int(atom.get("auth_seq_id") or atom.get("label_seq_id") or serial)
        seq = max(1, min(seq, 9999))
        occ = atom.get("occupancy")
        occ_v = 1.0 if occ is None else float(occ)
        biso = atom.get("b_iso")
        b_v = 0.0 if biso is None else float(biso)
        elem = str(atom.get("element") or name[:1] or "C")[:2].rjust(2)
        line = (
            f"{group}{serial:5d} {name_field} {res} {chain}{seq:4d}    "
            f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}"
            f"{occ_v:6.2f}{b_v:6.2f}          {elem}"
        )
        lines.append(line)
    if serial == 0:
        raise StructureError("No finite atomic coordinates for molecular surface.", "INVALID_INPUT")
    lines.append("END")
    return "\n".join(lines) + "\n"


def parse_ply_ascii(text: str) -> dict:
    """Parse an ASCII PLY mesh. Does not invent vertices.

    Args:
        text: PLY file body.

    Returns:
        Dict vertices, faces, n_vertices, n_faces.

    Raises:
        StructureError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    blob = str(text or "")
    if not blob.lstrip().startswith("ply"):
        raise StructureError("Surface output is not an ASCII PLY mesh.", "PARSING_ERROR")
    lines = blob.splitlines()
    n_vertices = 0
    n_faces = 0
    header_end = None
    fmt = ""
    for index, line in enumerate(lines):
        raw = line.strip()
        if raw.startswith("format"):
            fmt = raw.lower()
        elif raw.startswith("element vertex"):
            parts = raw.split()
            n_vertices = int(parts[-1])
        elif raw.startswith("element face"):
            parts = raw.split()
            n_faces = int(parts[-1])
        elif raw == "end_header":
            header_end = index
            break
    if header_end is None:
        raise StructureError("PLY header has no end_header.", "PARSING_ERROR")
    if "ascii" not in fmt:
        raise StructureError("Only ASCII PLY surface meshes are parsed.", "PARSING_ERROR")
    if n_vertices <= 0 or n_vertices > MAX_SURFACE_VERTICES:
        raise StructureError(
            f"PLY vertex count {n_vertices} is outside 1..{MAX_SURFACE_VERTICES}.",
            "RESOURCE_LIMIT" if n_vertices > MAX_SURFACE_VERTICES else "PARSING_ERROR",
        )
    if n_faces < 0 or n_faces > MAX_SURFACE_FACES:
        raise StructureError(
            f"PLY face count {n_faces} exceeds {MAX_SURFACE_FACES}.",
            "RESOURCE_LIMIT",
        )
    body = lines[header_end + 1 :]
    vertices: List[dict] = []
    for row in body[:n_vertices]:
        parts = row.split()
        if len(parts) < 3:
            raise StructureError("PLY vertex row is truncated.", "PARSING_ERROR")
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        if not scientific_checks.atom_coordinate_is_finite(x, y, z):
            raise StructureError("PLY vertex has non-finite coordinates.", "PARSING_ERROR")
        item = {"x": x, "y": y, "z": z}
        if len(parts) >= 6:
            nx, ny, nz = float(parts[3]), float(parts[4]), float(parts[5])
            if scientific_checks.atom_coordinate_is_finite(nx, ny, nz):
                item["nx"], item["ny"], item["nz"] = nx, ny, nz
        vertices.append(item)
    if len(vertices) != n_vertices:
        raise StructureError("PLY vertex block is shorter than the header count.", "PARSING_ERROR")
    faces: List[list[int]] = []
    for row in body[n_vertices : n_vertices + n_faces]:
        parts = row.split()
        if not parts:
            continue
        count = int(float(parts[0]))
        idx = [int(item) for item in parts[1 : 1 + count]]
        if len(idx) < 3:
            raise StructureError("PLY face has fewer than 3 indices.", "PARSING_ERROR")
        if any(item < 0 or item >= n_vertices for item in idx):
            raise StructureError("PLY face index is out of range.", "PARSING_ERROR")
        faces.append(idx)
    if len(faces) != n_faces:
        raise StructureError("PLY face block is shorter than the header count.", "PARSING_ERROR")
    validate_surface_mesh(vertices)
    return {
        "vertices": vertices,
        "faces": faces,
        "n_vertices": n_vertices,
        "n_faces": n_faces,
        "format": "ply_ascii",
    }


def atoms_to_xyzr(atoms: Sequence[Mapping[str, Any]]) -> str:
    """Write MSMS xyzr input from deposited coordinates and Bondi radii.

    Args:
        atoms: Parsed atom dicts with finite x, y, z.

    Returns:
        Text with one `x y z radius` row per atom.

    Raises:
        StructureError: INVALID_INPUT se nao houver atomos com coordenadas.

    Nota biologica:
        Raios sao Bondi 1964 (os mesmos de Shrake-Rupley neste modulo), nao a
        tabela atmtypenumbers do pdb_to_xyzr da distribuicao MSMS. O motor SES
        continua a ser o MSMS oficial.
    """
    lines: List[str] = []
    serial = 0
    for atom in atoms:
        if not scientific_checks.atom_coordinate_is_finite(atom.get("x"), atom.get("y"), atom.get("z")):
            continue
        serial += 1
        if serial > MAX_ATOMS:
            raise StructureError(
                f"MSMS input exceeds {MAX_ATOMS:,} atoms.",
                "RESOURCE_LIMIT",
            )
        elem = str(atom.get("element") or atom.get("atom_name") or "C").strip().upper()
        if elem[:2] in BONDI_RADII_ANGSTROM:
            radius = float(BONDI_RADII_ANGSTROM[elem[:2]])
        elif elem[:1] in BONDI_RADII_ANGSTROM:
            radius = float(BONDI_RADII_ANGSTROM[elem[:1]])
        else:
            radius = float(BONDI_RADII_ANGSTROM["C"])
        lines.append(
            f"{float(atom['x']):9.5f} {float(atom['y']):9.5f} "
            f"{float(atom['z']):9.5f} {radius:6.3f}"
        )
    if serial == 0:
        raise StructureError("No finite atomic coordinates for MSMS xyzr input.", "INVALID_INPUT")
    return "\n".join(lines) + "\n"


def parse_msms_vert(text: str) -> list[dict]:
    """Parse an MSMS .vert file. Does not invent vertices.

    Args:
        text: Body of the .vert file.

    Returns:
        Lista de dicts x,y,z e normais se presentes.

    Raises:
        StructureError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    vertices: List[dict] = []
    for raw in str(text or "").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        offset = 0
        if len(parts) >= 7:
            try:
                int(parts[0])
                float(parts[1])
                offset = 1
            except ValueError:
                offset = 0
        try:
            x = float(parts[offset])
            y = float(parts[offset + 1])
            z = float(parts[offset + 2])
        except (TypeError, ValueError, IndexError) as exc:
            raise StructureError("MSMS vertex row is not numeric.", "PARSING_ERROR") from exc
        if not scientific_checks.atom_coordinate_is_finite(x, y, z):
            raise StructureError("MSMS vertex has non-finite coordinates.", "PARSING_ERROR")
        item: dict[str, Any] = {"x": x, "y": y, "z": z}
        if len(parts) >= offset + 6:
            try:
                nx = float(parts[offset + 3])
                ny = float(parts[offset + 4])
                nz = float(parts[offset + 5])
            except (TypeError, ValueError):
                nx = ny = nz = None
            if nx is not None and scientific_checks.atom_coordinate_is_finite(nx, ny, nz):
                item["nx"], item["ny"], item["nz"] = nx, ny, nz
        vertices.append(item)
        if len(vertices) > MAX_SURFACE_VERTICES:
            raise StructureError(
                f"MSMS vertex count exceeds {MAX_SURFACE_VERTICES}.",
                "RESOURCE_LIMIT",
            )
    if not vertices:
        raise StructureError("MSMS .vert file has no vertex rows.", "PARSING_ERROR")
    validate_surface_mesh(vertices)
    return vertices


def parse_msms_face(text: str, n_vertices: int) -> list[list[int]]:
    """Parse an MSMS .face file. Indices 1-based in the official format.

    Args:
        text: Body of the .face file.
        n_vertices: Numero de vertices ja parseados.

    Returns:
        Lista de triangulos com indices 0-based.

    Raises:
        StructureError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    faces: List[list[int]] = []
    n_vert = int(n_vertices)
    for raw in str(text or "").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                float(parts[2])
                float(parts[3])
            except (TypeError, ValueError):
                pass
            else:
                if "." in parts[2] or "." in parts[3]:
                    continue
        nums: list[int] = []
        for token in parts[:6]:
            try:
                nums.append(int(float(token)))
            except (TypeError, ValueError):
                break
        if len(nums) < 3:
            continue
        if len(nums) >= 4 and nums[0] == 3:
            idx = nums[1:4]
        else:
            idx = nums[:3]
        if max(idx) <= n_vert and min(idx) >= 1:
            idx = [item - 1 for item in idx]
        if len(idx) < 3 or any(item < 0 or item >= n_vert for item in idx):
            raise StructureError("MSMS face index is out of range.", "PARSING_ERROR")
        faces.append(idx[:3])
        if len(faces) > MAX_SURFACE_FACES:
            raise StructureError(
                f"MSMS face count exceeds {MAX_SURFACE_FACES}.",
                "RESOURCE_LIMIT",
            )
    if not faces:
        raise StructureError("MSMS .face file has no triangle rows.", "PARSING_ERROR")
    return faces


def _run_edtsurf_surface(
    atoms: Sequence[Mapping[str, Any]],
    *,
    kind: str,
    structure_hash: str,
    surface_code: int,
    probe_radius: float,
    run_fn=None,
) -> dict:
    code = int(surface_code)
    if code not in EDTSURF_SURFACE_CODES:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={"reason": f"Unknown EDTSurf surface code {surface_code}.", "n_vertices": 0},
            status="ERROR",
            algorithm="EDTSurf",
            parameters={"surface_code": surface_code},
            source="HelixScope EDTSurf",
        ) | {"status": "ERROR", "available": False, "reason": "Unknown EDTSurf surface code."}
    try:
        probe = float(probe_radius)
    except (TypeError, ValueError):
        probe = -1.0
    if not math.isfinite(probe) or probe < 0.0 or probe > 2.0:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={"reason": "EDTSurf probe radius must be in [0, 2.0].", "n_vertices": 0},
            status="ERROR",
            algorithm="EDTSurf",
            parameters={"probe_radius": probe_radius},
            source="HelixScope EDTSurf",
        ) | {"status": "ERROR", "available": False, "reason": "Invalid EDTSurf probe radius."}
    definition = EDTSURF_SURFACE_CODES[code]
    executable = tool_detection.resolve_allowlisted_executable(
        EDTSURF_NAMES,
        extra_file_candidates=tool_paths.candidates_for(EDTSURF_NAMES),
        env_var=SURFACE_ENV,
    )
    if executable is None:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={"reason": "EDTSurf executable is not allowlisted.", "n_vertices": 0},
            status="UNAVAILABLE",
            algorithm="none",
            parameters={},
            source="HelixScope EDTSurf",
        ) | {"status": "UNAVAILABLE", "available": False, "reason": "EDTSurf executable is not allowlisted."}
    try:
        pdb_text = atoms_to_pdb_text(atoms)
    except StructureError as exc:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={"reason": str(exc), "n_vertices": 0},
            status="ERROR",
            algorithm="EDTSurf",
            parameters={},
            source="HelixScope EDTSurf",
        ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": exc.category}
    tmpdir = tempfile.mkdtemp(prefix="helixscope_edtsurf_")
    try:
        infile = os.path.join(tmpdir, "input.pdb")
        with open(infile, "w", encoding="ascii", errors="replace") as handle:
            handle.write(pdb_text)
        prefix = os.path.join(tmpdir, "out")
        runner = run_fn or subprocess.run
        try:
            completed = runner(
                [
                    executable,
                    "-i",
                    infile,
                    "-o",
                    prefix,
                    "-s",
                    str(code),
                    "-p",
                    f"{probe:.2f}",
                    "-t",
                    "2",
                    "-h",
                    "2",
                ],
                capture_output=True,
                text=True,
                timeout=SURFACE_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmpdir,
            )
        except subprocess.TimeoutExpired as exc:
            return provenance.analysis_envelope(
                module="PROTEIN surface",
                payload={"reason": f"EDTSurf exceeded {int(SURFACE_TIMEOUT_S)}s.", "n_vertices": 0},
                status="ERROR",
                algorithm="EDTSurf",
                parameters={"timeout_s": SURFACE_TIMEOUT_S},
                source="HelixScope EDTSurf",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": "TIMEOUT"}
        except OSError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN surface",
                payload={"reason": str(exc), "n_vertices": 0},
                status="ERROR",
                algorithm="EDTSurf",
                parameters={},
                source="HelixScope EDTSurf",
            ) | {"status": "ERROR", "available": False, "reason": str(exc)}
        ply_path = prefix + ".ply"
        if not os.path.isfile(ply_path):
            stderr = str(getattr(completed, "stderr", "") or "")[:400]
            return provenance.analysis_envelope(
                module="PROTEIN surface",
                payload={
                    "reason": "EDTSurf did not write a PLY mesh. No surface is invented.",
                    "stderr": stderr,
                    "n_vertices": 0,
                },
                status="ERROR",
                algorithm="EDTSurf",
                parameters={"returncode": getattr(completed, "returncode", None)},
                source="HelixScope EDTSurf",
            ) | {
                "status": "ERROR",
                "available": False,
                "reason": "EDTSurf produced no PLY mesh.",
                "category": "TOOL_FAILED",
            }
        ply_text = Path(ply_path).read_text(encoding="utf-8", errors="replace")
        try:
            mesh = parse_ply_ascii(ply_text)
        except StructureError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN surface",
                payload={"reason": str(exc), "n_vertices": 0},
                status="ERROR",
                algorithm="EDTSurf PLY parser",
                parameters={},
                source="HelixScope EDTSurf",
            ) | {
                "status": "ERROR",
                "available": False,
                "reason": str(exc),
                "category": exc.category,
            }
        status = "EXPERIMENTAL" if kind == "experimental" else "COMPUTED"
        if run_fn is None:
            engine_validation.record_live_validation(
                "EDTSurf",
                ok=True,
                version="zhanggroup.org Windows executable",
                details={
                    "n_vertices": mesh["n_vertices"],
                    "n_faces": mesh["n_faces"],
                    "surface_definition": definition,
                    "probe_radius": probe,
                },
            )
        envelope = provenance.analysis_envelope(
            module="PROTEIN surface",
            payload=mesh,
            status=status,
            algorithm="EDTSurf Euclidean distance transform surface",
            parameters={
                "backend": "EDTSurf",
                "surface_definition": definition,
                "surface_code": code,
                "probe_radius_angstrom": probe,
                "triangulation": "VCMC",
                "outer_only": True,
                "structure_hash": structure_hash,
                "source_structure": structure_hash,
            },
            source="HelixScope EDTSurf",
        ) | {
            "status": status,
            "available": True,
            "tool": "EDTSurf",
            "backend": "EDTSurf",
            "surface_definition": definition,
            "probe_radius_angstrom": probe,
            "structure_hash": structure_hash,
            "not_shrake_rupley": True,
            "reason": (
                f"EDTSurf {definition} mesh from deposited/predicted coordinates. "
                "This is not Shrake-Rupley SASA."
            ),
        }
        return envelope
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _run_msms_surface(
    atoms: Sequence[Mapping[str, Any]],
    *,
    kind: str,
    structure_hash: str,
    probe_radius: float,
    run_fn=None,
) -> dict:
    """Executa MSMS sobre xyzr derivado de atomos depositados.

    Args:
        atoms: Atomos parseados com coordenadas Cartesianas.
        kind: experimental ou predicted.
        structure_hash: Hash da estrutura de origem.
        probe_radius: Raio da sonda em Angstrom (-probe_radius).
        run_fn: subprocess.run injetavel.

    Returns:
        Envelope com malha validada, UNAVAILABLE ou ERROR.

    Raises:
        Nenhum.

    Nota biologica:
        MSMS (Sanner, Olson, Spehner 1996) calcula a superficie analitica
        solvent-excluded. HelixScope nao reimplementa o algoritmo; so converte
        atomos para xyzr com raios Bondi, executa o binario oficial e valida
        .vert/.face. Os raios Bondi nao sao as tabelas pdb_to_xyzr originais
        do MSMS.
    """
    avail = msms_availability()
    if not avail.get("available"):
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={
                "kind": kind,
                "n_vertices": 0,
                "structure_hash": structure_hash,
                "reason": avail.get("reason"),
                "tool": "MSMS",
            },
            status="UNAVAILABLE",
            algorithm="none",
            parameters={"kind": kind},
            source="HelixScope surface",
        )
    executable = tool_detection.resolve_allowlisted_executable(
        MSMS_NAMES,
        extra_file_candidates=tool_paths.candidates_for(MSMS_NAMES),
        env_var=MSMS_ENV,
    )
    if executable is None:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={
                "kind": kind,
                "n_vertices": 0,
                "structure_hash": structure_hash,
                "reason": "MSMS executable is not allowlisted.",
                "tool": "MSMS",
            },
            status="UNAVAILABLE",
            algorithm="none",
            parameters={"kind": kind},
            source="HelixScope surface",
        )
    try:
        xyzr_text = atoms_to_xyzr(atoms)
    except StructureError as exc:
        return {
            "status": "ERROR",
            "available": False,
            "reason": str(exc),
            "tool": "MSMS",
        }
    tmpdir = tempfile.mkdtemp(prefix="helixscope_msms_")
    try:
        xyzr_path = os.path.join(tmpdir, "protein.xyzr")
        with open(xyzr_path, "w", encoding="ascii", newline="\n") as handle:
            handle.write(xyzr_text)
        runner = run_fn or subprocess.run
        cmd = [
            executable,
            "-if",
            "protein.xyzr",
            "-of",
            "out",
            "-probe_radius",
            f"{float(probe_radius):.4f}",
        ]
        try:
            completed = runner(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
                shell=False,
                cwd=tmpdir,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return {
                "status": "ERROR",
                "available": False,
                "reason": f"MSMS execution failed: {exc}",
                "tool": "MSMS",
            }
        rc = int(getattr(completed, "returncode", 1) or 0)
        if rc != 0:
            stderr = str(getattr(completed, "stderr", "") or "")[:400]
            return {
                "status": "ERROR",
                "available": False,
                "reason": f"MSMS exited {rc}. {stderr}".strip(),
                "tool": "MSMS",
            }
        vert_path = os.path.join(tmpdir, "out.vert")
        face_path = os.path.join(tmpdir, "out.face")
        if not os.path.isfile(vert_path) or not os.path.isfile(face_path):
            return {
                "status": "ERROR",
                "available": False,
                "reason": (
                    "MSMS exited 0 but did not write .vert/.face. "
                    "No mesh is invented."
                ),
                "tool": "MSMS",
            }
        with open(vert_path, encoding="utf-8", errors="replace") as handle:
            vertices = parse_msms_vert(handle.read())
        with open(face_path, encoding="utf-8", errors="replace") as handle:
            faces = parse_msms_face(handle.read(), n_vertices=len(vertices))
        try:
            validate_surface_mesh(vertices)
        except StructureError as exc:
            return {
                "status": "ERROR",
                "available": False,
                "reason": str(exc),
                "tool": "MSMS",
            }
        mesh_status = "EXPERIMENTAL" if kind == "experimental" else "COMPUTED"
        if run_fn is None:
            blob = (
                f"{getattr(completed, 'stdout', '') or ''}\n"
                f"{getattr(completed, 'stderr', '') or ''}"
            )
            found = re.search(r"MSMS\s+([0-9]+(?:\.[0-9]+)*)", blob, flags=re.I)
            live_version = found.group(1) if found else str(avail.get("version") or "")
            engine_validation.record_live_validation(
                "MSMS",
                ok=True,
                version=live_version,
                details={
                    "n_vertices": len(vertices),
                    "n_faces": len(faces),
                    "probe_radius": float(probe_radius),
                    "radii": "Bondi (not MSMS pdb_to_xyzr tables)",
                },
            )
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={
                "kind": kind,
                "vertices": vertices,
                "faces": faces,
                "n_vertices": len(vertices),
                "n_faces": len(faces),
                "structure_hash": structure_hash,
                "backend": "MSMS",
                "surface_definition": "SES",
                "not_shrake_rupley": True,
                "radii_source": "Bondi (HelixScope); not MSMS pdb_to_xyzr",
                "probe_radius": float(probe_radius),
                "reason": None,
                "tool": "MSMS",
            },
            status=mesh_status,
            algorithm="MSMS SES (Sanner 1996 binary)",
            parameters={
                "kind": kind,
                "probe_radius": float(probe_radius),
                "n_atoms": len(atoms),
            },
            source="HelixScope surface",
        ) | {
            "status": mesh_status,
            "available": True,
            "backend": "MSMS",
            "surface_definition": "SES",
            "not_shrake_rupley": True,
            "n_vertices": len(vertices),
            "n_faces": len(faces),
            "vertices": vertices,
            "faces": faces,
            "tool": "MSMS",
        }
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def generate_molecular_surface(
    *,
    kind: str,
    vertices: Optional[Sequence[Mapping[str, Any]]] = None,
    structure_hash: str = "",
    parsed: Optional[Mapping[str, Any]] = None,
    surface_code: int = 3,
    probe_radius: float = 1.4,
    run_fn=None,
) -> dict:
    """Gera superficie molecular so a partir de coordenadas EXPERIMENTAL/PREDICTED.

    Args:
        kind: experimental, predicted ou illustrative.
        vertices: Ignorado; a malha vem do motor, nunca de vertices inventados.
        structure_hash: Hash da estrutura de origem.
        parsed: Estrutura parseada com atoms reais.
        surface_code: EDTSurf -s: 1 VDW, 2 SAS, 3 MS, 4 SES.
        probe_radius: Raio de sonda em Angstrom (EDTSurf -p).
        run_fn: subprocess.run injetavel.

    Returns:
        Envelope com malha validada, UNAVAILABLE ou ERROR.

    Raises:
        StructureError: INVALID_INPUT se kind for illustrative.

    Nota biologica:
        Superficie EDTSurf (VDW/SAS/MS/SES) ou MSMS SES nao e o escalar
        Shrake-Rupley SASA. ILLUSTRATIVE nao recebe superficie molecular.
        Quando EDTSurf e MSMS estao ambos instalados, EDTSurf e o default.
    """
    key = str(kind or "").strip().lower()
    if key in {"illustrative", "unavailable", ""}:
        raise StructureError(
            "Molecular surface is not generated from ILLUSTRATIVE or missing coordinates.",
            "INVALID_INPUT",
        )
    if key not in {"experimental", "predicted"}:
        raise StructureError(
            "Molecular surface requires experimental or predicted deposited coordinates.",
            "INVALID_INPUT",
        )
    avail = surface_availability()
    if not avail.get("available"):
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={
                "kind": key,
                "n_vertices": 0,
                "structure_hash": structure_hash,
                "reason": avail.get("reason"),
                "tool": None,
            },
            status="UNAVAILABLE",
            algorithm="none",
            parameters={"kind": key},
            source="HelixScope surface",
        )
    _ = vertices
    atoms = list((parsed or {}).get("atoms") or [])
    if not atoms:
        return provenance.analysis_envelope(
            module="PROTEIN surface",
            payload={
                "kind": key,
                "n_vertices": 0,
                "structure_hash": structure_hash,
                "reason": (
                    "A surface engine is present but no atomic coordinates were "
                    "provided. No mesh is invented."
                ),
                "tool": avail.get("tool"),
            },
            status="UNAVAILABLE",
            algorithm="none",
            parameters={"kind": key},
            source="HelixScope surface",
        ) | {
            "status": "UNAVAILABLE",
            "available": False,
            "reason": "No atomic coordinates were provided for molecular surface.",
            "tool": avail.get("tool"),
        }
    tool_name = str(avail.get("tool") or "").lower()
    if "edtsurf" in tool_name:
        return _run_edtsurf_surface(
            atoms,
            kind=key,
            structure_hash=structure_hash,
            surface_code=surface_code,
            probe_radius=probe_radius,
            run_fn=run_fn,
        )
    if "msms" in tool_name:
        return _run_msms_surface(
            atoms,
            kind=key,
            structure_hash=structure_hash,
            probe_radius=probe_radius,
            run_fn=run_fn,
        )
    return provenance.analysis_envelope(
        module="PROTEIN surface",
        payload={
            "kind": key,
            "n_vertices": 0,
            "structure_hash": structure_hash,
            "reason": (
                f"Surface tool {avail.get('tool')} is detected but has no HelixScope "
                "execution adapter. No mesh is invented."
            ),
            "tool": avail.get("tool"),
        },
        status="UNAVAILABLE",
        algorithm="none",
        parameters={"kind": key},
        source="HelixScope surface",
    ) | {
        "status": "UNAVAILABLE",
        "available": False,
        "reason": "Detected surface tool has no execution adapter.",
        "tool": avail.get("tool"),
    }


def validate_surface_mesh(vertices: Sequence[Mapping[str, Any]]) -> dict:
    """Valida vertices de uma malha: finitude, coordenadas, topologia minima.

    Args:
        vertices: Lista de dicts x,y,z e opcionalmente nx,ny,nz.

    Returns:
        Dict ok, n_vertices, issues.

    Raises:
        StructureError: PARSING_ERROR se NaN/Inf ou lista vazia.
    """
    if not vertices:
        raise StructureError("Surface mesh has no vertices.", "PARSING_ERROR")
    issues: List[str] = []
    for index, item in enumerate(vertices):
        if not isinstance(item, Mapping):
            raise StructureError(f"Surface vertex {index} is not a mapping.", "PARSING_ERROR")
        x, y, z = item.get("x"), item.get("y"), item.get("z")
        if not scientific_checks.atom_coordinate_is_finite(x, y, z):
            raise StructureError(
                f"Surface vertex {index} has non-finite coordinates.",
                "PARSING_ERROR",
            )
        nx, ny, nz = item.get("nx"), item.get("ny"), item.get("nz")
        if any(value is not None for value in (nx, ny, nz)):
            if not scientific_checks.atom_coordinate_is_finite(nx, ny, nz):
                raise StructureError(
                    f"Surface vertex {index} has non-finite normals.",
                    "PARSING_ERROR",
                )
    if len(vertices) < 3:
        issues.append("Fewer than 3 vertices: topology cannot form a triangle.")
    return {
        "ok": not issues,
        "n_vertices": len(vertices),
        "issues": issues,
        "finite": True,
    }


def parse_dssp_legacy_output(text: str) -> dict:
    """Parseia o formato DSSP de colunas fixas. Nao inventa H/E/C.

    Args:
        text: Stdout do mkdssp --output-format=dssp.

    Returns:
        Dict residues (index, aa, code), n_helix, n_sheet, n_coil.

    Raises:
        StructureError: PARSING_ERROR se o texto nao tiver residuos DSSP.

    Nota biologica:
        DSSP e assignment geometrico de uma estrutura ja resolvida, nao
        predicao a partir da sequencia. Codigos: H/G/I helice, E/B folha,
        T/S/espaco coil no esquema Kabsch-Sander.
    """
    blob = str(text or "")
    lines = blob.splitlines()
    start = None
    for index, line in enumerate(lines):
        if "RESIDUE" in line and "STRUCTURE" in line:
            start = index + 1
            break
    if start is None:
        raise StructureError(
            "DSSP output has no RESIDUE/STRUCTURE header. Return code is not proof of assignment.",
            "PARSING_ERROR",
        )
    residues: List[dict] = []
    for line in lines[start:]:
        if len(line) < 17:
            continue
        code = line[16]
        aa = line[13] if len(line) > 13 else " "
        num_blob = line[5:11].strip()
        if not num_blob:
            continue
        try:
            seqnum = int(num_blob)
        except ValueError:
            continue
        residues.append(
            {
                "auth_seq_id": seqnum,
                "chain_id": (line[11].strip() if len(line) > 11 else ""),
                "aa": aa,
                "code": code if code.strip() else " ",
            }
        )
    if not residues:
        raise StructureError(
            "DSSP output header was present but no residue rows parsed.",
            "PARSING_ERROR",
        )
    helix = sum(1 for item in residues if item["code"] in {"H", "G", "I"})
    sheet = sum(1 for item in residues if item["code"] in {"E", "B"})
    coil = sum(1 for item in residues if item["code"] not in {"H", "G", "I", "E", "B"})
    return {
        "residues": residues,
        "n_residues": len(residues),
        "n_helix": helix,
        "n_sheet": sheet,
        "n_coil": coil,
        "format": "legacy_dssp_columns",
        "note": (
            "Counts follow Kabsch-Sander DSSP letters in the parsed file. "
            "This is assignment, not sequence prediction."
        ),
    }


def map_dssp_assignments(
    dssp_residues: Sequence[Mapping[str, Any]],
    structure_residues: Sequence[Mapping[str, Any]],
) -> dict:
    """Mapeia linhas DSSP a residuos da estrutura. Sem residuo, sem assignment.

    Args:
        dssp_residues: Saida de parse_dssp_legacy_output()['residues'].
        structure_residues: Residuos depositados (auth_seq_id, chain, coords).

    Returns:
        Dict rows, n_mapped, n_missing_structure, n_unmapped_dssp.

    Raises:
        Nenhum.

    Nota biologica:
        DSSP nao preenche gaps da estrutura. Residuo sem coordenadas ou sem
        linha DSSP fica sem codigo, nao como coil inventado.
    """
    by_key: Dict[tuple, dict] = {}
    for item in dssp_residues:
        chain = str(item.get("chain_id") or "").strip() or "A"
        seq_id = item.get("auth_seq_id")
        if seq_id is None:
            continue
        by_key[(chain, int(seq_id))] = dict(item)
    rows: List[dict] = []
    mapped = 0
    missing = 0
    for residue in structure_residues:
        chain = str(
            residue.get("chain_id")
            or residue.get("auth_asym_id")
            or residue.get("label_asym_id")
            or ""
        ).strip() or "A"
        seq_id = residue.get("auth_seq_id")
        has_coords = bool(residue.get("has_coordinates"))
        if seq_id is None or not has_coords:
            rows.append(
                {
                    "chain_id": chain,
                    "auth_seq_id": seq_id,
                    "structure_aa": residue.get("one_letter"),
                    "dssp_code": None,
                    "mapped": False,
                    "reason": "Residue has no deposited coordinates for DSSP mapping.",
                }
            )
            missing += 1
            continue
        hit = by_key.pop((chain, int(seq_id)), None)
        if hit is None:
            rows.append(
                {
                    "chain_id": chain,
                    "auth_seq_id": int(seq_id),
                    "structure_aa": residue.get("one_letter"),
                    "dssp_code": None,
                    "mapped": False,
                    "reason": "No DSSP row for this residue.",
                }
            )
            missing += 1
            continue
        rows.append(
            {
                "chain_id": chain,
                "auth_seq_id": int(seq_id),
                "structure_aa": residue.get("one_letter"),
                "dssp_aa": hit.get("aa"),
                "dssp_code": hit.get("code"),
                "mapped": True,
                "reason": "",
            }
        )
        mapped += 1
    unmapped_dssp = [
        {"chain_id": key[0], "auth_seq_id": key[1], "aa": item.get("aa"), "code": item.get("code")}
        for key, item in by_key.items()
    ]
    return {
        "status": "EXPERIMENTAL" if mapped else "UNMAPPED",
        "method": "DSSP residue map by chain + auth_seq_id",
        "n_mapped": mapped,
        "n_missing_structure": missing,
        "n_unmapped_dssp": len(unmapped_dssp),
        "rows": rows,
        "unmapped_dssp": unmapped_dssp,
        "disclaimer": (
            "Geometric DSSP assignment mapped onto deposited residues. "
            "Not a sequence secondary-structure predictor."
        ),
    }


def parse_stride_output(text: str) -> dict:
    """Parseia linhas ASG do STRIDE. Nao inventa H/E/C.

    Args:
        text: Stdout do binario STRIDE.

    Returns:
        Dict residues, n_helix, n_sheet, n_coil.

    Raises:
        StructureError: PARSING_ERROR se nao houver linhas ASG.

    Nota biologica:
        STRIDE (Frishman and Argos 1995) e assignment geometrico. Sem o
        executavel oficial o parser so corre sobre texto injetado.
    """
    residues: List[dict] = []
    for line in str(text or "").splitlines():
        if not line.startswith("ASG "):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        code = parts[5]
        try:
            seqnum = int(parts[3])
        except ValueError:
            continue
        residues.append(
            {
                "comp_id": parts[1],
                "chain_id": parts[2],
                "auth_seq_id": seqnum,
                "aa": residue_one_letter(parts[1]),
                "code": code,
            }
        )
    if not residues:
        raise StructureError(
            "STRIDE output has no ASG residue rows. Return code is not assignment.",
            "PARSING_ERROR",
        )
    helix = sum(1 for item in residues if item["code"] in {"H", "G", "I"})
    sheet = sum(1 for item in residues if item["code"] in {"E", "B", "b"})
    coil = sum(1 for item in residues if item["code"] not in {"H", "G", "I", "E", "B", "b"})
    return {
        "residues": residues,
        "n_residues": len(residues),
        "n_helix": helix,
        "n_sheet": sheet,
        "n_coil": coil,
        "format": "stride_asg",
        "note": "STRIDE ASG codes from the tool output, not a house geometric cartoon.",
    }


def assign_secondary_structure_stride(
    structure_text: str = "",
    *,
    run_fn=None,
    output_text: Optional[str] = None,
) -> dict:
    """Executa STRIDE real ou parseia output injetado. Ausente = UNAVAILABLE.

    Args:
        structure_text: PDB/mmCIF para o binario.
        run_fn: subprocess.run injetavel.
        output_text: Texto STRIDE para testes sem binario.

    Returns:
        Envelope EXPERIMENTAL, UNAVAILABLE ou ERROR.

    Raises:
        Nenhum.

    Nota biologica:
        Nao reproduz o algoritmo STRIDE internamente.
    """
    if output_text is not None:
        try:
            parsed = parse_stride_output(output_text)
        except StructureError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN STRIDE",
                payload={"reason": str(exc), "n_residues": 0},
                status="ERROR",
                algorithm="STRIDE output parser",
                parameters={"injected_output": True},
                source="injected STRIDE text",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": exc.category}
        return provenance.analysis_envelope(
            module="PROTEIN STRIDE",
            payload=parsed,
            status="EXPERIMENTAL",
            algorithm="STRIDE (Frishman and Argos 1995) assignment parser",
            parameters={"injected_output": True},
            source="injected STRIDE text",
        ) | {"status": "EXPERIMENTAL", "available": True, "tool": "STRIDE", "version": ""}

    avail = stride_availability()
    if not avail.get("available"):
        return provenance.analysis_envelope(
            module="PROTEIN STRIDE",
            payload={"reason": avail.get("reason"), "n_residues": 0},
            status="UNAVAILABLE",
            algorithm="none",
            parameters={},
            source="HelixScope STRIDE",
        ) | {
            "status": "UNAVAILABLE",
            "available": False,
            "tool": "STRIDE",
            "version": "",
            "reason": avail.get("reason"),
        }
    text = str(structure_text or "")
    if not text.strip():
        return provenance.analysis_envelope(
            module="PROTEIN STRIDE",
            payload={"reason": "No structure text to assign.", "n_residues": 0},
            status="ERROR",
            algorithm="STRIDE",
            parameters={},
            source="HelixScope STRIDE",
        ) | {"status": "ERROR", "available": False, "reason": "No structure text to assign."}
    executable = tool_detection.resolve_allowlisted_executable(
        STRIDE_NAMES,
        extra_file_candidates=tool_paths.candidates_for(STRIDE_NAMES),
        env_var=STRIDE_ENV,
    )
    tmpdir = tempfile.mkdtemp(prefix="helixscope_stride_")
    try:
        infile = tmpdir + "/input.pdb"
        with open(infile, "w", encoding="utf-8") as handle:
            handle.write(text)
        runner = run_fn or subprocess.run
        try:
            completed = runner(
                [executable, infile],
                capture_output=True,
                text=True,
                timeout=DSSP_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmpdir,
            )
        except subprocess.TimeoutExpired as exc:
            return provenance.analysis_envelope(
                module="PROTEIN STRIDE",
                payload={"reason": f"STRIDE exceeded {int(DSSP_TIMEOUT_S)}s.", "n_residues": 0},
                status="ERROR",
                algorithm="STRIDE",
                parameters={"timeout_s": DSSP_TIMEOUT_S},
                source="HelixScope STRIDE",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": "TIMEOUT"}
        except OSError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN STRIDE",
                payload={"reason": str(exc), "n_residues": 0},
                status="ERROR",
                algorithm="STRIDE",
                parameters={},
                source="HelixScope STRIDE",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": "TOOL_ERROR"}
        stdout = str(getattr(completed, "stdout", "") or "")
        try:
            parsed = parse_stride_output(stdout)
        except StructureError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN STRIDE",
                payload={"reason": str(exc), "n_residues": 0},
                status="ERROR",
                algorithm="STRIDE",
                parameters={},
                source="HelixScope STRIDE",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": "PARSING_ERROR"}
        if run_fn is None:
            engine_validation.record_live_validation(
                "STRIDE",
                ok=True,
                version=str(avail.get("version") or ""),
                details={"n_residues": int(parsed.get("n_residues") or 0)},
            )
        return provenance.analysis_envelope(
            module="PROTEIN STRIDE",
            payload=parsed,
            status="EXPERIMENTAL",
            algorithm="STRIDE assignment",
            parameters={"returncode": int(completed.returncode) if completed.returncode is not None else 1},
            source="HelixScope STRIDE",
        ) | {
            "status": "EXPERIMENTAL",
            "available": True,
            "tool": "STRIDE",
            "version": str(avail.get("version") or ""),
        }
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def assign_secondary_structure_dssp(
    structure_text: str = "",
    *,
    run_fn=None,
    output_text: Optional[str] = None,
) -> dict:
    """Executa DSSP de forma segura ou parseia output injetado.

    Args:
        structure_text: mmCIF/PDB para o binario (ficheiro temporario).
        run_fn: subprocess.run injetavel (testes).
        output_text: Se fornecido, so o parser corre (sem binario).

    Returns:
        Envelope EXPERIMENTAL assignment ou UNAVAILABLE/ERROR.

    Raises:
        Nenhum. Falhas voltam classificadas no dict.

    Nota biologica:
        Implementacao chama mkdssp/dssp real. Nao reproduz o algoritmo DSSP
        internamente. Output vazio com returncode 0 e PARSING_ERROR.
    """
    if output_text is not None:
        try:
            parsed = parse_dssp_legacy_output(output_text)
        except StructureError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN DSSP",
                payload={"reason": str(exc), "n_residues": 0},
                status="ERROR",
                algorithm="DSSP output parser",
                parameters={"injected_output": True},
                source="injected DSSP text",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": exc.category}
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload=parsed,
            status="EXPERIMENTAL",
            algorithm="DSSP (Kabsch and Sander) assignment parser",
            parameters={"injected_output": True},
            source="injected DSSP text",
        ) | {"status": "EXPERIMENTAL", "available": True, "tool": "DSSP", "version": ""}

    avail = dssp_availability()
    if not avail.get("available"):
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={"reason": avail.get("reason"), "n_residues": 0},
            status="UNAVAILABLE",
            algorithm="none",
            parameters={},
            source="HelixScope DSSP",
        ) | {
            "status": "UNAVAILABLE",
            "available": False,
            "tool": "DSSP",
            "version": "",
            "reason": avail.get("reason"),
        }
    text = str(structure_text or "")
    if not text.strip():
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={"reason": "No structure text to assign.", "n_residues": 0},
            status="ERROR",
            algorithm="DSSP",
            parameters={},
            source="HelixScope DSSP",
        ) | {"status": "ERROR", "available": False, "reason": "No structure text to assign."}

    structure_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    cache_identity = hashlib.sha256(
        f"{structure_sha}|dssp|{avail.get('version') or ''}|{provenance.HELIXSCOPE_VERSION}".encode("utf-8")
    ).hexdigest()
    cached = _DSSP_CACHE.get(cache_identity)
    if isinstance(cached, dict):
        copied = copy.deepcopy(cached)
        copied["cache_status"] = "cached"
        return copied

    executable = tool_detection.resolve_allowlisted_executable(
        DSSP_NAMES,
        extra_file_candidates=tool_paths.candidates_for(DSSP_NAMES),
        env_var=DSSP_ENV,
    )
    if executable is None:
        if run_fn is None:
            return provenance.analysis_envelope(
                module="PROTEIN DSSP",
                payload={"reason": avail.get("reason"), "n_residues": 0},
                status="UNAVAILABLE",
                algorithm="none",
                parameters={},
                source="HelixScope DSSP",
            ) | {
                "status": "UNAVAILABLE",
                "available": False,
                "tool": "DSSP",
                "version": "",
                "reason": "DSSP executable is not allowlisted.",
            }
        executable = "mkdssp"
    tmpdir = tempfile.mkdtemp(prefix="helixscope_dssp_")
    try:
        infile = tmpdir + "/input.cif"
        with open(infile, "w", encoding="utf-8") as handle:
            handle.write(text)
        runner = run_fn or subprocess.run
        try:
            completed = runner(
                [executable, "--output-format=dssp", infile],
                capture_output=True,
                text=True,
                timeout=DSSP_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmpdir,
            )
        except subprocess.TimeoutExpired as exc:
            return provenance.analysis_envelope(
                module="PROTEIN DSSP",
                payload={"reason": f"DSSP exceeded {int(DSSP_TIMEOUT_S)}s.", "n_residues": 0},
                status="ERROR",
                algorithm="DSSP",
                parameters={"timeout_s": DSSP_TIMEOUT_S},
                source="HelixScope DSSP",
            ) | {"status": "ERROR", "available": False, "reason": str(exc), "category": "TIMEOUT"}
        except OSError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN DSSP",
                payload={"reason": str(exc), "n_residues": 0},
                status="ERROR",
                algorithm="DSSP",
                parameters={},
                source="HelixScope DSSP",
            ) | {"status": "ERROR", "available": False, "reason": str(exc)}
        stdout = str(getattr(completed, "stdout", "") or "")
        returncode = int(completed.returncode) if completed.returncode is not None else 1
        try:
            parsed = parse_dssp_legacy_output(stdout)
        except StructureError as exc:
            return provenance.analysis_envelope(
                module="PROTEIN DSSP",
                payload={
                    "reason": str(exc),
                    "returncode": returncode,
                    "n_residues": 0,
                },
                status="ERROR",
                algorithm="DSSP",
                parameters={"returncode": returncode},
                source="HelixScope DSSP",
            ) | {
                "status": "ERROR",
                "available": False,
                "reason": str(exc),
                "category": "PARSING_ERROR",
                "returncode": returncode,
            }
        envelope = provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload=parsed,
            status="EXPERIMENTAL",
            algorithm="mkdssp/dssp assignment",
            parameters={
                "returncode": returncode,
                "structure_sha256": structure_sha,
            },
            source="HelixScope DSSP",
        ) | {
            "status": "EXPERIMENTAL",
            "available": True,
            "tool": "DSSP",
            "version": str(avail.get("version") or ""),
            "returncode": returncode,
            "cache_status": "live",
            "cache_key": cache_identity,
        }
        if run_fn is None:
            engine_validation.record_live_validation(
                "DSSP",
                ok=True,
                version=str(avail.get("version") or ""),
                details={"n_residues": int(parsed.get("n_residues") or 0)},
            )
        _DSSP_CACHE[cache_identity] = copy.deepcopy(envelope)
        if len(_DSSP_CACHE) > MAX_CACHE * 4:
            oldest = next(iter(_DSSP_CACHE))
            _DSSP_CACHE.pop(oldest, None)
        return envelope
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def assign_secondary_structure_dssp_remote(
    structure_text: str = "",
    *,
    urlopen_fn=None,
    coordinate_kind: str = "experimental",
) -> dict:
    """POST coordinates to the official PDB-REDO DSSP API and parse legacy output.

    Args:
        structure_text: mmCIF or PDB text of the structure being assigned.
        urlopen_fn: HTTP client injectavel (testes). None usa a rede.
        coordinate_kind: experimental ou predicted. Predicted coordinates are
            not labeled EXPERIMENTAL.

    Returns:
        Envelope EXPERIMENTAL/COMPUTED assignment, UNAVAILABLE, or ERROR.
        Never LIVE_VALIDATED. Never sequence prediction.

    Raises:
        Nenhum.

    Nota biologica:
        Implementacao envia coordenadas ao mkdssp do PDB-REDO. Nao reproduz o
        algoritmo DSSP internamente. Coordenadas predicted permanecem predicted.
    """
    text = str(structure_text or "")
    if not text.strip():
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={"reason": "No structure text to assign.", "n_residues": 0},
            status="ERROR",
            algorithm="PDB-REDO DSSP API",
            parameters={},
            source=DSSP_REMOTE_TOOL,
        ) | {"status": "ERROR", "available": False, "reason": "No structure text to assign."}
    encoded = text.encode("utf-8")
    if len(encoded) > MAX_STRUCTURE_BYTES:
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={
                "reason": f"Structure exceeds {MAX_STRUCTURE_BYTES:,} bytes.",
                "n_residues": 0,
            },
            status="ERROR",
            algorithm="PDB-REDO DSSP API",
            parameters={},
            source=DSSP_REMOTE_TOOL,
        ) | {
            "status": "ERROR",
            "available": False,
            "reason": f"Structure exceeds {MAX_STRUCTURE_BYTES:,} bytes.",
            "category": "RESOURCE_LIMIT",
        }
    kind = str(coordinate_kind or "experimental").strip().lower()
    if kind not in {"experimental", "predicted"}:
        kind = "experimental"
    assignment_status = "EXPERIMENTAL" if kind == "experimental" else "COMPUTED"
    structure_sha = hashlib.sha256(encoded).hexdigest()
    cache_identity = hashlib.sha256(
        f"{structure_sha}|dssp-pdb-redo|{kind}|{provenance.HELIXSCOPE_VERSION}".encode("utf-8")
    ).hexdigest()
    cached = _DSSP_CACHE.get(cache_identity)
    if isinstance(cached, dict):
        copied = copy.deepcopy(cached)
        copied["cache_status"] = "cached"
        return copied
    try:
        stdout = _http_text(
            PDB_REDO_DSSP_DO,
            method="POST",
            multipart_fields={
                "format": ("", "dssp"),
                "data": ("input.cif", text),
            },
            urlopen_fn=urlopen_fn,
            limit=MAX_STRUCTURE_BYTES,
            timeout_s=DSSP_REMOTE_TIMEOUT_S,
        )
    except StructureError as exc:
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={"reason": str(exc), "n_residues": 0},
            status="ERROR",
            algorithm="PDB-REDO DSSP API",
            parameters={"endpoint": PDB_REDO_DSSP_DO},
            source=DSSP_REMOTE_TOOL,
        ) | {
            "status": "ERROR",
            "available": False,
            "reason": str(exc),
            "category": exc.category,
            "engine_location": "remote",
            "tool": DSSP_REMOTE_TOOL,
        }
    try:
        parsed = parse_dssp_legacy_output(stdout)
    except StructureError as exc:
        return provenance.analysis_envelope(
            module="PROTEIN DSSP",
            payload={"reason": str(exc), "n_residues": 0},
            status="ERROR",
            algorithm="PDB-REDO DSSP API parser",
            parameters={"endpoint": PDB_REDO_DSSP_DO},
            source=DSSP_REMOTE_TOOL,
        ) | {
            "status": "ERROR",
            "available": False,
            "reason": str(exc),
            "category": exc.category,
            "engine_location": "remote",
            "tool": DSSP_REMOTE_TOOL,
        }
    envelope = provenance.analysis_envelope(
        module="PROTEIN DSSP",
        payload=parsed,
        status=assignment_status,
        algorithm="mkdssp assignment via PDB-REDO HTTP API",
        parameters={
            "endpoint": PDB_REDO_DSSP_DO,
            "structure_sha256": structure_sha,
            "coordinate_kind": kind,
            "output_format": "legacy_dssp",
        },
        source=DSSP_REMOTE_TOOL,
    ) | {
        "status": assignment_status,
        "available": True,
        "tool": DSSP_REMOTE_TOOL,
        "version": "mkdssp via PDB-REDO",
        "engine_location": "remote",
        "coordinate_kind": kind,
        "not_sequence_prediction": True,
        "not_local_mkdssp": True,
        "cache_status": "live",
        "cache_key": cache_identity,
        "reason": (
            "Secondary structure assigned by PDB-REDO mkdssp from the submitted "
            "atomic coordinates. This is not sequence prediction and not a "
            "local LIVE_VALIDATED mkdssp binary."
        ),
    }
    if urlopen_fn is None:
        engine_validation.record_remote_validation(
            DSSP_REMOTE_TOOL,
            ok=True,
            version="mkdssp via PDB-REDO",
            details={
                "n_residues": int(parsed.get("n_residues") or 0),
                "endpoint": PDB_REDO_DSSP_DO,
                "structure_sha256": structure_sha,
            },
        )
    _DSSP_CACHE[cache_identity] = copy.deepcopy(envelope)
    if len(_DSSP_CACHE) > MAX_CACHE * 4:
        oldest = next(iter(_DSSP_CACHE))
        _DSSP_CACHE.pop(oldest, None)
    return envelope


def ca_contact_map(
    result: Mapping[str, Any],
    *,
    threshold_angstrom: float = 8.0,
) -> dict:
    """Pares CA-CA com distancia Euclidiana <= limiar. Coordenadas depositadas.

    Args:
        result: Envelope estrutural com residue_mapping.
        threshold_angstrom: Distancia maxima em Angstrom (declarada).

    Returns:
        Dict pairs, threshold, method, n_residues_with_ca. Sem score
        proprietario.

    Raises:
        StructureError: INVALID_INPUT se o limiar nao for finito e positivo.

    Nota biologica:
        Contacto geometrico CA-CA, nao um potencial estatistico (e.g. CASP).
        8 A e um limiar comum de vizinhança, nao uma verdade biologica unica.
    """
    try:
        cutoff = float(threshold_angstrom)
    except (TypeError, ValueError) as exc:
        raise StructureError("Contact threshold must be a number.", "INVALID_INPUT") from exc
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise StructureError("Contact threshold must be a finite positive length.", "INVALID_INPUT")
    points: List[tuple[int, float, float, float]] = []
    for item in list(result.get("residue_mapping") or []):
        if item.get("query_index_0based") is None:
            continue
        ca = item.get("ca") or {}
        if ca.get("x") is None or ca.get("y") is None or ca.get("z") is None:
            continue
        points.append(
            (
                int(item["query_index_0based"]),
                float(ca["x"]),
                float(ca["y"]),
                float(ca["z"]),
            )
        )
    pairs: List[dict] = []
    for i, (idx_a, xa, ya, za) in enumerate(points):
        for idx_b, xb, yb, zb in points[i + 1 :]:
            dist = math.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2)
            if dist <= cutoff:
                pairs.append(
                    {
                        "i": idx_a,
                        "j": idx_b,
                        "distance_angstrom": dist,
                    }
                )
    return {
        "status": "COMPUTED" if points else "UNAVAILABLE",
        "method": "CA-CA Euclidean distance on deposited coordinates",
        "threshold_angstrom": cutoff,
        "n_residues_with_ca": len(points),
        "n_pairs": len(pairs),
        "pairs": pairs,
        "disclaimer": (
            "Geometric neighborhood, not a published statistical contact potential. "
            "Absence of CA coordinates yields no pair, not a fabricated score. "
            "A contact is not proof of a biological interaction."
        ),
    }


def interchain_contacts(
    parsed: Mapping[str, Any],
    *,
    threshold_angstrom: float = 8.0,
    atom_name: str = "CA",
) -> dict:
    """Contactos entre chains por distancia Euclidiana do atomo escolhido.

    Args:
        parsed: Estrutura parseada (atoms).
        threshold_angstrom: Distancia maxima declarada.
        atom_name: Nome do atomo (CA para proteina, P para acido nucleico).

    Returns:
        Dict pairs, threshold, method. 0 pares e COMPUTED, nao UNAVAILABLE.

    Raises:
        StructureError: INVALID_INPUT se o limiar nao for finito e positivo.

    Nota biologica:
        Contacto geometrico inter-chain nao prova interface biologica.
    """
    try:
        cutoff = float(threshold_angstrom)
    except (TypeError, ValueError) as exc:
        raise StructureError("Contact threshold must be a number.", "INVALID_INPUT") from exc
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise StructureError("Contact threshold must be a finite positive length.", "INVALID_INPUT")
    wanted = str(atom_name or "CA").strip().upper()
    points: List[tuple[str, int, str, float, float, float]] = []
    for atom in list(parsed.get("atoms") or []):
        if str(atom.get("group") or "") != "ATOM":
            continue
        if str(atom.get("atom_name") or "").strip().upper() != wanted:
            continue
        chain = str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        seq_id = atom.get("auth_seq_id")
        if seq_id is None or atom.get("x") is None:
            continue
        points.append(
            (
                chain,
                int(seq_id),
                str(atom.get("comp_id") or ""),
                float(atom["x"]),
                float(atom["y"]),
                float(atom["z"]),
            )
        )
    pairs: List[dict] = []
    for i, (chain_a, seq_a, comp_a, xa, ya, za) in enumerate(points):
        for chain_b, seq_b, comp_b, xb, yb, zb in points[i + 1 :]:
            if chain_a == chain_b:
                continue
            dist = math.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2)
            if dist <= cutoff:
                pairs.append(
                    {
                        "chain_a": chain_a,
                        "auth_seq_id_a": seq_a,
                        "comp_a": comp_a,
                        "chain_b": chain_b,
                        "auth_seq_id_b": seq_b,
                        "comp_b": comp_b,
                        "atom": wanted,
                        "distance_angstrom": dist,
                    }
                )
    return {
        "status": "COMPUTED",
        "method": f"{wanted}-{wanted} Euclidean distance between different chains",
        "threshold_angstrom": cutoff,
        "atom_name": wanted,
        "n_atoms": len(points),
        "n_pairs": len(pairs),
        "pairs": pairs,
        "disclaimer": (
            "Geometric inter-chain neighborhood. Not a biological interface "
            "proof and not a docking score."
        ),
    }


BONDI_RADII_ANGSTROM: MappingProxyType = MappingProxyType(
    {
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "S": 1.80,
        "P": 1.80,
        "SE": 1.90,
    }
)
SASA_PROBE_RADIUS: float = 1.4
SASA_SPHERE_POINTS: int = 92
MAX_SASA_ATOMS: int = 4_000


def _golden_sphere_points(count: int) -> List[Tuple[float, float, float]]:
    points: List[Tuple[float, float, float]] = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(int(count)):
        z = 1.0 - (2.0 * index + 1.0) / float(count)
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * index
        points.append((math.cos(theta) * radius, math.sin(theta) * radius, z))
    return points


def _sasa_neighbor_lists(atoms: Sequence[Mapping[str, Any]], probe: float) -> list[list[int]]:
    """Pares que podem ocluir um ponto da esfera (mesmo criterio geometrico)."""
    neighbors: list[list[int]] = [[] for _ in range(len(atoms))]
    count = len(atoms)
    for i in range(count):
        xi = float(atoms[i]["x"])
        yi = float(atoms[i]["y"])
        zi = float(atoms[i]["z"])
        ri = float(atoms[i]["radius"])
        for j in range(i + 1, count):
            dx = xi - float(atoms[j]["x"])
            dy = yi - float(atoms[j]["y"])
            dz = zi - float(atoms[j]["z"])
            limit = ri + float(atoms[j]["radius"]) + 2.0 * probe
            if dx * dx + dy * dy + dz * dz <= limit * limit:
                neighbors[i].append(j)
                neighbors[j].append(i)
    return neighbors


def shrake_rupley_sasa(
    parsed: Mapping[str, Any],
    *,
    probe_radius: float = SASA_PROBE_RADIUS,
    n_points: int = SASA_SPHERE_POINTS,
) -> dict:
    """Area acessivel ao solvente (Shrake e Rupley 1973) sobre atomos depositados.

    Args:
        parsed: Estrutura parseada com atoms.
        probe_radius: Raio da sonda (A), tipicamente 1.4 A (agua).
        n_points: Pontos na esfera por atomo (92 no paper original).

    Returns:
        Envelope COMPUTED com sasa por residuo, ou RESOURCE_LIMIT/UNAVAILABLE.

    Raises:
        StructureError: INVALID_INPUT se parametros nao forem finitos positivos.

    Nota biologica:
        SASA numerica nao e malha de superficie molecular (MSMS/SES). Nao e DSSP ACC.
        Raios atomicos: Bondi 1964. Implementacao do algoritmo publicado de
        amostragem esferica; nao e um wrapper de MSMS.
    """
    try:
        probe = float(probe_radius)
        samples = int(n_points)
    except (TypeError, ValueError) as exc:
        raise StructureError("SASA probe radius and point count must be numeric.", "INVALID_INPUT") from exc
    if not math.isfinite(probe) or probe <= 0 or samples < 4:
        raise StructureError("SASA probe radius must be positive and n_points >= 4.", "INVALID_INPUT")
    atoms: List[dict] = []
    for atom in list(parsed.get("atoms") or []):
        if str(atom.get("group") or "") != "ATOM":
            continue
        element = str(atom.get("element") or "").strip().upper()
        if not element:
            name = str(atom.get("atom_name") or "").strip().upper()
            element = name[:2] if name[:2] in BONDI_RADII_ANGSTROM else name[:1]
        if element == "H":
            continue
        radius = BONDI_RADII_ANGSTROM.get(element)
        if radius is None:
            continue
        if atom.get("x") is None:
            continue
        atoms.append(
            {
                "element": element,
                "radius": float(radius),
                "x": float(atom["x"]),
                "y": float(atom["y"]),
                "z": float(atom["z"]),
                "chain_id": str(atom.get("auth_asym_id") or atom.get("label_asym_id") or ""),
                "auth_seq_id": atom.get("auth_seq_id"),
                "atom_name": str(atom.get("atom_name") or ""),
            }
        )
    if not atoms:
        return provenance.analysis_envelope(
            module="PROTEIN SASA",
            payload={"reason": "No heavy ATOM coordinates for Shrake-Rupley.", "n_atoms": 0},
            status="UNAVAILABLE",
            algorithm="none",
            parameters={"probe_radius": probe, "n_points": samples},
            source="HelixScope Shrake-Rupley",
        ) | {"status": "UNAVAILABLE", "sasa_angstrom2": None, "n_atoms": 0}
    if len(atoms) > MAX_SASA_ATOMS:
        raise StructureError(
            f"SASA refuses more than {MAX_SASA_ATOMS} heavy atoms (RESOURCE_LIMIT).",
            "RESOURCE_LIMIT",
        )
    structure_hash = str(parsed.get("content_hash") or parsed.get("structure_hash") or "")
    cache_identity = hashlib.sha256(
        f"{structure_hash}|{probe}|{samples}|{SASA_ALGO_VERSION}|{len(atoms)}".encode("utf-8")
    ).hexdigest()
    cached = _SASA_CACHE.get(cache_identity)
    if isinstance(cached, dict) and structure_hash:
        copied = copy.deepcopy(cached)
        copied["cache_status"] = "cached"
        return copied
    sphere = _golden_sphere_points(samples)
    neighbor_lists = _sasa_neighbor_lists(atoms, probe)
    per_residue: Dict[tuple, float] = {}
    total = 0.0
    for index, atom in enumerate(atoms):
        extended = atom["radius"] + probe
        area_per_point = 4.0 * math.pi * extended * extended / float(samples)
        accessible = 0
        for sx, sy, sz in sphere:
            px = atom["x"] + sx * extended
            py = atom["y"] + sy * extended
            pz = atom["z"] + sz * extended
            buried = False
            for other_index in neighbor_lists[index]:
                other = atoms[other_index]
                limit = other["radius"] + probe
                dx = px - other["x"]
                dy = py - other["y"]
                dz = pz - other["z"]
                if dx * dx + dy * dy + dz * dz <= limit * limit:
                    buried = True
                    break
            if not buried:
                accessible += 1
        sasa = area_per_point * accessible
        total += sasa
        key = (atom["chain_id"], atom["auth_seq_id"])
        per_residue[key] = per_residue.get(key, 0.0) + sasa
    residues = [
        {
            "chain_id": chain,
            "auth_seq_id": seq_id,
            "sasa_angstrom2": round(value, 3),
        }
        for (chain, seq_id), value in sorted(per_residue.items(), key=lambda item: (item[0][0], item[0][1] or 0))
    ]
    payload = {
        "sasa_angstrom2": round(total, 3),
        "n_atoms": len(atoms),
        "n_residues": len(residues),
        "residues": residues,
        "probe_radius_angstrom": probe,
        "n_points": samples,
        "radii": "Bondi 1964",
        "algorithm": "Shrake and Rupley 1973 sphere sampling",
        "neighbor_list": True,
        "algorithm_version": SASA_ALGO_VERSION,
        "not_msms": True,
        "not_ses_mesh": True,
        "disclaimer": (
            "Solvent-accessible surface area, not a triangulated molecular surface."
        ),
    }
    result = provenance.analysis_envelope(
        module="PROTEIN SASA",
        payload=payload,
        status="COMPUTED",
        algorithm="Shrake and Rupley 1973",
        parameters={
            "probe_radius_angstrom": probe,
            "n_points": samples,
            "radii": "Bondi 1964",
            "neighbor_list": True,
            "algorithm_version": SASA_ALGO_VERSION,
        },
        source="HelixScope Shrake-Rupley",
    ) | {"status": "COMPUTED", **payload, "cache_status": "live"}
    if structure_hash:
        _SASA_CACHE[cache_identity] = copy.deepcopy(result)
        if len(_SASA_CACHE) > MAX_CACHE * 4:
            oldest = next(iter(_SASA_CACHE))
            _SASA_CACHE.pop(oldest, None)
    return result


def classify_structure_identifier(text: str) -> dict:
    """Classifica PDB ID, UniProt accession ou identificador AlphaFold.

    Args:
        text: Identificador bruto.

    Returns:
        Dict kind, canonical, uniprot (quando extraivel).

    Raises:
        Nenhum. Identificadores invalidos devolvem kind=invalid.
    """
    raw = str(text or "").strip()
    if not raw or len(raw) > 40:
        return {"kind": "invalid", "canonical": "", "uniprot": ""}
    if "://" in raw or "/" in raw or "\\" in raw or ".." in raw:
        return {"kind": "invalid", "canonical": "", "uniprot": ""}
    compact = raw.replace(" ", "")
    af = ALPHAFOLD_ID_PATTERN.match(compact)
    if af:
        uniprot = af.group(1).upper()
        if not UNIPROT_PATTERN.match(uniprot):
            return {"kind": "invalid", "canonical": "", "uniprot": ""}
        fragment = compact.split("-F")[-1]
        return {
            "kind": "alphafold",
            "canonical": f"AF-{uniprot}-F{fragment.upper()}",
            "uniprot": uniprot,
        }
    if PDB_ID_PATTERN.match(compact):
        return {"kind": "pdb", "canonical": compact.upper(), "uniprot": ""}
    accession = compact.split(".")[0]
    if UNIPROT_PATTERN.match(accession):
        return {"kind": "uniprot", "canonical": accession.upper(), "uniprot": accession.upper()}
    return {"kind": "invalid", "canonical": "", "uniprot": ""}


def url_is_allowed(url: str) -> bool:
    """Confere host e prefixo de caminho contra a allowlist estrutural.

    Args:
        url: URL absoluta.

    Returns:
        True somente para HTTPS nos endpoints documentados.

    Raises:
        Nenhum.
    """
    parsed = urlparse(str(url or ""))
    if parsed.scheme != "https":
        return False
    host = (parsed.hostname or "").lower()
    path = parsed.path or ""
    for allowed_host, prefix in ALLOWED_URLS:
        if host == allowed_host and path.startswith(prefix):
            return True
    return False


def remote_source_status() -> dict:
    """Declara as fontes estruturais integradas. Nao testa a rede.

    Args:
        Nenhum.

    Returns:
        Dict por fonte com endpoint e kind esperado.

    Raises:
        Nenhum.
    """
    return {
        "rcsb_pdb": {
            "available": True,
            "source": "RCSB PDB",
            "kind_when_experimental": "experimental",
            "endpoints": [RCSB_DATA, RCSB_SEARCH, RCSB_FILES],
            "method": "RCSB Data API, Search API (MMseqs2), Files API mmCIF",
        },
        "alphafold_db": {
            "available": True,
            "source": "AlphaFold DB",
            "kind": "predicted",
            "endpoints": [ALPHAFOLD_API],
            "method": "AlphaFold DB prediction API and mmCIF files",
        },
        "uniprot": {
            "available": True,
            "source": "UniProtKB",
            "kind": "cross-reference",
            "endpoints": [UNIPROT_REST],
            "method": "UniProt REST uniprotkb JSON cross-references",
        },
        "local_3d_renderer": {
            "available": True,
            "engine": "plotly",
            "library": "plotly.graph_objects.Scatter3d",
            "reason": (
                "Plotly 3D visualizes validated structure_3d_input coordinates after "
                "an explicit View 3D Structure action. It does not fetch PDB, "
                "AlphaFold, NCBI or UniProt data and does not invent coordinates."
            ),
        },
        "pdb_redo_dssp": {
            "available": True,
            "source": "PDB-REDO DSSP",
            "kind": "remote_assignment",
            "endpoints": [PDB_REDO_DSSP_DO, PDB_REDO_DSSP_DB],
            "method": (
                "Official PDB-REDO mkdssp HTTP API. Geometric assignment of "
                "submitted coordinates. Not sequence prediction. Not local mkdssp."
            ),
        },
    }


def source_disclaimer(kind: str, source: str) -> str:
    """Texto especifico da fonte, em cima de provenance.structure_disclaimer.

    Args:
        kind: experimental, predicted, retrieved, illustrative ou unavailable.
        source: Nome da base.

    Returns:
        Frases em ingles.

    Raises:
        ValueError: Se kind nao for reconhecido.
    """
    base = provenance.structure_disclaimer(kind)
    src = str(source or "").strip()
    extra = ""
    if src == "RCSB PDB" and kind == "experimental":
        extra = " Experimental structure retrieved from RCSB PDB."
    elif src == "AlphaFold DB" and kind == "predicted":
        extra = " Predicted structure from AlphaFold DB. pLDDT is not a probability that the fold is correct."
    elif src:
        extra = f" Source: {src}."
    return base + extra


def cache_key(*, source: str, identifier: str, chain: str, parameters: Mapping[str, Any]) -> str:
    """Chave de cache de estrutura.

    Args:
        source: Banco.
        identifier: PDB/AF/UniProt.
        chain: Chain selecionada ou vazio.
        parameters: Parametros efetivos.

    Returns:
        SHA-256 hexadecimal.

    Raises:
        Nenhum.
    """
    params = dict(parameters)
    payload = json.dumps(
        {
            "source": source,
            "identifier": identifier,
            "chain": chain,
            "molecule": str(params.get("molecule") or params.get("cache_molecule") or ""),
            "parameters": params,
        },
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def mark_cached_result(result: Mapping[str, Any]) -> dict:
    """Copia um resultado marcando cache, preservando o timestamp original.

    Args:
        result: Envelope anterior.

    Returns:
        Copia com cache_status=cached.

    Raises:
        Nenhum.
    """
    copied = dict(result)
    copied["cache_status"] = "cached"
    copied["served_from_cache_at_utc"] = provenance.utc_now()
    return copied


def parse_mmcif(text: str) -> dict:
    """Interpreta mmCIF o suficiente para atomos, chains e metadados de metodo.

    Args:
        text: Conteudo mmCIF.

    Returns:
        Dict format, entry_id, method, resolution, chains, atoms, models.

    Raises:
        StructureError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    if not isinstance(text, str) or not text.strip():
        raise StructureError("mmCIF content is empty.", "PARSING_ERROR")
    encoded = text.encode("utf-8")
    if len(encoded) > MAX_STRUCTURE_BYTES:
        raise StructureError(
            f"Structure file exceeds {MAX_STRUCTURE_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    try:
        block = _parse_cif_data(text)
    except StructureError:
        raise
    except Exception as exc:
        raise StructureError(f"mmCIF parser failed: {exc}.", "PARSING_ERROR") from exc
    return _structure_from_cif_block(block, encoded)


def parse_pdb(text: str) -> dict:
    """Interpreta ATOM/HETATM/MODEL do formato PDB legado.

    Args:
        text: Conteudo PDB.

    Returns:
        Dict analogo ao mmCIF.

    Raises:
        StructureError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    if not isinstance(text, str) or not text.strip():
        raise StructureError("PDB content is empty.", "PARSING_ERROR")
    encoded = text.encode("utf-8")
    if len(encoded) > MAX_STRUCTURE_BYTES:
        raise StructureError(
            f"Structure file exceeds {MAX_STRUCTURE_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    atoms: List[dict] = []
    model = 1
    for raw_line in text.splitlines():
        record = raw_line[:6]
        if record.startswith("MODEL"):
            try:
                model = int(raw_line[10:14])
            except ValueError:
                model = 1
            continue
        if record.startswith("ENDMDL"):
            continue
        if not (record.startswith("ATOM") or record.startswith("HETATM")):
            continue
        if len(raw_line) < 54:
            raise StructureError("PDB ATOM line is shorter than coordinate columns.", "PARSING_ERROR")
        try:
            atom = {
                "group": "ATOM" if record.startswith("ATOM") else "HETATM",
                "atom_name": raw_line[12:16].strip(),
                "alt_id": raw_line[16:17].strip(),
                "comp_id": raw_line[17:20].strip(),
                "label_asym_id": raw_line[21:22].strip() or "A",
                "auth_asym_id": raw_line[21:22].strip() or "A",
                "label_seq_id": _optional_int(raw_line[22:26].strip()),
                "auth_seq_id": _optional_int(raw_line[22:26].strip()),
                "insertion_code": raw_line[26:27].strip(),
                "x": float(raw_line[30:38]),
                "y": float(raw_line[38:46]),
                "z": float(raw_line[46:54]),
                "occupancy": _optional_float(raw_line[54:60].strip() if len(raw_line) >= 60 else ""),
                "b_iso": _optional_float(raw_line[60:66].strip() if len(raw_line) >= 66 else ""),
                "element": raw_line[76:78].strip() if len(raw_line) >= 78 else "",
                "model": model,
                "entity_id": "1",
            }
        except (TypeError, ValueError) as exc:
            raise StructureError("PDB coordinates are not numeric.", "PARSING_ERROR") from exc
        if not scientific_checks.atom_coordinate_is_finite(atom["x"], atom["y"], atom["z"]):
            raise StructureError("PDB coordinates are not finite numbers.", "PARSING_ERROR")
        atoms.append(atom)
        if len(atoms) > MAX_ATOMS:
            raise StructureError(
                f"Structure exceeds {MAX_ATOMS:,} atoms (RESOURCE_LIMIT).",
                "RESOURCE_LIMIT",
            )
    if not atoms:
        raise StructureError("PDB file has no ATOM/HETATM records.", "PARSING_ERROR")
    return _structure_from_atoms(
        atoms,
        entry_id="",
        method="",
        resolution=None,
        polymer_sequences={},
        polymer_types={},
        content_hash=hashlib.sha256(encoded).hexdigest(),
        fmt="pdb",
    )


def parse_structure_text(text: str, *, hint: str = "") -> dict:
    """Escolhe mmCIF ou PDB pelo conteudo. Nao corrige o ficheiro.

    Args:
        text: Texto estrutural.
        hint: mmcif, pdb ou vazio.

    Returns:
        Dict parseado.

    Raises:
        StructureError.
    """
    body = str(text or "")
    kind = str(hint or "").lower()
    stripped = body.lstrip()
    if kind == "pdb" or stripped.startswith(("HEADER", "ATOM", "HETATM", "MODEL", "REMARK", "CRYST1", "TITLE")):
        if stripped.startswith("data_"):
            return parse_mmcif(body)
        return parse_pdb(body)
    if stripped.startswith("data_") or "_atom_site." in body or kind == "mmcif":
        return parse_mmcif(body)
    if "ATOM" in body[:4000] or "HETATM" in body[:4000]:
        return parse_pdb(body)
    raise StructureError("Structure text is neither mmCIF nor PDB.", "PARSING_ERROR")


def parse_rcsb_entry_metadata(payload: Mapping[str, Any]) -> dict:
    """Extrai metodo, resolucao e kind do JSON da Data API do RCSB.

    Args:
        payload: JSON de /rest/v1/core/entry/{id}.

    Returns:
        Dict structure_id, method, resolution_angstrom, kind, release_date.

    Raises:
        StructureError: PARSING_ERROR.
    """
    if not isinstance(payload, Mapping):
        raise StructureError("RCSB entry JSON is not an object.", "PARSING_ERROR")
    entry = payload.get("entry") if isinstance(payload.get("entry"), Mapping) else {}
    structure_id = str((entry or {}).get("id") or payload.get("rcsb_id") or "").upper()
    if not structure_id:
        raise StructureError("RCSB entry JSON has no entry id.", "PARSING_ERROR")
    info = payload.get("rcsb_entry_info") if isinstance(payload.get("rcsb_entry_info"), Mapping) else {}
    methodology = str(info.get("structure_determination_methodology") or "").strip().lower()
    methods: List[str] = []
    exptl = payload.get("exptl")
    if isinstance(exptl, list):
        for item in exptl:
            if isinstance(item, Mapping) and item.get("method"):
                methods.append(str(item.get("method")))
    method = methods[0] if methods else str(info.get("experimental_method") or "")
    resolution = _first_finite(
        (info.get("resolution_combined") or [None])[0]
        if isinstance(info.get("resolution_combined"), list)
        else info.get("resolution_combined")
    )
    if resolution is None:
        refine = payload.get("refine")
        if isinstance(refine, list) and refine and isinstance(refine[0], Mapping):
            resolution = _first_finite(refine[0].get("ls_d_res_high"))
    accession = payload.get("rcsb_accession_info") if isinstance(payload.get("rcsb_accession_info"), Mapping) else {}
    release = str(accession.get("initial_release_date") or "")
    if methodology in {"computational", "computational model"} or structure_id.startswith(("AF_", "MA_")):
        kind = "predicted"
    elif methodology == "experimental" or methods:
        kind = "experimental"
    else:
        kind = "predicted" if not methods else "experimental"
    return {
        "structure_id": structure_id,
        "source": "RCSB PDB",
        "kind": kind,
        "method": method,
        "resolution_angstrom": resolution,
        "release_date": release,
        "methodology": methodology,
        "deposited_atom_count": info.get("deposited_atom_count"),
        "deposited_model_count": info.get("deposited_model_count"),
    }


def parse_alphafold_prediction(payload: Any) -> dict:
    """Interpreta o JSON /api/prediction/{uniprot} da AlphaFold DB.

    Args:
        payload: Lista JSON da API.

    Returns:
        Dict predicted com URLs e metricas originais.

    Raises:
        StructureError.
    """
    rows = payload if isinstance(payload, list) else [payload]
    if not rows or not isinstance(rows[0], Mapping):
        raise StructureError("AlphaFold DB prediction JSON is empty.", "PARSING_ERROR")
    item = rows[0]
    entry_id = str(item.get("entryId") or item.get("modelEntityId") or "")
    uniprot = str(item.get("uniprotAccession") or "")
    if not entry_id or not uniprot:
        raise StructureError("AlphaFold DB JSON lacks entryId or uniprotAccession.", "PARSING_ERROR")
    cif_url = str(item.get("cifUrl") or "")
    if cif_url and not url_is_allowed(cif_url):
        raise StructureError("AlphaFold DB cifUrl is not on the allowlist.", "INVALID_INPUT")
    plddt_url = str(item.get("plddtDocUrl") or "")
    if plddt_url and not url_is_allowed(plddt_url):
        plddt_url = ""
    global_metric = _first_finite(item.get("globalMetricValue"))
    sequence = str(item.get("uniprotSequence") or item.get("sequence") or "")
    return {
        "structure_id": entry_id,
        "source": "AlphaFold DB",
        "kind": "predicted",
        "method": str(item.get("toolUsed") or "AlphaFold"),
        "uniprot": uniprot,
        "organism": str(item.get("organismScientificName") or ""),
        "latest_version": item.get("latestVersion"),
        "model_created_date": str(item.get("modelCreatedDate") or ""),
        "global_metric_value": global_metric,
        "global_metric_name": "AlphaFold DB globalMetricValue (mean pLDDT)",
        "sequence": sequence,
        "sequence_start": item.get("sequenceStart") or item.get("uniprotStart"),
        "sequence_end": item.get("sequenceEnd") or item.get("uniprotEnd"),
        "cif_url": cif_url,
        "pdb_url": str(item.get("pdbUrl") or ""),
        "plddt_doc_url": plddt_url,
        "chain_id": str(item.get("chainId") or "A"),
    }


def parse_uniprot_structure_xrefs(payload: Mapping[str, Any]) -> dict:
    """Extrai PDB e AlphaFoldDB a partir do JSON UniProtKB.

    Args:
        payload: JSON uniprotkb.

    Returns:
        Dict accession, sequence, pdb, alphafold, truncated.

    Raises:
        StructureError.
    """
    if not isinstance(payload, Mapping):
        raise StructureError("UniProt JSON is not an object.", "PARSING_ERROR")
    accession = str(payload.get("primaryAccession") or "")
    if not accession:
        raise StructureError("UniProt JSON has no primaryAccession.", "PARSING_ERROR")
    sequence_block = payload.get("sequence") if isinstance(payload.get("sequence"), Mapping) else {}
    sequence = str((sequence_block or {}).get("value") or "")
    organism = ""
    org = payload.get("organism")
    if isinstance(org, Mapping):
        organism = str(org.get("scientificName") or "")
    pdb_hits: List[dict] = []
    alphafold_hits: List[dict] = []
    xrefs = payload.get("uniProtKBCrossReferences")
    if not isinstance(xrefs, list):
        xrefs = []
    n_pdb = 0
    n_af = 0
    truncated = False
    for item in xrefs:
        if not isinstance(item, Mapping):
            continue
        database = str(item.get("database") or "")
        identifier = str(item.get("id") or "")
        properties = {
            str(prop.get("key")): str(prop.get("value"))
            for prop in (item.get("properties") or [])
            if isinstance(prop, Mapping)
        }
        if database == "PDB":
            n_pdb += 1
            if len(pdb_hits) < MAX_UNIPROT_XREFS:
                resolution = _parse_resolution_text(properties.get("Resolution") or "")
                pdb_hits.append(
                    {
                        "source": "RCSB PDB",
                        "structure_id": identifier.upper(),
                        "kind": "experimental",
                        "method": properties.get("Method") or "",
                        "resolution_angstrom": resolution,
                        "chains": properties.get("Chains") or "",
                        "match_basis": "UniProtKB PDB cross-reference",
                    }
                )
            else:
                truncated = True
        elif database == "AlphaFoldDB":
            n_af += 1
            if len(alphafold_hits) < MAX_UNIPROT_XREFS:
                alphafold_hits.append(
                    {
                        "source": "AlphaFold DB",
                        "structure_id": f"AF-{identifier}-F1" if not identifier.startswith("AF-") else identifier,
                        "kind": "predicted",
                        "method": "AlphaFold",
                        "uniprot": identifier,
                        "match_basis": "UniProtKB AlphaFoldDB cross-reference",
                    }
                )
            else:
                truncated = True
    return {
        "accession": accession,
        "organism": organism,
        "sequence": sequence,
        "pdb": pdb_hits,
        "alphafold": alphafold_hits,
        "n_pdb_reported": n_pdb,
        "n_alphafold_reported": n_af,
        "truncated": truncated,
        "source": "UniProtKB",
    }


def _labels_are_polymer_indices(residues: Sequence[Mapping[str, Any]], sequence_length: int) -> bool:
    """True se label_seq_id cabe em 1..len(polymer), o indice mmCIF usual."""
    if sequence_length <= 0:
        return False
    labels: List[int] = []
    for item in residues:
        lid = item.get("label_seq_id")
        if lid is None:
            continue
        try:
            labels.append(int(lid))
        except (TypeError, ValueError):
            return False
    if not labels:
        return False
    return all(1 <= value <= sequence_length for value in labels)


def _ordered_polymer_residues(chain: Mapping[str, Any]) -> List[Optional[dict]]:
    """Residuos na ordem do polimero, nao na numeracao auth do PDB.

    mmCIF: label_seq_id em 1..n e o indice do polimero.
    PDB legado: auth_seq_id frequentemente nao comeca em 1; nesse caso a
    ordem e a dos residuos observados, sem inventar label i+1.
    """
    sequence = str(chain.get("sequence") or "")
    residues = list(chain.get("residues") or [])
    length = len(sequence)
    if length <= 0:
        return []
    if _labels_are_polymer_indices(residues, length):
        by_label: Dict[int, dict] = {}
        for item in residues:
            lid = item.get("label_seq_id")
            if lid is None:
                continue
            key = int(lid)
            previous = by_label.get(key)
            if previous is None or (item.get("has_coordinates") and not previous.get("has_coordinates")):
                by_label[key] = item
        return [by_label.get(index + 1) for index in range(length)]
    observed = [item for item in residues if item.get("has_coordinates")]
    observed.sort(
        key=lambda item: (
            int(item.get("auth_seq_id") or 0),
            str(item.get("insertion_code") or ""),
            int(item.get("label_seq_id") or 0),
        )
    )
    if len(observed) == length:
        return list(observed)
    slots: List[Optional[dict]] = [None] * length
    cursor = 0
    for item in observed:
        letter = str(item.get("one_letter") or "")
        while cursor < length and sequence[cursor] != letter:
            cursor += 1
        if cursor >= length:
            break
        slots[cursor] = item
        cursor += 1
    return slots


def classify_mapping_status(
    *,
    sequences_identical: bool,
    polymer_index_mode: bool,
    observed_coordinate_residues: int,
    polymer_length: int,
) -> dict:
    """Classifica o mapping sem promover best-effort a exact.

    Args:
        sequences_identical: Query igual a sequencia do polimero.
        polymer_index_mode: label_seq_id em 1..n do polimero.
        observed_coordinate_residues: Residuos com coordenadas ATOM.
        polymer_length: Comprimento da sequencia do polimero.

    Returns:
        Dict mapping_status e mapping_status_note.

    Raises:
        Nenhum.

    Nota biologica:
        Sem indices de polimero completos (_entity_poly_seq / label_seq_id
        1..n), o mapping segue a ordem dos ATOM e nao e um mapeamento exacto
        da sequencia. Residuos sem coordenadas nao alteram essa classificacao.
    """
    if polymer_length <= 0:
        return {
            "mapping_status": "UNCERTAIN",
            "mapping_status_note": (
                "Polymer length is zero or unknown. Mapping status cannot be proven."
            ),
        }
    if not sequences_identical:
        return {
            "mapping_status": "ALIGNED",
            "mapping_status_note": (
                "Sequences differ. Mapping used Needleman-Wunsch (BLOSUM62). "
                "This is not exact sequence mapping."
            ),
        }
    if polymer_index_mode:
        return {
            "mapping_status": "EXACT",
            "mapping_status_note": (
                "Identical sequences and polymer indices from label_seq_id in "
                "1..n (mmCIF polymer order). Deposited auth_seq_id is kept separate."
            ),
        }
    if observed_coordinate_residues == polymer_length:
        return {
            "mapping_status": "PARTIAL",
            "mapping_status_note": (
                "Identical letter sequence, but residue identifiers are not "
                "polymer indices (complete _entity_poly_seq / label_seq_id 1..n "
                "unavailable). Mapping follows observed ATOM residue order. "
                "This is not exact sequence mapping."
            ),
        }
    return {
        "mapping_status": "BEST_EFFORT",
        "mapping_status_note": (
            "Identical letter sequence, but the polymer is longer than observed "
            "ATOM residues and complete polymer indices are unavailable. Mapping "
            "is a letter walk along ATOM residues. This is not exact sequence "
            "mapping."
        ),
    }


def map_query_to_chain(
    query_sequence: str,
    chain: Mapping[str, Any],
    *,
    molecule: str = "PROTEIN",
) -> dict:
    """Mapeia a sequencia de analise a uma chain da estrutura.

    Usa identidade exata quando as sequencias coincidem; caso contrario
    Needleman-Wunsch global (BLOSUM62 para proteina; identidade para DNA/RNA).
    PDB residue number nao e assumido igual ao indice da sequencia.

    Args:
        query_sequence: Polimero de analise.
        chain: Chain parseada.
        molecule: PROTEIN, DNA ou RNA.

    Returns:
        Dict mapping, identity, coverage, method.

    Raises:
        StructureError: INVALID_INPUT.
    """
    kind = str(molecule or "PROTEIN").strip().upper()
    if kind not in {"PROTEIN", "DNA", "RNA"}:
        raise StructureError("Mapping molecule must be PROTEIN, DNA or RNA.", "INVALID_INPUT")
    query_info = dna_analysis.validate_for_molecule(query_sequence, kind)
    if not query_info["is_valid"]:
        raise StructureError(
            str(query_info.get("rejection_reason") or f"Query is not a standard {kind}."),
            "INVALID_INPUT",
        )
    query = str(query_info["sequence"])
    chain_seq = str(chain.get("sequence") or "")
    if not chain_seq:
        raise StructureError("Selected chain has no polymer sequence.", "PARSING_ERROR")
    polymer = _ordered_polymer_residues(chain)
    if query == chain_seq:
        aligned_q, aligned_c = query, chain_seq
        method = (
            "Identical sequences. Query index i (0-based) maps to polymer index i "
            "on the selected chain. Deposited label_seq_id and auth_seq_id come "
            "from the file and are not assumed equal to i+1."
        )
        identity = 100.0
    else:
        aligned = alignment.pairwise_global(query, chain_seq)
        aligned_q = str(aligned.get("aligned_seq1") or "")
        aligned_c = str(aligned.get("aligned_seq2") or "")
        identity = aligned.get("ungapped_identity_pct")
        if kind == "PROTEIN":
            method = (
                "Needleman-Wunsch global alignment (BLOSUM62) between the analysis "
                "sequence and the structure polymer sequence. Residue numbers are "
                "taken from the matched polymer residue, not from i+1."
            )
        else:
            method = (
                "Needleman-Wunsch global alignment (nucleotide identity) between "
                "the analysis sequence and the deposited polymer. Residue numbers "
                "come from the file, not from i+1."
            )
    mapping: List[dict] = []
    q_index = 0
    c_index = 0
    for aa_q, aa_c in zip(aligned_q, aligned_c):
        query_pos = None
        chain_pos = None
        if aa_q != "-":
            query_pos = q_index
            q_index += 1
        if aa_c != "-":
            chain_pos = c_index
            c_index += 1
        if query_pos is None:
            continue
        residue = None
        if chain_pos is not None and 0 <= chain_pos < len(polymer):
            residue = polymer[chain_pos]
        has_coords = bool(residue and residue.get("has_coordinates"))
        row_status = "UNMAPPED" if residue is None or chain_pos is None else None
        mapping.append(
            {
                "query_index_0based": query_pos,
                "query_residue": query[query_pos],
                "chain_index_0based": chain_pos,
                "chain_residue": chain_seq[chain_pos] if chain_pos is not None else None,
                "aligned_query": aa_q,
                "aligned_chain": aa_c,
                "match": aa_q != "-" and aa_c != "-" and aa_q == aa_c,
                "chain_id": chain.get("chain_id"),
                "label_seq_id": None if residue is None else residue.get("label_seq_id"),
                "auth_seq_id": None if residue is None else residue.get("auth_seq_id"),
                "insertion_code": "" if residue is None else residue.get("insertion_code") or "",
                "comp_id": None if residue is None else residue.get("comp_id"),
                "has_coordinates": has_coords,
                "coordinate_status": "available" if has_coords else "unavailable",
                "sequence_residue_exists": True,
                "mapping_row_status": row_status,
                "ca": None if residue is None else residue.get("ca"),
                "p": None if residue is None else residue.get("p"),
                "c1": None if residue is None else residue.get("c1"),
                "ca_b_iso": None if residue is None else residue.get("ca_b_iso"),
                "model": None if residue is None else residue.get("model"),
                "alt_ids": [] if residue is None else list(residue.get("alt_ids") or []),
            }
        )
    residues = list(chain.get("residues") or [])
    polymer_index_mode = _labels_are_polymer_indices(residues, len(chain_seq))
    observed_n = sum(1 for item in residues if item.get("has_coordinates"))
    n_letter_matches = sum(1 for item in mapping if item.get("match"))
    mapped_query_indices = [
        int(item["query_index_0based"])
        for item in mapping
        if item.get("chain_index_0based") is not None
    ]
    unmapped_query_indices = [
        int(item["query_index_0based"])
        for item in mapping
        if item.get("chain_index_0based") is None
    ]
    classified = classify_mapping_status(
        sequences_identical=query == chain_seq,
        polymer_index_mode=polymer_index_mode,
        observed_coordinate_residues=observed_n,
        polymer_length=len(chain_seq),
    )
    if n_letter_matches == 0:
        classified = {
            "mapping_status": "UNMAPPED",
            "mapping_status_note": (
                "No identical aligned letters between the analysis sequence and "
                "the deposited polymer. HelixScope does not promote a forced "
                "best alignment as a valid mapping."
            ),
        }
    elif unmapped_query_indices and query != chain_seq:
        classified = {
            "mapping_status": "PARTIAL",
            "mapping_status_note": (
                "Only a subset of analysis positions align to polymer residues. "
                "mapped_query_indices and unmapped_query_indices list the extent."
            ),
        }
    for item in mapping:
        if item.get("mapping_row_status") is None:
            item["mapping_row_status"] = classified["mapping_status"]
        if classified["mapping_status"] == "UNMAPPED":
            item["mapping_row_status"] = "UNMAPPED"
    n_with_coords = sum(1 for item in mapping if item["has_coordinates"])
    coverage = coverage_from_mapping(len(query), n_with_coords)
    return {
        "mapping": mapping,
        "identity_ungapped_pct": identity,
        "identity_definition": (
            "Ungapped identity percent from identical sequences or from "
            "Needleman-Wunsch (BLOSUM62) between query and polymer sequence."
        ),
        "coverage_query_with_coordinates": coverage,
        "coverage_scale": "fraction_0_1",
        "coverage_definition": "n_query_residues_with_coordinates / len(analysis_sequence)",
        "n_query": len(query),
        "n_with_coordinates": n_with_coords,
        "n_without_coordinates": len(query) - n_with_coords,
        "n_letter_matches": n_letter_matches,
        "mapped_query_indices": mapped_query_indices,
        "unmapped_query_indices": unmapped_query_indices,
        "chain_length": len(chain_seq),
        "mapping_method": method,
        "mapping_status": classified["mapping_status"],
        "mapping_status_note": classified["mapping_status_note"],
        "polymer_index_mode": polymer_index_mode,
        "query_hash": provenance.sequence_digest(query),
        "chain_hash": provenance.sequence_digest(chain_seq),
        "sequences_identical": query == chain_seq,
        "coordinate_system_query": "0-based",
        "coordinate_system_structure": (
            "Polymer index 0-based; deposited label_seq_id and auth_seq_id as in the file"
        ),
    }


def coverage_from_mapping(query_length: int, n_with_coordinates: int) -> Optional[float]:
    """coverage = residuos da query com coordenadas / comprimento da query.

    Args:
        query_length: Comprimento da proteina de analise.
        n_with_coordinates: Quantos residuos da query tem atomos.

    Returns:
        Fracao em [0, 1], ou None se o comprimento for zero.

    Raises:
        Nenhum.
    """
    if query_length <= 0:
        return None
    return n_with_coordinates / float(query_length)


def _mapping_rows(result: Mapping[str, Any]) -> List[dict]:
    """Linhas de mapping do envelope cientifico ou do structure_3d_input."""
    rows = result.get("residue_mapping")
    if rows is None:
        rows = result.get("mappings")
    return list(rows or [])


def query_position_to_residue(result: Mapping[str, Any], position: int) -> Optional[dict]:
    """Devolve o mapeamento de uma posicao 0-based da query.

    Args:
        result: Envelope com residue_mapping ou mappings.
        position: Indice 0-based.

    Returns:
        Item de mapping ou None.

    Raises:
        StructureError: INVALID_INPUT se fora do intervalo.
    """
    sequence = str(result.get("sequence") or "")
    rows = _mapping_rows(result)
    if sequence:
        if position < 0 or position >= len(sequence):
            raise StructureError("Position is outside the analysis protein.", "INVALID_INPUT")
    elif rows:
        max_index = max(int(item.get("query_index_0based") or 0) for item in rows)
        if position < 0 or position > max_index:
            raise StructureError("Position is outside the analysis protein.", "INVALID_INPUT")
    else:
        raise StructureError("Position is outside the analysis protein.", "INVALID_INPUT")
    for item in rows:
        if int(item.get("query_index_0based")) == position:
            return dict(item)
    return None


def structure_residue_to_query(
    result: Mapping[str, Any],
    *,
    chain_id: str,
    label_seq_id: Optional[int] = None,
    auth_seq_id: Optional[int] = None,
    insertion_code: str = "",
    model: Optional[int] = None,
) -> Optional[dict]:
    """Devolve a posicao da query para um residuo estrutural. Nao adivinha.

    Args:
        result: Envelope cientifico ou structure_3d_input.
        chain_id: Chain depositada.
        label_seq_id: label_seq_id mmCIF, se conhecido.
        auth_seq_id: auth_seq_id PDB, se conhecido.
        insertion_code: Codigo de insercao depositado.
        model: Modelo NMR, se relevante.

    Returns:
        Linha de mapping ou None se nao houver correspondencia.

    Raises:
        StructureError: INVALID_INPUT se chain_id estiver vazio.

    Nota biologica:
        label_seq_id e auth_seq_id nao sao o mesmo numero. A funcao nao
        converte um no outro e nao aplica +1/-1.
    """
    if not str(chain_id or "").strip():
        raise StructureError("Chain ID is required for reverse mapping.", "INVALID_INPUT")
    ins = str(insertion_code or "")
    matches: List[dict] = []
    for item in _mapping_rows(result):
        if str(item.get("chain_id") or "") != str(chain_id):
            continue
        if model is not None and item.get("model") is not None and int(item.get("model")) != int(model):
            continue
        if str(item.get("insertion_code") or "") != ins:
            continue
        if label_seq_id is not None:
            if item.get("label_seq_id") is None or int(item.get("label_seq_id")) != int(label_seq_id):
                continue
        if auth_seq_id is not None:
            if item.get("auth_seq_id") is None or int(item.get("auth_seq_id")) != int(auth_seq_id):
                continue
        if label_seq_id is None and auth_seq_id is None:
            continue
        matches.append(dict(item))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        return None
    unique_query = {item.get("query_index_0based") for item in matches}
    if len(unique_query) == 1:
        return matches[0]
    return None


def map_msa_column_to_structure_residue(
    coordinate_map_row: Sequence[Optional[int]],
    column: int,
    result: Mapping[str, Any],
) -> Optional[dict]:
    """Coluna de MSA -> posicao original -> residuo estrutural. Nao assume igualdade.

    Args:
        coordinate_map_row: maps[seq] do MSA.
        column: Coluna 0-based.
        result: Envelope estrutural.

    Returns:
        Mapping do residuo ou None se gap.

    Raises:
        StructureError: INVALID_INPUT.
    """
    if column < 0 or column >= len(coordinate_map_row):
        raise StructureError("MSA column is outside the coordinate map.", "INVALID_INPUT")
    query_pos = coordinate_map_row[column]
    if query_pos is None:
        return None
    return query_position_to_residue(result, int(query_pos))


def residues_for_sequence_span(
    result: Mapping[str, Any], start: int, end: int
) -> List[dict]:
    """Residuos estruturais de um intervalo da sequencia de analise [start, end).

    Args:
        result: Envelope.
        start: 0-based inclusivo.
        end: 0-based exclusivo.

    Returns:
        Itens de mapping no intervalo.

    Raises:
        StructureError.
    """
    sequence = str(result.get("sequence") or "")
    if start < 0 or end > len(sequence) or start >= end:
        raise StructureError("Sequence span is invalid.", "INVALID_INPUT")
    rows = []
    for item in list(result.get("residue_mapping") or result.get("mappings") or []):
        index = int(item.get("query_index_0based"))
        if start <= index < end:
            rows.append(dict(item))
    return rows


def read_bundled_mmcif(filename: str) -> str:
    """Le uma copia mmCIF depositada em tests/fixtures. Sem fetch.

    Args:
        filename: Nome do ficheiro, so basename .cif.

    Returns:
        Texto mmCIF.

    Raises:
        StructureError: INVALID_INPUT ou NOT_FOUND.
    """
    name = Path(str(filename or "")).name
    if name != str(filename or "").strip() or not re.match(r"^[0-9A-Za-z][0-9A-Za-z_.-]*\.cif$", name):
        raise StructureError("Bundled fixture name is invalid.", "INVALID_INPUT")
    path = BUNDLED_MMCIF_DIR / name
    if not path.is_file():
        raise StructureError("Bundled deposited mmCIF is not present.", "NOT_FOUND")
    return path.read_text(encoding="utf-8")


def bundled_deposited_hit(structure_id: str, query: str, molecule: str) -> Optional[dict]:
    """Hit de busca a partir de mmCIF local se a sequencia do polimero coincidir.

    Args:
        structure_id: PDB ID.
        query: Sequencia de analise ja validada.
        molecule: PROTEIN, DNA ou RNA.

    Returns:
        Hit experimental ou None. None nao e estrutura.

    Raises:
        Nenhum. Parse falhado devolve None.

    Nota biologica:
        A copia em tests/fixtures e o ficheiro RCSB, nao geometria procedural.
        Nao substitui um fetch live; o source declara bundled deposited mmCIF copy.
    """
    pdb_id = str(structure_id or "").strip().upper()
    mol = str(molecule or "PROTEIN").strip().upper()
    if mol not in {"PROTEIN", "DNA", "RNA"} or not pdb_id:
        return None
    cif_path = BUNDLED_MMCIF_DIR / f"{pdb_id}.cif"
    meta_path = BUNDLED_MMCIF_DIR / f"rcsb_entry_{pdb_id}.json"
    if not cif_path.is_file() or not meta_path.is_file():
        return None
    needle = {
        "PROTEIN": "polypeptide",
        "DNA": "polydeoxyribonucleotide",
        "RNA": "polyribonucleotide",
    }[mol]
    query_u = str(query or "").strip().upper()
    try:
        parsed = parse_mmcif(cif_path.read_text(encoding="utf-8"))
        raw_meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (StructureError, OSError, json.JSONDecodeError, TypeError, ValueError):
        return None
    matched = False
    for chain in list(parsed.get("chains") or []):
        ptype = str(chain.get("polymer_type") or "").lower()
        seq = str(chain.get("sequence") or "").upper()
        if needle in ptype and seq == query_u:
            matched = True
            break
    if not matched:
        return None
    if isinstance(raw_meta, Mapping) and (
        "rcsb_entry_info" in raw_meta or isinstance(raw_meta.get("entry"), Mapping)
    ):
        try:
            meta = parse_rcsb_entry_metadata(raw_meta)
        except StructureError:
            meta = {}
    else:
        meta = dict(raw_meta) if isinstance(raw_meta, Mapping) else {}
    return {
        **meta,
        "structure_id": pdb_id,
        "source": "RCSB PDB",
        "kind": "experimental",
        "match_basis": (
            f"Bundled deposited mmCIF copy ({cif_path.name}); polymer sequence "
            "identical to the analysis sequence. Not a live RCSB fetch."
        ),
        "bundled_fixture": cif_path.name,
        "cache_status": "bundled",
    }


def search_structure_sources(
    query_sequence: str,
    *,
    pdb_id: str = "",
    uniprot: str = "",
    alphafold_id: str = "",
    search_pdb_by_sequence: bool = False,
    lookup_uniprot: bool = False,
    lookup_alphafold: bool = False,
    identity_cutoff: float = DEFAULT_IDENTITY_CUTOFF,
    molecule: str = "PROTEIN",
    urlopen_fn=None,
) -> dict:
    """Procura fontes reais. Nao escolhe uma estrutura 'melhor'.

    Args:
        query_sequence: Polimero de analise.
        pdb_id: PDB ID opcional.
        uniprot: Accession UniProt opcional.
        alphafold_id: AF-... ou UniProt para AlphaFold DB.
        search_pdb_by_sequence: Se True, RCSB Search API (MMseqs2).
        lookup_uniprot: Se True, cruzamentos UniProt (proteina).
        lookup_alphafold: Se True, metadados AlphaFold DB (proteina).
        identity_cutoff: Corte da Search API, 0-1, nao um score caseiro.
        molecule: PROTEIN, DNA ou RNA.
        urlopen_fn: Cliente HTTP injetavel.

    Returns:
        Dict status, hits, provenance. Hits vazios sao NO_STRUCTURE, nao fake.

    Raises:
        StructureError.
    """
    molecule_kind = str(molecule or "PROTEIN").strip().upper()
    if molecule_kind not in {"PROTEIN", "DNA", "RNA"}:
        raise StructureError("Structure search molecule must be PROTEIN, DNA or RNA.", "INVALID_INPUT")
    query_info = dna_analysis.validate_for_molecule(query_sequence, molecule_kind)
    if not query_info["is_valid"]:
        raise StructureError(
            str(query_info.get("rejection_reason") or f"Query is not a standard {molecule_kind}."),
            "INVALID_INPUT",
        )
    query = str(query_info["sequence"])
    hits: List[dict] = []
    errors: List[dict] = []
    parameters = {
        "identity_cutoff": identity_cutoff,
        "evalue_cutoff": DEFAULT_EVALUE_CUTOFF,
        "search_pdb_by_sequence": bool(search_pdb_by_sequence),
        "lookup_uniprot": bool(lookup_uniprot),
        "lookup_alphafold": bool(lookup_alphafold),
        "molecule": molecule_kind,
        "sequence_search_engine": "RCSB Search API sequence service (MMseqs2)",
    }
    if pdb_id:
        classified = classify_structure_identifier(pdb_id)
        if classified["kind"] != "pdb":
            raise StructureError("PDB identifier is invalid.", "INVALID_INPUT")
        bundled = bundled_deposited_hit(classified["canonical"], query, molecule_kind)
        if bundled:
            hits.append(bundled)
        try:
            meta = fetch_rcsb_entry(classified["canonical"], urlopen_fn=urlopen_fn)
            hits.append({**meta, "match_basis": "PDB ID supplied by the user"})
        except StructureError as exc:
            errors.append({"source": "RCSB PDB", "category": exc.category, "message": str(exc)})
    if search_pdb_by_sequence:
        try:
            hits.extend(
                search_rcsb_by_sequence(
                    query,
                    identity_cutoff=identity_cutoff,
                    sequence_type={"PROTEIN": "protein", "DNA": "dna", "RNA": "rna"}[molecule_kind],
                    urlopen_fn=urlopen_fn,
                )
            )
        except StructureError as exc:
            errors.append({"source": "RCSB PDB Search", "category": exc.category, "message": str(exc)})
    if molecule_kind != "PROTEIN" and (lookup_uniprot or lookup_alphafold):
        errors.append(
            {
                "source": "AlphaFold/UniProt",
                "category": "UNAVAILABLE",
                "message": "AlphaFold DB and UniProt structure xrefs are protein resources.",
            }
        )
    uniprot_acc = ""
    af_uniprot = ""
    if molecule_kind == "PROTEIN":
        if uniprot:
            classified = classify_structure_identifier(uniprot)
            if classified["kind"] not in {"uniprot", "alphafold"}:
                raise StructureError("UniProt accession is invalid.", "INVALID_INPUT")
            uniprot_acc = classified["uniprot"]
        if alphafold_id:
            classified = classify_structure_identifier(alphafold_id)
            if classified["kind"] == "alphafold":
                af_uniprot = classified["uniprot"]
            elif classified["kind"] == "uniprot":
                af_uniprot = classified["uniprot"]
            else:
                raise StructureError("AlphaFold identifier is invalid.", "INVALID_INPUT")
        if lookup_uniprot:
            accession = uniprot_acc or af_uniprot
            if not accession:
                raise StructureError("UniProt lookup requires a UniProt accession.", "INVALID_INPUT")
            try:
                xref = fetch_uniprot_entry(accession, urlopen_fn=urlopen_fn)
                hits.extend(xref["pdb"])
                hits.extend(xref["alphafold"])
            except StructureError as exc:
                errors.append({"source": "UniProtKB", "category": exc.category, "message": str(exc)})
        if lookup_alphafold:
            accession = af_uniprot or uniprot_acc
            if not accession:
                raise StructureError("AlphaFold lookup requires a UniProt accession.", "INVALID_INPUT")
            try:
                prediction = fetch_alphafold_prediction(accession, urlopen_fn=urlopen_fn)
                hits.append(
                    {
                        "source": prediction["source"],
                        "structure_id": prediction["structure_id"],
                        "kind": "predicted",
                        "method": prediction["method"],
                        "uniprot": prediction["uniprot"],
                        "global_metric_value": prediction["global_metric_value"],
                        "global_metric_name": prediction["global_metric_name"],
                        "latest_version": prediction["latest_version"],
                        "cif_url": prediction["cif_url"],
                        "match_basis": "AlphaFold DB prediction API",
                    }
                )
            except StructureError as exc:
                errors.append({"source": "AlphaFold DB", "category": exc.category, "message": str(exc)})
    hits = _dedupe_hits(hits)[:MAX_SEARCH_HITS]
    if hits:
        status = "RETRIEVED"
        category = ""
        message = "Structure records were retrieved. Select one; HelixScope does not pick a best structure."
    elif errors and not hits:
        first = errors[0]
        raise StructureError(
            f"Structure retrieval failed. {first['source']}: {first['message']}",
            str(first["category"]),
        )
    else:
        status = "UNAVAILABLE"
        category = "NO_STRUCTURE"
        message = "No validated structure found."
    envelope = provenance.analysis_envelope(
        module="PROTEIN structure search",
        payload={"n_hits": len(hits), "n_errors": len(errors)},
        status=status,
        algorithm="RCSB/AlphaFold/UniProt lookup",
        parameters=parameters,
        sequence=query,
        source="structure source search",
    )
    envelope.update(
        {
            "status": status,
            "category": category,
            "message": message,
            "hits": hits,
            "errors": errors,
            "sequence": query,
            "sequence_hash": provenance.sequence_digest(query),
            "disclaimer": (
                "A database hit is not sequence identity. Experimental and predicted "
                "records are listed separately. Structure is not function."
            ),
        }
    )
    return envelope


def search_rcsb_by_sequence(
    sequence: str,
    *,
    identity_cutoff: float = DEFAULT_IDENTITY_CUTOFF,
    sequence_type: str = "protein",
    urlopen_fn=None,
) -> List[dict]:
    """Busca RCSB Search API (MMseqs2). Nao calcula um score estrutural proprio.

    Args:
        sequence: Polimero da busca.
        identity_cutoff: Parametro da API, 0-1.
        sequence_type: protein, dna ou rna (parametro oficial da Search API).
        urlopen_fn: HTTP injetavel.

    Returns:
        Lista de hits polymer_entity.

    Raises:
        StructureError.
    """
    kind = str(sequence_type or "protein").strip().lower()
    if kind not in {"protein", "dna", "rna"}:
        raise StructureError(
            "RCSB sequence_type must be protein, dna or rna.",
            "INVALID_INPUT",
        )
    if len(sequence) < MIN_SEARCH_SEQUENCE:
        raise StructureError(
            f"RCSB sequence search requires at least {MIN_SEARCH_SEQUENCE} residues.",
            "INVALID_INPUT",
        )
    if len(sequence) > MAX_SEARCH_SEQUENCE:
        raise StructureError(
            f"RCSB sequence search is limited to {MAX_SEARCH_SEQUENCE} residues.",
            "RESOURCE_LIMIT",
        )
    try:
        cutoff = float(identity_cutoff)
    except (TypeError, ValueError) as exc:
        raise StructureError("identity_cutoff must be a number.", "INVALID_INPUT") from exc
    if cutoff < 0.0 or cutoff > 1.0:
        raise StructureError("identity_cutoff must be in [0, 1] (RCSB Search API units).", "INVALID_INPUT")
    body = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": DEFAULT_EVALUE_CUTOFF,
                "identity_cutoff": cutoff,
                "sequence_type": kind,
                "value": sequence,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": MAX_SEARCH_HITS},
            "results_content_type": ["experimental"],
            "scoring_strategy": "sequence",
        },
    }
    payload = _http_json(
        RCSB_SEARCH,
        method="POST",
        json_body=body,
        urlopen_fn=urlopen_fn,
        empty_ok=True,
    )
    if payload in (None, "", {}):
        return []
    if not isinstance(payload, Mapping):
        raise StructureError("RCSB Search API JSON is malformed.", "PARSING_ERROR")
    result_set = payload.get("result_set")
    if not isinstance(result_set, list):
        return []
    hits: List[dict] = []
    for item in result_set:
        if not isinstance(item, Mapping):
            continue
        identifier = str(item.get("identifier") or "")
        pdb_id, entity_id = _split_polymer_entity(identifier)
        if not pdb_id:
            continue
        score = _first_finite(item.get("score"))
        hits.append(
            {
                "source": "RCSB PDB",
                "structure_id": pdb_id,
                "entity_id": entity_id,
                "kind": "experimental",
                "rcsb_score": score,
                "identity_cutoff_requested": cutoff,
                "match_basis": "RCSB Search API sequence service (MMseqs2)",
                "sequence_type": kind,
            }
        )
    return hits


def fetch_rcsb_entry(pdb_id: str, *, urlopen_fn=None) -> dict:
    """GET Data API /rest/v1/core/entry/{id}.

    Args:
        pdb_id: PDB ID de 4 caracteres.
        urlopen_fn: HTTP injetavel.

    Returns:
        Metadados parseados.

    Raises:
        StructureError.
    """
    classified = classify_structure_identifier(pdb_id)
    if classified["kind"] != "pdb":
        raise StructureError("PDB identifier is invalid.", "INVALID_INPUT")
    url = f"{RCSB_DATA}/rest/v1/core/entry/{classified['canonical']}"
    payload = _http_json(url, method="GET", urlopen_fn=urlopen_fn)
    if not isinstance(payload, Mapping):
        raise StructureError("RCSB entry response is not JSON object.", "PARSING_ERROR")
    return parse_rcsb_entry_metadata(payload)


def fetch_alphafold_prediction(uniprot: str, *, urlopen_fn=None) -> dict:
    """GET AlphaFold DB /api/prediction/{accession}.

    Args:
        uniprot: Accession UniProt.
        urlopen_fn: HTTP injetavel.

    Returns:
        Metadados predicted.

    Raises:
        StructureError.
    """
    classified = classify_structure_identifier(uniprot)
    accession = classified["uniprot"]
    if not accession:
        raise StructureError("AlphaFold lookup requires a UniProt accession.", "INVALID_INPUT")
    url = f"{ALPHAFOLD_API}/api/prediction/{accession}"
    payload = _http_json(url, method="GET", urlopen_fn=urlopen_fn)
    return parse_alphafold_prediction(payload)


def fetch_uniprot_entry(uniprot: str, *, urlopen_fn=None) -> dict:
    """GET UniProtKB JSON (accession + sequence + xrefs estruturais).

    Args:
        uniprot: Accession.
        urlopen_fn: HTTP injetavel.

    Returns:
        parse_uniprot_structure_xrefs.

    Raises:
        StructureError.
    """
    classified = classify_structure_identifier(uniprot)
    if not classified["uniprot"]:
        raise StructureError("UniProt accession is invalid.", "INVALID_INPUT")
    url = (
        f"{UNIPROT_REST}/uniprotkb/{classified['uniprot']}.json"
        "?fields=accession,id,organism_name,sequence,xref_pdb,xref_alphafolddb"
    )
    payload = _http_json(url, method="GET", urlopen_fn=urlopen_fn)
    if not isinstance(payload, Mapping):
        raise StructureError("UniProt response is not JSON object.", "PARSING_ERROR")
    return parse_uniprot_structure_xrefs(payload)


def structure_cache_root() -> Path:
    """Project-local cache for authenticated RCSB/AlphaFold file retrievals.

    Args:
        Nenhum.

    Returns:
        Directory path. The directory may not exist yet.

    Raises:
        Nenhum.
    """
    env = (os.environ.get(STRUCTURE_CACHE_ENV) or "").strip().strip('"')
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[1] / "data" / "structure_cache"


def structure_file_cache_identity(url: str) -> str:
    """SHA-256 of an allowlisted structure file URL.

    Args:
        url: Absolute HTTPS URL already accepted by url_is_allowed.

    Returns:
        Hex digest.

    Raises:
        Nenhum.
    """
    return hashlib.sha256(str(url or "").encode("utf-8")).hexdigest()


def read_cached_structure_file(url: str) -> Optional[tuple[str, dict]]:
    """Return cached mmCIF/PDB text if the sidecar identity still matches.

    Args:
        url: Allowlisted download URL.

    Returns:
        (text, provenance) or None. Never invents coordinates.

    Raises:
        Nenhum.
    """
    if not url_is_allowed(url):
        return None
    identity = structure_file_cache_identity(url)
    root = structure_cache_root()
    meta_path = root / identity[:2] / f"{identity}.json"
    body_path = root / identity[:2] / f"{identity}.cif"
    if not meta_path.is_file() or not body_path.is_file():
        return None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        raw = body_path.read_bytes()
    except (OSError, json.JSONDecodeError, UnicodeError):
        return None
    if not isinstance(meta, dict):
        return None
    expected = str(meta.get("sha256") or "")
    actual = hashlib.sha256(raw).hexdigest()
    if expected != actual or str(meta.get("url") or "") != url:
        return None
    try:
        text = raw.decode("utf-8")
    except UnicodeError:
        return None
    meta = dict(meta)
    meta["cache_hit"] = True
    return text, meta


def write_cached_structure_file(url: str, text: str) -> dict:
    """Atomically store a retrieved structure file with identity provenance.

    Args:
        url: Allowlisted source URL.
        text: Validated decoded file body.

    Returns:
        Provenance dict written beside the file.

    Raises:
        Nenhum. Cache failures are silent; retrieval already succeeded.
    """
    if not url_is_allowed(url):
        return {}
    encoded = str(text or "").encode("utf-8")
    identity = structure_file_cache_identity(url)
    parsed = urlparse(url)
    record = {
        "url": url,
        "host": (parsed.hostname or "").lower(),
        "path": parsed.path,
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "n_bytes": len(encoded),
        "retrieved_at": provenance.utc_now(),
        "source": "authoritative HTTPS retrieval",
    }
    root = structure_cache_root() / identity[:2]
    body_path = root / f"{identity}.cif"
    meta_path = root / f"{identity}.json"
    tmp_body = root / f"{identity}.cif.partial"
    tmp_meta = root / f"{identity}.json.partial"
    try:
        root.mkdir(parents=True, exist_ok=True)
        tmp_body.write_bytes(encoded)
        tmp_meta.write_text(json.dumps(record, indent=2), encoding="utf-8")
        os.replace(tmp_body, body_path)
        os.replace(tmp_meta, meta_path)
    except OSError:
        try:
            tmp_body.unlink(missing_ok=True)
            tmp_meta.unlink(missing_ok=True)
        except OSError:
            pass
        return {}
    return record


def fetch_structure_file(url: str, *, urlopen_fn=None) -> str:
    """Descarrega mmCIF/PDB de um URL allowlisted.

    Args:
        url: HTTPS na allowlist.
        urlopen_fn: HTTP injetavel.

    Returns:
        Texto UTF-8. Pode ser um cache local da mesma URL com SHA-256.

    Raises:
        StructureError.
    """
    if urlopen_fn is None:
        cached = read_cached_structure_file(url)
        if cached is not None:
            return cached[0]
    text = _http_text(
        url,
        method="GET",
        urlopen_fn=urlopen_fn,
        limit=MAX_STRUCTURE_BYTES,
        timeout_s=STRUCTURE_FILE_TIMEOUT_S,
    )
    if urlopen_fn is None:
        write_cached_structure_file(url, text)
    return text


def load_protein_structure(
    query_sequence: str,
    *,
    source: str,
    structure_id: str,
    chain_id: str = "",
    model_number: Optional[int] = None,
    cif_url: str = "",
    urlopen_fn=None,
    structure_text: str = "",
    metadata: Optional[Mapping[str, Any]] = None,
    cache: Optional[dict] = None,
    molecule: str = "PROTEIN",
) -> dict:
    """Recupera, valida, mapeia e empacota uma estrutura. Sem coordenadas inventadas.

    Args:
        query_sequence: Polimero de analise.
        source: RCSB PDB ou AlphaFold DB.
        structure_id: Identificador.
        chain_id: Chain a mapear; vazio usa a primeira do tipo pedido.
        model_number: Modelo NMR; None lista todos e mapeia o primeiro declarado.
        cif_url: URL allowlisted, tipicamente da AlphaFold DB.
        urlopen_fn: HTTP injetavel.
        structure_text: mmCIF/PDB ja obtido (testes).
        metadata: Metadados ja parseados.
        cache: Dict injetavel.
        molecule: PROTEIN, DNA ou RNA.

    Returns:
        ProteinStructureResult (tambem usado para DNA/RNA depositados).

    Raises:
        StructureError.
    """
    kind_mol = str(molecule or "PROTEIN").strip().upper()
    if kind_mol not in {"PROTEIN", "DNA", "RNA"}:
        raise StructureError("Structure molecule must be PROTEIN, DNA or RNA.", "INVALID_INPUT")
    if kind_mol != "PROTEIN" and str(source or "") == "AlphaFold DB":
        raise StructureError(
            "AlphaFold DB models are protein predictions, not DNA/RNA 3D.",
            "INVALID_INPUT",
        )
    query_info = dna_analysis.validate_for_molecule(query_sequence, kind_mol)
    if not query_info["is_valid"]:
        raise StructureError(
            str(query_info.get("rejection_reason") or f"Query is not a standard {kind_mol}."),
            "INVALID_INPUT",
        )
    query = str(query_info["sequence"])
    parameters = {
        "chain_id": chain_id,
        "model_number": model_number,
        "file_preference": "mmCIF",
        "molecule": kind_mol,
    }
    key = cache_key(
        source=source,
        identifier=structure_id,
        chain=chain_id,
        parameters={**parameters, "query_hash": provenance.sequence_digest(query)},
    )
    store = cache if cache is not None else _CACHE
    cached = store.get(key)
    if isinstance(cached, dict):
        return mark_cached_result(cached)
    meta = dict(metadata or {})
    src = str(source or "")
    if src == "RCSB PDB" and not meta:
        meta = fetch_rcsb_entry(structure_id, urlopen_fn=urlopen_fn)
    if src == "AlphaFold DB":
        needs_af = not str(cif_url or "").strip() and not str(meta.get("cif_url") or "").strip()
        if not meta or (needs_af and not structure_text):
            classified = classify_structure_identifier(structure_id)
            fetched = fetch_alphafold_prediction(
                classified["uniprot"] or str(meta.get("uniprot") or structure_id),
                urlopen_fn=urlopen_fn,
            )
            meta = {**meta, **fetched}
    text = structure_text
    if not text:
        url = cif_url or str(meta.get("cif_url") or "")
        if not url and src == "RCSB PDB":
            url = f"{RCSB_FILES}/download/{str(structure_id).upper()}.cif"
        if not url:
            raise StructureError("No allowlisted structure file URL is available.", "INVALID_INPUT")
        text = fetch_structure_file(url, urlopen_fn=urlopen_fn)
    parsed = parse_structure_text(text, hint="mmcif")
    kind = str(meta.get("kind") or parsed.get("kind") or "unavailable")
    if src == "AlphaFold DB":
        kind = "predicted"
    method = str(meta.get("method") or parsed.get("method") or "")
    resolution = meta.get("resolution_angstrom")
    if resolution is None:
        resolution = parsed.get("resolution_angstrom")
    chains = list(parsed.get("chains") or [])
    if not chains:
        raise StructureError("Structure has no polymer chains.", "PARSING_ERROR")
    selected = _select_chain(chains, chain_id, molecule=kind_mol)
    models = list(parsed.get("models") or [1])
    declared_model = model_number if model_number is not None else models[0]
    if declared_model not in models:
        raise StructureError(
            f"Model {declared_model} is not present. Available models: {models}. "
            "NMR ensembles are not silently reduced to model 1.",
            "INVALID_INPUT",
        )
    selected = _chain_for_model(selected, declared_model)
    mapped = map_query_to_chain(query, selected, molecule=kind_mol)
    plddt = _plddt_from_chain(selected, kind)
    result = _build_structure_result(
        query=query,
        parsed=parsed,
        metadata=meta,
        source=src,
        kind=kind,
        method=method,
        resolution=resolution,
        selected=selected,
        mapping=mapped,
        model_number=declared_model,
        models=models,
        plddt=plddt,
        raw_text=text,
        parameters=parameters,
    )
    _remember(store, key, result)
    return result


def load_deposited_macromolecule(
    *,
    source: str,
    structure_id: str,
    structure_text: str = "",
    metadata: Optional[Mapping[str, Any]] = None,
    cif_url: str = "",
    urlopen_fn=None,
    cache: Optional[dict] = None,
) -> dict:
    """Carrega um mmCIF/PDB depositado sem exigir uma query proteica.

    Args:
        source: RCSB PDB.
        structure_id: PDB ID.
        structure_text: Ficheiro ja obtido (testes).
        metadata: Metadados da entry.
        cif_url: URL allowlisted.
        urlopen_fn: HTTP injetavel.
        cache: Cache injetavel.

    Returns:
        Dict parsed chains/atoms. Sem mapping de query. Sem XYZ inventados.

    Raises:
        StructureError.
    """
    parameters = {"file_preference": "mmCIF", "no_query": True, "cache_molecule": "COMPLEX"}
    key = cache_key(
        source=source,
        identifier=structure_id,
        chain="",
        parameters=parameters,
    )
    store = cache if cache is not None else _CACHE
    cached = store.get(key)
    if isinstance(cached, dict):
        return mark_cached_result(cached)
    meta = dict(metadata or {})
    src = str(source or "")
    if src == "RCSB PDB" and not meta:
        meta = fetch_rcsb_entry(structure_id, urlopen_fn=urlopen_fn)
    text = structure_text
    if not text:
        url = cif_url or str(meta.get("cif_url") or "")
        if not url and src == "RCSB PDB":
            url = f"{RCSB_FILES}/download/{str(structure_id).upper()}.cif"
        if not url:
            raise StructureError("No allowlisted structure file URL is available.", "INVALID_INPUT")
        text = fetch_structure_file(url, urlopen_fn=urlopen_fn)
    parsed = parse_structure_text(text, hint="mmcif")
    kind = str(meta.get("kind") or parsed.get("kind") or "experimental")
    parsed["kind"] = kind
    parsed["source"] = src
    parsed["structure_id"] = str(structure_id).upper()
    parsed["metadata"] = meta
    parsed["method"] = str(meta.get("method") or parsed.get("method") or "")
    if meta.get("resolution_angstrom") is not None:
        parsed["resolution_angstrom"] = meta.get("resolution_angstrom")
    parsed["retrieved_at"] = provenance.utc_now()
    parsed["structure_hash"] = parsed.get("content_hash")
    first_seq = ""
    chains = list(parsed.get("chains") or [])
    if chains:
        first_seq = str(chains[0].get("sequence") or "")
        parsed["chain_id"] = chains[0].get("chain_id")
    parsed["sequence_hash"] = provenance.sequence_digest(first_seq) if first_seq else ""
    _remember(store, key, parsed)
    return parsed


def _without_network_fields(value: Any) -> Any:
    """Remove chaves de rede/credencial do envelope enviado ao renderer."""
    blocked = ("url", "token", "api_key", "password", "secret", "http")
    if isinstance(value, Mapping):
        cleaned: Dict[str, Any] = {}
        for key, item in value.items():
            low = str(key).lower()
            if any(part in low for part in blocked):
                continue
            cleaned[str(key)] = _without_network_fields(item)
        return cleaned
    if isinstance(value, list):
        return [_without_network_fields(item) for item in value]
    if isinstance(value, tuple):
        return [_without_network_fields(item) for item in value]
    return value


def structure_kind_label(kind: object) -> str:
    """Rotulo cientifico do kind. Nao inventa experimental vs predicted.

    Args:
        kind: experimental, predicted, illustrative, unavailable ou outro.

    Returns:
        Rotulo em ingles, ou N/A se kind estiver vazio.

    Raises:
        Nenhum.
    """
    key = str(kind or "").strip().lower()
    if key == "experimental":
        return "Experimental Structure"
    if key == "predicted":
        return "Predicted Structure"
    if key == "illustrative":
        return "Illustrative Structure"
    if key == "unavailable":
        return "Coordinates Unavailable"
    if not key:
        return "N/A"
    return key.replace("_", " ")


def structure_3d_input(result: Mapping[str, Any]) -> dict:
    """Envelope isolado para o renderer 3D. Nao renderiza e nao faz fetch.

    Args:
        result: Envelope estrutural validado.

    Returns:
        Structure3DInput (copia profunda). O renderer nao deve alterar o
        envelope original. O campo renderer e data-only: este contrato nao
        desenha atomos.

    Raises:
        Nenhum.

    Nota biologica:
        Coordenadas, mapping e proveniencia sao copiados do resultado
        validado. Nenhuma geometria e inventada aqui.
    """
    atoms = copy.deepcopy(list(result.get("atoms") or []))
    selected = result.get("selected_chain") or {}
    mappings = copy.deepcopy(list(result.get("residue_mapping") or []))
    missing = [
        int(item.get("query_index_0based"))
        for item in mappings
        if not item.get("has_coordinates")
    ]
    scene = {
        "source": result.get("source"),
        "kind": result.get("kind"),
        "kind_label": structure_kind_label(result.get("kind")),
        "structure_id": result.get("structure_id"),
        "sequence_hash": result.get("sequence_hash"),
        "structure_hash": result.get("content_hash"),
        "chain_id": selected.get("chain_id"),
        "model_number": result.get("model_number"),
        "models": copy.deepcopy(result.get("models")),
        "chains": copy.deepcopy(result.get("chains")),
        "n_atoms": result.get("n_atoms"),
        "n_residues": len(list(selected.get("residues") or [])),
        "atoms": atoms,
        "molecule_type": (
            result.get("molecule")
            or result.get("molecule_type")
            or (result.get("parameters") or {}).get("molecule")
        ),
        "residues": copy.deepcopy(selected.get("residues")),
        "coordinates": [
            {
                "atom": item.get("atom_name"),
                "element": item.get("element"),
                "residue": item.get("comp_id"),
                "chain": item.get("auth_asym_id") or item.get("label_asym_id"),
                "label_seq_id": item.get("label_seq_id"),
                "auth_seq_id": item.get("auth_seq_id"),
                "insertion_code": item.get("insertion_code") or "",
                "alt_id": item.get("alt_id") or "",
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
        "mapping_status": result.get("mapping_status"),
        "mapping_status_note": result.get("mapping_status_note"),
        "coverage_query_with_coordinates": result.get("coverage_query_with_coordinates"),
        "coverage_formula": result.get("coverage_formula"),
        "n_with_coordinates": result.get("n_with_coordinates"),
        "n_without_coordinates": result.get("n_without_coordinates"),
        "residues_without_coordinates": missing,
        "confidence_data": copy.deepcopy(result.get("confidence_data")),
        "metadata": {
            "method": result.get("method"),
            "resolution_angstrom": result.get("resolution_angstrom"),
            "models": copy.deepcopy(result.get("models")),
            "model_note": result.get("model_note"),
            "alt_location_policy": result.get("alt_location_policy"),
            "disclaimer": result.get("disclaimer"),
        },
        "provenance": copy.deepcopy(
            result.get("provenance")
            or {
                "source": result.get("source"),
                "structure_id": result.get("structure_id"),
                "sequence_hash": result.get("sequence_hash"),
                "kind": result.get("kind"),
                "retrieved_at": result.get("retrieved_at"),
            }
        ),
        "renderer": "data-only",
    }
    return _without_network_fields(scene)


def freeze_structure_3d_input(envelope: Mapping[str, Any]) -> MappingProxyType:
    """Vista somente-leitura do envelope 3D. Mutacao levanta TypeError.

    Args:
        envelope: Saida de structure_3d_input.

    Returns:
        MappingProxyType de uma copia profunda.

    Raises:
        Nenhum.

    Nota:
        O contrato Structure3DInput continua a ser um dict serializavel na
        sessao. freeze_structure_3d_input e a verificacao de imutabilidade.
    """
    return MappingProxyType(copy.deepcopy(dict(envelope)))


def export_mapping_csv(result: Mapping[str, Any]) -> str:
    """CSV do mapping. N/A para ausente; 0 permanece 0.

    Args:
        result: Envelope.

    Returns:
        Texto CSV.

    Raises:
        Nenhum.
    """
    import csv

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=[
            "query_index_0based",
            "query_residue",
            "chain_id",
            "label_seq_id",
            "auth_seq_id",
            "insertion_code",
            "has_coordinates",
            "coordinate_status",
            "mapping_row_status",
            "mapping_status",
            "coverage_query_with_coordinates",
            "identity_ungapped_pct",
            "structure_id",
            "source",
            "kind",
            "method",
            "resolution_angstrom",
            "sequence_hash",
            "status",
        ],
    )
    writer.writeheader()
    coverage = provenance.csv_cell(result.get("coverage_query_with_coordinates"))
    identity = provenance.csv_cell(result.get("identity_ungapped_pct"))
    resolution = provenance.csv_cell(result.get("resolution_angstrom"))
    rows = list(result.get("residue_mapping") or [])
    mapping_status = result.get("mapping_status") or "N/A"
    if not rows:
        writer.writerow(
            {
                "query_index_0based": "N/A",
                "query_residue": "N/A",
                "chain_id": "N/A",
                "label_seq_id": "N/A",
                "auth_seq_id": "N/A",
                "insertion_code": "N/A",
                "has_coordinates": "N/A",
                "coordinate_status": "N/A",
                "mapping_row_status": "N/A",
                "mapping_status": mapping_status,
                "coverage_query_with_coordinates": coverage,
                "identity_ungapped_pct": identity,
                "structure_id": result.get("structure_id"),
                "source": result.get("source"),
                "kind": result.get("kind"),
                "method": result.get("method"),
                "resolution_angstrom": resolution,
                "sequence_hash": result.get("sequence_hash"),
                "status": result.get("status"),
            }
        )
        return buffer.getvalue()
    for item in rows:
        writer.writerow(
            {
                "query_index_0based": item.get("query_index_0based"),
                "query_residue": item.get("query_residue"),
                "chain_id": item.get("chain_id"),
                "label_seq_id": provenance.csv_cell(item.get("label_seq_id")),
                "auth_seq_id": provenance.csv_cell(item.get("auth_seq_id")),
                "insertion_code": item.get("insertion_code") or "",
                "has_coordinates": item.get("has_coordinates"),
                "coordinate_status": item.get("coordinate_status"),
                "mapping_row_status": item.get("mapping_row_status") or mapping_status,
                "mapping_status": mapping_status,
                "coverage_query_with_coordinates": coverage,
                "identity_ungapped_pct": identity,
                "structure_id": result.get("structure_id"),
                "source": result.get("source"),
                "kind": result.get("kind"),
                "method": result.get("method"),
                "resolution_angstrom": resolution,
                "sequence_hash": result.get("sequence_hash"),
                "status": result.get("status"),
            }
        )
    return buffer.getvalue()


def export_result_bundle(result: Mapping[str, Any]) -> dict:
    """JSON-safe. NaN nao vira 0.

    Args:
        result: Envelope.

    Returns:
        Dict serializavel sem o texto bruto gigante duplicado alem de content_hash.

    Raises:
        Nenhum.
    """
    copied = dict(result)
    copied.pop("raw_structure_text", None)
    return provenance.json_safe(copied)


def _build_structure_result(
    *,
    query: str,
    parsed: Mapping[str, Any],
    metadata: Mapping[str, Any],
    source: str,
    kind: str,
    method: str,
    resolution: Optional[float],
    selected: Mapping[str, Any],
    mapping: Mapping[str, Any],
    model_number: int,
    models: Sequence[int],
    plddt: Optional[dict],
    raw_text: str,
    parameters: Mapping[str, Any],
) -> dict:
    structure_id = str(metadata.get("structure_id") or parsed.get("entry_id") or "")
    sequence_hash = provenance.sequence_digest(query)
    chain_seq = str(selected.get("sequence") or "")
    status = "RETRIEVED"
    envelope = provenance.analysis_envelope(
        module="PROTEIN structure",
        payload={
            "structure_id": structure_id,
            "kind": kind,
            "method": method,
            "n_atoms": parsed.get("n_atoms"),
            "coverage_query_with_coordinates": mapping.get("coverage_query_with_coordinates"),
        },
        status=status,
        algorithm=f"{source} {method}".strip(),
        parameters=dict(parameters),
        sequence=query,
        source=source,
        accession=str(metadata.get("uniprot") or structure_id),
        organism=str(metadata.get("organism") or ""),
        input_identifier=structure_id,
    )
    envelope.update(
        {
            "status": status,
            "kind": kind,
            "structure_kind": kind,
            "source": source,
            "structure_source": source,
            "sequence_source": "analysis sequence",
            "structure_id": structure_id,
            "method": method,
            "resolution_angstrom": resolution,
            "tool": source,
            "tool_version": "" if metadata.get("latest_version") is None else str(metadata.get("latest_version")),
            "database": source,
            "release_date": metadata.get("release_date") or metadata.get("model_created_date") or "",
            "sequence": query,
            "sequence_hash": sequence_hash,
            "structure_sequence": chain_seq,
            "structure_sequence_hash": provenance.sequence_digest(chain_seq) if chain_seq else "",
            "sequences_identical": mapping.get("sequences_identical"),
            "chains": parsed.get("chains_summary"),
            "selected_chain": {
                "chain_id": selected.get("chain_id"),
                "label_asym_id": selected.get("label_asym_id"),
                "auth_asym_id": selected.get("auth_asym_id"),
                "length": len(chain_seq),
                "sequence": chain_seq,
                "polymer_type": selected.get("polymer_type"),
                "residues": selected.get("residues"),
            },
            "models": list(models),
            "model_number": int(model_number),
            "model_note": (
                "Multiple models are present. The displayed model is selected explicitly; "
                "model 1 is not treated as the structure."
                if len(models) > 1
                else "Single model in this file."
            ),
            "residue_mapping": mapping.get("mapping"),
            "mapping_method": mapping.get("mapping_method"),
            "mapping_status": mapping.get("mapping_status"),
            "mapping_status_note": mapping.get("mapping_status_note"),
            "mapped_query_indices": list(mapping.get("mapped_query_indices") or []),
            "unmapped_query_indices": list(mapping.get("unmapped_query_indices") or []),
            "n_letter_matches": mapping.get("n_letter_matches"),
            "polymer_index_mode": mapping.get("polymer_index_mode"),
            "identity_ungapped_pct": mapping.get("identity_ungapped_pct"),
            "coverage_query_with_coordinates": mapping.get("coverage_query_with_coordinates"),
            "coverage_formula": "n_query_residues_with_coordinates / len(analysis_sequence)",
            "n_with_coordinates": mapping.get("n_with_coordinates"),
            "n_without_coordinates": mapping.get("n_without_coordinates"),
            "atoms": parsed.get("atoms"),
            "n_atoms": parsed.get("n_atoms"),
            "content_hash": parsed.get("content_hash"),
            "structure_hash": parsed.get("content_hash"),
            "deposited_secondary_structure": parsed.get("deposited_secondary_structure"),
            "confidence_data": plddt,
            "alt_location_policy": (
                "Alternate locations are retained. The residue CA used for mapping is the "
                "highest occupancy among recorded alts; this is not a claim of biological truth."
            ),
            "retrieved_at": provenance.utc_now(),
            "cache_status": "live",
            "coordinate_system": "0-based query index; structure uses deposited auth_seq_id",
            "disclaimer": source_disclaimer(kind, source),
            "raw_structure_text": raw_text,
        }
    )
    return envelope


def _plddt_from_chain(chain: Mapping[str, Any], kind: str) -> Optional[dict]:
    if kind != "predicted":
        return None
    values: List[float] = []
    for residue in list(chain.get("residues") or []):
        b_iso = residue.get("ca_b_iso")
        number = _first_finite(b_iso)
        if number is not None:
            values.append(number)
    if not values:
        return None
    return {
        "metric": "pLDDT stored in B-factor (AlphaFold convention when kind is predicted)",
        "per_residue": values,
        "mean": sum(values) / len(values),
        "n": len(values),
        "range": "0-100 as deposited, not a probability that the structure is correct",
    }


def _select_chain(chains: Sequence[dict], chain_id: str, *, molecule: str = "PROTEIN") -> dict:
    mol = str(molecule or "PROTEIN").strip().upper()
    typed: List[dict] = []
    for item in chains:
        ptype = str(item.get("polymer_type") or "").lower()
        if mol == "DNA" and ("deoxy" in ptype or ptype == "dna"):
            typed.append(item)
        elif mol == "RNA" and (("ribo" in ptype and "deoxy" not in ptype) or ptype == "rna"):
            typed.append(item)
        elif mol == "PROTEIN" and ptype.startswith("polypeptide"):
            typed.append(item)
    pool = typed if typed else ([] if mol in {"DNA", "RNA"} else list(chains))
    if mol in {"DNA", "RNA"} and not pool:
        raise StructureError(
            f"This deposited file has no {mol} polymer chain to map.",
            "NO_STRUCTURE",
        )
    if not pool:
        pool = list(chains)
    if chain_id:
        for item in pool:
            if str(item.get("chain_id")) == chain_id or str(item.get("auth_asym_id")) == chain_id:
                return item
        raise StructureError(f"Chain {chain_id} is not present.", "INVALID_INPUT")
    return pool[0]


def _chain_for_model(chain: Mapping[str, Any], model_number: int) -> dict:
    residues = [
        item for item in list(chain.get("residues") or []) if int(item.get("model") or 1) == int(model_number)
    ]
    copied = dict(chain)
    copied["residues"] = residues
    return copied


def _dedupe_hits(hits: Sequence[dict]) -> List[dict]:
    seen: set[tuple] = set()
    unique: List[dict] = []
    for item in hits:
        key = (item.get("source"), str(item.get("structure_id") or "").upper(), item.get("kind"))
        if key in seen:
            continue
        seen.add(key)
        unique.append(dict(item))
    return unique


def _split_polymer_entity(identifier: str) -> Tuple[str, str]:
    text = str(identifier or "")
    if "_" in text:
        pdb_id, entity = text.split("_", 1)
        classified = classify_structure_identifier(pdb_id)
        if classified["kind"] == "pdb":
            return classified["canonical"], entity
    classified = classify_structure_identifier(text)
    if classified["kind"] == "pdb":
        return classified["canonical"], ""
    return "", ""


def _remember(store: dict, key: str, result: dict) -> None:
    store[key] = result
    extra = list(store.keys())[:-MAX_CACHE] if len(store) > MAX_CACHE else []
    for old in extra:
        store.pop(old, None)


def _structure_from_cif_block(block: Mapping[str, Any], encoded: bytes) -> dict:
    items = dict(block.get("items") or {})
    loops = list(block.get("loops") or [])
    entry_id = str(items.get("_entry.id") or "").strip()
    method = str(items.get("_exptl.method") or "").strip().strip("'\"")
    resolution = _first_finite(items.get("_refine.ls_d_res_high"))
    polymer_sequences: Dict[str, str] = {}
    polymer_types: Dict[str, str] = {}
    entity_descriptions: Dict[str, str] = {}
    entity_id = str(items.get("_entity_poly.entity_id") or "1")
    seq = _strip_cif_sequence(str(items.get("_entity_poly.pdbx_seq_one_letter_code_can") or items.get("_entity_poly.pdbx_seq_one_letter_code") or ""))
    ptype = str(items.get("_entity_poly.type") or "")
    strand = str(items.get("_entity_poly.pdbx_strand_id") or "")
    if seq:
        polymer_sequences[entity_id] = seq
        polymer_types[entity_id] = ptype
        if strand:
            for chain_id in strand.replace(" ", "").split(","):
                if chain_id:
                    polymer_sequences[f"strand:{chain_id}"] = seq
    for loop in loops:
        tags = list(loop.get("tags") or [])
        if tags and tags[0].startswith("_entity_poly."):
            for row in loop.get("rows") or []:
                rec = dict(zip(tags, row))
                eid = str(rec.get("_entity_poly.entity_id") or "")
                one = _strip_cif_sequence(
                    str(rec.get("_entity_poly.pdbx_seq_one_letter_code_can") or rec.get("_entity_poly.pdbx_seq_one_letter_code") or "")
                )
                if eid and one:
                    polymer_sequences[eid] = one
                    polymer_types[eid] = str(rec.get("_entity_poly.type") or "")
        if tags and any(tag == "_entity.id" or tag.startswith("_entity.") for tag in tags) and "_entity.id" in tags:
            for row in loop.get("rows") or []:
                rec = dict(zip(tags, row))
                eid = str(rec.get("_entity.id") or "")
                desc = str(rec.get("_entity.pdbx_description") or "").strip().strip("'\"")
                if eid and desc:
                    entity_descriptions[eid] = desc
    atom_rows: List[dict] = []
    for loop in loops:
        tags = [str(tag) for tag in (loop.get("tags") or [])]
        if not any(tag.startswith("_atom_site.") for tag in tags):
            continue
        for values in loop.get("rows") or []:
            rec = dict(zip(tags, values))
            atom_rows.append(rec)
    if not atom_rows:
        raise StructureError("mmCIF has no _atom_site loop.", "PARSING_ERROR")
    atoms: List[dict] = []
    for rec in atom_rows:
        group = str(rec.get("_atom_site.group_PDB") or "ATOM")
        try:
            atom = {
                "group": group,
                "atom_name": str(rec.get("_atom_site.label_atom_id") or rec.get("_atom_site.auth_atom_id") or ""),
                "alt_id": _blank_cif(rec.get("_atom_site.label_alt_id")),
                "comp_id": str(rec.get("_atom_site.label_comp_id") or rec.get("_atom_site.auth_comp_id") or ""),
                "label_asym_id": str(rec.get("_atom_site.label_asym_id") or ""),
                "auth_asym_id": str(rec.get("_atom_site.auth_asym_id") or rec.get("_atom_site.label_asym_id") or ""),
                "entity_id": str(rec.get("_atom_site.label_entity_id") or "1"),
                "label_seq_id": _optional_int(rec.get("_atom_site.label_seq_id")),
                "auth_seq_id": _optional_int(rec.get("_atom_site.auth_seq_id")),
                "insertion_code": _blank_cif(rec.get("_atom_site.pdbx_PDB_ins_code")),
                "x": float(rec.get("_atom_site.Cartn_x")),
                "y": float(rec.get("_atom_site.Cartn_y")),
                "z": float(rec.get("_atom_site.Cartn_z")),
                "occupancy": _optional_float(rec.get("_atom_site.occupancy")),
                "b_iso": _optional_float(rec.get("_atom_site.B_iso_or_equiv")),
                "element": str(rec.get("_atom_site.type_symbol") or ""),
                "model": int(rec.get("_atom_site.pdbx_PDB_model_num") or 1),
            }
        except (TypeError, ValueError) as exc:
            raise StructureError("mmCIF atom coordinates are not numeric.", "PARSING_ERROR") from exc
        if not scientific_checks.atom_coordinate_is_finite(atom["x"], atom["y"], atom["z"]):
            raise StructureError("mmCIF atom coordinates are not finite.", "PARSING_ERROR")
        atoms.append(atom)
        if len(atoms) > MAX_ATOMS:
            raise StructureError(
                f"Structure exceeds {MAX_ATOMS:,} atoms (RESOURCE_LIMIT).",
                "RESOURCE_LIMIT",
            )
    return _structure_from_atoms(
        atoms,
        entry_id=entry_id,
        method=method,
        resolution=resolution,
        polymer_sequences=polymer_sequences,
        polymer_types=polymer_types,
        entity_descriptions=entity_descriptions,
        content_hash=hashlib.sha256(encoded).hexdigest(),
        fmt="mmcif",
        deposited_secondary_structure=_extract_deposited_secondary_structure(loops),
    )


def _structure_from_atoms(
    atoms: Sequence[dict],
    *,
    entry_id: str,
    method: str,
    resolution: Optional[float],
    polymer_sequences: Mapping[str, str],
    polymer_types: Mapping[str, str],
    content_hash: str,
    fmt: str,
    entity_descriptions: Optional[Mapping[str, str]] = None,
    deposited_secondary_structure: Optional[Mapping[str, Any]] = None,
) -> dict:
    models = sorted({int(atom.get("model") or 1) for atom in atoms})
    chains_map: Dict[tuple, dict] = {}
    for atom in atoms:
        if str(atom.get("group")) != "ATOM":
            continue
        key = (str(atom.get("auth_asym_id") or atom.get("label_asym_id")), int(atom.get("model") or 1))
        chain = chains_map.setdefault(
            key,
            {
                "chain_id": key[0],
                "auth_asym_id": str(atom.get("auth_asym_id") or key[0]),
                "label_asym_id": str(atom.get("label_asym_id") or key[0]),
                "entity_id": str(atom.get("entity_id") or "1"),
                "model": key[1],
                "residues_index": {},
            },
        )
        label_seq = atom.get("label_seq_id")
        auth_seq = atom.get("auth_seq_id")
        ins = atom.get("insertion_code") or ""
        res_key = (label_seq, auth_seq, ins, atom.get("alt_id") or "")
        residue = chain["residues_index"].setdefault(
            (label_seq, auth_seq, ins),
            {
                "label_seq_id": label_seq,
                "auth_seq_id": auth_seq,
                "insertion_code": ins,
                "comp_id": atom.get("comp_id"),
                "one_letter": residue_one_letter(str(atom.get("comp_id") or "")),
                "has_coordinates": True,
                "model": key[1],
                "alt_ids": [],
                "ca": None,
                "p": None,
                "c1": None,
                "ca_b_iso": None,
                "occupancies": [],
            },
        )
        alt = atom.get("alt_id") or ""
        if alt not in residue["alt_ids"]:
            residue["alt_ids"].append(alt)
        if atom.get("occupancy") is not None:
            residue["occupancies"].append(float(atom["occupancy"]))
        if str(atom.get("atom_name")) in {"CA", "C"} and residue["ca"] is None:
            residue["ca"] = {"x": atom["x"], "y": atom["y"], "z": atom["z"]}
            residue["ca_b_iso"] = atom.get("b_iso")
        elif str(atom.get("atom_name")) == "CA":
            current_occ = max(residue["occupancies"]) if residue["occupancies"] else None
            this_occ = atom.get("occupancy")
            if this_occ is not None and (current_occ is None or float(this_occ) >= current_occ):
                residue["ca"] = {"x": atom["x"], "y": atom["y"], "z": atom["z"]}
                residue["ca_b_iso"] = atom.get("b_iso")
        if str(atom.get("atom_name")) == "P" and residue.get("p") is None:
            residue["p"] = {"x": atom["x"], "y": atom["y"], "z": atom["z"]}
        if str(atom.get("atom_name")) in {"C1'", "C1*"} and residue.get("c1") is None:
            residue["c1"] = {"x": atom["x"], "y": atom["y"], "z": atom["z"]}
        _ = res_key
    chains: List[dict] = []
    for (_cid, _model), chain in sorted(chains_map.items(), key=lambda item: (item[0][1], item[0][0])):
        residues = sorted(
            chain["residues_index"].values(),
            key=lambda item: (
                int(item["label_seq_id"] or 0),
                int(item["auth_seq_id"] or 0),
                str(item["insertion_code"] or ""),
            ),
        )
        if len(residues) > MAX_RESIDUES:
            raise StructureError(
                f"Structure exceeds {MAX_RESIDUES:,} residues (RESOURCE_LIMIT).",
                "RESOURCE_LIMIT",
            )
        entity_id = str(chain.get("entity_id") or "1")
        sequence = str(polymer_sequences.get(entity_id) or polymer_sequences.get(f"strand:{chain['chain_id']}") or "")
        if not sequence:
            sequence = "".join(str(item.get("one_letter") or "X") for item in residues)
        if sequence and _labels_are_polymer_indices(residues, len(sequence)):
            present = {
                int(item["label_seq_id"])
                for item in residues
                if item.get("label_seq_id") is not None
            }
            for index, aa in enumerate(sequence):
                label = index + 1
                if label not in present:
                    residues.append(
                        {
                            "label_seq_id": label,
                            "auth_seq_id": None,
                            "insertion_code": "",
                            "comp_id": "",
                            "one_letter": aa,
                            "has_coordinates": False,
                            "model": chain["model"],
                            "alt_ids": [],
                            "ca": None,
                            "p": None,
                            "c1": None,
                            "ca_b_iso": None,
                            "occupancies": [],
                        }
                    )
            residues.sort(key=lambda item: int(item.get("label_seq_id") or 0))
        polymer = _ordered_polymer_residues({"sequence": sequence, "residues": residues})
        for index, item in enumerate(polymer):
            if item is not None:
                item["polymer_index_0based"] = index
        chains.append(
            {
                "chain_id": chain["chain_id"],
                "auth_asym_id": chain["auth_asym_id"],
                "label_asym_id": chain["label_asym_id"],
                "entity_id": entity_id,
                "model": chain["model"],
                "sequence": sequence,
                "polymer_type": polymer_types.get(entity_id, "polypeptide(L)"),
                "description": str((entity_descriptions or {}).get(entity_id) or ""),
                "residues": residues,
                "length": len(sequence),
            }
        )
    summary = [
        {
            "chain_id": item["chain_id"],
            "length": item["length"],
            "model": item["model"],
            "polymer_type": item["polymer_type"],
        }
        for item in chains
    ]
    kind = "experimental" if method and "THEORETICAL" not in method.upper() else "retrieved"
    return {
        "format": fmt,
        "entry_id": entry_id,
        "method": method,
        "resolution_angstrom": resolution,
        "kind": kind,
        "atoms": list(atoms),
        "n_atoms": len(atoms),
        "models": models,
        "chains": chains,
        "chains_summary": summary,
        "content_hash": content_hash,
        "entity_descriptions": dict(entity_descriptions or {}),
        "deposited_secondary_structure": dict(deposited_secondary_structure or {}),
    }


def _parse_cif_data(text: str) -> dict:
    tokens = _cif_tokens(text)
    items: Dict[str, str] = {}
    loops: List[dict] = []
    index = 0
    while index < len(tokens):
        token, quoted = tokens[index]
        if token.startswith("data_") and not quoted:
            index += 1
            continue
        if token == "loop_" and not quoted:
            index += 1
            tags: List[str] = []
            while index < len(tokens) and _is_cif_tag(tokens[index]):
                tags.append(tokens[index][0])
                index += 1
            values: List[str] = []
            while index < len(tokens) and not _is_cif_tag(tokens[index]) and not _is_cif_block(tokens[index]):
                values.append(tokens[index][0])
                index += 1
            if not tags:
                raise StructureError("mmCIF loop has no tags.", "PARSING_ERROR")
            if len(values) % len(tags) != 0:
                raise StructureError("mmCIF loop row length does not match tags.", "PARSING_ERROR")
            rows = [
                values[offset : offset + len(tags)]
                for offset in range(0, len(values), len(tags))
            ]
            loops.append({"tags": tags, "rows": rows})
            continue
        if _is_cif_tag((token, quoted)):
            if index + 1 >= len(tokens):
                raise StructureError("mmCIF tag is missing a value.", "PARSING_ERROR")
            items[token] = tokens[index + 1][0]
            index += 2
            continue
        index += 1
    return {"items": items, "loops": loops}


def _is_cif_tag(token: Tuple[str, bool]) -> bool:
    text, quoted = token
    return (not quoted) and text.startswith("_") and "." in text


def _is_cif_block(token: Tuple[str, bool]) -> bool:
    text, quoted = token
    if quoted:
        return False
    return text in {"loop_", "global_"} or text.startswith("data_") or text.startswith("save_")


def _cif_tokens(text: str) -> List[Tuple[str, bool]]:
    tokens: List[Tuple[str, bool]] = []
    length = len(text)
    index = 0
    line_start = True
    while index < length:
        char = text[index]
        if char in " \t\r":
            index += 1
            continue
        if char == "\n":
            index += 1
            line_start = True
            continue
        if char == "#" and line_start:
            while index < length and text[index] != "\n":
                index += 1
            continue
        if char == ";" and line_start:
            index += 1
            start = index
            while index < length:
                if text[index] == "\n" and index + 1 < length and text[index + 1] == ";":
                    tokens.append((text[start:index], True))
                    index += 2
                    line_start = False
                    break
                index += 1
            else:
                raise StructureError("Unterminated mmCIF semicolon text field.", "PARSING_ERROR")
            continue
        if char in {'"', "'"}:
            quote = char
            index += 1
            start = index
            while index < length and text[index] != quote:
                index += 1
            if index >= length:
                raise StructureError("Unterminated mmCIF quoted string.", "PARSING_ERROR")
            tokens.append((text[start:index], True))
            index += 1
            line_start = False
            continue
        start = index
        while index < length and text[index] not in " \t\r\n":
            index += 1
        tokens.append((text[start:index], False))
        line_start = False
    return tokens


class _AllowlistRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not url_is_allowed(newurl):
            raise StructureError(
                "Refusing HTTP redirect off the structure allowlist.",
                "INVALID_INPUT",
            )
        return HTTPRedirectHandler.redirect_request(self, req, fp, code, msg, headers, newurl)


def _http_json(
    url: str,
    *,
    method: str,
    json_body: Optional[Mapping[str, Any]] = None,
    urlopen_fn=None,
    empty_ok: bool = False,
) -> Any:
    text = _http_text(
        url,
        method=method,
        json_body=json_body,
        urlopen_fn=urlopen_fn,
        limit=MAX_METADATA_BYTES,
        empty_ok=empty_ok,
    )
    if empty_ok and not text:
        return {}
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise StructureError("Remote JSON is malformed.", "PARSING_ERROR") from exc


def _encode_multipart(fields: Mapping[str, tuple[str, str]]) -> tuple[bytes, str]:
    """Encode a small multipart/form-data body. Filenames are not user-controlled.

    Args:
        fields: name -> (filename or empty, value). Empty filename is a text field.

    Returns:
        (body_bytes, Content-Type header value).

    Raises:
        StructureError: INVALID_INPUT se o conjunto estiver vazio.
    """
    if not fields:
        raise StructureError("Multipart body is empty.", "INVALID_INPUT")
    boundary = "----HelixScope" + secrets.token_hex(16)
    chunks: list[bytes] = []
    for name, (filename, value) in fields.items():
        safe_name = str(name or "").replace('"', "")
        safe_file = str(filename or "").replace('"', "").replace("\\", "").replace("/", "")
        payload = str(value or "").encode("utf-8")
        disposition = f'Content-Disposition: form-data; name="{safe_name}"'
        if safe_file:
            disposition += f'; filename="{safe_file}"'
        header = f"--{boundary}\r\n{disposition}\r\n"
        if safe_file:
            header += "Content-Type: application/octet-stream\r\n"
        header += "\r\n"
        chunks.append(header.encode("utf-8"))
        chunks.append(payload)
        chunks.append(b"\r\n")
    chunks.append(f"--{boundary}--\r\n".encode("utf-8"))
    return b"".join(chunks), f"multipart/form-data; boundary={boundary}"


def _http_text(
    url: str,
    *,
    method: str,
    json_body: Optional[Mapping[str, Any]] = None,
    form_body: Optional[Mapping[str, str]] = None,
    multipart_fields: Optional[Mapping[str, tuple[str, str]]] = None,
    urlopen_fn=None,
    limit: int,
    empty_ok: bool = False,
    timeout_s: Optional[float] = None,
) -> str:
    bodies = [item for item in (json_body, form_body, multipart_fields) if item is not None]
    if len(bodies) > 1:
        raise StructureError("HTTP body cannot mix JSON, form, and multipart.", "INVALID_INPUT")
    if not url_is_allowed(url):
        raise StructureError("Refusing a non-allowlisted structure URL.", "INVALID_INPUT")
    data = None
    headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if json_body is not None:
        data = json.dumps(json_body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    elif form_body is not None:
        data = urlencode(dict(form_body)).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    elif multipart_fields is not None:
        data, content_type = _encode_multipart(multipart_fields)
        headers["Content-Type"] = content_type
    request = Request(url, data=data, method=method, headers=headers)
    opener = urlopen_fn
    if opener is None:
        opener = build_opener(_AllowlistRedirect()).open
    wait = float(HTTP_TIMEOUT_S if timeout_s is None else timeout_s)
    try:
        with opener(request, timeout=wait) as handle:
            raw = handle.read()
            final = str(getattr(handle, "geturl", lambda: url)() or url)
    except HTTPError as exc:
        raise _http_error(exc) from exc
    except StructureError:
        raise
    except (URLError, TimeoutError, socket.timeout, OSError) as exc:
        message = str(exc).lower()
        if "timed out" in message or isinstance(exc, (TimeoutError, socket.timeout)):
            raise StructureError(
                "Structure request timed out. This is TIMEOUT, not a missing structure.",
                "TIMEOUT",
            ) from exc
        raise StructureError(
            f"Structure network error: {exc}. This is NETWORK_ERROR, not NOT_FOUND.",
            "NETWORK_ERROR",
        ) from exc
    if final and not url_is_allowed(final):
        raise StructureError("Refusing a non-allowlisted structure URL.", "INVALID_INPUT")
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    if raw[:2] == b"\x1f\x8b":
        try:
            raw = gzip.decompress(raw)
        except OSError as exc:
            raise StructureError("Gzip structure payload is malformed.", "PARSING_ERROR") from exc
    if len(raw) > limit:
        raise StructureError(
            f"Structure response exceeds {limit:,} bytes.",
            "RESOURCE_LIMIT",
        )
    text = raw.decode("utf-8", errors="replace")
    if not text.strip() and not empty_ok:
        raise StructureError("Structure response is empty.", "PARSING_ERROR")
    return text


def _http_error(exc: HTTPError) -> StructureError:
    code = int(getattr(exc, "code", 0) or 0)
    if code == 204:
        return StructureError("No validated structure found.", "NO_STRUCTURE")
    if code == 404:
        return StructureError("Structure identifier was not found.", "NOT_FOUND")
    if code == 429:
        return StructureError("Structure service rate limited the request (HTTP 429).", "RATE_LIMITED")
    if code >= 500:
        return StructureError(f"Structure service unavailable (HTTP {code}).", "SERVICE_UNAVAILABLE")
    if code in {400, 415}:
        return StructureError(f"Structure service rejected the request (HTTP {code}).", "INVALID_INPUT")
    return StructureError(f"Structure HTTP {code}: {exc}.", "SERVICE_UNAVAILABLE")


def _optional_int(value: Any) -> Optional[int]:
    text = _blank_cif(value)
    if text == "":
        return None
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return None


def _optional_float(value: Any) -> Optional[float]:
    text = _blank_cif(value)
    if text == "":
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def _first_finite(value: Any) -> Optional[float]:
    return _optional_float(value)


def _blank_cif(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if text in {".", "?", ""}:
        return ""
    return text


def _strip_cif_sequence(value: str) -> str:
    return re.sub(r"\s+", "", _blank_cif(value))


def _parse_resolution_text(text: str) -> Optional[float]:
    cleaned = str(text or "").replace("A", " ").replace("Å", " ").strip()
    if cleaned in {"", "-", "N/A", "NA"}:
        return None
    return _first_finite(cleaned.split()[0] if cleaned else "")


def _extract_deposited_secondary_structure(loops: Sequence[Mapping[str, Any]]) -> dict:
    """Le _struct_conf e _struct_sheet_range do mmCIF. Nao e DSSP.

    Args:
        loops: Loops CIF ja tokenizados.

    Returns:
        Dict helices, sheets, source. Assignment do depositor, nao predicao.

    Raises:
        Nenhum.
    """
    helices: List[dict] = []
    sheets: List[dict] = []
    for loop in loops:
        tags = [str(tag) for tag in (loop.get("tags") or [])]
        if "_struct_conf.conf_type_id" in tags:
            for values in loop.get("rows") or []:
                rec = dict(zip(tags, values))
                conf = str(rec.get("_struct_conf.conf_type_id") or "")
                beg = _optional_int(rec.get("_struct_conf.beg_label_seq_id"))
                end = _optional_int(rec.get("_struct_conf.end_label_seq_id"))
                chain = str(rec.get("_struct_conf.beg_label_asym_id") or "")
                if beg is None or end is None:
                    continue
                kind = "helix" if conf.upper().startswith("HELX") else "other"
                helices.append(
                    {
                        "conf_type_id": conf,
                        "kind": kind,
                        "chain_id": chain,
                        "beg_label_seq_id": beg,
                        "end_label_seq_id": end,
                        "helix_id": str(rec.get("_struct_conf.pdbx_PDB_helix_id") or ""),
                    }
                )
        if "_struct_sheet_range.sheet_id" in tags:
            for values in loop.get("rows") or []:
                rec = dict(zip(tags, values))
                beg = _optional_int(rec.get("_struct_sheet_range.beg_label_seq_id"))
                end = _optional_int(rec.get("_struct_sheet_range.end_label_seq_id"))
                if beg is None or end is None:
                    continue
                sheets.append(
                    {
                        "sheet_id": str(rec.get("_struct_sheet_range.sheet_id") or ""),
                        "chain_id": str(rec.get("_struct_sheet_range.beg_label_asym_id") or ""),
                        "beg_label_seq_id": beg,
                        "end_label_seq_id": end,
                    }
                )
    return {
        "source": "mmCIF _struct_conf / _struct_sheet_range",
        "method": "depositor annotation in the coordinate file",
        "not_dssp": True,
        "not_sequence_prediction": True,
        "helices": helices,
        "sheets": sheets,
        "n_helix_records": len(helices),
        "n_sheet_records": len(sheets),
    }


def parse_rcsb_polymer_entity(payload: Mapping[str, Any]) -> dict:
    """Extrai tipo, sequencia, descricao e chain IDs do JSON Data API.

    Args:
        payload: JSON de /rest/v1/core/polymer_entity/{pdb}/{entity}.

    Returns:
        Dict entity fields. Sem coordenadas.

    Raises:
        StructureError: PARSING_ERROR.
    """
    if not isinstance(payload, Mapping):
        raise StructureError("Polymer entity JSON is not an object.", "PARSING_ERROR")
    poly = payload.get("entity_poly") if isinstance(payload.get("entity_poly"), Mapping) else {}
    ent = payload.get("rcsb_polymer_entity") if isinstance(payload.get("rcsb_polymer_entity"), Mapping) else {}
    ids = (
        payload.get("rcsb_polymer_entity_container_identifiers")
        if isinstance(payload.get("rcsb_polymer_entity_container_identifiers"), Mapping)
        else {}
    )
    seq = _strip_cif_sequence(
        str(poly.get("pdbx_seq_one_letter_code_can") or poly.get("pdbx_seq_one_letter_code") or "")
    )
    entry_id = str(ids.get("entry_id") or "").upper()
    entity_id = str(ids.get("entity_id") or "")
    if not seq and not str(poly.get("type") or ""):
        raise StructureError("Polymer entity JSON has no sequence or type.", "PARSING_ERROR")
    return {
        "entry_id": entry_id,
        "entity_id": entity_id,
        "polymer_type": str(poly.get("type") or ""),
        "rcsb_entity_polymer_type": str(poly.get("rcsb_entity_polymer_type") or ""),
        "sequence": seq,
        "length": len(seq),
        "description": str(ent.get("pdbx_description") or ""),
        "asym_ids": list(ids.get("asym_ids") or []),
        "auth_asym_ids": list(ids.get("auth_asym_ids") or []),
        "strand_id": str(poly.get("pdbx_strand_id") or ""),
        "source": "RCSB PDB Data API polymer_entity",
    }


def fetch_rcsb_polymer_entity(pdb_id: str, entity_id: str, *, urlopen_fn=None) -> dict:
    """GET Data API /rest/v1/core/polymer_entity/{pdb}/{entity}.

    Args:
        pdb_id: PDB ID.
        entity_id: Numero da entidade.
        urlopen_fn: HTTP injetavel.

    Returns:
        parse_rcsb_polymer_entity.

    Raises:
        StructureError.
    """
    classified = classify_structure_identifier(pdb_id)
    if classified["kind"] != "pdb":
        raise StructureError("PDB identifier is invalid.", "INVALID_INPUT")
    eid = str(entity_id or "").strip()
    if not eid.isdigit():
        raise StructureError("Polymer entity id must be a positive integer.", "INVALID_INPUT")
    url = f"{RCSB_DATA}/rest/v1/core/polymer_entity/{classified['canonical']}/{eid}"
    payload = _http_json(url, method="GET", urlopen_fn=urlopen_fn)
    if not isinstance(payload, Mapping):
        raise StructureError("RCSB polymer entity response is not JSON object.", "PARSING_ERROR")
    parsed = parse_rcsb_polymer_entity(payload)
    parsed["retrieved_at"] = provenance.utc_now()
    return parsed


def experimental_structure_provenance(result: Mapping[str, Any]) -> dict:
    """Oito campos de proveniencia estrutural preenchidos a partir do envelope.

    Args:
        result: Saida de load_protein_structure / load_deposited_macromolecule.

    Returns:
        Dict pdb_id, source, chain, method, resolution, sequence_hash,
        structure_hash, timestamp.

    Raises:
        StructureError: INVALID_INPUT se um campo obrigatorio estiver vazio.
    """
    selected = result.get("selected_chain") or {}
    record = {
        "pdb_id": str(result.get("structure_id") or result.get("entry_id") or ""),
        "source": str(result.get("source") or ""),
        "chain": str(selected.get("chain_id") or result.get("chain_id") or ""),
        "method": str(result.get("method") or ""),
        "resolution_angstrom": result.get("resolution_angstrom"),
        "sequence_hash": str(result.get("sequence_hash") or ""),
        "structure_hash": str(result.get("structure_hash") or result.get("content_hash") or ""),
        "timestamp": str(result.get("retrieved_at") or result.get("retrieval_timestamp") or ""),
    }
    missing = [
        key
        for key in ("pdb_id", "source", "chain", "method", "sequence_hash", "structure_hash", "timestamp")
        if not record.get(key)
    ]
    if missing:
        raise StructureError(
            "Structure provenance is incomplete: " + ", ".join(missing) + ".",
            "INVALID_INPUT",
        )
    return record


def check_structure_integrity(
    parsed: Mapping[str, Any],
    *,
    molecule: str,
    chain_id: str = "",
    expected_sequence: str = "",
) -> dict:
    """Seis checagens isoladas: coords, polimero, chain, ordem, comprimento, sequencia.

    Args:
        parsed: parse_mmcif/parse_pdb.
        molecule: DNA, RNA ou PROTEIN.
        chain_id: Chain exigida.
        expected_sequence: Sequencia esperada da chain (opcional).

    Returns:
        Dict ok, issues (lista vazia se ok).

    Raises:
        StructureError: PARSING_ERROR/INVALID_INPUT no primeiro defeito grave.
    """
    mol = str(molecule or "PROTEIN").strip().upper()
    issues: List[str] = []
    for atom in list(parsed.get("atoms") or []):
        if not scientific_checks.atom_coordinate_is_finite(atom.get("x"), atom.get("y"), atom.get("z")):
            raise StructureError("Atom coordinates are not finite.", "PARSING_ERROR")
    chains = list(parsed.get("chains") or [])
    if chain_id:
        found = [
            item
            for item in chains
            if str(item.get("chain_id")) == chain_id or str(item.get("auth_asym_id")) == chain_id
        ]
        if not found:
            raise StructureError(f"Chain {chain_id} is not present.", "INVALID_INPUT")
        selected = found[0]
    else:
        if not chains:
            raise StructureError("Structure has no polymer chains.", "PARSING_ERROR")
        selected = chains[0]
    ptype = str(selected.get("polymer_type") or "").lower()
    if mol == "DNA" and not ("deoxy" in ptype or ptype == "dna"):
        raise StructureError("Selected chain is not a DNA polymer.", "INVALID_INPUT")
    if mol == "RNA" and not (("ribo" in ptype and "deoxy" not in ptype) or ptype == "rna"):
        raise StructureError("Selected chain is not an RNA polymer.", "INVALID_INPUT")
    if mol == "PROTEIN" and ptype and not ptype.startswith("polypeptide"):
        raise StructureError("Selected chain is not a polypeptide.", "INVALID_INPUT")
    residues = list(selected.get("residues") or [])
    labels = [item.get("label_seq_id") for item in residues if item.get("label_seq_id") is not None]
    if labels != sorted(labels):
        raise StructureError("Residue label_seq_id order is not monotonic.", "PARSING_ERROR")
    sequence = str(selected.get("sequence") or "")
    if int(selected.get("length") or 0) != len(sequence):
        raise StructureError("Chain length does not match polymer sequence length.", "PARSING_ERROR")
    if expected_sequence and sequence != expected_sequence:
        raise StructureError("Deposited polymer sequence does not match the expected sequence.", "INVALID_INPUT")
    if mol == "DNA" and any(letter not in "ACGT" for letter in sequence):
        raise StructureError("DNA polymer sequence contains an invalid character.", "INVALID_INPUT")
    if mol == "RNA" and any(letter not in "ACGU" for letter in sequence):
        raise StructureError("RNA polymer sequence contains an invalid character.", "INVALID_INPUT")
    return {"ok": True, "issues": issues, "chain_id": selected.get("chain_id"), "length": len(sequence)}
