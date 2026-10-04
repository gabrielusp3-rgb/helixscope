"""US-align local structural alignment (proteins and nucleic acids).

US-align (Zhang et al., Nature Methods 2022; Nature Protocols 2025) alinha
proteinas, acidos nucleicos e complexos. Binario oficial Windows:
USalignWin64.zip em https://zhanggroup.org/US-align/bin/module/USalignWin64.zip
(o pacote Win64 observado nesta maquina reporta Version 20241108; o historico
do site pode listar datas posteriores no codigo-fonte).

HelixScope nao descarrega o binario em runtime. Sem executavel allowlisted:
NOT_INSTALLED. DETECTED nao e LIVE_VALIDATED. A API RCSB permanece um backend
separado; US-align nunca a substitui em silencio.

Nenhuma funcao aqui importa Streamlit. Subprocesso usa lista de argumentos
com shell remaining False.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
import subprocess
import tempfile
from typing import Any, Mapping, Optional, Sequence

from . import engine_validation, provenance, tool_detection, tool_paths

USALIGN_NAMES: tuple[str, ...] = (
    "USalign",
    "USalign.exe",
    "usalign",
    "usalign.exe",
)
USALIGN_ENV: str = "HELIXSCOPE_USALIGN"
OFFICIAL_SOURCE: str = "https://zhanggroup.org/US-align/"
OFFICIAL_WIN64_ZIP: str = "https://zhanggroup.org/US-align/bin/module/USalignWin64.zip"
CITATION: str = (
    "Zhang, Shine, Pyle, Zhang. Nature Methods 19:1109-1115 (2022). "
    "Zhang, Freddolino, Zhang. Nat Protoc 21:517-541 (2025)."
)
VERSION_OBSERVED_IN_RESEARCH: str = "20260527"
VERSION_RE: re.Pattern[str] = re.compile(r"US-align\s*\(Version\s+([0-9]+)\)", re.I)
ALIGNED_RE: re.Pattern[str] = re.compile(
    r"Aligned length=\s*(?P<n>\d+)\s*,\s*RMSD=\s*(?P<rmsd>[-+0-9.]+)\s*,"
    r"\s*Seq_ID=n_identical/n_aligned=\s*(?P<seqid>[-+0-9.]+)",
    re.I,
)
TM_RE: re.Pattern[str] = re.compile(
    r"TM-score=\s*(?P<tm>[-+0-9.]+)\s*\(normalized by length of Structure_(?P<which>[12])",
    re.I,
)
LENGTH_RE: re.Pattern[str] = re.compile(
    r"Length of Structure_(?P<which>[12]):\s*(?P<n>\d+)\s*residues",
    re.I,
)
NAME_RE: re.Pattern[str] = re.compile(
    r"Name of Structure_(?P<which>[12]):\s*(?P<name>.+)$",
    re.I,
)
MATRIX_ROW_RE: re.Pattern[str] = re.compile(
    r"^\s*([012])\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s*$"
)
LOCAL_TIMEOUT_S: float = 60.0
MAX_OUTPUT_CHARS: int = 200_000
ALLOWED_MOL: frozenset[str] = frozenset({"auto", "prot", "RNA"})
# PDB auth_asym_id is one character. mmCIF label_asym_id in deposited
# entries used here is a short ASCII token, not a command-line option.
_CHAIN_ID_RE: re.Pattern[str] = re.compile(r"[A-Za-z0-9]{1,4}\Z")
_STRUCTURE_NAME_RE: re.Pattern[str] = re.compile(
    r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}\.(?:cif|pdb|ent)\Z",
    re.IGNORECASE,
)


class USAlignError(RuntimeError):
    """Falha classificada do US-align.

    Attributes:
        category: TOOL_NOT_INSTALLED, INVALID_INPUT, TIMEOUT, TOOL_FAILED,
            PARSING_ERROR, RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "TOOL_FAILED")


def detect_usalign(which_fn: Optional[tool_detection.WhichFn] = None) -> dict:
    """Procura USalign allowlisted e le a versao do banner.

    Args:
        which_fn: shutil.which injetavel.

    Returns:
        Dict available, status, version, reason, source.

    Raises:
        Nenhum.
    """
    executable = tool_detection.resolve_allowlisted_executable(
        USALIGN_NAMES,
        extra_file_candidates=tool_paths.candidates_for(USALIGN_NAMES),
        env_var=USALIGN_ENV,
        which_fn=which_fn,
    )
    if not executable:
        return {
            "available": False,
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "path": "",
            "version": "",
            "source": OFFICIAL_SOURCE,
            "citation": CITATION,
            "research_version_note": VERSION_OBSERVED_IN_RESEARCH,
            "reason": (
                "US-align is not installed. HelixScope does not download unofficial "
                f"binaries. Official source: {OFFICIAL_SOURCE}. Nucleic-acid and "
                "complex structural comparison therefore remain UNAVAILABLE."
            ),
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
    version = _probe_version(executable)
    return {
        "available": True,
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
        "path": tool_detection.sanitize_tool_path(executable),
        "version": version,
        "source": OFFICIAL_SOURCE,
        "citation": CITATION,
        "reason": "US-align executable detected. DETECTED is not LIVE_VALIDATED.",
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def _probe_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = f"{completed.stdout or ''}\n{completed.stderr or ''}"
    match = VERSION_RE.search(text)
    if match:
        return match.group(1)
    return ""


def nucleic_structural_comparison_availability(which_fn: Optional[tool_detection.WhichFn] = None) -> dict:
    """DNA/RNA estrutural so com engine que suporte acidos nucleicos.

    Args:
        which_fn: which injetavel.

    Returns:
        Dict available False e status UNAVAILABLE se US-align ausente.
        A RCSB Alignment API desta versao e documentada para proteinas
        (C-alpha, minimo 10 residuos).

    Raises:
        Nenhum.
    """
    detected = detect_usalign(which_fn=which_fn)
    live = engine_validation.live_record("US-align")
    nucleic_ok = bool(
        isinstance(live, dict)
        and live.get("ok")
        and str((live.get("details") or {}).get("nucleic_validated") or "") == "true"
    )
    if not detected.get("available"):
        return {
            "available": False,
            "status": "UNAVAILABLE",
            "engine": "US-align",
            "reason": detected.get("reason"),
            "rcsb_alignment_api": (
                "The RCSB Alignment API fits C-alpha atoms and requires at least "
                "10 CA residues. It is not used here as a DNA/RNA aligner."
            ),
        }
    if nucleic_ok:
        return {
            "available": True,
            "status": "LIVE_VALIDATED",
            "engine": "US-align",
            "version": detected.get("version") or (live or {}).get("version"),
            "reason": (
                "US-align nucleic-acid alignment was live-validated on real coordinate "
                "fixtures. RCSB Alignment API remains a separate protein backend."
            ),
        }
    return {
        "available": True,
        "status": "DETECTED",
        "engine": "US-align",
        "version": detected.get("version"),
        "reason": (
            "US-align is DETECTED but nucleic comparison is not LIVE_VALIDATED until "
            "a real nucleic-acid pair is run and the mapping is checked."
        ),
    }


def parse_usalign_stdout(text: str) -> dict:
    """Extrai TM-score, RMSD, identidade e alinhamento impresso.

    Args:
        text: stdout/stderr combinados do executavel.

    Returns:
        Dict com campos numericos ou None quando ausentes (nunca 0 inventado).

    Raises:
        USAlignError: PARSING_ERROR se o bloco essencial estiver ausente.
    """
    blob = str(text or "")
    if len(blob) > MAX_OUTPUT_CHARS:
        raise USAlignError("US-align output exceeds the size cap.", "RESOURCE_LIMIT")
    version = ""
    match_v = VERSION_RE.search(blob)
    if match_v:
        version = match_v.group(1)
    aligned = ALIGNED_RE.search(blob)
    if not aligned:
        raise USAlignError(
            "US-align output did not contain Aligned length / RMSD / Seq_ID.",
            "PARSING_ERROR",
        )
    n_aligned = int(aligned.group("n"))
    rmsd = _finite_float(aligned.group("rmsd"))
    seq_id = _finite_float(aligned.group("seqid"))
    tm_scores: list[dict] = []
    for match in TM_RE.finditer(blob):
        tm_scores.append(
            {
                "type": f"TM-score_structure_{match.group('which')}",
                "value": _finite_float(match.group("tm")),
                "unit": "dimensionless_0_1",
                "source": "US-align stdout",
                "normalized_by": f"Structure_{match.group('which')}",
            }
        )
    lengths: dict[str, int] = {}
    for match in LENGTH_RE.finditer(blob):
        lengths[match.group("which")] = int(match.group("n"))
    names: dict[str, str] = {}
    for match in NAME_RE.finditer(blob):
        names[match.group("which")] = match.group("name").strip()[:240]
    seq1, seq2, colon = _parse_alignment_block(blob)
    return {
        "version": version,
        "n_aligned": n_aligned,
        "rmsd_angstrom": rmsd,
        "sequence_identity": seq_id,
        "tm_scores": tm_scores,
        "length_1": lengths.get("1"),
        "length_2": lengths.get("2"),
        "name_1": names.get("1", ""),
        "name_2": names.get("2", ""),
        "aligned_seq1": seq1,
        "aligned_seq2": seq2,
        "alignment_mark": colon,
        "raw_excerpt": blob[:4000],
    }


def _parse_alignment_block(blob: str) -> tuple[str, str, str]:
    lines = [line.rstrip("\n") for line in blob.splitlines()]
    for index, line in enumerate(lines):
        if "denotes residue pairs" in line.lower():
            seq_lines = [item for item in lines[index + 1 :] if item.strip() and not item.startswith("#")]
            if len(seq_lines) >= 3:
                return seq_lines[0].strip(), seq_lines[2].strip(), seq_lines[1]
    return "", "", ""


def parse_usalign_matrix(text: str) -> Optional[list[float]]:
    """Le a matriz 3x3 + translacao do ficheiro -m para 4x4 column-major.

    Args:
        text: Conteudo do ficheiro escrito por US-align -m.

    Returns:
        16 floats column-major, ou None se o parse falhar. Nao inventa identidade.

    Raises:
        Nenhum.
    """
    rows: dict[int, tuple[float, float, float, float]] = {}
    for line in str(text or "").splitlines():
        match = MATRIX_ROW_RE.match(line)
        if not match:
            continue
        index = int(match.group(1))
        t_i = _finite_float(match.group(2))
        u0 = _finite_float(match.group(3))
        u1 = _finite_float(match.group(4))
        u2 = _finite_float(match.group(5))
        if None in {t_i, u0, u1, u2}:
            return None
        rows[index] = (float(t_i), float(u0), float(u1), float(u2))
    if set(rows) != {0, 1, 2}:
        return None
    matrix = [0.0] * 16
    for i in range(3):
        t_i, u0, u1, u2 = rows[i]
        matrix[i] = u0
        matrix[4 + i] = u1
        matrix[8 + i] = u2
        matrix[12 + i] = t_i
    matrix[15] = 1.0
    return matrix


def _finite_float(value: object) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def rotation_translation_to_column_major(u: Sequence[Sequence[float]], t: Sequence[float]) -> list[float]:
    """Converte u 3x3 e t 3 para 4x4 column-major (convencao RCSB/HelixScope).

    Args:
        u: Matriz de rotacao linha-maior como no stdout US-align.
        t: Translacao.

    Returns:
        16 floats.

    Raises:
        USAlignError: INVALID_INPUT se as dimensoes falharem.
    """
    if len(u) != 3 or any(len(row) != 3 for row in u) or len(t) != 3:
        raise USAlignError("US-align rotation/translation has the wrong shape.", "PARSING_ERROR")
    matrix = [0.0] * 16
    for i in range(3):
        matrix[i] = float(u[i][0])
        matrix[4 + i] = float(u[i][1])
        matrix[8 + i] = float(u[i][2])
        matrix[12 + i] = float(t[i])
    matrix[15] = 1.0
    return matrix


def normalize_chain_id(chain_id: str) -> str:
    """Return a US-align chain token, or empty when no chain filter is requested.

    Args:
        chain_id: Caller-supplied chain label. Empty means the whole structure.

    Returns:
        The same ASCII letters or digits, length 1 to 4.

    Raises:
        USAlignError: INVALID_INPUT when the label is not a structure chain id.

    Nota biologica:
        O token e o auth_asym_id ou label_asym_id da estrutura, nao uma opcao
        do US-align. Um hifen inicial seria lido como outra flag.
    """
    if chain_id is None:
        return ""
    if not isinstance(chain_id, str):
        raise USAlignError("US-align chain id must be text.", "INVALID_INPUT")
    if chain_id == "":
        return ""
    if chain_id != chain_id.strip() or any(ord(char) < 32 or ord(char) == 127 for char in chain_id):
        raise USAlignError("US-align chain id contains whitespace or controls.", "INVALID_INPUT")
    matched = _CHAIN_ID_RE.fullmatch(chain_id)
    if matched is None:
        raise USAlignError(
            "US-align chain id must be 1 to 4 ASCII letters or digits.",
            "INVALID_INPUT",
        )
    return matched.group(0)


def _existing_structure_file(path: str) -> str:
    """Resolve a caller path to an existing structure file.

    Args:
        path: PDB or mmCIF path.

    Returns:
        Real path of an existing .cif, .pdb, or .ent file.

    Raises:
        USAlignError: INVALID_INPUT when the path is not that file.
    """
    if not isinstance(path, str) or path == "" or "\x00" in path:
        raise USAlignError("US-align structure path is empty or contains NUL.", "INVALID_INPUT")
    if any(ord(char) < 32 for char in path):
        raise USAlignError("US-align structure path contains a control character.", "INVALID_INPUT")
    resolved = os.path.realpath(path)
    if not os.path.isfile(resolved):
        raise USAlignError("US-align input structure file is missing.", "INVALID_INPUT")
    if not resolved.lower().endswith((".cif", ".pdb", ".ent")):
        raise USAlignError("US-align input must be .cif, .pdb, or .ent.", "INVALID_INPUT")
    return resolved


def _stage_structure(source: str, directory: str, which: str, chain_id: str) -> str:
    """Copy a structure into the work directory under a fixed name.

    Args:
        source: Real path of the caller file. It is read, not passed through.
        directory: Temporary directory created by this module.
        which: "1" or "2".
        chain_id: Empty keeps every chain. Otherwise only that chain is written.

    Returns:
        Path of the staged file. The name is chosen here, not by the caller.

    Raises:
        USAlignError: INVALID_INPUT when the chain is absent or the file
            cannot be reduced to that chain.

    Nota biologica:
        O US-align recebe o arquivo estagiado. O identificador de cadeia nao
        entra na linha de comando. A selecao, quando pedida, acontece antes,
        no conteudo da estrutura.
    """
    if which == "1":
        name = "structure1.pdb" if source.lower().endswith((".pdb", ".ent")) else "structure1.cif"
    elif which == "2":
        name = "structure2.pdb" if source.lower().endswith((".pdb", ".ent")) else "structure2.cif"
    else:
        raise USAlignError("US-align stage index is not 1 or 2.", "INVALID_INPUT")
    dest = os.path.join(directory, name)
    if chain_id == "":
        shutil.copyfile(source, dest)
        return dest
    _write_single_chain(source, dest, chain_id)
    return dest


def _write_single_chain(source: str, dest: str, chain_id: str) -> None:
    """Write one chain from a PDB or mmCIF file.

    Args:
        source: Input structure path.
        dest: Output path inside the temporary directory.
        chain_id: ASCII chain token already validated.

    Returns:
        None.

    Raises:
        USAlignError: INVALID_INPUT when that chain is not in the file.
    """
    from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser, Select

    class _OnlyChain(Select):
        def accept_chain(self, chain: object) -> bool:
            return str(getattr(chain, "id", "")) == chain_id

    pdb = source.lower().endswith((".pdb", ".ent"))
    try:
        if pdb:
            structure = PDBParser(QUIET=True).get_structure("helixscope", source)
        else:
            structure = MMCIFParser(QUIET=True).get_structure("helixscope", source)
    except Exception as exc:
        raise USAlignError(
            "US-align could not read the structure while selecting a chain.",
            "INVALID_INPUT",
        ) from exc
    present = {
        str(chain.id)
        for model in structure
        for chain in model
    }
    if chain_id not in present:
        raise USAlignError(
            f"Chain {chain_id} is not in the structure.",
            "INVALID_INPUT",
        )
    if pdb:
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(dest, select=_OnlyChain())
    else:
        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(dest, select=_OnlyChain())


def align_structure_files(
    path1: str,
    path2: str,
    *,
    mol: str = "auto",
    chain1: str = "",
    chain2: str = "",
    mm: int = 0,
    timeout_s: float = LOCAL_TIMEOUT_S,
    record_validation: bool = False,
) -> dict:
    """Executa USalign sobre dois ficheiros PDB/mmCIF locais.

    Args:
        path1: Estrutura 1 (US-align roda esta sobre a estrutura 2).
        path2: Estrutura 2.
        mol: auto, prot ou RNA (RNA inclui DNA segundo o help upstream).
        chain1: Cadeia opcional (-chain1).
        chain2: Cadeia opcional (-chain2).
        mm: Modo multimeric. 0 = monomeros (suportado nesta interface).
        timeout_s: Timeout do subprocesso.
        record_validation: Se True e o parse passar, grava LIVE_VALIDATED.

    Returns:
        Registro de comparacao estrutural com engine US-align.

    Raises:
        USAlignError.

    Nota biologica:
        TM-score e RMSD sao os valores impressos pelo US-align. HelixScope nao
        os renormaliza. A API RCSB nao e consultada neste caminho.
    """
    molecule = str(mol or "auto").strip()
    if molecule not in ALLOWED_MOL:
        raise USAlignError(f"US-align -mol '{mol}' is not offered.", "INVALID_INPUT")
    if int(mm) != 0:
        raise USAlignError(
            "HelixScope US-align workflow is monomer alignment (-mm 0). "
            "Oligomer/complex modes are not exposed until a dedicated interface is validated.",
            "INVALID_INPUT",
        )
    safe_chain1 = normalize_chain_id(chain1)
    safe_chain2 = normalize_chain_id(chain2)
    if molecule == "auto":
        mol_flag = "auto"
    elif molecule == "prot":
        mol_flag = "prot"
    elif molecule == "RNA":
        mol_flag = "RNA"
    else:
        raise USAlignError(f"US-align -mol '{mol}' is not offered.", "INVALID_INPUT")
    source1 = _existing_structure_file(path1)
    source2 = _existing_structure_file(path2)
    detected = detect_usalign()
    if not detected.get("available"):
        raise USAlignError(str(detected.get("reason") or "US-align is not installed."), "TOOL_NOT_INSTALLED")
    executable = tool_detection.resolve_allowlisted_executable(
        USALIGN_NAMES,
        extra_file_candidates=tool_paths.candidates_for(USALIGN_NAMES),
        env_var=USALIGN_ENV,
    )
    if not executable:
        raise USAlignError("US-align executable UNAVAILABLE.", "TOOL_NOT_INSTALLED")
    try:
        timeout = float(timeout_s)
    except (TypeError, ValueError) as exc:
        raise USAlignError("US-align timeout must be numeric.", "INVALID_INPUT") from exc
    if not math.isfinite(timeout) or timeout <= 0:
        raise USAlignError("US-align timeout must be positive.", "INVALID_INPUT")
    with tempfile.TemporaryDirectory(prefix="helixscope_usalign_") as tmp:
        file1 = _stage_structure(source1, tmp, "1", safe_chain1)
        file2 = _stage_structure(source2, tmp, "2", safe_chain2)
        matrix_path = os.path.join(tmp, "matrix.txt")
        argv = [executable, file1, file2, "-ter", "2", "-mol", mol_flag, "-m", matrix_path]
        try:
            completed = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
                shell=False,
                cwd=tmp,
            )
        except subprocess.TimeoutExpired as exc:
            raise USAlignError(
                f"US-align exceeded {timeout:.0f}s and was stopped.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise USAlignError(f"US-align could not be executed: {exc}", "TOOL_FAILED") from exc
        stdout = str(completed.stdout or "")
        stderr = str(completed.stderr or "")
        blob = stdout + "\n" + stderr
        if len(blob) > MAX_OUTPUT_CHARS:
            raise USAlignError("US-align output exceeds the size cap.", "RESOURCE_LIMIT")
        parsed = parse_usalign_stdout(blob)
        matrix = None
        if os.path.isfile(matrix_path):
            matrix = parse_usalign_matrix(_read_capped(matrix_path))
        record = build_alignment_record(
            parsed,
            path1=source1,
            path2=source2,
            mol=molecule,
            chain1=safe_chain1,
            chain2=safe_chain2,
            matrix=matrix,
            executable_version=str(parsed.get("version") or detected.get("version") or ""),
        )
        if record_validation:
            nucleic = molecule == "RNA" or _looks_nucleic(parsed)
            previous = engine_validation.live_record("US-align") or {}
            details = dict(previous.get("details") or {})
            details.update(
                {
                    "fixture_1": os.path.basename(source1),
                    "fixture_2": os.path.basename(source2),
                    "n_aligned": record.get("n_aligned_residue_pairs"),
                    "rmsd_angstrom": record.get("rmsd_global_angstrom"),
                    "mol": molecule,
                }
            )
            if nucleic:
                details["nucleic_validated"] = "true"
            else:
                details.setdefault("nucleic_validated", "false")
            engine_validation.record_live_validation(
                "US-align",
                ok=True,
                version=str(record.get("engine_version") or ""),
                details=details,
            )
        return record


def _read_capped(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read(MAX_OUTPUT_CHARS)


def _looks_nucleic(parsed: Mapping[str, Any]) -> bool:
    seq = str(parsed.get("aligned_seq1") or "").upper()
    letters = {char for char in seq if char.isalpha()}
    if not letters:
        return False
    return letters <= set("ACGTUacgtu")


def build_alignment_record(
    parsed: Mapping[str, Any],
    *,
    path1: str,
    path2: str,
    mol: str,
    chain1: str,
    chain2: str,
    matrix: Optional[Sequence[float]],
    executable_version: str,
) -> dict:
    """Monta o objeto de comparacao no contrato usado pela UI Compare.

    Args:
        parsed: Saida de parse_usalign_stdout.
        path1, path2: Ficheiros de entrada.
        mol: Flag -mol.
        chain1, chain2: Cadeias pedidas.
        matrix: 4x4 column-major da estrutura 1, ou None.
        executable_version: Banner US-align.

    Returns:
        Dict alignment-like. RMSD global e de bloco coincidem no monomero -mm 0
        porque ha um unico bloco; os campos permanecem separados.

    Raises:
        Nenhum.
    """
    n_pairs = parsed.get("n_aligned")
    rmsd = parsed.get("rmsd_angstrom")
    tm_scores = list(parsed.get("tm_scores") or [])
    identity = parsed.get("sequence_identity")
    atoms = "CA" if mol == "prot" else ("C3'" if mol == "RNA" else "CA or C3' (US-align auto)")
    identity_hash = hashlib.sha256(
        f"{os.path.basename(path1)}|{os.path.basename(path2)}|{mol}|{n_pairs}|{rmsd}".encode("utf-8")
    ).hexdigest()
    pairs = _pairs_from_alignment(
        str(parsed.get("aligned_seq1") or ""),
        str(parsed.get("aligned_seq2") or ""),
        chain1 or "A",
        chain2 or "A",
        atoms,
    )
    transform = list(matrix) if matrix is not None and len(matrix) == 16 else []
    identity_4 = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    block = {
        "block_index": 0,
        "n_pairs": len(pairs) if pairs else n_pairs,
        "pairs": pairs,
        "reference_transform_column_major_4x4": transform or identity_4,
        "target_transform_column_major_4x4": identity_4,
        "transform_convention": (
            "US-align rotates Structure_1 onto Structure_2. The 4x4 is stored on "
            "the reference (structure 1). HelixScope does not apply this matrix to "
            "the RCSB target-transform overlay without a dedicated US-align viewer path."
        ),
        "rmsd_angstrom": rmsd,
        "rmsd_scope": "US-align aligned residues of the monomer comparison",
    }
    len1 = parsed.get("length_1")
    len2 = parsed.get("length_2")
    coverage: list[float] = []
    if isinstance(n_pairs, int) and isinstance(len1, int) and len1 > 0:
        coverage.append(round(100.0 * n_pairs / len1, 1))
    if isinstance(n_pairs, int) and isinstance(len2, int) and len2 > 0:
        coverage.append(round(100.0 * n_pairs / len2, 1))
    record = {
        "kind": "structural_comparison",
        "status": "COMPUTED",
        "evidence_status": "COMPUTED",
        "source": "US-align local executable",
        "docs": OFFICIAL_SOURCE,
        "engine": "US-align",
        "engine_location": "local",
        "engine_version": executable_version,
        "method": "US-align",
        "method_display": f"US-align -mol {mol} -mm 0",
        "method_kind": "rigid",
        "alignment_mode": "pairwise_monomer",
        "atoms_fitted": atoms,
        "atoms_fitted_note": (
            "US-align default representation is CA for protein and C3' for RNA/DNA "
            "when -mol RNA. auto may mix polymer types present in the files."
        ),
        "reference": {
            "entry_id": os.path.splitext(os.path.basename(path1))[0],
            "path": os.path.basename(path1),
            "asym_id": chain1 or "",
            "kind_label": "EXPERIMENTAL",
            "n_residues": len1,
        },
        "target": {
            "entry_id": os.path.splitext(os.path.basename(path2))[0],
            "path": os.path.basename(path2),
            "asym_id": chain2 or "",
            "kind_label": "EXPERIMENTAL",
            "n_residues": len2,
        },
        "n_blocks": 1,
        "blocks": [block],
        "residue_pairs": pairs,
        "n_aligned_residue_pairs": n_pairs,
        "mapping_policy": (
            "Pairs are 1-based indices along the ungapped US-align printed chains, "
            "not auth_seq_id. Residues absent from the alignment are UNMAPPED."
        ),
        "global_scores": [
            {"type": "RMSD", "value": rmsd, "unit": "angstrom", "source": "US-align"},
            {"type": "sequence-identity", "value": identity, "unit": "fraction_0_1", "source": "US-align"},
            *tm_scores,
        ],
        "rmsd_global_angstrom": rmsd,
        "rmsd_block0_angstrom": rmsd,
        "rmsd_note": (
            "For monomer -mm 0, US-align prints one RMSD. HelixScope stores it as "
            "both global and block-0 fields without averaging additional blocks."
        ),
        "tm_scores": tm_scores,
        "tm_score_note": (
            "Two TM-scores are reported (normalized by each structure length). "
            "They are not fused. Absence is Unavailable, not 0."
        ),
        "sequence_identity": identity,
        "aln_coverage_percent": coverage,
        "aligned_seq1": parsed.get("aligned_seq1"),
        "aligned_seq2": parsed.get("aligned_seq2"),
        "superposition_visual": "UNAVAILABLE",
        "superposition_visual_note": (
            "US-align rotates structure 1 onto structure 2. The Compare 3D overlay "
            "path is built for RCSB target-transform matrices. Metrics above are "
            "from US-align; the overlay is not silently reused from RCSB."
        ),
        "identity_hash": identity_hash,
        "citation": CITATION,
        "disclaimer": (
            "LOCAL · US-align. This backend is not the RCSB Alignment API. "
            "TM-score >= 0.5 (proteins) or 0.45 (RNA) is the upstream topology "
            "guideline, not a proof of equivalent function."
        ),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
    return record


def _pairs_from_alignment(
    seq1: str,
    seq2: str,
    chain1: str,
    chain2: str,
    atoms: str,
) -> list[dict]:
    if not seq1 or not seq2 or len(seq1) != len(seq2):
        return []
    pairs: list[dict] = []
    i = 0
    j = 0
    for a, b in zip(seq1, seq2):
        gap_a = a in "-."
        gap_b = b in "-."
        if not gap_a:
            i += 1
        if not gap_b:
            j += 1
        if gap_a or gap_b:
            continue
        pairs.append(
            {
                "reference_asym_id": chain1,
                "target_asym_id": chain2,
                "reference_label_seq_id": i,
                "target_label_seq_id": j,
                "numbering": "1-based index along US-align printed ungapped chain, not auth_seq_id",
                "atoms": atoms,
            }
        )
    return pairs


def write_structure_text(text: str, directory: str, filename: str) -> str:
    """Grava mmCIF/PDB num diretorio temporario controlado.

    Args:
        text: Coordenadas ja carregadas.
        directory: Pasta temporaria do caller.
        filename: Basename allowlisted (.cif/.pdb).

    Returns:
        Caminho absoluto escrito.

    Raises:
        USAlignError: INVALID_INPUT se o texto ou o nome forem invalidos.
    """
    body = str(text or "")
    if not body.strip() or "\x00" in body:
        raise USAlignError("US-align cannot run without structure text.", "INVALID_INPUT")
    if not isinstance(directory, str) or directory == "" or "\x00" in directory:
        raise USAlignError("US-align output directory is not usable.", "INVALID_INPUT")
    base = os.path.basename(str(filename or ""))
    matched = _STRUCTURE_NAME_RE.fullmatch(base)
    if matched is None or ".." in base:
        raise USAlignError("US-align output filename is not a basename.", "INVALID_INPUT")
    base = matched.group(0)
    root = os.path.realpath(directory)
    path = os.path.realpath(os.path.join(root, base))
    try:
        shared = os.path.commonpath([root, path])
    except ValueError as exc:
        raise USAlignError("US-align output path leaves the temporary directory.", "INVALID_INPUT") from exc
    if os.path.normcase(shared) != os.path.normcase(root):
        raise USAlignError("US-align output path leaves the temporary directory.", "INVALID_INPUT")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(body if body.endswith("\n") else body + "\n")
    return path
