"""Predicao de estrutura secundaria de RNA via ViennaRNA real.

Backends: modulo Python `RNA` (RNAlib) ou executavel allowlisted `RNAfold`.
Nenhum algoritmo de folding e reimplementado neste modulo. Sem a ferramenta,
o estado e UNAVAILABLE; falha de execucao nao vira estrutura ilustrativa.

Nenhuma funcao aqui importa Streamlit. Subprocesso usa lista de argumentos,
nunca shell=True. Arquivos temporarios sao do sistema e sempre removidos.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from . import dna_analysis, engine_validation, provenance, scientific_checks

MAX_FOLD_NT: int = 600
"""Teto de nucleotidos submetidos ao folding global (protecao O(n^3), nao regra biologica)."""

MAX_LOCAL_RNA_NT: int = 150_000
"""Teto da analise local RNAplfold. Nao autoriza MFE global nessa escala."""

LOCAL_WINDOW_DEFAULT: int = 70
"""Janela default documentada de RNAplfold (-W)."""

LOCAL_ANALYSIS_TIMEOUT_S: float = 120.0
"""Timeout da analise local. Estouro e TIMEOUT, nao estrutura vazia."""

LOCAL_PAIR_CUTOFF: float = 0.01
"""Pares abaixo deste valor sao omitidos pelo proprio RNAplfold, nao por corte da sequencia."""

LOCAL_PREVIEW_PAIRS: int = 20
"""Pares mostrados no preview. O total e calculado sobre a lista completa."""

MIN_FOLD_NT: int = 1
"""Comprimento minimo aceito; ViennaRNA pode devolver so pontos."""

LOCAL_TIMEOUT_S: float = 20.0
"""Timeout do folding (subprocesso ou chamada RNAlib em thread)."""

MAX_OUTPUT_BYTES: int = 100_000
"""Teto do texto devolvido por RNAfold."""

VIENNA_DEFAULT_TEMPERATURE_C: float = 37.0
"""Default documentado de RNAfold --temp."""

VIENNA_DEFAULT_SALT_M: float = 1.021
"""Default documentado de RNAfold --salt."""

RNAFOLD_NAMES: Tuple[str, ...] = ("RNAfold", "RNAfold.exe", "rnafold", "rnafold.exe")
"""Basenames allowlisted. Nenhum path arbitrario do usuario."""

CANONICAL_FOLD_PAIRS: frozenset[Tuple[str, str]] = frozenset(
    {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
)
"""Pares permitidos no modelo Turner usado por ViennaRNA (Watson-Crick + wobble)."""

DOT_BRACKET_CHARS: frozenset[str] = frozenset(".()")
"""Alfabeto MFE sem pseudoknots. Colchetes nao sao aceitos neste pipeline."""

ENERGY_TAIL: re.Pattern[str] = re.compile(
    r"^\s*(?P<structure>[.()]+)\s*(?:\(\s*(?P<energy>[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*(?:kcal/mol)?\))?\s*$"
)
"""Linha de estrutura RNAfold: notacao de colchetes e MFE opcional em kcal/mol."""


class FoldingError(Exception):
    """Falha classificada de folding. Nunca deve ser reescrita como estrutura falsa.

    Attributes:
        category: TOOL_NOT_INSTALLED, TOOL_FAILED, TIMEOUT, RESOURCE_LIMIT,
            PARSING_ERROR, INVALID_INPUT, JOB_FAILED.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "TOOL_FAILED")


def tool_availability() -> dict:
    """Detecta ViennaRNA real (RNAlib Python ou RNAfold) sem fingir versao.

    Args:
        Nenhum.

    Returns:
        Dict com available, backends, version (vazia se nao detectada) e reason.

    Raises:
        Nenhum.

    Nota biologica:
        A predicao MFE de ViennaRNA e o algoritmo de programacao dinamica de
        Zuker com parametros de Turner. Sem RNAlib/RNAfold, nao ha predicao.
    """
    python = detect_python_bindings()
    binary = detect_rnafold_executable()
    backends = [python, binary]
    any_available = bool(python["available"] or binary["available"])
    if any_available:
        reason = (
            "RNA secondary structure uses ViennaRNA (RNAlib and/or RNAfold). "
            "The result is a predicted MFE structure, not an experimental structure."
        )
    else:
        reason = (
            "RNA folding unavailable in this environment. ViennaRNA RNAlib "
            "(Python module RNA) and the RNAfold executable were not detected. "
            "HelixScope will not draw an invented structure."
        )
    return {
        "available": any_available,
        "any_available": any_available,
        "backends": backends,
        "python_rna": python,
        "rnafold": binary,
        "reason": reason,
        "recommended_tools": [
            "ViennaRNA RNAfold / RNAlib (Lorenz et al. 2011)",
            "mfold / UNAfold",
            "RNAstructure",
        ],
        "summary_reason": reason,
        "kind": "predicted",
        "experimental": False,
        "partition_function": False,
        "pseudoknots": False,
    }


def detect_python_bindings() -> dict:
    """Tenta importar o modulo RNA oficial do ViennaRNA.

    Args:
        Nenhum.

    Returns:
        Dict available, version (detectada ou vazia), method.

    Raises:
        Nenhum.
    """
    try:
        import RNA  # type: ignore
    except ImportError:
        return {
            "id": "viennarna_python",
            "tool": "ViennaRNA RNAlib",
            "method": "Python module RNA",
            "available": False,
            "version": "",
            "version_status": "unavailable",
            "reason": "Python module RNA is not installed.",
        }
    version = _module_version(RNA)
    return {
        "id": "viennarna_python",
        "tool": "ViennaRNA RNAlib",
        "method": "RNA.fold (MFE)",
        "available": True,
        "version": version,
        "version_status": "detected" if version else "unavailable",
        "reason": "ViennaRNA Python bindings imported.",
        "module": "RNA",
    }


def detect_rnafold_executable() -> dict:
    """Procura apenas basenames allowlisted via shutil.which.

    Args:
        Nenhum.

    Returns:
        Dict available, version, path.

    Raises:
        Nenhum.
    """
    executable = _resolve_allowlisted_executable(RNAFOLD_NAMES)
    if executable is None:
        return {
            "id": "rnafold_exe",
            "tool": "ViennaRNA RNAfold",
            "method": "RNAfold executable",
            "available": False,
            "version": "",
            "version_status": "unavailable",
            "reason": "RNAfold executable not installed.",
            "path": "",
        }
    version = _probe_rnafold_version(executable)
    return {
        "id": "rnafold_exe",
        "tool": "ViennaRNA RNAfold",
        "method": "RNAfold --noPS --noconv",
        "available": True,
        "version": version,
        "version_status": "detected" if version else "unavailable",
        "reason": f"Found allowlisted executable {os.path.basename(executable)}.",
        "path": executable,
    }


def available_backend_ids() -> List[str]:
    """IDs de backend que realmente podem correr.

    Args:
        Nenhum.

    Returns:
        Lista possivelmente vazia.

    Raises:
        Nenhum.
    """
    ids: List[str] = []
    if detect_python_bindings()["available"]:
        ids.append("viennarna_python")
    if detect_rnafold_executable()["available"]:
        ids.append("rnafold_exe")
    return ids


def default_fold_parameters(region_mode: str = "full sequence") -> dict:
    """Parametros efetivos gravados na proveniencia (defaults documentados do RNAfold).

    Args:
        region_mode: full sequence ou selected subsequence.

    Returns:
        Dict serializavel.

    Raises:
        Nenhum.
    """
    return _default_parameters(region_mode)


def validate_rna_for_folding(sequence: str) -> dict:
    """Valida RNA canonico ACGU. Nao converte T em U.

    Args:
        sequence: Texto bruto.

    Returns:
        Dict com sequence (maiuscula, sem espaco), hash, length, normalization.

    Raises:
        FoldingError: INVALID_INPUT ou RESOURCE_LIMIT.
    """
    if not isinstance(sequence, str):
        raise FoldingError("RNA sequence must be a string.", "INVALID_INPUT")
    stripped = "".join(sequence.split())
    if not stripped:
        raise FoldingError("RNA sequence is empty.", "INVALID_INPUT")
    upper = stripped.upper()
    if any(char == "T" for char in upper):
        raise FoldingError(
            "Sequence contains T. DNA is not sent to RNA folding and T is not "
            "silently converted to U.",
            "INVALID_INPUT",
        )
    info = dna_analysis.validate_for_molecule(upper, "RNA")
    if not info["is_valid"]:
        raise FoldingError(
            str(info.get("rejection_reason") or "Sequence is not canonical RNA (A, C, G, U)."),
            "INVALID_INPUT",
        )
    residues = str(info["sequence"])
    ambiguous = sorted({char for char in residues if char not in dna_analysis.CANONICAL_RNA_ALPHABET})
    if ambiguous:
        raise FoldingError(
            "Ambiguous RNA symbols are not substituted. Rejected: "
            + ", ".join(ambiguous)
            + ". ViennaRNA is only invoked on A, C, G and U.",
            "INVALID_INPUT",
        )
    if len(residues) < MIN_FOLD_NT:
        raise FoldingError("RNA is shorter than the folding minimum.", "INVALID_INPUT")
    if len(residues) > MAX_FOLD_NT:
        raise FoldingError(
            f"RNA folding is limited to {MAX_FOLD_NT} nt (RESOURCE_LIMIT, not "
            "an empty structure).",
            "RESOURCE_LIMIT",
        )
    return {
        "sequence": residues,
        "length": len(residues),
        "hash": provenance.sequence_digest(residues),
        "normalization": {
            "whitespace_removed": "".join(sequence.split()) != sequence,
            "uppercased": any(char.islower() for char in stripped),
            "nucleotides_changed": False,
        },
        "molecule": "RNA",
    }


def parse_rnafold_output(body: str, expected_sequence: str) -> dict:
    """Interpreta stdout de RNAfold. Nao corrige sequencia nem energia.

    Args:
        body: Texto da ferramenta.
        expected_sequence: RNA submetido, ja validado.

    Returns:
        Dict com structure e mfe_kcal_mol.

    Raises:
        FoldingError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    if body is None:
        raise FoldingError("RNAfold output is empty.", "PARSING_ERROR")
    raw = body if isinstance(body, str) else body.decode("utf-8", errors="replace")
    if len(raw.encode("utf-8")) > MAX_OUTPUT_BYTES:
        raise FoldingError(
            f"RNAfold output exceeds {MAX_OUTPUT_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if not lines:
        raise FoldingError("RNAfold output is empty.", "PARSING_ERROR")
    sequence_line = ""
    structure_line = ""
    for line in lines:
        if line.startswith(">"):
            continue
        if not sequence_line:
            candidate = line.replace(" ", "").upper()
            if re.fullmatch(r"[ACGU]+", candidate):
                sequence_line = candidate
                continue
        if sequence_line and not structure_line:
            structure_line = line
            break
    if not sequence_line or not structure_line:
        raise FoldingError(
            "RNAfold output did not contain a sequence line and a structure line.",
            "PARSING_ERROR",
        )
    if sequence_line != expected_sequence:
        raise FoldingError(
            "RNAfold reported a sequence that does not match the submitted RNA. "
            "The output was not repaired.",
            "PARSING_ERROR",
        )
    parsed = _parse_structure_line(structure_line)
    if len(parsed["structure"]) != len(expected_sequence):
        raise FoldingError(
            "Dot-bracket length does not match the RNA sequence length.",
            "PARSING_ERROR",
        )
    if parsed["mfe_kcal_mol"] is None:
        raise FoldingError(
            "RNAfold did not report a numeric MFE. Energy was not invented.",
            "PARSING_ERROR",
        )
    return parsed


def parse_dot_bracket(structure: str) -> List[Tuple[int, int]]:
    """Extrai pares 0-based da notacao de colchetes MFE.

    Args:
        structure: String de '.', '(' e ')'.

    Returns:
        Lista de (i, j) com i < j.

    Raises:
        FoldingError: PARSING_ERROR se desbalanceada ou com simbolos extra.
    """
    text = str(structure or "")
    if not text:
        raise FoldingError("Dot-bracket structure is empty.", "PARSING_ERROR")
    if any(char not in DOT_BRACKET_CHARS for char in text):
        raise FoldingError(
            "Dot-bracket contains characters other than '.', '(' and ')'. "
            "Pseudoknot symbols are not accepted in this MFE pipeline.",
            "PARSING_ERROR",
        )
    if not scientific_checks.rna_parentheses_are_balanced(text):
        raise FoldingError("Dot-bracket parentheses are not balanced.", "PARSING_ERROR")
    stack: List[int] = []
    pairs: List[Tuple[int, int]] = []
    for index, char in enumerate(text):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                raise FoldingError("Dot-bracket has a closing pair without an opening pair.", "PARSING_ERROR")
            left = stack.pop()
            pairs.append((left, index))
    if stack:
        raise FoldingError("Dot-bracket has unmatched opening pairs.", "PARSING_ERROR")
    pairs.sort()
    return pairs


def reconstruct_dot_bracket(length: int, pairs: Sequence[Tuple[int, int]]) -> str:
    """Reconstroi a notacao de colchetes a partir dos pares 0-based.

    Args:
        length: Comprimento da sequencia.
        pairs: Pares (i, j).

    Returns:
        String de comprimento `length`.

    Raises:
        FoldingError: PARSING_ERROR se os pares sairem do intervalo.
    """
    if length < 0:
        raise FoldingError("Structure length is negative.", "PARSING_ERROR")
    chars = ["."] * length
    for item in pairs:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise FoldingError("A base pair is not an (i, j) tuple.", "PARSING_ERROR")
        left, right = int(item[0]), int(item[1])
        if left < 0 or right <= left or right >= length:
            raise FoldingError("A base pair is outside the sequence coordinates.", "PARSING_ERROR")
        if chars[left] != "." or chars[right] != ".":
            raise FoldingError("Base pairs overlap in the reconstructed structure.", "PARSING_ERROR")
        chars[left] = "("
        chars[right] = ")"
    return "".join(chars)


def validate_pairs_against_sequence(
    sequence: str, pairs: Sequence[Tuple[int, int]]
) -> None:
    """Recusa pares que o modelo Turner padrao nao permite.

    Args:
        sequence: RNA ACGU.
        pairs: Pares 0-based.

    Returns:
        None.

    Raises:
        FoldingError: PARSING_ERROR se um par for impossivel.
    """
    residues = str(sequence or "")
    for left, right in pairs:
        if left < 0 or right >= len(residues) or left >= right:
            raise FoldingError("A parsed pair is outside the RNA sequence.", "PARSING_ERROR")
        bases = (residues[left], residues[right])
        if bases not in CANONICAL_FOLD_PAIRS:
            raise FoldingError(
                f"Pair {left}:{residues[left]}-{right}:{residues[right]} is not a "
                "Watson-Crick or GU wobble pair. The structure was not kept.",
                "PARSING_ERROR",
            )


def classify_structural_elements(
    structure: str, pairs: Sequence[Tuple[int, int]]
) -> List[dict]:
    """Classifica stems, hairpins, bulges e loops internos a partir do dot-bracket.

    Args:
        structure: Notacao validada.
        pairs: Pares 0-based.

    Returns:
        Lista de elementos topologicos. Nao e anotacao funcional.

    Raises:
        Nenhum.

    Nota biologica:
        Hairpin/stem aqui sao classes da notacao de colchetes, nao elementos
        regulatorios nem evidencia experimental.
    """
    length = len(structure)
    pair_map = {left: right for left, right in pairs}
    pair_map.update({right: left for left, right in pairs})
    elements: List[dict] = []
    stacked: List[Tuple[int, int]] = []
    for left, right in sorted(pairs):
        prev = stacked[-1] if stacked else None
        if prev and left == prev[0] + 1 and right == prev[1] - 1:
            stacked.append((left, right))
            continue
        if stacked:
            elements.append(_stem_record(stacked))
        stacked = [(left, right)]
    if stacked:
        elements.append(_stem_record(stacked))
    for left, right in pairs:
        inner = structure[left + 1 : right]
        if inner and all(char == "." for char in inner):
            elements.append(
                {
                    "type": "hairpin_loop",
                    "start": left + 1,
                    "end": right,
                    "length": right - left - 1,
                    "status": "COMPUTED",
                    "label": "hairpin loop (dot-bracket)",
                }
            )
    unpaired_runs = _unpaired_runs(structure)
    paired_set = {left for left, _right in pairs} | {right for _left, right in pairs}
    for start, end in unpaired_runs:
        if start == 0 or end == length:
            elements.append(
                {
                    "type": "exterior_loop",
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "status": "COMPUTED",
                    "label": "exterior unpaired (dot-bracket)",
                }
            )
            continue
        left_partner = pair_map.get(start - 1)
        right_partner = pair_map.get(end)
        if left_partner is None or right_partner is None:
            continue
        if start - 1 in paired_set and end in paired_set:
            other_unpaired = abs(left_partner - right_partner) - 1
            kind = "bulge" if other_unpaired <= 0 else "internal_loop"
            elements.append(
                {
                    "type": kind,
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "status": "COMPUTED",
                    "label": f"{kind.replace('_', ' ')} (dot-bracket)",
                }
            )
    return elements


def cache_key(
    *,
    sequence_hash: str,
    backend: str,
    tool_version: str,
    parameters: Mapping[str, Any],
) -> str:
    """Chave de cache: hash + ferramenta + versao + parametros.

    Args:
        sequence_hash: SHA-256 da RNA submetida.
        backend: viennarna_python ou rnafold_exe.
        tool_version: Versao detectada; vazia entra como unknown.
        parameters: Parametros efetivos.

    Returns:
        Digest SHA-256.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [
            "rnafold-v1",
            str(backend),
            tool_version or "unknown",
            str(sequence_hash or ""),
            repr(sorted((str(k), str(v)) for k, v in dict(parameters).items())),
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def mark_cached_result(result: Mapping[str, Any]) -> dict:
    """Marca resultado de cache sem fingir nova execucao.

    Args:
        result: Envelope de folding.

    Returns:
        Copia com cache_status=cached e timestamp original preservado.

    Raises:
        Nenhum.
    """
    copied = dict(result)
    copied["cache_status"] = "cached"
    copied["served_from_cache_at_utc"] = provenance.utc_now()
    return copied


def fold_rna(
    sequence: str,
    *,
    source: str = "user input",
    identifier: str = "",
    accession: str = "",
    version: str = "",
    organism: str = "",
    database: str = "",
    retrieved_at: str = "",
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
    full_sequence: str = "",
    backend: str = "auto",
    fold_fn: Optional[Callable[[str], Tuple[str, float]]] = None,
) -> dict:
    """Prediz a estrutura MFE com ViennaRNA. Sem ferramenta, TOOL_NOT_INSTALLED.

    Args:
        sequence: RNA a dobrar, ou a sequencia completa se houver regiao.
        source: Origem da sequencia (user input, workspace, NCBI Entrez).
        identifier, accession, version, organism, database, retrieved_at:
            Proveniencia da sequencia, nao da estrutura.
        region_start, region_end: Intervalo 0-based [start, end) se o folding
            for de uma subsequencia. None dobra a sequencia inteira.
        full_sequence: RNA completo quando uma regiao e dobrada.
        backend: auto, viennarna_python ou rnafold_exe.
        fold_fn: Motor injetavel para testes. Deve devolver (dot_bracket, mfe).

    Returns:
        Envelope RNASecondaryStructureResult.

    Raises:
        FoldingError.
    """
    parent = str(full_sequence or sequence)
    if region_start is not None or region_end is not None:
        folded = _extract_region(parent, region_start, region_end)
        region_mode = "selected subsequence"
    else:
        folded = sequence
        region_mode = "full sequence"
        region_start = 0
        region_end = None
    validated = validate_rna_for_folding(folded)
    residues = validated["sequence"]
    if region_end is None:
        region_end = (region_start or 0) + len(residues)
    parameters = _default_parameters(region_mode)
    backend_id, tool, method, tool_version = _select_backend(backend, fold_fn is not None)
    if fold_fn is not None:
        structure, mfe = _call_injected(fold_fn, residues)
        tool = tool or "injected ViennaRNA-compatible engine"
        method = method or "injected MFE"
        backend_id = backend_id or "injected"
    elif backend_id == "viennarna_python":
        structure, mfe = _fold_python(residues)
    elif backend_id == "rnafold_exe":
        structure, mfe = _fold_executable(residues)
    else:
        raise FoldingError(
            "RNA folding unavailable in this environment. ViennaRNA was not detected.",
            "TOOL_NOT_INSTALLED",
        )
    result = build_fold_result(
        sequence=residues,
        structure=structure,
        mfe_kcal_mol=mfe,
        tool=tool,
        tool_version=tool_version,
        method=method,
        backend=backend_id,
        parameters=parameters,
        source=source,
        identifier=identifier,
        accession=accession,
        version=version,
        organism=organism,
        database=database,
        retrieved_at=retrieved_at,
        region_mode=region_mode,
        region_start=int(region_start or 0),
        region_end=int(region_end),
        full_sequence_hash=provenance.sequence_digest(parent) if parent else validated["hash"],
        normalization=validated["normalization"],
    )
    if backend_id == "viennarna_python":
        engine_validation.record_live_validation(
            "ViennaRNA Python",
            ok=True,
            version=str(result.get("tool_version") or tool_version or ""),
            details={
                "backend": backend_id,
                "n_nt": len(residues),
                "mfe_kcal_mol": result.get("mfe_kcal_mol"),
            },
        )
    elif backend_id == "rnafold_exe":
        engine_validation.record_live_validation(
            "RNAfold",
            ok=True,
            version=str(result.get("tool_version") or tool_version or ""),
            details={"backend": backend_id, "n_nt": len(residues)},
        )
    return result


def build_fold_result(
    *,
    sequence: str,
    structure: str,
    mfe_kcal_mol: float,
    tool: str,
    tool_version: str,
    method: str,
    backend: str,
    parameters: Mapping[str, Any],
    source: str,
    identifier: str = "",
    accession: str = "",
    version: str = "",
    organism: str = "",
    database: str = "",
    retrieved_at: str = "",
    region_mode: str = "full sequence",
    region_start: int = 0,
    region_end: Optional[int] = None,
    full_sequence_hash: str = "",
    normalization: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Valida output do motor e monta o envelope tipado.

    Args:
        sequence: RNA submetido.
        structure: Dot-bracket do motor.
        mfe_kcal_mol: Energia reportada, kcal/mol.
        tool, tool_version, method, backend, parameters: proveniencia do folding.
        source: Origem da sequencia.
        identifier, accession, version, organism, database, retrieved_at:
            Metadados da sequencia.
        region_mode, region_start, region_end: Escopo dobrado.
        full_sequence_hash: Hash da molecula completa.
        normalization: Operacoes de texto aplicadas antes do folding.

    Returns:
        Envelope com status PREDICTED e kind predicted.

    Raises:
        FoldingError: PARSING_ERROR se a estrutura nao for coerente.
    """
    if not scientific_checks.rna_structure_length_matches(sequence, structure):
        raise FoldingError(
            "Dot-bracket length does not match the RNA sequence length.",
            "PARSING_ERROR",
        )
    if not scientific_checks.mfe_kcal_mol_is_valid(mfe_kcal_mol):
        raise FoldingError(
            "MFE is not a finite number. Energy was not replaced with zero.",
            "PARSING_ERROR",
        )
    pairs = parse_dot_bracket(structure)
    rebuilt = reconstruct_dot_bracket(len(sequence), pairs)
    if rebuilt != structure:
        raise FoldingError(
            "Reconstructed dot-bracket does not match the engine output.",
            "PARSING_ERROR",
        )
    validate_pairs_against_sequence(sequence, pairs)
    pair_records = [
        {
            "position_0based": left,
            "base": sequence[left],
            "paired_position_0based": right,
            "paired_base": sequence[right],
            "coordinate_system": "0-based",
        }
        for left, right in pairs
    ]
    elements = classify_structural_elements(structure, pairs)
    sequence_hash = provenance.sequence_digest(sequence)
    envelope = provenance.analysis_envelope(
        module="RNA folding",
        payload={
            "tool": tool,
            "tool_version": tool_version or "",
            "method": method,
            "mfe_kcal_mol": float(mfe_kcal_mol),
            "n_pairs": len(pairs),
            "length": len(sequence),
        },
        status="PREDICTED",
        algorithm=f"{tool} {method}".strip(),
        parameters=dict(parameters),
        source=source or "user input",
        sequence=sequence,
        accession=accession,
        organism=organism,
        input_identifier=identifier or accession or sequence_hash[:12],
    )
    envelope["status"] = "PREDICTED"
    envelope["kind"] = "predicted"
    envelope["structure_kind"] = "predicted"
    envelope["tool"] = tool
    envelope["tool_version"] = tool_version or ""
    envelope["version_status"] = "detected" if tool_version else "unavailable"
    envelope["method"] = method
    envelope["backend"] = backend
    envelope["parameters"] = dict(parameters)
    envelope["sequence"] = sequence
    envelope["sequence_hash"] = sequence_hash
    envelope["full_sequence_hash"] = full_sequence_hash or sequence_hash
    envelope["region_mode"] = region_mode
    envelope["region_start"] = int(region_start)
    envelope["region_end"] = int(region_end if region_end is not None else len(sequence))
    envelope["identifier"] = identifier
    envelope["accession"] = accession
    envelope["version"] = version
    envelope["organism"] = organism
    envelope["database"] = database
    envelope["retrieved_at"] = retrieved_at
    envelope["sequence_source"] = source
    envelope["structure_source"] = tool
    envelope["structure"] = structure
    envelope["dot_bracket"] = structure
    envelope["base_pairs"] = pair_records
    envelope["n_pairs"] = len(pair_records)
    envelope["mfe"] = float(mfe_kcal_mol)
    envelope["mfe_kcal_mol"] = float(mfe_kcal_mol)
    envelope["mfe_unit"] = "kcal/mol"
    envelope["elements"] = elements
    envelope["normalization"] = dict(normalization or {})
    envelope["cache_status"] = "live"
    envelope["computed_at_utc"] = provenance.utc_now()
    envelope["coordinate_system"] = "0-based"
    envelope["disclaimer"] = (
        "Predicted RNA secondary structure (minimum free energy). "
        "This is a computational prediction, not an experimental structure, "
        "not a 3D model, and not a functional annotation. "
        + (
            "Folding was performed on a selected subsequence; this is not the "
            "native local structure inside the full RNA. "
            if region_mode == "selected subsequence"
            else ""
        )
        + provenance.structure_disclaimer("predicted")
    )
    return envelope


def pair_at_position(result: Mapping[str, Any], position: int) -> Optional[dict]:
    """Devolve o par que envolve uma posicao 0-based, ou None se nao pareada.

    Args:
        result: Envelope validado.
        position: Coordenada 0-based.

    Returns:
        Registro de par ou None.

    Raises:
        FoldingError: INVALID_INPUT se a posicao for invalida.
    """
    sequence = str(result.get("sequence") or "")
    if position < 0 or position >= len(sequence):
        raise FoldingError("Position is outside the folded RNA.", "INVALID_INPUT")
    for item in list(result.get("base_pairs") or []):
        left = int(item.get("position_0based"))
        right = int(item.get("paired_position_0based"))
        if position in {left, right}:
            return dict(item)
    return None


def spans_for_position(result: Mapping[str, Any], position: int) -> dict:
    """Span de 1 nt para motif/MSA/3D futuros.

    Args:
        result: Envelope validado.
        position: Coordenada 0-based na RNA dobrada.

    Returns:
        scientific_checks.make_span.

    Raises:
        FoldingError: INVALID_INPUT.
    """
    sequence = str(result.get("sequence") or "")
    if position < 0 or position >= len(sequence):
        raise FoldingError("Position is outside the folded RNA.", "INVALID_INPUT")
    pair = pair_at_position(result, position)
    label = "paired" if pair else "unpaired"
    return scientific_checks.make_span(
        start=position,
        end=position + 1,
        sequence=sequence,
        strand="+",
        source="ViennaRNA predicted secondary structure",
        label=label,
        kind="rna_residue",
        status="PREDICTED",
    )


def map_msa_column_to_fold_position(
    coordinate_map_row: Sequence[Optional[int]], column: int
) -> Optional[int]:
    """Traduz coluna de MSA para coordenada original da RNA. Nao assume igualdade.

    Args:
        coordinate_map_row: maps[seq] do MSA.
        column: Coluna alinhada 0-based.

    Returns:
        Posicao original ou None se a coluna for gap.

    Raises:
        FoldingError: INVALID_INPUT se a coluna for invalida.
    """
    if column < 0 or column >= len(coordinate_map_row):
        raise FoldingError("MSA column is outside the coordinate map.", "INVALID_INPUT")
    value = coordinate_map_row[column]
    if value is None:
        return None
    return int(value)


def arc_paths(
    length: int, pairs: Sequence[Mapping[str, Any]]
) -> List[dict]:
    """Geometria de arcos no eixo da sequencia, derivada so dos pares reais.

    Args:
        length: Comprimento da RNA.
        pairs: Registros de base_pairs.

    Returns:
        Lista de paths com x, y, i, j.

    Raises:
        FoldingError: INVALID_INPUT se length for invalido.
    """
    if length < 1:
        raise FoldingError("Cannot build an arc diagram for an empty RNA.", "INVALID_INPUT")
    paths: List[dict] = []
    for item in pairs:
        left = int(item["position_0based"])
        right = int(item["paired_position_0based"])
        span = max(1, right - left)
        xs: List[float] = []
        ys: List[float] = []
        steps = max(8, min(40, span * 2))
        for step in range(steps + 1):
            t = step / steps
            xs.append(left + (right - left) * t)
            ys.append(math.sin(math.pi * t) * span)
        paths.append(
            {
                "x": xs,
                "y": ys,
                "i": left,
                "j": right,
                "base_i": str(item.get("base") or ""),
                "base_j": str(item.get("paired_base") or ""),
            }
        )
    return paths


def circular_layout(
    sequence: str, pairs: Sequence[Mapping[str, Any]]
) -> dict:
    """Layout circular: nucleotidos no circulo, cordas = pares reais.

    Args:
        sequence: RNA dobrada.
        pairs: Registros de base_pairs.

    Returns:
        Dict com x, y, bases, chords.

    Raises:
        FoldingError: INVALID_INPUT se vazia.

    Nota biologica:
        Isto nao e RNAPlot/NAView nem estrutura 3D. E um diagrama da topologia
        de pareamento predita.
    """
    residues = str(sequence or "")
    n = len(residues)
    if n < 1:
        raise FoldingError("Cannot build a circular layout for an empty RNA.", "INVALID_INPUT")
    xs = [math.cos(2.0 * math.pi * i / n) for i in range(n)]
    ys = [math.sin(2.0 * math.pi * i / n) for i in range(n)]
    chords = []
    for item in pairs:
        left = int(item["position_0based"])
        right = int(item["paired_position_0based"])
        chords.append(
            {
                "x": [xs[left], xs[right]],
                "y": [ys[left], ys[right]],
                "i": left,
                "j": right,
                "base_i": residues[left],
                "base_j": residues[right],
            }
        )
    return {"x": xs, "y": ys, "bases": list(residues), "chords": chords, "n": n}


def export_dot_bracket(result: Mapping[str, Any]) -> str:
    """Texto DBN com proveniencia em comentarios.

    Args:
        result: Envelope validado.

    Returns:
        Texto.

    Raises:
        Nenhum.
    """
    header = (
        f"# HelixScope predicted RNA secondary structure\n"
        f"# status={result.get('status')}\n"
        f"# kind={result.get('kind')}\n"
        f"# tool={result.get('tool')}\n"
        f"# tool_version={result.get('tool_version') or 'unavailable'}\n"
        f"# method={result.get('method')}\n"
        f"# parameters={result.get('parameters')}\n"
        f"# sequence_hash={result.get('sequence_hash')}\n"
        f"# sequence_source={result.get('sequence_source')}\n"
        f"# structure_source={result.get('structure_source')}\n"
        f"# mfe_kcal_mol={result.get('mfe_kcal_mol')}\n"
        f"# computed_at_utc={result.get('computed_at_utc')}\n"
        f"# {result.get('disclaimer')}\n"
        f">{result.get('identifier') or result.get('accession') or 'rna'}\n"
    )
    return header + str(result.get("sequence") or "") + "\n" + str(result.get("dot_bracket") or "") + "\n"


def export_pairs_csv(result: Mapping[str, Any]) -> str:
    """CSV dos pares. Zero real permanece 0; ausente vira N/A.

    Args:
        result: Envelope validado.

    Returns:
        Texto CSV.

    Raises:
        Nenhum.
    """
    import csv
    import io as io_module

    buffer = io_module.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=[
            "position_0based",
            "base",
            "paired_position_0based",
            "paired_base",
            "mfe_kcal_mol",
            "tool",
            "tool_version",
            "sequence_hash",
            "status",
        ],
    )
    writer.writeheader()
    version = result.get("tool_version") or "N/A"
    mfe = provenance.csv_cell(result.get("mfe_kcal_mol"))
    rows = list(result.get("base_pairs") or [])
    if not rows:
        writer.writerow(
            {
                "position_0based": "N/A",
                "base": "N/A",
                "paired_position_0based": "N/A",
                "paired_base": "N/A",
                "mfe_kcal_mol": mfe,
                "tool": result.get("tool"),
                "tool_version": version,
                "sequence_hash": result.get("sequence_hash"),
                "status": result.get("status"),
            }
        )
        return buffer.getvalue()
    for item in rows:
        writer.writerow(
            {
                "position_0based": item.get("position_0based"),
                "base": item.get("base"),
                "paired_position_0based": item.get("paired_position_0based"),
                "paired_base": item.get("paired_base"),
                "mfe_kcal_mol": mfe,
                "tool": result.get("tool"),
                "tool_version": version,
                "sequence_hash": result.get("sequence_hash"),
                "status": result.get("status"),
            }
        )
    return buffer.getvalue()


def export_result_bundle(result: Mapping[str, Any]) -> dict:
    """Copia JSON-safe. NaN nao vira 0.

    Args:
        result: Envelope.

    Returns:
        Dict serializavel.

    Raises:
        Nenhum.
    """
    return provenance.json_safe(dict(result))


def _default_parameters(region_mode: str) -> dict:
    return {
        "algorithm": "MFE (ViennaRNA Zuker-style dynamic programming)",
        "temperature_c": VIENNA_DEFAULT_TEMPERATURE_C,
        "temperature_source": (
            "ViennaRNA documented default. RNAfold is invoked with --temp="
            f"{VIENNA_DEFAULT_TEMPERATURE_C}. Python RNA.cvar.temperature is "
            "set when that attribute exists."
        ),
        "temperature_passed_by_helixscope": True,
        "salt_M": VIENNA_DEFAULT_SALT_M,
        "salt_source": (
            "ViennaRNA documented RNAfold default. HelixScope does not pass "
            "--salt; the library value is used."
        ),
        "salt_passed_by_helixscope": False,
        "energy_parameters": (
            "ViennaRNA default nearest-neighbor set (Turner 2004 in unmodified "
            "ViennaRNA 2 distributions). HelixScope does not load a custom file."
        ),
        "noconv": True,
        "partition_function": False,
        "pseudoknots": False,
        "noPS": True,
        "region_mode": region_mode,
    }


def _select_backend(requested: str, injected: bool) -> Tuple[str, str, str, str]:
    if injected:
        return "injected", "injected engine", "injected MFE", ""
    python = detect_python_bindings()
    binary = detect_rnafold_executable()
    choice = str(requested or "auto").strip()
    if choice == "auto":
        if python["available"]:
            choice = "viennarna_python"
        elif binary["available"]:
            choice = "rnafold_exe"
        else:
            return "", "", "", ""
    if choice == "viennarna_python":
        if not python["available"]:
            raise FoldingError(
                "ViennaRNA Python module RNA is not installed.",
                "TOOL_NOT_INSTALLED",
            )
        return (
            "viennarna_python",
            "ViennaRNA RNAlib",
            "RNA.fold MFE",
            str(python.get("version") or ""),
        )
    if choice == "rnafold_exe":
        if not binary["available"]:
            raise FoldingError(
                "RNAfold executable is not installed.",
                "TOOL_NOT_INSTALLED",
            )
        return (
            "rnafold_exe",
            "ViennaRNA RNAfold",
            "RNAfold --noPS --noconv MFE",
            str(binary.get("version") or ""),
        )
    raise FoldingError(f"Folding backend '{choice}' is not offered.", "INVALID_INPUT")


def _call_injected(
    fold_fn: Callable[[str], Tuple[str, float]], sequence: str
) -> Tuple[str, float]:
    try:
        structure, mfe = fold_fn(sequence)
    except FoldingError:
        raise
    except Exception as exc:
        raise FoldingError(f"Injected folding engine failed: {exc}.", "TOOL_FAILED") from exc
    return str(structure), float(mfe)


def _apply_library_temperature(module: Any) -> None:
    """Aplica o default documentado de 37 C quando a API Python expoe o campo.

    Args:
        module: Modulo RNA importado.

    Returns:
        None.

    Raises:
        Nenhum. A ausencia do atributo nao e fabricacao de energia.
    """
    cvar = getattr(module, "cvar", None)
    if cvar is not None and hasattr(cvar, "temperature"):
        cvar.temperature = VIENNA_DEFAULT_TEMPERATURE_C


def _fold_python(sequence: str) -> Tuple[str, float]:
    try:
        import RNA  # type: ignore
    except ImportError as exc:
        raise FoldingError(
            "ViennaRNA Python module RNA is not installed.",
            "TOOL_NOT_INSTALLED",
        ) from exc
    if not hasattr(RNA, "fold"):
        raise FoldingError(
            "The imported RNA module does not expose fold(). This is TOOL_FAILED, "
            "not a computed structure.",
            "TOOL_FAILED",
        )

    def _run() -> Tuple[str, float]:
        _apply_library_temperature(RNA)
        structure, mfe = RNA.fold(sequence)
        return str(structure), float(mfe)

    try:
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(_run)
            return future.result(timeout=LOCAL_TIMEOUT_S)
    except FutureTimeout as exc:
        raise FoldingError(
            f"ViennaRNA RNA.fold exceeded {int(LOCAL_TIMEOUT_S)}s. This is TIMEOUT, "
            "not an empty structure.",
            "TIMEOUT",
        ) from exc
    except FoldingError:
        raise
    except Exception as exc:
        raise FoldingError(f"ViennaRNA RNA.fold failed: {exc}.", "TOOL_FAILED") from exc


def _fold_executable(sequence: str) -> Tuple[str, float]:
    info = detect_rnafold_executable()
    executable = str(info.get("path") or "")
    if not executable or not info.get("available"):
        raise FoldingError("RNAfold executable is not installed.", "TOOL_NOT_INSTALLED")
    tmpdir = tempfile.mkdtemp(prefix="helixscope_rnafold_")
    try:
        args = [
            executable,
            "--noPS",
            "--noconv",
            f"--temp={VIENNA_DEFAULT_TEMPERATURE_C}",
        ]
        try:
            completed = subprocess.run(
                args,
                input=sequence + "\n",
                capture_output=True,
                timeout=LOCAL_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmpdir,
                text=True,
            )
        except subprocess.TimeoutExpired as exc:
            raise FoldingError(
                f"RNAfold exceeded {int(LOCAL_TIMEOUT_S)}s. This is TIMEOUT, "
                "not an empty structure.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise FoldingError(f"RNAfold failed to start: {exc}.", "TOOL_FAILED") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()[:300]
            raise FoldingError(
                f"RNAfold exited {completed.returncode}. {detail}",
                "TOOL_FAILED",
            )
        parsed = parse_rnafold_output(completed.stdout or "", sequence)
        return str(parsed["structure"]), float(parsed["mfe_kcal_mol"])
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _parse_structure_line(line: str) -> dict:
    match = ENERGY_TAIL.match(line)
    if not match:
        raise FoldingError(
            "RNAfold structure line is malformed.",
            "PARSING_ERROR",
        )
    energy_text = match.group("energy")
    if energy_text is None:
        mfe: Optional[float] = None
    else:
        mfe = float(energy_text)
    return {"structure": match.group("structure"), "mfe_kcal_mol": mfe}


def _extract_region(sequence: str, start: Optional[int], end: Optional[int]) -> str:
    text = str(sequence or "")
    if start is None:
        start = 0
    if end is None:
        end = len(text)
    try:
        start_i = int(start)
        end_i = int(end)
    except (TypeError, ValueError) as exc:
        raise FoldingError("Region coordinates must be integers.", "INVALID_INPUT") from exc
    if start_i < 0 or end_i > len(text) or start_i >= end_i:
        raise FoldingError(
            "Selected region is outside the RNA sequence.",
            "INVALID_INPUT",
        )
    return text[start_i:end_i]


def _stem_record(stacked: Sequence[Tuple[int, int]]) -> dict:
    lefts = [item[0] for item in stacked]
    rights = [item[1] for item in stacked]
    return {
        "type": "stem",
        "start": min(lefts),
        "end": max(rights) + 1,
        "length": len(stacked),
        "status": "COMPUTED",
        "label": "stem (dot-bracket)",
        "outer_pair": [min(lefts), max(rights)],
    }


def _unpaired_runs(structure: str) -> List[Tuple[int, int]]:
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for index, char in enumerate(structure):
        if char == ".":
            if start is None:
                start = index
        elif start is not None:
            runs.append((start, index))
            start = None
    if start is not None:
        runs.append((start, len(structure)))
    return runs


def _module_version(module: Any) -> str:
    for attr in ("__version__", "VERSION", "version"):
        value = getattr(module, attr, None)
        if value:
            text = str(value).strip()
            if text:
                return text[:80]
    getter = getattr(module, "get_version", None)
    if callable(getter):
        try:
            text = str(getter()).strip()
        except Exception:
            text = ""
        if text:
            return text[:80]
    return ""


def _probe_rnafold_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            timeout=5,
            check=False,
            shell=False,
            text=True,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = ((completed.stdout or "") + (completed.stderr or "")).strip()
    line = text.splitlines()[0] if text else ""
    return line[:80]


def _normalized_tool_basename(path: str) -> str:
    base = os.path.basename(str(path or "")).lower()
    for suffix in (".exe", ".bat", ".cmd"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def _allowlisted_basenames(names: Sequence[str]) -> set[str]:
    allowed: set[str] = set()
    for item in names:
        low = str(item).lower()
        allowed.add(low)
        allowed.add(_normalized_tool_basename(low))
    return allowed


def folding_admission(length: int) -> dict:
    """Decide o modo de estrutura sem executar o engine.

    Args:
        length: Comprimento da RNA ja validada.

    Returns:
        Dict com global_fold, local_long_rna e status. Sequencias acima de
        MAX_FOLD_NT nao recebem folding global. Acima de MAX_LOCAL_RNA_NT o
        status e RESOURCE_LIMIT.

    Raises:
        ValueError: Se length nao for um inteiro nao negativo.

    Nota biologica:
        Folding global de MFE e cubico no comprimento. Analise local em janela
        nao e a estrutura global da molecula.
    """
    if isinstance(length, bool) or not isinstance(length, int) or length < 0:
        raise ValueError("length must be a non-negative integer.")
    if length > MAX_LOCAL_RNA_NT:
        return {
            "global_fold": "RESOURCE_LIMIT",
            "local_long_rna": "RESOURCE_LIMIT",
            "status": "RESOURCE_LIMIT",
            "reason": (
                f"RNA length {length} exceeds the local-analysis limit of "
                f"{MAX_LOCAL_RNA_NT} nt."
            ),
        }
    global_state = "AVAILABLE" if length <= MAX_FOLD_NT else "RESOURCE_LIMIT"
    return {
        "global_fold": global_state,
        "local_long_rna": "AVAILABLE",
        "status": "AVAILABLE",
        "reason": (
            "Global MFE is limited to "
            f"{MAX_FOLD_NT} nt. Local pair probabilities use ViennaRNA pfl_fold "
            "and are not a global structure."
        ),
    }


def ensemble_metrics(sequence: str) -> dict:
    """MFE, ensemble, centroide, MEA e diversidade para RNA curta.

    Usa RNA.fold_compound da ViennaRNA instalada. Nao estima valores se a
    chamada falhar.

    Args:
        sequence: RNA canonica A, C, G, U.

    Returns:
        Dict PREDICTED com estruturas, energias, diversidade, comprimentos e
        proveniencia do engine.

    Raises:
        FoldingError: INVALID_INPUT, RESOURCE_LIMIT, TOOL_NOT_INSTALLED,
            TOOL_FAILED ou TIMEOUT.

    Nota biologica:
        A estrutura centroide e a MEA sao estruturas representativas do
        ensemble, nao a estrutura experimental e nao um modelo 3D. A
        diversidade e a distancia media de pares no ensemble.
    """
    checked = validate_rna_for_folding(sequence)
    residues = str(checked["sequence"])
    if len(residues) > MAX_FOLD_NT:
        raise FoldingError(
            f"Ensemble metrics use the global partition function and are limited "
            f"to {MAX_FOLD_NT} nt (RESOURCE_LIMIT). Use local long-RNA analysis "
            "for longer sequences. No ensemble value was invented.",
            "RESOURCE_LIMIT",
        )
    try:
        import RNA  # type: ignore
    except ImportError as exc:
        raise FoldingError(
            "ViennaRNA Python module RNA is not installed.",
            "TOOL_NOT_INSTALLED",
        ) from exc

    def _run() -> dict:
        _apply_library_temperature(RNA)
        compound = RNA.fold_compound(residues)
        mfe_structure, mfe = compound.mfe()
        partition_structure, ensemble_energy = compound.pf()
        centroid_structure, centroid_distance = compound.centroid()
        mea_structure, mea_value = compound.MEA()
        diversity = compound.mean_bp_distance()
        return {
            "mfe_structure": str(mfe_structure),
            "mfe_kcal_mol": float(mfe),
            "partition_structure": str(partition_structure),
            "ensemble_free_energy_kcal_mol": float(ensemble_energy),
            "centroid_structure": str(centroid_structure),
            "centroid_ensemble_distance": float(centroid_distance),
            "mea_structure": str(mea_structure),
            "mea_value": float(mea_value),
            "ensemble_diversity": float(diversity),
        }

    try:
        with ThreadPoolExecutor(max_workers=1) as pool:
            payload = pool.submit(_run).result(timeout=LOCAL_TIMEOUT_S)
    except FutureTimeout as exc:
        raise FoldingError(
            f"ViennaRNA ensemble metrics exceeded {int(LOCAL_TIMEOUT_S)}s. "
            "This is TIMEOUT, not an empty structure.",
            "TIMEOUT",
        ) from exc
    except FoldingError:
        raise
    except Exception as exc:
        raise FoldingError(
            f"ViennaRNA ensemble metrics failed: {exc}.",
            "TOOL_FAILED",
        ) from exc
    version = _module_version(RNA)
    return {
        "status": "PREDICTED",
        "mode": "GLOBAL_ENSEMBLE",
        "input_length": len(residues),
        "analyzed_length": len(residues),
        "visualized_length": len(residues),
        "method": "ViennaRNA fold_compound MFE, partition function, centroid, MEA",
        "engine": "ViennaRNA Python",
        "engine_version": version,
        "sequence_hash": checked["hash"],
        "disclaimer": (
            "Predicted ensemble statistics are not an experimental structure "
            "and not an RNA 3D model."
        ),
        **payload,
    }


def local_long_rna_analysis(
    sequence: str,
    *,
    window: int = LOCAL_WINDOW_DEFAULT,
    max_pair_span: int | None = None,
    cutoff: float = LOCAL_PAIR_CUTOFF,
) -> dict:
    """Probabilidades de par local via ViennaRNA pfl_fold.

    Nao chama RNA.fold nem a funcao de particao global. O resultado e uma
    analise local, mesmo quando a sequencia e curta.

    Args:
        sequence: RNA canonica A, C, G, U, ate MAX_LOCAL_RNA_NT.
        window: Janela deslizante. Default 70, o -W documentado de RNAplfold.
        max_pair_span: Distancia maxima do par. None usa a janela. Nao pode
            exceder a janela.
        cutoff: Probabilidade minima devolvida por pfl_fold.

    Returns:
        Dict PREDICTED com comprimentos, parametros, contagem de pares na
        lista completa, preview rotulado e resumo de probabilidade desemparelhada.

    Raises:
        FoldingError: INVALID_INPUT, RESOURCE_LIMIT, TOOL_NOT_INSTALLED,
            TOOL_FAILED ou TIMEOUT.

    Nota biologica:
        pfl_fold estima pares dentro de janelas. Nao e a estrutura secundaria
        global da RNA e nao e um modelo 3D.
    """
    if isinstance(window, bool) or not isinstance(window, int) or window < 2:
        raise FoldingError("Local window must be an integer >= 2.", "INVALID_INPUT")
    span = window if max_pair_span is None else max_pair_span
    if isinstance(span, bool) or not isinstance(span, int) or span < 1 or span > window:
        raise FoldingError(
            "max_pair_span must be an integer between 1 and the window size.",
            "INVALID_INPUT",
        )
    if not isinstance(cutoff, (int, float)) or isinstance(cutoff, bool):
        raise FoldingError("cutoff must be a number in (0, 1].", "INVALID_INPUT")
    cutoff_value = float(cutoff)
    if not 0.0 < cutoff_value <= 1.0:
        raise FoldingError("cutoff must be a number in (0, 1].", "INVALID_INPUT")
    stripped = "".join(sequence.split()).upper() if isinstance(sequence, str) else ""
    if len(stripped) > MAX_LOCAL_RNA_NT:
        raise FoldingError(
            f"Local RNA analysis is limited to {MAX_LOCAL_RNA_NT} nt "
            "(RESOURCE_LIMIT). The sequence was not shortened.",
            "RESOURCE_LIMIT",
        )
    if len(stripped) > dna_analysis.MAX_INPUT_RESIDUES:
        raise FoldingError(
            "Sequence exceeds the technical input limit.",
            "RESOURCE_LIMIT",
        )
    checked = validate_rna_for_folding(sequence) if len(stripped) <= MAX_FOLD_NT else None
    if checked is None:
        if not isinstance(sequence, str):
            raise FoldingError("RNA sequence must be a string.", "INVALID_INPUT")
        if not stripped:
            raise FoldingError("RNA sequence is empty.", "INVALID_INPUT")
        if "T" in stripped:
            raise FoldingError(
                "Sequence contains T. DNA is not sent to RNA folding and T is not "
                "silently converted to U.",
                "INVALID_INPUT",
            )
        if any(base not in dna_analysis.CANONICAL_RNA_ALPHABET for base in stripped):
            raise FoldingError(
                "Local RNA analysis accepts only A, C, G and U.",
                "INVALID_INPUT",
            )
        residues = stripped
        sequence_hash = provenance.sequence_digest(residues)
    else:
        residues = str(checked["sequence"])
        sequence_hash = str(checked["hash"])
    try:
        import RNA  # type: ignore
    except ImportError as exc:
        raise FoldingError(
            "ViennaRNA Python module RNA is not installed.",
            "TOOL_NOT_INSTALLED",
        ) from exc

    def _run() -> tuple[list, list]:
        _apply_library_temperature(RNA)
        raw_pairs = RNA.pfl_fold(residues, window, span, cutoff_value)
        unpaired = RNA.pfl_fold_up(residues, 1, window, span)
        pairs: list[dict] = []
        for item in list(raw_pairs or []):
            left = int(item.i)
            right = int(item.j)
            if left <= 0 or right <= 0:
                continue
            pairs.append(
                {
                    "i": left,
                    "j": right,
                    "probability": float(item.p),
                }
            )
        access: list[float] = []
        matrix = list(unpaired or [])
        for position in range(1, len(residues) + 1):
            if position >= len(matrix) or len(matrix[position]) < 2:
                break
            access.append(float(matrix[position][1]))
        return pairs, access

    try:
        with ThreadPoolExecutor(max_workers=1) as pool:
            pairs, access = pool.submit(_run).result(timeout=LOCAL_ANALYSIS_TIMEOUT_S)
    except FutureTimeout as exc:
        raise FoldingError(
            f"ViennaRNA pfl_fold exceeded {int(LOCAL_ANALYSIS_TIMEOUT_S)}s. "
            "This is TIMEOUT, not an empty structure.",
            "TIMEOUT",
        ) from exc
    except FoldingError:
        raise
    except Exception as exc:
        raise FoldingError(f"ViennaRNA pfl_fold failed: {exc}.", "TOOL_FAILED") from exc
    access_status = "COMPUTED" if len(access) == len(residues) else "PARTIAL"
    if access and access_status == "COMPUTED":
        access_mean = sum(access) / len(access)
        access_min = min(access)
        access_max = max(access)
    else:
        access_mean = None
        access_min = None
        access_max = None
    ordered = sorted(pairs, key=lambda item: item["probability"], reverse=True)
    preview = ordered[:LOCAL_PREVIEW_PAIRS]
    probabilities = [item["probability"] for item in pairs]
    return {
        "status": "PREDICTED",
        "mode": "LOCAL_LONG_RNA",
        "input_length": len(residues),
        "analyzed_length": len(residues),
        "visualized_length": len(preview),
        "visualization": "PREVIEW",
        "window_size": window,
        "max_pair_span": span,
        "cutoff": cutoff_value,
        "method": "ViennaRNA RNA.pfl_fold local pair probabilities",
        "engine": "ViennaRNA Python",
        "engine_version": _module_version(RNA),
        "sequence_hash": sequence_hash,
        "pair_count": len(pairs),
        "max_pair_probability": max(probabilities) if probabilities else None,
        "mean_pair_probability": (
            sum(probabilities) / len(probabilities) if probabilities else None
        ),
        "pairs_preview": preview,
        "pairs_omitted_from_preview": max(0, len(pairs) - len(preview)),
        "preview_is_complete_catalogue": len(preview) == len(pairs),
        "unpaired_probability_status": access_status,
        "unpaired_probability_mean": access_mean,
        "unpaired_probability_min": access_min,
        "unpaired_probability_max": access_max,
        "global_structure": False,
        "disclaimer": (
            "Local pair probabilities are not the global secondary structure, "
            "not an experimental structure, and not an RNA 3D model."
        ),
    }


def _resolve_allowlisted_executable(names: Sequence[str]) -> Optional[str]:
    from . import tool_detection

    from . import tool_paths

    extra = list(tool_detection.VIENNARNA_WINDOWS_CANDIDATES)
    extra.extend(tool_paths.candidates_for(names))
    return tool_detection.resolve_allowlisted_executable(
        names,
        extra_file_candidates=tuple(extra),
        env_var="HELIXSCOPE_RNAFOLD",
        which_fn=shutil.which,
    )
