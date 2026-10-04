"""Inferencia filogenetica a partir de um MSA validado.

Metodos locais: Neighbor-Joining (Saitou e Nei 1987) e UPGMA via Biopython
DistanceTreeConstructor, sobre uma matriz de distancia calculada aqui
(p-distance, Jukes-Cantor 1969 ou Kimura 1980). FastTree e IQ-TREE so correm
se o executavel allowlisted estiver instalado; a ausencia e UNAVAILABLE, nunca
uma reimplementacao rotulada com o nome da ferramenta.

A topologia e os comprimentos de ramo saem do algoritmo. Nao ha arvore
manual, bootstrap inventado, relacao de especies inventada nem escala
arbitraria. Distancia computacional nao e verdade evolutiva. Arvore
filogenetica nao e arvore taxonomica.

Nenhuma funcao importa Streamlit. Subprocesso usa argv estruturado e
shell=False.
"""

from __future__ import annotations

import copy
import hashlib
import io
import math
import os
import random
import re
import subprocess
import tempfile
import time
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from . import provenance, tool_detection, tool_paths

try:
    from Bio import __version__ as BIOPYTHON_VERSION
except Exception:  # pragma: no cover - Bio sempre presente neste projeto
    BIOPYTHON_VERSION = ""

RunFn = Callable[..., subprocess.CompletedProcess]

GAP_CHARS: str = "-.~"
DNA_CANONICAL: str = "ACGT"
RNA_CANONICAL: str = "ACGU"
PROTEIN_CANONICAL: str = "ACDEFGHIKLMNPQRSTVWY"
PURINES: frozenset[str] = frozenset({"A", "G"})
PYRIMIDINES: frozenset[str] = frozenset({"C", "T", "U"})
NEWICK_UNSAFE: re.Pattern[str] = re.compile(r"[():,;\[\]\s'\"]+")

MAX_SEQUENCES: int = 30
MAX_ALIGNMENT_COLUMNS: int = 8000
MAX_BOOTSTRAP_REPLICATES: int = 50
MAX_NEWICK_CHARS: int = 200_000
MAX_ML_OUTPUT_CHARS: int = 200_000
MIN_SEQUENCES: int = 2
DEFAULT_BOOTSTRAP_SEED: int = 1
ML_TIMEOUT_S: float = 90.0
MIN_UFBOOT_REPLICATES: int = 1000
MIN_ALRT_REPLICATES: int = 1000
MAX_UFBOOT_REPLICATES: int = 1000
MAX_ALRT_REPLICATES: int = 1000
MAX_IQTREE_REPORT_CHARS: int = 200_000
IQTREE_MODEL_USER: str = "user"
IQTREE_MODEL_MFP: str = "MFP"
SUPPORT_FELSENSTEIN: str = "FELSENSTEIN_BOOTSTRAP"
SUPPORT_UFBOOT: str = "ULTRAFAST_BOOTSTRAP"
SUPPORT_SH_ALRT: str = "SH_ALRT"
IQTREE_SUPPORT_TOKEN: re.Pattern[str] = re.compile(
    r"^[0-9]+(?:\.[0-9]+)?(?:/[0-9]+(?:\.[0-9]+)?)?$"
)
BEST_FIT_MODEL_RE: re.Pattern[str] = re.compile(
    r"Best-fit model according to (AIC|AICc|BIC):\s*(\S+)",
    re.I,
)

DISTANCE_P: str = "p_distance"
DISTANCE_JC69: str = "jukes_cantor_1969"
DISTANCE_K2P: str = "kimura_1980"
DISTANCE_BLOSUM62: str = "blosum62_biopython"

METHOD_NJ: str = "neighbor_joining"
METHOD_UPGMA: str = "upgma"
METHOD_FASTTREE: str = "fasttree_ml"
METHOD_IQTREE: str = "iqtree_ml"

ROOTING_UNROOTED: str = "unrooted"
ROOTING_UPGMA: str = "upgma_rooted"
ROOTING_MIDPOINT: str = "midpoint"
ROOTING_OUTGROUP: str = "outgroup"

FASTTREE_NAMES: tuple[str, ...] = (
    "fasttree",
    "FastTree",
    "fasttree.exe",
    "FastTree.exe",
    "FastTreeMP",
    "FastTreeMP.exe",
)
IQTREE_NAMES: tuple[str, ...] = (
    "iqtree",
    "iqtree2",
    "iqtree3",
    "iqtree.exe",
    "iqtree2.exe",
    "iqtree3.exe",
)
FASTTREE_ENV: str = "HELIXSCOPE_FASTTREE"
IQTREE_ENV: str = "HELIXSCOPE_IQTREE"

SUPPORTED_METHODS: tuple[str, ...] = (
    METHOD_NJ,
    METHOD_UPGMA,
    METHOD_FASTTREE,
    METHOD_IQTREE,
)
NUCLEOTIDE_DISTANCES: tuple[str, ...] = (DISTANCE_P, DISTANCE_JC69, DISTANCE_K2P)
PROTEIN_DISTANCES: tuple[str, ...] = (DISTANCE_P, DISTANCE_BLOSUM62)


class PhylogenyError(Exception):
    """Falha classificada da pipeline filogenetica.

    Attributes:
        category: INVALID_INPUT, PARSING_ERROR, RESOURCE_LIMIT,
            TOOL_NOT_INSTALLED, TIMEOUT, UNAVAILABLE ou INTERNAL.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "INTERNAL")


def tool_availability() -> dict:
    """Declara metodos filogeneticos reais e ferramentas opcionais.

    Args:
        Nenhum.

    Returns:
        Dict com neighbor_joining, upgma, fasttree, iqtree, limites e razoes.
        FastTree/IQ-TREE available e True so se o binario allowlisted existir
        e a versao, se obtida, nao for inventada.

    Raises:
        Nenhum.

    Nota biologica:
        NJ e UPGMA sao algoritmos publicados executados pelo Biopython.
        FastTree (Price et al.) e IQ-TREE (Nguyen/Minh et al.) nao sao
        reimplementados. UPGMA assume relogio molecular / ultrametricidade.
    """
    fasttree = detect_fasttree()
    iqtree = detect_iqtree()
    return {
        "neighbor_joining": True,
        "upgma": True,
        "distance_methods": {
            "DNA": list(NUCLEOTIDE_DISTANCES),
            "RNA": list(NUCLEOTIDE_DISTANCES),
            "PROTEIN": list(PROTEIN_DISTANCES),
        },
        "fasttree": fasttree,
        "iqtree": iqtree,
        "maximum_likelihood_local": bool(fasttree["available"] or iqtree["available"]),
        "bootstrap": {
            "method": "Felsenstein 1985 column resampling",
            "available_for": [METHOD_NJ, METHOD_UPGMA],
            "max_replicates": MAX_BOOTSTRAP_REPLICATES,
            "not_available_for_ml_unless_tool_emits_support": True,
        },
        "iqtree_support": {
            "ufboot": {
                "method": SUPPORT_UFBOOT,
                "available_when_iqtree_installed": bool(iqtree["available"]),
                "min_replicates": MIN_UFBOOT_REPLICATES,
                "not_felsenstein_bootstrap": True,
            },
            "sh_alrt": {
                "method": SUPPORT_SH_ALRT,
                "available_when_iqtree_installed": bool(iqtree["available"]),
                "min_replicates": MIN_ALRT_REPLICATES,
                "not_ufboot": True,
            },
            "modelfinder": {
                "keyword": "MFP",
                "available_when_iqtree_installed": bool(iqtree["available"]),
                "criterion_default": "BIC",
            },
        },
        "limits": {
            "max_sequences": MAX_SEQUENCES,
            "max_alignment_columns": MAX_ALIGNMENT_COLUMNS,
            "max_bootstrap_replicates": MAX_BOOTSTRAP_REPLICATES,
        },
        "biopython_version": BIOPYTHON_VERSION or "",
        "summary_reason": (
            "Phylogenetic inference uses Biopython Neighbor-Joining or UPGMA "
            "on a published distance, or FastTree/IQ-TREE when installed. "
            "A distance is not phylogenetic truth. The tree is not taxonomy."
        ),
    }


def detect_fasttree() -> dict:
    """Detecta FastTree allowlisted e a versao reportada pelo binario.

    Args:
        Nenhum.

    Returns:
        Dict available, path sanitizado, version, status, reason.

    Raises:
        Nenhum.
    """
    return _detect_ml_binary(
        names=FASTTREE_NAMES,
        env_var=FASTTREE_ENV,
        tool="FastTree",
        tool_id="fasttree",
        version_args=("-expert",),
        version_pattern=re.compile(
            r"FastTree(?:\s+Version)?\s+([0-9][0-9.\w-]*)",
            re.I,
        ),
    )


def detect_iqtree() -> dict:
    """Detecta IQ-TREE allowlisted e a versao reportada pelo binario.

    Args:
        Nenhum.

    Returns:
        Dict available, path sanitizado, version, status, reason.

    Raises:
        Nenhum.
    """
    return _detect_ml_binary(
        names=IQTREE_NAMES,
        env_var=IQTREE_ENV,
        tool="IQ-TREE",
        tool_id="iqtree",
        version_args=("-version",),
        version_pattern=re.compile(r"IQ-TREE[^\n]*?version\s+([0-9][0-9.\w-]*)", re.I),
    )


def alignment_identity_hash(msa_result: Mapping[str, Any]) -> str:
    """Hash SHA-256 da identidade do MSA usado como entrada da arvore.

    Args:
        msa_result: Envelope de MSA (linhas alinhadas, hashes, molecula).

    Returns:
        Hex digest. Nao e o hash de uma sequencia isolada.

    Raises:
        PhylogenyError: INVALID_INPUT se o envelope nao tiver linhas.

    Nota biologica:
        A identidade do alinhamento inclui ordem, IDs, hashes e colunas
        alinhadas. Metadados de organismo nao entram neste hash.
    """
    rows = list(msa_result.get("rows") or [])
    if not rows:
        raise PhylogenyError("MSA has no aligned rows.", "INVALID_INPUT")
    payload = hashlib.sha256()
    payload.update(str(msa_result.get("molecule") or "").encode("utf-8"))
    payload.update(b"\n")
    payload.update(str(msa_result.get("alignment_length") or "").encode("utf-8"))
    payload.update(b"\n")
    for digest in list(msa_result.get("input_hashes") or []):
        payload.update(str(digest).encode("utf-8"))
        payload.update(b"\n")
    for row in rows:
        payload.update(str(row.get("identifier") or "").encode("utf-8"))
        payload.update(b"\n")
        payload.update(str(row.get("aligned") or "").encode("utf-8"))
        payload.update(b"\n")
        payload.update(str(row.get("hash") or "").encode("utf-8"))
        payload.update(b"\n")
    return payload.hexdigest()


def unique_leaf_labels(
    identifiers: Sequence[str],
    hashes: Sequence[str],
) -> Tuple[List[dict], List[dict]]:
    """Gera IDs unicos e seguros para Newick sem apagar sequencias.

    Args:
        identifiers: Identificadores originais na ordem do MSA.
        hashes: Hashes das sequencias na mesma ordem.

    Returns:
        (leaf_records, transformations). transformations registra cada
        normalizacao. Duplicatas de sequencia nao sao removidas.

    Raises:
        PhylogenyError: INVALID_INPUT se as listas tiverem tamanhos diferentes
            ou estiverem vazias.
    """
    if len(identifiers) != len(hashes) or not identifiers:
        raise PhylogenyError(
            "Leaf identifiers and sequence hashes must be non-empty and aligned.",
            "INVALID_INPUT",
        )
    used: set[str] = set()
    records: List[dict] = []
    transformations: List[dict] = []
    for index, (raw_id, digest) in enumerate(zip(identifiers, hashes)):
        original = str(raw_id or "").strip() or f"seq_{index + 1}"
        safe = NEWICK_UNSAFE.sub("_", original).strip("_") or f"seq_{index + 1}"
        candidate = safe
        suffix = 1
        while candidate.lower() in used:
            candidate = f"{safe}_{str(digest)[:8]}_{suffix}"
            suffix += 1
        used.add(candidate.lower())
        changed = candidate != original
        if changed:
            transformations.append(
                {
                    "original_id": original,
                    "tree_id": candidate,
                    "sequence_hash": str(digest),
                    "reason": "Newick-safe unique leaf identifier",
                }
            )
        records.append(
            {
                "original_id": original,
                "tree_id": candidate,
                "sequence_hash": str(digest),
                "index": index,
            }
        )
    return records, transformations


def validate_msa_for_phylogeny(msa_result: Mapping[str, Any]) -> dict:
    """Valida um MSA COMPLETED antes da inferencia.

    Args:
        msa_result: Envelope produzido por modules.msa.

    Returns:
        Dict com molecule, n_sequences, alignment_length, leaf_records,
        transformations, identical_groups, alignment_hash, aligned_rows.

    Raises:
        PhylogenyError: INVALID_INPUT, RESOURCE_LIMIT.

    Nota biologica:
        DNA, RNA e proteina nao sao misturados. Duplicatas de sequencia sao
        identificadas e preservadas. Sem fallback para FASTA bruto.
    """
    status = str(msa_result.get("status") or "")
    if status not in {"COMPLETED", "COMPUTED"}:
        raise PhylogenyError(
            "Phylogeny requires a completed MSA. Unaligned FASTA is not an MSA.",
            "INVALID_INPUT",
        )
    rows = list(msa_result.get("rows") or [])
    n_seq = len(rows)
    if n_seq < MIN_SEQUENCES:
        raise PhylogenyError(
            f"Phylogeny requires at least {MIN_SEQUENCES} aligned sequences.",
            "INVALID_INPUT",
        )
    if n_seq > MAX_SEQUENCES:
        raise PhylogenyError(
            f"Phylogeny is limited to {MAX_SEQUENCES} sequences (RESOURCE_LIMIT).",
            "RESOURCE_LIMIT",
        )
    molecule = str(msa_result.get("molecule") or "").upper()
    if molecule not in {"DNA", "RNA", "PROTEIN"}:
        raise PhylogenyError("MSA molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")
    lengths = {len(str(row.get("aligned") or "")) for row in rows}
    if len(lengths) != 1:
        raise PhylogenyError(
            "Aligned sequences must share one alignment length.",
            "INVALID_INPUT",
        )
    alignment_length = lengths.pop()
    if alignment_length <= 0:
        raise PhylogenyError("Alignment length must be positive.", "INVALID_INPUT")
    if alignment_length > MAX_ALIGNMENT_COLUMNS:
        raise PhylogenyError(
            f"Alignment exceeds {MAX_ALIGNMENT_COLUMNS} columns (RESOURCE_LIMIT).",
            "RESOURCE_LIMIT",
        )
    identifiers = [str(row.get("identifier") or "") for row in rows]
    hashes = [str(row.get("hash") or "") for row in rows]
    if any(not item for item in hashes):
        raise PhylogenyError("Every MSA row must carry a sequence hash.", "INVALID_INPUT")
    leaf_records, transformations = unique_leaf_labels(identifiers, hashes)
    aligned_rows = [str(row.get("aligned") or "") for row in rows]
    return {
        "molecule": molecule,
        "n_sequences": n_seq,
        "alignment_length": alignment_length,
        "leaf_records": leaf_records,
        "id_transformations": transformations,
        "identical_groups": list(msa_result.get("identical_groups") or []),
        "alignment_hash": alignment_identity_hash(msa_result),
        "aligned_rows": aligned_rows,
        "rows": rows,
        "exclusions": [],
    }


def pairwise_p_distance(
    seq_a: str,
    seq_b: str,
    molecule: str,
) -> dict:
    """p-distance em sitios comparaveis (sem gap e sem ambiguo).

    Args:
        seq_a: Sequencia alinhada A.
        seq_b: Sequencia alinhada B.
        molecule: DNA, RNA ou PROTEIN.

    Returns:
        Dict p, n_sites, n_differences, method. p e fracao, nao percentual.

    Raises:
        PhylogenyError: INVALID_INPUT se comprimentos diferirem ou n_sites==0.

    Nota biologica:
        Implementacao da p-distance classica (proporcao de sitios distintos).
        Gaps sao pairwise deletion, nao mismatch. Ambiguos nao sao sorteados.
    """
    if len(seq_a) != len(seq_b):
        raise PhylogenyError("Pairwise sequences must have equal aligned length.", "INVALID_INPUT")
    alphabet = _canonical_alphabet(molecule)
    differences = 0
    sites = 0
    for raw_a, raw_b in zip(seq_a.upper(), seq_b.upper()):
        if raw_a in GAP_CHARS or raw_b in GAP_CHARS:
            continue
        if raw_a not in alphabet or raw_b not in alphabet:
            continue
        sites += 1
        if raw_a != raw_b:
            differences += 1
    if sites == 0:
        raise PhylogenyError(
            "A pair has zero comparable sites after removing gaps and ambiguities.",
            "INVALID_INPUT",
        )
    return {
        "p": differences / sites,
        "n_sites": sites,
        "n_differences": differences,
        "method": "p-distance (pairwise deletion of gaps and non-canonical symbols)",
        "status": "COMPUTED",
    }


def jukes_cantor_distance(p: float) -> dict:
    """Correcao Jukes-Cantor 1969 para quatro estados.

    Args:
        p: p-distance em [0, 1).

    Returns:
        Dict d, p, formula, status.

    Raises:
        PhylogenyError: INVALID_INPUT se p >= 0.75 ou p nao for finito.
            Nao ha fallback silencioso para p rotulado como JC.

    Nota biologica:
        Implementacao da formula publicada d = -3/4 ln(1 - 4p/3).
        Nao se aplica a proteinas.
    """
    if not math.isfinite(p) or p < 0:
        raise PhylogenyError("Jukes-Cantor p must be a finite non-negative fraction.", "INVALID_INPUT")
    if p >= 0.75:
        raise PhylogenyError(
            "Jukes-Cantor is undefined for p >= 0.75. HelixScope does not fall back to p-distance under the JC label.",
            "INVALID_INPUT",
        )
    distance = -0.75 * math.log(1.0 - (4.0 * p / 3.0))
    return {
        "d": distance,
        "p": p,
        "formula": "d = -3/4 ln(1 - 4p/3)",
        "citation": "Jukes and Cantor 1969",
        "status": "COMPUTED",
    }


def kimura_two_parameter_distance(seq_a: str, seq_b: str, molecule: str) -> dict:
    """Distancia Kimura 1980 (K2P) para DNA/RNA.

    Args:
        seq_a: Sequencia alinhada A.
        seq_b: Sequencia alinhada B.
        molecule: DNA ou RNA.

    Returns:
        Dict d, P (transicoes), Q (transversoes), n_sites.

    Raises:
        PhylogenyError: INVALID_INPUT se molecula for proteina, n_sites==0
            ou a formula for indefinida (1-2P-Q<=0 ou 1-2Q<=0).

    Nota biologica:
        Implementacao da formula publicada. Em RNA, U e pirimidina como T.
        Ambiguos e gaps sao omitidos; nao sao substituidos.
    """
    kind = str(molecule or "").upper()
    if kind not in {"DNA", "RNA"}:
        raise PhylogenyError("Kimura 1980 is defined for DNA/RNA, not protein.", "INVALID_INPUT")
    if len(seq_a) != len(seq_b):
        raise PhylogenyError("Pairwise sequences must have equal aligned length.", "INVALID_INPUT")
    alphabet = _canonical_alphabet(kind)
    transitions = 0
    transversions = 0
    sites = 0
    for raw_a, raw_b in zip(seq_a.upper(), seq_b.upper()):
        a_sym = "T" if kind == "RNA" and raw_a == "U" else raw_a
        b_sym = "T" if kind == "RNA" and raw_b == "U" else raw_b
        if raw_a in GAP_CHARS or raw_b in GAP_CHARS:
            continue
        src = raw_a if kind == "DNA" else raw_a
        dst = raw_b if kind == "DNA" else raw_b
        if src not in alphabet or dst not in alphabet:
            continue
        a_class = "T" if a_sym == "U" else a_sym
        b_class = "T" if b_sym == "U" else b_sym
        sites += 1
        if a_class == b_class:
            continue
        a_pur = a_class in PURINES
        b_pur = b_class in PURINES
        if a_pur == b_pur:
            transitions += 1
        else:
            transversions += 1
    if sites == 0:
        raise PhylogenyError(
            "A pair has zero comparable sites for Kimura 1980.",
            "INVALID_INPUT",
        )
    p_ts = transitions / sites
    q_tv = transversions / sites
    term_a = 1.0 - 2.0 * p_ts - q_tv
    term_b = 1.0 - 2.0 * q_tv
    if term_a <= 0 or term_b <= 0:
        raise PhylogenyError(
            "Kimura 1980 is undefined for this pair (1-2P-Q<=0 or 1-2Q<=0). "
            "HelixScope does not invent a substitute distance.",
            "INVALID_INPUT",
        )
    distance = -0.5 * math.log(term_a) - 0.25 * math.log(term_b)
    return {
        "d": distance,
        "P": p_ts,
        "Q": q_tv,
        "n_sites": sites,
        "n_transitions": transitions,
        "n_transversions": transversions,
        "formula": "d = -1/2 ln(1-2P-Q) - 1/4 ln(1-2Q)",
        "citation": "Kimura 1980",
        "rna_note": (
            "RNA U is treated as a pyrimidine, equivalent to T in the K2P class."
            if kind == "RNA"
            else ""
        ),
        "status": "COMPUTED",
    }


def distance_matrix_from_msa(
    aligned_rows: Sequence[str],
    leaf_records: Sequence[Mapping[str, Any]],
    molecule: str,
    distance_model: str,
) -> dict:
    """Matriz de distancia simetrica a partir do MSA.

    Args:
        aligned_rows: Linhas alinhadas na ordem dos leaf_records.
        leaf_records: Metadados de folhas (tree_id, hash).
        molecule: DNA, RNA ou PROTEIN.
        distance_model: p_distance, jukes_cantor_1969, kimura_1980 ou
            blosum62_biopython.

    Returns:
        Dict names, matrix (cheia), lower_triangle (formato Biopython),
        pairwise (detalhe), method, warnings.

    Raises:
        PhylogenyError: INVALID_INPUT se o modelo nao se aplicar a molecula
            ou algum par for indefinido.

    Nota biologica:
        A matriz e distancia computacional. Nao e uma relacao evolutiva
        comprovada. JC e K2P nao sao aplicados a proteina.
    """
    n_seq = len(aligned_rows)
    if n_seq != len(leaf_records) or n_seq < MIN_SEQUENCES:
        raise PhylogenyError("Distance matrix requires matching aligned rows and leaves.", "INVALID_INPUT")
    model = str(distance_model or "").strip()
    kind = str(molecule or "").upper()
    _assert_distance_allowed(kind, model)
    names = [str(item["tree_id"]) for item in leaf_records]
    pairwise: List[dict] = []
    warnings: List[str] = []
    full: List[List[float]] = [[0.0] * n_seq for _ in range(n_seq)]
    if model == DISTANCE_BLOSUM62:
        values = _blosum62_distances(aligned_rows, names)
        for i in range(n_seq):
            for j in range(i + 1, n_seq):
                full[i][j] = values[i][j]
                full[j][i] = values[i][j]
                pairwise.append(
                    {
                        "a": names[i],
                        "b": names[j],
                        "d": values[i][j],
                        "method": DISTANCE_BLOSUM62,
                    }
                )
        method_label = (
            "Biopython DistanceCalculator blosum62 scoring-matrix distance. "
            "Not WAG, not LG, not a phylogenetic substitution model."
        )
    else:
        for i in range(n_seq):
            for j in range(i + 1, n_seq):
                pair = _pair_distance(aligned_rows[i], aligned_rows[j], kind, model)
                if int(pair.get("n_sites") or 0) < 10:
                    warnings.append(
                        f"Pair {names[i]}-{names[j]} used {pair.get('n_sites')} comparable sites."
                    )
                full[i][j] = float(pair["d"])
                full[j][i] = float(pair["d"])
                pairwise.append({"a": names[i], "b": names[j], **pair})
        method_label = {
            DISTANCE_P: "p-distance (pairwise deletion)",
            DISTANCE_JC69: "Jukes-Cantor 1969 on p-distance (DNA/RNA 4-state)",
            DISTANCE_K2P: "Kimura 1980 two-parameter (DNA/RNA)",
        }[model]
    lower: List[List[float]] = []
    for i in range(n_seq):
        row = [float(full[i][j]) for j in range(i)]
        row.append(0.0)
        lower.append(row)
    return {
        "names": names,
        "matrix": full,
        "lower_triangle": lower,
        "pairwise": pairwise,
        "method": method_label,
        "distance_model": model,
        "molecule": kind,
        "warnings": warnings,
        "status": "COMPUTED",
        "disclaimer": (
            "These values are computational distances on this alignment. "
            "They are not phylogenetic truth and not taxonomic rank."
        ),
    }


def infer_phylogeny(
    msa_result: Mapping[str, Any],
    *,
    method: str = METHOD_NJ,
    distance_model: str = DISTANCE_JC69,
    rooting: str = ROOTING_UNROOTED,
    outgroup_tree_id: str = "",
    bootstrap_replicates: int = 0,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
    run_fn: Optional[RunFn] = None,
    cache: Optional[Dict[str, dict]] = None,
    iqtree_model_mode: str = IQTREE_MODEL_USER,
    ufboot_replicates: int = 0,
    alrt_replicates: int = 0,
) -> dict:
    """Constroi uma arvore filogenetica real a partir de um MSA validado.

    Args:
        msa_result: Envelope MSA COMPLETED.
        method: neighbor_joining, upgma, fasttree_ml ou iqtree_ml.
        distance_model: Modelo de distancia para NJ/UPGMA; JC/K2P/WAG para IQ-TREE user mode.
        rooting: unrooted, midpoint, outgroup ou upgma_rooted.
        outgroup_tree_id: Folha usada como outgroup quando rooting=outgroup.
        bootstrap_replicates: 0 = sem suporte. Maximo MAX_BOOTSTRAP_REPLICATES.
            So NJ/UPGMA. Truncado explicitamente se o pedido exceder o teto.
        bootstrap_seed: Semente do resampling (reproduzivel).
        run_fn: subprocess.run injetavel para FastTree/IQ-TREE.
        cache: Cache opcional chaveado por cache_key.
        iqtree_model_mode: user (modelo mapeado) ou MFP (ModelFinder Plus).
        ufboot_replicates: 0 ou >=1000. So IQ-TREE. Nao e bootstrap de Felsenstein.
        alrt_replicates: 0 ou >=1000. So IQ-TREE. Nao e UFBoot.

    Returns:
        Envelope PhylogeneticTree com status COMPUTED ou UNAVAILABLE.

    Raises:
        PhylogenyError: categorias classificadas. Nunca devolve uma arvore
            fabricada.

    Nota biologica:
        Inferencia filogenetica dependente das sequencias, do alinhamento,
        do metodo e do modelo. Nao e a historia evolutiva verdadeira.
        UPGMA pressupoe relogio molecular. Midpoint rooting nao e raiz
        biologica. BLAST similarity nao e este resultado. FastTree e ML
        aproximado, nao IQ-TREE. UFBoot nao e bootstrap classico.
    """
    started = time.perf_counter()
    validated = validate_msa_for_phylogeny(msa_result)
    chosen_method = str(method or "").strip()
    if chosen_method not in SUPPORTED_METHODS:
        raise PhylogenyError(f"Unsupported phylogenetic method: {method}.", "INVALID_INPUT")
    molecule = validated["molecule"]
    distance = str(distance_model or "").strip()
    if chosen_method in {METHOD_NJ, METHOD_UPGMA}:
        _assert_distance_allowed(molecule, distance)
    requested_boot = int(bootstrap_replicates)
    if requested_boot < 0:
        raise PhylogenyError("Bootstrap replicates cannot be negative.", "INVALID_INPUT")
    truncated_boot = requested_boot > MAX_BOOTSTRAP_REPLICATES
    used_boot = min(requested_boot, MAX_BOOTSTRAP_REPLICATES)
    if chosen_method in {METHOD_FASTTREE, METHOD_IQTREE} and used_boot > 0:
        raise PhylogenyError(
            "HelixScope does not attach Felsenstein column-bootstrap to FastTree/IQ-TREE output. "
            "Set bootstrap_replicates=0 for ML tools, or use Neighbor-Joining/UPGMA for resampling support.",
            "INVALID_INPUT",
        )
    model_mode = str(iqtree_model_mode or IQTREE_MODEL_USER).strip()
    if model_mode not in {IQTREE_MODEL_USER, IQTREE_MODEL_MFP}:
        raise PhylogenyError(
            "iqtree_model_mode must be user or MFP (ModelFinder Plus).",
            "INVALID_INPUT",
        )
    requested_ufboot = int(ufboot_replicates)
    requested_alrt = int(alrt_replicates)
    if requested_ufboot < 0 or requested_alrt < 0:
        raise PhylogenyError("UFBoot and SH-aLRT replicates cannot be negative.", "INVALID_INPUT")
    if chosen_method != METHOD_IQTREE:
        if model_mode == IQTREE_MODEL_MFP:
            raise PhylogenyError("ModelFinder Plus is an IQ-TREE option, not FastTree/NJ/UPGMA.", "INVALID_INPUT")
        if requested_ufboot or requested_alrt:
            raise PhylogenyError(
                "UFBoot and SH-aLRT are IQ-TREE branch tests. They are not attached to NJ, UPGMA or FastTree.",
                "INVALID_INPUT",
            )
    if requested_ufboot and requested_ufboot < MIN_UFBOOT_REPLICATES:
        raise PhylogenyError(
            f"IQ-TREE UFBoot requires 0 or at least {MIN_UFBOOT_REPLICATES} replicates. "
            "A smaller value is not silently raised.",
            "INVALID_INPUT",
        )
    if requested_alrt and requested_alrt < MIN_ALRT_REPLICATES:
        raise PhylogenyError(
            f"IQ-TREE SH-aLRT requires 0 or at least {MIN_ALRT_REPLICATES} replicates. "
            "A smaller value is not silently raised.",
            "INVALID_INPUT",
        )
    truncated_ufboot = requested_ufboot > MAX_UFBOOT_REPLICATES
    truncated_alrt = requested_alrt > MAX_ALRT_REPLICATES
    used_ufboot = min(requested_ufboot, MAX_UFBOOT_REPLICATES) if requested_ufboot else 0
    used_alrt = min(requested_alrt, MAX_ALRT_REPLICATES) if requested_alrt else 0
    rooting_method = str(rooting or ROOTING_UNROOTED).strip()
    if chosen_method == METHOD_UPGMA and rooting_method == ROOTING_UNROOTED:
        rooting_method = ROOTING_UPGMA
    if chosen_method == METHOD_NJ and rooting_method == ROOTING_UPGMA:
        raise PhylogenyError(
            "UPGMA rooting label is reserved for UPGMA trees.",
            "INVALID_INPUT",
        )
    ml_distance_token = ""
    if chosen_method == METHOD_IQTREE:
        ml_distance_token = f"{model_mode}:{distance}:{used_ufboot}:{used_alrt}"
    elif chosen_method == METHOD_FASTTREE:
        ml_distance_token = "fasttree"
    key = cache_key(
        alignment_hash=validated["alignment_hash"],
        method=chosen_method,
        distance_model=distance if chosen_method in {METHOD_NJ, METHOD_UPGMA} else ml_distance_token,
        rooting=rooting_method,
        outgroup_tree_id=outgroup_tree_id,
        bootstrap_replicates=used_boot,
        bootstrap_seed=int(bootstrap_seed) if used_boot else 0,
        tool_version=_method_tool_version(chosen_method),
    )
    if cache is not None and key in cache:
        cached = dict(cache[key])
        cached["cache_status"] = "cached"
        return cached

    if chosen_method in {METHOD_FASTTREE, METHOD_IQTREE}:
        bio_tree, ml_meta = _infer_ml_tree(
            validated=validated,
            method=chosen_method,
            molecule=molecule,
            distance_model=distance,
            run_fn=run_fn,
            iqtree_model_mode=model_mode,
            ufboot_replicates=used_ufboot,
            alrt_replicates=used_alrt,
        )
        distance_payload = {
            "method": ml_meta.get("model"),
            "distance_model": "",
            "status": "UNAVAILABLE",
            "disclaimer": (
                "Maximum-likelihood tools do not use the HelixScope NJ distance matrix."
            ),
            "matrix": [],
            "names": [item["tree_id"] for item in validated["leaf_records"]],
            "warnings": [],
        }
        tool_name = ml_meta["tool"]
        tool_version = ml_meta["version"]
        algorithm = ml_meta["algorithm"]
    else:
        distance_payload = distance_matrix_from_msa(
            validated["aligned_rows"],
            validated["leaf_records"],
            molecule,
            distance,
        )
        bio_tree = _distance_tree(distance_payload, chosen_method)
        ml_meta = {}
        tool_name = "Biopython DistanceTreeConstructor"
        tool_version = BIOPYTHON_VERSION or ""
        algorithm = (
            "Saitou and Nei 1987 Neighbor-Joining"
            if chosen_method == METHOD_NJ
            else "Sokal and Michener 1958 / Sneath and Sokal UPGMA"
        )

    rooted_tree, rooting_record = _apply_rooting(
        bio_tree,
        method=chosen_method,
        rooting=rooting_method,
        outgroup_tree_id=outgroup_tree_id,
        leaf_records=validated["leaf_records"],
    )
    support_payload = {
        "status": "UNAVAILABLE",
        "method": "",
        "n_replicates": 0,
        "n_replicates_requested": requested_boot,
        "truncated": truncated_boot,
        "seed": None,
        "values": {},
        "methods": [],
        "disclaimer": "No support values. Absence is N/A, not 0.",
    }
    if used_boot > 0:
        support_payload = _bootstrap_support(
            validated=validated,
            distance_model=distance,
            tree_method=chosen_method,
            original=rooted_tree,
            n_replicates=used_boot,
            seed=int(bootstrap_seed),
            truncated=truncated_boot,
        )
        support_payload["methods"] = [
            {
                "method": SUPPORT_FELSENSTEIN,
                "n_replicates": used_boot,
                "n_replicates_requested": requested_boot,
                "truncated": truncated_boot,
            }
        ]
        _assign_clade_support(rooted_tree, support_payload["values"])
    elif chosen_method == METHOD_IQTREE and (used_ufboot or used_alrt):
        _annotate_iqtree_support_from_labels(
            rooted_tree,
            ufboot=bool(used_ufboot),
            alrt=bool(used_alrt),
        )
        support_payload = _iqtree_support_payload(
            ufboot_replicates=used_ufboot,
            ufboot_requested=requested_ufboot,
            ufboot_truncated=truncated_ufboot,
            alrt_replicates=used_alrt,
            alrt_requested=requested_alrt,
            alrt_truncated=truncated_alrt,
        )

    serialized = serialize_bio_tree(
        rooted_tree,
        leaf_records=validated["leaf_records"],
        msa_result=msa_result,
        support_values=support_payload.get("values") or {},
    )
    _validate_serialized_tree(serialized, validated["leaf_records"])
    newick = export_newick(rooted_tree)
    parse_newick(newick)
    tree_digest = tree_hash(
        newick=newick,
        method=chosen_method,
        distance_model=distance if chosen_method in {METHOD_NJ, METHOD_UPGMA} else str(ml_meta.get("model") or ""),
        alignment_hash=validated["alignment_hash"],
        rooting=rooting_record["method"],
    )
    layouts = {
        "rectangular": layout_phylogram(serialized, mode="rectangular"),
        "radial": layout_phylogram(serialized, mode="radial"),
    }
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    envelope = provenance.analysis_envelope(
        module="PHYLOGENY",
        payload={
            "method": chosen_method,
            "distance_model": distance if chosen_method in {METHOD_NJ, METHOD_UPGMA} else "",
            "n_leaves": serialized["n_leaves"],
            "alignment_hash": validated["alignment_hash"],
            "tree_hash": tree_digest,
        },
        status="COMPUTED",
        algorithm=algorithm,
        parameters={
            "method": chosen_method,
            "distance_model": distance if chosen_method in {METHOD_NJ, METHOD_UPGMA} else str(ml_meta.get("model") or ""),
            "rooting": rooting_record["method"],
            "outgroup_tree_id": rooting_record.get("outgroup_tree_id") or "",
            "bootstrap_replicates": used_boot,
            "bootstrap_replicates_requested": requested_boot,
            "bootstrap_truncated": truncated_boot,
            "bootstrap_seed": int(bootstrap_seed) if used_boot else 0,
            "iqtree_model_mode": model_mode if chosen_method == METHOD_IQTREE else "",
            "ufboot_replicates": used_ufboot,
            "ufboot_replicates_requested": requested_ufboot,
            "ufboot_truncated": truncated_ufboot,
            "alrt_replicates": used_alrt,
            "alrt_replicates_requested": requested_alrt,
            "alrt_truncated": truncated_alrt,
            "selected_model": str(ml_meta.get("selected_model") or ml_meta.get("model") or ""),
            "model_criterion": str(ml_meta.get("model_criterion") or ""),
        },
        source=tool_name,
        input_identifier=validated["alignment_hash"],
    )
    envelope.update(serialized)
    envelope["status"] = "COMPUTED"
    envelope["kind"] = "phylogenetic_inference"
    envelope["not_taxonomic_tree"] = True
    envelope["not_true_evolutionary_history"] = True
    envelope["method"] = chosen_method
    envelope["method_label"] = _method_label(chosen_method)
    envelope["model"] = (
        distance
        if chosen_method in {METHOD_NJ, METHOD_UPGMA}
        else str(ml_meta.get("model") or "")
    )
    envelope["model_source"] = (
        distance_payload.get("method")
        if chosen_method in {METHOD_NJ, METHOD_UPGMA}
        else str(ml_meta.get("model_source") or "")
    )
    envelope["tool"] = tool_name
    envelope["tool_version"] = tool_version
    envelope["molecule"] = molecule
    envelope["alignment_hash"] = validated["alignment_hash"]
    envelope["tree_hash"] = tree_digest
    envelope["msa_n_sequences"] = validated["n_sequences"]
    envelope["msa_alignment_length"] = validated["alignment_length"]
    envelope["id_transformations"] = validated["id_transformations"]
    envelope["identical_groups"] = validated["identical_groups"]
    envelope["exclusions"] = validated["exclusions"]
    envelope["distance"] = {
        "method": distance_payload.get("method"),
        "distance_model": distance_payload.get("distance_model"),
        "matrix": distance_payload.get("matrix"),
        "names": distance_payload.get("names"),
        "warnings": distance_payload.get("warnings") or [],
        "disclaimer": distance_payload.get("disclaimer"),
        "status": distance_payload.get("status"),
    }
    envelope["rooting"] = rooting_record
    envelope["support"] = support_payload
    envelope["newick"] = newick
    envelope["layouts"] = layouts
    envelope["cache_key"] = key
    envelope["cache_status"] = "live"
    envelope["elapsed_ms"] = elapsed_ms
    envelope["retrieved_at_utc"] = provenance.utc_now()
    envelope["disclaimer"] = (
        "Phylogenetic inference from this MSA, method and model. "
        "Not a true evolutionary tree, not a taxonomic classification, "
        "and not BLAST similarity. Branch lengths are algorithm output. "
        "Support is N/A when no test was run, never 0. Felsenstein bootstrap, "
        "UFBoot and SH-aLRT are distinct methods and are not averaged."
    )
    envelope["selected_model"] = str(ml_meta.get("selected_model") or "")
    envelope["model_criterion"] = str(ml_meta.get("model_criterion") or "")
    envelope["iqtree_model_mode"] = model_mode if chosen_method == METHOD_IQTREE else ""
    if run_fn is None and chosen_method in {METHOD_IQTREE, METHOD_FASTTREE}:
        from . import engine_validation

        engine_validation.record_live_validation(
            "IQ-TREE" if chosen_method == METHOD_IQTREE else "FastTree",
            ok=True,
            version=str(ml_meta.get("version") or tool_version or ""),
            details={
                "n_leaves": serialized["n_leaves"],
                "model": str(ml_meta.get("model") or ""),
                "selected_model": str(ml_meta.get("selected_model") or ""),
                "method": chosen_method,
                "elapsed_ms": elapsed_ms,
            },
        )
    if cache is not None:
        cache[key] = envelope
    return envelope


def cache_key(
    *,
    alignment_hash: str,
    method: str,
    distance_model: str,
    rooting: str,
    outgroup_tree_id: str,
    bootstrap_replicates: int,
    bootstrap_seed: int,
    tool_version: str,
) -> str:
    """Chave de cache da arvore. Taxonomia nao entra.

    Args:
        alignment_hash: Hash do MSA.
        method: Metodo filogenetico.
        distance_model: Modelo de distancia ou vazio.
        rooting: Metodo de rooting.
        outgroup_tree_id: Outgroup ou vazio.
        bootstrap_replicates: Replicates usados.
        bootstrap_seed: Semente.
        tool_version: Versao detectada do software.

    Returns:
        Hex SHA-256.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [
            str(alignment_hash),
            str(method),
            str(distance_model),
            str(rooting),
            str(outgroup_tree_id),
            str(int(bootstrap_replicates)),
            str(int(bootstrap_seed)),
            str(tool_version),
            provenance.HELIXSCOPE_VERSION,
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def tree_hash(
    *,
    newick: str,
    method: str,
    distance_model: str,
    alignment_hash: str,
    rooting: str,
) -> str:
    """Hash deterministico da representacao da arvore.

    Args:
        newick: Newick exportado.
        method: Metodo.
        distance_model: Modelo.
        alignment_hash: Hash do MSA.
        rooting: Rooting.

    Returns:
        Hex SHA-256.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [newick.strip(), method, distance_model, alignment_hash, rooting]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def export_newick(tree: Tree) -> str:
    """Exporta a arvore Biopython exatamente como Newick.

    Args:
        tree: Bio.Phylo Tree.

    Returns:
        String Newick com ponto e virgula final.

    Raises:
        PhylogenyError: PARSING_ERROR se a escrita falhar.
    """
    handle = io.StringIO()
    try:
        Phylo.write(tree, handle, "newick")
    except Exception as exc:
        raise PhylogenyError(f"Failed to export Newick: {exc}", "PARSING_ERROR") from exc
    text = handle.getvalue().strip()
    if not text.endswith(";"):
        text += ";"
    if len(text) > MAX_NEWICK_CHARS:
        raise PhylogenyError("Newick exceeds the export size limit.", "RESOURCE_LIMIT")
    return text


def parse_newick(text: str) -> Tree:
    """Parseia Newick. Falha e PARSING_ERROR, nunca um placeholder.

    Args:
        text: String Newick.

    Returns:
        Bio.Phylo Tree.

    Raises:
        PhylogenyError: PARSING_ERROR ou RESOURCE_LIMIT.
    """
    raw = str(text or "").strip()
    if not raw:
        raise PhylogenyError("Newick text is empty.", "PARSING_ERROR")
    if len(raw) > MAX_NEWICK_CHARS:
        raise PhylogenyError("Newick exceeds the parse size limit.", "RESOURCE_LIMIT")
    try:
        tree = Phylo.read(io.StringIO(raw), "newick")
    except Exception as exc:
        raise PhylogenyError(f"Malformed Newick: {exc}", "PARSING_ERROR") from exc
    if tree is None:
        raise PhylogenyError("Malformed Newick: parser returned no tree.", "PARSING_ERROR")
    return tree


def newick_roundtrip(text: str) -> dict:
    """Parse -> export -> parse. Preserva o conjunto de folhas.

    Args:
        text: Newick de entrada.

    Returns:
        Dict input_leaves, output_leaves, newick.

    Raises:
        PhylogenyError: PARSING_ERROR se a estrutura nao round-tripar.
    """
    first = parse_newick(text)
    exported = export_newick(first)
    second = parse_newick(exported)
    first_leaves = _terminal_names(first)
    second_leaves = _terminal_names(second)
    if first_leaves != second_leaves:
        raise PhylogenyError(
            "Newick round-trip changed the leaf set.",
            "PARSING_ERROR",
        )
    return {
        "input_leaves": sorted(first_leaves),
        "output_leaves": sorted(second_leaves),
        "newick": exported,
        "status": "COMPUTED",
    }


def serialize_bio_tree(
    tree: Tree,
    *,
    leaf_records: Sequence[Mapping[str, Any]],
    msa_result: Mapping[str, Any],
    support_values: Optional[Mapping[str, float]] = None,
) -> dict:
    """Converte Bio.Phylo Tree no objeto PhylogeneticTree do HelixScope.

    Args:
        tree: Arvore Biopython ja enraizada ou nao conforme o metodo.
        leaf_records: Mapeamento tree_id -> hash e ID original.
        msa_result: Envelope MSA para metadados de accession/organismo.
        support_values: Mapa node_id -> bootstrap percent, se houver.

    Returns:
        Dict nodes, edges, leaves, n_leaves, n_nodes, branch_lengths_present.

    Raises:
        PhylogenyError: PARSING_ERROR se uma folha nao corresponder ao MSA.
    """
    by_tree_id = {str(item["tree_id"]): dict(item) for item in leaf_records}
    rows = list(msa_result.get("rows") or [])
    row_by_hash = {str(row.get("hash") or ""): row for row in rows}
    nodes: List[dict] = []
    edges: List[dict] = []
    leaves: List[dict] = []
    counter = {"n": 0}

    def next_id() -> str:
        counter["n"] += 1
        return f"n{counter['n']}"

    def visit(clade: Clade, parent_id: Optional[str]) -> str:
        this_id = next_id()
        is_leaf = bool(clade.is_terminal())
        label = str(clade.name or "") if is_leaf else str(clade.name or this_id)
        length = clade.branch_length
        length_value: Optional[float]
        if length is None:
            length_value = None
        else:
            length_value = float(length)
        support = None
        support_by_method: Dict[str, float] = {}
        if support_values and this_id in support_values:
            support = float(support_values[this_id])
            support_by_method[SUPPORT_FELSENSTEIN] = support
        elif getattr(clade, "iqtree_ufboot", None) is not None or getattr(clade, "iqtree_sh_alrt", None) is not None:
            uf_val = getattr(clade, "iqtree_ufboot", None)
            sh_val = getattr(clade, "iqtree_sh_alrt", None)
            if sh_val is not None:
                support_by_method[SUPPORT_SH_ALRT] = float(sh_val)
            if uf_val is not None:
                support_by_method[SUPPORT_UFBOOT] = float(uf_val)
            if len(support_by_method) == 1:
                support = next(iter(support_by_method.values()))
        elif getattr(clade, "confidence", None) is not None and not is_leaf:
            try:
                support = float(clade.confidence)
            except (TypeError, ValueError):
                support = None
        node = {
            "id": this_id,
            "label": label,
            "is_leaf": is_leaf,
            "branch_length": length_value,
            "support": support,
            "support_by_method": support_by_method,
            "support_method": (
                next(iter(support_by_method.keys())) if len(support_by_method) == 1 else ""
            ),
        }
        nodes.append(node)
        if parent_id is not None:
            edges.append(
                {
                    "parent": parent_id,
                    "child": this_id,
                    "length": length_value,
                }
            )
        if is_leaf:
            record = by_tree_id.get(label)
            if record is None:
                raise PhylogenyError(
                    f"Tree leaf '{label}' is not in the MSA leaf map.",
                    "PARSING_ERROR",
                )
            row = row_by_hash.get(str(record["sequence_hash"])) or {}
            organism = str(row.get("organism") or "").strip()
            source = str(row.get("source") or "").strip()
            leaves.append(
                {
                    "node_id": this_id,
                    "tree_id": label,
                    "original_id": record["original_id"],
                    "sequence_hash": record["sequence_hash"],
                    "index": record["index"],
                    "accession": str(row.get("accession") or ""),
                    "version": str(row.get("version") or ""),
                    "organism_declared": organism,
                    "organism_display": organism if organism else "unknown / user-provided",
                    "metadata_source": source or "user-provided",
                    "taxonomy": None,
                }
            )
        for child in list(clade.clades or []):
            visit(child, this_id)
        return this_id

    root_id = visit(tree.root, None)
    lengths = [edge["length"] for edge in edges]
    present = any(item is not None for item in lengths)
    negatives = [
        item for item in lengths if isinstance(item, float) and math.isfinite(item) and item < 0
    ]
    return {
        "nodes": nodes,
        "edges": edges,
        "leaves": leaves,
        "root_id": root_id,
        "n_nodes": len(nodes),
        "n_leaves": len(leaves),
        "n_edges": len(edges),
        "branch_lengths_present": present,
        "negative_branch_lengths": len(negatives),
        "rooted_flag": bool(getattr(tree, "rooted", False)),
    }


def layout_phylogram(serialized: Mapping[str, Any], mode: str = "rectangular") -> dict:
    """Coordenadas de desenho a partir dos edges reais. Layout nao muda topologia.

    Args:
        serialized: Saida de serialize_bio_tree.
        mode: rectangular ou radial.

    Returns:
        Dict nodes (id, x, y, ...), edges (x0,y0,x1,y1,length), scale_bar.

    Raises:
        PhylogenyError: INVALID_INPUT se faltar um no de um edge.

    Nota biologica:
        O layout e visual. Radial/circular nao altera a inferencia.
        Escala so existe se houver branch lengths reais.
    """
    chosen = str(mode or "rectangular").strip().lower()
    if chosen not in {"rectangular", "radial"}:
        raise PhylogenyError("Tree layout must be rectangular or radial.", "INVALID_INPUT")
    nodes = {str(item["id"]): dict(item) for item in list(serialized.get("nodes") or [])}
    edges = list(serialized.get("edges") or [])
    children: Dict[str, List[str]] = {key: [] for key in nodes}
    for edge in edges:
        parent = str(edge["parent"])
        child = str(edge["child"])
        if parent not in nodes or child not in nodes:
            raise PhylogenyError("Tree edge references a missing node.", "PARSING_ERROR")
        children[parent].append(child)
    leaf_order: List[str] = []

    def collect(node_id: str) -> None:
        kids = children.get(node_id) or []
        if not kids:
            leaf_order.append(node_id)
            return
        for child_id in kids:
            collect(child_id)

    root_id = str(serialized.get("root_id") or "")
    collect(root_id)
    y_of: Dict[str, float] = {}
    for index, leaf_id in enumerate(leaf_order):
        y_of[leaf_id] = float(index)

    def place_y(node_id: str) -> float:
        if node_id in y_of:
            return y_of[node_id]
        kids = children.get(node_id) or []
        values = [place_y(child_id) for child_id in kids]
        y_of[node_id] = sum(values) / len(values) if values else 0.0
        return y_of[node_id]

    place_y(root_id)
    x_of: Dict[str, float] = {root_id: 0.0}
    drawn_edges: List[dict] = []
    for edge in edges:
        parent = str(edge["parent"])
        child = str(edge["child"])
        length = edge.get("length")
        delta = 0.0 if length is None else float(length)
        x_of[child] = x_of[parent] + delta
        drawn_edges.append(
            {
                "parent": parent,
                "child": child,
                "length": None if length is None else float(length),
                "x0": x_of[parent],
                "y0": y_of[parent],
                "x1": x_of[child],
                "y1": y_of[child],
            }
        )
    max_x = max(x_of.values()) if x_of else 0.0
    placed: List[dict] = []
    for node_id, node in nodes.items():
        x_val = x_of[node_id]
        y_val = y_of[node_id]
        if chosen == "radial":
            angle = 0.0 if not leaf_order else (2.0 * math.pi * y_val / max(len(leaf_order), 1))
            px = x_val * math.cos(angle)
            py = x_val * math.sin(angle)
        else:
            px, py = x_val, y_val
        placed.append(
            {
                "id": node_id,
                "label": node.get("label"),
                "is_leaf": node.get("is_leaf"),
                "support": node.get("support"),
                "branch_length": node.get("branch_length"),
                "x": px,
                "y": py,
                "distance_from_root": x_val,
            }
        )
    if chosen == "radial":
        for edge in drawn_edges:
            parent_node = next(item for item in placed if item["id"] == edge["parent"])
            child_node = next(item for item in placed if item["id"] == edge["child"])
            edge["x0"] = parent_node["x"]
            edge["y0"] = parent_node["y"]
            edge["x1"] = child_node["x"]
            edge["y1"] = child_node["y"]
    scale = None
    if bool(serialized.get("branch_lengths_present")) and max_x > 0:
        magnitude = 10 ** math.floor(math.log10(max_x / 4.0)) if max_x / 4.0 > 0 else max_x
        if magnitude <= 0:
            magnitude = max_x
        scale = {
            "length": magnitude,
            "label": f"{magnitude:g} substitutions/site (algorithm units)",
            "present": True,
        }
    else:
        scale = {"length": None, "label": "", "present": False}
    return {
        "mode": chosen,
        "nodes": placed,
        "edges": drawn_edges,
        "scale_bar": scale,
        "topology_unchanged": True,
    }


def scientific_report(tree_result: Mapping[str, Any]) -> dict:
    """Relatorio reconstruivel da inferencia.

    Args:
        tree_result: Envelope de infer_phylogeny.

    Returns:
        Dict JSON-serializavel com metodo, modelo, MSA, rooting, suporte,
        exclusoes, taxonomia anexada e limitacoes.

    Raises:
        PhylogenyError: INVALID_INPUT se o envelope estiver vazio.
    """
    if not tree_result:
        raise PhylogenyError("No phylogenetic result to report.", "INVALID_INPUT")
    leaves = []
    for leaf in list(tree_result.get("leaves") or []):
        taxonomy = leaf.get("taxonomy") or {}
        leaves.append(
            {
                "tree_id": leaf.get("tree_id"),
                "original_id": leaf.get("original_id"),
                "sequence_hash": leaf.get("sequence_hash"),
                "accession": leaf.get("accession") or "",
                "version": leaf.get("version") or "",
                "organism_declared": leaf.get("organism_declared") or "",
                "taxon_id": (taxonomy or {}).get("taxon_id") or "",
                "scientific_name": (taxonomy or {}).get("scientific_name") or "",
                "taxonomy_source": (taxonomy or {}).get("source") or "",
            }
        )
    support = tree_result.get("support") or {}
    return {
        "kind": "phylogenetic_inference",
        "status": tree_result.get("status"),
        "method": tree_result.get("method"),
        "method_label": tree_result.get("method_label"),
        "model": tree_result.get("model"),
        "model_source": tree_result.get("model_source"),
        "tool": tree_result.get("tool"),
        "tool_version": tree_result.get("tool_version") or "",
        "molecule": tree_result.get("molecule"),
        "alignment_hash": tree_result.get("alignment_hash"),
        "tree_hash": tree_result.get("tree_hash"),
        "newick": tree_result.get("newick"),
        "rooting": tree_result.get("rooting"),
        "branch_lengths_present": tree_result.get("branch_lengths_present"),
        "negative_branch_lengths": tree_result.get("negative_branch_lengths"),
        "support": {
            "status": support.get("status"),
            "method": support.get("method") or "",
            "n_replicates": support.get("n_replicates"),
            "truncated": support.get("truncated"),
        },
        "distance_method": (tree_result.get("distance") or {}).get("method"),
        "id_transformations": tree_result.get("id_transformations") or [],
        "exclusions": tree_result.get("exclusions") or [],
        "identical_groups": tree_result.get("identical_groups") or [],
        "leaves": leaves,
        "elapsed_ms": tree_result.get("elapsed_ms"),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "limitations": [
            "Phylogenetic inference, not true evolutionary history.",
            "Not a taxonomic tree.",
            "Sampling, alignment, method and model dependent.",
            "Distance is computational, not phylogenetic proof.",
            "Midpoint rooting is not a biological root.",
            "UPGMA assumes an ultrametric / molecular-clock process.",
            "BLAST hits used as MSA members are not evolutionary proof.",
        ],
        "disclaimer": tree_result.get("disclaimer"),
        "retrieved_at_utc": tree_result.get("retrieved_at_utc"),
    }


def export_leaf_csv(tree_result: Mapping[str, Any]) -> str:
    """CSV de metadados das folhas. Nao e a topologia.

    Args:
        tree_result: Envelope da arvore.

    Returns:
        Texto CSV.

    Raises:
        PhylogenyError: INVALID_INPUT se nao houver folhas.
    """
    leaves = list(tree_result.get("leaves") or [])
    if not leaves:
        raise PhylogenyError("No leaves to export.", "INVALID_INPUT")
    lines = [
        "tree_id,original_id,sequence_hash,accession,version,organism_declared,"
        "taxon_id,scientific_name,rank,taxonomy_source"
    ]
    for leaf in leaves:
        taxonomy = leaf.get("taxonomy") or {}
        cells = [
            leaf.get("tree_id"),
            leaf.get("original_id"),
            leaf.get("sequence_hash"),
            leaf.get("accession"),
            leaf.get("version"),
            leaf.get("organism_declared"),
            (taxonomy or {}).get("taxon_id"),
            (taxonomy or {}).get("scientific_name"),
            (taxonomy or {}).get("rank"),
            (taxonomy or {}).get("source"),
        ]
        lines.append(",".join(_csv_cell(item) for item in cells))
    return "\n".join(lines) + "\n"


def leaf_by_sequence_hash(tree_result: Mapping[str, Any], sequence_hash: str) -> Optional[dict]:
    """Localiza a folha cujo hash coincide com a sequencia.

    Args:
        tree_result: Envelope da arvore.
        sequence_hash: Hash SHA-256 da sequencia.

    Returns:
        Dict da folha ou None. Multiplas folhas do mesmo hash devolvem a primeira.

    Raises:
        Nenhum.
    """
    digest = str(sequence_hash or "")
    for leaf in list(tree_result.get("leaves") or []):
        if str(leaf.get("sequence_hash") or "") == digest:
            return dict(leaf)
    return None


def leaf_by_tree_id(tree_result: Mapping[str, Any], tree_id: str) -> Optional[dict]:
    """Localiza a folha pelo ID da arvore.

    Args:
        tree_result: Envelope da arvore.
        tree_id: Identificador Newick.

    Returns:
        Dict da folha ou None.

    Raises:
        Nenhum.
    """
    label = str(tree_id or "")
    for leaf in list(tree_result.get("leaves") or []):
        if str(leaf.get("tree_id") or "") == label:
            return dict(leaf)
    return None


def _canonical_alphabet(molecule: str) -> str:
    kind = str(molecule or "").upper()
    if kind == "DNA":
        return DNA_CANONICAL
    if kind == "RNA":
        return RNA_CANONICAL
    if kind == "PROTEIN":
        return PROTEIN_CANONICAL
    raise PhylogenyError("Molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")


def _assert_distance_allowed(molecule: str, model: str) -> None:
    kind = str(molecule or "").upper()
    chosen = str(model or "").strip()
    if kind in {"DNA", "RNA"}:
        if chosen not in NUCLEOTIDE_DISTANCES:
            raise PhylogenyError(
                f"Distance {chosen} is not valid for {kind}. JC/K2P/p-distance only.",
                "INVALID_INPUT",
            )
        return
    if kind == "PROTEIN":
        if chosen not in PROTEIN_DISTANCES:
            raise PhylogenyError(
                "Protein phylogeny does not use DNA models (JC/K2P).",
                "INVALID_INPUT",
            )
        return
    raise PhylogenyError("Molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")


def _pair_distance(seq_a: str, seq_b: str, molecule: str, model: str) -> dict:
    p_info = pairwise_p_distance(seq_a, seq_b, molecule)
    if model == DISTANCE_P:
        return {"d": float(p_info["p"]), **p_info}
    if model == DISTANCE_JC69:
        jc = jukes_cantor_distance(float(p_info["p"]))
        return {"d": float(jc["d"]), **p_info, **jc}
    if model == DISTANCE_K2P:
        k2p = kimura_two_parameter_distance(seq_a, seq_b, molecule)
        return {**k2p, **p_info, "d": float(k2p["d"])}
    raise PhylogenyError(f"Unknown distance model: {model}.", "INVALID_INPUT")


def _blosum62_distances(aligned_rows: Sequence[str], names: Sequence[str]) -> List[List[float]]:
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = [
        SeqRecord(Seq(row), id=name, name=name, description="")
        for row, name in zip(aligned_rows, names)
    ]
    alignment = MultipleSeqAlignment(records)
    calculator = DistanceCalculator("blosum62")
    dm = calculator.get_distance(alignment)
    n_seq = len(names)
    full = [[0.0] * n_seq for _ in range(n_seq)]
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            if i == j:
                continue
            full[i][j] = float(dm[left, right])
    return full


def _distance_tree(distance_payload: Mapping[str, Any], method: str) -> Tree:
    names = list(distance_payload["names"])
    lower = list(distance_payload["lower_triangle"])
    try:
        matrix = DistanceMatrix(names, lower)
    except Exception as exc:
        raise PhylogenyError(f"Invalid distance matrix: {exc}", "INVALID_INPUT") from exc
    constructor = DistanceTreeConstructor()
    try:
        if method == METHOD_NJ:
            return constructor.nj(matrix)
        if method == METHOD_UPGMA:
            return constructor.upgma(matrix)
    except Exception as exc:
        raise PhylogenyError(f"Tree construction failed: {exc}", "INTERNAL") from exc
    raise PhylogenyError(f"Unknown distance method: {method}.", "INVALID_INPUT")


def _apply_rooting(
    tree: Tree,
    *,
    method: str,
    rooting: str,
    outgroup_tree_id: str,
    leaf_records: Sequence[Mapping[str, Any]],
) -> Tuple[Tree, dict]:
    working = copy.deepcopy(tree)
    leaves = {str(item["tree_id"]) for item in leaf_records}
    if rooting == ROOTING_OUTGROUP:
        target = str(outgroup_tree_id or "").strip()
        if target not in leaves:
            raise PhylogenyError(
                "Outgroup must be a leaf identifier present in this tree.",
                "INVALID_INPUT",
            )
        try:
            working.root_with_outgroup(target)
        except Exception as exc:
            raise PhylogenyError(f"Outgroup rooting failed: {exc}", "INVALID_INPUT") from exc
        working.rooted = True
        return working, {
            "method": ROOTING_OUTGROUP,
            "rooted": True,
            "biological_root": False,
            "outgroup_tree_id": target,
            "note": "Outgroup-rooted display. The outgroup choice is a user parameter, not an inferred ancestor.",
        }
    if rooting == ROOTING_MIDPOINT:
        try:
            working.root_at_midpoint()
        except Exception as exc:
            raise PhylogenyError(f"Midpoint rooting failed: {exc}", "INVALID_INPUT") from exc
        working.rooted = True
        return working, {
            "method": ROOTING_MIDPOINT,
            "rooted": True,
            "biological_root": False,
            "outgroup_tree_id": "",
            "note": "Midpoint rooted. This is not a biological root.",
        }
    if method == METHOD_UPGMA or rooting == ROOTING_UPGMA:
        working.rooted = True
        return working, {
            "method": ROOTING_UPGMA,
            "rooted": True,
            "biological_root": False,
            "outgroup_tree_id": "",
            "ultrametric_assumption": True,
            "note": "UPGMA produces a rooted ultrametric tree and assumes a molecular clock. Not a universal method.",
        }
    working.rooted = False
    return working, {
        "method": ROOTING_UNROOTED,
        "rooted": False,
        "biological_root": False,
        "outgroup_tree_id": "",
        "display_root": "algorithm_node",
        "note": "Unrooted inference. The drawing root is an algorithm node, not an evolutionary root.",
    }


def _bootstrap_support(
    *,
    validated: Mapping[str, Any],
    distance_model: str,
    tree_method: str,
    original: Tree,
    n_replicates: int,
    seed: int,
    truncated: bool,
) -> dict:
    all_leaves = _terminal_names(original)
    original_splits = _node_splits(original, all_leaves)
    rng = random.Random(int(seed))
    counts: Dict[str, int] = {node_id: 0 for node_id in original_splits}
    aligned = list(validated["aligned_rows"])
    n_col = len(aligned[0])
    for _ in range(int(n_replicates)):
        indices = [rng.randrange(n_col) for _ in range(n_col)]
        resampled = ["".join(row[index] for index in indices) for row in aligned]
        distance_payload = distance_matrix_from_msa(
            resampled,
            validated["leaf_records"],
            str(validated["molecule"]),
            distance_model,
        )
        replica = _distance_tree(distance_payload, tree_method)
        replica_splits = _unrooted_splits(replica, all_leaves)
        for node_id, split in original_splits.items():
            if split in replica_splits:
                counts[node_id] += 1
    values = {
        node_id: 100.0 * counts[node_id] / float(n_replicates)
        for node_id in counts
    }
    return {
        "status": "COMPUTED",
        "method": "Felsenstein 1985 bootstrap (column resampling)",
        "n_replicates": int(n_replicates),
        "n_replicates_requested": int(n_replicates),
        "truncated": bool(truncated),
        "seed": int(seed),
        "values": values,
        "label": "bootstrap support",
        "not_confidence": True,
        "methods": [],
        "disclaimer": (
            "Bootstrap support (percent of resampled trees containing the split). "
            "Not a posterior probability and not named confidence."
            + (" Replicate count was truncated to the resource limit." if truncated else "")
        ),
    }


def _assign_clade_support(tree: Tree, values: Mapping[str, float]) -> None:
    counter = {"n": 0}

    def next_id() -> str:
        counter["n"] += 1
        return f"n{counter['n']}"

    def walk(clade: Clade) -> None:
        this_id = next_id()
        if this_id in values and not clade.is_terminal():
            clade.confidence = float(values[this_id])
        for child in list(clade.clades or []):
            walk(child)

    walk(tree.root)


def _terminal_names(tree: Tree) -> frozenset[str]:
    return frozenset(str(clade.name) for clade in tree.get_terminals() if clade.name)


def _leaf_set(clade: Clade) -> frozenset[str]:
    if clade.is_terminal():
        return frozenset([str(clade.name or "")])
    names: set[str] = set()
    for child in list(clade.clades or []):
        names |= set(_leaf_set(child))
    return frozenset(names)


def _unrooted_splits(tree: Tree, all_leaves: frozenset[str]) -> set[frozenset[frozenset[str]]]:
    splits: set[frozenset[frozenset[str]]] = set()

    def walk(clade: Clade) -> frozenset[str]:
        if clade.is_terminal():
            return _leaf_set(clade)
        inside: set[str] = set()
        for child in list(clade.clades or []):
            child_set = walk(child)
            inside |= set(child_set)
            if 0 < len(child_set) < len(all_leaves):
                other = all_leaves - child_set
                splits.add(frozenset((child_set, other)))
        return frozenset(inside)

    walk(tree.root)
    return splits


def _node_splits(tree: Tree, all_leaves: frozenset[str]) -> Dict[str, frozenset[frozenset[str]]]:
    mapping: Dict[str, frozenset[frozenset[str]]] = {}
    counter = {"n": 0}

    def next_id() -> str:
        counter["n"] += 1
        return f"n{counter['n']}"

    def walk(clade: Clade) -> frozenset[str]:
        this_id = next_id()
        if clade.is_terminal():
            return _leaf_set(clade)
        inside: set[str] = set()
        for child in list(clade.clades or []):
            inside |= set(walk(child))
        child_set = frozenset(inside)
        if 0 < len(child_set) < len(all_leaves):
            other = all_leaves - child_set
            mapping[this_id] = frozenset((child_set, other))
        return child_set

    walk(tree.root)
    return mapping


def _validate_serialized_tree(
    serialized: Mapping[str, Any],
    leaf_records: Sequence[Mapping[str, Any]],
) -> None:
    leaves = list(serialized.get("leaves") or [])
    expected = [str(item["tree_id"]) for item in leaf_records]
    got = [str(item["tree_id"]) for item in leaves]
    if sorted(got) != sorted(expected):
        raise PhylogenyError(
            "Tree leaves do not match the MSA sequences. No silent exclusions.",
            "PARSING_ERROR",
        )
    if len(got) != len(set(got)):
        raise PhylogenyError("Tree contains duplicated leaf identifiers.", "PARSING_ERROR")
    if len(leaves) != len(leaf_records):
        raise PhylogenyError("Leaf count does not match the MSA.", "PARSING_ERROR")
    for edge in list(serialized.get("edges") or []):
        length = edge.get("length")
        if length is None:
            continue
        if not isinstance(length, (int, float)) or not math.isfinite(float(length)):
            raise PhylogenyError("A branch length is not a finite number.", "PARSING_ERROR")


def _method_label(method: str) -> str:
    return {
        METHOD_NJ: "Neighbor-Joining (Saitou and Nei 1987)",
        METHOD_UPGMA: "UPGMA (ultrametric / molecular-clock assumption)",
        METHOD_FASTTREE: "FastTree approximate maximum likelihood",
        METHOD_IQTREE: "IQ-TREE maximum likelihood",
    }.get(method, method)


def _method_tool_version(method: str) -> str:
    if method in {METHOD_NJ, METHOD_UPGMA}:
        return BIOPYTHON_VERSION or ""
    if method == METHOD_FASTTREE:
        return str(detect_fasttree().get("version") or "")
    if method == METHOD_IQTREE:
        return str(detect_iqtree().get("version") or "")
    return ""


def _detect_ml_binary(
    *,
    names: Sequence[str],
    env_var: str,
    tool: str,
    tool_id: str,
    version_args: Sequence[str],
    version_pattern: re.Pattern[str],
) -> dict:
    executable = tool_detection.resolve_allowlisted_executable(
        names,
        extra_file_candidates=tool_paths.candidates_for(names),
        env_var=env_var,
    )
    if executable is None:
        return {
            "id": tool_id,
            "tool": tool,
            "available": False,
            "version": "",
            "version_status": "unavailable",
            "status": tool_detection.TOOL_STATUS_NOT_INSTALLED,
            "path": "",
            "reason": (
                f"{tool} executable not installed. HelixScope does not reimplement it. "
                "Neighbor-Joining and UPGMA remain available."
            ),
        }
    version = ""
    try:
        completed = subprocess.run(
            [executable, *version_args],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
        blob = f"{completed.stdout or ''}\n{completed.stderr or ''}"
        match = version_pattern.search(blob)
        if match:
            version = match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        version = ""
    return {
        "id": tool_id,
        "tool": tool,
        "available": True,
        "version": version,
        "version_status": "detected" if version else tool_detection.TOOL_STATUS_VERSION_UNAVAILABLE,
        "status": tool_detection.classify_tool_record(
            available=True, path=executable, version=version
        ),
        "path": tool_detection.sanitize_tool_path(executable),
        "reason": f"{tool} executable detected.",
    }


def _parse_iqtree_support_token(
    raw: str,
    *,
    ufboot: bool,
    alrt: bool,
) -> Tuple[Optional[float], Optional[float]]:
    text = str(raw or "").strip()
    if not IQTREE_SUPPORT_TOKEN.match(text):
        return None, None
    if "/" in text:
        left, right = text.split("/", 1)
        sh_val = float(left)
        uf_val = float(right)
        if alrt and ufboot:
            return sh_val, uf_val
        if alrt and not ufboot:
            return sh_val, None
        if ufboot and not alrt:
            return None, uf_val
        return None, None
    single = float(text)
    if alrt and ufboot:
        return None, None
    if alrt:
        return single, None
    if ufboot:
        return None, single
    return None, None


def _annotate_iqtree_support_from_labels(
    tree: Tree,
    *,
    ufboot: bool,
    alrt: bool,
) -> None:
    """Le suporte IQ-TREE em labels internos (ex: 95.4/98). Nao mistura metodos."""
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        raw = str(clade.name or "")
        if not IQTREE_SUPPORT_TOKEN.match(raw):
            continue
        sh_val, uf_val = _parse_iqtree_support_token(raw, ufboot=ufboot, alrt=alrt)
        clade.name = None
        clade.iqtree_sh_alrt = sh_val
        clade.iqtree_ufboot = uf_val


def _iqtree_support_payload(
    *,
    ufboot_replicates: int,
    ufboot_requested: int,
    ufboot_truncated: bool,
    alrt_replicates: int,
    alrt_requested: int,
    alrt_truncated: bool,
) -> dict:
    methods: List[dict] = []
    if alrt_replicates:
        methods.append(
            {
                "method": SUPPORT_SH_ALRT,
                "label": "SH-aLRT",
                "n_replicates": int(alrt_replicates),
                "n_replicates_requested": int(alrt_requested),
                "truncated": bool(alrt_truncated),
                "not_ufboot": True,
                "not_felsenstein_bootstrap": True,
            }
        )
    if ufboot_replicates:
        methods.append(
            {
                "method": SUPPORT_UFBOOT,
                "label": "Ultrafast bootstrap (UFBoot)",
                "n_replicates": int(ufboot_replicates),
                "n_replicates_requested": int(ufboot_requested),
                "truncated": bool(ufboot_truncated),
                "not_felsenstein_bootstrap": True,
                "not_sh_alrt": True,
            }
        )
    labels = [str(item["label"]) for item in methods]
    return {
        "status": "COMPUTED" if methods else "UNAVAILABLE",
        "method": "; ".join(str(item["method"]) for item in methods),
        "label": " + ".join(labels) if labels else "",
        "n_replicates": 0,
        "n_replicates_requested": 0,
        "truncated": bool(ufboot_truncated or alrt_truncated),
        "seed": None,
        "values": {},
        "methods": methods,
        "disclaimer": (
            "IQ-TREE branch tests as requested. UFBoot is not Felsenstein 1985 "
            "bootstrap. SH-aLRT is not UFBoot. Values are N/A when the tool did "
            "not label an internal node, never 0."
        ),
    }


def _parse_iqtree_report(text: str) -> dict:
    blob = str(text or "")
    if len(blob) > MAX_IQTREE_REPORT_CHARS:
        raise PhylogenyError("IQ-TREE report exceeds the size limit.", "RESOURCE_LIMIT")
    selected = ""
    criterion = ""
    for match in BEST_FIT_MODEL_RE.finditer(blob):
        criterion = match.group(1).upper()
        selected = match.group(2)
        if criterion == "BIC":
            break
    return {
        "selected_model": selected,
        "model_criterion": criterion,
    }


def _infer_ml_tree(
    *,
    validated: Mapping[str, Any],
    method: str,
    molecule: str,
    distance_model: str,
    run_fn: Optional[RunFn],
    iqtree_model_mode: str = IQTREE_MODEL_USER,
    ufboot_replicates: int = 0,
    alrt_replicates: int = 0,
) -> Tuple[Tree, dict]:
    if method == METHOD_FASTTREE:
        info = detect_fasttree()
        if not info["available"]:
            raise PhylogenyError(info["reason"], "TOOL_NOT_INSTALLED")
        executable = tool_detection.resolve_allowlisted_executable(
            FASTTREE_NAMES,
            extra_file_candidates=tool_paths.candidates_for(FASTTREE_NAMES),
            env_var=FASTTREE_ENV,
        )
        if not executable:
            raise PhylogenyError(info["reason"], "TOOL_NOT_INSTALLED")
        newick, version = _run_fasttree(
            executable=executable,
            aligned_rows=validated["aligned_rows"],
            leaf_records=validated["leaf_records"],
            molecule=molecule,
            run_fn=run_fn,
        )
        tree = parse_newick(newick)
        model = "FastTree default approximate ML (nucleotide JC-like / protein CAT as implemented by FastTree)"
        if molecule == "PROTEIN" and distance_model in {DISTANCE_JC69, DISTANCE_K2P}:
            raise PhylogenyError("DNA distance models cannot be applied to protein FastTree runs.", "INVALID_INPUT")
        return tree, {
            "tool": "FastTree",
            "version": version or str(info.get("version") or ""),
            "algorithm": "FastTree approximate maximum likelihood",
            "model": model,
            "model_source": "FastTree executable; not the HelixScope NJ distance matrix",
        }
    if method == METHOD_IQTREE:
        info = detect_iqtree()
        if not info["available"]:
            raise PhylogenyError(info["reason"], "TOOL_NOT_INSTALLED")
        executable = tool_detection.resolve_allowlisted_executable(
            IQTREE_NAMES,
            extra_file_candidates=tool_paths.candidates_for(IQTREE_NAMES),
            env_var=IQTREE_ENV,
        )
        if not executable:
            raise PhylogenyError(info["reason"], "TOOL_NOT_INSTALLED")
        if iqtree_model_mode == IQTREE_MODEL_MFP:
            iq_model = IQTREE_MODEL_MFP
        else:
            iq_model = _iqtree_model(molecule, distance_model)
        newick, version, report_meta = _run_iqtree(
            executable=executable,
            aligned_rows=validated["aligned_rows"],
            leaf_records=validated["leaf_records"],
            iq_model=iq_model,
            ufboot_replicates=int(ufboot_replicates),
            alrt_replicates=int(alrt_replicates),
            run_fn=run_fn,
        )
        tree = parse_newick(newick)
        selected = str(report_meta.get("selected_model") or "")
        criterion = str(report_meta.get("model_criterion") or "")
        displayed_model = selected or iq_model
        model_source = (
            "IQ-TREE ModelFinder Plus (-m MFP); selected model from the .iqtree report"
            if iqtree_model_mode == IQTREE_MODEL_MFP
            else "IQ-TREE -m as requested; version detected from the binary"
        )
        if iqtree_model_mode == IQTREE_MODEL_MFP and not selected:
            model_source += (
                " (report lacked a Best-fit model line; requested keyword remains MFP)"
            )
        return tree, {
            "tool": "IQ-TREE",
            "version": version or str(info.get("version") or ""),
            "algorithm": "IQ-TREE maximum likelihood",
            "model": displayed_model,
            "selected_model": selected,
            "model_criterion": criterion,
            "model_source": model_source,
        }
    raise PhylogenyError(f"Unknown ML method: {method}.", "INVALID_INPUT")


def _iqtree_model(molecule: str, distance_model: str) -> str:
    kind = str(molecule or "").upper()
    if kind == "PROTEIN":
        if distance_model in {DISTANCE_JC69, DISTANCE_K2P}:
            raise PhylogenyError("IQ-TREE protein runs cannot use DNA models.", "INVALID_INPUT")
        return "WAG"
    if kind in {"DNA", "RNA"}:
        mapping = {
            DISTANCE_P: "JC",
            DISTANCE_JC69: "JC",
            DISTANCE_K2P: "K2P",
            "": "JC",
        }
        if distance_model == DISTANCE_P:
            raise PhylogenyError(
                "IQ-TREE does not run under p-distance. Choose Jukes-Cantor or Kimura 1980.",
                "INVALID_INPUT",
            )
        if distance_model not in mapping:
            raise PhylogenyError(
                f"IQ-TREE nucleotide model mapping does not include {distance_model}.",
                "INVALID_INPUT",
            )
        return mapping[distance_model]
    raise PhylogenyError("IQ-TREE molecule must be DNA, RNA or PROTEIN.", "INVALID_INPUT")


def _write_aligned_fasta(
    directory: str,
    aligned_rows: Sequence[str],
    leaf_records: Sequence[Mapping[str, Any]],
) -> str:
    path = os.path.join(directory, "alignment.fa")
    lines: List[str] = []
    for row, record in zip(aligned_rows, leaf_records):
        lines.append(f">{record['tree_id']}")
        lines.append(str(row))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _run_fasttree(
    *,
    executable: str,
    aligned_rows: Sequence[str],
    leaf_records: Sequence[Mapping[str, Any]],
    molecule: str,
    run_fn: Optional[RunFn],
) -> Tuple[str, str]:
    runner = run_fn or subprocess.run
    with tempfile.TemporaryDirectory(prefix="helixscope_fasttree_") as tmp:
        fasta = _write_aligned_fasta(tmp, aligned_rows, leaf_records)
        argv = [executable]
        if molecule in {"DNA", "RNA"}:
            argv.append("-nt")
        argv.append(fasta)
        try:
            completed = runner(
                argv,
                capture_output=True,
                text=True,
                timeout=ML_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmp,
            )
        except subprocess.TimeoutExpired as exc:
            raise PhylogenyError(
                f"FastTree exceeded {ML_TIMEOUT_S:.0f}s and was stopped.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise PhylogenyError(f"FastTree could not be executed: {exc}", "INTERNAL") from exc
        stdout = str(getattr(completed, "stdout", "") or "")
        stderr = str(getattr(completed, "stderr", "") or "")
        if len(stdout) > MAX_ML_OUTPUT_CHARS:
            raise PhylogenyError("FastTree output exceeds the size limit.", "RESOURCE_LIMIT")
        if completed.returncode not in {0} or not stdout.strip():
            snippet = stderr.strip()[:300] or "no Newick on stdout"
            raise PhylogenyError(f"FastTree did not return a tree: {snippet}", "PARSING_ERROR")
        version = ""
        match = re.search(r"FastTree\s+Version\s+([0-9][0-9.\w-]*)", stderr, re.I)
        if match:
            version = match.group(1)
        return stdout.strip(), version


def _run_iqtree(
    *,
    executable: str,
    aligned_rows: Sequence[str],
    leaf_records: Sequence[Mapping[str, Any]],
    iq_model: str,
    ufboot_replicates: int = 0,
    alrt_replicates: int = 0,
    run_fn: Optional[RunFn],
) -> Tuple[str, str, dict]:
    runner = run_fn or subprocess.run
    with tempfile.TemporaryDirectory(prefix="helixscope_iqtree_") as tmp:
        fasta = _write_aligned_fasta(tmp, aligned_rows, leaf_records)
        prefix = os.path.join(tmp, "run")
        argv = [
            executable,
            "-s",
            fasta,
            "-m",
            iq_model,
            "-nt",
            "1",
            "--prefix",
            prefix,
            "-redo",
            "-quiet",
        ]
        if ufboot_replicates:
            argv.extend(["-B", str(int(ufboot_replicates))])
        if alrt_replicates:
            argv.extend(["-alrt", str(int(alrt_replicates))])
        try:
            completed = runner(
                argv,
                capture_output=True,
                text=True,
                timeout=ML_TIMEOUT_S,
                check=False,
                shell=False,
                cwd=tmp,
            )
        except subprocess.TimeoutExpired as exc:
            raise PhylogenyError(
                f"IQ-TREE exceeded {ML_TIMEOUT_S:.0f}s and was stopped.",
                "TIMEOUT",
            ) from exc
        except OSError as exc:
            raise PhylogenyError(f"IQ-TREE could not be executed: {exc}", "TOOL_ERROR") from exc
        returncode = int(getattr(completed, "returncode", 1) or 0)
        stderr = str(getattr(completed, "stderr", "") or "")
        if len(stderr) > MAX_ML_OUTPUT_CHARS:
            raise PhylogenyError("IQ-TREE stderr exceeds the size limit.", "RESOURCE_LIMIT")
        treefile = prefix + ".treefile"
        if not os.path.isfile(treefile):
            category = "TOOL_ERROR" if returncode not in {0} else "PARSING_ERROR"
            raise PhylogenyError(
                f"IQ-TREE did not write a treefile. {stderr.strip()[:300] or 'No tree output.'}",
                category,
            )
        with open(treefile, "r", encoding="utf-8") as handle:
            newick = handle.read()
        if len(newick) > MAX_ML_OUTPUT_CHARS:
            raise PhylogenyError("IQ-TREE treefile exceeds the size limit.", "RESOURCE_LIMIT")
        report_meta = {"selected_model": "", "model_criterion": ""}
        report_path = prefix + ".iqtree"
        if os.path.isfile(report_path):
            with open(report_path, "r", encoding="utf-8", errors="replace") as handle:
                report_meta = _parse_iqtree_report(handle.read())
        version = ""
        blob = f"{getattr(completed, 'stdout', '')}\n{stderr}"
        match = re.search(r"version\s+([0-9][0-9.\w-]*)", blob, re.I)
        if match:
            version = match.group(1)
        return newick.strip(), version, report_meta


def _csv_cell(value: object) -> str:
    text = "" if value is None else str(value)
    if any(char in text for char in {",", '"', "\n"}):
        return '"' + text.replace('"', '""') + '"'
    return text
