"""Alinhamento de sequencias e matriz de dotplot.

Implementa alinhamento global (Needleman-Wunsch) e local (Smith-Waterman) por
meio do PairwiseAligner do Biopython, alem de uma matriz de dotplot binaria.

Nenhuma funcao aqui importa Streamlit; o Biopython e importado de forma tardia.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from . import dna_analysis

MATCH_SCORE: float = 2.0
"""Pontuacao atribuida a um par de residuos identicos em DNA/RNA."""

MISMATCH_SCORE: float = -1.0
"""Pontuacao atribuida a um par de residuos diferentes em DNA/RNA."""

OPEN_GAP_SCORE: float = -5.0
"""Penalizacao por abrir um gap."""

EXTEND_GAP_SCORE: float = -0.5
"""Penalizacao por estender um gap ja aberto."""

PROTEIN_SUBSTITUTION_MATRIX: str = "BLOSUM62"
"""Matriz de substituicao usada apenas em alinhamentos de proteina."""

MAX_ALIGN_LENGTH: int = 8_000
"""Teto tecnico por sequencia no Needleman-Wunsch / Smith-Waterman (O(n*m))."""

MAX_DOTPLOT_LENGTH: int = 2_000
"""Teto tecnico por sequencia na matriz de dotplot (O(n*m) em memoria)."""


def _ensure_pair_length(seq1: str, seq2: str, max_length: int, task: str) -> None:
    """Recusa pares que estourariam memoria ou CPU no algoritmo pedido."""
    if len(seq1) > max_length or len(seq2) > max_length:
        raise ValueError(
            f"{task} is limited to {max_length:,} residues per sequence. "
            "Paste a shorter region of interest."
        )


def _prepare_pair(seq1: str, seq2: str) -> Tuple[str, str, str]:
    """Normaliza e exige duas sequencias do mesmo tipo molecular (interno).

    Args:
        seq1: Primeira sequencia bruta.
        seq2: Segunda sequencia bruta.

    Returns:
        Tupla (seq1 canonica, seq2 canonica, tipo molecular).

    Raises:
        ValueError: Se alguma sequencia for vazia, invalida, ou se os tipos
            moleculares forem diferentes.

    Nota biologica:
        Alinhar DNA com proteina ou DNA com RNA nao tem interpretacao de
        homologia posicao a posicao sob a mesma matriz de identidade.
    """
    info1 = dna_analysis.validate_sequence(seq1)
    info2 = dna_analysis.validate_sequence(seq2)
    if not info1["is_valid"] or not info2["is_valid"]:
        reason1 = info1.get("rejection_reason") or "sequence 1 is invalid"
        reason2 = info2.get("rejection_reason") or "sequence 2 is invalid"
        raise ValueError(
            "Both sequences must be valid canonical DNA, RNA or protein. "
            f"Sequence 1: {reason1} Sequence 2: {reason2}"
        )
    if info1["type"] != info2["type"]:
        raise ValueError(
            f"Cannot align {info1['type']} with {info2['type']}. "
            "Both sequences must be the same molecule type."
        )
    return str(info1["sequence"]), str(info2["sequence"]), str(info1["type"])


def _load_blosum62():
    """Carrega BLOSUM62 do Biopython (uso interno)."""
    try:
        from Bio.Align import substitution_matrices
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    return substitution_matrices.load(PROTEIN_SUBSTITUTION_MATRIX)


def _build_aligner(mode: str, molecule: str):
    """Cria um PairwiseAligner configurado (uso interno).

    Args:
        mode: Modo de alinhamento, "global" ou "local".
        molecule: Tipo molecular canonico ("DNA", "RNA" ou "PROTEIN").

    Returns:
        Tupla (aligner, nome da pontuacao, matriz ou None).

    Raises:
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        DNA e RNA usam pontuacao de identidade. Proteinas usam BLOSUM62, que
        pontua substituicoes observadas em blocos conservados, nao identidade
        bruta.
    """
    try:
        from Bio.Align import PairwiseAligner
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = OPEN_GAP_SCORE
    aligner.extend_gap_score = EXTEND_GAP_SCORE
    aligner.open_end_gap_score = OPEN_GAP_SCORE
    aligner.extend_end_gap_score = EXTEND_GAP_SCORE
    if molecule == "PROTEIN":
        matrix = _load_blosum62()
        aligner.substitution_matrix = matrix
        return aligner, PROTEIN_SUBSTITUTION_MATRIX, matrix
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_SCORE
    return aligner, "identity match=2 mismatch=-1", None


def _empty_alignment(
    molecule: str, scoring: str, query_length: int = 0, subject_length: int = 0
) -> Dict[str, object]:
    """Resultado vazio de alinhamento local sem subalinhamento positivo."""
    return {
        "score": 0.0,
        "aligned_seq1": "",
        "aligned_seq2": "",
        "identity_pct": 0.0,
        "ungapped_identity_pct": float("nan"),
        "similarity_pct": 0.0,
        "matches": 0,
        "mismatches": 0,
        "gaps": 0,
        "length": 0,
        "coverage_seq1_pct": 0.0,
        "coverage_seq2_pct": 0.0,
        "query_length": query_length,
        "subject_length": subject_length,
        "molecule": molecule,
        "scoring": scoring,
        "method": "Smith-Waterman",
        "status": "COMPUTED",
        "parameters": {
            "mode": "local",
            "match": MATCH_SCORE if molecule != "PROTEIN" else None,
            "mismatch": MISMATCH_SCORE if molecule != "PROTEIN" else None,
            "substitution": PROTEIN_SUBSTITUTION_MATRIX if molecule == "PROTEIN" else None,
            "open_gap": OPEN_GAP_SCORE,
            "extend_gap": EXTEND_GAP_SCORE,
        },
    }


def _summarize(
    alignment,
    score: float,
    *,
    molecule: str,
    scoring: str,
    query_length: int,
    subject_length: int,
    mode: str,
    substitution_matrix: Optional[object] = None,
) -> Dict[str, object]:
    """Resume um alinhamento em um dicionario padronizado (uso interno).

    Args:
        alignment: Objeto de alinhamento do Biopython.
        score: Pontuacao otima do alinhamento.
        molecule: Tipo molecular alinhado.
        scoring: Nome do esquema de pontuacao efetivamente usado.
        query_length: Comprimento da primeira sequencia de entrada.
        subject_length: Comprimento da segunda sequencia de entrada.
        mode: "global" ou "local".
        substitution_matrix: Matriz BLOSUM/PAM, ou None para identidade.

    Returns:
        Dicionario com score, sequencias alinhadas, identidade (incluindo gaps
        no denominador), identidade sem gaps, cobertura, similaridade,
        contagens e metadados de metodo.

    Nota biologica:
        identity_pct = matches / alignment_length, com colunas de gap no
        denominador. ungapped_identity_pct = matches / (matches + mismatches).
        coverage_seq1_pct = residuos nao-gap da linha 1 / comprimento original 1.
    """
    aligned_seq1 = str(alignment[0])
    aligned_seq2 = str(alignment[1])
    length = len(aligned_seq1)

    matches = 0
    mismatches = 0
    gaps = 0
    positives = 0
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-" or b == "-":
            gaps += 1
            continue
        if a == b:
            matches += 1
            positives += 1
            continue
        mismatches += 1
        if substitution_matrix is not None:
            try:
                if float(substitution_matrix[a, b]) > 0:
                    positives += 1
            except (KeyError, ValueError, TypeError, IndexError):
                pass
    identity_pct = round((matches / length) * 100.0, 2) if length else 0.0
    ungapped = matches + mismatches
    ungapped_identity_pct = (
        round((matches / ungapped) * 100.0, 2) if ungapped else float("nan")
    )
    if substitution_matrix is None:
        similarity_pct = identity_pct
    else:
        similarity_pct = round((positives / length) * 100.0, 2) if length else 0.0
    nongap1 = sum(1 for residue in aligned_seq1 if residue != "-")
    nongap2 = sum(1 for residue in aligned_seq2 if residue != "-")
    coverage_seq1_pct = (
        round((nongap1 / query_length) * 100.0, 2) if query_length else float("nan")
    )
    coverage_seq2_pct = (
        round((nongap2 / subject_length) * 100.0, 2) if subject_length else float("nan")
    )

    return {
        "score": round(float(score), 2),
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "identity_pct": identity_pct,
        "ungapped_identity_pct": ungapped_identity_pct,
        "similarity_pct": similarity_pct,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "length": length,
        "coverage_seq1_pct": coverage_seq1_pct,
        "coverage_seq2_pct": coverage_seq2_pct,
        "identity_definition": "matches / alignment_length including gap columns, percent",
        "ungapped_identity_definition": "matches / (matches + mismatches), percent",
        "coverage_seq1_definition": "non-gap residues in aligned seq1 / original seq1 length, percent",
        "coverage_seq2_definition": "non-gap residues in aligned seq2 / original seq2 length, percent",
        "query_length": query_length,
        "subject_length": subject_length,
        "molecule": molecule,
        "scoring": scoring,
        "method": (
            "Needleman-Wunsch" if mode == "global" else "Smith-Waterman"
        ),
        "status": "COMPUTED",
        "parameters": {
            "mode": mode,
            "match": MATCH_SCORE if molecule != "PROTEIN" else None,
            "mismatch": MISMATCH_SCORE if molecule != "PROTEIN" else None,
            "substitution": PROTEIN_SUBSTITUTION_MATRIX if molecule == "PROTEIN" else None,
            "open_gap": OPEN_GAP_SCORE,
            "extend_gap": EXTEND_GAP_SCORE,
        },
    }


def pairwise_global(seq1: str, seq2: str) -> Dict[str, object]:
    """Realiza o alinhamento global de Needleman-Wunsch entre duas sequencias.

    Usa o PairwiseAligner do Biopython em modo global com match=2, mismatch=-1,
    abertura de gap=-5.0 e extensao de gap=-0.5.

    Args:
        seq1: Primeira sequencia (sera normalizada).
        seq2: Segunda sequencia (sera normalizada).

    Returns:
        Dicionario com "score" (float), "aligned_seq1" (str), "aligned_seq2"
        (str), "identity_pct" (float), "similarity_pct" (float), "matches"
        (int), "mismatches" (int), "gaps" (int), "length" (int), "molecule"
        (str), "scoring" (str) e "status" (str).

    Raises:
        ValueError: Se alguma sequencia for vazia, invalida, se os tipos
            moleculares forem diferentes, ou se alguma exceder MAX_ALIGN_LENGTH.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        O alinhamento global (Needleman-Wunsch com gaps afins de Gotoh) e
        adequado quando se espera homologia ao longo de toda a extensao das
        sequencias, como genes ortologos completos. Gaps terminais recebem a
        mesma penalidade dos internos. DNA/RNA usam identidade (match=2,
        mismatch=-1). Proteinas usam BLOSUM62.
    """
    a, b, molecule = _prepare_pair(seq1, seq2)
    _ensure_pair_length(a, b, MAX_ALIGN_LENGTH, "Pairwise alignment")
    aligner, scoring, matrix = _build_aligner("global", molecule)
    alignments = aligner.align(a, b)
    try:
        best = alignments[0]
    except IndexError as exc:
        raise ValueError("Nenhum alinhamento global foi produzido.") from exc
    return _summarize(
        best,
        best.score,
        molecule=molecule,
        scoring=scoring,
        query_length=len(a),
        subject_length=len(b),
        mode="global",
        substitution_matrix=matrix,
    )


def pairwise_local(seq1: str, seq2: str) -> Dict[str, object]:
    """Realiza o alinhamento local de Smith-Waterman entre duas sequencias.

    Usa o PairwiseAligner do Biopython em modo local com os mesmos parametros de
    score do alinhamento global (match=2, mismatch=-1, abertura=-5.0,
    extensao=-0.5).

    Args:
        seq1: Primeira sequencia (sera normalizada).
        seq2: Segunda sequencia (sera normalizada).

    Returns:
        Dicionario com "score" (float), "aligned_seq1" (str), "aligned_seq2"
        (str), "identity_pct" (float), "similarity_pct" (float), "matches"
        (int), "mismatches" (int), "gaps" (int), "length" (int), "molecule"
        (str), "scoring" (str) e "status" (str).

    Raises:
        ValueError: Se alguma sequencia for vazia, invalida, se os tipos
            moleculares forem diferentes, ou se alguma exceder MAX_ALIGN_LENGTH.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        O alinhamento local identifica as sub-regioes de maior semelhanca, util
        para detectar dominios conservados em sequencias divergentes. DNA/RNA
        usam identidade; proteinas usam BLOSUM62.
    """
    a, b, molecule = _prepare_pair(seq1, seq2)
    _ensure_pair_length(a, b, MAX_ALIGN_LENGTH, "Pairwise alignment")
    aligner, scoring, matrix = _build_aligner("local", molecule)
    alignments = aligner.align(a, b)
    try:
        best = alignments[0]
    except IndexError:
        return _empty_alignment(molecule, scoring, len(a), len(b))
    return _summarize(
        best,
        best.score,
        molecule=molecule,
        scoring=scoring,
        query_length=len(a),
        subject_length=len(b),
        mode="local",
        substitution_matrix=matrix,
    )


def dotplot_matrix(seq1: str, seq2: str, window: int = 1) -> np.ndarray:
    """Constroi uma matriz de dotplot binaria entre duas sequencias.

    Cada celula (i, j) recebe 1 quando a subsequencia de tamanho window iniciada
    em i na primeira sequencia e identica a iniciada em j na segunda, e 0 caso
    contrario.

    Args:
        seq1: Primeira sequencia (sera normalizada).
        seq2: Segunda sequencia (sera normalizada).
        window: Tamanho da janela de comparacao; deve ser positivo e nao maior
            que o comprimento da menor sequencia.

    Returns:
        Matriz numpy de dtype uint8 com forma
        (len(seq1) - window + 1, len(seq2) - window + 1).

    Raises:
        ValueError: Se window nao for positivo, se window for maior que
            min(len(seq1), len(seq2)), ou se alguma sequencia exceder
            MAX_DOTPLOT_LENGTH.
        RuntimeError: Nunca; nao depende do Biopython.

    Nota biologica:
        O dotplot revela visualmente regioes de similaridade, repeticoes e
        inversoes entre duas sequencias; janelas maiores reduzem o ruido de
        coincidencias curtas.
    """
    a, b, _molecule = _prepare_pair(seq1, seq2)
    _ensure_pair_length(a, b, MAX_DOTPLOT_LENGTH, "Dot plot")
    if window <= 0:
        raise ValueError("window deve ser positivo.")
    if window > min(len(a), len(b)):
        raise ValueError("window nao pode ser maior que a menor sequencia.")

    rows = len(a) - window + 1
    cols = len(b) - window + 1
    if window == 1:
        a_arr = np.frombuffer(a.encode("ascii"), dtype=np.uint8)
        b_arr = np.frombuffer(b.encode("ascii"), dtype=np.uint8)
        matrix = (a_arr[:, None] == b_arr[None, :]).astype(np.uint8)
        return matrix
    matrix = np.zeros((rows, cols), dtype=np.uint8)
    a_windows = [a[i : i + window] for i in range(rows)]
    b_windows = [b[j : j + window] for j in range(cols)]
    for i, wa in enumerate(a_windows):
        for j, wb in enumerate(b_windows):
            if wa == wb:
                matrix[i, j] = 1
    return matrix


def classify_alignment_columns(aligned_seq1: str, aligned_seq2: str) -> List[str]:
    """Classifica cada coluna do alinhamento como match, mismatch ou gap.

    Args:
        aligned_seq1: Primeira sequencia alinhada.
        aligned_seq2: Segunda sequencia alinhada.

    Returns:
        Lista do mesmo comprimento, com "match", "mismatch" ou "gap".

    Raises:
        ValueError: Se os comprimentos diferirem ou alguma linha for vazia.
    """
    if not aligned_seq1 or not aligned_seq2:
        raise ValueError("Alignment is empty.")
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("Aligned sequences have different lengths.")
    classes: List[str] = []
    for residue_a, residue_b in zip(aligned_seq1, aligned_seq2):
        if residue_a == "-" or residue_b == "-":
            classes.append("gap")
        elif residue_a == residue_b:
            classes.append("match")
        else:
            classes.append("mismatch")
    return classes


def format_pairwise_export(
    result: Dict[str, object],
    name1: str = "seq1",
    name2: str = "seq2",
    width: int = 60,
) -> str:
    """Formata um alinhamento pairwise como texto verificavel (estilo CLUSTAL).

    Args:
        result: Dicionario retornado por pairwise_global ou pairwise_local.
        name1: Identificador da primeira sequencia.
        name2: Identificador da segunda sequencia.
        width: Residuos por linha; deve ser positivo.

    Returns:
        Texto com cabecalho de metodo, blocos alinhados e linha de identidade.

    Raises:
        ValueError: Se o resultado nao contiver sequencias alinhadas do mesmo
            comprimento, ou se width nao for positivo.

    Nota biologica:
        O texto e uma representacao do alinhamento ja calculado; nao realinha.
    """
    if width <= 0:
        raise ValueError("width deve ser positivo.")
    aligned1 = str(result.get("aligned_seq1") or "")
    aligned2 = str(result.get("aligned_seq2") or "")
    if not aligned1 or not aligned2:
        raise ValueError("Alignment is empty; nothing to export.")
    if len(aligned1) != len(aligned2):
        raise ValueError("Aligned sequences have different lengths.")

    match_line = []
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            match_line.append(" ")
        elif a == b:
            match_line.append("*")
        else:
            match_line.append(" ")
    matches = "".join(match_line)

    label1 = (name1 or "seq1")[:10].ljust(10)
    label2 = (name2 or "seq2")[:10].ljust(10)
    marker = " " * 10
    lines = [
        f"HelixScope pairwise alignment",
        f"molecule={result.get('molecule', 'unknown')}",
        f"scoring={result.get('scoring', 'unknown')}",
        f"status={result.get('status', 'COMPUTED')}",
        f"score={result.get('score')}",
        f"identity_pct={result.get('identity_pct')}",
        f"ungapped_identity_pct={result.get('ungapped_identity_pct')}",
        f"similarity_pct={result.get('similarity_pct')}",
        f"coverage_seq1_pct={result.get('coverage_seq1_pct')}",
        f"coverage_seq2_pct={result.get('coverage_seq2_pct')}",
        f"gaps={result.get('gaps')} length={result.get('length')}",
        "",
    ]
    for start in range(0, len(aligned1), width):
        end = start + width
        lines.append(f"{label1} {aligned1[start:end]}")
        lines.append(f"{marker} {matches[start:end]}")
        lines.append(f"{label2} {aligned2[start:end]}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def multiple_alignment_availability() -> dict:
    """Disponibilidade de MSA real; pairwise continua nao sendo MSA.

    Args:
        Nenhum.

    Returns:
        Dict com available (True se houver backend real), reason enfatizando
        que Needleman-Wunsch nao e MSA, e backends.

    Raises:
        Nenhum.

    Nota biologica:
        Alinhamento multiplo e um algoritmo proprio (Clustal Omega, MUSCLE,
        MAFFT). Simular um MSA a partir de alinhamentos pairwise seria
        fabricacao.
    """
    from . import msa

    tools = msa.tool_availability()
    return {
        "available": bool(tools.get("any_available")),
        "reason": (
            "Pairwise Needleman-Wunsch and Smith-Waterman are not MSA. "
            + str(tools.get("summary_reason") or "")
        ),
        "recommended_tools": [
            "Clustal Omega (EMBL-EBI Job Dispatcher)",
            "MAFFT",
            "MUSCLE",
        ],
        "backends": tools.get("backends"),
        "blast_is_not_msa": True,
    }


def blast_search_availability() -> dict:
    """Declara que BLAST / QBLAST nao esta implementado neste projeto.

    Args:
        Nenhum.

    Returns:
        Dicionario com "available" (False), "reason" e "recommended_tools".

    Raises:
        Nenhum.

    Nota biologica:
        BLAST e um algoritmo e um servico especificos. Um score pairwise local
        nao e um resultado BLAST. Inventar um "BLAST-like score" seria
        fabricacao.
    """
    return {
        "available": False,
        "reason": (
            "BLAST and QBLAST are not implemented. Pairwise Smith-Waterman in "
            "HelixScope is not a BLAST search and is not reported as one."
        ),
        "recommended_tools": [
            "NCBI BLAST",
            "BLAST+",
        ],
    }
