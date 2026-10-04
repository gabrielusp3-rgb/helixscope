"""Verificacoes de sanidade sobre resultados cientificos ja calculados.

Nao substitui o algoritmo original: apenas recusa inconsistencias numericas ou
coordenadas impossiveis antes da visualizacao. Nenhuma funcao inventa um valor
corrigido em silencio.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence


def percent_is_valid(value: float) -> bool:
    """Aceita NaN (desconhecido) ou um percentual no intervalo [0, 100].

    Args:
        value: Percentual ou NaN.

    Returns:
        True se o valor puder ser exibido; False se estiver fora de [0, 100].

    Raises:
        Nenhum.
    """
    if isinstance(value, float) and math.isnan(value):
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return 0.0 <= number <= 100.0


def composition_counts_match_length(counts: Mapping[str, int], length: int) -> bool:
    """Confere se a soma das contagens iguala o comprimento da sequencia.

    Args:
        counts: Mapeamento simbolo -> contagem.
        length: Comprimento esperado.

    Returns:
        True quando as somas coincidem e length e nao negativo.

    Raises:
        Nenhum.
    """
    if length < 0:
        return False
    total = 0
    for value in counts.values():
        try:
            total += int(value)
        except (TypeError, ValueError):
            return False
    return total == length


def span_is_valid(start: int, end: int, sequence_length: int, *, base: int = 0) -> bool:
    """Valida um intervalo [start, end) ou 1-based inclusivo segundo base.

    Args:
        start: Inicio do intervalo.
        end: Fim do intervalo (exclusivo se base=0).
        sequence_length: Comprimento da sequencia de referencia.
        base: 0 para semiaberto 0-based; 1 para inclusivo 1-based.

    Returns:
        True se o intervalo cabe na sequencia.

    Raises:
        ValueError: Se base nao for 0 nem 1.
    """
    if base not in {0, 1}:
        raise ValueError("base must be 0 or 1.")
    if sequence_length <= 0:
        return False
    if base == 0:
        return 0 <= start < end <= sequence_length
    return 1 <= start <= end <= sequence_length


def to_1based_inclusive(start_0: int, end_0_exclusive: int) -> tuple[int, int]:
    """Converte [start, end) 0-based em start/end 1-based inclusivos.

    Args:
        start_0: Inicio 0-based.
        end_0_exclusive: Fim exclusivo 0-based.

    Returns:
        Par (start_1, end_1) inclusivo.

    Raises:
        ValueError: Se o intervalo for vazio ou invertido.
    """
    if end_0_exclusive <= start_0:
        raise ValueError("Half-open interval must have end > start.")
    return start_0 + 1, end_0_exclusive


def to_0based_half_open(start_1: int, end_1_inclusive: int) -> tuple[int, int]:
    """Converte start/end 1-based inclusivos em [start, end) 0-based.

    Args:
        start_1: Inicio 1-based inclusivo.
        end_1_inclusive: Fim 1-based inclusivo.

    Returns:
        Par (start_0, end_0_exclusive).

    Raises:
        ValueError: Se as coordenadas 1-based forem invalidas.
    """
    if start_1 < 1 or end_1_inclusive < start_1:
        raise ValueError("1-based coordinates must satisfy 1 <= start <= end.")
    return start_1 - 1, end_1_inclusive


def alignment_rows_are_consistent(aligned_seq1: str, aligned_seq2: str) -> bool:
    """Exige duas linhas de alinhamento do mesmo comprimento.

    Args:
        aligned_seq1: Primeira sequencia alinhada, com gaps.
        aligned_seq2: Segunda sequencia alinhada, com gaps.

    Returns:
        True se ambas forem nao vazias e tiverem o mesmo comprimento.

    Raises:
        Nenhum.
    """
    if not aligned_seq1 or not aligned_seq2:
        return False
    return len(aligned_seq1) == len(aligned_seq2)


def motif_hit_matches_sequence(sequence: str, start: int, end: int, match: str) -> bool:
    """Confere se o trecho reportado coincide com a sequencia na posicao dada.

    Args:
        sequence: Sequencia na mesma capitalizacao do hit.
        start: Inicio 0-based.
        end: Fim exclusivo.
        match: Trecho reportado na fita +.

    Returns:
        True se [start:end] == match. Hits da fita - nao sao confrontados com
        o trecho direto, porque o match e lido na fita reversa.

    Raises:
        Nenhum.
    """
    if start < 0 or end > len(sequence) or start >= end:
        return False
    return sequence[start:end] == match


def crispr_guide_in_sequence(sequence: str, guide: Mapping[str, object]) -> bool:
    """Confere se um candidato CRISPR cabe na sequencia e cita o protoespacador.

    Args:
        sequence: DNA alvo ja normalizado.
        guide: Dict com ao menos guide_sequence, position (1-based) e strand.

    Returns:
        True se a posicao 1-based e o comprimento do guia cabem na sequencia e
        o protoespacador aparece na fita correspondente. Nao valida o PAM.

    Raises:
        Nenhum.

    Nota biologica:
        Um guia cuja coordenada cai fora da sequencia analisada nao pode ser
        apresentado como sítio daquela entrada.
    """
    spacer = str(guide.get("guide_sequence") or "")
    try:
        position = int(guide.get("position"))
    except (TypeError, ValueError):
        return False
    strand = str(guide.get("strand") or "+")
    if not spacer or position < 1:
        return False
    end = position + len(spacer) - 1
    if end > len(sequence):
        return False
    window = sequence[position - 1 : end]
    if strand == "+":
        return window == spacer
    try:
        from . import dna_analysis

        return dna_analysis.reverse_complement(window) == spacer
    except ValueError:
        return False


def chart_xy_lengths_match(x: Sequence[object], y: Sequence[object]) -> bool:
    """Recusa graficos cujo eixo X e Y tem comprimentos diferentes.

    Args:
        x: Coordenadas do eixo X.
        y: Valores do eixo Y.

    Returns:
        True se len(x) == len(y) e ambos sao nao vazios.

    Raises:
        Nenhum.
    """
    return len(x) == len(y) and len(x) > 0


def msa_row_length_matches(aligned: str, alignment_length: int) -> bool:
    """Confere se uma linha alinhada tem exatamente alignment_length colunas.

    Args:
        aligned: Sequencia alinhada, com gaps.
        alignment_length: Comprimento esperado.

    Returns:
        True se len(aligned) == alignment_length e alignment_length > 0.

    Raises:
        Nenhum.
    """
    return alignment_length > 0 and len(aligned) == alignment_length


def conservation_score_is_valid(value: float) -> bool:
    """Aceita NaN (coluna sem residuos) ou um score no intervalo [0, 1].

    Args:
        value: Score de conservacao Shannon ou NaN.

    Returns:
        True se puder ser exibido.

    Raises:
        Nenhum.
    """
    if isinstance(value, float) and math.isnan(value):
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return 0.0 <= number <= 1.0


def rna_structure_length_matches(sequence: str, structure: str) -> bool:
    """Confere se a notacao de colchetes tem o mesmo comprimento da RNA.

    Args:
        sequence: RNA submetida ao folding.
        structure: Dot-bracket.

    Returns:
        True se os comprimentos coincidem e nao sao vazios.

    Raises:
        Nenhum.
    """
    return bool(sequence) and bool(structure) and len(sequence) == len(structure)


def rna_parentheses_are_balanced(structure: str) -> bool:
    """Confere balanceamento de '(' e ')' na notacao MFE, sem pseudoknots.

    Args:
        structure: Dot-bracket.

    Returns:
        True se o saldo nunca fica negativo e termina em zero.

    Raises:
        Nenhum.
    """
    depth = 0
    for char in str(structure or ""):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def atom_coordinate_is_finite(x: object, y: object, z: object) -> bool:
    """Confere se as coordenadas cartesianas sao numeros finitos.

    Args:
        x, y, z: Coordenadas de um atomo.

    Returns:
        True se os tres valores sao float/int finitos. Zero e valido.

    Raises:
        Nenhum.
    """
    try:
        coords = (float(x), float(y), float(z))
    except (TypeError, ValueError):
        return False
    return not any(math.isnan(value) or math.isinf(value) for value in coords)


def mapping_status_is_declared(status: object) -> bool:
    """Aceita apenas os estados de mapping estrutural declarados.

    Args:
        status: Rotulo proposto.

    Returns:
        True se status for EXACT, ALIGNED, PARTIAL, BEST_EFFORT, UNMAPPED
        ou UNCERTAIN.

    Raises:
        Nenhum.
    """
    allowed = {
        "EXACT",
        "ALIGNED",
        "PARTIAL",
        "BEST_EFFORT",
        "UNMAPPED",
        "UNCERTAIN",
    }
    return str(status or "").strip().upper() in allowed


def protein_mapping_covers_query(n_mapped: int, query_length: int) -> bool:
    """Confere se o mapping cobre todas as posicoes da query 0-based.

    Args:
        n_mapped: Numero de linhas de mapping (uma por residuo da query).
        query_length: Comprimento da proteina de analise.

    Returns:
        True se n_mapped == query_length e query_length > 0.

    Raises:
        Nenhum.
    """
    return query_length > 0 and n_mapped == query_length


def mfe_kcal_mol_is_valid(value: float) -> bool:
    """Aceita um MFE finito em kcal/mol. Recusa NaN e infinitos.

    Args:
        value: Energia reportada pela ferramenta.

    Returns:
        True se value e um float/int finito. Zero e valido.

    Raises:
        Nenhum.
    """
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    if isinstance(number, float) and (math.isnan(number) or math.isinf(number)):
        return False
    return True


def blast_hsp_is_valid(
    hsp: Mapping[str, object],
    *,
    query_len: int = 0,
    hit_len: int = 0,
) -> bool:
    """Confere invariantes de um HSP NCBI antes da renderizacao.

    Args:
        hsp: Dict com identities, alignment_length, coordenadas e alinhamento.
        query_len: Comprimento da query (0 ignora o teto da query).
        hit_len: Comprimento do sujeito (0 ignora o teto do sujeito).

    Returns:
        True se identities <= alignment_length, percentuais em [0, 100],
        coordenadas 1-based e linhas de alinhamento consistentes.

    Raises:
        Nenhum.

    Nota biologica:
        Estes numeros vem do XML BlastOutput do NCBI. A funcao nao calcula
        E-value nem bit score.
    """
    try:
        identities = int(hsp.get("identities"))
        align_len = int(hsp.get("alignment_length"))
        query_from = int(hsp.get("query_from"))
        query_to = int(hsp.get("query_to"))
        hit_from = int(hsp.get("hit_from"))
        hit_to = int(hsp.get("hit_to"))
    except (TypeError, ValueError):
        return False
    if align_len <= 0 or identities < 0 or identities > align_len:
        return False
    positives = hsp.get("positives")
    if positives is not None:
        try:
            if int(positives) < 0 or int(positives) > align_len:
                return False
        except (TypeError, ValueError):
            return False
    gaps = hsp.get("gaps")
    if gaps is not None:
        try:
            if int(gaps) < 0 or int(gaps) > align_len:
                return False
        except (TypeError, ValueError):
            return False
    identity_pct = hsp.get("identity_pct")
    coverage = hsp.get("query_coverage_pct")
    if identity_pct is not None and not percent_is_valid(float(identity_pct)):
        return False
    if coverage is not None:
        try:
            cov_value = float(coverage)
        except (TypeError, ValueError):
            return False
        if not percent_is_valid(cov_value):
            return False
    if query_from <= 0 or query_to <= 0 or hit_from <= 0 or hit_to <= 0:
        return False
    if query_len > 0 and max(query_from, query_to) > query_len:
        return False
    if hit_len > 0 and max(hit_from, hit_to) > hit_len:
        return False
    qseq = str(hsp.get("qseq") or "")
    hseq = str(hsp.get("hseq") or "")
    midline = str(hsp.get("midline") or "")
    if qseq and hseq and not alignment_rows_are_consistent(qseq, hseq):
        return False
    if qseq and len(qseq) != align_len:
        return False
    if midline and qseq and len(midline) != len(qseq):
        return False
    accession = str(hsp.get("accession") or "")
    hit_id = str(hsp.get("hit_id") or "")
    if "accession" in hsp and not accession and not hit_id:
        return False
    return True


def protein_composition_matches_length(counts: Mapping[str, int], length: int) -> bool:
    """Confere soma de aminoacidos contra o comprimento da cadeia.

    Args:
        counts: Mapeamento residuo -> contagem.
        length: Comprimento esperado.

    Returns:
        True quando as somas coincidem e length e nao negativo.
    """
    return composition_counts_match_length(counts, length)


def make_span(
    *,
    start: int,
    end: int,
    sequence: str,
    strand: str = "+",
    source: str = "computed",
    label: str = "",
    kind: str = "interval",
    status: str = "COMPUTED",
) -> dict:
    """Representacao de um intervalo para sincronizacao futura 2D/3D.

    Args:
        start: Inicio 0-based.
        end: Fim exclusivo.
        sequence: Sequencia de referencia.
        strand: "+" ou "-".
        source: Origem do intervalo.
        label: Rotulo opcional (Predicted ORF, motif, etc.).
        kind: Tipo do intervalo (interval, orf, motif, guide).
        status: Estado de evidencia; nao promove ORF a gene.

    Returns:
        Dict com start, end, strand, length, subsequence (fita +), source,
        label, type e status.

    Raises:
        ValueError: Se o intervalo for invalido.
    """
    if strand not in {"+", "-"}:
        raise ValueError('strand must be "+" or "-".')
    if not span_is_valid(start, end, len(sequence), base=0):
        raise ValueError("Span is not valid on the reference sequence.")
    return {
        "start": start,
        "end": end,
        "strand": strand,
        "length": end - start,
        "sequence": sequence[start:end],
        "source": source,
        "label": label,
        "type": kind,
        "kind": kind,
        "status": status,
    }


def offtarget_hit_coordinates_valid(
    hit: Mapping[str, object], contig_length: int
) -> bool:
    """Confere se o intervalo do hit e do PAM cabem no contig.

    Args:
        hit: Off-target com start_0based, end_0based, pam_start_0based,
            pam_end_0based e strand.
        contig_length: Comprimento do contig.

    Returns:
        True se os intervalos 0-based semiabertos forem validos e a fita for
        "+" ou "-".

    Raises:
        Nenhum.
    """
    strand = str(hit.get("strand") or "")
    if strand not in {"+", "-"}:
        return False
    try:
        start = int(hit.get("start_0based"))
        end = int(hit.get("end_0based"))
        pam_start = int(hit.get("pam_start_0based"))
        pam_end = int(hit.get("pam_end_0based"))
    except (TypeError, ValueError):
        return False
    return span_is_valid(start, end, contig_length, base=0) and span_is_valid(
        pam_start, pam_end, contig_length, base=0
    )


def phylogenetic_leaf_set_matches_msa(
    tree_leaf_ids: Sequence[object],
    msa_tree_ids: Sequence[object],
) -> bool:
    """Confere se as folhas da arvore sao exatamente as sequencias do MSA.

    Args:
        tree_leaf_ids: IDs das folhas serializadas.
        msa_tree_ids: IDs esperados do MSA (ja normalizados).

    Returns:
        True se os conjuntos coincidirem e nao houver duplicata.

    Raises:
        Nenhum.
    """
    tree_ids = [str(item) for item in tree_leaf_ids]
    msa_ids = [str(item) for item in msa_tree_ids]
    if len(tree_ids) != len(set(tree_ids)):
        return False
    if len(msa_ids) != len(set(msa_ids)):
        return False
    return sorted(tree_ids) == sorted(msa_ids)


def branch_length_is_valid(value: object) -> bool:
    """Aceita None (ausente) ou um numero finito. Nao substitui ausente por 1.0.

    Args:
        value: Comprimento de ramo ou None.

    Returns:
        True se puder ser exibido. False se for NaN/inf ou nao numerico.

    Raises:
        Nenhum.
    """
    if value is None:
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number)


def bootstrap_support_is_valid(value: object, *, replicates: int) -> bool:
    """Suporte bootstrap: N/A se nao houver replicates; senao [0, 100].

    Args:
        value: Percentual ou None.
        replicates: Numero de replicates executados.

    Returns:
        True se a semantica estiver correta. 0 so e valido com replicates>0.

    Raises:
        Nenhum.
    """
    if int(replicates or 0) <= 0:
        return value is None
    if value is None:
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    if not math.isfinite(number):
        return False
    return 0.0 <= number <= 100.0
