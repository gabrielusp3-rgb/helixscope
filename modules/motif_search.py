"""Busca de motivos em DNA com suporte a codigos IUPAC nas duas fitas.

Procura um padrao (incluindo codigos de ambiguidade IUPAC) na fita informada e
no seu complemento reverso, retornando as posicoes em coordenadas da fita
informada. Os dados externos (expansao IUPAC para classes de regex e mapa de
complementaridade) estao embutidos como constantes. Nenhuma funcao aqui importa
Streamlit.
"""

from __future__ import annotations

import re
from typing import Dict, List

IUPAC_TO_REGEX: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
"""Expansao de cada simbolo IUPAC para a classe de caracteres regex equivalente."""

COMPLEMENT_IUPAC: Dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
"""Mapa de complementaridade de bases com suporte a codigos IUPAC."""

MAX_MOTIF_PATTERN_LENGTH: int = 256
"""Teto tecnico do padrao de motivo; padroes maiores nao tem uso biologico aqui."""

MAX_MOTIF_HITS: int = 10_000
"""Teto de hits sobrepostos. Acima disso a busca e abortada, nao truncada."""

MAX_MOTIF_PATTERNS: int = 20
"""Numero maximo de padroes numa busca multipla."""


def _reverse_complement(seq: str) -> str:
    """Calcula o complemento reverso de uma sequencia com suporte IUPAC.

    Args:
        seq: Sequencia de DNA ja normalizada em maiusculas.

    Returns:
        Complemento reverso da sequencia.

    Raises:
        ValueError: Se a sequencia contiver simbolos sem complemento definido.

    Nota biologica:
        Motivos podem ocorrer na fita oposta; o complemento reverso permite
        localiza-los lendo a fita antisentido no sentido 5'->3'.
    """
    out: List[str] = []
    for base in reversed(seq):
        if base not in COMPLEMENT_IUPAC:
            raise ValueError(f"Simbolo sem complemento definido: {base}.")
        out.append(COMPLEMENT_IUPAC[base])
    return "".join(out)


def _pattern_to_regex(pattern: str) -> str:
    """Converte um padrao com codigos IUPAC em uma expressao regular.

    Args:
        pattern: Padrao de motivo (sera normalizado para maiusculas).

    Returns:
        String de expressao regular equivalente ao padrao.

    Raises:
        ValueError: Se o padrao for vazio ou exceder MAX_MOTIF_PATTERN_LENGTH.

    Nota biologica:
        Padroes degenerados (IUPAC) capturam a variabilidade natural de sitios
        funcionais sem exigir uma sequencia exata.
    """
    cleaned = "".join(pattern.split()).upper()
    if not cleaned:
        raise ValueError("O padrao de motivo esta vazio.")
    if len(cleaned) > MAX_MOTIF_PATTERN_LENGTH:
        raise ValueError(
            f"O padrao de motivo excede {MAX_MOTIF_PATTERN_LENGTH} simbolos."
        )
    parts: List[str] = []
    for symbol in cleaned:
        parts.append(IUPAC_TO_REGEX.get(symbol, re.escape(symbol)))
    return "".join(parts)


def find_motif(seq: str, pattern: str) -> List[dict]:
    """Localiza um motivo (com codigos IUPAC) nas duas fitas de uma sequencia.

    Procura o padrao na fita informada (strand "+") e no complemento reverso
    (strand "-"), de forma insensivel a maiusculas e capturando ocorrencias
    sobrepostas. As posicoes sao reportadas em coordenadas da fita informada.

    Args:
        seq: Sequencia de DNA (sera normalizada para maiusculas).
        pattern: Padrao de motivo, aceitando codigos de ambiguidade IUPAC.

    Returns:
        Lista de dicionarios com "start" (int, base 0), "end" (int, exclusivo),
        "strand" (str: "+" ou "-") e "match" (str, trecho correspondente lido na
        respectiva fita), ordenada por "start".

    Raises:
        ValueError: Se a sequencia ou o padrao forem vazios.

    Nota biologica:
        Motivos costumam corresponder a sitios de ligacao de fatores de
        transcricao, sitios de restricao ou elementos regulatorios, que podem
        residir em qualquer das duas fitas.
    """
    cleaned = "".join(seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia esta vazia.")
    regex = _pattern_to_regex(pattern)
    compiled = re.compile(f"(?=({regex}))", re.IGNORECASE)

    results: List[dict] = []

    def _add(hit: dict) -> None:
        if len(results) >= MAX_MOTIF_HITS:
            raise ValueError(
                f"Motif search exceeded {MAX_MOTIF_HITS:,} hits. Narrow the "
                "pattern or search a shorter region. No partial hit list is returned."
            )
        results.append(hit)

    for match in compiled.finditer(cleaned):
        text = match.group(1)
        start = match.start()
        _add(
            {
                "start": start,
                "end": start + len(text),
                "strand": "+",
                "match": text,
            }
        )

    rc = _reverse_complement(cleaned)
    n = len(cleaned)
    for match in compiled.finditer(rc):
        text = match.group(1)
        rc_start = match.start()
        forward_start = n - (rc_start + len(text))
        _add(
            {
                "start": forward_start,
                "end": forward_start + len(text),
                "strand": "-",
                "match": text,
            }
        )

    results.sort(key=lambda item: item["start"])
    return results


IUPAC_TO_REGEX_RNA: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "U": "U",
    "R": "[AG]", "Y": "[CU]", "S": "[GC]", "W": "[AU]",
    "K": "[GU]", "M": "[AC]", "B": "[CGU]", "D": "[AGU]",
    "H": "[ACU]", "V": "[ACG]", "N": "[ACGU]",
}
"""Expansao IUPAC para RNA (U no lugar de T)."""

RNA_COMPLEMENT_IUPAC: Dict[str, str] = {
    "A": "U", "U": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
"""Complementaridade de RNA com IUPAC."""

PROTEIN_MOTIF_SYMBOLS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWYX")
"""Aminoacidos padrao mais X (qualquer residuo)."""


def _reverse_complement_rna(seq: str) -> str:
    """Complemento reverso de RNA com IUPAC (uso interno)."""
    out: List[str] = []
    for base in reversed(seq):
        if base not in RNA_COMPLEMENT_IUPAC:
            raise ValueError(f"Simbolo sem complemento definido: {base}.")
        out.append(RNA_COMPLEMENT_IUPAC[base])
    return "".join(out)


def _pattern_to_regex_rna(pattern: str) -> str:
    """Converte um padrao IUPAC de RNA em regex (uso interno)."""
    cleaned = "".join(pattern.split()).upper().replace("T", "U")
    if not cleaned:
        raise ValueError("O padrao de motivo esta vazio.")
    if len(cleaned) > MAX_MOTIF_PATTERN_LENGTH:
        raise ValueError(
            f"O padrao de motivo excede {MAX_MOTIF_PATTERN_LENGTH} simbolos."
        )
    parts: List[str] = []
    for symbol in cleaned:
        parts.append(IUPAC_TO_REGEX_RNA.get(symbol, re.escape(symbol)))
    return "".join(parts)


def find_motif_rna(seq: str, pattern: str) -> List[dict]:
    """Localiza um motivo IUPAC nas duas fitas de uma sequencia de RNA.

    Args:
        seq: Sequencia de RNA (A/C/G/U; T no padrao e tratado como U).
        pattern: Padrao IUPAC de RNA.

    Returns:
        Lista de hits com "start" (base 0), "end" (exclusivo), "strand" e
        "match", ordenada por start.

    Raises:
        ValueError: Se a sequencia ou o padrao forem vazios.

    Nota biologica:
        Motivos em RNA (por exemplo sitios de splicing ou seed de miRNA) podem
        ocorrer na fita transcrita ou, em dsRNA, na complementar.
    """
    cleaned = "".join(seq.split()).upper().replace("T", "U")
    if not cleaned:
        raise ValueError("A sequencia esta vazia.")
    regex = _pattern_to_regex_rna(pattern)
    compiled = re.compile(f"(?=({regex}))")
    results: List[dict] = []

    def _add(hit: dict) -> None:
        if len(results) >= MAX_MOTIF_HITS:
            raise ValueError(
                f"Motif search exceeded {MAX_MOTIF_HITS:,} hits. Narrow the "
                "pattern or search a shorter region. No partial hit list is returned."
            )
        results.append(hit)

    for match in compiled.finditer(cleaned):
        text = match.group(1)
        start = match.start()
        _add(
            {
                "start": start,
                "end": start + len(text),
                "strand": "+",
                "match": text,
            }
        )
    rc = _reverse_complement_rna(cleaned)
    n = len(cleaned)
    for match in compiled.finditer(rc):
        text = match.group(1)
        rc_start = match.start()
        forward_start = n - (rc_start + len(text))
        _add(
            {
                "start": forward_start,
                "end": forward_start + len(text),
                "strand": "-",
                "match": text,
            }
        )
    results.sort(key=lambda item: item["start"])
    return results


def find_motif_protein(seq: str, pattern: str) -> List[dict]:
    """Localiza um motivo proteico (aminoacidos padrao e X = qualquer).

    Nao calcula complemento reverso: proteinas nao tem fita complementar.

    Args:
        seq: Sequencia proteica.
        pattern: Padrao de aminoacidos; X casa qualquer um dos 20 padrao.

    Returns:
        Lista de hits com "start" (base 0), "end" (exclusivo), "strand" "+" e
        "match".

    Raises:
        ValueError: Se a sequencia ou o padrao forem vazios, ou se o padrao
            contiver simbolos fora dos 20 aminoacidos padrao e X.

    Nota biologica:
        Isto e busca exata/degenerada de sequencia, nao um PWM nem um dominio
        Pfam. Um hit nao significa dominio anotado.
    """
    cleaned = "".join(seq.split()).upper()
    if not cleaned:
        raise ValueError("A sequencia esta vazia.")
    motif = "".join(pattern.split()).upper()
    if not motif:
        raise ValueError("O padrao de motivo esta vazio.")
    if len(motif) > MAX_MOTIF_PATTERN_LENGTH:
        raise ValueError(
            f"O padrao de motivo excede {MAX_MOTIF_PATTERN_LENGTH} simbolos."
        )
    invalid = sorted({symbol for symbol in motif if symbol not in PROTEIN_MOTIF_SYMBOLS})
    if invalid:
        raise ValueError(
            "Protein motif symbols must be standard amino acids or X. Invalid: "
            + ", ".join(invalid)
        )
    parts: List[str] = []
    for symbol in motif:
        if symbol == "X":
            parts.append("[ACDEFGHIKLMNPQRSTVWY]")
        else:
            parts.append(re.escape(symbol))
    compiled = re.compile(f"(?=({''.join(parts)}))")
    results: List[dict] = []
    for match in compiled.finditer(cleaned):
        if len(results) >= MAX_MOTIF_HITS:
            raise ValueError(
                f"Motif search exceeded {MAX_MOTIF_HITS:,} hits. Narrow the "
                "pattern or search a shorter region. No partial hit list is returned."
            )
        text = match.group(1)
        start = match.start()
        results.append(
            {
                "start": start,
                "end": start + len(text),
                "strand": "+",
                "match": text,
            }
        )
    return results


def parse_motif_patterns(text: str) -> List[str]:
    """Separa padroes por virgula, ponto-e-virgula ou quebra de linha.

    Args:
        text: Texto com um ou mais padroes.

    Returns:
        Lista de padroes nao vazios, na ordem informada, sem duplicar entradas
        identicas consecutivas.

    Raises:
        ValueError: Se nenhum padrao restar ou se houver mais de MAX_MOTIF_PATTERNS.
    """
    raw = text or ""
    parts: List[str] = []
    for chunk in raw.replace(";", ",").split(","):
        for line in chunk.splitlines():
            item = line.strip()
            if item:
                parts.append(item)
    if not parts:
        raise ValueError("Provide at least one motif pattern.")
    if len(parts) > MAX_MOTIF_PATTERNS:
        raise ValueError(
            f"At most {MAX_MOTIF_PATTERNS} motif patterns are searched at once."
        )
    return parts


def filter_motif_hits(
    hits: List[dict],
    *,
    strand: str = "",
    min_start: int | None = None,
    max_end: int | None = None,
) -> List[dict]:
    """Filtra hits ja calculados sem repetir a busca.

    Args:
        hits: Lista no formato de find_motif.
        strand: "+" ou "-" para restringir; vazio mantem ambas.
        min_start: Inicio minimo 0-based, inclusive.
        max_end: Fim maximo exclusivo.

    Returns:
        Nova lista filtrada.

    Raises:
        ValueError: Se strand nao for vazio, "+" ou "-".
    """
    wanted = str(strand or "").strip()
    if wanted not in {"", "+", "-"}:
        raise ValueError('strand must be "+", "-" or empty.')
    selected: List[dict] = []
    for hit in hits:
        if wanted and str(hit.get("strand") or "") != wanted:
            continue
        start = int(hit["start"])
        end = int(hit["end"])
        if min_start is not None and start < min_start:
            continue
        if max_end is not None and end > max_end:
            continue
        selected.append(hit)
    return selected


def find_motifs(seq: str, patterns: List[str], molecule: str) -> List[dict]:
    """Busca varios padroes, rotulando cada hit com o padrao de origem.

    Args:
        seq: Sequencia alvo.
        patterns: Padroes ja separados por parse_motif_patterns.
        molecule: "DNA", "RNA" ou "PROTEIN".

    Returns:
        Hits combinados, cada um com "pattern", ordenados por start e padrao.
        O teto MAX_MOTIF_HITS aplica-se ao total combinado.

    Raises:
        ValueError: Se molecule for invalido, se patterns estiver vazio, ou se
            o total de hits exceder o teto.
    """
    key = str(molecule or "").strip().upper()
    if key not in {"DNA", "RNA", "PROTEIN"}:
        raise ValueError("molecule must be DNA, RNA or PROTEIN.")
    if not patterns:
        raise ValueError("Provide at least one motif pattern.")
    combined: List[dict] = []
    searchers = {
        "DNA": find_motif,
        "RNA": find_motif_rna,
        "PROTEIN": find_motif_protein,
    }
    search = searchers[key]
    for pattern in patterns:
        for hit in search(seq, pattern):
            if len(combined) >= MAX_MOTIF_HITS:
                raise ValueError(
                    f"Motif search exceeded {MAX_MOTIF_HITS:,} hits. Narrow the "
                    "pattern or search a shorter region. No partial hit list is returned."
                )
            item = dict(hit)
            item["pattern"] = "".join(pattern.split()).upper()
            combined.append(item)
    combined.sort(key=lambda item: (int(item["start"]), str(item.get("pattern") or "")))
    return combined
