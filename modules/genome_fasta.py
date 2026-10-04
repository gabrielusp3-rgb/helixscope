"""Indice FASTA (FAI-like) e acesso aleatorio sem carregar o genoma em RAM.

Formato de indice: cinco colunas tab-separadas no estilo samtools faidx
(name, length, offset, linebases, linewidth). O offset e o byte do primeiro
nucleotido apos o cabecalho. Coordenadas internas: 0-based half-open.
Exibicao: 1-based inclusive, rotulada explicitamente.

Nenhuma funcao importa Streamlit. Nenhum contig e renomeado (NC_000001.11
nao vira chr1).
"""

from __future__ import annotations

import hashlib
import json
import os
from typing import Any, Dict, Iterable, Mapping, Optional

from . import dna_analysis, provenance

FAI_LINEBASES_DEFAULT: int = 60
MAX_CONTIG_NAME: int = 256
MAX_FETCH_NT: int = 1_000_000
"""Teto de uma extracao de subsequencia (viewer/hit), nao do genoma."""

MappingLike = Mapping[str, Any]


class GenomeFastaError(ValueError):
    """Falha de indice ou de acesso aleatorio.

    Attributes:
        category: INVALID_INPUT, PARSING_ERROR, RESOURCE_LIMIT, STALE, ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def display_interval_1based_inclusive(start_0based: int, end_0based: int) -> dict:
    """Converte um intervalo interno 0-based half-open para 1-based inclusive.

    Args:
        start_0based: Inicio interno (incluido).
        end_0based: Fim interno (excluido).

    Returns:
        Dict start_1based, end_1based, convention.

    Raises:
        GenomeFastaError: INVALID_INPUT se o intervalo for impossivel.
    """
    start = int(start_0based)
    end = int(end_0based)
    if start < 0 or end < start:
        raise GenomeFastaError(
            "Interval must be 0-based half-open with end >= start >= 0.",
            "INVALID_INPUT",
        )
    if start == end:
        return {
            "start_1based": None,
            "end_1based": None,
            "empty": True,
            "convention": "1-based inclusive (empty interval has no display bases)",
        }
    return {
        "start_1based": start + 1,
        "end_1based": end,
        "empty": False,
        "convention": "1-based inclusive genomic display; internal is 0-based half-open",
    }


def parse_fai_text(text: str) -> dict:
    """Interpreta um ficheiro .fai.

    Args:
        text: Conteudo do indice.

    Returns:
        Dict contigs (name -> row), n_contigs.

    Raises:
        GenomeFastaError: PARSING_ERROR.
    """
    rows: Dict[str, dict] = {}
    for line_no, raw in enumerate(str(text or "").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 5:
            raise GenomeFastaError(
                f"FAI line {line_no} must have 5 tab-separated fields.",
                "PARSING_ERROR",
            )
        name = parts[0]
        if not name or len(name) > MAX_CONTIG_NAME:
            raise GenomeFastaError(
                f"FAI line {line_no} has an invalid contig name.",
                "PARSING_ERROR",
            )
        try:
            length = int(parts[1])
            offset = int(parts[2])
            linebases = int(parts[3])
            linewidth = int(parts[4])
        except ValueError as exc:
            raise GenomeFastaError(
                f"FAI line {line_no} has a non-integer numeric field.",
                "PARSING_ERROR",
            ) from exc
        if length < 0 or offset < 0 or linebases < 1 or linewidth < linebases:
            raise GenomeFastaError(
                f"FAI line {line_no} has impossible index numbers.",
                "PARSING_ERROR",
            )
        if name in rows:
            raise GenomeFastaError(
                f"FAI repeats contig name {name}. Names are not aliased.",
                "PARSING_ERROR",
            )
        rows[name] = {
            "name": name,
            "length": length,
            "offset": offset,
            "linebases": linebases,
            "linewidth": linewidth,
        }
    return {"contigs": rows, "n_contigs": len(rows)}


def load_fai(path: str) -> dict:
    """Le um .fai do disco.

    Args:
        path: Caminho do indice.

    Returns:
        parse_fai_text resultado mais path.

    Raises:
        GenomeFastaError: INVALID_INPUT se o ficheiro nao existir.
    """
    if not os.path.isfile(path):
        raise GenomeFastaError("FAI index file does not exist.", "INVALID_INPUT")
    with open(path, "r", encoding="utf-8") as handle:
        parsed = parse_fai_text(handle.read())
    parsed["path"] = os.path.abspath(path)
    return parsed


def write_fai(path: str, rows: Iterable[MappingLike]) -> str:
    """Grava um .fai.

    Args:
        path: Destino.
        rows: Linhas com name, length, offset, linebases, linewidth.

    Returns:
        Caminho absoluto.

    Raises:
        GenomeFastaError: INVALID_INPUT.
    """
    lines = []
    for row in rows:
        name = str(row.get("name") or "")
        if not name:
            raise GenomeFastaError("FAI row is missing contig name.", "INVALID_INPUT")
        lines.append(
            "\t".join(
                [
                    name,
                    str(int(row["length"])),
                    str(int(row["offset"])),
                    str(int(row["linebases"])),
                    str(int(row["linewidth"])),
                ]
            )
        )
    target = os.path.abspath(path)
    with open(target, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + ("\n" if lines else ""))
    return target


def build_fai_for_fasta(
    fasta_path: str,
    *,
    fai_path: str = "",
    reference_sha256: str = "",
) -> dict:
    """Constroi um indice FAI em streaming. Nao carrega o FASTA inteiro.

    Args:
        fasta_path: FASTA existente.
        fai_path: Destino; default fasta_path + '.fai'.
        reference_sha256: Hash do ficheiro FASTA ja conhecido; se vazio,
            calcula SHA-256 em streaming na mesma passagem.

    Returns:
        Dict fai_path, n_contigs, file_sha256, rows, algorithm.

    Raises:
        GenomeFastaError: PARSING_ERROR ou INVALID_INPUT.

    Nota:
        O formato FAI guarda uma unica largura de linha por contig, por isso o
        acesso aleatorio calcula o deslocamento de uma base a partir dessa
        largura. Um FASTA com linhas de comprimento irregular no meio de um
        contig e recusado: aceitar geraria um indice que devolve a sequencia
        errada em silencio, tal como o `samtools faidx` recusa o mesmo ficheiro.
    """
    source = os.path.abspath(fasta_path)
    if not os.path.isfile(source):
        raise GenomeFastaError("FASTA path does not exist.", "INVALID_INPUT")
    dest = os.path.abspath(fai_path or source + ".fai")
    hasher = hashlib.sha256()
    rows: list[dict] = []
    current: Optional[dict] = None
    first_seq_line_bases: Optional[int] = None
    first_seq_line_bytes: Optional[int] = None
    short_line_seen = False
    line_no = 0
    seq_len = 0
    offset = 0
    with open(source, "rb") as handle:
        while True:
            raw = handle.readline()
            if not raw:
                break
            hasher.update(raw)
            line_no += 1
            if raw.startswith(b">"):
                if current is not None:
                    current["length"] = seq_len
                    if first_seq_line_bases is None:
                        raise GenomeFastaError(
                            f"Contig {current['name']} has no sequence.",
                            "PARSING_ERROR",
                        )
                    current["linebases"] = first_seq_line_bases
                    current["linewidth"] = first_seq_line_bytes or (first_seq_line_bases + 1)
                    rows.append(current)
                header = raw.decode("utf-8", errors="replace").strip()[1:]
                name = header.split()[0] if header.split() else ""
                if not name:
                    raise GenomeFastaError("FASTA header has an empty contig name.", "PARSING_ERROR")
                current = {
                    "name": name,
                    "header": header,
                    "offset": handle.tell(),
                    "length": 0,
                    "linebases": 0,
                    "linewidth": 0,
                }
                first_seq_line_bases = None
                first_seq_line_bytes = None
                short_line_seen = False
                seq_len = 0
                continue
            if current is None:
                raise GenomeFastaError("FASTA sequence appeared before a header.", "PARSING_ERROR")
            bases = raw.strip().replace(b"\r", b"")
            n_bases = len(bases)
            if first_seq_line_bases is None:
                first_seq_line_bases = n_bases if n_bases else 1
                first_seq_line_bytes = len(raw)
            else:
                if short_line_seen:
                    raise GenomeFastaError(
                        f"FASTA line {line_no} continues contig {current['name']} after a "
                        f"line shorter than {first_seq_line_bases} bases. A FAI index "
                        "assumes a single line width per contig, so this file would "
                        "yield wrong offsets and return the wrong sequence.",
                        "PARSING_ERROR",
                    )
                if n_bases > first_seq_line_bases:
                    raise GenomeFastaError(
                        f"FASTA line {line_no} of contig {current['name']} has "
                        f"{n_bases} bases but the contig started with "
                        f"{first_seq_line_bases}. Line width must be constant for "
                        "random access to be correct.",
                        "PARSING_ERROR",
                    )
                if n_bases == first_seq_line_bases and len(raw) != first_seq_line_bytes:
                    raise GenomeFastaError(
                        f"FASTA line {line_no} of contig {current['name']} mixes line "
                        "terminators. A FAI index stores one line width per contig, so "
                        "mixed CRLF and LF would shift every downstream offset.",
                        "PARSING_ERROR",
                    )
                if n_bases < first_seq_line_bases:
                    short_line_seen = True
            seq_len += n_bases
            offset = handle.tell()
        if current is not None:
            if first_seq_line_bases is None:
                raise GenomeFastaError(
                    f"Contig {current['name']} has no sequence.",
                    "PARSING_ERROR",
                )
            current["length"] = seq_len
            current["linebases"] = first_seq_line_bases
            current["linewidth"] = first_seq_line_bytes or (first_seq_line_bases + 1)
            rows.append(current)
    if not rows:
        raise GenomeFastaError("FASTA contains no contigs.", "PARSING_ERROR")
    file_sha = hasher.hexdigest()
    expected = str(reference_sha256 or "").strip()
    if expected and expected != file_sha:
        raise GenomeFastaError(
            "FASTA SHA-256 does not match the hash recorded for this index.",
            "STALE",
        )
    write_fai(dest, rows)
    write_fai_provenance(
        dest,
        reference_sha256=file_sha,
        fasta_path=source,
        n_contigs=len(rows),
    )
    return {
        "fai_path": dest,
        "fasta_path": source,
        "n_contigs": len(rows),
        "file_sha256": file_sha,
        "rows": rows,
        "algorithm": "samtools-faidx-compatible 5-column index",
        "reference_hash": file_sha,
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def fai_provenance_path(fai_path: str) -> str:
    """Caminho do sidecar JSON do indice (hash da referencia).

    Args:
        fai_path: Ficheiro .fai.

    Returns:
        Caminho .fai.meta.json.

    Raises:
        Nenhum.
    """
    return os.path.abspath(str(fai_path) + ".meta.json")


def write_fai_provenance(
    fai_path: str,
    *,
    reference_sha256: str,
    fasta_path: str,
    n_contigs: int,
) -> dict:
    """Grava proveniencia do indice. Sem o hash, o FAI e tratado como stale.

    Args:
        fai_path: Indice.
        reference_sha256: SHA-256 do FASTA.
        fasta_path: FASTA indexado.
        n_contigs: Numero de contigs.

    Returns:
        Dict gravado.

    Raises:
        GenomeFastaError: INVALID_INPUT.
    """
    digest = str(reference_sha256 or "").strip()
    if not digest:
        raise GenomeFastaError("FAI provenance requires the FASTA SHA-256.", "INVALID_INPUT")
    source = os.path.abspath(fasta_path)
    payload = {
        "reference_sha256": digest,
        "fasta_path_basename": os.path.basename(source),
        "fasta_size": os.path.getsize(source),
        "fasta_mtime_ns": os.stat(source).st_mtime_ns,
        "n_contigs": int(n_contigs),
        "fai_path_basename": os.path.basename(os.path.abspath(fai_path)),
        "algorithm": "samtools-faidx-compatible 5-column index",
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
    dest = fai_provenance_path(fai_path)
    with open(dest, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return payload


def read_fai_provenance(fai_path: str) -> Optional[dict]:
    """Le o sidecar do indice, ou None.

    Args:
        fai_path: .fai.

    Returns:
        dict ou None.

    Raises:
        Nenhum.
    """
    path = fai_provenance_path(fai_path)
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def fai_is_current(fasta_path: str, fai_path: str, reference_sha256: str) -> bool:
    """True se o indice existe, o sidecar casa com o hash, e o FASTA nao mudou.

    Nao rehash o FASTA inteiro (genomas multi-GB). Mudanca de tamanho ou mtime
    invalida o indice. Hash vazio nunca e corrente.

    Args:
        fasta_path: FASTA.
        fai_path: Indice.
        reference_sha256: Hash esperado do FASTA (manifesto).

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    expected = str(reference_sha256 or "").strip()
    if not expected:
        return False
    if not os.path.isfile(fasta_path) or not os.path.isfile(fai_path):
        return False
    meta = read_fai_provenance(fai_path)
    if not meta:
        return False
    if str(meta.get("reference_sha256") or "") != expected:
        return False
    try:
        size = os.path.getsize(fasta_path)
        mtime_ns = os.stat(fasta_path).st_mtime_ns
    except OSError:
        return False
    if int(meta.get("fasta_size") or -1) != int(size):
        return False
    if int(meta.get("fasta_mtime_ns") or -1) != int(mtime_ns):
        return False
    return True


def file_sha256(path: str, *, chunk: int = 1024 * 1024) -> str:
    """SHA-256 de um ficheiro em blocos.

    Args:
        path: Ficheiro.
        chunk: Tamanho do bloco.

    Returns:
        Hex digest.

    Raises:
        GenomeFastaError: INVALID_INPUT.
    """
    if not os.path.isfile(path):
        raise GenomeFastaError("Cannot hash a missing file.", "INVALID_INPUT")
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(int(chunk))
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def file_md5(path: str, *, chunk: int = 1024 * 1024) -> str:
    """MD5 de um ficheiro (so para comparar com md5checksums.txt da NCBI).

    Args:
        path: Ficheiro.
        chunk: Tamanho do bloco.

    Returns:
        Hex digest.

    Raises:
        GenomeFastaError: INVALID_INPUT.
    """
    if not os.path.isfile(path):
        raise GenomeFastaError("Cannot hash a missing file.", "INVALID_INPUT")
    hasher = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(int(chunk))
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def fetch_sequence(
    fasta_path: str,
    fai_row: MappingLike,
    start_0based: int,
    end_0based: int,
    *,
    strand: str = "+",
) -> str:
    """Extrai [start, end) do contig via seek. Nao le o genoma inteiro.

    Args:
        fasta_path: FASTA indexado.
        fai_row: Linha FAI do contig.
        start_0based: Inicio 0-based incluido.
        end_0based: Fim 0-based excluido.
        strand: '+' sentido do ficheiro; '-' complemento reverso.

    Returns:
        Sequencia DNA maiuscula.

    Raises:
        GenomeFastaError: INVALID_INPUT ou RESOURCE_LIMIT.
    """
    start = int(start_0based)
    end = int(end_0based)
    length = int(fai_row["length"])
    if start < 0 or end < start or end > length:
        raise GenomeFastaError(
            "Requested interval is outside the contig (0-based half-open).",
            "INVALID_INPUT",
        )
    span = end - start
    if span > MAX_FETCH_NT:
        raise GenomeFastaError(
            f"Subsequence fetch exceeds {MAX_FETCH_NT} nt.",
            "RESOURCE_LIMIT",
        )
    if span == 0:
        return ""
    linebases = int(fai_row["linebases"])
    linewidth = int(fai_row["linewidth"])
    offset = int(fai_row["offset"])
    line_index = start // linebases
    column = start % linebases
    byte_pos = offset + line_index * linewidth + column
    collected: list[str] = []
    remaining = span
    with open(fasta_path, "rb") as handle:
        handle.seek(byte_pos)
        while remaining > 0:
            raw = handle.read(min(65536, remaining + remaining // linebases + 8))
            if not raw:
                break
            for byte in raw:
                if byte in (10, 13):
                    continue
                collected.append(chr(byte).upper())
                remaining -= 1
                if remaining <= 0:
                    break
    if len(collected) != span:
        raise GenomeFastaError(
            "FASTA fetch ended before the requested interval. Index may be stale.",
            "STALE",
        )
    sequence = "".join(collected)
    if strand == "-":
        sequence = dna_analysis.reverse_complement(sequence)
    elif strand not in {"+", ""}:
        raise GenomeFastaError("strand must be '+' or '-'.", "INVALID_INPUT")
    return sequence


def region_hash(
    fasta_path: str,
    fai_row: MappingLike,
    start_0based: int,
    end_0based: int,
    *,
    strand: str = "+",
) -> str:
    """SHA-256 da subsequencia extraida.

    Args:
        fasta_path: FASTA.
        fai_row: Linha FAI.
        start_0based: Inicio.
        end_0based: Fim.
        strand: Fita.

    Returns:
        Hex digest.

    Raises:
        GenomeFastaError: repassado de fetch_sequence.
    """
    sequence = fetch_sequence(
        fasta_path, fai_row, start_0based, end_0based, strand=strand
    )
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()
