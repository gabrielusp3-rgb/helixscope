"""Referencia genomica fornecida pelo usuario para busca CRISPR.

Nao baixa assemblies humanos/murinos automaticamente. Nao interpreta o nome do
arquivo (hg38.fa, genome.fa) como identidade de um genoma publicado. A
identidade cientifica e o conteudo FASTA (contigs, alfabeto, hash SHA-256) mais
os metadados que o usuario ou um registro NCBI realmente forneceram.

Nenhuma funcao aqui importa Streamlit. Nenhuma funcao aceita URL arbitraria nem
caminho de ficheiro como mecanismo de execucao.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from typing import Any, Mapping, Optional

from . import dna_analysis, provenance

MAX_REFERENCE_RAW_CHARS: int = 250_000
"""Teto de caracteres brutos do FASTA de referencia (limite tecnico)."""

MAX_REFERENCE_NT: int = 200_000
"""Teto da soma dos comprimentos dos contigs, em nucleotidos."""

MAX_REFERENCE_CONTIGS: int = 50
"""Teto de registros FASTA na referencia fornecida."""

MAX_UNCOMPRESSED_BYTES: int = 1_048_576
"""Teto de expansao gzip (protecao contra bomba de descompressao)."""

MAX_COMPRESSED_BYTES: int = 1_048_576
"""Teto do payload compactado ou em bruto antes de decodificar."""

ALLOWED_REFERENCE_ALPHABET: frozenset[str] = dna_analysis.DNA_ALPHABET
"""DNA canonico mais IUPAC. U, aminoacidos e simbolos livres sao recusados."""

SEARCH_SCOPE_PASTED: str = "pasted_sequence"
SEARCH_SCOPE_PROVIDED: str = "provided_reference"
SEARCH_SCOPE_NCBI_RECORD: str = "ncbi_nucleotide_record"
SEARCH_SCOPE_WORKSPACE: str = "workspace_dna"
"""Rotulos de escopo. Nenhum deles e genome-wide."""


class ReferenceError(ValueError):
    """Falha ao carregar ou validar uma referencia CRISPR.

    Attributes:
        category: INVALID_INPUT, PARSING_ERROR, RESOURCE_LIMIT ou ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def decode_reference_payload(data: bytes) -> str:
    """Decodifica bytes de upload para texto FASTA, com protecao gzip.

    Args:
        data: Conteudo do ficheiro enviado (possivelmente gzip).

    Returns:
        Texto UTF-8 do FASTA ou da sequencia crua.

    Raises:
        ReferenceError: RESOURCE_LIMIT se o payload ou a expansao excederem os
            tetos; INVALID_INPUT se o formato compactado nao for gzip de um
            unico membro; PARSING_ERROR se o texto nao for UTF-8.

    Nota biologica:
        Um ficheiro chamado hg38.fa.gz nao e o assembly hg38. Apenas o texto
        descomprimido entra na validacao posterior.
    """
    if not isinstance(data, (bytes, bytearray)):
        raise ReferenceError(
            "Reference payload must be bytes.",
            "INVALID_INPUT",
        )
    payload = bytes(data)
    if len(payload) > MAX_COMPRESSED_BYTES:
        raise ReferenceError(
            f"Reference file exceeds {MAX_COMPRESSED_BYTES:,} bytes. "
            "Paste a gene, plasmid or contig FASTA, not a chromosome-scale assembly.",
            "RESOURCE_LIMIT",
        )
    if payload[:2] == b"PK":
        raise ReferenceError(
            "ZIP archives are not accepted as CRISPR references. "
            "Upload a FASTA or a single-member .gz FASTA.",
            "INVALID_INPUT",
        )
    if payload[:2] == b"\x1f\x8b":
        payload = _decompress_gzip_limited(payload)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ReferenceError(
            "Reference FASTA is not valid UTF-8. No sequence was invented.",
            "PARSING_ERROR",
        ) from exc
    if len(text) > MAX_REFERENCE_RAW_CHARS:
        raise ReferenceError(
            f"Reference FASTA exceeds {MAX_REFERENCE_RAW_CHARS:,} characters. "
            "Paste a smaller contig or locus.",
            "RESOURCE_LIMIT",
        )
    return text


def _decompress_gzip_limited(payload: bytes) -> bytes:
    """Descomprime gzip recusando expansao acima de MAX_UNCOMPRESSED_BYTES."""
    try:
        handle = gzip.GzipFile(fileobj=io.BytesIO(payload), mode="rb")
    except OSError as exc:
        raise ReferenceError(
            "gzip header is present but the archive could not be opened.",
            "PARSING_ERROR",
        ) from exc
    chunks: list[bytes] = []
    total = 0
    try:
        while True:
            chunk = handle.read(65_536)
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_UNCOMPRESSED_BYTES:
                raise ReferenceError(
                    f"gzip expansion exceeded {MAX_UNCOMPRESSED_BYTES:,} bytes. "
                    "Treated as a decompression bomb; no sequence was loaded.",
                    "RESOURCE_LIMIT",
                )
            chunks.append(chunk)
    except ReferenceError:
        raise
    except OSError as exc:
        raise ReferenceError(
            "gzip FASTA could not be decompressed (corrupt or truncated).",
            "PARSING_ERROR",
        ) from exc
    finally:
        handle.close()
    return b"".join(chunks)


def load_reference_from_text(
    text: str,
    *,
    source: str,
    organism: str = "",
    assembly: str = "",
    version: str = "",
    accession: str = "",
    upload_filename: str = "",
) -> dict:
    """Interpreta FASTA ou DNA cru como referencia CRISPR fornecida.

    Args:
        text: FASTA (um ou varios registros) ou uma sequencia DNA sem cabecalho.
        source: Origem declarada (user FASTA, workspace DNA, NCBI record).
        organism: Organismo declarado pelo usuario ou pelo registro NCBI.
        assembly: Assembly declarado pelo usuario; nao e inferido do filename.
        version: Versao/release declarada, se houver.
        accession: Accession NCBI quando a sequencia veio de Entrez.
        upload_filename: Nome do ficheiro enviado, so para proveniencia; nunca
            usado como assembly ou organismo.

    Returns:
        Dicionario com contigs validados, hashes, metadados declarados, escopo
        e avisos. identity_hash e o SHA-256 do FASTA canonico.

    Raises:
        ReferenceError: INVALID_INPUT, PARSING_ERROR ou RESOURCE_LIMIT.

    Nota biologica:
        Off-targets dependem da referencia. Um plasmideo rotulado pelo usuario
        como hg38 continua a ser o plasmideo: o hash do conteudo e a identidade
        cientifica; o rotulo de assembly e declaracao do usuario.
    """
    if not isinstance(text, str):
        raise ReferenceError("Reference text must be a string.", "INVALID_INPUT")
    if len(text) > MAX_REFERENCE_RAW_CHARS:
        raise ReferenceError(
            f"Reference FASTA exceeds {MAX_REFERENCE_RAW_CHARS:,} characters.",
            "RESOURCE_LIMIT",
        )
    stripped = text.strip()
    if not stripped:
        raise ReferenceError(
            "Reference FASTA is empty. No genome was assumed.",
            "INVALID_INPUT",
        )
    try:
        parsed = dna_analysis.parse_fasta_records(
            text, max_records=MAX_REFERENCE_CONTIGS
        )
    except ValueError as exc:
        message = str(exc)
        category = "RESOURCE_LIMIT" if "limit" in message.lower() else "PARSING_ERROR"
        if "technical limit" in message.lower() or "exceed" in message.lower():
            category = "RESOURCE_LIMIT"
        raise ReferenceError(message, category) from exc

    raw_records = list(parsed.get("records") or [])
    if not raw_records:
        raise ReferenceError(
            "FASTA header found but no sequence records could be parsed.",
            "PARSING_ERROR",
        )

    is_raw = str(parsed.get("format") or "") == "raw"
    contigs: list[dict] = []
    seen_ids: dict[str, int] = {}
    total_nt = 0
    for index, record in enumerate(raw_records, start=1):
        identifier = str(record.get("identifier") or "").strip()
        identifier_source = "fasta_header"
        if not identifier:
            if is_raw and len(raw_records) == 1:
                identifier = "unnamed_contig_1"
                identifier_source = "assigned"
            else:
                raise ReferenceError(
                    f"FASTA record {index} has an empty identifier. "
                    "Empty headers are not replaced silently.",
                    "PARSING_ERROR",
                )
        if identifier in seen_ids:
            raise ReferenceError(
                f"Duplicate contig identifier '{identifier}' in the reference "
                f"FASTA (records {seen_ids[identifier]} and {index}). "
                "Duplicates were not merged.",
                "PARSING_ERROR",
            )
        seen_ids[identifier] = index
        sequence = _normalize_contig_sequence(str(record.get("sequence") or ""))
        if not sequence:
            raise ReferenceError(
                f"Contig '{identifier}' has an empty sequence.",
                "PARSING_ERROR",
            )
        extra = set(sequence) - ALLOWED_REFERENCE_ALPHABET
        if extra:
            shown = "".join(sorted(extra)[:12])
            raise ReferenceError(
                f"Contig '{identifier}' contains non-DNA symbols ({shown}). "
                "RNA (U) and protein alphabets are not accepted as CRISPR DNA "
                "references.",
                "INVALID_INPUT",
            )
        total_nt += len(sequence)
        if total_nt > MAX_REFERENCE_NT:
            raise ReferenceError(
                f"Reference exceeds {MAX_REFERENCE_NT:,} nt across contigs. "
                "This is not a genome-wide human/mouse search. Provide a locus, "
                "plasmid or small contig set.",
                "RESOURCE_LIMIT",
            )
        contigs.append(
            {
                "identifier": identifier,
                "identifier_source": identifier_source,
                "sequence": sequence,
                "length": len(sequence),
                "sequence_hash": provenance.sequence_digest(sequence),
            }
        )

    canonical = _canonical_fasta(contigs)
    identity_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    declared_organism = str(organism or "").strip()
    declared_assembly = str(assembly or "").strip()
    declared_version = str(version or "").strip()
    filename = str(upload_filename or "").strip()
    warnings: list[str] = []
    if filename and not declared_assembly:
        warnings.append(
            "The upload filename is recorded for provenance only and is not "
            "treated as an assembly name. Declare organism/assembly/version "
            "explicitly if they are known."
        )
    if declared_assembly and filename:
        stem = filename.lower().replace(".gz", "").replace(".fasta", "").replace(".fa", "")
        if stem and stem in declared_assembly.lower().replace(" ", ""):
            warnings.append(
                "Assembly was typed by the user. HelixScope did not verify that "
                "the FASTA content matches that public assembly."
            )
    if not declared_assembly:
        warnings.append(
            "No assembly/version was declared. The scientific identity of this "
            "reference is the FASTA content hash, not a public genome name."
        )

    scope = SEARCH_SCOPE_PROVIDED
    source_key = str(source or "user FASTA").strip() or "user FASTA"
    if source_key == "NCBI Entrez" or accession:
        scope = SEARCH_SCOPE_NCBI_RECORD
    elif source_key.lower().startswith("workspace"):
        scope = SEARCH_SCOPE_WORKSPACE

    packed = {
        "status": "READY",
        "source": source_key,
        "scope": scope,
        "scope_label": _scope_label(scope),
        "genome_wide": False,
        "organism_declared": declared_organism,
        "assembly_declared": declared_assembly,
        "version_declared": declared_version,
        "accession": str(accession or "").strip(),
        "upload_filename": filename,
        "filename_trusted_as_assembly": False,
        "verified_assembly": "",
        "assembly_identity_status": "user_declared_unverified",
        "contigs": contigs,
        "n_contigs": len(contigs),
        "total_nt": total_nt,
        "identity_hash": identity_hash,
        "software_version": provenance.HELIXSCOPE_VERSION,
        "timestamp": provenance.utc_now(),
        "warnings": warnings,
        "limits": {
            "max_reference_nt": MAX_REFERENCE_NT,
            "max_contigs": MAX_REFERENCE_CONTIGS,
        },
    }
    from . import crispr_assemblies

    return crispr_assemblies.annotate_reference(packed)


def load_reference_from_bytes(
    data: bytes,
    *,
    source: str = "user FASTA upload",
    organism: str = "",
    assembly: str = "",
    version: str = "",
    accession: str = "",
    upload_filename: str = "",
) -> dict:
    """Decodifica bytes de upload e valida o FASTA como referencia.

    Args:
        data: Bytes do ficheiro.
        source: Origem declarada.
        organism: Organismo declarado.
        assembly: Assembly declarado.
        version: Versao declarada.
        accession: Accession, se conhecido.
        upload_filename: Nome do ficheiro; nao e a identidade do genoma.

    Returns:
        O mesmo dicionario de load_reference_from_text.

    Raises:
        ReferenceError: Como decode_reference_payload e load_reference_from_text.
    """
    text = decode_reference_payload(data)
    return load_reference_from_text(
        text,
        source=source,
        organism=organism,
        assembly=assembly,
        version=version,
        accession=accession,
        upload_filename=upload_filename,
    )


def reference_from_ncbi_record(record: Mapping[str, Any]) -> dict:
    """Constroi uma referencia a partir de um registro NCBI ja recuperado.

    Args:
        record: Dicionario de ncbi_fetch.fetch_by_accession, ja na sessao.

    Returns:
        Referencia com um contig, accession, organismo e versao do registro.

    Raises:
        ReferenceError: INVALID_INPUT se nao houver sequencia nucleotidica
            recuperada; RESOURCE_LIMIT se exceder os tetos.

    Nota biologica:
        Isto nao e um assembly completo. E o registro Entrez que o usuario
        buscou. Nao e genome-wide.
    """
    if not isinstance(record, Mapping):
        raise ReferenceError(
            "NCBI record is missing. Fetch a nucleotide accession first.",
            "INVALID_INPUT",
        )
    if not bool(record.get("sequence_available")):
        raise ReferenceError(
            "This NCBI record has no retrieved nucleotide sequence. "
            "WGS masters and gene records cannot be used as CRISPR references.",
            "INVALID_INPUT",
        )
    database = str(record.get("database") or "").strip().lower()
    molecule = str(record.get("molecule") or "").strip().lower()
    if database == "protein" or molecule in {"aa", "protein"}:
        raise ReferenceError(
            "Protein NCBI records are not DNA CRISPR references.",
            "INVALID_INPUT",
        )
    sequence = str(record.get("sequence") or "")
    if not sequence:
        raise ReferenceError(
            "NCBI record sequence is empty. No bases were invented.",
            "INVALID_INPUT",
        )
    accession = str(record.get("accession") or record.get("version") or "NCBI")
    fasta = f">{accession}\n{sequence}\n"
    return load_reference_from_text(
        fasta,
        source="NCBI Entrez",
        organism=str(record.get("organism") or ""),
        assembly="",
        version=str(record.get("version") or accession),
        accession=accession,
        upload_filename="",
    )


def reference_from_workspace_dna(
    sequence: str,
    *,
    identifier: str = "",
    organism: str = "",
    accession: str = "",
    version: str = "",
) -> dict:
    """Usa a DNA do workspace como referencia fornecida, nao como genoma publico.

    Args:
        sequence: DNA ja normalizado do workspace.
        identifier: Identificador FASTA ou accession, se houver.
        organism: Organismo do snapshot, se houver.
        accession: Accession NCBI copiado para o workspace.
        version: Versao NCBI, se houver.

    Returns:
        Referencia de um contig.

    Raises:
        ReferenceError: Se a sequencia for vazia ou invalida.
    """
    seq = str(sequence or "").strip()
    if not seq:
        raise ReferenceError(
            "Workspace DNA is empty. Copy a nucleotide sequence first.",
            "INVALID_INPUT",
        )
    name = str(identifier or accession or "workspace_dna").strip() or "workspace_dna"
    fasta = f">{name}\n{seq}\n"
    return load_reference_from_text(
        fasta,
        source="workspace DNA",
        organism=organism,
        assembly="",
        version=version,
        accession=accession,
        upload_filename="",
    )


def contig_by_id(reference: Mapping[str, Any], contig_id: str) -> dict:
    """Devolve o contig cujo identifier casa exatamente.

    Args:
        reference: Saida de load_reference_from_text.
        contig_id: Identificador do contig.

    Returns:
        O dicionario do contig.

    Raises:
        ReferenceError: Se o contig nao existir. Nao escolhe o primeiro.
    """
    wanted = str(contig_id or "")
    for contig in list(reference.get("contigs") or []):
        if str(contig.get("identifier") or "") == wanted:
            return dict(contig)
    raise ReferenceError(
        f"Contig '{wanted}' is not in the loaded reference.",
        "INVALID_INPUT",
    )


def extract_sense_window(
    contig_sequence: str,
    start_0based: int,
    end_0based: int,
) -> str:
    """Extrai [start, end) da fita sense do contig.

    Args:
        contig_sequence: Sequencia do contig ja normalizada.
        start_0based: Inicio inclusivo 0-based.
        end_0based: Fim exclusivo 0-based.

    Returns:
        Subsequencia sense.

    Raises:
        ReferenceError: Se o intervalo cair fora do contig.
    """
    seq = str(contig_sequence or "")
    if start_0based < 0 or end_0based > len(seq) or start_0based >= end_0based:
        raise ReferenceError(
            "Requested window falls outside the contig. No sequence was invented.",
            "INVALID_INPUT",
        )
    return seq[start_0based:end_0based]


def reference_cache_key(reference: Mapping[str, Any]) -> str:
    """Identidade de cache: hash do conteudo + rotulos declarados + versao.

    Args:
        reference: Referencia validada.

    Returns:
        SHA-256 hexadecimal. Assembly A e assembly B nao compartilham cache
        mesmo se o FASTA for identico, porque o rotulo faz parte da identidade
        do resultado.

    Raises:
        ReferenceError: Se identity_hash estiver ausente.
    """
    content = str(reference.get("identity_hash") or "").strip()
    if not content:
        raise ReferenceError(
            "Reference has no identity_hash; it cannot be cached.",
            "ERROR",
        )
    payload = "|".join(
        [
            content,
            str(reference.get("organism_declared") or ""),
            str(reference.get("assembly_declared") or ""),
            str(reference.get("version_declared") or ""),
            str(reference.get("accession") or ""),
            str(reference.get("source") or ""),
            str(reference.get("software_version") or provenance.HELIXSCOPE_VERSION),
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def provenance_record(reference: Mapping[str, Any]) -> dict:
    """Monta o bloco de proveniencia da referencia.

    Args:
        reference: Referencia validada.

    Returns:
        Dict JSON-safe com fonte, organismo, assembly, versao, hash e avisos.
    """
    return provenance.json_safe(
        {
            "source": reference.get("source"),
            "organism_declared": reference.get("organism_declared") or "",
            "assembly_declared": reference.get("assembly_declared") or "",
            "version_declared": reference.get("version_declared") or "",
            "accession": reference.get("accession") or "",
            "identity_hash": reference.get("identity_hash"),
            "n_contigs": reference.get("n_contigs"),
            "total_nt": reference.get("total_nt"),
            "genome_wide": False,
            "scope": reference.get("scope"),
            "filename_trusted_as_assembly": False,
            "upload_filename": reference.get("upload_filename") or "",
            "software_version": reference.get("software_version"),
            "timestamp": reference.get("timestamp"),
        }
    )


def _normalize_contig_sequence(sequence: str) -> str:
    """Remove espacos e padroniza maiusculas; U nao e convertido em T."""
    return "".join(str(sequence or "").split()).upper()


def to_fasta_text(reference: Mapping[str, Any]) -> str:
    """Serializa os contigs validados para FASTA canonico.

    Args:
        reference: Referencia carregada por load_reference_from_text.

    Returns:
        Texto FASTA na ordem original. Nao e um genoma publico.

    Raises:
        ReferenceError: INVALID_INPUT se nao houver contigs.
    """
    contigs = list(reference.get("contigs") or [])
    if not contigs:
        raise ReferenceError("Reference has no contigs to serialize.", "INVALID_INPUT")
    return _canonical_fasta(contigs)


def _canonical_fasta(contigs: list[dict]) -> str:
    """FASTA canonico na ordem original dos registros, um contig por bloco."""
    parts: list[str] = []
    for contig in contigs:
        identifier = str(contig.get("identifier") or "")
        sequence = str(contig.get("sequence") or "")
        parts.append(f">{identifier}\n{sequence}\n")
    return "".join(parts)


def _scope_label(scope: str) -> str:
    """Texto de interface para o escopo; nunca diz genome-wide."""
    labels = {
        SEARCH_SCOPE_PASTED: (
            "Pasted target sequence only. Not genome-wide."
        ),
        SEARCH_SCOPE_PROVIDED: (
            "User-provided FASTA contigs. Not a public assembly and not genome-wide."
        ),
        SEARCH_SCOPE_NCBI_RECORD: (
            "Single NCBI nucleotide record already in this session. Not genome-wide."
        ),
        SEARCH_SCOPE_WORKSPACE: (
            "Workspace DNA currently loaded. Not genome-wide."
        ),
    }
    return labels.get(scope, "Provided sequences. Not genome-wide.")
