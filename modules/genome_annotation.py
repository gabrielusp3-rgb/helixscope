"""Anotacao GFF3/GTF tratada como input nao confiavel.

So liga anotacao a uma assembly quando o accession/id coincide. Checksum
oficial do .gff.gz e exigido para READY de anotacao publica. Overlaps usam
coordenadas reais. Hit em exao nao e efeito funcional.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import gzip
import json
import os
from typing import Any, Mapping, Optional, Sequence

from . import crispr_assemblies, genome_coordinates, genome_fasta, genome_store, provenance

MAX_GFF_BYTES: int = 512 * 1024 * 1024
MAX_GFF_FEATURES: int = 2_000_000
MAX_GFF_LINE: int = 32_768
MAX_ATTR_CHARS: int = 4_096
ANNOTATION_NAME: str = "annotation.gff3"
MANIFEST_NAME: str = "annotation_manifest.json"

CODING_TYPES: frozenset[str] = frozenset({"cds"})
UTR_TYPES: frozenset[str] = frozenset(
    {"five_prime_utr", "three_prime_utr", "utr", "3utr", "5utr"}
)
EXON_TYPES: frozenset[str] = frozenset({"exon"})
GENE_TYPES: frozenset[str] = frozenset({"gene", "ncrna_gene", "pseudogene"})


class AnnotationError(ValueError):
    """Falha de anotacao.

    Attributes:
        category: INVALID_INPUT, PARSING_ERROR, RESOURCE_LIMIT, CORRUPT.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def _safe_seqid(seqid: str) -> str:
    ident = str(seqid or "").strip()
    if not ident or len(ident) > 256:
        raise AnnotationError("GFF seqid is missing or too long.", "INVALID_INPUT")
    if any(ch in ident for ch in "/\\:*?\"<>|\t") or ident in {".", ".."}:
        raise AnnotationError(
            "GFF seqid contains path or delimiter characters. Refused.",
            "INVALID_INPUT",
        )
    if any(ord(ch) < 0x20 or ord(ch) == 0x7F for ch in ident):
        raise AnnotationError(
            "GFF seqid contains a control character. The GFF3 specification "
            "requires such bytes to be percent-encoded, and a NUL byte can "
            "truncate a path in a downstream C library.",
            "INVALID_INPUT",
        )
    return ident


def _parse_attributes(raw: str) -> dict:
    text = str(raw or "")[:MAX_ATTR_CHARS]
    attrs: dict[str, str] = {}
    if not text or text == ".":
        return attrs
    parts = text.split(";") if ";" in text else text.split()
    if "=" in text:
        for part in text.split(";"):
            item = part.strip()
            if not item:
                continue
            if "=" not in item:
                continue
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')[:512]
        return attrs
    for part in parts:
        item = part.strip()
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key.strip()] = value.strip().strip('"')[:512]
    return attrs


def parse_gff_text(text: str, *, declared_assembly: str) -> dict:
    """Interpreta GFF3 ou GTF limitado. Coordenadas internas 0-based half-open.

    Args:
        text: Conteudo (nao caminho).
        declared_assembly: Assembly declarada pelo ficheiro ou pelo catalogo.

    Returns:
        Dict features, n_features, declared_assembly.

    Raises:
        AnnotationError: PARSING_ERROR / RESOURCE_LIMIT / INVALID_INPUT.
    """
    if not str(declared_assembly or "").strip():
        raise AnnotationError(
            "Annotation must declare an assembly. Empty assembly is refused.",
            "INVALID_INPUT",
        )
    payload = str(text or "")
    if len(payload.encode("utf-8")) > MAX_GFF_BYTES:
        raise AnnotationError(
            "GFF/GTF exceeds the byte cap. Treated as untrusted oversized input.",
            "RESOURCE_LIMIT",
        )
    features: list[dict] = []
    for line_no, raw in enumerate(payload.splitlines(), start=1):
        if len(raw) > MAX_GFF_LINE:
            raise AnnotationError(
                f"GFF line {line_no} exceeds {MAX_GFF_LINE} characters.",
                "RESOURCE_LIMIT",
            )
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            raise AnnotationError(
                f"GFF line {line_no} does not have 8+ tab-separated columns.",
                "PARSING_ERROR",
            )
        seqid = _safe_seqid(cols[0])
        source = cols[1][:80]
        ftype = cols[2].strip().lower()
        if not ftype:
            raise AnnotationError(f"GFF line {line_no} has an empty type.", "PARSING_ERROR")
        try:
            start_1 = int(cols[3])
            end_1 = int(cols[4])
        except ValueError as exc:
            raise AnnotationError(
                f"GFF line {line_no} start/end are not integers.",
                "PARSING_ERROR",
            ) from exc
        if start_1 < 1 or end_1 < start_1:
            raise AnnotationError(
                f"GFF line {line_no} has an invalid 1-based interval.",
                "PARSING_ERROR",
            )
        strand = cols[6].strip()
        if strand not in {"+", "-", "."}:
            raise AnnotationError(
                f"GFF line {line_no} strand {strand!r} is not +, - or .",
                "PARSING_ERROR",
            )
        phase = cols[7].strip()
        if phase not in {"", ".", "0", "1", "2"}:
            raise AnnotationError(
                f"GFF line {line_no} has phase {phase!r}. The GFF3 specification "
                "allows only 0, 1, 2 or '.'; any other value would place the "
                "reading frame on the wrong base and translate the wrong protein.",
                "PARSING_ERROR",
            )
        attrs = _parse_attributes(cols[8] if len(cols) > 8 else "")
        internal = genome_coordinates.display_to_internal(start_1, end_1)
        gene_id = str(
            attrs.get("gene_id")
            or attrs.get("geneID")
            or attrs.get("gene")
            or attrs.get("ID")
            or ""
        )[:120]
        gene_name = str(attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene_symbol") or "")[:80]
        transcript_id = str(
            attrs.get("transcript_id")
            or attrs.get("transcript_id")
            or attrs.get("Parent")
            or ""
        )[:120]
        features.append(
            {
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start_0based": internal["start_0based"],
                "end_0based": internal["end_0based"],
                "start_1based": start_1,
                "end_1based": end_1,
                "strand": strand,
                "phase": phase,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "transcript_id": transcript_id,
                "attributes": attrs,
                "line_no": line_no,
            }
        )
        if len(features) > MAX_GFF_FEATURES:
            raise AnnotationError(
                "GFF feature cap reached. Remaining lines were not parsed.",
                "RESOURCE_LIMIT",
            )
    return {
        "declared_assembly": str(declared_assembly).strip(),
        "n_features": len(features),
        "features": features,
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def parse_gff_path(path: str, *, declared_assembly: str) -> dict:
    """Le um GFF/GTF ou .gz no store/fixtures.

    Args:
        path: Ficheiro.
        declared_assembly: Assembly declarada.

    Returns:
        Saida de parse_gff_text.

    Raises:
        AnnotationError.
    """
    abs_path = os.path.abspath(path)
    if not genome_store.path_is_in_store(abs_path):
        raise AnnotationError(
            "Annotation file is outside the HelixScope reference store/fixtures.",
            "INVALID_INPUT",
        )
    if not os.path.isfile(abs_path):
        raise AnnotationError("Annotation file does not exist.", "INVALID_INPUT")
    size = os.path.getsize(abs_path)
    if size > MAX_GFF_BYTES:
        raise AnnotationError("Annotation file exceeds the byte cap.", "RESOURCE_LIMIT")
    if abs_path.lower().endswith(".gz"):
        try:
            with gzip.open(abs_path, "rt", encoding="utf-8", errors="replace") as handle:
                text = handle.read(MAX_GFF_BYTES + 1)
        except OSError as exc:
            raise AnnotationError(f"gzip annotation could not be read: {exc}", "PARSING_ERROR") from exc
    else:
        with open(abs_path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read(MAX_GFF_BYTES + 1)
    if len(text.encode("utf-8")) > MAX_GFF_BYTES:
        raise AnnotationError("Decompressed annotation exceeds the byte cap.", "RESOURCE_LIMIT")
    return parse_gff_text(text, declared_assembly=declared_assembly)


def bind_annotation_to_assembly(assembly_id: str, declared_assembly: str) -> None:
    """Recusa anotacao de outro assembly.

    Args:
        assembly_id: id READY da referencia.
        declared_assembly: texto da anotacao.

    Returns:
        None.

    Raises:
        AnnotationError: CORRUPT se nao coincidir.
    """
    if not crispr_assemblies.gtf_matches_assembly(assembly_id, declared_assembly):
        raise AnnotationError(
            f"Annotation declared {declared_assembly!r} does not match "
            f"reference {assembly_id}. Wrong-assembly tracks are refused.",
            "CORRUPT",
        )


def build_interval_index(features: Sequence[Mapping[str, Any]]) -> dict:
    """Indice por contig: listas ordenadas por start interno.

    Args:
        features: Features parseadas.

    Returns:
        Dict contig -> lista de features.

    Raises:
        Nenhum.
    """
    buckets: dict[str, list[dict]] = {}
    for feat in features:
        seqid = str(feat.get("seqid") or "")
        buckets.setdefault(seqid, []).append(dict(feat))
    for seqid, rows in buckets.items():
        rows.sort(key=lambda item: (int(item["start_0based"]), int(item["end_0based"])))
        buckets[seqid] = rows
    return {
        "by_contig": buckets,
        "n_contigs": len(buckets),
        "n_features": sum(len(v) for v in buckets.values()),
    }


def overlapping_features(
    index: Mapping[str, Any],
    contig: str,
    start_0based: int,
    end_0based: int,
) -> list[dict]:
    """Features que intersectam [start, end) no contig.

    Args:
        index: Saida de build_interval_index.
        contig: seqid.
        start_0based: inicio interno.
        end_0based: fim exclusivo.

    Returns:
        Lista de features (copia rasa).

    Raises:
        AnnotationError: INVALID_INPUT.
    """
    ident = _safe_seqid(contig)
    try:
        start = int(start_0based)
        end = int(end_0based)
    except (TypeError, ValueError) as exc:
        raise AnnotationError("Overlap query coordinates must be integers.") from exc
    if start < 0 or end < start:
        raise AnnotationError("Overlap query interval is invalid.", "INVALID_INPUT")
    rows = list((index.get("by_contig") or {}).get(ident) or [])
    hits: list[dict] = []
    for feat in rows:
        f_start = int(feat["start_0based"])
        f_end = int(feat["end_0based"])
        if f_end <= start:
            continue
        if f_start >= end:
            break
        hits.append(dict(feat))
    return hits


def _gene_exons(features: Sequence[Mapping[str, Any]], gene_id: str) -> list[dict]:
    rows = []
    for feat in features:
        if str(feat.get("gene_id") or "") != gene_id:
            continue
        if str(feat.get("type") or "") in EXON_TYPES | CODING_TYPES | UTR_TYPES:
            rows.append(feat)
    return rows


def classify_hit_context(
    index: Mapping[str, Any],
    *,
    contig: str,
    start_0based: int,
    end_0based: int,
) -> dict:
    """Contexto genomico do intervalo. Sem claim de efeito.

    Args:
        index: Indice.
        contig: seqid.
        start_0based: inicio do spacer/hit.
        end_0based: fim exclusivo.

    Returns:
        Dict primary_label, overlaps, genes, transcripts, note.

    Raises:
        AnnotationError.
    """
    overlaps = overlapping_features(index, contig, start_0based, end_0based)
    types = {str(item.get("type") or "") for item in overlaps}
    genes: list[dict] = []
    seen_genes: set[tuple[str, str]] = set()
    transcripts: list[str] = []
    for item in overlaps:
        gene_id = str(item.get("gene_id") or "")
        gene_name = str(item.get("gene_name") or "")
        key = (gene_id, gene_name)
        if key != ("", "") and key not in seen_genes:
            seen_genes.add(key)
            genes.append(
                {
                    "gene_id": gene_id,
                    "symbol": gene_name,
                    "strand": str(item.get("strand") or "."),
                    "feature_type": str(item.get("type") or ""),
                }
            )
        tx = str(item.get("transcript_id") or "")
        if tx and tx not in transcripts:
            transcripts.append(tx)

    intron = False
    contig_feats = list((index.get("by_contig") or {}).get(contig) or [])
    for gene in genes:
        gid = str(gene.get("gene_id") or "")
        if not gid:
            continue
        gene_spans = [
            item
            for item in contig_feats
            if str(item.get("gene_id") or "") == gid and str(item.get("type") or "") in GENE_TYPES
        ]
        exons = _gene_exons(contig_feats, gid)
        if not gene_spans or not exons:
            continue
        in_gene = any(
            int(g["start_0based"]) < end_0based and int(g["end_0based"]) > start_0based
            for g in gene_spans
        )
        in_exon = any(
            int(ex["start_0based"]) < end_0based and int(ex["end_0based"]) > start_0based
            for ex in exons
        )
        if in_gene and not in_exon:
            intron = True

    if types & CODING_TYPES:
        primary = "CDS"
    elif types & UTR_TYPES:
        primary = "UTR"
    elif types & EXON_TYPES:
        primary = "exon"
    elif intron:
        primary = "intron"
    elif types & GENE_TYPES:
        primary = "gene"
    elif overlaps:
        primary = "feature"
    else:
        primary = "intergenic"

    coding = primary in {"CDS", "exon"}
    return {
        "primary_label": primary,
        "coding": coding,
        "noncoding": not coding and primary != "intergenic",
        "intergenic": primary == "intergenic",
        "overlaps": overlaps,
        "genes": genes,
        "transcripts": transcripts,
        "transcript_selection": "all_overlapping",
        "effect": None,
        "effect_note": (
            "Overlap with CDS/exon/intron/gene is coordinate context only. "
            "It is not deleterious, pathogenic or functional evidence."
        ),
        "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
        "display": genome_coordinates.internal_to_display(start_0based, end_0based),
    }


def annotate_hits(
    hits: Sequence[Mapping[str, Any]],
    index: Optional[Mapping[str, Any]],
    *,
    assembly_id: str,
    annotation_sha256: str = "",
) -> list[dict]:
    """Anexa contexto. Nao altera scores nem o resultado cru da busca.

    Args:
        hits: Hits verificados.
        index: Indice ou None (UNKNOWN).
        assembly_id: assembly da busca.
        annotation_sha256: identidade da anotacao.

    Returns:
        Nova lista de hits.

    Raises:
        Nenhum.
    """
    rows: list[dict] = []
    for hit in hits:
        packed = dict(hit)
        if index is None:
            packed["genomic_context"] = {
                "primary_label": "UNKNOWN",
                "genes": [],
                "transcripts": [],
                "effect": None,
                "effect_note": "No matching GFF/GTF is loaded for this assembly.",
                "annotation_sha256": "",
            }
            packed["feature_annotation"] = "UNKNOWN"
            rows.append(packed)
            continue
        contig = str(hit.get("chromosome_or_contig") or "")
        try:
            context = classify_hit_context(
                index,
                contig=contig,
                start_0based=int(hit.get("start_0based") or 0),
                end_0based=int(hit.get("end_0based") or 0),
            )
        except AnnotationError as exc:
            packed["genomic_context"] = {
                "primary_label": "UNKNOWN",
                "reason": str(exc),
                "effect": None,
                "annotation_sha256": annotation_sha256,
            }
            packed["feature_annotation"] = "UNKNOWN"
            rows.append(packed)
            continue
        context["annotation_sha256"] = annotation_sha256
        context["assembly_id"] = assembly_id
        packed["genomic_context"] = context
        packed["feature_annotation"] = context["primary_label"]
        rows.append(packed)
    return rows


def context_cache_identity(*, search_identity: str, annotation_sha256: str) -> str:
    """Chave de cache de contexto. Separada do resultado cru.

    Args:
        search_identity: identidade da busca Cas-OFFinder.
        annotation_sha256: hash da anotacao.

    Returns:
        SHA-256 hex.

    Raises:
        Nenhum.
    """
    import hashlib

    payload = f"{search_identity}|ann|{annotation_sha256}|{provenance.HELIXSCOPE_VERSION}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def annotation_catalog_urls(assembly_id: str) -> dict:
    """URLs oficiais NCBI GFF/GTF para o mesmo accession.

    Args:
        assembly_id: id de catalogo.

    Returns:
        Dict urls e MD5 oficiais, ou vazio.

    Raises:
        Nenhum.
    """
    row = crispr_assemblies.catalog_by_id(assembly_id) or {}
    return {
        "gff_url": str(row.get("official_gff_url") or ""),
        "gff_md5": str(row.get("official_gff_md5") or ""),
        "gtf_url": str(row.get("official_gtf_url") or ""),
        "gtf_md5": str(row.get("official_gtf_md5") or ""),
        "accession": str(row.get("accession") or ""),
    }


def load_ready_annotation(assembly_id: str) -> Optional[dict]:
    """Indice READY no store, ou None.

    Args:
        assembly_id: id.

    Returns:
        Dict index, sha256, declared_assembly, ou None.

    Raises:
        Nenhum.
    """
    directory = genome_store.assembly_dir(assembly_id)
    manifest_path = os.path.join(directory, MANIFEST_NAME)
    gff_path = os.path.join(directory, ANNOTATION_NAME)
    if not os.path.isfile(manifest_path) or not os.path.isfile(gff_path):
        return None
    try:
        with open(manifest_path, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(manifest, dict):
        return None
    if str(manifest.get("status") or "") != "READY":
        return None
    expected = str(manifest.get("file_sha256") or "")
    try:
        observed = genome_fasta.file_sha256(gff_path)
    except genome_fasta.GenomeFastaError:
        return None
    if expected and observed.lower() != expected.lower():
        return None
    try:
        parsed = parse_gff_path(
            gff_path,
            declared_assembly=str(manifest.get("declared_assembly") or assembly_id),
        )
        bind_annotation_to_assembly(assembly_id, parsed["declared_assembly"])
    except AnnotationError:
        return None
    return {
        "index": build_interval_index(parsed["features"]),
        "file_sha256": observed,
        "declared_assembly": parsed["declared_assembly"],
        "n_features": parsed["n_features"],
        "path": gff_path,
        "manifest": manifest,
    }


def install_annotation_from_text(
    assembly_id: str,
    text: str,
    *,
    declared_assembly: str,
) -> dict:
    """Grava anotacao local apos match de assembly. Uso de testes e fixtures.

    Args:
        assembly_id: id da referencia.
        text: GFF3/GTF.
        declared_assembly: deve coincidir.

    Returns:
        Manifesto READY da anotacao.

    Raises:
        AnnotationError.
    """
    bind_annotation_to_assembly(assembly_id, declared_assembly)
    parsed = parse_gff_text(text, declared_assembly=declared_assembly)
    directory = genome_store.assembly_dir(assembly_id)
    os.makedirs(directory, exist_ok=True)
    dest = os.path.join(directory, ANNOTATION_NAME)
    with open(dest, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    digest = genome_fasta.file_sha256(dest)
    manifest = {
        "status": "READY",
        "assembly_id": assembly_id,
        "declared_assembly": declared_assembly,
        "file": dest,
        "file_sha256": digest,
        "n_features": parsed["n_features"],
        "source": "local_text",
        "verified_at": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
    path = os.path.join(directory, MANIFEST_NAME)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest
