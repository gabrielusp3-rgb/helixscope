"""Catalogo de assemblies publicos para rotulagem honesta de referencias CRISPR.

Nao contem sequencias. Nao descarrega genomas no import. Nao trata o nome de
um ficheiro (hg38.fa) como prova de identidade. verified_assembly so e
preenchido quando o accession casa exactamente com o catalogo; READY so apos
checksum do ficheiro local.

Checksums oficiais: NCBI md5checksums.txt recuperados em 2026-08-28 para o
ficheiro *_genomic.fna.gz de cada assembly. Nao sao placeholders.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional

from . import provenance

CATALOG_SNAPSHOT: str = "2026-08-28"
CATALOG_SOURCE: str = (
    "NCBI Assembly / RefSeq FTP md5checksums.txt retrieved 2026-08-28: "
    "GRCh38.p14 GCF_000001405.40; T2T-CHM13v2.0 GCF_009914755.1; "
    "GRCm39 GCF_000001635.27. Official MD5 applies to genomic.fna.gz only."
)

NCBI_FTP_HTTPS: str = "https://ftp.ncbi.nlm.nih.gov"

ASSEMBLY_CATALOG: tuple[dict, ...] = (
    {
        "id": "GRCh38.p14",
        "kind": "public_assembly",
        "organism": "Homo sapiens",
        "assembly": "GRCh38.p14",
        "assembly_level": "Chromosome",
        "release": "patch 14; Genome Reference Consortium 2022-02-03",
        "release_date": "2022-02-03",
        "source": "NCBI RefSeq",
        "accession": "GCF_000001405.40",
        "genbank_accession": "GCA_000001405.29",
        "nicknames": ("hg38", "grch38", "grch38.p14", "gcf_000001405.40"),
        "checksum": "",
        "official_compressed_md5": "c30471567037b2b2389d43c908c653e1",
        "official_checksum_algo": "MD5",
        "official_checksum_file": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/md5checksums.txt"
        ),
        "official_fasta_name": "GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        "download_url": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
        ),
        "estimated_compressed_bytes": 1_000_000_000,
        "estimated_size_label": "about 1 GB compressed genomic.fna.gz (NCBI FTP listing order)",
        "checksum_status": "official_ncbi_md5checksums_retrieved",
        "downloaded": False,
        "status": "catalog_metadata_only",
        "license_note": "NCBI RefSeq sequence data; cite the assembly accession.",
        "notes": (
            "Official human reference patch 14. HelixScope does not ship this "
            "assembly. Typing hg38 does not verify FASTA content. READY requires "
            "a completed download or local genomic.fna.gz whose MD5 matches."
        ),
        "reference_scope": (
            "NCBI RefSeq genomic.fna.gz for GCF_000001405.40 after decompress "
            "(all sequences in that file: chromosomes, unlocalized, unplaced, "
            "and any alternate loci NCBI included). HelixScope does not silently "
            "subset to primary assembly, patches-only, or decoys-only. Not "
            "GRCh37, not hg19, not UCSC hg38.fa."
        ),
        "official_gff_name": "GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
        "official_gff_md5": "24b731562b9d4cae9e37b23404b5be16",
        "official_gff_url": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
        ),
        "official_gtf_name": "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
        "official_gtf_md5": "a561524a1ac438a2a95aed54d00e490e",
        "official_gtf_url": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
        ),
        "annotation_checksum_file": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/405/"
            "GCF_000001405.40_GRCh38.p14/md5checksums.txt"
        ),
    },
    {
        "id": "T2T-CHM13v2.0",
        "kind": "public_assembly",
        "organism": "Homo sapiens",
        "assembly": "T2T-CHM13v2.0",
        "assembly_level": "Complete Genome",
        "release": "T2T Consortium 2022-01-24",
        "release_date": "2022-01-24",
        "source": "NCBI RefSeq",
        "accession": "GCF_009914755.1",
        "genbank_accession": "GCA_009914755.4",
        "nicknames": ("t2t", "chm13", "t2t-chm13v2.0", "hs1", "gcf_009914755.1"),
        "checksum": "",
        "official_compressed_md5": "9e6bf6b586bc8954208d1cc1d5f2fc99",
        "official_checksum_algo": "MD5",
        "official_checksum_file": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/009/914/755/"
            "GCF_009914755.1_T2T-CHM13v2.0/md5checksums.txt"
        ),
        "official_fasta_name": "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz",
        "download_url": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/009/914/755/"
            "GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
        ),
        "estimated_compressed_bytes": 889_000_000,
        "estimated_size_label": "about 889 MB compressed genomic.fna.gz (NCBI FTP listing)",
        "checksum_status": "official_ncbi_md5checksums_retrieved",
        "downloaded": False,
        "status": "catalog_metadata_only",
        "license_note": "NCBI RefSeq / T2T Consortium; cite GCF_009914755.1.",
        "notes": (
            "Telomere-to-telomere alternate human assembly. Not GRCh38. Not "
            "downloaded at startup."
        ),
    },
    {
        "id": "GRCm39",
        "kind": "public_assembly",
        "organism": "Mus musculus",
        "assembly": "GRCm39",
        "assembly_level": "Chromosome",
        "release": "Genome Reference Consortium 2020-06-24",
        "release_date": "2020-06-24",
        "source": "NCBI RefSeq",
        "accession": "GCF_000001635.27",
        "genbank_accession": "GCA_000001635.9",
        "nicknames": ("mm39", "grcm39", "gcf_000001635.27"),
        "checksum": "",
        "official_compressed_md5": "c0b0c4c3f54d2b480efe68a18bf7e42b",
        "official_checksum_algo": "MD5",
        "official_checksum_file": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/635/"
            "GCF_000001635.27_GRCm39/md5checksums.txt"
        ),
        "official_fasta_name": "GCF_000001635.27_GRCm39_genomic.fna.gz",
        "download_url": (
            f"{NCBI_FTP_HTTPS}/genomes/all/GCF/000/001/635/"
            "GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
        ),
        "estimated_compressed_bytes": 800_000_000,
        "estimated_size_label": "about 0.8 GB compressed genomic.fna.gz",
        "checksum_status": "official_ncbi_md5checksums_retrieved",
        "downloaded": False,
        "status": "catalog_metadata_only",
        "license_note": "NCBI RefSeq mouse reference; cite GCF_000001635.27.",
        "notes": (
            "Current major mouse reference (C57BL/6J). Not mm10/GRCm38. Not "
            "downloaded at startup."
        ),
    },
)
"""Metadados oficiais. Nenhuma sequencia esta embutida."""

TEST_REFERENCE_ID: str = "HELIXSCOPE_TEST_REF"


def list_catalog(*, include_test: bool = False) -> list[dict]:
    """Copia do catalogo com o snapshot de pesquisa.

    Args:
        include_test: Se True, inclui o descritor do FASTA de teste (nunca
            apresentado como genoma publico).

    Returns:
        Lista de dicts (sem sequencias).

    Raises:
        Nenhum.
    """
    rows = []
    for item in ASSEMBLY_CATALOG:
        row = dict(item)
        row["catalog_snapshot"] = CATALOG_SNAPSHOT
        row["catalog_source"] = CATALOG_SOURCE
        rows.append(row)
    if include_test:
        rows.append(test_reference_descriptor())
    return rows


def catalog_by_id(assembly_id: str) -> Optional[dict]:
    """Uma entrada do catalogo publico, ou None.

    Args:
        assembly_id: id interno (GRCh38.p14, ...).

    Returns:
        Copia do dict ou None.

    Raises:
        Nenhum.
    """
    key = str(assembly_id or "").strip()
    if key == TEST_REFERENCE_ID:
        return test_reference_descriptor()
    for item in ASSEMBLY_CATALOG:
        if str(item.get("id") or "") == key:
            row = dict(item)
            row["catalog_snapshot"] = CATALOG_SNAPSHOT
            row["catalog_source"] = CATALOG_SOURCE
            return row
    return None


def test_reference_descriptor() -> dict:
    """Descritor do FASTA sintetico de testes. Nao e um assembly NCBI.

    Args:
        Nenhum.

    Returns:
        Dict kind=synthetic_test_reference.

    Raises:
        Nenhum.
    """
    return {
        "id": TEST_REFERENCE_ID,
        "kind": "synthetic_test_reference",
        "organism": "HelixScope test construct (not an organism)",
        "assembly": TEST_REFERENCE_ID,
        "assembly_level": "test_contigs",
        "source": "HelixScope tests/fixtures/helixscope_test_reference.fa",
        "accession": "",
        "nicknames": (),
        "checksum": "",
        "official_compressed_md5": "",
        "checksum_status": "local_fixture_sha256_after_install",
        "downloaded": False,
        "status": "test_fixture",
        "estimated_compressed_bytes": 2048,
        "estimated_size_label": "< 2 KB",
        "notes": (
            "Synthetic multi-contig FASTA for unit/live Cas-OFFinder tests. "
            "It is not GRCh38, not T2T, not GRCm39, and must never be labelled "
            "genome-wide for a public assembly."
        ),
        "reference_scope": (
            "Synthetic two-contig FASTA (HS_TEST_1, HS_TEST_2) only. "
            "Completed search status is COMPLETED_TEST_REFERENCE, never "
            "COMPLETED_FULL_REFERENCE."
        ),
        "catalog_snapshot": CATALOG_SNAPSHOT,
        "not_a_public_assembly": True,
    }


def lookup_declaration(
    declared_assembly: str = "",
    *,
    declared_organism: str = "",
    accession: str = "",
) -> dict:
    """Compara a declaracao do usuario com o catalogo, sem verificar FASTA.

    Args:
        declared_assembly: Texto digitado (hg38, GRCh38.p14, ...).
        declared_organism: Organismo declarado.
        accession: Accession se o usuario ou NCBI o forneceu.

    Returns:
        Dict com declared_assembly, catalog_match, verified_assembly (so se o
        accession casar), status e avisos. Nunca afirma que o FASTA e esse
        assembly. checksum (conteudo local) permanece vazio.

    Raises:
        Nenhum.
    """
    declared = str(declared_assembly or "").strip()
    org = str(declared_organism or "").strip()
    acc = str(accession or "").strip()
    key = _normalize(declared)
    acc_key = _normalize(acc)
    match: Optional[dict] = None
    match_via = ""
    if acc_key:
        for item in ASSEMBLY_CATALOG:
            if acc_key == _normalize(str(item.get("accession") or "")):
                match = dict(item)
                match_via = "accession"
                break
            if acc_key == _normalize(str(item.get("genbank_accession") or "")):
                match = dict(item)
                match_via = "genbank_accession"
                break
    if match is None and key:
        for item in ASSEMBLY_CATALOG:
            nicknames = {_normalize(str(item.get("assembly") or ""))}
            nicknames.update(_normalize(n) for n in item.get("nicknames") or ())
            if key in nicknames:
                match = dict(item)
                match_via = "declared_name_or_nickname"
                break
    verified = ""
    identity_status = "user_declared_unverified"
    warnings: list[str] = []
    if match is None:
        if declared or acc:
            warnings.append(
                "Declaration did not match the HelixScope catalog snapshot. "
                "That does not prove the label is wrong; it was not verified "
                "against FASTA content."
            )
        else:
            identity_status = "no_assembly_declared"
    else:
        if match_via in {"accession", "genbank_accession"}:
            verified = str(match.get("assembly") or "")
            identity_status = "accession_matches_catalog_not_fasta"
            warnings.append(
                "Accession matches a catalog entry. FASTA bases were not "
                "checksum-verified against that public assembly unless a local "
                "READY reference exists. verified_assembly is the catalog name."
            )
        else:
            identity_status = "nickname_matches_catalog_not_verified"
            warnings.append(
                f"'{declared}' is a known nickname/name for catalog assembly "
                f"{match.get('assembly')} ({match.get('accession')}). The "
                "uploaded or pasted FASTA was not verified to be that genome."
            )
        if org and match.get("organism") and _normalize(org) not in _normalize(
            str(match.get("organism") or "")
        ):
            warnings.append(
                "Declared organism does not match the catalog organism for "
                "that assembly name. Neither label was promoted to fact."
            )
    return provenance.json_safe(
        {
            "declared_assembly": declared,
            "declared_organism": org,
            "declared_accession": acc,
            "catalog_match": match,
            "match_via": match_via,
            "verified_assembly": verified,
            "assembly_identity_status": identity_status,
            "filename_trusted_as_assembly": False,
            "downloaded": False,
            "checksum": "",
            "official_compressed_md5": (match or {}).get("official_compressed_md5") or "",
            "catalog_snapshot": CATALOG_SNAPSHOT,
            "warnings": warnings,
        }
    )


def annotate_reference(reference: Mapping[str, Any]) -> dict:
    """Acrescenta lookup de catalogo a uma referencia ja carregada.

    Args:
        reference: Saida de crispr_reference.load_reference_from_text.

    Returns:
        Copia com assembly_catalog e verified_assembly. genome_wide permanece
        False para FASTA fornecido na sessao.

    Raises:
        Nenhum.
    """
    packed = dict(reference)
    lookup = lookup_declaration(
        str(reference.get("assembly_declared") or ""),
        declared_organism=str(reference.get("organism_declared") or ""),
        accession=str(reference.get("accession") or ""),
    )
    packed["assembly_catalog"] = lookup
    packed["verified_assembly"] = lookup.get("verified_assembly") or ""
    packed["assembly_identity_status"] = lookup.get("assembly_identity_status")
    packed["filename_trusted_as_assembly"] = False
    packed["genome_wide"] = False
    warnings = list(packed.get("warnings") or [])
    warnings.extend(list(lookup.get("warnings") or []))
    packed["warnings"] = warnings
    return packed


def gtf_matches_assembly(assembly_id: str, annotation_declared_assembly: str) -> bool:
    """Recusa GTF/GFF rotulado para outro assembly (ex. GRCh37 em GRCh38).

    Args:
        assembly_id: id do catalogo da referencia READY.
        annotation_declared_assembly: Texto associado ao ficheiro de anotacao.

    Returns:
        True so se os identificadores normalizados coincidirem.

    Raises:
        Nenhum.
    """
    row = catalog_by_id(assembly_id)
    if row is None:
        return False
    declared = _normalize(annotation_declared_assembly)
    if row.get("kind") == "synthetic_test_reference":
        return declared in {
            _normalize(TEST_REFERENCE_ID),
            "helixscopetestref",
        }
    if not declared:
        return False
    tokens = {
        _normalize(str(row.get("id") or "")),
        _normalize(str(row.get("assembly") or "")),
        _normalize(str(row.get("accession") or "")),
        _normalize(str(row.get("genbank_accession") or "")),
    }
    tokens.update(_normalize(n) for n in row.get("nicknames") or ())
    if assembly_id == "GRCh38.p14" and declared in {"grch37", "hg19", "b37"}:
        return False
    return declared in tokens


def _normalize(text: str) -> str:
    return "".join(ch for ch in str(text or "").lower() if ch.isalnum())
