"""Recuperacao de registros de sequencias e genes no NCBI via Entrez (Biopython).

Recupera um registro por identificador usando as E-utilities: Entrez.efetch para
nucleotide e protein (formato GenBank/GenPept) e Entrez.esummary/elink para gene.
Tambem classifica o identificador para escolher o banco apropriado e interpreta
as localizacoes das features CDS. Aplica limitacao de taxa conforme as regras do
NCBI. Nenhuma funcao aqui importa Streamlit; o Biopython e importado de forma
tardia.
"""

from __future__ import annotations

import os
import re
import time
import urllib.error
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple
import socket

SLEEP_WITHOUT_KEY: float = 0.34
"""Atraso (s) entre requisicoes sem API key (limite de ~3 req/s do NCBI)."""

SLEEP_WITH_KEY: float = 0.11
"""Atraso (s) entre requisicoes com API key (limite de ~10 req/s do NCBI)."""

ENTREZ_TOOL: str = "HelixScope"
"""Valor de Entrez.tool enviado em toda requisicao, conforme a politica do NCBI."""

MAX_SEQUENCE_FETCH_NT: int = 200_000
"""Limite de bases ou residuos baixados por efetch. Scaffolds WGS cromossomicos
nao cabem na interface Streamlit e nao devem ser lidos por inteiro."""

MAX_IDENTIFIER_LENGTH: int = 64
"""Tamanho maximo de um accession ou GeneID aceito antes da consulta Entrez."""

MAX_GENE_ID_DIGITS: int = 12
"""Numero maximo de digitos em um GeneID; IDs reais do NCBI cabem folgadamente."""

MAX_EMAIL_LENGTH: int = 254
"""Comprimento maximo de um e-mail de contato Entrez (RFC 5321)."""

ENTREZ_TIMEOUT_S: float = 45.0
"""Timeout em segundos de cada requisicao HTTP ao NCBI."""

SUPPORTED_DATABASES: tuple = ("nucleotide", "protein", "gene", "pubmed")
"""Bancos Entrez consultados por fetch_by_accession (pubmed e PMID)."""

SEARCHABLE_DATABASES: tuple = ("nucleotide", "protein", "gene", "pubmed")
"""Bancos Entrez consultados por search_records (esearch + esummary)."""

MAX_SEARCH_TERM_LENGTH: int = 200
"""Tamanho maximo do termo de busca Entrez."""

MAX_SEARCH_RETMAX: int = 20
"""Numero maximo de resumos devolvidos por busca."""

SEARCH_RETRY_ATTEMPTS: int = 3
"""Tentativas de esearch/esummary diante de timeout, HTTP 429 ou 5xx."""

LOCATION_SEGMENT_PATTERN: re.Pattern[str] = re.compile(
    r"\[<?(\d+)\:>?(\d+)\](?:\(([+?-])\))?"
)
"""Padrao de um segmento na representacao textual de uma localizacao Biopython,
como "[0:1234](+)"; aceita limites incertos ("<" e ">") e fita ausente."""

PROTEIN_REFSEQ_PREFIXES: frozenset[str] = frozenset(
    {"NP_", "XP_", "YP_", "WP_", "AP_"}
)
"""Prefixos RefSeq de registros de proteina."""

NUCLEOTIDE_REFSEQ_PREFIXES: frozenset[str] = frozenset(
    {
        "NM_",
        "NR_",
        "XM_",
        "XR_",
        "NC_",
        "NG_",
        "NT_",
        "NW_",
        "NZ_",
        "AC_",
        "AE_",
        "CP_",
        "CH_",
    }
)
"""Prefixos RefSeq e de cromossomo/contig frequentemente depositados em nucleotide."""

ASSEMBLY_PREFIXES: tuple = ("GCA_", "GCF_")
"""Prefixos de accession do banco Assembly, que nao e consultado aqui."""

REFSEQ_PATTERN: re.Pattern[str] = re.compile(
    r"^([A-Z]{2})_\d{5,12}(?:\.\d+)?$", re.IGNORECASE
)
"""Accession RefSeq do tipo PREFIXO_digitos, com versao opcional."""

WGS_FOUR_PATTERN: re.Pattern[str] = re.compile(
    r"^[A-Z]{4}\d{8,10}(?:\.\d+)?$", re.IGNORECASE
)
"""Projeto WGS classico: quatro letras e oito a dez digitos (ex.: AGTM000000000)."""

WGS_SIX_PATTERN: re.Pattern[str] = re.compile(
    r"^[A-Z]{6}\d{8,11}(?:\.\d+)?$", re.IGNORECASE
)
"""Projeto WGS de seis letras (ex.: JBFSEQ000000000)."""

GENBANK_NUC_SHORT: re.Pattern[str] = re.compile(
    r"^[A-Z]\d{5}(?:\.\d+)?$", re.IGNORECASE
)
"""Accession GenBank classico de nucleotideo: uma letra e cinco digitos."""

GENBANK_NUC_LONG: re.Pattern[str] = re.compile(
    r"^[A-Z]{2}\d{6,8}(?:\.\d+)?$", re.IGNORECASE
)
"""Accession GenBank de nucleotideo: duas letras e seis a oito digitos."""

GENBANK_PROTEIN: re.Pattern[str] = re.compile(
    r"^[A-Z]{3}\d{5,7}(?:\.\d+)?$", re.IGNORECASE
)
"""Accession GenBank de proteina: tres letras e cinco a sete digitos."""

GENE_ID_PATTERN: re.Pattern[str] = re.compile(r"^\d+$")
"""GeneID numerico do banco Gene."""

TAXON_ID_PATTERN: re.Pattern[str] = re.compile(r"^\d{1,12}$")
"""Taxon ID NCBI: apenas digitos, ate 12 caracteres."""

MAX_TAXONOMY_QUERY_LENGTH: int = 200
"""Tamanho maximo de um nome cientifico consultado no banco Taxonomy."""

RELATED_LINK_LIMIT: int = 8
"""Numero maximo de identificadores relacionados retornados por tipo de elink."""

LOCUS_LENGTH_PATTERN: re.Pattern[str] = re.compile(
    r"^LOCUS\s+\S+\s+(\d+)\s+(bp|aa|nt)\b",
    re.IGNORECASE | re.MULTILINE,
)
"""Linha LOCUS de um registro GenBank, da qual se le o comprimento oficial."""


class NCBIQueryError(RuntimeError):
    """Falha classificada ao consultar o NCBI.

    Attributes:
        category: Uma de "invalid_input", "not_found", "wrong_database",
            "source_unavailable", "configuration" ou "internal".
            source_unavailable cobre timeout, HTTP 429, HTTP 5xx e rede.
        failure_kind: Distincao fina para a interface: o mesmo que category,
            ou "timeout" / "rate_limited" quando essa for a causa real.
    """

    def __init__(
        self,
        message: str,
        category: str,
        failure_kind: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.category = category
        self.failure_kind = failure_kind or category


def api_key_from_environment() -> Optional[str]:
    """Le NCBI_API_KEY do ambiente, sem registrá-la.

    Args:
        Nenhum.

    Returns:
        Chave nao vazia ou None. Nunca devolve espacos isolados.

    Raises:
        Nenhum.
    """
    key = (os.environ.get("NCBI_API_KEY") or os.environ.get("ENTREZ_API_KEY") or "").strip()
    return key or None


def classify_identifier(query: str) -> Dict[str, object]:
    """Classifica um identificador NCBI pelo formato, sem consultar a rede.

    Distingue GeneID numerico, accessions de nucleotideo e proteina, projetos WGS,
    accessions de assembly, nomes de organismo e texto que nao e accession.

    Args:
        query: Texto informado pelo usuario (accession, GeneID ou nome).

    Returns:
        Dicionario com "raw" (str), "normalized" (str), "kind" (str),
        "suggested_db" (str ou None), "has_version" (bool) e "reason" (str).
        kind e um de: "empty", "gene_id", "nucleotide", "protein", "wgs_master",
        "assembly", "organism_name", "unknown", "invalid".

    Raises:
        Nenhum.

    Nota biologica:
        O prefixo do accession indica o banco Entrez correto. Um GeneID numerico
        nao deve ser enviado a nucleotide, porque nesse banco o mesmo numero e
        interpretado como GI e pode devolver um registro de outro organismo.
    """
    raw = query or ""
    stripped = raw.strip()
    if not stripped:
        return {
            "raw": raw,
            "normalized": "",
            "kind": "empty",
            "suggested_db": None,
            "has_version": False,
            "reason": "O identificador esta vazio.",
        }
    if len(stripped) > MAX_IDENTIFIER_LENGTH:
        return {
            "raw": raw,
            "normalized": "",
            "kind": "invalid",
            "suggested_db": None,
            "has_version": False,
            "reason": (
                f"O identificador excede {MAX_IDENTIFIER_LENGTH} caracteres e "
                "nao sera enviado ao NCBI."
            ),
        }

    if any(ch.isspace() for ch in stripped):
        return {
            "raw": raw,
            "normalized": stripped,
            "kind": "organism_name",
            "suggested_db": None,
            "has_version": False,
            "reason": (
                "O texto parece um nome de organismo ou frase, nao um accession "
                "NCBI nem um GeneID."
            ),
        }

    normalized = stripped.upper()
    has_version = "." in normalized and normalized.rsplit(".", 1)[-1].isdigit()

    if GENE_ID_PATTERN.fullmatch(stripped):
        if len(stripped) > MAX_GENE_ID_DIGITS:
            return {
                "raw": raw,
                "normalized": stripped,
                "kind": "invalid",
                "suggested_db": None,
                "has_version": False,
                "reason": (
                    f"GeneID numerico excede {MAX_GENE_ID_DIGITS} digitos e nao "
                    "sera enviado ao NCBI."
                ),
            }
        return {
            "raw": raw,
            "normalized": stripped,
            "kind": "gene_id",
            "suggested_db": "gene",
            "has_version": False,
            "reason": "Identificador numerico interpretado como NCBI GeneID.",
        }

    if any(normalized.startswith(prefix) for prefix in ASSEMBLY_PREFIXES):
        return {
            "raw": raw,
            "normalized": normalized,
            "kind": "assembly",
            "suggested_db": None,
            "has_version": has_version,
            "reason": (
                "Accession de assembly (GCA_/GCF_). O HelixScope consulta "
                "nucleotide, protein e gene, nao o banco Assembly."
            ),
        }

    refseq = REFSEQ_PATTERN.fullmatch(normalized)
    if refseq:
        prefix = refseq.group(1).upper() + "_"
        if prefix in PROTEIN_REFSEQ_PREFIXES:
            return {
                "raw": raw,
                "normalized": normalized,
                "kind": "protein",
                "suggested_db": "protein",
                "has_version": has_version,
                "reason": f"Accession RefSeq de proteina ({prefix}).",
            }
        if prefix in NUCLEOTIDE_REFSEQ_PREFIXES:
            return {
                "raw": raw,
                "normalized": normalized,
                "kind": "nucleotide",
                "suggested_db": "nucleotide",
                "has_version": has_version,
                "reason": f"Accession RefSeq de nucleotideo ({prefix}).",
            }

    if WGS_FOUR_PATTERN.fullmatch(normalized) or WGS_SIX_PATTERN.fullmatch(normalized):
        return {
            "raw": raw,
            "normalized": normalized,
            "kind": "wgs_master",
            "suggested_db": "nucleotide",
            "has_version": has_version,
            "reason": "Accession de projeto WGS no banco nucleotide.",
        }

    if GENBANK_PROTEIN.fullmatch(normalized):
        return {
            "raw": raw,
            "normalized": normalized,
            "kind": "protein",
            "suggested_db": "protein",
            "has_version": has_version,
            "reason": "Accession GenBank de proteina (tres letras + digitos).",
        }

    if GENBANK_NUC_SHORT.fullmatch(normalized) or GENBANK_NUC_LONG.fullmatch(normalized):
        return {
            "raw": raw,
            "normalized": normalized,
            "kind": "nucleotide",
            "suggested_db": "nucleotide",
            "has_version": has_version,
            "reason": "Accession GenBank de nucleotideo.",
        }

    if re.fullmatch(r"[A-Z][A-Z0-9_.]+", normalized):
        return {
            "raw": raw,
            "normalized": normalized,
            "kind": "unknown",
            "suggested_db": None,
            "has_version": has_version,
            "reason": (
                "Formato nao classificado com certeza; o banco selecionado na "
                "interface sera usado."
            ),
        }

    return {
        "raw": raw,
        "normalized": normalized,
        "kind": "invalid",
        "suggested_db": None,
        "has_version": has_version,
        "reason": "O texto nao corresponde a um accession NCBI nem a um GeneID.",
    }


def resolve_database(query: str, requested_db: str) -> Tuple[str, Dict[str, object], List[str]]:
    """Escolhe o banco Entrez adequado ao identificador.

    Args:
        query: Accession, GeneID ou texto informado.
        requested_db: Banco pedido pela interface ("nucleotide", "protein",
            "gene" ou "pubmed").

    Returns:
        Tupla (banco resolvido, classificacao de classify_identifier, avisos).

    Raises:
        NCBIQueryError: Se o identificador for vazio, invalido, for um nome de
            organismo, for um accession de assembly, ou se requested_db nao for
            um dos bancos suportados.

    Nota biologica:
        Consultar o banco errado nao e neutro: um GeneID em nucleotide e lido
        como GI e pode devolver uma sequencia de outra especie.
    """
    requested = (requested_db or "nucleotide").strip().lower()
    if requested not in SUPPORTED_DATABASES:
        raise NCBIQueryError(
            f"Banco '{requested_db}' nao e suportado. Use nucleotide, protein, gene ou pubmed.",
            "wrong_database",
        )

    info = classify_identifier(query)
    kind = str(info["kind"])
    warnings: List[str] = []

    if kind == "empty":
        raise NCBIQueryError("O numero de acesso nao pode ser vazio.", "invalid_input")
    if kind == "organism_name":
        raise NCBIQueryError(
            f"'{query.strip()}' nao e um accession. Informe um accession NCBI "
            "(por exemplo NM_000518.5, NP_000509.1 ou AGTM000000000.1) ou um "
            "GeneID numerico, nao o nome da especie.",
            "invalid_input",
        )
    if kind == "assembly":
        raise NCBIQueryError(
            f"{info['normalized']} e um accession de assembly. O HelixScope nao "
            "consulta o banco Assembly. Use o accession nucleotide do genoma ou "
            "de um contig/scaffold WGS.",
            "wrong_database",
        )
    if kind == "invalid":
        raise NCBIQueryError(
            f"'{query.strip()}' nao e um accession NCBI reconhecivel nem um "
            "GeneID. Nenhum registro foi consultado.",
            "invalid_input",
        )

    if requested == "pubmed":
        if kind == "gene_id":
            pmid_info = dict(info)
            pmid_info["kind"] = "pmid"
            pmid_info["suggested_db"] = "pubmed"
            pmid_info["reason"] = (
                "Identificador numerico interpretado como PMID porque o banco "
                "pedido e pubmed, nao Gene."
            )
            return "pubmed", pmid_info, []
        raise NCBIQueryError(
            "PubMed fetch espera um PMID numerico. Accessions de nucleotideo "
            "ou proteina nao sao registros PubMed.",
            "wrong_database",
        )

    suggested = str(info["suggested_db"] or requested)
    if info["suggested_db"] and suggested != requested:
        warnings.append(
            f"Identificador classificado como {kind}; consultando '{suggested}' "
            f"em vez de '{requested}'. {info['reason']}"
        )
    return suggested, info, warnings


def fetch_by_accession(
    accession: str,
    email: str,
    db: str = "nucleotide",
    api_key: Optional[str] = None,
) -> Dict[str, object]:
    """Recupera um registro do NCBI pelo identificador e o interpreta.

    Classifica o identificador, escolhe o banco Entrez correto e consulta
    nucleotide/protein via efetch (GenBank), gene via esummary/elink, ou
    pubmed via esummary (e abstract via efetch quando o NCBI indica abstract).
    Accessions sem versao sao aceitos: o NCBI devolve a versao corrente.

    Args:
        accession: Numero de acesso, GeneID ou identificador do registro.
        email: E-mail valido exigido pelo Entrez para identificar o usuario.
        db: Banco pedido pela interface; e substituido pelo banco inferido
            quando o formato do identificador indica outro.
        api_key: Chave de API opcional do NCBI; aumenta o limite de requisicoes.

    Returns:
        Dicionario com "accession", "description", "organism", "length",
        "sequence", "features", alem de "database", "record_kind",
        "sequence_available", "identifier_kind", "warnings", "related_records"
        e, para gene, "gene". Cada feature tem "type", "location" e "qualifiers".

    Raises:
        ValueError: Se o e-mail for vazio.
        NCBIQueryError: Se o identificador for invalido, o registro nao existir,
            o banco for inadequado, a fonte estiver indisponivel ou o parse
            falhar.
        RuntimeError: Se o Biopython nao estiver instalado.

    Nota biologica:
        O numero de acesso identifica de forma estavel e versionada uma
        sequencia depositada. GeneIDs identificam o gene, nao uma sequencia:
        a sequencia correspondente vem de accessions nucleotide/protein ligados.
    """
    if not email or not email.strip():
        raise ValueError("O e-mail e exigido pelo NCBI Entrez.")
    if not _valid_entrez_email(email):
        raise ValueError(
            "Informe um e-mail de contato valido para o NCBI Entrez "
            "(formato local@dominio)."
        )

    resolved_db, info, warnings = resolve_database(accession, db)
    identifier = str(info["normalized"])
    _configure_entrez(email, api_key)

    if resolved_db == "gene":
        record = _fetch_gene_record(identifier)
    elif resolved_db == "pubmed":
        record = _fetch_pubmed_record(identifier)
    else:
        record = _fetch_sequence_record(identifier, resolved_db)

    record["database"] = resolved_db
    record["identifier_kind"] = info["kind"]
    record["query"] = accession.strip()
    record["source"] = "NCBI Entrez"
    record["retrieved_at_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    record["warnings"] = list(record.get("warnings") or []) + warnings
    record.setdefault("related_records", [])
    record.setdefault("gene", None)
    record.setdefault("genbank_date", "")
    record.setdefault("taxonomy", [])
    record.setdefault("molecule", "")
    record.setdefault("topology", "")
    record.setdefault("references", [])
    record.setdefault("comment", "")
    return record


def search_records(
    term: str,
    email: str,
    db: str = "nucleotide",
    api_key: Optional[str] = None,
    retmax: int = 15,
    retstart: int = 0,
) -> Dict[str, object]:
    """Busca registros reais no NCBI via Entrez esearch + esummary.

    Nao inventa hits. Lista vazia significa zero IDs devolvidos pelo NCBI para
    aquele termo, nao falha de rede. Falhas de rede, timeout e HTTP 429 sao
    NCBIQueryError.

    Args:
        term: Consulta Entrez (por exemplo 'BRCA1 AND human[orgn]').
        email: E-mail de contato exigido pelo NCBI.
        db: Banco Entrez: nucleotide, protein, gene ou pubmed.
        api_key: Chave opcional do NCBI.
        retmax: Maximo de resumos a recuperar (1 a MAX_SEARCH_RETMAX).
        retstart: Offset 0-based no resultado Entrez para paginacao.

    Returns:
        Dicionario com "query", "database", "count" (total NCBI), "records"
        (lista de resumos realmente devolvidos), "retstart", "retmax",
        "source", "retrieved_at_utc" e "status" ("retrieved" ou "no_records_found").

    Raises:
        ValueError: Se o e-mail for invalido.
        NCBIQueryError: Se o termo for vazio/longo, o banco for invalido, a
            resposta for malformada, ou a fonte estiver indisponivel.

    Nota biologica:
        esearch devolve IDs; esummary devolve metadados. A sequencia completa
        so e obtida por fetch_by_accession. BLAST nao e esta funcao.
    """
    if not email or not email.strip():
        raise ValueError("O e-mail e exigido pelo NCBI Entrez.")
    if not _valid_entrez_email(email):
        raise ValueError(
            "Informe um e-mail de contato valido para o NCBI Entrez "
            "(formato local@dominio)."
        )
    cleaned = (term or "").strip()
    if not cleaned:
        raise NCBIQueryError("The search query cannot be empty.", "invalid_input")
    if len(cleaned) > MAX_SEARCH_TERM_LENGTH:
        raise NCBIQueryError(
            f"Search query exceeds {MAX_SEARCH_TERM_LENGTH} characters.",
            "invalid_input",
        )
    database = (db or "nucleotide").strip().lower()
    if database not in SEARCHABLE_DATABASES:
        raise NCBIQueryError(
            f"Database '{db}' is not searchable here. Use nucleotide, protein, "
            "gene or pubmed.",
            "wrong_database",
        )
    try:
        limit = int(retmax)
    except (TypeError, ValueError) as exc:
        raise NCBIQueryError("retmax must be an integer.", "invalid_input") from exc
    limit = max(1, min(limit, MAX_SEARCH_RETMAX))
    try:
        start = int(retstart)
    except (TypeError, ValueError) as exc:
        raise NCBIQueryError("retstart must be an integer.", "invalid_input") from exc
    if start < 0:
        raise NCBIQueryError("retstart cannot be negative.", "invalid_input")
    _configure_entrez(email, api_key)

    from Bio import Entrez

    def _esearch() -> dict:
        handle = Entrez.esearch(
            db=database,
            term=cleaned,
            retmax=limit,
            retstart=start,
            retmode="xml",
        )
        try:
            payload = Entrez.read(handle)
        finally:
            handle.close()
        if not isinstance(payload, dict):
            raise NCBIQueryError(
                "Malformed NCBI search response (esearch).",
                "internal",
            )
        return payload

    payload = _entrez_retry(_esearch, cleaned, database)
    try:
        count = int(payload.get("Count", 0))
    except (TypeError, ValueError):
        count = 0
    raw_ids = payload.get("IdList") or []
    if isinstance(raw_ids, str):
        id_list = [raw_ids] if raw_ids else []
    else:
        id_list = [str(item) for item in raw_ids if str(item).strip()]

    retrieved_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if not id_list:
        return {
            "query": cleaned,
            "database": database,
            "count": count,
            "records": [],
            "retstart": start,
            "retmax": limit,
            "source": "NCBI Entrez",
            "retrieved_at_utc": retrieved_at,
            "status": "no_records_found",
        }

    def _esummary() -> object:
        handle = Entrez.esummary(db=database, id=",".join(id_list), retmode="xml")
        try:
            return Entrez.read(handle)
        finally:
            handle.close()

    summaries = _entrez_retry(_esummary, cleaned, database)
    parsed = list(_iter_document_summaries(summaries))
    records: List[dict] = []
    for index, item_id in enumerate(id_list):
        if index < len(parsed):
            records.append(
                _search_hit_from_summary(parsed[index], database, item_id)
            )
            continue
        records.append(
            {
                "id": item_id,
                "accession": item_id,
                "title": "Not available",
                "organism": "Not available",
                "length": None,
                "database": database,
                "extra": "esummary did not return a parseable document for this ID.",
            }
        )
    return {
        "query": cleaned,
        "database": database,
        "count": count,
        "records": records,
        "retstart": start,
        "retmax": limit,
        "source": "NCBI Entrez",
        "retrieved_at_utc": retrieved_at,
        "status": "retrieved",
    }


def fetch_taxonomy(
    query: str,
    email: str,
    api_key: Optional[str] = None,
) -> Dict[str, object]:
    """Recupera um taxon real do NCBI Taxonomy (esearch + efetch XML).

    Args:
        query: Taxon ID numerico ou nome cientifico. Nao e URL.
        email: E-mail de contato exigido pelo Entrez.
        api_key: Chave NCBI opcional.

    Returns:
        Dict com taxon_id, scientific_name, rank, lineage (lista de
        {taxon_id, scientific_name, rank}), common_name somente se o NCBI
        fornecer, source, retrieved_at_utc, status RETRIEVED.

    Raises:
        ValueError: E-mail invalido.
        NCBIQueryError: Query invalida, taxon inexistente, XML malformado
            ou NCBI indisponivel.
        RuntimeError: Biopython ausente.

    Nota biologica:
        Isto e classificacao taxonomica NCBI, nao inferencia filogenetica
        a partir de sequencias. Lineage vem de LineageEx quando presente.
    """
    if not email or not email.strip():
        raise ValueError("O e-mail e exigido pelo NCBI Entrez.")
    if not _valid_entrez_email(email):
        raise ValueError(
            "Informe um e-mail de contato valido para o NCBI Entrez "
            "(formato local@dominio)."
        )
    cleaned = (query or "").strip()
    if not cleaned:
        raise NCBIQueryError("Taxonomy query cannot be empty.", "invalid_input")
    if len(cleaned) > MAX_TAXONOMY_QUERY_LENGTH:
        raise NCBIQueryError(
            f"Taxonomy query exceeds {MAX_TAXONOMY_QUERY_LENGTH} characters.",
            "invalid_input",
        )
    lowered = cleaned.lower()
    if "://" in lowered or ".." in cleaned or "<" in cleaned or ">" in cleaned:
        raise NCBIQueryError(
            "Taxonomy query cannot be a URL or markup.",
            "invalid_input",
        )
    _configure_entrez(email, api_key)
    from Bio import Entrez

    tax_id = cleaned if TAXON_ID_PATTERN.fullmatch(cleaned) else ""
    if not tax_id:
        def _esearch() -> dict:
            handle = Entrez.esearch(
                db="taxonomy",
                term=cleaned,
                retmax=1,
                retmode="xml",
            )
            try:
                payload = Entrez.read(handle)
            finally:
                handle.close()
            if not isinstance(payload, dict):
                raise NCBIQueryError(
                    "Malformed NCBI taxonomy search response.",
                    "internal",
                )
            return payload

        payload = _entrez_retry(_esearch, cleaned, "taxonomy")
        raw_ids = payload.get("IdList") or []
        if isinstance(raw_ids, str):
            id_list = [raw_ids] if raw_ids else []
        else:
            id_list = [str(item) for item in raw_ids if str(item).strip()]
        if not id_list:
            raise NCBIQueryError(
                f"No NCBI Taxonomy record for '{cleaned}'.",
                "not_found",
            )
        tax_id = id_list[0]
        if not TAXON_ID_PATTERN.fullmatch(tax_id):
            raise NCBIQueryError(
                "NCBI Taxonomy returned a non-numeric TaxId.",
                "internal",
            )

    def _efetch() -> object:
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        try:
            return Entrez.read(handle)
        finally:
            handle.close()

    raw = _entrez_retry(_efetch, tax_id, "taxonomy")
    record = _first_taxonomy_record(raw)
    parsed = _parse_taxonomy_record(record)
    parsed["query"] = cleaned
    parsed["source"] = "NCBI Taxonomy"
    parsed["database"] = "taxonomy"
    parsed["retrieved_at_utc"] = datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    parsed["status"] = "RETRIEVED"
    return parsed


def homology_search_availability() -> dict:
    """Disponibilidade de BLAST: NCBI Common URL API, sem BLAST+ local.

    Args:
        Nenhum.

    Returns:
        Dict de blast_search.ncbi_blast_availability().

    Raises:
        Nenhum.
    """
    from . import blast_search

    return blast_search.ncbi_blast_availability()


def parse_genbank_location(location: str) -> List[dict]:
    """Interpreta a representacao textual de uma localizacao de feature GenBank.

    Reconhece localizacoes simples ("[0:1234](+)") e compostas
    ("join{[0:100](+), [200:300](+)}"), incluindo limites incertos marcados por
    "<" e ">". As coordenadas sao mantidas como no Biopython: base 0, intervalo
    semiaberto e sempre referidas a fita direta (plus).

    Args:
        location: String de localizacao, tal como armazenada no campo "location"
            de cada feature retornada por fetch_by_accession().

    Returns:
        Lista de dicionarios com "start" (int), "end" (int) e "strand" (str, "+"
        ou "-"), na ordem em que os segmentos aparecem na localizacao. Lista vazia
        quando nenhum segmento pode ser reconhecido.

    Raises:
        Nenhum.

    Nota biologica:
        Genes eucarioticos e alguns genes virais sao descritos por localizacoes
        compostas porque o CDS maduro resulta da juncao de exons; ignorar os
        segmentos e usar apenas o intervalo externo incluiria introns e
        corromperia a leitura em codons.
    """
    segments: List[dict] = []
    for match in LOCATION_SEGMENT_PATTERN.finditer(location or ""):
        start = int(match.group(1))
        end = int(match.group(2))
        strand_symbol = match.group(3)
        if end <= start:
            continue
        segments.append(
            {
                "start": start,
                "end": end,
                "strand": "-" if strand_symbol == "-" else "+",
            }
        )
    return segments


def extract_cds_regions(record: Dict[str, object]) -> List[dict]:
    """Extrai as regioes codificadoras (CDS) anotadas de um registro do NCBI.

    Percorre as features do registro, seleciona as de tipo "CDS" e converte cada
    localizacao em segmentos de coordenadas, preservando gene, produto,
    identificador de proteina e o deslocamento de leitura (codon_start).

    Args:
        record: Dicionario retornado por fetch_by_accession(), contendo a chave
            "features".

    Returns:
        Lista de dicionarios, um por CDS, com "gene" (str), "product" (str),
        "protein_id" (str), "strand" (str), "segments" (list de tuplas
        (start, end) em coordenadas da fita direta, base 0 e semiabertas),
        "start" (int), "end" (int), "length_nt" (int) e "codon_start" (int, 1 a
        3). Lista vazia quando o registro nao possui CDS anotados.

    Raises:
        Nenhum.

    Nota biologica:
        A anotacao oficial e a fonte mais confiavel de regioes codificadoras
        porque resulta de evidencia experimental e curadoria, e nao de predicao
        estatistica. Um mesmo genoma pode conter CDS sobrepostos e em ambas as
        fitas, algo que nenhuma regra simples de predicao reproduz.
    """
    features = record.get("features") or []
    if not isinstance(features, list):
        return []

    regions: List[dict] = []
    for feature in features:
        if not isinstance(feature, dict):
            continue
        if str(feature.get("type", "")).upper() != "CDS":
            continue
        segments = parse_genbank_location(str(feature.get("location", "")))
        if not segments:
            continue

        qualifiers = feature.get("qualifiers") or {}
        if not isinstance(qualifiers, dict):
            qualifiers = {}

        codon_start = 1
        raw_codon_start = qualifiers.get("codon_start")
        if isinstance(raw_codon_start, list) and raw_codon_start:
            raw_codon_start = raw_codon_start[0]
        try:
            parsed_codon_start = int(str(raw_codon_start))
            if 1 <= parsed_codon_start <= 3:
                codon_start = parsed_codon_start
        except (TypeError, ValueError):
            codon_start = 1

        spans = [(segment["start"], segment["end"]) for segment in segments]
        strand = "-" if any(segment["strand"] == "-" for segment in segments) else "+"
        regions.append(
            {
                "gene": _first_qualifier(qualifiers, "gene"),
                "product": _first_qualifier(qualifiers, "product"),
                "protein_id": _first_qualifier(qualifiers, "protein_id"),
                "strand": strand,
                "segments": spans,
                "start": min(span[0] for span in spans),
                "end": max(span[1] for span in spans),
                "length_nt": sum(span[1] - span[0] for span in spans),
                "codon_start": codon_start,
            }
        )
    return regions


TRACKED_FEATURE_TYPES: frozenset[str] = frozenset(
    {"gene", "cds", "exon", "intron", "misc_feature"}
)
"""Tipos GenBank inspecionados pelo explorador de features. Outros tipos nao
sao promovidos a regioes de analise."""


def split_accession_version(accession: str) -> Tuple[str, str]:
    """Separa o accession da versao NCBI quando o sufixo e numerico.

    Args:
        accession: Identificador, com ou sem versao (ex.: NM_007294.4).

    Returns:
        Tupla (accession_sem_versao, accession_versionado). A versao e string
        vazia quando o identificador nao traz sufixo .N.

    Raises:
        Nenhum.
    """
    text = str(accession or "").strip()
    if not text:
        return "", ""
    if "." in text:
        base, suffix = text.rsplit(".", 1)
        if base and suffix.isdigit():
            return base, text
    return text, ""


def record_explorer_fields(record: Dict[str, object]) -> dict:
    """Campos de um registro NCBI que realmente existem neste payload.

    Nao inclui chaves vazias. Nao inventa taxonomia, referencias ou sequencia.

    Args:
        record: Dicionario retornado por fetch_by_accession().

    Returns:
        Dicionario possivelmente com accession, version, definition, organism,
        taxonomy, molecule, topology, length, sequence, sequence_available,
        features, references e links.

    Raises:
        Nenhum.
    """
    accession_full = str(record.get("accession") or "")
    acc, version = split_accession_version(accession_full)
    fields: Dict[str, object] = {}
    if acc:
        fields["accession"] = acc
    if version:
        fields["version"] = version
    definition = str(record.get("description") or "")
    if definition:
        fields["definition"] = definition
    organism = str(record.get("organism") or "")
    if organism:
        fields["organism"] = organism
    taxonomy = record.get("taxonomy") or []
    if isinstance(taxonomy, list):
        cleaned = [str(item).strip() for item in taxonomy if str(item).strip()]
        if cleaned:
            fields["taxonomy"] = cleaned
    molecule = str(record.get("molecule") or "")
    if molecule:
        fields["molecule"] = molecule
    topology = str(record.get("topology") or "")
    if topology:
        fields["topology"] = topology
    try:
        length = int(record.get("length") or 0)
    except (TypeError, ValueError):
        length = 0
    if length > 0:
        fields["length"] = length
    sequence = str(record.get("sequence") or "")
    available = bool(record.get("sequence_available") and sequence)
    fields["sequence_available"] = available
    if available:
        fields["sequence"] = sequence
    features = record.get("features")
    if isinstance(features, list) and features:
        fields["features"] = features
    references = record.get("references")
    if isinstance(references, list) and references:
        fields["references"] = references
    links = record.get("related_records")
    if isinstance(links, list) and links:
        fields["links"] = links
    return fields


def extract_annotated_features(record: Dict[str, object]) -> List[dict]:
    """Extrai gene, CDS, exon, intron e misc_feature com coordenadas e fita.

    Join e complement sao preservados. Se a localizacao nao puder ser lida ou
    cair fora da sequencia recuperada, o status e UNAVAILABLE e nenhuma
    sequencia e inventada.

    Args:
        record: Dicionario retornado por fetch_by_accession().

    Returns:
        Lista de dicionarios com type, location, strand, is_join,
        is_complement, segments, start, end, length_nt, gene, product, status,
        unavailable_reason e region_sequence (fita +, intervalos concatenados
        quando recuperaveis).

    Raises:
        Nenhum.

    Nota biologica:
        region_sequence e a concatenacao dos intervalos na fita direta do
        registro, nao o transcrito maduro nem o reverso-complemento. Nao
        interpreta splicing.
    """
    features = record.get("features") or []
    if not isinstance(features, list):
        return []
    sequence = str(record.get("sequence") or "")
    if "sequence_available" in record:
        sequence_available = bool(record.get("sequence_available") and sequence)
    else:
        sequence_available = bool(sequence)
    parsed: List[dict] = []
    for feature in features:
        if not isinstance(feature, dict):
            continue
        ftype = str(feature.get("type") or "").strip()
        if ftype.lower() not in TRACKED_FEATURE_TYPES:
            continue
        location = str(feature.get("location") or "")
        segments = parse_genbank_location(location)
        is_join = "join" in location.lower() or len(segments) > 1
        is_complement = (
            "complement" in location.lower()
            or any(segment["strand"] == "-" for segment in segments)
        )
        status = "AVAILABLE"
        reason = ""
        if not segments:
            status = "UNAVAILABLE"
            reason = (
                "Feature location could not be parsed. No sequence region "
                "was invented."
            )
        elif sequence_available:
            for segment in segments:
                start = int(segment["start"])
                end = int(segment["end"])
                if not (0 <= start < end <= len(sequence)):
                    status = "UNAVAILABLE"
                    reason = (
                        "Feature coordinates fall outside the retrieved "
                        "sequence. The region was not reconstructed."
                    )
                    break
        else:
            status = "UNAVAILABLE"
            reason = (
                "Sequence bases were not retrieved for this record, so the "
                "feature region cannot be recovered."
            )
        qualifiers = feature.get("qualifiers") or {}
        if not isinstance(qualifiers, dict):
            qualifiers = {}
        spans = [(int(item["start"]), int(item["end"])) for item in segments]
        strand = "-" if is_complement else "+"
        region_sequence = ""
        if status == "AVAILABLE" and sequence_available:
            region_sequence = "".join(sequence[start:end] for start, end in spans)
        parsed.append(
            {
                "type": ftype,
                "location": location,
                "strand": strand if segments else "",
                "is_join": is_join,
                "is_complement": is_complement,
                "segments": spans,
                "start": min((span[0] for span in spans), default=None),
                "end": max((span[1] for span in spans), default=None),
                "length_nt": sum(span[1] - span[0] for span in spans) if spans else None,
                "gene": _first_qualifier(qualifiers, "gene"),
                "product": _first_qualifier(qualifiers, "product"),
                "locus_tag": _first_qualifier(qualifiers, "locus_tag"),
                "note": _first_qualifier(qualifiers, "note"),
                "status": status,
                "unavailable_reason": reason,
                "region_sequence": region_sequence,
                "qualifiers": qualifiers,
            }
        )
    return parsed


def sequence_for_workspace(record: Dict[str, object]) -> dict:
    """Copia a sequencia NCBI para o workspace sem alterar o conteudo.

    Args:
        record: Dicionario retornado por fetch_by_accession().

    Returns:
        Dict com sequence (verbatim), molecule, accession, version, source,
        organism, hash, timestamp e length.

    Raises:
        ValueError: Se a sequencia nao tiver sido recuperada.
    """
    from . import provenance

    sequence = str(record.get("sequence") or "")
    if not record.get("sequence_available") or not sequence:
        raise ValueError(
            "This NCBI record has no retrieved sequence. Nothing was copied "
            "to the workspace and no bases were invented."
        )
    accession_full = str(record.get("accession") or "")
    acc, version = split_accession_version(accession_full)
    database = str(record.get("database") or "").strip().lower()
    if database == "protein":
        molecule = "PROTEIN"
    else:
        molecule = "DNA"
    retrieved = str(record.get("retrieved_at_utc") or provenance.utc_now())
    return {
        "sequence": sequence,
        "molecule": molecule,
        "accession": acc or accession_full,
        "version": version or accession_full,
        "source": str(record.get("source") or "NCBI Entrez"),
        "organism": str(record.get("organism") or ""),
        "hash": provenance.sequence_digest(sequence),
        "timestamp": retrieved,
        "length": len(sequence),
    }


def feature_spans_for_analysis(
    feature: Dict[str, object], sequence: str
) -> List[dict]:
    """Converte segmentos recuperaveis em spans para CRISPR/3D posteriores.

    Args:
        feature: Item de extract_annotated_features().
        sequence: Sequencia do registro (fita +), inalterada.

    Returns:
        Lista de dicts scientific_checks.make_span, um por segmento.

    Raises:
        ValueError: Se a feature estiver UNAVAILABLE ou as coordenadas
            nao couberem na sequencia.
    """
    from . import scientific_checks

    if str(feature.get("status") or "") != "AVAILABLE":
        raise ValueError(
            str(feature.get("unavailable_reason") or "Feature region unavailable.")
        )
    spans = feature.get("segments") or []
    if not spans:
        raise ValueError("Feature has no recoverable coordinates.")
    built: List[dict] = []
    strand = str(feature.get("strand") or "+")
    if strand not in {"+", "-"}:
        strand = "+"
    label = str(feature.get("gene") or feature.get("type") or "feature")
    for start, end in spans:
        built.append(
            scientific_checks.make_span(
                start=int(start),
                end=int(end),
                sequence=sequence,
                strand=strand,
                source="NCBI GenBank feature",
                label=label,
                kind=str(feature.get("type") or "feature"),
                status="RETRIEVED",
            )
        )
    return built


def _first_qualifier(qualifiers: Dict[str, object], name: str) -> str:
    """Retorna o primeiro valor textual de um qualifier GenBank (uso interno).

    Args:
        qualifiers: Dicionario de qualifiers de uma feature.
        name: Nome do qualifier desejado.

    Returns:
        Primeiro valor como string, ou string vazia quando ausente.

    Nota biologica:
        Qualifiers GenBank sao listas porque uma feature pode acumular varios
        sinonimos ou anotacoes historicas para o mesmo campo.
    """
    value = qualifiers.get(name)
    if isinstance(value, list):
        return str(value[0]) if value else ""
    if value is None:
        return ""
    return str(value)


def _valid_entrez_email(email: str) -> bool:
    """Verifica se o e-mail de contato tem forma minima aceitavel (interno)."""
    text = (email or "").strip()
    if not text or len(text) > MAX_EMAIL_LENGTH:
        return False
    if text.count("@") != 1:
        return False
    local, domain = text.split("@", 1)
    if not local or not domain or "." not in domain:
        return False
    if any(ch.isspace() for ch in text):
        return False
    if domain.startswith(".") or domain.endswith("."):
        return False
    return True


def _entrez_urlopen(url, data=None, timeout=None, *, context=None):
    """urlopen do Entrez com timeout explicito (uso interno)."""
    from urllib.request import urlopen

    if timeout is None:
        timeout = ENTREZ_TIMEOUT_S
    return urlopen(url, data=data, timeout=timeout, context=context)


def _configure_entrez(email: str, api_key: Optional[str]) -> None:
    """Configura Entrez.email, Entrez.tool e a pausa de taxa (uso interno)."""
    try:
        from Bio import Entrez
    except ImportError as exc:
        raise RuntimeError(
            "Biopython nao esta instalado; instale a dependencia 'biopython'."
        ) from exc

    Entrez.email = email.strip()
    Entrez.tool = ENTREZ_TOOL
    Entrez.urlopen = _entrez_urlopen
    if api_key:
        Entrez.api_key = api_key.strip()
        time.sleep(SLEEP_WITH_KEY)
    else:
        time.sleep(SLEEP_WITHOUT_KEY)


def _summary_length(summary: Dict[str, object]) -> int:
    """Extrai o comprimento numerico de um DocumentSummary Entrez (interno)."""
    raw = summary.get("Length")
    if raw is None:
        return 0
    try:
        return int(raw)
    except (TypeError, ValueError):
        return 0


def _locus_length_from_genbank_text(raw: str) -> int:
    """Le o comprimento declarado na linha LOCUS de um texto GenBank (interno)."""
    match = LOCUS_LENGTH_PATTERN.search(raw or "")
    if not match:
        return 0
    return int(match.group(1))


def _peek_reported_length(accession: str, db: str) -> int:
    """Consulta so o cabecalho (1 resíduo) para obter o comprimento LOCUS.

    Evita baixar scaffolds WGS inteiros quando o esummary falha. Nao inventa
    bases: se o peek falhar, devolve 0 e o chamador decide se o efetch completo
    ainda e seguro.
    """
    from Bio import Entrez

    try:
        time.sleep(SLEEP_WITHOUT_KEY)
        handle = Entrez.efetch(
            db=db,
            id=accession,
            rettype="gb",
            retmode="text",
            seq_start=1,
            seq_stop=1,
        )
        try:
            raw = handle.read()
        finally:
            handle.close()
    except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError):
        return 0
    return _locus_length_from_genbank_text(str(raw))


def _fetch_sequence_summary(accession: str, db: str) -> Optional[dict]:
    """Obtem o resumo Entrez (comprimento e titulo) sem baixar a sequencia."""
    from Bio import Entrez

    try:
        handle = Entrez.esummary(db=db, id=accession, retmode="xml")
        try:
            payload = Entrez.read(handle)
        finally:
            handle.close()
    except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError, OSError, ValueError):
        return None
    return _first_gene_summary(payload)


def _large_record_payload(
    accession: str,
    db: str,
    summary: dict,
    reported_length: int,
) -> Dict[str, object]:
    """Monta metadados de um registro grande demais para download completo."""
    from Bio import Entrez, SeqIO
    from io import StringIO

    title = str(summary.get("Title") or summary.get("Caption") or accession)
    version = str(summary.get("AccessionVersion") or accession)
    organism = ""
    genbank_date = ""
    features: List[dict] = []
    try:
        time.sleep(SLEEP_WITHOUT_KEY)
        handle = Entrez.efetch(
            db=db,
            id=accession,
            rettype="gb",
            retmode="text",
            seq_start=1,
            seq_stop=1,
        )
        try:
            raw = handle.read()
        finally:
            handle.close()
        if str(raw).strip():
            header = SeqIO.read(StringIO(str(raw)), "genbank")
            organism = str(header.annotations.get("organism") or "")
            genbank_date = str(header.annotations.get("date") or "")
            title = header.description or title
            version = header.id or version
            for feature in header.features:
                features.append(
                    {
                        "type": feature.type,
                        "location": str(feature.location),
                        "qualifiers": dict(feature.qualifiers),
                    }
                )
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError):
        pass

    unit = "residues" if db == "protein" else "bp"
    return {
        "accession": version,
        "description": title,
        "organism": organism,
        "taxonomy": [],
        "molecule": "",
        "topology": "",
        "length": reported_length,
        "sequence": "",
        "features": features,
        "references": [],
        "record_kind": "large_sequence",
        "sequence_available": False,
        "warnings": [
            f"O registro existe no NCBI ({reported_length:,} {unit}), mas o "
            f"HelixScope nao baixa sequencias acima de {MAX_SEQUENCE_FETCH_NT:,} "
            f"{unit}. Nenhum nucleotideo foi inventado. Use um accession de "
            "mRNA/CDS (por exemplo NM_*) ou um contig menor."
        ],
        "related_records": [],
        "gene": None,
        "genbank_date": genbank_date,
        "comment": "",
    }


def _fetch_sequence_record(accession: str, db: str) -> Dict[str, object]:
    """Busca um registro GenBank/GenPept e o converte em dicionario (interno)."""
    from Bio import Entrez, SeqIO

    summary = _fetch_sequence_summary(accession, db)
    reported_length = _summary_length(summary) if summary is not None else 0
    if reported_length <= 0:
        reported_length = _peek_reported_length(accession, db)
    if reported_length > MAX_SEQUENCE_FETCH_NT:
        return _large_record_payload(
            accession, db, summary or {}, reported_length
        )

    try:
        UndefinedSequenceError = _undefined_sequence_error_type()
        time.sleep(SLEEP_WITHOUT_KEY)
        handle = Entrez.efetch(
            db=db, id=accession, rettype="gb", retmode="text"
        )
        try:
            raw = handle.read()
        finally:
            handle.close()
        if not str(raw).strip():
            raise NCBIQueryError(
                f"O NCBI devolveu uma resposta vazia para '{accession}' em {db}.",
                "source_unavailable",
            )
        from io import StringIO

        record = SeqIO.read(StringIO(str(raw)), "genbank")
    except NCBIQueryError:
        raise
    except urllib.error.HTTPError as exc:
        raise _http_error(accession, db, exc) from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise NCBIQueryError(
            f"NCBI indisponivel ao buscar '{accession}' em {db}: {exc}.",
            "source_unavailable",
        ) from exc
    except ValueError as exc:
        message = str(exc)
        if "No records found" in message:
            raise NCBIQueryError(
                f"Nenhum registro encontrado para '{accession}' no banco {db}. "
                "O identificador nao existe nesse banco ou foi retirado.",
                "not_found",
            ) from exc
        raise NCBIQueryError(
            f"Nao foi possivel interpretar o registro '{accession}' em {db}: {exc}.",
            "internal",
        ) from exc

    features: List[dict] = []
    for feature in record.features:
        features.append(
            {
                "type": feature.type,
                "location": str(feature.location),
                "qualifiers": dict(feature.qualifiers),
            }
        )

    sequence, available = _sequence_payload(record.seq, UndefinedSequenceError)
    warnings: List[str] = []
    wgs = record.annotations.get("wgs")
    keywords = record.annotations.get("keywords") or []
    is_wgs = bool(wgs) or any(str(item).upper() == "WGS" for item in keywords)
    if not available:
        warnings.append(
            "Este registro existe no NCBI, mas o conteudo da sequencia nao esta "
            "depositado neste accession (comum em projetos WGS mestres). "
            "Organismo, definicao e features foram recuperados."
        )
        if isinstance(wgs, (list, tuple)) and len(wgs) >= 2:
            warnings.append(
                f"Contigs WGS associados: {wgs[0]} a {wgs[1]}."
            )
        elif wgs:
            warnings.append(f"Anotacao WGS: {wgs}.")

    organism = record.annotations.get("organism", "")
    taxonomy = record.annotations.get("taxonomy") or []
    if not isinstance(taxonomy, list):
        taxonomy = []
    references: List[dict] = []
    for ref in list(record.annotations.get("references") or [])[:10]:
        pubmed = str(getattr(ref, "pubmed_id", "") or "")
        references.append(
            {
                "title": str(getattr(ref, "title", "") or ""),
                "authors": str(getattr(ref, "authors", "") or ""),
                "journal": str(getattr(ref, "journal", "") or ""),
                "pubmed": pubmed,
            }
        )
    return {
        "accession": record.id,
        "description": record.description,
        "organism": organism,
        "taxonomy": [str(item) for item in taxonomy if str(item).strip()],
        "molecule": str(record.annotations.get("molecule_type") or ""),
        "topology": str(record.annotations.get("topology") or ""),
        "length": len(sequence) if available else reported_length,
        "sequence": sequence,
        "features": features,
        "references": references,
        "record_kind": "wgs_master" if is_wgs or not available else "sequence",
        "sequence_available": available,
        "warnings": warnings,
        "related_records": [],
        "gene": None,
        "genbank_date": str(record.annotations.get("date") or ""),
        "comment": str(record.annotations.get("comment") or ""),
    }


def _pubmed_authors(summary: Dict[str, object]) -> List[str]:
    authors = summary.get("AuthorList") or summary.get("Authors") or []
    if isinstance(authors, str):
        return [authors] if authors.strip() else []
    if not isinstance(authors, list):
        return []
    names: List[str] = []
    for item in authors:
        if isinstance(item, dict):
            name = str(item.get("Name") or item.get("CollectiveName") or "").strip()
        else:
            name = str(item).strip()
        if name:
            names.append(name)
    return names


def _fetch_pubmed_record(pmid: str) -> Dict[str, object]:
    """Busca metadados PubMed via esummary e abstract via efetch (interno).

    Args:
        pmid: PMID numerico ja validado pelo resolve_database.

    Returns:
        Registro HelixScope com record_kind pubmed. sequence_available e False.

    Raises:
        NCBIQueryError: not_found, source_unavailable ou parse falho.
    """
    from Bio import Entrez

    def _esummary() -> object:
        handle = Entrez.esummary(db="pubmed", id=pmid, retmode="xml")
        try:
            return Entrez.read(handle)
        finally:
            handle.close()

    raw = _entrez_retry(_esummary, pmid, "pubmed")
    summaries = list(_iter_document_summaries(raw))
    if not summaries:
        raise NCBIQueryError(
            f"PMID {pmid} nao existe no banco PubMed, ou o NCBI nao devolveu um resumo.",
            "not_found",
        )
    summary = summaries[0] if isinstance(summaries[0], dict) else {}
    title = _summary_field(summary, "Title")
    journal = _summary_field(summary, "FullJournalName", "Source")
    pubdate = _summary_field(summary, "PubDate", "EPubDate", "SortPubDate")
    doi = _summary_field(summary, "DOI")
    volume = _summary_field(summary, "Volume")
    issue = _summary_field(summary, "Issue")
    pages = _summary_field(summary, "Pages")
    authors = _pubmed_authors(summary)
    has_abstract_raw = str(summary.get("HasAbstract") or "").strip()
    abstract_available = has_abstract_raw in {"1", "Y", "Yes", "true", "True"}
    year = ""
    if pubdate not in {"", "Not available"}:
        year = pubdate[:4] if pubdate[:4].isdigit() else ""
    abstract = ""
    abstract_status = "UNAVAILABLE"
    if abstract_available:
        def _efetch_abstract() -> object:
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                rettype="abstract",
                retmode="text",
            )
            try:
                return handle.read()
            finally:
                handle.close()

        try:
            text = str(_entrez_retry(_efetch_abstract, pmid, "pubmed") or "").strip()
            if text:
                abstract = text
                abstract_status = "RETRIEVED"
        except NCBIQueryError:
            abstract_status = "UNAVAILABLE"
    return {
        "accession": pmid,
        "pmid": pmid,
        "description": title if title != "Not available" else f"PubMed {pmid}",
        "title": title,
        "organism": "Not available",
        "taxonomy": [],
        "molecule": "",
        "topology": "",
        "length": 0,
        "sequence": "",
        "features": [],
        "references": [],
        "record_kind": "pubmed",
        "sequence_available": False,
        "authors": authors,
        "journal": journal,
        "pubdate": pubdate,
        "publication_year": year,
        "volume": volume if volume != "Not available" else "",
        "issue": issue if issue != "Not available" else "",
        "pages": pages if pages != "Not available" else "",
        "doi": doi if doi != "Not available" else "",
        "abstract": abstract,
        "abstract_available": bool(abstract),
        "abstract_status": abstract_status,
        "has_abstract_flag": abstract_available,
        "warnings": [
            "PubMed records do not contain a nucleotide or protein sequence. "
            "HelixScope retrieved bibliographic metadata, not a full paper."
        ],
        "comment": "",
    }


def _fetch_gene_record(gene_id: str) -> Dict[str, object]:
    """Busca metadados de um GeneID via esummary e elink (uso interno)."""
    from Bio import Entrez

    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        try:
            payload = Entrez.read(handle)
        finally:
            handle.close()
    except urllib.error.HTTPError as exc:
        raise _http_error(gene_id, "gene", exc) from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise NCBIQueryError(
            f"NCBI indisponivel ao buscar GeneID {gene_id}: {exc}.",
            "source_unavailable",
        ) from exc
    except RuntimeError as exc:
        raise NCBIQueryError(
            f"GeneID {gene_id} nao existe no banco Gene, ou o NCBI nao devolveu "
            f"um resumo: {exc}.",
            "not_found",
        ) from exc

    summary = _first_gene_summary(payload)
    if summary is None:
        raise NCBIQueryError(
            f"GeneID {gene_id} nao existe no banco Gene.",
            "not_found",
        )

    name = str(summary.get("NomenclatureName") or summary.get("Description") or "")
    symbol = str(summary.get("NomenclatureSymbol") or summary.get("Name") or "")
    organism_block = summary.get("Organism") or {}
    if not isinstance(organism_block, dict):
        organism_block = {}
    organism = str(organism_block.get("ScientificName") or "")
    tax_id = str(organism_block.get("TaxID") or "")
    aliases = [
        item.strip()
        for item in str(summary.get("OtherAliases") or "").split(",")
        if item.strip()
    ]
    omim = summary.get("Mim") or []
    if not isinstance(omim, list):
        omim = [omim] if omim else []

    genomic_info = summary.get("GenomicInfo") or []
    related: List[dict] = []
    if isinstance(genomic_info, list):
        for item in genomic_info:
            if not isinstance(item, dict):
                continue
            chrom_acc = str(item.get("ChrAccVer") or "")
            if chrom_acc:
                related.append(
                    {
                        "database": "nucleotide",
                        "id": chrom_acc,
                        "role": "genomic location",
                        "detail": (
                            f"chr {item.get('ChrLoc', '')} "
                            f"{item.get('ChrStart', '')}-{item.get('ChrStop', '')}"
                        ).strip(),
                    }
                )

    related.extend(_gene_related_records(gene_id))

    description = name or str(summary.get("Description") or symbol or f"Gene {gene_id}")
    gene_meta = {
        "id": gene_id,
        "symbol": symbol,
        "name": name or description,
        "aliases": aliases,
        "chromosome": str(summary.get("Chromosome") or ""),
        "map_location": str(summary.get("MapLocation") or ""),
        "summary": str(summary.get("Summary") or ""),
        "tax_id": tax_id,
        "nomenclature_status": str(summary.get("NomenclatureStatus") or ""),
        "omim": [str(item) for item in omim if str(item)],
        "genetic_source": str(summary.get("GeneticSource") or ""),
    }
    return {
        "accession": gene_id,
        "description": description,
        "organism": organism,
        "taxonomy": [],
        "molecule": "",
        "topology": "",
        "length": 0,
        "sequence": "",
        "features": [],
        "references": [],
        "record_kind": "gene",
        "sequence_available": False,
        "warnings": [
            "Registros do banco Gene nao contem a sequencia; use um accession "
            "nucleotide ou protein relacionado para obter bases ou residuos."
        ],
        "related_records": related,
        "gene": gene_meta,
        "genbank_date": "",
        "comment": "",
    }


def _gene_related_records(gene_id: str) -> List[dict]:
    """Obtem nucleotide, protein e PubMed ligados ao GeneID (uso interno)."""
    from Bio import Entrez

    related: List[dict] = []
    queries = (
        ("nucleotide", "gene_nuccore_refseqrna", "RefSeq RNA"),
        ("protein", "gene_protein_refseq", "RefSeq protein"),
        ("pubmed", "gene_pubmed", "PubMed"),
    )
    for database, link_name, role in queries:
        try:
            time.sleep(SLEEP_WITHOUT_KEY)
            handle = Entrez.elink(dbfrom="gene", db=database, id=gene_id)
            try:
                payload = Entrez.read(handle)
            finally:
                handle.close()
        except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError, OSError):
            continue
        ids = _elink_ids(payload, link_name)[:RELATED_LINK_LIMIT]
        accession_map: Dict[str, str] = {}
        if database in {"nucleotide", "protein"} and ids:
            accession_map = _accessions_from_gis(database, ids)
        for identifier in ids:
            related.append(
                {
                    "database": database,
                    "id": accession_map.get(identifier, identifier),
                    "role": role,
                    "detail": "",
                }
            )
    return related


def _elink_ids(payload: object, preferred_link: str) -> List[str]:
    """Extrai IDs de um payload elink, preferindo um LinkName (uso interno)."""
    if not payload:
        return []
    block = payload[0] if isinstance(payload, list) else payload
    if not isinstance(block, dict):
        return []
    link_sets = block.get("LinkSetDb") or []
    chosen = None
    for item in link_sets:
        if isinstance(item, dict) and item.get("LinkName") == preferred_link:
            chosen = item
            break
    if chosen is None and link_sets:
        first = link_sets[0]
        chosen = first if isinstance(first, dict) else None
    if not chosen:
        return []
    ids: List[str] = []
    for link in chosen.get("Link") or []:
        if isinstance(link, dict) and link.get("Id"):
            ids.append(str(link["Id"]))
    return ids


def _accessions_from_gis(db: str, gis: List[str]) -> Dict[str, str]:
    """Converte uma lista de GIs em accession.version via um esummary (interno)."""
    mapping: Dict[str, str] = {gi: gi for gi in gis}
    if not gis:
        return mapping
    from Bio import Entrez

    try:
        time.sleep(SLEEP_WITHOUT_KEY)
        handle = Entrez.esummary(db=db, id=",".join(gis), retmode="xml")
        try:
            payload = Entrez.read(handle)
        finally:
            handle.close()
    except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError, OSError, ValueError):
        return mapping

    documents: List[dict] = []
    if isinstance(payload, list):
        documents = [item for item in payload if isinstance(item, dict)]
    else:
        index = 0
        while True:
            try:
                item = payload[index]  # type: ignore[index]
            except (IndexError, KeyError, TypeError):
                break
            if isinstance(item, dict):
                documents.append(item)
            index += 1
            if index > len(gis) + 2:
                break

    for document in documents:
        uid = str(document.get("Id") or "")
        if not uid and hasattr(document, "attributes"):
            uid = str(document.attributes.get("uid") or "")
        accession = ""
        for key in ("AccessionVersion", "Caption"):
            value = document.get(key)
            if value:
                accession = str(value).strip().split()[0]
                break
        if uid and accession:
            mapping[uid] = accession
    return mapping


def _accession_from_gi(db: str, gi: str) -> Optional[str]:
    """Converte um GI em accession.version via esummary, se possivel (interno)."""
    mapped = _accessions_from_gis(db, [gi]).get(gi)
    if mapped and mapped != gi:
        return mapped
    if mapped == gi:
        return None
    return mapped


def _first_gene_summary(payload: object) -> Optional[dict]:
    """Devolve o primeiro DocumentSummary de um esummary (uso interno)."""
    if payload is None:
        return None
    try:
        first = payload[0]  # type: ignore[index]
        if isinstance(first, dict) and (
            "Name" in first
            or "NomenclatureSymbol" in first
            or "Description" in first
            or "Caption" in first
            or "AccessionVersion" in first
        ):
            return first
    except (KeyError, IndexError, TypeError):
        pass
    if isinstance(payload, dict):
        if "Name" in payload or "NomenclatureSymbol" in payload or "Caption" in payload:
            return payload
        nested = payload.get("DocumentSummarySet") or {}
        if isinstance(nested, dict):
            docs = nested.get("DocumentSummary") or []
            if isinstance(docs, list) and docs and isinstance(docs[0], dict):
                return docs[0]
            if isinstance(docs, dict):
                return docs
    return None


def _sequence_payload(seq_obj: object, undefined_error: type) -> Tuple[str, bool]:
    """Extrai o texto da sequencia ou declara que ela e indefinida (interno)."""
    if getattr(seq_obj, "defined", True) is False:
        return "", False
    try:
        text = str(seq_obj)
    except undefined_error:
        return "", False
    except (TypeError, ValueError, AttributeError):
        return "", False
    if not text or set(text) == {"?"}:
        return "", False
    return text, True


def _undefined_sequence_error_type() -> type:
    """Retorna UndefinedSequenceError do Biopython, ou Exception (interno)."""
    try:
        from Bio.Seq import UndefinedSequenceError

        return UndefinedSequenceError
    except ImportError:
        return Exception


def _entrez_retry(operation, identifier: str, db: str):
    """Executa uma chamada Entrez com backoff diante de timeout/429/5xx."""
    delay = 0.5
    last_error: Optional[Exception] = None
    for attempt in range(SEARCH_RETRY_ATTEMPTS):
        try:
            return operation()
        except urllib.error.HTTPError as exc:
            last_error = exc
            code = int(getattr(exc, "code", 0) or 0)
            if code == 429 or code >= 500:
                if attempt < SEARCH_RETRY_ATTEMPTS - 1:
                    time.sleep(delay)
                    delay *= 2
                    continue
            raise _http_error(identifier, db, exc) from exc
        except (TimeoutError, socket.timeout) as exc:
            last_error = exc
            if attempt < SEARCH_RETRY_ATTEMPTS - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise NCBIQueryError(
                f"NCBI request timed out while querying '{identifier}' in {db}.",
                "source_unavailable",
                failure_kind="timeout",
            ) from exc
        except urllib.error.URLError as exc:
            last_error = exc
            reason = str(getattr(exc, "reason", exc)).lower()
            retryable = "timed out" in reason or "timeout" in reason
            if retryable and attempt < SEARCH_RETRY_ATTEMPTS - 1:
                time.sleep(delay)
                delay *= 2
                continue
            if retryable:
                raise NCBIQueryError(
                    f"NCBI request timed out while querying '{identifier}' in {db}.",
                    "source_unavailable",
                    failure_kind="timeout",
                ) from exc
            raise NCBIQueryError(
                f"NCBI service temporarily unavailable while querying "
                f"'{identifier}' in {db}.",
                "source_unavailable",
            ) from exc
        except (RuntimeError, OSError, ValueError) as exc:
            last_error = exc
            message = str(exc).lower()
            if "timeout" in message and attempt < SEARCH_RETRY_ATTEMPTS - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise NCBIQueryError(
                f"NCBI request failed while querying '{identifier}' in {db}: {exc}",
                "source_unavailable",
            ) from exc
    raise NCBIQueryError(
        f"NCBI request failed while querying '{identifier}' in {db}: {last_error}",
        "source_unavailable",
    )


def _iter_document_summaries(payload: object):
    """Itera DocumentSummary reais de um esummary Entrez (uso interno)."""
    if payload is None:
        return
    if isinstance(payload, list):
        for item in payload:
            if isinstance(item, dict):
                yield item
        return
    if not isinstance(payload, dict):
        return
    nested = payload.get("DocumentSummarySet")
    if isinstance(nested, dict):
        docs = nested.get("DocumentSummary")
        if isinstance(docs, list):
            for item in docs:
                if isinstance(item, dict):
                    yield item
            return
        if isinstance(docs, dict):
            yield docs
            return
    if any(
        key in payload
        for key in ("Title", "Caption", "Name", "NomenclatureSymbol", "Id")
    ):
        yield payload


def _summary_field(summary: Dict[str, object], *names: str) -> str:
    """Primeiro campo textual nao vazio de um DocumentSummary."""
    for name in names:
        value = summary.get(name)
        if value is None:
            continue
        if isinstance(value, list):
            text = " ".join(str(item).strip() for item in value if str(item).strip())
        else:
            text = str(value).strip()
        if text:
            return text
    return "Not available"


def _search_hit_from_summary(
    summary: Dict[str, object], database: str, fallback_id: str
) -> Dict[str, object]:
    """Extrai metadados realmente presentes em um esummary (uso interno)."""
    item_id = _summary_field(summary, "Id", "Gi", "UID")
    if item_id == "Not available":
        item_id = fallback_id or "Not available"
    accession = _summary_field(
        summary, "AccessionVersion", "Caption", "Accession", "Name"
    )
    if accession == "Not available":
        accession = item_id
    title = _summary_field(
        summary, "Title", "Description", "NomenclatureName", "Summary"
    )
    organism = _summary_field(summary, "Organism", "ScientificName")
    source = _summary_field(summary, "Source")
    length_value = _summary_length(summary)
    extra_parts = []
    if source != "Not available":
        extra_parts.append(f"source={source}")
    create_date = _summary_field(summary, "CreateDate", "PubDate", "SortDate")
    if create_date != "Not available":
        extra_parts.append(f"date={create_date}")
    chromosome = _summary_field(summary, "Chromosome")
    if chromosome != "Not available":
        extra_parts.append(f"chromosome={chromosome}")
    tax_id = _summary_field(summary, "TaxId")
    if tax_id != "Not available":
        extra_parts.append(f"taxid={tax_id}")
    authors = _pubmed_authors(summary) if database == "pubmed" else []
    journal = ""
    pubdate = ""
    doi = ""
    has_abstract: Optional[bool] = None
    if database == "pubmed":
        journal = _summary_field(summary, "FullJournalName", "Source")
        pubdate = _summary_field(summary, "PubDate", "EPubDate")
        doi = _summary_field(summary, "DOI")
        flag = str(summary.get("HasAbstract") or "").strip()
        has_abstract = flag in {"1", "Y", "Yes", "true", "True"} if flag else None
        if authors:
            extra_parts.append("authors=" + ", ".join(authors[:8]))
        if journal != "Not available":
            extra_parts.append(f"journal={journal}")
    return {
        "id": item_id,
        "accession": accession,
        "title": title,
        "organism": organism,
        "length": length_value if length_value > 0 else None,
        "database": database,
        "extra": "; ".join(extra_parts) if extra_parts else "Not available",
        "authors": authors,
        "journal": journal if journal != "Not available" else "",
        "pubdate": pubdate if pubdate != "Not available" else "",
        "doi": doi if doi != "Not available" else "",
        "has_abstract": has_abstract,
    }


def _first_taxonomy_record(raw: object) -> object:
    """Extrai o primeiro taxon de uma resposta Entrez.read (interno)."""
    if raw is None:
        raise NCBIQueryError("Empty NCBI Taxonomy response.", "internal")
    if isinstance(raw, list):
        if not raw:
            raise NCBIQueryError("NCBI Taxonomy returned an empty taxon list.", "not_found")
        return raw[0]
    if isinstance(raw, dict):
        for key in ("Taxon", "TaxaSet"):
            nested = raw.get(key)
            if nested is None:
                continue
            return _first_taxonomy_record(nested)
        if raw.get("TaxId") or raw.get("ScientificName"):
            return raw
    raise NCBIQueryError("Malformed NCBI Taxonomy XML payload.", "internal")


def _parse_taxonomy_record(record: object) -> Dict[str, object]:
    """Normaliza um taxon Entrez em campos HelixScope (interno)."""
    if not isinstance(record, dict):
        try:
            record = dict(record)
        except (TypeError, ValueError) as exc:
            raise NCBIQueryError(
                "Malformed NCBI Taxonomy record.",
                "internal",
            ) from exc
    tax_id = str(record.get("TaxId") or record.get("taxid") or "").strip()
    if not TAXON_ID_PATTERN.fullmatch(tax_id):
        raise NCBIQueryError("NCBI Taxonomy record has no numeric TaxId.", "internal")
    scientific = str(record.get("ScientificName") or "").strip()
    if not scientific:
        raise NCBIQueryError("NCBI Taxonomy record has no scientific name.", "internal")
    rank = str(record.get("Rank") or "").strip()
    lineage: List[dict] = []
    lineage_ex = record.get("LineageEx") or record.get("Lineage") or []
    if isinstance(lineage_ex, dict):
        lineage_ex = lineage_ex.get("Taxon") or []
    if isinstance(lineage_ex, str):
        for name in [part.strip() for part in lineage_ex.split(";") if part.strip()]:
            lineage.append(
                {"taxon_id": "", "scientific_name": name, "rank": ""}
            )
    elif isinstance(lineage_ex, list):
        for item in lineage_ex:
            if not isinstance(item, dict):
                try:
                    item = dict(item)
                except (TypeError, ValueError):
                    continue
            lineage.append(
                {
                    "taxon_id": str(item.get("TaxId") or "").strip(),
                    "scientific_name": str(item.get("ScientificName") or "").strip(),
                    "rank": str(item.get("Rank") or "").strip(),
                }
            )
    other = record.get("OtherNames") or {}
    if not isinstance(other, dict):
        try:
            other = dict(other)
        except (TypeError, ValueError):
            other = {}
    common = str(
        other.get("GenbankCommonName")
        or other.get("CommonName")
        or record.get("CommonName")
        or ""
    ).strip()
    if isinstance(other.get("CommonName"), list):
        names = [str(item).strip() for item in other.get("CommonName") or [] if str(item).strip()]
        common = names[0] if names else common
    return {
        "taxon_id": tax_id,
        "scientific_name": scientific,
        "rank": rank,
        "lineage": lineage,
        "common_name": common,
    }


def _http_error(identifier: str, db: str, exc: urllib.error.HTTPError) -> NCBIQueryError:
    """Converte HTTPError do NCBI em NCBIQueryError classificado (interno)."""
    code = int(getattr(exc, "code", 0) or 0)
    if code == 429:
        return NCBIQueryError(
            f"NCBI recusou a consulta a '{identifier}' por limite de taxa "
            f"(HTTP 429) no banco {db}. Aguarde e tente de novo.",
            "source_unavailable",
            failure_kind="rate_limited",
        )
    if code in {400, 404}:
        return NCBIQueryError(
            f"Nenhum registro valido para '{identifier}' no banco {db} "
            f"(HTTP {code}).",
            "not_found",
        )
    if code >= 500:
        return NCBIQueryError(
            f"NCBI indisponivel (HTTP {code}) ao buscar '{identifier}' em {db}.",
            "source_unavailable",
        )
    return NCBIQueryError(
        f"Erro HTTP {code} do NCBI ao buscar '{identifier}' em {db}: {exc}.",
        "source_unavailable",
    )
