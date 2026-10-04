"""Dominios proteicos do InterPro e anotacao de referencia do UniProt.

HelixScope nao prediz dominios. As entradas vem do InterPro (Blum et al., Nucleic
Acids Research), que integra Pfam, PROSITE, SMART, CDD, PANTHER e outros, e a
anotacao de proteina vem do UniProtKB. Cada intervalo devolvido preserva a base
de dados de origem: um dominio recuperado via InterPro cuja evidencia e Pfam
continua identificado como Pfam.

Endpoints usados:
    - GET https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{acc}
    - GET https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}
    - GET https://rest.uniprot.org/uniprotkb/{acc}.json

Limite cientifico explicito: estar dentro de um dominio anotado NAO prova efeito
funcional de uma variante. O dominio e contexto estrutural/evolutivo, nao
mecanismo.

PRIVACIDADE DE REDE: estas consultas enviam o accession UniProt para o EMBL-EBI.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import json
import re
import socket
import time
import urllib.parse
from typing import Any, Mapping, Optional
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

from . import provenance

INTERPRO_HOST: str = "www.ebi.ac.uk"
UNIPROT_HOST: str = "rest.uniprot.org"

ALLOWED_URLS: tuple[tuple[str, str], ...] = (
    (INTERPRO_HOST, "/interpro/api/"),
    (UNIPROT_HOST, "/uniprotkb/"),
)
"""Pares host/prefixo permitidos. Qualquer outro destino e recusado."""

USER_AGENT: str = (
    f"HelixScope/{provenance.HELIXSCOPE_VERSION} (bioinformatics; InterPro and UniProt client)"
)

HTTP_TIMEOUT_S: float = 60.0
MAX_RESPONSE_BYTES: int = 8 * 1024 * 1024
MAX_RETRIES: int = 3
BACKOFF_BASE_S: float = 1.5
MAX_BACKOFF_S: float = 20.0
MAX_PAGES: int = 10

UNIPROT_ACCESSION_PATTERN: re.Pattern[str] = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:-\d{1,3})?$"
)
"""Padrao oficial de accession UniProtKB, com sufixo de isoforma opcional."""

STRUCTURAL_ENTRY_TYPES: frozenset[str] = frozenset(
    {"domain", "family", "homologous_superfamily", "repeat", "conserved_site",
     "active_site", "binding_site", "ptm"}
)
"""Tipos de entrada InterPro considerados anotacao posicional utilizavel."""


class DomainError(RuntimeError):
    """Falha ao consultar InterPro ou UniProt.

    Attributes:
        category: NOT_FOUND, TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE,
            PARSING_ERROR, NETWORK_ERROR, INVALID_INPUT ou RESOURCE_LIMIT.
        retry_after_s: Segundos sugeridos pelo servico, quando informado.
    """

    def __init__(self, message: str, category: str, retry_after_s: Optional[float] = None) -> None:
        super().__init__(message)
        self.category = str(category or "SERVICE_UNAVAILABLE").strip().upper()
        self.retry_after_s = retry_after_s


class _AllowlistRedirect(HTTPRedirectHandler):
    """Recusa redirects para fora da allowlist InterPro/UniProt."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        if not url_is_allowed(newurl):
            raise DomainError(
                f"Refusing a redirect to a non-allowlisted URL: {newurl}.",
                "INVALID_INPUT",
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def url_is_allowed(url: str) -> bool:
    """Valida esquema, host e prefixo de caminho contra a allowlist.

    Args:
        url: URL absoluta.

    Returns:
        True apenas para HTTPS nos endpoints InterPro/UniProt documentados.

    Raises:
        Nenhum.
    """
    parsed = urllib.parse.urlparse(str(url or ""))
    if parsed.scheme != "https":
        return False
    host = (parsed.hostname or "").lower()
    path = parsed.path or ""
    return any(host == allowed and path.startswith(prefix) for allowed, prefix in ALLOWED_URLS)


def network_disclosure() -> dict:
    """Descreve o que sai da maquina nas consultas de dominio.

    Args:
        Nenhum.

    Returns:
        Dict com hosts, dados enviados e limites de servico.

    Raises:
        Nenhum.
    """
    return {
        "remote": True,
        "services": [
            {
                "name": "InterPro API (EMBL-EBI)",
                "host": INTERPRO_HOST,
                "reference": "https://interpro-documentation.readthedocs.io/",
            },
            {
                "name": "UniProtKB REST",
                "host": UNIPROT_HOST,
                "reference": "https://www.uniprot.org/help/api",
            },
        ],
        "transport": "HTTPS",
        "data_sent": ["the UniProt accession you query"],
        "data_not_sent": [
            "variant coordinates",
            "any sequence loaded in HelixScope",
            "local file names or paths",
        ],
        "user_agent": USER_AGENT,
        "rate_limit_policy": (
            "InterPro returns HTTP 408 for queries longer than a minute and moves "
            "them to background; HelixScope treats 408 as TIMEOUT and does not "
            "hammer the endpoint. Pagination is bounded to "
            f"{MAX_PAGES} pages."
        ),
        "scientific_limit": (
            "A residue falling inside an annotated domain is structural and "
            "evolutionary context. It is not evidence of a functional effect."
        ),
    }


def validate_uniprot_accession(accession: str) -> str:
    """Valida um accession UniProtKB contra o padrao oficial.

    Args:
        accession: Accession candidato, com ou sem sufixo de isoforma.

    Returns:
        Accession em maiusculas.

    Raises:
        DomainError: INVALID_INPUT se nao corresponder ao padrao UniProtKB.

    Nota biologica:
        VEP devolve accessions com versao (`P38398.282`) e isoformas
        (`P38398-1`). A versao e removida porque os endpoints aceitam o
        accession base; a isoforma e preservada porque identifica outra
        sequencia proteica.
    """
    text = str(accession or "").strip().upper()
    if "." in text:
        text = text.split(".", 1)[0]
    if not text:
        raise DomainError("UniProt accession is empty.", "INVALID_INPUT")
    if len(text) > 16 or not UNIPROT_ACCESSION_PATTERN.match(text):
        raise DomainError(
            f"'{accession}' is not a valid UniProtKB accession.", "INVALID_INPUT"
        )
    return text


def parse_interpro_entries(payload: Any, *, accession: str) -> list[dict]:
    """Extrai entradas e intervalos de uma pagina da API InterPro.

    Args:
        payload: JSON desserializado de `/entry/.../protein/uniprot/{acc}`.
        accession: Accession consultado, para proveniencia.

    Returns:
        Lista de dicts com entry_accession, name, type, source_database,
        member_databases, ranges e protein_length.

    Raises:
        DomainError: PARSING_ERROR se a estrutura nao tiver `results`.

    Nota biologica:
        Um dominio pode ser descontinuo: a API devolve `fragments` e cada
        fragmento e um intervalo 1-based inclusive na sequencia da proteina.
        Todos os fragmentos sao preservados; unir fragmentos num unico intervalo
        falsearia a arquitetura do dominio.
    """
    if not isinstance(payload, Mapping):
        raise DomainError("InterPro response is not a JSON object.", "PARSING_ERROR")
    results = payload.get("results")
    if results is None:
        raise DomainError("InterPro response has no results block.", "PARSING_ERROR")
    if not isinstance(results, list):
        raise DomainError("InterPro results block is not a list.", "PARSING_ERROR")
    rows: list[dict] = []
    for item in results:
        if not isinstance(item, Mapping):
            continue
        metadata = item.get("metadata")
        if not isinstance(metadata, Mapping):
            continue
        members: list[dict] = []
        raw_members = metadata.get("member_databases")
        if isinstance(raw_members, Mapping):
            for database, entries in raw_members.items():
                if isinstance(entries, Mapping):
                    for member_accession, member_name in entries.items():
                        members.append(
                            {
                                "source_database": str(database),
                                "accession": str(member_accession),
                                "name": str(member_name),
                            }
                        )
        protein_length: Optional[int] = None
        ranges: list[dict] = []
        for protein in item.get("proteins") or []:
            if not isinstance(protein, Mapping):
                continue
            length = _optional_int(protein.get("protein_length"))
            if length is not None:
                protein_length = length
            for location in protein.get("entry_protein_locations") or []:
                if not isinstance(location, Mapping):
                    continue
                fragments = []
                for fragment in location.get("fragments") or []:
                    if not isinstance(fragment, Mapping):
                        continue
                    start = _optional_int(fragment.get("start"))
                    end = _optional_int(fragment.get("end"))
                    if start is None or end is None:
                        continue
                    fragments.append(
                        {
                            "start_1based": start,
                            "end_1based": end,
                            "continuity": str(fragment.get("dc-status") or ""),
                        }
                    )
                if fragments:
                    ranges.append(
                        {
                            "fragments": fragments,
                            "representative": bool(location.get("representative")),
                            "model": str(location.get("model") or ""),
                            "score": _optional_float(location.get("score")),
                        }
                    )
        go_terms = []
        for term in metadata.get("go_terms") or []:
            if isinstance(term, Mapping):
                go_terms.append(
                    {
                        "identifier": str(term.get("identifier") or ""),
                        "name": str(term.get("name") or ""),
                        "category": str(
                            (term.get("category") or {}).get("code")
                            if isinstance(term.get("category"), Mapping)
                            else term.get("category") or ""
                        ),
                    }
                )
        rows.append(
            {
                "entry_accession": str(metadata.get("accession") or "").upper(),
                "name": str(metadata.get("name") or ""),
                "type": str(metadata.get("type") or ""),
                "source_database": str(metadata.get("source_database") or ""),
                "integrated": str(metadata.get("integrated") or ""),
                "member_databases": members,
                "go_terms": go_terms,
                "ranges": ranges,
                "protein_length": protein_length,
                "uniprot_accession": str(accession or "").upper(),
                "evidence_status": "RETRIEVED",
                "source": "InterPro (EMBL-EBI)",
                "url": (
                    "https://www.ebi.ac.uk/interpro/entry/"
                    f"{str(metadata.get('source_database') or 'InterPro')}/"
                    f"{str(metadata.get('accession') or '')}/"
                ),
            }
        )
    return rows


def fetch_protein_domains(
    *,
    accession: str,
    urlopen_fn=None,
    include_member_databases: bool = False,
) -> dict:
    """Recupera as entradas InterPro que anotam uma proteina. FAZ REDE.

    Args:
        accession: Accession UniProtKB.
        urlopen_fn: Injecao para teste.
        include_member_databases: True consulta `/entry/all/` (InterPro mais
            entradas de bases membro); False consulta apenas entradas InterPro
            integradas.

    Returns:
        Dict status (AVAILABLE ou NOT_FOUND), entries, count, accession,
        proveniencia e limite cientifico.

    Raises:
        DomainError: INVALID_INPUT, TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE,
            PARSING_ERROR ou NETWORK_ERROR.

    Nota biologica:
        HTTP 204 do InterPro significa "sem conteudo para esta combinacao", ou
        seja proteina sem entradas, e nao falha de servico.
    """
    uniprot = validate_uniprot_accession(accession)
    scope = "all" if include_member_databases else "interpro"
    url = (
        f"https://{INTERPRO_HOST}/interpro/api/entry/{scope}/protein/uniprot/"
        f"{urllib.parse.quote(uniprot, safe='')}?page_size=100"
    )
    entries: list[dict] = []
    pages = 0
    total: Optional[int] = None
    while url and pages < MAX_PAGES:
        pages += 1
        try:
            payload = _http_json(url, urlopen_fn=urlopen_fn)
        except DomainError as exc:
            if exc.category == "NOT_FOUND" and pages == 1:
                return {
                    "status": "NOT_FOUND",
                    "source": "InterPro (EMBL-EBI)",
                    "evidence_status": "RETRIEVED",
                    "uniprot_accession": uniprot,
                    "entries": [],
                    "count": 0,
                    "message": (
                        f"InterPro has no entry annotating {uniprot}. Absence of "
                        "an annotated domain is not evidence that the region is "
                        "unstructured or unimportant."
                    ),
                    "scientific_limit": network_disclosure()["scientific_limit"],
                    "retrieved_at_utc": provenance.utc_now(),
                    "software_version": provenance.HELIXSCOPE_VERSION,
                }
            raise
        if not payload:
            break
        if total is None:
            total = _optional_int(payload.get("count") if isinstance(payload, Mapping) else None)
        entries.extend(parse_interpro_entries(payload, accession=uniprot))
        next_url = payload.get("next") if isinstance(payload, Mapping) else None
        url = str(next_url) if next_url else ""
        if url and not url_is_allowed(url):
            raise DomainError(
                "InterPro pagination pointed outside the allowlist.", "INVALID_INPUT"
            )
        if url:
            time.sleep(0.35)
    positional = [row for row in entries if row.get("ranges")]
    return {
        "status": "AVAILABLE" if entries else "NOT_FOUND",
        "source": "InterPro (EMBL-EBI)",
        "evidence_status": "RETRIEVED",
        "uniprot_accession": uniprot,
        "entries": entries,
        "count": len(entries),
        "positional_count": len(positional),
        "reported_total": total,
        "pages_fetched": pages,
        "truncated": bool(url),
        "scope": scope,
        "message": (
            f"{len(entries)} InterPro entry/entries retrieved for {uniprot}; "
            f"{len(positional)} carry explicit residue ranges."
            if entries
            else (
                f"InterPro returned no entry for {uniprot}. Absence of an "
                "annotated domain is not evidence of absence of structure."
            )
        ),
        "scientific_limit": network_disclosure()["scientific_limit"],
        "method": (
            "InterPro integrated signatures retrieved through the official "
            "InterPro API. HelixScope does not run HMM searches locally and does "
            "not predict domains."
        ),
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def domains_covering_residue(
    entries: Any,
    *,
    residue_1based: int,
) -> dict:
    """Seleciona as entradas cujos intervalos contem um residuo.

    Args:
        entries: Lista de entradas de `fetch_protein_domains` ou o proprio dict.
        residue_1based: Posicao do residuo na proteina, 1-based.

    Returns:
        Dict residue_1based, covering (entradas com o fragmento que casa),
        covering_count e note.

    Raises:
        DomainError: INVALID_INPUT se o residuo nao for um inteiro positivo.

    Nota biologica:
        Cobertura e uma afirmacao posicional. Um residuo dentro de um dominio
        BRCT nao passa a ser patogenico por isso; apenas fica descrito no
        contexto daquele dominio.
    """
    try:
        residue = int(residue_1based)
    except (TypeError, ValueError) as exc:
        raise DomainError("residue_1based must be an integer.", "INVALID_INPUT") from exc
    if residue < 1:
        raise DomainError("residue_1based must be >= 1.", "INVALID_INPUT")
    rows = entries.get("entries") if isinstance(entries, Mapping) else entries
    if rows is None:
        rows = []
    covering: list[dict] = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        hits = []
        for interval in row.get("ranges") or []:
            for fragment in (interval or {}).get("fragments") or []:
                start = fragment.get("start_1based")
                end = fragment.get("end_1based")
                if start is None or end is None:
                    continue
                if int(start) <= residue <= int(end):
                    hits.append(
                        {
                            "start_1based": int(start),
                            "end_1based": int(end),
                            "continuity": str(fragment.get("continuity") or ""),
                            "representative": bool(interval.get("representative")),
                        }
                    )
        if hits:
            covering.append(
                {
                    "entry_accession": str(row.get("entry_accession") or ""),
                    "name": str(row.get("name") or ""),
                    "type": str(row.get("type") or ""),
                    "source_database": str(row.get("source_database") or ""),
                    "member_databases": list(row.get("member_databases") or []),
                    "matching_fragments": hits,
                    "url": str(row.get("url") or ""),
                    "evidence_status": "RETRIEVED",
                    "source": str(row.get("source") or "InterPro (EMBL-EBI)"),
                }
            )
    return {
        "residue_1based": residue,
        "covering": covering,
        "covering_count": len(covering),
        "note": (
            f"{len(covering)} annotated InterPro entry/entries span residue "
            f"{residue}. Overlap is positional context retrieved from InterPro; "
            "it does not demonstrate a functional effect at this residue."
            if covering
            else (
                f"No annotated InterPro entry spans residue {residue}. This is "
                "absence of annotation, not evidence that the residue is "
                "unimportant."
            )
        ),
        "retrieved_at_utc": provenance.utc_now(),
    }


def fetch_uniprot_protein(*, accession: str, urlopen_fn=None) -> dict:
    """Recupera anotacao de referencia de uma proteina no UniProtKB. FAZ REDE.

    Args:
        accession: Accession UniProtKB.
        urlopen_fn: Injecao para teste.

    Returns:
        Dict com accession, entry_name, protein_name, gene_names, organism,
        sequence_length, sequence_checksum, existence, reviewed e xrefs de PDB e
        AlphaFold quando presentes.

    Raises:
        DomainError: INVALID_INPUT, NOT_FOUND e restantes categorias de rede.

    Nota biologica:
        O comprimento e o checksum CRC64 da sequencia permitem verificar que a
        posicao de um residuo faz sentido nesta isoforma antes de a mapear numa
        estrutura.
    """
    uniprot = validate_uniprot_accession(accession)
    url = f"https://{UNIPROT_HOST}/uniprotkb/{urllib.parse.quote(uniprot, safe='')}.json"
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    if not isinstance(payload, Mapping):
        raise DomainError("UniProt response is not a JSON object.", "PARSING_ERROR")
    description = payload.get("proteinDescription")
    protein_name = ""
    if isinstance(description, Mapping):
        recommended = description.get("recommendedName")
        if isinstance(recommended, Mapping):
            full = recommended.get("fullName")
            if isinstance(full, Mapping):
                protein_name = str(full.get("value") or "")
    genes = []
    for gene in payload.get("genes") or []:
        if isinstance(gene, Mapping):
            name = gene.get("geneName")
            if isinstance(name, Mapping) and name.get("value"):
                genes.append(str(name["value"]))
    organism = payload.get("organism")
    organism_name = ""
    taxon_id: Optional[int] = None
    if isinstance(organism, Mapping):
        organism_name = str(organism.get("scientificName") or "")
        taxon_id = _optional_int(organism.get("taxonId"))
    sequence = payload.get("sequence")
    sequence_length: Optional[int] = None
    checksum = ""
    sequence_residues = ""
    mol_weight_da: Optional[float] = None
    if isinstance(sequence, Mapping):
        sequence_length = _optional_int(sequence.get("length"))
        checksum = str(sequence.get("crc64") or "")
        sequence_residues = str(sequence.get("value") or "")
        mol_weight_da = _optional_float(sequence.get("molWeight"))
    pdb_ids: list[str] = []
    alphafold_ids: list[str] = []
    cross_refs: list[dict[str, str]] = []
    for xref in payload.get("uniProtKBCrossReferences") or []:
        if not isinstance(xref, Mapping):
            continue
        database = str(xref.get("database") or "").upper()
        identifier = str(xref.get("id") or "")
        if not identifier:
            continue
        if database == "PDB":
            pdb_ids.append(identifier)
        elif database == "ALPHAFOLDDB":
            alphafold_ids.append(identifier)
        if len(cross_refs) < 80:
            cross_refs.append({"database": database, "id": identifier})
    comments, function_comments = _uniprot_comments(payload)
    features = _uniprot_features(payload)
    literature = _uniprot_literature(payload)
    entry_audit = _uniprot_entry_audit(payload)
    return {
        "status": "AVAILABLE",
        "source": "UniProtKB",
        "evidence_status": "RETRIEVED",
        "accession": str(payload.get("primaryAccession") or uniprot).upper(),
        "entry_name": str(payload.get("uniProtkbId") or ""),
        "protein_name": protein_name,
        "gene_names": genes,
        "organism": organism_name,
        "taxon_id": taxon_id,
        "sequence_length": sequence_length,
        "sequence_checksum_crc64": checksum,
        "protein_existence": str(payload.get("proteinExistence") or ""),
        "reviewed": str(payload.get("entryType") or "").lower().startswith("uniprotkb reviewed"),
        "entry_type": str(payload.get("entryType") or ""),
        "pdb_cross_references": pdb_ids,
        "alphafold_cross_references": alphafold_ids,
        "url": f"https://www.uniprot.org/uniprotkb/{uniprot}/entry",
        "method": (
            "UniProtKB entry retrieved through the official REST API. HelixScope "
            "does not edit or infer UniProt annotation."
        ),
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "retrieved_sequence": sequence_residues,
        "retrieved_sequence_status": (
            "RETRIEVED" if sequence_residues else "UNAVAILABLE"
        ),
        "molecular_weight_da": mol_weight_da,
        "function_comments": function_comments,
        "comments": comments,
        "features": features,
        "literature": literature,
        "cross_references": cross_refs,
        "entry_audit": entry_audit,
        "identity_scope": "database_retrieved",
        "identity_note": (
            "These fields are UniProtKB annotation for the requested accession. "
            "They were not computed from an anonymous pasted sequence."
        ),
    }


def _http_json(url: str, *, urlopen_fn=None) -> Any:
    if not url_is_allowed(url):
        raise DomainError(
            "Refusing a non-allowlisted InterPro/UniProt URL.", "INVALID_INPUT"
        )
    request = Request(
        url,
        method="GET",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    opener = urlopen_fn or build_opener(_AllowlistRedirect()).open
    attempt = 0
    while True:
        attempt += 1
        try:
            with opener(request, timeout=HTTP_TIMEOUT_S) as handle:
                raw = handle.read()
                status = int(getattr(handle, "status", 0) or getattr(handle, "code", 0) or 0)
                final = str(getattr(handle, "geturl", lambda: url)() or url)
            break
        except HTTPError as exc:
            error = _http_error(exc)
            if error.category in {"RATE_LIMITED", "SERVICE_UNAVAILABLE"} and attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, error.retry_after_s))
                continue
            raise error from exc
        except DomainError:
            raise
        except (URLError, TimeoutError, socket.timeout, OSError) as exc:
            message = str(exc).lower()
            timed_out = isinstance(exc, (TimeoutError, socket.timeout)) or "timed out" in message
            error = DomainError(
                (
                    "InterPro/UniProt request timed out. This is TIMEOUT, not an "
                    "absence of domains."
                )
                if timed_out
                else (
                    f"InterPro/UniProt network error: {exc}. This is "
                    "NETWORK_ERROR, not NOT_FOUND."
                ),
                "TIMEOUT" if timed_out else "NETWORK_ERROR",
            )
            if attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, None))
                continue
            raise error from exc
    if final and not url_is_allowed(final):
        raise DomainError(
            "InterPro/UniProt answered from a non-allowlisted URL.", "INVALID_INPUT"
        )
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    if len(raw) > MAX_RESPONSE_BYTES:
        raise DomainError(
            f"Response exceeds {MAX_RESPONSE_BYTES:,} bytes.", "RESOURCE_LIMIT"
        )
    text = raw.decode("utf-8", errors="replace")
    if status == 204 or not text.strip():
        return {}
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise DomainError(
            "InterPro/UniProt response is not valid JSON.", "PARSING_ERROR"
        ) from exc


def _backoff_seconds(attempt: int, retry_after_s: Optional[float]) -> float:
    if retry_after_s is not None and retry_after_s >= 0:
        return min(float(retry_after_s), MAX_BACKOFF_S)
    return min(BACKOFF_BASE_S ** attempt, MAX_BACKOFF_S)


def _http_error(exc: HTTPError) -> DomainError:
    code = int(getattr(exc, "code", 0) or 0)
    retry_after: Optional[float] = None
    headers = getattr(exc, "headers", None)
    if headers is not None:
        raw_retry = headers.get("Retry-After")
        if raw_retry:
            try:
                retry_after = float(str(raw_retry).strip())
            except ValueError:
                retry_after = None
    if code == 404:
        return DomainError("InterPro/UniProt has no record for this accession.", "NOT_FOUND")
    if code == 410:
        return DomainError(
            "This InterPro entry was retired (HTTP 410).", "NOT_FOUND"
        )
    if code == 408:
        return DomainError(
            "InterPro moved the query to background processing (HTTP 408). This "
            "is TIMEOUT; retry later instead of hammering the endpoint.",
            "TIMEOUT",
            retry_after,
        )
    if code == 429:
        return DomainError(
            "InterPro/UniProt rate limited the request (HTTP 429).",
            "RATE_LIMITED",
            retry_after,
        )
    if code == 400:
        return DomainError(
            f"InterPro/UniProt rejected the request (HTTP 400): {exc}.", "INVALID_INPUT"
        )
    if code >= 500:
        return DomainError(
            f"InterPro/UniProt service unavailable (HTTP {code}).",
            "SERVICE_UNAVAILABLE",
            retry_after,
        )
    return DomainError(f"InterPro/UniProt HTTP {code}: {exc}.", "SERVICE_UNAVAILABLE", retry_after)


def _optional_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _uniprot_comments(payload: Mapping[str, Any]) -> tuple[list[dict[str, str]], list[str]]:
    comments: list[dict[str, str]] = []
    function_comments: list[str] = []
    for item in payload.get("comments") or []:
        if not isinstance(item, Mapping):
            continue
        comment_type = str(item.get("commentType") or "")
        texts: list[str] = []
        for text in item.get("texts") or []:
            if isinstance(text, Mapping) and text.get("value"):
                texts.append(str(text["value"]))
        joined = " ".join(texts).strip()
        if not comment_type and not joined:
            continue
        comments.append({"comment_type": comment_type, "text": joined})
        if comment_type.upper() == "FUNCTION" and joined:
            function_comments.append(joined)
    return comments, function_comments


def _uniprot_features(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    features: list[dict[str, Any]] = []
    for item in payload.get("features") or []:
        if not isinstance(item, Mapping):
            continue
        location = item.get("location") if isinstance(item.get("location"), Mapping) else {}
        start = location.get("start") if isinstance(location.get("start"), Mapping) else {}
        end = location.get("end") if isinstance(location.get("end"), Mapping) else {}
        features.append(
            {
                "type": str(item.get("type") or ""),
                "description": str(item.get("description") or ""),
                "start": _optional_int(start.get("value")),
                "end": _optional_int(end.get("value")),
            }
        )
    return features


def _uniprot_literature(payload: Mapping[str, Any]) -> list[dict[str, str]]:
    literature: list[dict[str, str]] = []
    for item in payload.get("references") or []:
        if not isinstance(item, Mapping):
            continue
        citation = item.get("citation") if isinstance(item.get("citation"), Mapping) else {}
        pmid = ""
        for xref in citation.get("citationCrossReferences") or []:
            if isinstance(xref, Mapping) and str(xref.get("database") or "").upper() == "PUBMED":
                pmid = str(xref.get("id") or "")
                break
        literature.append(
            {
                "title": str(citation.get("title") or ""),
                "journal": str(citation.get("journal") or ""),
                "publication_date": str(citation.get("publicationDate") or ""),
                "pmid": pmid,
            }
        )
    return literature


def _uniprot_entry_audit(payload: Mapping[str, Any]) -> dict[str, Any]:
    audit = payload.get("entryAudit")
    if not isinstance(audit, Mapping):
        return {}
    return {
        "first_public_date": str(audit.get("firstPublicDate") or ""),
        "last_annotation_update_date": str(audit.get("lastAnnotationUpdateDate") or ""),
        "entry_version": _optional_int(audit.get("entryVersion")),
        "sequence_version": _optional_int(audit.get("sequenceVersion")),
    }


def _optional_float(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in {float("inf"), float("-inf")}:
        return None
    return number
