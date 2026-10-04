"""Evidencia clinica RECUPERADA do NCBI ClinVar. Nao e um preditor HelixScope.

ClinVar e um arquivo publico de interpretacoes SUBMETIDAS por laboratorios e
consorcios. A classificacao de um registo submetido (SCV) e do submissor; a
classificacao agregada (VCV/RCV) e calculada pelo NCBI a partir dos SCV segundo
o review status. HelixScope nao classifica variantes: apenas recupera, preserva
identificadores e mostra divergencias entre submissoes.

Consequencias diretas desta politica:
    - Ausencia de registo devolve NOT_FOUND, nunca "benign".
    - Submissoes divergentes sao mostradas como conflito, nao reduzidas a uma
      conclusao unica.
    - `retrieved_at_utc` acompanha todo resultado: ClinVar muda continuamente.
    - Nada aqui constitui diagnostico ou recomendacao clinica.

Endpoints E-utilities usados (documentados em
https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/):
    - esearch.fcgi?db=clinvar  resolve termo/identificador para UID
    - esummary.fcgi?db=clinvar&retmode=json  resumo do registo agregado

PRIVACIDADE DE REDE: consultar ClinVar ENVIA o termo de busca (coordenada,
notacao HGVS, rsID ou accession) para o NCBI.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import json
import re
import socket
import time
import urllib.parse
from typing import Any, Mapping, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

from . import ncbi_fetch, provenance

EUTILS_HOST: str = "eutils.ncbi.nlm.nih.gov"
EUTILS_BASE: str = f"https://{EUTILS_HOST}/entrez/eutils"

ALLOWED_HOSTS: frozenset[str] = frozenset({EUTILS_HOST})
ALLOWED_PATH_PREFIX: str = "/entrez/eutils/"

USER_AGENT: str = f"HelixScope/{provenance.HELIXSCOPE_VERSION} (bioinformatics; ClinVar E-utilities client)"

HTTP_TIMEOUT_S: float = 45.0
MAX_RESPONSE_BYTES: int = 4 * 1024 * 1024
MAX_RETRIES: int = 3
BACKOFF_BASE_S: float = 1.5
MAX_BACKOFF_S: float = 15.0
MAX_TERM_LENGTH: int = 300
MAX_UIDS: int = 20

CONFLICT_MARKERS: tuple[str, ...] = (
    "conflicting",
    "conflict",
)
"""Marcadores textuais que o NCBI usa para sinalizar divergencia agregada."""

REVIEW_STATUS_ORDER: tuple[str, ...] = (
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
    "criteria provided, conflicting classifications",
    "criteria provided, single submitter",
    "no assertion criteria provided",
    "no classification provided",
    "no classification for the single variant",
)
"""Ordem de precedencia de review status publicada pelo ClinVar.

Usada apenas para ORDENAR a apresentacao. HelixScope nao recalcula a
classificacao agregada a partir desta ordem.
"""

CONTIG_TOKEN_PATTERN = re.compile(r"^[A-Za-z0-9._-]{1,64}$")
"""Token de cromossoma aceito numa busca por coordenada.

Restringir aqui nao e cosmetico. O termo Entrez `17[chr] AND 43093557[chrpos]`
aceita operadores; um token como `17[chr] OR 13` alargaria a busca a outro
cromossoma e as submissoes recuperadas seriam atribuidas a variante errada.
"""

ACCESSION_PATTERN = re.compile(r"^(VCV|RCV|SCV)\d{6,12}(\.\d+)?$", re.IGNORECASE)
"""Accession ClinVar: VCV/RCV/SCV seguido de digitos e versao opcional."""

RSID_PATTERN = re.compile(r"^rs\d{1,12}$", re.IGNORECASE)
"""rsID dbSNP."""

HGVS_FORBIDDEN_CHARS: str = "[]()\"' \t\n\r"
"""Caracteres recusados numa notacao HGVS usada como termo Entrez.

Os parenteses e parenteses retos sao operadores de campo e agrupamento do
Entrez, por isso uma notacao que os contenha deixaria de identificar uma unica
variante.
"""

CLINICAL_DISCLAIMER: str = (
    "ClinVar records are clinical interpretations submitted by external "
    "laboratories and expert panels, retrieved here unchanged. They are not a "
    "HelixScope assessment, not a diagnosis and not a clinical recommendation. "
    "Interpretation of a variant for any clinical purpose requires a qualified "
    "professional and the primary submitted evidence."
)


class ClinVarError(RuntimeError):
    """Falha ao consultar ou interpretar ClinVar.

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
    """Recusa redirects para fora do host E-utilities."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        if not url_is_allowed(newurl):
            raise ClinVarError(
                f"Refusing a ClinVar redirect to a non-allowlisted URL: {newurl}.",
                "INVALID_INPUT",
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def url_is_allowed(url: str) -> bool:
    """Valida esquema, host e prefixo de caminho contra a allowlist E-utilities.

    Args:
        url: URL absoluta.

    Returns:
        True apenas para HTTPS em eutils.ncbi.nlm.nih.gov sob /entrez/eutils/.

    Raises:
        Nenhum.
    """
    parsed = urllib.parse.urlparse(str(url or ""))
    if parsed.scheme != "https":
        return False
    if (parsed.hostname or "").lower() not in ALLOWED_HOSTS:
        return False
    return (parsed.path or "").startswith(ALLOWED_PATH_PREFIX)


def network_disclosure() -> dict:
    """Descreve o que sai da maquina numa consulta ClinVar.

    Args:
        Nenhum.

    Returns:
        Dict com host, dados enviados, dados nao enviados e politica de uso.

    Raises:
        Nenhum.
    """
    return {
        "remote": True,
        "service": "NCBI ClinVar via E-utilities",
        "host": EUTILS_HOST,
        "transport": "HTTPS",
        "data_sent": [
            "the ClinVar search term you submit (coordinate, HGVS notation, rsID or accession)",
            "contact e-mail required by the NCBI usage policy",
            "NCBI API key when one is configured in the environment",
        ],
        "data_not_sent": [
            "any other sequence loaded in HelixScope",
            "local file names or paths",
            "reference genome files",
        ],
        "user_agent": USER_AGENT,
        "rate_limit_policy": (
            "NCBI allows 3 requests/second without an API key and 10 "
            "requests/second with one. HelixScope serialises ClinVar calls and "
            "backs off on HTTP 429."
        ),
        "terms": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        "disclaimer": CLINICAL_DISCLAIMER,
    }


def build_search_term(
    *,
    assembly: str = "",
    contig: str = "",
    position_1based: Optional[int] = None,
    hgvs: str = "",
    rsid: str = "",
    accession: str = "",
) -> str:
    """Monta um termo de busca ClinVar a partir da identidade disponivel.

    Args:
        assembly: Assembly declarada (usada apenas em busca por coordenada).
        contig: Cromossoma sem prefixo chr.
        position_1based: Posicao 1-based.
        hgvs: Notacao HGVS completa, preferida quando existe.
        rsid: rsID.
        accession: Accession VCV ou RCV.

    Returns:
        Termo de busca no idioma do Entrez.

    Raises:
        ClinVarError: INVALID_INPUT quando nada identificavel foi fornecido,
            quando a coordenada e dada sem assembly, ou quando um identificador
            contem operadores da linguagem de busca do Entrez.

    Nota biologica:
        A busca por accession e a mais precisa. A busca por coordenada e a menos
        precisa: uma posicao pode conter varias variantes diferentes, por isso o
        resultado devolve todas em vez de escolher uma.

    Nota de seguranca:
        Cada identificador e validado contra um padrao fechado. Um token com
        operadores Entrez (`OR`, `[chr]`, parenteses) alargaria a busca e as
        submissoes clinicas recuperadas seriam atribuidas a variante errada.
    """
    acc = str(accession or "").strip()
    if acc:
        if not ACCESSION_PATTERN.match(acc):
            raise ClinVarError(
                f"'{acc}' is not a ClinVar accession. Expected VCV, RCV or SCV "
                "followed by digits and an optional version.",
                "INVALID_INPUT",
            )
        return acc.upper()
    rs = str(rsid or "").strip()
    if rs:
        if not RSID_PATTERN.match(rs):
            raise ClinVarError(
                f"'{rs}' is not an rsID. Expected 'rs' followed by digits.",
                "INVALID_INPUT",
            )
        return rs.lower()
    notation = str(hgvs or "").strip()
    if notation:
        if len(notation) > MAX_TERM_LENGTH:
            raise ClinVarError(
                f"HGVS notation exceeds {MAX_TERM_LENGTH} characters.",
                "RESOURCE_LIMIT",
            )
        if any(char in notation for char in HGVS_FORBIDDEN_CHARS):
            raise ClinVarError(
                "HGVS notation contains characters that are Entrez search "
                "operators and would broaden the query beyond this variant.",
                "INVALID_INPUT",
            )
        return notation
    chrom = str(contig or "").strip()
    if chrom and position_1based is not None:
        if not str(assembly or "").strip():
            raise ClinVarError(
                "A ClinVar coordinate search needs the assembly: the same "
                "chromosomal position is a different variant in GRCh37 and "
                "GRCh38.",
                "INVALID_INPUT",
            )
        if not CONTIG_TOKEN_PATTERN.match(chrom):
            raise ClinVarError(
                f"'{chrom}' is not a safe chromosome token for an Entrez query.",
                "INVALID_INPUT",
            )
        position = int(position_1based)
        if position < 1:
            raise ClinVarError("position_1based must be >= 1.", "INVALID_INPUT")
        return f"{chrom}[chr] AND {position}[chrpos]"
    raise ClinVarError(
        "No usable ClinVar query. Provide an accession, rsID, HGVS notation or "
        "assembly plus coordinate.",
        "INVALID_INPUT",
    )


def _eutils_params(email: str) -> dict:
    contact = str(email or "").strip()
    if not ncbi_fetch._valid_entrez_email(contact):  # noqa: SLF001 - shared validation
        raise ClinVarError(
            "NCBI requires a valid contact e-mail for E-utilities requests.",
            "INVALID_INPUT",
        )
    params = {"tool": "HelixScope", "email": contact}
    key = ncbi_fetch.api_key_from_environment()
    if key:
        params["api_key"] = key
    return params


def search_clinvar(
    *,
    term: str,
    email: str,
    retmax: int = 10,
    urlopen_fn=None,
) -> dict:
    """Resolve um termo para UIDs ClinVar. FAZ REDE.

    Args:
        term: Termo de busca Entrez.
        email: E-mail de contacto exigido pela politica do NCBI.
        retmax: Numero maximo de UIDs.
        urlopen_fn: Injecao para teste.

    Returns:
        Dict uids, count, term, translation, retrieved_at_utc.

    Raises:
        ClinVarError: INVALID_INPUT, TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE,
            PARSING_ERROR ou NETWORK_ERROR. Lista vazia NAO e erro.

    Nota biologica:
        Zero UIDs significa "nenhuma submissao recuperada para este termo", nao
        "variante benigna" e nao "variante inexistente".
    """
    text = str(term or "").strip()
    if not text:
        raise ClinVarError("ClinVar search term is empty.", "INVALID_INPUT")
    if len(text) > MAX_TERM_LENGTH:
        raise ClinVarError(
            f"ClinVar search term exceeds {MAX_TERM_LENGTH} characters.",
            "RESOURCE_LIMIT",
        )
    limit = max(1, min(int(retmax), MAX_UIDS))
    params = _eutils_params(email)
    params.update({"db": "clinvar", "retmode": "json", "term": text, "retmax": str(limit)})
    url = f"{EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    result = payload.get("esearchresult") if isinstance(payload, Mapping) else None
    if not isinstance(result, Mapping):
        raise ClinVarError("ClinVar esearch response is malformed.", "PARSING_ERROR")
    uids = [str(item) for item in (result.get("idlist") or []) if str(item).strip().isdigit()]
    return {
        "uids": uids,
        "count": int(str(result.get("count") or "0") or 0),
        "term": text,
        "translation": str(result.get("querytranslation") or ""),
        "retrieved_at_utc": provenance.utc_now(),
    }


def _classification_block(payload: Any, label: str) -> Optional[dict]:
    if not isinstance(payload, Mapping):
        return None
    description = str(payload.get("description") or "").strip()
    review = str(payload.get("review_status") or "").strip()
    if not description and not review:
        return None
    traits = []
    for item in payload.get("trait_set") or []:
        if not isinstance(item, Mapping):
            continue
        xrefs = [
            {
                "db_source": str((ref or {}).get("db_source") or ""),
                "db_id": str((ref or {}).get("db_id") or ""),
            }
            for ref in (item.get("trait_xrefs") or [])
            if isinstance(ref, Mapping)
        ]
        traits.append(
            {"trait_name": str(item.get("trait_name") or ""), "trait_xrefs": xrefs}
        )
    last_evaluated = str(payload.get("last_evaluated") or "").strip()
    if last_evaluated.startswith("1/01/01"):
        last_evaluated = ""
    lowered = f"{description} {review}".lower()
    return {
        "classification_type": label,
        "description": description,
        "review_status": review,
        "review_status_rank": (
            REVIEW_STATUS_ORDER.index(review.lower())
            if review.lower() in REVIEW_STATUS_ORDER
            else len(REVIEW_STATUS_ORDER)
        ),
        "last_evaluated": last_evaluated,
        "fda_recognized_database": str(payload.get("fda_recognized_database") or ""),
        "conditions": traits,
        "has_conflict": any(marker in lowered for marker in CONFLICT_MARKERS),
        "aggregated_by": (
            "NCBI ClinVar aggregates submitted SCV classifications into this "
            "value using submission review status. HelixScope does not "
            "recompute it."
        ),
    }


def parse_clinvar_summary(payload: Any, *, uid: str = "") -> dict:
    """Converte um esummary ClinVar em evidencia estruturada e rastreavel.

    Args:
        payload: JSON desserializado do esummary.
        uid: UID esperado; se vazio usa o primeiro do resultado.

    Returns:
        Dict com accession, titulo, classificacoes (germline, clinical impact,
        oncogenicity), condicoes, localizacoes por assembly, submissoes de
        suporte, deteccao de conflito, disclaimer e proveniencia.

    Raises:
        ClinVarError: PARSING_ERROR estrutura inesperada; NOT_FOUND se o UID nao
            estiver no resultado.

    Nota biologica:
        As localizacoes vem por assembly e podem citar um accession de patch
        diferente do usado localmente (por exemplo GCF_000001405.38 em vez de
        .40). Isso e preservado como veio, para que a divergencia fique visivel
        em vez de ser normalizada em silencio.
    """
    if not isinstance(payload, Mapping):
        raise ClinVarError("ClinVar esummary response is not a JSON object.", "PARSING_ERROR")
    result = payload.get("result")
    if not isinstance(result, Mapping):
        raise ClinVarError("ClinVar esummary response has no result block.", "PARSING_ERROR")
    key = str(uid or "").strip()
    if not key:
        candidates = [str(item) for item in (result.get("uids") or [])]
        if not candidates:
            raise ClinVarError("ClinVar esummary returned no UID.", "NOT_FOUND")
        key = candidates[0]
    record = result.get(key)
    if not isinstance(record, Mapping):
        raise ClinVarError(
            f"ClinVar esummary has no record for UID {key}.", "NOT_FOUND"
        )
    classifications = []
    for field, label in (
        ("germline_classification", "germline"),
        ("clinical_impact_classification", "somatic clinical impact"),
        ("oncogenicity_classification", "oncogenicity"),
    ):
        block = _classification_block(record.get(field), label)
        if block:
            classifications.append(block)
    classifications.sort(key=lambda item: item["review_status_rank"])
    locations = []
    variation_names = []
    cdna_changes = []
    for item in record.get("variation_set") or []:
        if not isinstance(item, Mapping):
            continue
        name = str(item.get("variation_name") or "")
        if name:
            variation_names.append(name)
        cdna = str(item.get("cdna_change") or "")
        if cdna:
            cdna_changes.append(cdna)
        for loc in item.get("variation_loc") or []:
            if not isinstance(loc, Mapping):
                continue
            locations.append(
                {
                    "status": str(loc.get("status") or ""),
                    "assembly_name": str(loc.get("assembly_name") or ""),
                    "assembly_accession": str(loc.get("assembly_acc_ver") or ""),
                    "chromosome": str(loc.get("chr") or ""),
                    "band": str(loc.get("band") or ""),
                    "start_1based": _optional_int(loc.get("start")),
                    "stop_1based": _optional_int(loc.get("stop")),
                    "ref": str(loc.get("ref") or ""),
                    "alt": str(loc.get("alt") or ""),
                }
            )
    xref_rows = []
    for item in record.get("variation_set") or []:
        if not isinstance(item, Mapping):
            continue
        for ref in item.get("variation_xrefs") or []:
            if isinstance(ref, Mapping):
                xref_rows.append(
                    {
                        "db_source": str(ref.get("db_source") or ""),
                        "db_id": str(ref.get("db_id") or ""),
                    }
                )
    submissions = record.get("supporting_submissions")
    scv = []
    rcv = []
    if isinstance(submissions, Mapping):
        scv = [str(item) for item in (submissions.get("scv") or [])]
        rcv = [str(item) for item in (submissions.get("rcv") or [])]
    genes = [
        {
            "symbol": str((item or {}).get("symbol") or ""),
            "gene_id": str((item or {}).get("GeneID") or (item or {}).get("geneid") or ""),
        }
        for item in (record.get("genes") or [])
        if isinstance(item, Mapping)
    ]
    accession = str(record.get("accession") or "")
    conflict = any(block.get("has_conflict") for block in classifications)
    distinct_descriptions = sorted(
        {
            block["description"]
            for block in classifications
            if block["classification_type"] == "germline" and block["description"]
        }
    )
    return {
        "source": "NCBI ClinVar",
        "evidence_status": "RETRIEVED",
        "role": (
            "External submitted clinical evidence. Not a HelixScope "
            "pathogenicity prediction."
        ),
        "uid": key,
        "accession": accession,
        "accession_version": str(record.get("accession_version") or ""),
        "title": str(record.get("title") or ""),
        "object_type": str(record.get("obj_type") or ""),
        "record_status": str(record.get("record_status") or ""),
        "protein_change": str(record.get("protein_change") or ""),
        "molecular_consequences": [
            str(item) for item in (record.get("molecular_consequence_list") or [])
        ],
        "genes": genes,
        "variation_names": [name for name in variation_names if name],
        "cdna_changes": cdna_changes,
        "classifications": classifications,
        "germline_descriptions": distinct_descriptions,
        "has_conflict": conflict,
        "conflict_note": (
            "ClinVar reports conflicting submitted classifications for this "
            "variant. HelixScope shows the conflict and does not resolve it into "
            "a single conclusion."
            if conflict
            else "ClinVar did not flag a conflict among submitted classifications."
        ),
        "locations": locations,
        "cross_references": xref_rows,
        "supporting_scv": scv,
        "supporting_rcv": rcv,
        "submission_count": len(scv),
        "url": (
            f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{key}/" if key else ""
        ),
        "disclaimer": CLINICAL_DISCLAIMER,
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def record_matches_variant(
    record: Mapping[str, Any],
    *,
    assembly: str,
    contig: str,
    position_1based: int,
    ref: str = "",
    alt: str = "",
) -> dict:
    """Confere se um registo ClinVar cai na coordenada consultada.

    Args:
        record: Resultado de `parse_clinvar_summary`.
        assembly: Assembly da variante (ex. GRCh38.p14).
        contig: Cromossoma sem prefixo chr.
        position_1based: Posicao 1-based da variante normalizada.
        ref: Alelo de referencia consultado, quando conhecido.
        alt: Alelo alternativo consultado, quando conhecido.

    Returns:
        Dict identity_match (posicao), allele_match (True, False ou None quando
        ClinVar nao publica os alelos no esummary), matched_location e reason.

    Raises:
        Nenhum.

    Nota biologica:
        A busca textual do Entrez e difusa: pesquisar uma notacao HGVS pode
        devolver registos de OUTRAS variantes do mesmo gene. Sem esta
        verificacao, a interface mostraria a classificacao clinica da variante
        errada, que e o pior erro possivel nesta camada.

        A comparacao usa o nome da assembly (GRCh38/GRCh37) e nao o accession de
        patch: ClinVar anota GRCh38 em GCF_000001405.38 enquanto HelixScope pode
        ter GCF_000001405.40 instalado. Patches acrescentam sequencias novas mas
        nao deslocam coordenadas dos cromossomas primarios, por isso os nomes de
        assembly sao comparaveis e a diferenca de patch e reportada como nota.

        Posicao NAO e identidade: numa mesma base podem existir uma delecao
        patogenica, uma substituicao benigna e outra com classificacoes
        conflitantes. Quando o esummary nao publica ref/alt, `allele_match` fica
        None e o resultado avisa que os registos na posicao podem ser de alelos
        diferentes.
    """
    target_assembly = str(assembly or "").strip().upper().replace(" ", "")
    base = ""
    for candidate in ("GRCH38", "GRCH37"):
        if target_assembly.startswith(candidate):
            base = candidate
            break
    chrom = str(contig or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    try:
        position = int(position_1based)
    except (TypeError, ValueError):
        return {
            "identity_match": False,
            "matched_location": None,
            "reason": "No comparable 1-based position was supplied.",
        }
    if not base:
        return {
            "identity_match": False,
            "allele_match": None,
            "matched_location": None,
            "reason": (
                f"Assembly {assembly} is not comparable with the GRCh37/GRCh38 "
                "coordinates published by ClinVar."
            ),
        }
    wanted_ref = str(ref or "").strip().upper()
    wanted_alt = str(alt or "").strip().upper()
    for location in record.get("locations") or []:
        if not isinstance(location, Mapping):
            continue
        loc_assembly = str(location.get("assembly_name") or "").strip().upper()
        if loc_assembly != base:
            continue
        if str(location.get("chromosome") or "").strip() != chrom:
            continue
        start = location.get("start_1based")
        stop = location.get("stop_1based")
        if start is None:
            continue
        end = stop if stop is not None else start
        if int(start) <= position <= int(end):
            patch_note = ""
            loc_accession = str(location.get("assembly_accession") or "")
            if loc_accession:
                patch_note = (
                    f" ClinVar annotated this position on {loc_accession}; patch "
                    "level may differ from the locally installed assembly "
                    "without changing primary-chromosome coordinates."
                )
            loc_ref = str(location.get("ref") or "").strip().upper()
            loc_alt = str(location.get("alt") or "").strip().upper()
            allele_match: Optional[bool] = None
            allele_note = (
                " ClinVar's summary does not publish explicit ref/alt alleles "
                "for this location, so allele identity is UNVERIFIED: distinct "
                "variants at the same base can carry different classifications."
            )
            if loc_ref and loc_alt and wanted_ref and wanted_alt:
                allele_match = loc_ref == wanted_ref and loc_alt == wanted_alt
                allele_note = (
                    f" ClinVar alleles {loc_ref}>{loc_alt} "
                    f"{'match' if allele_match else 'do not match'} the queried "
                    f"{wanted_ref}>{wanted_alt}."
                )
            return {
                "identity_match": True,
                "allele_match": allele_match,
                "matched_location": dict(location),
                "reason": (
                    f"ClinVar reports this record at {base} chromosome {chrom} "
                    f"{start}-{end}, which contains the queried position "
                    f"{position}.{patch_note}{allele_note}"
                ),
            }
    return {
        "identity_match": False,
        "allele_match": None,
        "matched_location": None,
        "reason": (
            f"No {base} location in this ClinVar record contains chromosome "
            f"{chrom} position {position}. The Entrez text search returned it as "
            "a related record, not as evidence for the queried variant."
        ),
    }


def fetch_clinvar_records(
    *,
    term: str,
    email: str,
    retmax: int = 5,
    urlopen_fn=None,
    assembly: str = "",
    contig: str = "",
    position_1based: Optional[int] = None,
    ref: str = "",
    alt: str = "",
) -> dict:
    """Recupera os registos ClinVar de um termo. FAZ REDE.

    Args:
        term: Termo de busca Entrez.
        email: E-mail de contacto.
        retmax: Numero maximo de registos a resumir.
        urlopen_fn: Injecao para teste.
        assembly: Assembly da variante consultada, para verificar identidade.
        contig: Cromossoma da variante consultada.
        position_1based: Posicao 1-based da variante consultada.

    Returns:
        Dict status (AVAILABLE ou NOT_FOUND), records, count, term, conflito
        agregado e disclaimer. Quando assembly/contig/position sao fornecidos,
        cada registo recebe `identity_match` e o resultado traz
        `identity_matched_count`.

    Raises:
        ClinVarError: categorias de rede e parsing. NOT_FOUND e devolvido como
            status, nao como excecao, para nao destruir as outras camadas.

    Nota biologica:
        Uma busca por coordenada pode devolver varias variantes distintas na
        mesma posicao, e uma busca textual pode devolver variantes de outra
        posicao do mesmo gene. Todos os registos sao devolvidos e rotulados;
        nenhum e descartado e nenhum e promovido a evidencia da variante
        consultada sem casar a coordenada.
    """
    search = search_clinvar(term=term, email=email, retmax=retmax, urlopen_fn=urlopen_fn)
    uids = search["uids"]
    if not uids:
        return {
            "status": "NOT_FOUND",
            "source": "NCBI ClinVar",
            "evidence_status": "RETRIEVED",
            "records": [],
            "count": 0,
            "term": search["term"],
            "query_translation": search["translation"],
            "has_conflict": False,
            "message": (
                "No ClinVar record retrieved for this query. Absence of a "
                "submission is not evidence that the variant is benign."
            ),
            "disclaimer": CLINICAL_DISCLAIMER,
            "retrieved_at_utc": provenance.utc_now(),
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
    params = _eutils_params(email)
    params.update({"db": "clinvar", "retmode": "json", "id": ",".join(uids)})
    url = f"{EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    records = []
    for uid in uids:
        try:
            records.append(parse_clinvar_summary(payload, uid=uid))
        except ClinVarError as exc:
            if exc.category == "NOT_FOUND":
                continue
            raise
    if not records:
        raise ClinVarError(
            "ClinVar returned UIDs but no summary could be parsed for them.",
            "PARSING_ERROR",
        )
    verified = position_1based is not None and bool(contig) and bool(assembly)
    matched = 0
    for record in records:
        if not verified:
            record["identity_match"] = None
            record["allele_match"] = None
            record["identity_match_reason"] = (
                "No queried coordinate was supplied, so this record was not "
                "verified against a variant position. Treat it as a search hit, "
                "not as confirmed evidence for a specific variant."
            )
            continue
        outcome = record_matches_variant(
            record,
            assembly=assembly,
            contig=contig,
            position_1based=int(position_1based),
            ref=ref,
            alt=alt,
        )
        record["identity_match"] = bool(outcome["identity_match"])
        record["allele_match"] = outcome["allele_match"]
        record["identity_match_reason"] = str(outcome["reason"])
        record["matched_location"] = outcome["matched_location"]
        if record["identity_match"]:
            matched += 1
    matched_records = [record for record in records if record.get("identity_match")]
    return {
        "status": "AVAILABLE",
        "identity_verified": verified,
        "identity_matched_count": matched if verified else None,
        "identity_note": (
            (
                f"{matched} of {len(records)} retrieved record(s) fall on the "
                "queried coordinate. Records that do not match are shown as "
                "related search hits and are not evidence for this variant."
                + (
                    " More than one distinct variant is recorded at this "
                    "position, so each classification belongs to its own "
                    "record and they must not be merged."
                    if matched > 1
                    else ""
                )
            )
            if verified
            else (
                "Records were not coordinate-verified because no normalized "
                "variant position was supplied."
            )
        ),
        "allele_verified_count": (
            sum(1 for record in matched_records if record.get("allele_match") is True)
            if verified
            else None
        ),
        "allele_unverified_count": (
            sum(1 for record in matched_records if record.get("allele_match") is None)
            if verified
            else None
        ),
        "has_conflict_at_variant": any(
            record.get("has_conflict") for record in matched_records
        )
        if verified
        else None,
        "source": "NCBI ClinVar",
        "evidence_status": "RETRIEVED",
        "records": records,
        "count": len(records),
        "total_matches": search["count"],
        "term": search["term"],
        "query_translation": search["translation"],
        "has_conflict": any(record.get("has_conflict") for record in records),
        "message": (
            f"{len(records)} ClinVar record(s) retrieved. Classifications belong "
            "to the submitters listed in each record."
        ),
        "disclaimer": CLINICAL_DISCLAIMER,
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def _http_json(url: str, *, urlopen_fn=None) -> Any:
    if not url_is_allowed(url):
        raise ClinVarError(
            "Refusing a non-allowlisted ClinVar URL.", "INVALID_INPUT"
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
                final = str(getattr(handle, "geturl", lambda: url)() or url)
            break
        except HTTPError as exc:
            error = _http_error(exc)
            if error.category in {"RATE_LIMITED", "SERVICE_UNAVAILABLE"} and attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, error.retry_after_s))
                continue
            raise error from exc
        except ClinVarError:
            raise
        except (URLError, TimeoutError, socket.timeout, OSError) as exc:
            message = str(exc).lower()
            timed_out = isinstance(exc, (TimeoutError, socket.timeout)) or "timed out" in message
            error = ClinVarError(
                (
                    "ClinVar request timed out. This is TIMEOUT, not an absence "
                    "of clinical submissions."
                )
                if timed_out
                else (
                    f"ClinVar network error: {exc}. This is NETWORK_ERROR, not "
                    "NOT_FOUND."
                ),
                "TIMEOUT" if timed_out else "NETWORK_ERROR",
            )
            if attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, None))
                continue
            raise error from exc
    if final and not url_is_allowed(final):
        raise ClinVarError(
            "ClinVar answered from a non-allowlisted URL.", "INVALID_INPUT"
        )
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    if len(raw) > MAX_RESPONSE_BYTES:
        raise ClinVarError(
            f"ClinVar response exceeds {MAX_RESPONSE_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    text = raw.decode("utf-8", errors="replace")
    if not text.strip():
        raise ClinVarError("ClinVar response is empty.", "PARSING_ERROR")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ClinVarError("ClinVar response is not valid JSON.", "PARSING_ERROR") from exc


def _backoff_seconds(attempt: int, retry_after_s: Optional[float]) -> float:
    if retry_after_s is not None and retry_after_s >= 0:
        return min(float(retry_after_s), MAX_BACKOFF_S)
    return min(BACKOFF_BASE_S ** attempt, MAX_BACKOFF_S)


def _http_error(exc: HTTPError) -> ClinVarError:
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
        return ClinVarError("ClinVar endpoint returned HTTP 404.", "NOT_FOUND")
    if code == 429:
        return ClinVarError(
            "NCBI rate limited the ClinVar request (HTTP 429). Configure "
            "NCBI_API_KEY to raise the allowance.",
            "RATE_LIMITED",
            retry_after,
        )
    if code == 400:
        return ClinVarError(
            f"NCBI rejected the ClinVar request (HTTP 400): {exc}.", "INVALID_INPUT"
        )
    if code >= 500:
        return ClinVarError(
            f"ClinVar service unavailable (HTTP {code}).", "SERVICE_UNAVAILABLE", retry_after
        )
    return ClinVarError(f"ClinVar HTTP {code}: {exc}.", "SERVICE_UNAVAILABLE", retry_after)


def _optional_int(value: Any) -> Optional[int]:
    text = str(value if value is not None else "").strip()
    if not text:
        return None
    try:
        return int(text)
    except (TypeError, ValueError):
        return None
