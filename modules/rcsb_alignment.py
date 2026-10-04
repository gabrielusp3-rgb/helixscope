"""Cliente real da RCSB PDB Structure Alignment API.

Alinhamento estrutural pairwise e calculado pelo servico oficial
(https://alignment.rcsb.org). HelixScope nao implementa TM-align nem FATCAT:
submete a consulta, espera o ticket assincrono e interpreta o JSON.

Documentacao revalidada 2026-08-28:
    User Guide: https://alignment.rcsb.org/
    API Reference: https://alignment.rcsb.org/api-reference.html
    Web APIs Overview (rate limit 429): https://www.rcsb.org/docs/programmatic-access/web-apis-overview
    Publicacao: Bittrich et al., Bioinformatics 2024, 40:btae370.

Metodos expostos (nomes de API oficiais no campo method.name):
    fatcat-rigid, fatcat-flexible, tm-align, ce, ce-cp.

RMSD da API e sobre pares de C-alpha estruturalmente equivalentes do bloco.
HelixScope nao funde RMSD de bloco com RMSD do summary global, nem TM-scores
com normalizacoes diferentes.

Nenhuma funcao aqui importa Streamlit. Nao se enviam URLs controladas pelo
utilizador ao servico (SSRF via terceiro).
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

ALIGNMENT_HOST: str = "alignment.rcsb.org"
API_BASE: str = "https://alignment.rcsb.org/api/v1/structures"
SUBMIT_URL: str = f"{API_BASE}/submit"
RESULTS_URL: str = f"{API_BASE}/results"
DOCS_URL: str = "https://alignment.rcsb.org/"
RATE_LIMIT_DOCS: str = (
    "https://www.rcsb.org/docs/programmatic-access/web-apis-overview"
)

ALLOWED_HOSTS: frozenset[str] = frozenset({ALIGNMENT_HOST})
ALLOWED_PATH_PREFIXES: tuple[str, ...] = ("/api/v1/structures",)

USER_AGENT: str = (
    f"HelixScope/{provenance.HELIXSCOPE_VERSION} (bioinformatics; RCSB Alignment API client)"
)
HTTP_TIMEOUT_S: float = 60.0
POLL_TIMEOUT_S: float = 180.0
POLL_INTERVAL_S: float = 1.5
MAX_RESPONSE_BYTES: int = 8 * 1024 * 1024
MAX_RETRIES: int = 3
BACKOFF_BASE_S: float = 1.5
MAX_BACKOFF_S: float = 20.0

MODE_PAIRWISE: str = "pairwise"

METHOD_FATCAT_RIGID: str = "fatcat-rigid"
METHOD_FATCAT_FLEXIBLE: str = "fatcat-flexible"
METHOD_TMALIGN: str = "tm-align"
METHOD_CE: str = "ce"
METHOD_CE_CP: str = "ce-cp"

RIGID_METHODS: frozenset[str] = frozenset(
    {METHOD_FATCAT_RIGID, METHOD_TMALIGN, METHOD_CE}
)
FLEXIBLE_METHODS: frozenset[str] = frozenset(
    {METHOD_FATCAT_FLEXIBLE, METHOD_CE_CP}
)

METHOD_CATALOG: tuple[dict, ...] = (
    {
        "api_name": METHOD_FATCAT_RIGID,
        "display_name": "jFATCAT-rigid",
        "kind": "rigid",
        "citation": "Ye and Godzik 2003 FATCAT; Li et al. 2020 Java port (jFATCAT)",
        "selectable": True,
    },
    {
        "api_name": METHOD_FATCAT_FLEXIBLE,
        "display_name": "jFATCAT-flexible",
        "kind": "flexible",
        "citation": "Ye and Godzik 2003 FATCAT flexible (twists between AFPs)",
        "selectable": True,
    },
    {
        "api_name": METHOD_TMALIGN,
        "display_name": "TM-align",
        "kind": "rigid",
        "citation": "Zhang and Skolnick 2005, Nucleic Acids Res 33:2302",
        "selectable": True,
    },
    {
        "api_name": METHOD_CE,
        "display_name": "jCE",
        "kind": "rigid",
        "citation": "Shindyalov and Bourne 1998, Protein Eng 11:739",
        "selectable": True,
    },
    {
        "api_name": METHOD_CE_CP,
        "display_name": "jCE-CP",
        "kind": "flexible",
        "citation": "CE with circular permutation; may return multiple blocks",
        "selectable": True,
    },
)

ENTRY_ID_PATTERN = re.compile(
    r"^(?:[0-9][A-Za-z0-9]{3}|AF_[A-Za-z0-9._-]{3,40}|MA_[A-Za-z0-9._-]{3,40})$"
)
ALPHAFOLD_HYPHEN_PATTERN = re.compile(r"^AF-([A-Z0-9]+)-F([0-9]+)$", re.IGNORECASE)
RCSB_CSM_AF_PATTERN = re.compile(r"^AF_AF([A-Z0-9]+)F([0-9]+)$", re.IGNORECASE)
UNIPROT_ACCESSION_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3}[0-9])?)$",
    re.IGNORECASE,
)
ASYM_ID_PATTERN = re.compile(r"^[A-Za-z0-9]{1,4}$")
UUID_PATTERN = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)


class AlignmentApiError(RuntimeError):
    """Falha ao consultar a Alignment API.

    Attributes:
        category: INVALID_INPUT, TIMEOUT, RATE_LIMITED, SERVICE_UNAVAILABLE,
            NETWORK_ERROR, PARSING_ERROR, NOT_FOUND, RESOURCE_LIMIT ou ERROR.
        retry_after_s: Segundos sugeridos pelo servico, quando informados.
    """

    def __init__(self, message: str, category: str, retry_after_s: Optional[float] = None) -> None:
        super().__init__(message)
        self.category = str(category or "SERVICE_UNAVAILABLE").strip().upper()
        self.retry_after_s = retry_after_s


class _AllowlistRedirect(HTTPRedirectHandler):
    """Recusa redirects para fora de alignment.rcsb.org."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        if not url_is_allowed(newurl):
            raise AlignmentApiError(
                f"Refusing an Alignment API redirect to a non-allowlisted URL: {newurl}.",
                "INVALID_INPUT",
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def url_is_allowed(url: str) -> bool:
    """Valida HTTPS, host e prefixo de caminho da Alignment API.

    Args:
        url: URL absoluta.

    Returns:
        True apenas para https://alignment.rcsb.org/api/v1/structures...

    Raises:
        Nenhum.
    """
    parsed = urllib.parse.urlparse(str(url or ""))
    if parsed.scheme != "https":
        return False
    host = (parsed.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        return False
    path = parsed.path or ""
    return any(path.startswith(prefix) for prefix in ALLOWED_PATH_PREFIXES)


def supported_methods() -> list[dict]:
    """Metodos que o selector pode mostrar. Nenhum algoritmo inventado.

    Args:
        Nenhum.

    Returns:
        Lista de dicts com api_name, display_name, kind, citation.

    Raises:
        Nenhum.
    """
    return [dict(row) for row in METHOD_CATALOG if row.get("selectable")]


def method_kind(api_name: str) -> str:
    """Devolve rigid ou flexible para o nome de API.

    Args:
        api_name: fatcat-rigid, tm-align, etc.

    Returns:
        "rigid" ou "flexible".

    Raises:
        AlignmentApiError: INVALID_INPUT se o metodo nao estiver no catalogo.
    """
    name = str(api_name or "").strip()
    for row in METHOD_CATALOG:
        if row["api_name"] == name:
            return str(row["kind"])
    raise AlignmentApiError(
        f"Unsupported alignment method {api_name!r}. HelixScope only exposes "
        "methods documented by the RCSB Alignment API.",
        "INVALID_INPUT",
    )


def display_name_for(api_name: str) -> str:
    """Nome de apresentacao (jFATCAT-rigid, TM-align, ...).

    Args:
        api_name: Nome enviado no JSON.

    Returns:
        display_name do catalogo.

    Raises:
        AlignmentApiError: INVALID_INPUT se desconhecido.
    """
    name = str(api_name or "").strip()
    for row in METHOD_CATALOG:
        if row["api_name"] == name:
            return str(row["display_name"])
    raise AlignmentApiError(f"Unsupported alignment method {api_name!r}.", "INVALID_INPUT")


def network_disclosure() -> dict:
    """O que sai da maquina ao pedir um alinhamento RCSB.

    Args:
        Nenhum.

    Returns:
        Host, endpoints, dados enviados, politica de rate limit.

    Raises:
        Nenhum.
    """
    return {
        "host": ALIGNMENT_HOST,
        "endpoints": [SUBMIT_URL, RESULTS_URL],
        "docs": DOCS_URL,
        "data_sent": [
            "PDB or RCSB CSM entry identifiers",
            "polymer chain asym_id",
            "alignment method name",
        ],
        "not_sent": [
            "arbitrary user URLs (refused; would be SSRF via a third party)",
            "genomic sequences",
            "email",
        ],
        "rate_limit_policy": (
            "RCSB returns HTTP 429 when the shared IP exceeds the API budget. "
            f"HelixScope retries at most {MAX_RETRIES} times with exponential "
            f"backoff. See {RATE_LIMIT_DOCS}."
        ),
        "async": True,
        "engine_location": "remote RCSB Alignment API, not a local binary",
    }


def rcsb_csm_id_from_uniprot(accession: str, fragment: int = 1) -> str:
    """Converte UniProt + fragmento no CSM ID da RCSB.

    Args:
        accession: Accession UniProtKB (ex. P01308).
        fragment: Fragmento AlphaFold, habitualmente 1.

    Returns:
        Identificador AF_AF{accession}F{fragment} (ex. AF_AFP01308F1).

    Raises:
        AlignmentApiError: INVALID_INPUT.

    Nota:
        Documentacao RCSB (Computed Structure Models): AF-B3EWR1-F1 na
        AlphaFold DB corresponde a AF_AFB3EWR1F1 na RCSB. HelixScope nao
        inventa o mapeamento; aplica essa regra publicada.
        https://www.rcsb.org/docs/general-help/computed-structure-models-and-rcsborg
    """
    text = str(accession or "").strip().upper()
    if not UNIPROT_ACCESSION_PATTERN.match(text):
        raise AlignmentApiError(
            "UniProt accession is invalid; gene symbols are not protein identity.",
            "INVALID_INPUT",
        )
    frag = int(fragment)
    if frag < 1 or frag > 99:
        raise AlignmentApiError("AlphaFold fragment number is out of range.", "INVALID_INPUT")
    return f"AF_AF{text}F{frag}"


def parse_rcsb_csm_uniprot(entry_id: str) -> Optional[tuple[str, int]]:
    """Extrai UniProt e fragmento de um CSM AF_AF...F n.

    Args:
        entry_id: Identificador RCSB.

    Returns:
        (accession, fragment) ou None.

    Raises:
        Nenhum.
    """
    match = RCSB_CSM_AF_PATTERN.match(str(entry_id or "").strip())
    if not match:
        return None
    return match.group(1).upper(), int(match.group(2))


def validate_entry_id(entry_id: str) -> str:
    """Aceita PDB ID de 4 caracteres, CSM RCSB (AF_/MA_) ou AF-P01308-F1.

    Args:
        entry_id: Identificador.

    Returns:
        Identificador canonico (PDB em maiusculas; AF hyphen vira AF_AF...Fn).

    Raises:
        AlignmentApiError: INVALID_INPUT.
    """
    text = str(entry_id or "").strip()
    hyphen = ALPHAFOLD_HYPHEN_PATTERN.match(text)
    if hyphen:
        return rcsb_csm_id_from_uniprot(hyphen.group(1), int(hyphen.group(2)))
    if not ENTRY_ID_PATTERN.match(text):
        raise AlignmentApiError(
            "Structure identifier must be a 4-character PDB ID, an RCSB "
            "computed-structure-model id (AF_ / MA_), or an AlphaFold DB id "
            "(AF-P01308-F1). Arbitrary URLs and gene symbols are refused.",
            "INVALID_INPUT",
        )
    if len(text) == 4:
        return text.upper()
    return text


def validate_asym_id(asym_id: str) -> str:
    """Valida o identificador de cadeia (mmCIF label_asym_id / auth).

    Args:
        asym_id: Cadeia.

    Returns:
        Cadeia sem espacos.

    Raises:
        AlignmentApiError: INVALID_INPUT.
    """
    text = str(asym_id or "").strip()
    if not ASYM_ID_PATTERN.match(text):
        raise AlignmentApiError(
            "Chain identifier must be 1-4 alphanumeric characters (mmCIF asym_id).",
            "INVALID_INPUT",
        )
    return text


def build_pairwise_query(
    *,
    reference_entry: str,
    target_entry: str,
    reference_chain: str,
    target_chain: str,
    method: str,
) -> dict:
    """Monta o JSON de consulta documentado pela RCSB (campo query.context).

    Args:
        reference_entry: PDB/CSM de referencia (estrutura 0; as outras sao
            superpostas a esta).
        target_entry: PDB/CSM alvo.
        reference_chain: asym_id da referencia.
        target_chain: asym_id do alvo.
        method: Nome de API (fatcat-rigid, tm-align, ...).

    Returns:
        Dict com context.mode pairwise.

    Raises:
        AlignmentApiError: INVALID_INPUT.

    Nota biologica:
        A API superpoe a segunda estrutura sobre a primeira. RMSD e TM-score
        referem-se a essa correspondencia, nao a um alinhamento de sequencia.
    """
    kind = method_kind(method)
    ref = validate_entry_id(reference_entry)
    tgt = validate_entry_id(target_entry)
    if ref == tgt and validate_asym_id(reference_chain) == validate_asym_id(target_chain):
        raise AlignmentApiError(
            "Reference and target are the same entry and chain. A self-alignment "
            "is not submitted.",
            "INVALID_INPUT",
        )
    return {
        "context": {
            "mode": MODE_PAIRWISE,
            "method": {"name": str(method).strip()},
            "structures": [
                {
                    "entry_id": ref,
                    "selection": {"asym_id": validate_asym_id(reference_chain)},
                },
                {
                    "entry_id": tgt,
                    "selection": {"asym_id": validate_asym_id(target_chain)},
                },
            ],
        },
        "helixscope_method_kind": kind,
    }


def submit_pairwise(
    *,
    reference_entry: str,
    target_entry: str,
    reference_chain: str = "A",
    target_chain: str = "A",
    method: str = METHOD_TMALIGN,
    urlopen_fn=None,
    sleep_fn=None,
) -> dict:
    """Submete o job e devolve o ticket UUID.

    Args:
        reference_entry: PDB/CSM de referencia.
        target_entry: PDB/CSM alvo.
        reference_chain: Cadeia referencia.
        target_chain: Cadeia alvo.
        method: Nome de API.
        urlopen_fn: Cliente HTTP injetavel.
        sleep_fn: time.sleep injetavel.

    Returns:
        Dict ticket, query, method, retrieved_at_utc.

    Raises:
        AlignmentApiError: categorias de rede/entrada.
    """
    del sleep_fn
    packed = build_pairwise_query(
        reference_entry=reference_entry,
        target_entry=target_entry,
        reference_chain=reference_chain,
        target_chain=target_chain,
        method=method,
    )
    query = {"context": packed["context"]}
    body = urllib.parse.urlencode({"query": json.dumps(query, separators=(",", ":"))}).encode("utf-8")
    raw = _request_bytes(
        SUBMIT_URL,
        method="POST",
        body=body,
        content_type="application/x-www-form-urlencoded",
        urlopen_fn=urlopen_fn,
    )
    try:
        payload: Any = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        payload = raw.decode("utf-8", errors="replace").strip().strip('"')
    ticket = _parse_ticket(payload)
    return {
        "ticket": ticket,
        "query": query,
        "method": str(method).strip(),
        "method_display": display_name_for(method),
        "method_kind": packed["helixscope_method_kind"],
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "endpoint": SUBMIT_URL,
    }


def fetch_results(
    ticket: str,
    *,
    urlopen_fn=None,
    sleep_fn=None,
    poll: bool = True,
) -> dict:
    """Obtem o JSON de resultados, fazendo poll enquanto RUNNING.

    Args:
        ticket: UUID devolvido por /submit.
        urlopen_fn: HTTP injetavel.
        sleep_fn: time.sleep injetavel.
        poll: Se False, uma unica leitura.

    Returns:
        Payload JSON (status COMPLETE ou o estado corrente).

    Raises:
        AlignmentApiError: TIMEOUT se o poll exceder POLL_TIMEOUT_S; ERROR se
            o servico marcar o job como ERROR.
    """
    uuid = _validate_ticket(ticket)
    sleeper = sleep_fn or time.sleep
    started = time.perf_counter()
    url = f"{RESULTS_URL}?uuid={urllib.parse.quote(uuid)}"
    while True:
        try:
            payload = _request_json(url, method="GET", urlopen_fn=urlopen_fn)
        except AlignmentApiError as exc:
            if poll and exc.category == "NOT_FOUND":
                if time.perf_counter() - started > POLL_TIMEOUT_S:
                    raise AlignmentApiError(
                        f"Alignment job {uuid} was not found after {POLL_TIMEOUT_S:.0f}s.",
                        "NOT_FOUND",
                    ) from exc
                sleeper(POLL_INTERVAL_S)
                continue
            raise
        if not isinstance(payload, Mapping):
            raise AlignmentApiError("Alignment results payload is not a JSON object.", "PARSING_ERROR")
        status = str(((payload.get("info") or {}) if isinstance(payload, Mapping) else {}).get("status") or "").upper()
        if status == "COMPLETE":
            return payload
        if status == "ERROR":
            message = str((payload.get("info") or {}).get("message") or payload.get("message") or "Alignment job ERROR.")
            raise AlignmentApiError(message, "ERROR")
        if not poll:
            return payload
        if time.perf_counter() - started > POLL_TIMEOUT_S:
            raise AlignmentApiError(
                f"Alignment job {uuid} stayed RUNNING longer than {POLL_TIMEOUT_S:.0f}s.",
                "TIMEOUT",
            )
        sleeper(POLL_INTERVAL_S)


def align_pairwise(
    *,
    reference_entry: str,
    target_entry: str,
    reference_chain: str = "A",
    target_chain: str = "A",
    method: str = METHOD_TMALIGN,
    urlopen_fn=None,
    sleep_fn=None,
) -> dict:
    """Submete e espera COMPLETE. Nao interpreta RMSD (ver structure_alignment).

    Args:
        reference_entry: PDB/CSM referencia.
        target_entry: PDB/CSM alvo.
        reference_chain: Cadeia referencia.
        target_chain: Cadeia alvo.
        method: Nome de API.
        urlopen_fn: HTTP injetavel.
        sleep_fn: sleep injetavel.

    Returns:
        Dict ticket, method, raw (JSON COMPLETE), network.

    Raises:
        AlignmentApiError: como submit/fetch.
    """
    submitted = submit_pairwise(
        reference_entry=reference_entry,
        target_entry=target_entry,
        reference_chain=reference_chain,
        target_chain=target_chain,
        method=method,
        urlopen_fn=urlopen_fn,
        sleep_fn=sleep_fn,
    )
    raw = fetch_results(
        submitted["ticket"],
        urlopen_fn=urlopen_fn,
        sleep_fn=sleep_fn,
        poll=True,
    )
    return {
        "ticket": submitted["ticket"],
        "method": submitted["method"],
        "method_display": submitted["method_display"],
        "method_kind": submitted["method_kind"],
        "query": submitted["query"],
        "raw": raw,
        "network": network_disclosure(),
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "engine": "RCSB PDB Structure Alignment API",
        "engine_location": "remote",
        "docs": DOCS_URL,
    }


def _validate_ticket(ticket: str) -> str:
    text = str(ticket or "").strip()
    if not UUID_PATTERN.match(text):
        raise AlignmentApiError(
            "Alignment ticket must be a UUID. Path-like values are refused.",
            "INVALID_INPUT",
        )
    return text


def _parse_ticket(payload: Any) -> str:
    if isinstance(payload, Mapping):
        for key in ("ticket", "uuid", "id"):
            value = payload.get(key)
            if value:
                return _validate_ticket(str(value))
        info = payload.get("info")
        if isinstance(info, Mapping) and info.get("uuid"):
            return _validate_ticket(str(info["uuid"]))
        # some responses are a bare string inside JSON
        if isinstance(payload.get("message"), str) and UUID_PATTERN.match(payload["message"].strip()):
            return _validate_ticket(payload["message"].strip())
    if isinstance(payload, str) and UUID_PATTERN.match(payload.strip()):
        return _validate_ticket(payload.strip())
    text = str(payload or "").strip().strip('"')
    if UUID_PATTERN.match(text):
        return _validate_ticket(text)
    raise AlignmentApiError(
        "Alignment API submit response did not contain a UUID ticket.",
        "PARSING_ERROR",
    )


def _request_json(url: str, *, method: str, urlopen_fn=None, body: Optional[bytes] = None, content_type: str = "") -> Any:
    raw = _request_bytes(
        url,
        method=method,
        body=body,
        content_type=content_type,
        urlopen_fn=urlopen_fn,
    )
    if not raw:
        raise AlignmentApiError("Alignment API returned an empty body.", "PARSING_ERROR")
    try:
        return json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        # submit may return a quoted UUID
        text = raw.decode("utf-8", errors="replace").strip().strip('"')
        if UUID_PATTERN.match(text):
            return text
        raise AlignmentApiError("Alignment API JSON is malformed.", "PARSING_ERROR") from exc


def _request_bytes(
    url: str,
    *,
    method: str,
    body: Optional[bytes] = None,
    content_type: str = "",
    urlopen_fn=None,
) -> bytes:
    if not url_is_allowed(url):
        raise AlignmentApiError(
            "Refusing a non-allowlisted Alignment API URL.",
            "INVALID_INPUT",
        )
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json, text/plain"}
    if content_type:
        headers["Content-Type"] = content_type
    last_error: Optional[AlignmentApiError] = None
    for attempt in range(MAX_RETRIES):
        request = Request(url, data=body, method=method, headers=headers)
        opener = urlopen_fn
        if opener is None:
            opener = build_opener(_AllowlistRedirect()).open
        try:
            with opener(request, timeout=HTTP_TIMEOUT_S) as handle:
                raw = handle.read()
                final = str(getattr(handle, "geturl", lambda: url)() or url)
        except HTTPError as exc:
            last_error = _http_error(exc)
            if last_error.category == "RATE_LIMITED" and attempt + 1 < MAX_RETRIES:
                delay = min(MAX_BACKOFF_S, BACKOFF_BASE_S ** (attempt + 1))
                if last_error.retry_after_s:
                    delay = max(delay, float(last_error.retry_after_s))
                time.sleep(delay)
                continue
            raise last_error from exc
        except AlignmentApiError:
            raise
        except (URLError, TimeoutError, socket.timeout, OSError) as exc:
            message = str(exc).lower()
            if "timed out" in message or isinstance(exc, (TimeoutError, socket.timeout)):
                raise AlignmentApiError(
                    "Alignment API request timed out.",
                    "TIMEOUT",
                ) from exc
            raise AlignmentApiError(
                f"Alignment API network error: {exc}.",
                "NETWORK_ERROR",
            ) from exc
        if final and not url_is_allowed(final):
            raise AlignmentApiError(
                "Refusing a non-allowlisted Alignment API URL after redirect.",
                "INVALID_INPUT",
            )
        if isinstance(raw, str):
            raw = raw.encode("utf-8")
        if len(raw) > MAX_RESPONSE_BYTES:
            raise AlignmentApiError(
                "Alignment API response exceeded the size cap.",
                "RESOURCE_LIMIT",
            )
        return raw
    if last_error is not None:
        raise last_error
    raise AlignmentApiError("Alignment API request failed.", "NETWORK_ERROR")


def _http_error(exc: HTTPError) -> AlignmentApiError:
    code = int(getattr(exc, "code", 0) or 0)
    retry_after = None
    headers = getattr(exc, "headers", None)
    if headers is not None:
        raw = headers.get("Retry-After")
        if raw:
            try:
                retry_after = float(raw)
            except (TypeError, ValueError):
                retry_after = None
    if code == 429:
        return AlignmentApiError(
            "RCSB Alignment API rate limited the request (HTTP 429).",
            "RATE_LIMITED",
            retry_after_s=retry_after,
        )
    if code == 404:
        return AlignmentApiError("Alignment job or endpoint was not found.", "NOT_FOUND")
    if code == 400:
        return AlignmentApiError(
            "Alignment API rejected the query (HTTP 400).",
            "INVALID_INPUT",
        )
    if code >= 500:
        return AlignmentApiError(
            f"Alignment API server error HTTP {code}.",
            "SERVICE_UNAVAILABLE",
        )
    return AlignmentApiError(f"Alignment API HTTP {code}.", "NETWORK_ERROR")
