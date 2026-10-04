"""Cliente real do Ensembl Variant Effect Predictor (VEP) via REST API.

HelixScope nao implementa predicao de consequencia. As consequencias moleculares
vem do VEP oficial (McLaren et al., Genome Biology 2016, 17:122), consultado no
endpoint REST publico do Ensembl. Os termos devolvidos sao Sequence Ontology e
sao preservados literalmente: `missense_variant` nunca e reescrito como
"mudanca de aminoacido" nem traduzido para um rotulo caseiro.

Endpoints usados (documentados em https://rest.ensembl.org):
    - GET /vep/:species/region/:region/:allele  consequencias por coordenada
    - GET /vep/:species/hgvs/:hgvs_notation     consequencias por notacao HGVS
    - GET /vep/:species/id/:id                  consequencias por rsID
    - GET /info/data                            releases disponiveis
    - GET /info/ping                            disponibilidade do servico

Assembly: rest.ensembl.org serve apenas a assembly corrente (GRCh38). Coordenadas
GRCh37 exigem o host de arquivo grch37.rest.ensembl.org. O host e escolhido pela
assembly declarada na variante e nunca por defeito silencioso.

PRIVACIDADE DE REDE: consultar o VEP remoto ENVIA a coordenada e os alelos da
variante para o servico do EMBL-EBI. `network_disclosure()` descreve exatamente o
que sai da maquina e deve ser mostrado ao utilizador antes da consulta.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import json
import socket
import time
import urllib.parse
from typing import Any, Mapping, Optional, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

from . import provenance

VEP_HOST_GRCH38: str = "rest.ensembl.org"
VEP_HOST_GRCH37: str = "grch37.rest.ensembl.org"

ALLOWED_HOSTS: frozenset[str] = frozenset({VEP_HOST_GRCH38, VEP_HOST_GRCH37})
"""Hosts Ensembl permitidos. Qualquer outro host e recusado antes do socket."""

ALLOWED_PATH_PREFIXES: tuple[str, ...] = ("/vep/", "/info/")
"""Prefixos de caminho permitidos nos hosts Ensembl."""

USER_AGENT: str = f"HelixScope/{provenance.HELIXSCOPE_VERSION} (bioinformatics; VEP REST client)"

HTTP_TIMEOUT_S: float = 60.0
MAX_RESPONSE_BYTES: int = 8 * 1024 * 1024
MAX_RETRIES: int = 3
BACKOFF_BASE_S: float = 1.5
MAX_BACKOFF_S: float = 20.0

SPECIES_HUMAN: str = "homo_sapiens"

MAX_ALLELE_LENGTH: int = 1000
"""Comprimento maximo de alelo aceito pelo endpoint de regiao.

O endpoint `/vep/:species/region` foi desenhado para variantes curtas. Eventos
maiores exigem um fluxo de variante estrutural com metodo proprio, nao uma
consulta de regiao alargada.
"""

ASSEMBLY_HOSTS: Mapping[str, str] = {
    "GRCH38": VEP_HOST_GRCH38,
    "GRCH37": VEP_HOST_GRCH37,
    "HG38": VEP_HOST_GRCH38,
    "HG19": VEP_HOST_GRCH37,
}
"""Mapa assembly declarada -> host REST. GRCh38.p14 casa pelo prefixo GRCH38."""

CONSEQUENCE_CLASSES: Mapping[str, str] = {
    "transcript_ablation": "loss_of_transcript",
    "splice_acceptor_variant": "splice_related",
    "splice_donor_variant": "splice_related",
    "stop_gained": "stop_gained",
    "frameshift_variant": "frameshift",
    "stop_lost": "stop_lost",
    "start_lost": "start_lost",
    "transcript_amplification": "amplification",
    "feature_elongation": "feature_length_change",
    "feature_truncation": "feature_length_change",
    "inframe_insertion": "inframe_indel",
    "inframe_deletion": "inframe_indel",
    "missense_variant": "missense",
    "protein_altering_variant": "protein_altering",
    "splice_donor_5th_base_variant": "splice_related",
    "splice_region_variant": "splice_related",
    "splice_donor_region_variant": "splice_related",
    "splice_polypyrimidine_tract_variant": "splice_related",
    "incomplete_terminal_codon_variant": "incomplete_terminal_codon",
    "start_retained_variant": "synonymous",
    "stop_retained_variant": "synonymous",
    "synonymous_variant": "synonymous",
    "coding_sequence_variant": "coding_unknown",
    "mature_miRNA_variant": "non_coding",
    "5_prime_UTR_variant": "UTR",
    "3_prime_UTR_variant": "UTR",
    "non_coding_transcript_exon_variant": "non_coding",
    "intron_variant": "intronic",
    "NMD_transcript_variant": "non_coding",
    "non_coding_transcript_variant": "non_coding",
    "coding_transcript_variant": "coding_unknown",
    "upstream_gene_variant": "intergenic",
    "downstream_gene_variant": "intergenic",
    "TFBS_ablation": "regulatory",
    "TFBS_amplification": "regulatory",
    "TF_binding_site_variant": "regulatory",
    "regulatory_region_ablation": "regulatory",
    "regulatory_region_amplification": "regulatory",
    "regulatory_region_variant": "regulatory",
    "intergenic_variant": "intergenic",
    "sequence_variant": "unspecified",
}
"""Agrupamento apenas VISUAL de termos Sequence Ontology.

A classe existe para ordenar a tabela da interface. O termo SO original e sempre
preservado em `consequence_terms`; a classe nunca o substitui.
"""


class VepError(RuntimeError):
    """Falha ao consultar o Ensembl VEP.

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
    """Recusa redirects para fora da allowlist Ensembl (defesa contra SSRF)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        if not url_is_allowed(newurl):
            raise VepError(
                f"Refusing an Ensembl redirect to a non-allowlisted URL: {newurl}.",
                "INVALID_INPUT",
            )
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def url_is_allowed(url: str) -> bool:
    """Valida esquema, host e prefixo de caminho contra a allowlist Ensembl.

    Args:
        url: URL absoluta.

    Returns:
        True apenas para HTTPS num host Ensembl permitido e caminho /vep/ ou /info/.

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


def host_for_assembly(assembly: str) -> str:
    """Escolhe o host REST correspondente a assembly declarada.

    Args:
        assembly: Nome da assembly (ex. GRCh38.p14, GRCh37).

    Returns:
        Hostname do servico Ensembl adequado.

    Raises:
        VepError: INVALID_INPUT quando a assembly nao e humana suportada pelo
            REST publico.

    Nota biologica:
        Nao existe defeito seguro aqui. Servir coordenadas GRCh37 ao host GRCh38
        devolveria a consequencia da base errada, e o erro seria silencioso.
    """
    text = str(assembly or "").strip().upper().replace(" ", "")
    if not text:
        raise VepError(
            "Assembly is required to choose the Ensembl REST host. GRCh38 and "
            "GRCh37 are served by different endpoints.",
            "INVALID_INPUT",
        )
    for prefix, host in ASSEMBLY_HOSTS.items():
        if text.startswith(prefix):
            return host
    raise VepError(
        f"No Ensembl REST host is configured for assembly {assembly}. Supported: "
        "GRCh38 (rest.ensembl.org) and GRCh37 (grch37.rest.ensembl.org).",
        "INVALID_INPUT",
    )


def network_disclosure(assembly: str = "GRCh38") -> dict:
    """Descreve exatamente o que sai da maquina numa consulta VEP.

    Args:
        assembly: Assembly que definira o host.

    Returns:
        Dict com host, dados enviados, dados nao enviados e politica de uso.

    Raises:
        Nenhum.

    Nota biologica:
        Coordenada e alelos de uma variante humana podem ser dados sensiveis.
        O utilizador precisa de saber que a consulta e remota antes de a fazer.
    """
    try:
        host = host_for_assembly(assembly)
    except VepError:
        host = VEP_HOST_GRCH38
    return {
        "remote": True,
        "service": "Ensembl VEP REST API (EMBL-EBI / EMBL Ensembl)",
        "host": host,
        "transport": "HTTPS",
        "data_sent": [
            "chromosome/contig name",
            "genomic position",
            "reference and alternate alleles (or the HGVS notation / rsID you submit)",
            "species (homo_sapiens)",
        ],
        "data_not_sent": [
            "any other sequence loaded in HelixScope",
            "file names or local paths",
            "reference genome files",
            "user identity beyond the HTTP User-Agent",
        ],
        "user_agent": USER_AGENT,
        "rate_limit_policy": (
            "Ensembl REST publishes X-RateLimit-* headers and returns HTTP 429 "
            "with Retry-After when exceeded. HelixScope honours Retry-After and "
            "applies bounded exponential backoff."
        ),
        "terms": "https://www.ensembl.org/info/about/legal/disclaimer.html",
        "note": (
            "This query is remote. Nothing is sent until you trigger the "
            "annotation explicitly."
        ),
    }


def service_availability(*, urlopen_fn=None, assembly: str = "GRCh38") -> dict:
    """Verifica ping e releases do Ensembl REST. Faz rede.

    Args:
        urlopen_fn: Injecao para teste; None usa o opener com allowlist.
        assembly: Assembly usada para escolher o host.

    Returns:
        Dict available, host, releases, release, checked_at_utc, reason.

    Raises:
        Nenhum. Falhas sao devolvidas como available False com categoria.

    Nota biologica:
        Servico indisponivel nao e ausencia de consequencia. A distincao evita
        que um timeout seja lido como "variante sem efeito".
    """
    host = VEP_HOST_GRCH38
    try:
        host = host_for_assembly(assembly)
    except VepError as exc:
        return {
            "available": False,
            "host": "",
            "releases": [],
            "release": None,
            "category": exc.category,
            "checked_at_utc": provenance.utc_now(),
            "reason": str(exc),
        }
    try:
        ping = _http_json(f"https://{host}/info/ping", urlopen_fn=urlopen_fn)
        data = _http_json(f"https://{host}/info/data", urlopen_fn=urlopen_fn)
    except VepError as exc:
        return {
            "available": False,
            "host": host,
            "releases": [],
            "release": None,
            "category": exc.category,
            "checked_at_utc": provenance.utc_now(),
            "reason": (
                f"Ensembl REST at {host} did not answer: {exc}. This is "
                f"{exc.category}, not an absence of consequences."
            ),
        }
    releases = [int(item) for item in (data.get("releases") or []) if str(item).isdigit()]
    alive = int((ping or {}).get("ping") or 0) == 1
    return {
        "available": alive,
        "host": host,
        "releases": releases,
        "release": max(releases) if releases else None,
        "category": "AVAILABLE" if alive else "SERVICE_UNAVAILABLE",
        "checked_at_utc": provenance.utc_now(),
        "reason": (
            f"Ensembl REST {host} responded to /info/ping with releases "
            f"{releases}."
            if alive
            else f"Ensembl REST {host} did not confirm ping."
        ),
    }


def classify_consequence(term: str) -> str:
    """Devolve a classe visual de um termo Sequence Ontology.

    Args:
        term: Termo SO exatamente como devolvido pelo VEP.

    Returns:
        Nome da classe para agrupamento na interface, ou `unclassified`.

    Raises:
        Nenhum.

    Nota biologica:
        A classe e conveniencia de apresentacao. O termo SO original e o dado
        cientifico e nunca e descartado.
    """
    return CONSEQUENCE_CLASSES.get(str(term or "").strip(), "unclassified")


def _consequence_row(payload: Mapping[str, Any], *, feature_type: str) -> dict:
    terms = [str(item) for item in (payload.get("consequence_terms") or []) if str(item).strip()]
    classes = sorted({classify_consequence(term) for term in terms})
    protein_start = payload.get("protein_start")
    protein_end = payload.get("protein_end")
    amino_acids = str(payload.get("amino_acids") or "")
    reference_aa = ""
    variant_aa = ""
    if "/" in amino_acids:
        reference_aa, _, variant_aa = amino_acids.partition("/")
    elif amino_acids:
        reference_aa = amino_acids
        variant_aa = amino_acids
    codons = str(payload.get("codons") or "")
    reference_codon = ""
    variant_codon = ""
    if "/" in codons:
        reference_codon, _, variant_codon = codons.partition("/")
    swissprot = [str(item) for item in (payload.get("swissprot") or [])]
    uniprot_isoform = [str(item) for item in (payload.get("uniprot_isoform") or [])]
    domains = [
        {
            "database": str((item or {}).get("db") or ""),
            "name": str((item or {}).get("name") or ""),
        }
        for item in (payload.get("domains") or [])
        if isinstance(item, Mapping)
    ]
    return {
        "feature_type": feature_type,
        "transcript_id": str(payload.get("transcript_id") or ""),
        "gene_id": str(payload.get("gene_id") or ""),
        "gene_symbol": str(payload.get("gene_symbol") or ""),
        "gene_symbol_source": str(payload.get("gene_symbol_source") or ""),
        "hgnc_id": str(payload.get("hgnc_id") or ""),
        "biotype": str(payload.get("biotype") or ""),
        "strand": _optional_int(payload.get("strand")),
        "canonical": bool(payload.get("canonical")),
        "consequence_terms": terms,
        "consequence_classes": classes,
        "impact": str(payload.get("impact") or ""),
        "cdna_start": _optional_int(payload.get("cdna_start")),
        "cdna_end": _optional_int(payload.get("cdna_end")),
        "cds_start": _optional_int(payload.get("cds_start")),
        "cds_end": _optional_int(payload.get("cds_end")),
        "protein_start": _optional_int(protein_start),
        "protein_end": _optional_int(protein_end),
        "protein_id": str(payload.get("protein_id") or ""),
        "amino_acids": amino_acids,
        "reference_amino_acid": reference_aa,
        "variant_amino_acid": variant_aa,
        "codons": codons,
        "reference_codon": reference_codon,
        "variant_codon": variant_codon,
        "hgvsc": str(payload.get("hgvsc") or ""),
        "hgvsp": str(payload.get("hgvsp") or ""),
        "distance": _optional_int(payload.get("distance")),
        "swissprot": swissprot,
        "uniprot_isoform": uniprot_isoform,
        "trembl": [str(item) for item in (payload.get("trembl") or [])],
        "uniparc": [str(item) for item in (payload.get("uniparc") or [])],
        "vep_domains": domains,
        "sift_score": _optional_float(payload.get("sift_score")),
        "sift_prediction": str(payload.get("sift_prediction") or ""),
        "polyphen_score": _optional_float(payload.get("polyphen_score")),
        "polyphen_prediction": str(payload.get("polyphen_prediction") or ""),
        "prediction_note": (
            "SIFT and PolyPhen-2 values, when present, are computed by those "
            "external tools and retrieved through VEP. They are PREDICTED "
            "tolerance scores, not a HelixScope pathogenicity call and not a "
            "clinical classification."
        ),
    }


def parse_vep_response(
    payload: Any,
    *,
    host: str = VEP_HOST_GRCH38,
    query: str = "",
) -> dict:
    """Converte a resposta JSON do VEP sem distorcer termos nem perder transcritos.

    Args:
        payload: JSON desserializado do endpoint VEP.
        host: Host consultado, guardado na proveniencia.
        query: URL/consulta original, guardada na proveniencia.

    Returns:
        Dict com input_variant, assembly_name, most_severe_consequence,
        transcript_consequences (TODOS), regulatory/motif consequences,
        colocated_variants e proveniencia.

    Raises:
        VepError: PARSING_ERROR se a estrutura nao for a esperada; NOT_FOUND se a
            lista vier vazia.

    Nota biologica:
        Uma variante atinge normalmente varios transcritos com consequencias
        diferentes (por exemplo missense num transcrito e intron_variant noutro).
        Todos os transcritos sao preservados: escolher um silenciosamente
        esconderia a ambiguidade biologica real.
    """
    if isinstance(payload, Mapping):
        payload = [payload]
    if not isinstance(payload, Sequence) or isinstance(payload, (str, bytes)):
        raise VepError("Ensembl VEP response is not a JSON array.", "PARSING_ERROR")
    records = [item for item in payload if isinstance(item, Mapping)]
    if not records:
        raise VepError(
            "Ensembl VEP returned no record for this query. The coordinate may "
            "be outside the assembly or the identifier may be unknown.",
            "NOT_FOUND",
        )
    record = records[0]
    transcripts = [
        _consequence_row(item, feature_type="transcript")
        for item in (record.get("transcript_consequences") or [])
        if isinstance(item, Mapping)
    ]
    regulatory = [
        {
            "feature_type": "regulatory_feature",
            "regulatory_feature_id": str((item or {}).get("regulatory_feature_id") or ""),
            "biotype": str((item or {}).get("biotype") or ""),
            "impact": str((item or {}).get("impact") or ""),
            "consequence_terms": [
                str(term) for term in ((item or {}).get("consequence_terms") or [])
            ],
        }
        for item in (record.get("regulatory_feature_consequences") or [])
        if isinstance(item, Mapping)
    ]
    intergenic = [
        {
            "feature_type": "intergenic",
            "impact": str((item or {}).get("impact") or ""),
            "consequence_terms": [
                str(term) for term in ((item or {}).get("consequence_terms") or [])
            ],
        }
        for item in (record.get("intergenic_consequences") or [])
        if isinstance(item, Mapping)
    ]
    colocated = []
    for item in record.get("colocated_variants") or []:
        if not isinstance(item, Mapping):
            continue
        colocated.append(
            {
                "id": str(item.get("id") or ""),
                "allele_string": str(item.get("allele_string") or ""),
                "start": _optional_int(item.get("start")),
                "end": _optional_int(item.get("end")),
                "clin_sig": [str(term) for term in (item.get("clin_sig") or [])],
                "clin_sig_allele": str(item.get("clin_sig_allele") or ""),
                "source": (
                    "Ensembl colocated variant record (dbSNP/ClinVar cross-reference "
                    "reported by Ensembl, not a HelixScope assertion)."
                ),
            }
        )
    canonical_count = sum(1 for row in transcripts if row.get("canonical"))
    protein_rows = [row for row in transcripts if row.get("protein_start") is not None]
    return {
        "source": "Ensembl VEP REST API",
        "evidence_status": "RETRIEVED",
        "host": str(host or ""),
        "query": str(query or ""),
        "input_variant": str(record.get("input") or ""),
        "vep_id": str(record.get("id") or ""),
        "assembly_name": str(record.get("assembly_name") or ""),
        "seq_region_name": str(record.get("seq_region_name") or ""),
        "start": _optional_int(record.get("start")),
        "end": _optional_int(record.get("end")),
        "strand": _optional_int(record.get("strand")),
        "allele_string": str(record.get("allele_string") or ""),
        "most_severe_consequence": str(record.get("most_severe_consequence") or ""),
        "transcript_consequences": transcripts,
        "regulatory_consequences": regulatory,
        "intergenic_consequences": intergenic,
        "colocated_variants": colocated,
        "transcript_count": len(transcripts),
        "canonical_transcript_count": canonical_count,
        "protein_mapped_transcript_count": len(protein_rows),
        "transcript_ambiguity": len(transcripts) > 1,
        "transcript_note": (
            f"{len(transcripts)} transcript consequences were returned and all "
            "are kept. HelixScope does not silently collapse them to a single "
            "canonical transcript."
        ),
        "consequence_vocabulary": "Sequence Ontology (terms preserved verbatim)",
        "method": (
            "Ensembl Variant Effect Predictor (McLaren et al., Genome Biology "
            "2016, 17:122) executed by the Ensembl REST service. HelixScope "
            "retrieves and displays; it does not recompute consequences."
        ),
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "extra_records": len(records) - 1,
    }


def annotate_region(
    *,
    contig: str,
    position_1based: int,
    ref: str,
    alt: str,
    assembly: str = "GRCh38",
    species: str = SPECIES_HUMAN,
    urlopen_fn=None,
    include_domains: bool = True,
) -> dict:
    """Consulta consequencias por coordenada. FAZ REDE.

    Args:
        contig: Nome de cromossoma no estilo Ensembl (ex. 17, X, MT).
        position_1based: Posicao 1-based inclusive da primeira base de REF.
        ref: Alelo de referencia (vazio para insercao pura).
        alt: Alelo alternativo.
        assembly: Assembly declarada; escolhe o host REST.
        species: Especie Ensembl.
        urlopen_fn: Injecao para teste.
        include_domains: Pede sobreposicao de dominios proteicos ao VEP.

    Returns:
        Resultado de `parse_vep_response`.

    Raises:
        VepError: INVALID_INPUT, NOT_FOUND, TIMEOUT, RATE_LIMITED,
            SERVICE_UNAVAILABLE, PARSING_ERROR ou NETWORK_ERROR.

    Nota biologica:
        Para insercao pura o VEP espera start = end + 1, marcando o ponto entre
        duas bases. Para substituicoes e delecoes o intervalo cobre as bases de
        REF.
    """
    host = host_for_assembly(assembly)
    chrom = str(contig or "").strip()
    if not chrom or any(ch in chrom for ch in "/?&#%: "):
        raise VepError("Contig token is missing or unsafe for a VEP region query.", "INVALID_INPUT")
    try:
        start = int(position_1based)
    except (TypeError, ValueError) as exc:
        raise VepError("position_1based must be an integer.", "INVALID_INPUT") from exc
    if start < 1:
        raise VepError("position_1based must be >= 1.", "INVALID_INPUT")
    reference = str(ref or "").strip().upper()
    alternate = str(alt or "").strip().upper()
    if not alternate:
        raise VepError(
            "The VEP region endpoint needs an explicit alternate allele. A pure "
            "deletion must be submitted with its reference bases and a "
            "single-base anchor, or through the HGVS endpoint.",
            "INVALID_INPUT",
        )
    if set(alternate) - set("ACGT"):
        raise VepError("Alternate allele must be A/C/G/T for a region query.", "INVALID_INPUT")
    if set(reference) - set("ACGT"):
        raise VepError(
            "Reference allele must be A/C/G/T (or empty for a pure insertion). "
            "The region end is derived from its length, so an unchecked "
            "character would silently shift the queried interval.",
            "INVALID_INPUT",
        )
    if len(reference) > MAX_ALLELE_LENGTH or len(alternate) > MAX_ALLELE_LENGTH:
        raise VepError(
            f"Alleles longer than {MAX_ALLELE_LENGTH} nt are not submitted to the "
            "region endpoint. Large events need a structural-variant workflow, "
            "not a single-nucleotide region query.",
            "RESOURCE_LIMIT",
        )
    end = start + len(reference) - 1 if reference else start - 1
    region = f"{chrom}:{start}-{end}:1"
    params = {"content-type": "application/json", "canonical": "1", "hgvs": "1", "protein": "1", "uniprot": "1"}
    if include_domains:
        params["domains"] = "1"
    url = (
        f"https://{host}/vep/{urllib.parse.quote(str(species), safe='')}/region/"
        f"{urllib.parse.quote(region, safe='')}/{urllib.parse.quote(alternate, safe='')}"
        f"?{urllib.parse.urlencode(params)}"
    )
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    return parse_vep_response(payload, host=host, query=url)


def annotate_hgvs(
    *,
    notation: str,
    assembly: str = "GRCh38",
    species: str = SPECIES_HUMAN,
    urlopen_fn=None,
    include_domains: bool = True,
) -> dict:
    """Consulta consequencias por notacao HGVS. FAZ REDE.

    Args:
        notation: Notacao HGVS completa (g., c., n., p. conforme suportado pelo VEP).
        assembly: Assembly declarada; escolhe o host REST.
        species: Especie Ensembl.
        urlopen_fn: Injecao para teste.
        include_domains: Pede dominios proteicos.

    Returns:
        Resultado de `parse_vep_response`.

    Raises:
        VepError: INVALID_INPUT para notacao vazia/insegura; restantes categorias
            conforme a resposta do servico.

    Nota biologica:
        Este endpoint existe precisamente para nao escrever um parser HGVS local:
        o mapeamento de `c.` para coordenada genomica depende do modelo oficial de
        transcritos do Ensembl.
    """
    host = host_for_assembly(assembly)
    text = str(notation or "").strip()
    if not text:
        raise VepError("HGVS notation is empty.", "INVALID_INPUT")
    if len(text) > 256:
        raise VepError("HGVS notation exceeds 256 characters.", "RESOURCE_LIMIT")
    if any(ch in text for ch in "?&#% \t\n\r/"):
        raise VepError(
            "HGVS notation contains characters that are not valid in this "
            "notation and would change the request path.",
            "INVALID_INPUT",
        )
    params = {"content-type": "application/json", "canonical": "1", "hgvs": "1", "protein": "1", "uniprot": "1"}
    if include_domains:
        params["domains"] = "1"
    url = (
        f"https://{host}/vep/{urllib.parse.quote(str(species), safe='')}/hgvs/"
        f"{urllib.parse.quote(text, safe='')}?{urllib.parse.urlencode(params)}"
    )
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    return parse_vep_response(payload, host=host, query=url)


def annotate_identifier(
    *,
    identifier: str,
    assembly: str = "GRCh38",
    species: str = SPECIES_HUMAN,
    urlopen_fn=None,
    include_domains: bool = True,
) -> dict:
    """Consulta consequencias por rsID ou identificador de variante. FAZ REDE.

    Args:
        identifier: rsID (ex. rs80357382) ou outro id conhecido do Ensembl.
        assembly: Assembly declarada; escolhe o host REST.
        species: Especie Ensembl.
        urlopen_fn: Injecao para teste.
        include_domains: Pede dominios proteicos.

    Returns:
        Resultado de `parse_vep_response`.

    Raises:
        VepError: INVALID_INPUT, NOT_FOUND e restantes categorias de rede.

    Nota biologica:
        Resolver um rsID pelo servico oficial evita inventar a coordenada: o
        mesmo rsID tem posicoes diferentes em GRCh37 e GRCh38, e alguns rsIDs
        cobrem multiplos alelos alternativos.
    """
    host = host_for_assembly(assembly)
    text = str(identifier or "").strip()
    if not text or not text.replace("_", "").replace(".", "").isalnum():
        raise VepError(
            "Variant identifier must be alphanumeric (for example rs80357382).",
            "INVALID_INPUT",
        )
    if len(text) > 64:
        raise VepError("Variant identifier exceeds 64 characters.", "RESOURCE_LIMIT")
    params = {"content-type": "application/json", "canonical": "1", "hgvs": "1", "protein": "1", "uniprot": "1"}
    if include_domains:
        params["domains"] = "1"
    url = (
        f"https://{host}/vep/{urllib.parse.quote(str(species), safe='')}/id/"
        f"{urllib.parse.quote(text, safe='')}?{urllib.parse.urlencode(params)}"
    )
    payload = _http_json(url, urlopen_fn=urlopen_fn)
    return parse_vep_response(payload, host=host, query=url)


def _http_json(url: str, *, urlopen_fn=None) -> Any:
    if not url_is_allowed(url):
        raise VepError(
            f"Refusing a non-allowlisted Ensembl URL: {url}.", "INVALID_INPUT"
        )
    request = Request(
        url,
        method="GET",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    opener = urlopen_fn or build_opener(_AllowlistRedirect()).open
    attempt = 0
    last: Optional[VepError] = None
    while attempt < MAX_RETRIES:
        attempt += 1
        try:
            with opener(request, timeout=HTTP_TIMEOUT_S) as handle:
                raw = handle.read()
                final = str(getattr(handle, "geturl", lambda: url)() or url)
            break
        except HTTPError as exc:
            error = _http_error(exc)
            last = error
            if error.category in {"RATE_LIMITED", "SERVICE_UNAVAILABLE"} and attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, error.retry_after_s))
                continue
            raise error from exc
        except VepError:
            raise
        except (URLError, TimeoutError, socket.timeout, OSError) as exc:
            message = str(exc).lower()
            timed_out = isinstance(exc, (TimeoutError, socket.timeout)) or "timed out" in message
            error = VepError(
                (
                    "Ensembl VEP request timed out. This is TIMEOUT, not an "
                    "absence of consequences."
                )
                if timed_out
                else (
                    f"Ensembl VEP network error: {exc}. This is NETWORK_ERROR, "
                    "not NOT_FOUND."
                ),
                "TIMEOUT" if timed_out else "NETWORK_ERROR",
            )
            last = error
            if attempt < MAX_RETRIES:
                time.sleep(_backoff_seconds(attempt, None))
                continue
            raise error from exc
    else:  # pragma: no cover - defensive; loop always breaks or raises
        raise last or VepError("Ensembl VEP request failed.", "SERVICE_UNAVAILABLE")
    if final and not url_is_allowed(final):
        raise VepError(
            f"Ensembl VEP answered from a non-allowlisted URL: {final}.",
            "INVALID_INPUT",
        )
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    if len(raw) > MAX_RESPONSE_BYTES:
        raise VepError(
            f"Ensembl VEP response exceeds {MAX_RESPONSE_BYTES:,} bytes.",
            "RESOURCE_LIMIT",
        )
    text = raw.decode("utf-8", errors="replace")
    if not text.strip():
        raise VepError("Ensembl VEP response is empty.", "PARSING_ERROR")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise VepError("Ensembl VEP response is not valid JSON.", "PARSING_ERROR") from exc


def _backoff_seconds(attempt: int, retry_after_s: Optional[float]) -> float:
    if retry_after_s is not None and retry_after_s >= 0:
        return min(float(retry_after_s), MAX_BACKOFF_S)
    return min(BACKOFF_BASE_S ** attempt, MAX_BACKOFF_S)


def _http_error(exc: HTTPError) -> VepError:
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
        return VepError(
            "Ensembl VEP has no record for this variant query (HTTP 404).",
            "NOT_FOUND",
        )
    if code == 429:
        return VepError(
            "Ensembl VEP rate limited the request (HTTP 429). Respect the "
            "published rate limit before retrying.",
            "RATE_LIMITED",
            retry_after,
        )
    if code == 400:
        return VepError(
            f"Ensembl VEP rejected the request as invalid (HTTP 400): {exc}.",
            "INVALID_INPUT",
        )
    if code >= 500:
        return VepError(
            f"Ensembl VEP service unavailable (HTTP {code}).",
            "SERVICE_UNAVAILABLE",
            retry_after,
        )
    return VepError(f"Ensembl VEP HTTP {code}: {exc}.", "SERVICE_UNAVAILABLE", retry_after)


def _optional_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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
