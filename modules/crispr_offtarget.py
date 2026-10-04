"""Busca de off-targets CRISPR verificavel em referencia fornecida.

Algoritmo: indice de sitios PAM NGG e NAG (atividade reduzida, Hsu 2013) em
ambas as fitas de cada contig, comparacao do spacer de 20 nt, contagem real de
mismatches e extracao da sequencia a partir do contig. Sem bulges. Sem download
de assembly. Sem rotulo genome-wide.

Nenhuma funcao importa Streamlit. Nenhuma funcao aceita URL ou caminho de
executavel do usuario.
"""

from __future__ import annotations

import hashlib
import time
from typing import Any, Mapping, Optional, Sequence

from . import crispr, crispr_cfd, crispr_reference, crispr_specificity, provenance, scientific_checks

JOB_NOT_RUN: str = "NOT_RUN"
JOB_RUNNING: str = "RUNNING"
JOB_COMPLETED: str = "COMPLETED"
JOB_COMPLETED_FULL_REFERENCE: str = "COMPLETED_FULL_REFERENCE"
JOB_COMPLETED_TEST_REFERENCE: str = "COMPLETED_TEST_REFERENCE"
JOB_FAILED: str = "ERROR"
JOB_TIMEOUT: str = "TIMEOUT"
JOB_RESOURCE_LIMIT: str = "RESOURCE_LIMIT"
JOB_INVALID: str = "INVALID_INPUT"
JOB_INCOMPLETE: str = "INCOMPLETE"
JOB_UNAVAILABLE: str = "UNAVAILABLE"
JOB_QUEUED: str = "QUEUED"
JOB_CANCELLED: str = "CANCELLED"

COMPLETE_SEARCH_STATUSES: frozenset[str] = frozenset(
    {JOB_COMPLETED, JOB_COMPLETED_FULL_REFERENCE, JOB_COMPLETED_TEST_REFERENCE}
)

PAM_MODE_NGG_NAG: str = "ngg_or_nag"
"""PAMs aceitos na busca SpCas9: NGG (primario) e NAG (atividade reduzida)."""

SEED_LENGTH: int = 12
"""Semente PAM-proximal: 12 nt no 3' do spacer de 20 nt (indices 8:20).
Valor alinhado a crispr.SPCAS9_SEED_NT; literal para o import nao depender
de um modules.crispr ainda a meio do reload do Streamlit."""

DEFAULT_SEARCH_TIMEOUT_S: float = 20.0
"""Timeout padrao da busca, alinhado a crispr.REFERENCE_SEARCH_TIMEOUT_S."""

DEFAULT_MAX_HITS: int = 200
"""Teto padrao de hits, alinhado a crispr.MAX_REFERENCE_HITS."""


def build_pam_index(reference: Mapping[str, Any], *, cas_system: str = "SpCas9") -> dict:
    """Lista todos os sitios NGG/NAG com spacer de 20 nt em cada contig.

    Args:
        reference: Saida de crispr_reference.load_reference_from_text.
        cas_system: Apenas SpCas9 e suportado neste motor.

    Returns:
        Dict com sites (lista), index_hash, n_sites, method e tempo de construcao.

    Raises:
        crispr.CrisprError: INVALID_INPUT se o sistema nao for SpCas9.

    Nota biologica:
        O indice e a lista de PAMs reais extraidos da referencia. Nao e Bowtie
        nem Cas-OFFinder. A identidade do indice entra na proveniencia.
    """
    if cas_system != "SpCas9":
        raise crispr.CrisprError(
            "Provided-reference off-target search is implemented for SpCas9 only.",
            "UNAVAILABLE",
        )
    started = time.perf_counter()
    spacer_len = crispr.GUIDE_LENGTH_BY_CAS["SpCas9"]
    sites: list[dict] = []
    for contig in list(reference.get("contigs") or []):
        contig_id = str(contig.get("identifier") or "")
        seq = str(contig.get("sequence") or "")
        n = len(seq)
        for strand, strand_seq in (("+", seq), ("-", crispr._reverse_complement(seq))):
            limit = n - spacer_len - 3 + 1
            for i in range(0, max(0, limit)):
                pam = strand_seq[i + spacer_len : i + spacer_len + 3]
                if crispr._iupac_match(pam, "NGG"):
                    pam_class = "NGG"
                elif crispr._iupac_match(pam, "NAG"):
                    pam_class = "NAG_reduced_activity"
                else:
                    continue
                spacer = strand_seq[i : i + spacer_len]
                if set(spacer) - {"A", "C", "G", "T"}:
                    continue
                if strand == "+":
                    start_0based = i
                    end_0based = i + spacer_len
                    pam_start = i + spacer_len
                    pam_end = i + spacer_len + 3
                    position = i + 1
                else:
                    start_0based = n - (i + spacer_len)
                    end_0based = n - i
                    pam_start = start_0based - 3
                    pam_end = start_0based
                    position = start_0based + 1
                sites.append(
                    {
                        "contig_id": contig_id,
                        "strand": strand,
                        "position_1based": position,
                        "start_0based": start_0based,
                        "end_0based": end_0based,
                        "pam_start_0based": pam_start,
                        "pam_end_0based": pam_end,
                        "spacer": spacer,
                        "pam_sequence": pam,
                        "pam_class": pam_class,
                    }
                )
    canonical = "|".join(
        f"{s['contig_id']}:{s['strand']}:{s['start_0based']}:{s['pam_sequence']}"
        for s in sites
    )
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    return {
        "method": "pam_site_list_ngg_nag",
        "cas_system": cas_system,
        "spacer_nt": spacer_len,
        "seed_nt": SEED_LENGTH,
        "seed_orientation": "PAM-proximal 3' of spacer (positions 9-20 of 20 nt)",
        "bulges_supported": False,
        "n_sites": len(sites),
        "sites": sites,
        "index_hash": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "build_ms": round(elapsed_ms, 3),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def search_off_targets_in_reference(
    guide: Mapping[str, Any],
    reference: Mapping[str, Any],
    *,
    max_mismatches: int = 3,
    include_perfect: bool = True,
    timeout_s: float = DEFAULT_SEARCH_TIMEOUT_S,
    max_hits: int = DEFAULT_MAX_HITS,
    target_hash: str = "",
    index: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Busca off-targets verificados na referencia fornecida.

    Cada hit tem sequencia extraida do contig, PAM presente no locus, fita,
    coordenadas internas do contig, mismatches e posicoes calculados, score Hsu
    2013 por sitio quando o spacer tem 20 nt, e proveniencia. Zero hits apos
    busca completa e um resultado valido; NOT_RUN/ERROR/TIMEOUT nunca viram 0.

    Args:
        guide: Guia SpCas9 (spacer 20 nt + PAM NGG).
        reference: Referencia validada.
        max_mismatches: Maximo de mismatches no spacer (0 a 4 nesta fase).
        include_perfect: Inclui sitios 0-mismatch (padrao True na referencia).
        timeout_s: Timeout da comparacao, alem da construcao do indice.
        max_hits: Teto de hits; excesso devolve RESOURCE_LIMIT e truncated.
        target_hash: Hash da sequencia de desenho, para a chave de cache.
        index: Indice precomputado; se None, e construido agora.

    Returns:
        Envelope com status, hits, escopo, hashes, tempos e truncated.

    Raises:
        Nenhum. Falhas viram status ERROR/INVALID_INPUT/TIMEOUT/RESOURCE_LIMIT
        no envelope. Nunca uma lista vazia fingindo busca completa nesses casos.

    Nota biologica:
        Escopo = contigs do FASTA carregado. Nao e genome-wide. BLAST nao e
        especificidade CRISPR. Sobreposicao com uma feature NCBI nao e
        'gene disrupted'.
    """
    started = time.perf_counter()
    started_mono = time.monotonic()
    validation = crispr.validate_spcas9_guide(guide)
    if not validation["valid"]:
        return _empty_envelope(
            status=JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="; ".join(validation["errors"]),
            elapsed_ms=_elapsed_ms(started),
        )
    try:
        max_mm = int(max_mismatches)
    except (TypeError, ValueError):
        return _empty_envelope(
            status=JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="max_mismatches must be an integer.",
            elapsed_ms=_elapsed_ms(started),
        )
    if max_mm < 0 or max_mm > 4:
        return _empty_envelope(
            status=JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            reason="max_mismatches must be between 0 and 4 for this method.",
            elapsed_ms=_elapsed_ms(started),
        )
    try:
        hit_cap = int(max_hits)
        time_cap = float(timeout_s)
    except (TypeError, ValueError):
        return _empty_envelope(
            status=JOB_INVALID,
            guide=guide,
            reference=reference,
            max_mismatches=max_mm,
            include_perfect=include_perfect,
            reason="timeout_s and max_hits must be numeric.",
            elapsed_ms=_elapsed_ms(started),
        )
    if hit_cap < 1:
        return _empty_envelope(
            status=JOB_RESOURCE_LIMIT,
            guide=guide,
            reference=reference,
            max_mismatches=max_mm,
            include_perfect=include_perfect,
            reason="max_hits must be at least 1.",
            elapsed_ms=_elapsed_ms(started),
        )
    spacer = validation["guide_sequence"]
    try:
        pam_index = dict(index) if index is not None else build_pam_index(reference)
    except crispr.CrisprError as exc:
        return _empty_envelope(
            status=exc.category if exc.category else JOB_FAILED,
            guide=guide,
            reference=reference,
            max_mismatches=max_mm,
            include_perfect=include_perfect,
            reason=str(exc),
            elapsed_ms=_elapsed_ms(started),
        )

    hits: list[dict] = []
    truncated = False
    unverified = 0
    status = JOB_COMPLETED
    reason = ""
    sites = list(pam_index.get("sites") or [])
    search_started = time.perf_counter()
    for site_index, site in enumerate(sites):
        if time.monotonic() - started_mono > time_cap:
            status = JOB_TIMEOUT
            reason = (
                f"Off-target search exceeded {time_cap:.1f}s. Incomplete; "
                "not converted to zero hits."
            )
            truncated = True
            hits = []
            break
        mismatches, mismatch_positions = _mismatch_positions(spacer, str(site["spacer"]))
        if mismatches > max_mm:
            continue
        if mismatches == 0 and not include_perfect:
            continue
        try:
            hit = _verified_hit(
                guide=guide,
                spacer=spacer,
                site=site,
                reference=reference,
                mismatches=mismatches,
                mismatch_positions=mismatch_positions,
            )
        except (crispr.CrisprError, crispr_reference.ReferenceError):
            unverified += 1
            continue
        hits.append(hit)
        if len(hits) >= hit_cap:
            if site_index + 1 < len(sites):
                status = JOB_RESOURCE_LIMIT
                reason = (
                    f"Hit cap of {hit_cap} reached before the PAM index was "
                    "fully compared. Result is truncated and is not a complete "
                    "off-target catalogue."
                )
                truncated = True
            break

    search_ms = (time.perf_counter() - search_started) * 1000.0
    if status == JOB_COMPLETED and unverified > 0 and not hits:
        status = JOB_FAILED
        reason = (
            f"{unverified} PAM-index candidate(s) failed contig re-extraction. "
            "Not reported as zero off-targets."
        )
    elif status == JOB_COMPLETED and not hits:
        reason = (
            "0 verified hits under selected method/reference. This is not a "
            "safety claim and not 'no risk'."
        )
    envelope = build_search_envelope(
        status=status,
        reason=reason,
        guide=guide,
        reference=reference,
        hits=hits,
        max_mm=max_mm,
        include_perfect=include_perfect,
        truncated=truncated,
        hit_cap=hit_cap,
        time_cap=time_cap,
        target_hash=target_hash,
        elapsed_ms=_elapsed_ms(started),
        search_ms=round(search_ms, 3),
        index_hash=str(pam_index.get("index_hash") or ""),
        index_n_sites=pam_index.get("n_sites"),
        index_build_ms=pam_index.get("build_ms"),
        algorithm=crispr.ALGORITHM_SPCAS9_PAM_INDEX,
    )
    envelope["unverified_count"] = unverified
    return envelope


def annotate_verified_hits(
    hits: Sequence[Mapping[str, Any]],
    guide: Mapping[str, Any],
) -> list[dict]:
    """Acrescenta CFD por sitio, Hsu ponderado CRISPOR e site_role.

    Args:
        hits: Hits ja verificados contra o contig.
        guide: Guia com pam_sequence NGG.

    Returns:
        Lista anotada. Hsu por sitio permanece em score.

    Raises:
        Nenhum.
    """
    pam = str(guide.get("pam_sequence") or "")
    annotated: list[dict] = []
    for hit in hits:
        item = crispr_cfd.attach_cfd_to_hit(hit, pam)
        weighted = crispr_specificity.mit_weighted_hsu(item)
        item["mit_weighted_hsu"] = weighted
        item["mit_weighted_hsu_method"] = (
            "Hsu 2013 single-hit; CRISPOR multiplies by 0.2 when the PAM "
            "dinucleotide is not GG. CFD is not reweighted this way."
        )
        annotated.append(item)
    return crispr_specificity.assign_site_roles(annotated)


def build_search_envelope(
    *,
    status: str,
    reason: str,
    guide: Mapping[str, Any],
    reference: Mapping[str, Any],
    hits: Sequence[Mapping[str, Any]],
    max_mm: int,
    include_perfect: bool,
    truncated: bool,
    hit_cap: int,
    time_cap: float,
    target_hash: str,
    elapsed_ms: float,
    search_ms: float,
    index_hash: str,
    index_n_sites: Any,
    index_build_ms: Any,
    algorithm: str,
    tool: Optional[Mapping[str, Any]] = None,
    genome_wide: bool = False,
) -> dict:
    """Monta o envelope de busca com CFD, papeis de sitio e agregados.

    Args:
        status: COMPLETED, TIMEOUT, RESOURCE_LIMIT, ...
        reason: Texto honesto.
        guide: Guia pesquisado.
        reference: Referencia pesquisada.
        hits: Hits verificados (ignorados se TIMEOUT/ERROR).
        max_mm: Max mismatches.
        include_perfect: Se 0-mismatch entrou.
        truncated: True se incompleto.
        hit_cap: Teto.
        time_cap: Timeout s.
        target_hash: Hash do DNA de desenho.
        elapsed_ms: Tempo total.
        search_ms: Tempo da comparacao.
        index_hash: Hash do indice PAM, se houver.
        index_n_sites: N sitios indexados ou None.
        index_build_ms: Tempo do indice ou None.
        algorithm: Nome do motor.
        tool: Metadados de ferramenta externa, se houver.
        genome_wide: True so para COMPLETED_FULL_REFERENCE de assembly READY.

    Returns:
        Envelope JSON-safe. genome_wide so e True quando status e
        COMPLETED_FULL_REFERENCE e o parametro genome_wide e True.
        Agregados MIT/CFD so em buscas completas nao truncadas.

    Raises:
        Nenhum.
    """
    spacer = str(guide.get("guide_sequence") or "")
    ref_hash = str(reference.get("identity_hash") or "")
    assembly = str(reference.get("assembly_declared") or "")
    show_hits = status in {
        JOB_COMPLETED,
        JOB_COMPLETED_FULL_REFERENCE,
        JOB_COMPLETED_TEST_REFERENCE,
        JOB_RESOURCE_LIMIT,
    }
    raw_hits = [dict(hit) for hit in hits] if show_hits else []
    annotated = annotate_verified_hits(raw_hits, guide) if raw_hits else []
    result_id = crispr.offtarget_cache_key(
        guide_sequence=spacer,
        target_hash=target_hash,
        reference_hash=ref_hash,
        assembly_declared=assembly,
        algorithm=algorithm,
        model="Hsu2013_single_hit+Doench2016_CFD",
        model_version=(
            "CRISPOR-hitScoreM-20nt+" + crispr_cfd.CFD_VERSION
        ),
        max_mismatches=max_mm,
        include_perfect=include_perfect,
        pam_mode=PAM_MODE_NGG_NAG,
    )
    contig_ids = [
        str(contig.get("identifier") or "")
        for contig in list(reference.get("contigs") or [])
    ]
    envelope = {
        "status": status,
        "reason": reason,
        "algorithm": algorithm,
        "algorithm_disclosure": {
            "method": (
                "Cas-OFFinder subprocess"
                if algorithm != crispr.ALGORITHM_SPCAS9_PAM_INDEX
                else "PAM-indexed exact spacer comparison"
            ),
            "search_scope": reference.get("scope"),
            "mismatch_model": "Hamming distance on 20 nt spacer, positions stored",
            "pam_model": "NGG required as primary; NAG labelled reduced-activity",
            "seed_model": (
                f"{SEED_LENGTH} nt PAM-proximal (3' of spacer, positions 9-20)"
            ),
            "bulge_support": False,
            "score_model": (
                crispr.HSU_2013_HIT_METHOD + "; " + crispr_cfd.CFD_METHOD
            ),
        },
        "genome_wide": bool(genome_wide)
        and status == JOB_COMPLETED_FULL_REFERENCE
        and not truncated,
        "scope": reference.get("scope"),
        "scope_label": reference.get("scope_label"),
        "contigs_searched": contig_ids,
        "regions": "full contigs in the provided FASTA",
        "guide_sequence": spacer,
        "guide_id": guide.get("guide_id") or crispr._guide_id(guide),
        "guide_hash": provenance.sequence_digest(spacer),
        "pam_sequence": guide.get("pam_sequence"),
        "strand": guide.get("strand"),
        "reference_identity_hash": ref_hash,
        "organism_declared": reference.get("organism_declared") or "",
        "assembly_declared": assembly,
        "verified_assembly": reference.get("verified_assembly") or "",
        "assembly_identity_status": reference.get("assembly_identity_status") or "",
        "version_declared": reference.get("version_declared") or "",
        "accession": reference.get("accession") or "",
        "index_hash": index_hash,
        "index_n_sites": index_n_sites,
        "index_build_ms": index_build_ms,
        "search_ms": search_ms,
        "elapsed_ms": elapsed_ms,
        "max_mismatches": max_mm,
        "include_perfect": include_perfect,
        "pam_mode": PAM_MODE_NGG_NAG,
        "truncated": truncated,
        "hit_cap": hit_cap,
        "timeout_s": time_cap,
        "verified_hit_count": (
            len(annotated) if status in COMPLETE_SEARCH_STATUSES else None
        ),
        "hits": annotated,
        "result_identity": result_id,
        "execution": {
            "mode": "in_process",
            "distributed_workers": False,
            "authentication": False,
            "per_user_quotas": False,
        },
        "software_version": provenance.HELIXSCOPE_VERSION,
        "timestamp": provenance.utc_now(),
        "kind": "COMPUTED" if status in COMPLETE_SEARCH_STATUSES else "ERROR",
        "cfd_model": crispr_cfd.cfd_model_record(),
    }
    if tool is not None:
        envelope["tool"] = dict(tool)
    if status == JOB_TIMEOUT:
        envelope["hits"] = []
        envelope["verified_hit_count"] = None
        envelope["truncated"] = True
        envelope["kind"] = "ERROR"
    if status == JOB_RESOURCE_LIMIT:
        envelope["verified_hit_count"] = None
        envelope["kind"] = "ERROR"
    if status == JOB_UNAVAILABLE:
        envelope["hits"] = []
        envelope["verified_hit_count"] = None
        envelope["kind"] = "UNAVAILABLE"
    draft_for_scores = dict(envelope)
    envelope["guide_specificity"] = crispr_specificity.scores_from_search(draft_for_scores)
    return provenance.json_safe(envelope)


def validate_hit_against_reference(
    hit: Mapping[str, Any], reference: Mapping[str, Any]
) -> bool:
    """Confere se sequencia, PAM, coordenadas e mismatches batem no contig.

    Args:
        hit: Off-target candidato.
        reference: Referencia da qual o hit afirma ter sido extraido.

    Returns:
        True somente se a janela extraida do contig for a sequencia reportada,
        o PAM existir, a fita for valida e o numero de mismatches for o
        calculado. False nao e corrigido em silencio.

    Raises:
        Nenhum.
    """
    try:
        contig = crispr_reference.contig_by_id(
            reference, str(hit.get("chromosome_or_contig") or "")
        )
    except crispr_reference.ReferenceError:
        return False
    seq = str(contig.get("sequence") or "")
    try:
        start = int(hit.get("start_0based"))
        end = int(hit.get("end_0based"))
        pam_start = int(hit.get("pam_start_0based"))
        pam_end = int(hit.get("pam_end_0based"))
    except (TypeError, ValueError):
        return False
    strand = str(hit.get("strand") or "")
    if strand not in {"+", "-"}:
        return False
    if not scientific_checks.span_is_valid(start, end, len(seq), base=0):
        return False
    if not scientific_checks.span_is_valid(pam_start, pam_end, len(seq), base=0):
        return False
    sense_spacer = seq[start:end]
    sense_pam = seq[pam_start:pam_end]
    reported_spacer = str(hit.get("target_sequence") or "")
    reported_pam = str(hit.get("PAM") or hit.get("pam_sequence") or "")
    if strand == "+":
        if sense_spacer != reported_spacer or sense_pam != reported_pam:
            return False
    else:
        if crispr._reverse_complement(sense_spacer) != reported_spacer:
            return False
        if crispr._reverse_complement(sense_pam) != reported_pam:
            return False
    guide_seq = str(hit.get("guide_sequence") or "")
    mismatches, positions = _mismatch_positions(guide_seq, reported_spacer)
    try:
        reported_n = int(hit.get("mismatches"))
    except (TypeError, ValueError):
        return False
    if mismatches != reported_n:
        return False
    reported_pos = list(hit.get("mismatch_positions") or [])
    if reported_pos != positions:
        return False
    ref_hash = str(reference.get("identity_hash") or "")
    hit_hash = str(hit.get("reference_identity_hash") or "")
    if ref_hash and hit_hash and ref_hash != hit_hash:
        return False
    return True


def overlapping_ncbi_features(
    hit: Mapping[str, Any],
    features: Sequence[Mapping[str, Any]],
    *,
    contig_length: int,
) -> list[dict]:
    """Lista features NCBI cujo intervalo 0-based sobrepoe o hit.

    Args:
        hit: Hit com start_0based e end_0based no contig.
        features: Saida de ncbi_fetch.extract_annotated_features.
        contig_length: Comprimento do contig (para recusar coordenadas invalidas).

    Returns:
        Lista de {type, gene, product, start, end, overlap_nt}. Nunca afirma
        que o gene foi disrupted. Features UNAVAILABLE sao ignoradas.

    Raises:
        Nenhum.

    Nota biologica:
        Distancia e overlap usam intervalos 0-based semiabertos no mesmo contig.
        overlap_nt = comprimento da intersecao. Sem mapping validado, a lista
        fica vazia.
    """
    try:
        h_start = int(hit.get("start_0based"))
        h_end = int(hit.get("end_0based"))
    except (TypeError, ValueError):
        return []
    if not scientific_checks.span_is_valid(h_start, h_end, contig_length, base=0):
        return []
    found: list[dict] = []
    for feature in features:
        if str(feature.get("status") or "") != "AVAILABLE":
            continue
        try:
            f_start = int(feature.get("start"))
            f_end = int(feature.get("end"))
        except (TypeError, ValueError):
            continue
        if not scientific_checks.span_is_valid(f_start, f_end, contig_length, base=0):
            continue
        overlap = min(h_end, f_end) - max(h_start, f_start)
        if overlap <= 0:
            continue
        found.append(
            {
                "type": str(feature.get("type") or ""),
                "gene": str(feature.get("gene") or ""),
                "product": str(feature.get("product") or ""),
                "start": f_start,
                "end": f_end,
                "overlap_nt": overlap,
                "claim": (
                    "Interval overlap on this NCBI record. Not a demonstration "
                    "that the gene is disrupted."
                ),
            }
        )
    return found


def export_offtarget_rows(search: Mapping[str, Any]) -> list[dict]:
    """Linhas CSV/JSON de hits. Recusa exportar 0 hits quando a busca nao completou.

    Args:
        search: Envelope de search_off_targets_in_reference.

    Returns:
        Lista de dicts. Vazia quando status nao e COMPLETED (incluindo
        RESOURCE_LIMIT, TIMEOUT, NOT_RUN, ERROR). COMPLETED com zero hits
        devolve uma lista vazia acompanhada no envelope, nao uma linha '0 hits'
        para estados de falha.

    Raises:
        crispr.CrisprError: Se o chamador pedir exportacao de um estado que nao
            pode ser representado como catalogo de hits.
    """
    status = str(search.get("status") or JOB_NOT_RUN)
    if status in {JOB_NOT_RUN, JOB_FAILED, JOB_TIMEOUT, JOB_INCOMPLETE, JOB_INVALID, JOB_CANCELLED}:
        raise crispr.CrisprError(
            f"Off-target search status is {status}; hit rows are not exported "
            "as zero hits.",
            status if status in {"TIMEOUT", "INVALID_INPUT", "ERROR"} else "ERROR",
        )
    if status == JOB_RESOURCE_LIMIT:
        raise crispr.CrisprError(
            "Search hit the resource limit. A truncated list is not exported "
            "as a complete catalogue. Inspect the in-session result with "
            "truncated=true.",
            "RESOURCE_LIMIT",
        )
    if status not in COMPLETE_SEARCH_STATUSES:
        raise crispr.CrisprError(
            f"Cannot export hits for status {status}.",
            "ERROR",
        )
    rows: list[dict] = []
    for hit in list(search.get("hits") or []):
        rows.append(
            {
                "Guide": hit.get("guide_sequence"),
                "Chromosome_or_contig": hit.get("chromosome_or_contig"),
                "Position_1based": hit.get("position"),
                "Start_0based": hit.get("start_0based"),
                "End_0based": hit.get("end_0based"),
                "Strand": hit.get("strand"),
                "Target": hit.get("target_sequence"),
                "PAM": hit.get("PAM"),
                "PAM_class": hit.get("pam_class"),
                "Mismatches": hit.get("mismatches"),
                "Mismatch_positions": ",".join(
                    str(p) for p in list(hit.get("mismatch_positions") or [])
                ),
                "Hsu_hit_score": hit.get("score"),
                "Score_method": hit.get("score_method"),
                "Score_model": hit.get("score_model"),
                "Score_version": hit.get("score_version"),
                "Score_scale": "0-100",
                "CFD_score": hit.get("cfd_score"),
                "CFD_method": hit.get("cfd_method"),
                "CFD_model": hit.get("cfd_model"),
                "CFD_version": hit.get("cfd_version"),
                "CFD_scale": hit.get("cfd_scale"),
                "Site_role": hit.get("site_role"),
                "Reference_hash": hit.get("reference_identity_hash"),
                "Assembly_declared": hit.get("assembly_declared"),
                "Organism_declared": hit.get("organism_declared"),
                "Search_scope": search.get("scope"),
                "Algorithm": search.get("algorithm"),
                "Status": status,
                "Genome_wide": bool(search.get("genome_wide")),
            }
        )
    return rows


def not_run_envelope() -> dict:
    """Envelope NOT_RUN. verified_hit_count e None, nunca 0."""
    return provenance.json_safe(
        {
            "status": JOB_NOT_RUN,
            "reason": "Off-target search was not run.",
            "genome_wide": False,
            "verified_hit_count": None,
            "hits": [],
            "truncated": False,
            "kind": "UNAVAILABLE",
            "execution": {
                "mode": "in_process",
                "distributed_workers": False,
                "authentication": False,
                "per_user_quotas": False,
            },
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
    )


def _mismatch_positions(guide: str, window: str) -> tuple[int, list[int]]:
    """Hamming 1-based nas posicoes que diferem. Nao estima."""
    positions: list[int] = []
    for index, (a, b) in enumerate(zip(guide, window), start=1):
        if a != b:
            positions.append(index)
    extra = abs(len(guide) - len(window))
    if extra:
        raise crispr.CrisprError(
            "Guide and window lengths differ; mismatch count refused.",
            "INVALID_INPUT",
        )
    return len(positions), positions


def _verified_hit(
    *,
    guide: Mapping[str, Any],
    spacer: str,
    site: Mapping[str, Any],
    reference: Mapping[str, Any],
    mismatches: int,
    mismatch_positions: list[int],
) -> dict:
    """Monta um hit so depois de reextrair spacer e PAM do contig."""
    contig = crispr_reference.contig_by_id(reference, str(site["contig_id"]))
    seq = str(contig["sequence"])
    start = int(site["start_0based"])
    end = int(site["end_0based"])
    pam_start = int(site["pam_start_0based"])
    pam_end = int(site["pam_end_0based"])
    strand = str(site["strand"])
    sense_spacer = seq[start:end]
    sense_pam = seq[pam_start:pam_end]
    if strand == "+":
        target_seq = sense_spacer
        pam_seq = sense_pam
    else:
        target_seq = crispr._reverse_complement(sense_spacer)
        pam_seq = crispr._reverse_complement(sense_pam)
    if target_seq != str(site["spacer"]) or pam_seq != str(site["pam_sequence"]):
        raise crispr.CrisprError(
            "Index spacer/PAM does not match contig extraction.",
            "PARSING_ERROR",
        )
    if not (
        crispr._iupac_match(pam_seq, "NGG") or crispr._iupac_match(pam_seq, "NAG")
    ):
        raise crispr.CrisprError("Extracted PAM is not NGG or NAG.", "PARSING_ERROR")
    hsu = crispr.hsu_single_hit_score(spacer, target_seq)
    seed = spacer[-SEED_LENGTH:]
    seed_window = target_seq[-SEED_LENGTH:]
    seed_mismatches = sum(a != b for a, b in zip(seed, seed_window))
    hit = {
        "guide_sequence": spacer,
        "guide_hash": provenance.sequence_digest(spacer),
        "reference_identity_hash": reference.get("identity_hash"),
        "organism_declared": reference.get("organism_declared") or "",
        "assembly_declared": reference.get("assembly_declared") or "",
        "version_declared": reference.get("version_declared") or "",
        "chromosome_or_contig": contig["identifier"],
        "position": int(site["position_1based"]),
        "start_0based": start,
        "end_0based": end,
        "pam_start_0based": pam_start,
        "pam_end_0based": pam_end,
        "coordinate_system": "contig_sense_internal",
        "strand": strand,
        "target_sequence": target_seq,
        "PAM": pam_seq,
        "pam_sequence": pam_seq,
        "pam_class": site["pam_class"],
        "mismatches": mismatches,
        "mismatch_positions": mismatch_positions,
        "seed_nt": SEED_LENGTH,
        "seed_mismatches": seed_mismatches,
        "score": hsu["score"],
        "score_method": hsu["method"],
        "score_status": hsu["status"],
        "score_model": hsu["model"],
        "score_version": hsu["version"],
        "status": "VERIFIED",
        "bulges_supported": False,
        "site_class": "perfect_site" if mismatches == 0 else "mismatch_site",
        "provenance": {
            "timestamp": provenance.utc_now(),
            "algorithm": crispr.ALGORITHM_SPCAS9_PAM_INDEX,
            "software_version": provenance.HELIXSCOPE_VERSION,
        },
    }
    if not validate_hit_against_reference(hit, reference):
        raise crispr.CrisprError(
            "Hit failed post-extraction validation.",
            "PARSING_ERROR",
        )
    return hit


def _empty_envelope(
    *,
    status: str,
    guide: Mapping[str, Any],
    reference: Mapping[str, Any],
    max_mismatches: Any,
    include_perfect: bool,
    reason: str,
    elapsed_ms: float,
) -> dict:
    """Envelope de falha sem fingir 0 hits verificados."""
    payload = {
        "status": status,
        "reason": reason,
        "algorithm": crispr.ALGORITHM_SPCAS9_PAM_INDEX,
        "genome_wide": False,
        "scope": reference.get("scope"),
        "guide_sequence": guide.get("guide_sequence"),
        "reference_identity_hash": reference.get("identity_hash"),
        "assembly_declared": reference.get("assembly_declared") or "",
        "verified_hit_count": None,
        "hits": [],
        "truncated": status in {JOB_TIMEOUT, JOB_RESOURCE_LIMIT, JOB_INCOMPLETE},
        "max_mismatches": max_mismatches,
        "include_perfect": include_perfect,
        "elapsed_ms": elapsed_ms,
        "execution": {
            "mode": "in_process",
            "distributed_workers": False,
            "authentication": False,
            "per_user_quotas": False,
        },
        "software_version": provenance.HELIXSCOPE_VERSION,
        "timestamp": provenance.utc_now(),
        "kind": (
            "UNAVAILABLE"
            if status in {JOB_INVALID, JOB_UNAVAILABLE, "UNAVAILABLE"}
            else "ERROR"
        ),
    }
    payload["guide_specificity"] = crispr_specificity.scores_from_search(payload)
    return provenance.json_safe(payload)


def _elapsed_ms(started: float) -> float:
    """Milissegundos desde started (perf_counter)."""
    return round((time.perf_counter() - started) * 1000.0, 3)


def mismatch_summary_from_hits(
    hits: Sequence[Mapping[str, Any]], max_mismatches: int
) -> dict[int, int]:
    """Contagem real de hits por numero de mismatches, so para busca COMPLETED.

    Args:
        hits: Hits verificados.
        max_mismatches: Maior balde a incluir (0..max).

    Returns:
        Dict int -> contagem. Baldes sem hits sao 0 somente quando a busca
        correu e este resumo e pedido para hits reais.

    Raises:
        Nenhum.
    """
    summary = {key: 0 for key in range(0, max_mismatches + 1)}
    for hit in hits:
        try:
            mismatches = int(hit.get("mismatches"))
        except (TypeError, ValueError):
            continue
        if mismatches in summary:
            summary[mismatches] += 1
    return summary
