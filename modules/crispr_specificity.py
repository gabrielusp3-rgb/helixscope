"""Agregacao de especificidade CRISPR sobre buscas realmente completas.

Separa:

- Hsu 2013 single-hit (por sitio, ja calculado em crispr.hsu_single_hit_score);
- MIT Sguide (Hsu 2013 / crispr.mit.edu / CRISPOR calcMitGuideScore);
- CFD por sitio (Doench 2016);
- CFD de guia (agregacao GuideScan/CRISPOR sobre os CFD 0-1).

Nenhum agregado e emitido como genome-wide MIT Specificity ou CFD de assembly
completo a menos que a busca tenha processado um assembly declarado por inteiro.
NOT_RUN, TIMEOUT, RESOURCE_LIMIT, INCOMPLETE e ERROR nunca viram score 0.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from . import crispr_cfd, provenance

SITE_ROLE_INTENDED: str = "intended_on_target"
SITE_ROLE_ADDITIONAL_EXACT: str = "additional_exact_site"
SITE_ROLE_OFF_TARGET: str = "off_target"

MIT_SGUIDE_MODEL: str = "Hsu2013_MIT_Sguide_aggregate"
MIT_SGUIDE_VERSION: str = "CRISPOR-calcMitGuideScore"
MIT_SGUIDE_METHOD: str = (
    "MIT guide specificity Sguide = 100/(100+sum of Hsu 2013 single-hit scores) "
    "as documented at crispr.mit.edu/about and CRISPOR calcMitGuideScore "
    "(Haeussler et al. 2016). Scale 0-100. Not a per-site Hsu score. Not "
    "genome-wide unless the completed search processed a complete declared "
    "assembly."
)
MIT_SGUIDE_SCALE: str = "0-100"

CFD_SGUIDE_MODEL: str = "Doench2016_CFD_guide_specificity"
CFD_SGUIDE_VERSION: str = "GuideScan-CRISPOR-100-over-100-plus-sum"
CFD_SGUIDE_METHOD: str = (
    "Guide-level CFD specificity used by GuideScan and shown by CRISPOR since "
    "May 2019: 100/(100+sum) applied to per-hit CFD scores scaled 0-100, "
    "equivalent to 100/(1+sum of 0-1 CFD) and to FlashFry "
    "DoenchCFD_specificityscore * 100. Per-hit CFD remains 0-1. Not genome-wide "
    "unless the completed search processed a complete declared assembly. Not "
    "the Parkinson 10000/sum formula."
)
CFD_SGUIDE_SCALE: str = "0-100"

NAG_MIT_FACTOR: float = 0.2
"""CRISPOR: Hsu de PAM alternativo (ultimos 2 nt != GG) multiplicado por 0.2.

Heuristica de Haeussler/CRISPOR ('alternative PAMs ~10% of cleavage'), nao um
termo do paper Hsu 2013. CFD ja incorpora o peso de PAM AG~0.259 e nao recebe
este fator extra.
"""

COMPLETED: str = "COMPLETED"


def mit_sguide_from_hit_sum(hit_sum: float) -> dict:
    """Sguide a partir da soma dos Hsu single-hit (escala 0-100).

    Args:
        hit_sum: Soma dos scores Hsu 0-100 dos off-targets incluidos.

    Returns:
        Dict com score float 0-100, inteiro CRISPOR, formula e metadados.

    Raises:
        Nenhum.

    Nota biologica:
        CRISPOR devolve int(round(100 * 100/(100+hitSum))). HelixScope reporta
        o valor nao truncado e o inteiro CRISPOR em campos separados.
    """
    total = float(hit_sum)
    if total < 0.0:
        total = 0.0
    unit = 100.0 / (100.0 + total)
    score_0_100 = unit * 100.0
    return {
        "score": round(score_0_100, 10),
        "score_crispor_integer": int(round(score_0_100)),
        "hit_sum": round(total, 10),
        "formula": "100 * (100 / (100 + sum(Hsu_single_hit_0_100)))",
        "method": MIT_SGUIDE_METHOD,
        "model": MIT_SGUIDE_MODEL,
        "version": MIT_SGUIDE_VERSION,
        "scale": MIT_SGUIDE_SCALE,
        "status": "COMPUTED",
    }


def cfd_sguide_from_cfd_sum(cfd_sum_0_1: float) -> dict:
    """Especificidade CFD de guia a partir da soma dos CFD 0-1.

    Args:
        cfd_sum_0_1: Soma dos CFD por sitio (escala 0-1).

    Returns:
        Dict com score 0-100, soma e metadados.

    Raises:
        Nenhum.
    """
    total = float(cfd_sum_0_1)
    if total < 0.0:
        total = 0.0
    score_0_100 = 100.0 / (1.0 + total)
    return {
        "score": round(score_0_100, 10),
        "cfd_sum_0_1": round(total, 10),
        "formula": "100 / (1 + sum(CFD_0_1)) = MIT 100/(100+sum) after scaling CFD to 0-100",
        "method": CFD_SGUIDE_METHOD,
        "model": CFD_SGUIDE_MODEL,
        "version": CFD_SGUIDE_VERSION,
        "scale": CFD_SGUIDE_SCALE,
        "status": "COMPUTED",
    }


def mit_weighted_hsu(hit: Mapping[str, Any]) -> Optional[float]:
    """Hsu single-hit com o fator CRISPOR de PAM alternativo.

    Args:
        hit: Hit com score Hsu e PAM.

    Returns:
        Float 0-100 ou None se o Hsu nao foi computado.

    Raises:
        Nenhum.
    """
    raw = hit.get("score")
    if raw is None:
        raw = hit.get("hsu_hit_score")
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    pam = str(hit.get("PAM") or hit.get("pam_sequence") or "").upper().replace("U", "T")
    if len(pam) >= 2 and pam[-2:] != "GG":
        value *= NAG_MIT_FACTOR
    return value


def assign_site_roles(
    hits: Sequence[Mapping[str, Any]],
    *,
    designated: Optional[Mapping[str, Any]] = None,
) -> list[dict]:
    """Classifica on-target, copias exatas adicionais e off-targets.

    O primeiro sitio 0-mismatch, em ordem de contig/start/fita, e o on-target
    automatico (modo CRISPOR auto-ontarget) salvo se `designated` apontar
    contig+start_0based+strand. Copias 0-mismatch restantes sao
    additional_exact_site e entram no Sguide. Mismatch > 0 e off_target.

    Args:
        hits: Hits verificados.
        designated: Opcional {chromosome_or_contig, start_0based, strand}.

    Returns:
        Lista de hits com site_role. Nao esconde multiplicidade.

    Raises:
        Nenhum.
    """
    annotated = [dict(hit) for hit in hits]
    zero_keys: list[tuple[str, int, str]] = []
    for hit in annotated:
        try:
            mismatches = int(hit.get("mismatches"))
        except (TypeError, ValueError):
            hit["site_role"] = SITE_ROLE_OFF_TARGET
            continue
        if mismatches != 0:
            hit["site_role"] = SITE_ROLE_OFF_TARGET
            continue
        key = _hit_key(hit)
        zero_keys.append(key)

    intended: Optional[tuple[str, int, str]] = None
    if designated is not None:
        try:
            intended = (
                str(designated.get("chromosome_or_contig") or ""),
                int(designated.get("start_0based")),
                str(designated.get("strand") or ""),
            )
        except (TypeError, ValueError):
            intended = None
        if intended is not None and intended not in zero_keys:
            intended = None
    if intended is None and zero_keys:
        intended = sorted(zero_keys)[0]

    for hit in annotated:
        try:
            mismatches = int(hit.get("mismatches"))
        except (TypeError, ValueError):
            hit["site_role"] = SITE_ROLE_OFF_TARGET
            continue
        if mismatches != 0:
            hit["site_role"] = SITE_ROLE_OFF_TARGET
            continue
        if intended is not None and _hit_key(hit) == intended:
            hit["site_role"] = SITE_ROLE_INTENDED
        else:
            hit["site_role"] = SITE_ROLE_ADDITIONAL_EXACT
    return annotated


def scores_from_search(search: Mapping[str, Any]) -> dict:
    """MIT Sguide e CFD de guia a partir de um envelope de busca.

    Args:
        search: Envelope com status, hits, genome_wide e hashes.

    Returns:
        Dict com mit_sguide e cfd_sguide, cada um COMPUTED somente se status
        for COMPLETED. Caso contrario score e None e o status da busca e
        preservado (nunca 0).

    Raises:
        Nenhum.
    """
    status = str(search.get("status") or "NOT_RUN")
    truncated = bool(search.get("truncated"))
    genome_wide = (
        bool(search.get("genome_wide"))
        and status == "COMPLETED_FULL_REFERENCE"
        and not truncated
    )
    scope = str(search.get("scope") or "")
    base_unavailable = {
        "score": None,
        "status": status if status != COMPLETED else "UNAVAILABLE",
        "reason": (
            "Guide-level specificity is not computed unless the off-target "
            f"search completed. Current status is {status}. This is not 0."
        ),
        "genome_wide": False,
        "search_scope": scope,
        "n_sites_included": None,
        "n_intended_excluded": None,
        "n_additional_exact": None,
    }
    mit = {
        **base_unavailable,
        "method": MIT_SGUIDE_METHOD,
        "model": MIT_SGUIDE_MODEL,
        "version": MIT_SGUIDE_VERSION,
        "scale": MIT_SGUIDE_SCALE,
    }
    cfd = {
        **base_unavailable,
        "method": CFD_SGUIDE_METHOD,
        "model": CFD_SGUIDE_MODEL,
        "version": CFD_SGUIDE_VERSION,
        "scale": CFD_SGUIDE_SCALE,
    }
    if status not in {
        COMPLETED,
        "COMPLETED_FULL_REFERENCE",
        "COMPLETED_TEST_REFERENCE",
    } or truncated:
        return provenance.json_safe(
            {
                "mit_sguide": mit,
                "cfd_sguide": cfd,
                "genome_wide": False,
                "search_status": status,
            }
        )

    hits = assign_site_roles(list(search.get("hits") or []))
    included = [
        hit for hit in hits if hit.get("site_role") != SITE_ROLE_INTENDED
    ]
    intended_n = sum(1 for hit in hits if hit.get("site_role") == SITE_ROLE_INTENDED)
    extra_exact = sum(
        1 for hit in hits if hit.get("site_role") == SITE_ROLE_ADDITIONAL_EXACT
    )

    hsu_values: list[float] = []
    missing_hsu = False
    for hit in included:
        weighted = mit_weighted_hsu(hit)
        if weighted is None:
            missing_hsu = True
            break
        hsu_values.append(weighted)

    cfd_values: list[float] = []
    missing_cfd = False
    for hit in included:
        raw = hit.get("cfd_score")
        if raw is None:
            missing_cfd = True
            break
        try:
            cfd_values.append(float(raw))
        except (TypeError, ValueError):
            missing_cfd = True
            break

    scope_note = (
        "Computed over this completed search scope only. Not genome-wide MIT "
        "Specificity."
        if not genome_wide
        else (
            "CFD/MIT guide specificity over the completed READY assembly "
            f"{search.get('verified_assembly') or search.get('assembly_declared') or ''} "
            "search (COMPLETED_FULL_REFERENCE)."
        )
    )
    if missing_hsu:
        mit["status"] = "UNAVAILABLE"
        mit["reason"] = (
            "At least one included hit lacks a Hsu 2013 single-hit score. "
            "No aggregate is invented."
        )
    else:
        computed = mit_sguide_from_hit_sum(sum(hsu_values))
        mit.update(computed)
        mit["n_sites_included"] = len(included)
        mit["n_intended_excluded"] = intended_n
        mit["n_additional_exact"] = extra_exact
        mit["nag_factor"] = NAG_MIT_FACTOR
        mit["nag_factor_source"] = (
            "CRISPOR: Hsu of non-GG PAM dinucleotide multiplied by 0.2"
        )
        mit["genome_wide"] = genome_wide
        mit["search_scope"] = scope
        mit["reason"] = scope_note
        mit["reference_identity_hash"] = search.get("reference_identity_hash") or ""
        mit["assembly_declared"] = search.get("assembly_declared") or ""
        mit["algorithm"] = search.get("algorithm") or ""

    if missing_cfd:
        cfd["status"] = "UNAVAILABLE"
        cfd["reason"] = (
            "At least one included hit lacks a CFD 0-1 score. No aggregate "
            "is invented."
        )
    else:
        computed_cfd = cfd_sguide_from_cfd_sum(sum(cfd_values))
        cfd.update(computed_cfd)
        cfd["n_sites_included"] = len(included)
        cfd["n_intended_excluded"] = intended_n
        cfd["n_additional_exact"] = extra_exact
        cfd["genome_wide"] = genome_wide
        cfd["search_scope"] = scope
        cfd["reason"] = (
            "Computed over this completed search scope only. Not a genome-wide "
            "CFD specificity score."
            if not genome_wide
            else (
                "CFD guide specificity over the completed READY assembly "
                f"{search.get('verified_assembly') or search.get('assembly_declared') or ''} "
                "search (COMPLETED_FULL_REFERENCE)."
            )
        )
        cfd["per_hit_cfd_scale"] = crispr_cfd.CFD_SCALE
        cfd["per_hit_cfd_model"] = crispr_cfd.CFD_MODEL
        cfd["per_hit_cfd_version"] = crispr_cfd.CFD_VERSION
        cfd["reference_identity_hash"] = search.get("reference_identity_hash") or ""
        cfd["assembly_declared"] = search.get("assembly_declared") or ""
        cfd["algorithm"] = search.get("algorithm") or ""

    return provenance.json_safe(
        {
            "mit_sguide": mit,
            "cfd_sguide": cfd,
            "genome_wide": genome_wide,
            "search_status": status,
            "exact_site_count": intended_n + extra_exact,
            "intended_on_target_count": intended_n,
            "additional_exact_count": extra_exact,
        }
    )


def _hit_key(hit: Mapping[str, Any]) -> tuple[str, int, str]:
    try:
        start = int(hit.get("start_0based"))
    except (TypeError, ValueError):
        start = -1
    return (
        str(hit.get("chromosome_or_contig") or ""),
        start,
        str(hit.get("strand") or ""),
    )
