"""Pontuacao on-target CRISPR publicada: Rule Set 2 / Azimuth e DeepHF.

Este modulo NAO reimplementa os modelos. Sem pesos oficiais carregaveis de
forma segura neste Python, cada pedido devolve UNAVAILABLE com score None.
Zero nunca e substituto de modelo ausente.

Rule Set 2 (Doench et al. 2016) e o modelo Azimuth 2; nao sao dois motores
distintos. DeepHF (Wang et al. 2019) e uma rede neuronal separada.

Nenhuma funcao importa Streamlit. Nenhum pickle e carregado.
"""

from __future__ import annotations

from typing import Mapping, Optional

from . import provenance

RS2_CONTEXT_NT: int = 30
"""Janela Azimuth / Rule Set 2: 4 nt 5' + 20 spacer + 3 PAM NGG + 3 nt 3'."""

RS2_SPACER_NT: int = 20
RS2_NUCLEASE: str = "SpCas9"
RS2_PAPER: str = (
    "Doench JG et al. Nature Biotechnology 34:184-191 (2016). "
    "Azimuth 2 is the Rule Set 2 on-target model (Microsoft Research, archived 2024)."
)
DEEPHF_PAPER: str = (
    "Wang D et al. Nature Communications 10:5804 (2019). "
    "Optimized CRISPR guide RNA design for two high-fidelity Cas9 variants by deep learning."
)
DEEPHF_NUCLEASES: tuple[str, ...] = (
    "SpCas9",
    "SpCas9-HF1",
    "eSpCas9",
)


def ruleset2_availability() -> dict:
    """Rule Set 2 / Azimuth nao esta embutido.

    Args:
        Nenhum.

    Returns:
        available False, reason, requires, paper.

    Raises:
        Nenhum.
    """
    return {
        "available": False,
        "model": "Doench Rule Set 2 / Azimuth 2",
        "version": "",
        "paper": RS2_PAPER,
        "nuclease": RS2_NUCLEASE,
        "context_nt": RS2_CONTEXT_NT,
        "reason": (
            "Azimuth 2 (Rule Set 2) is an archived Microsoft Research gradient-boosting "
            "model whose official pickles target Python 2 / old scikit-learn. HelixScope "
            "does not load untrusted pickle weights and does not reimplement the trees. "
            "The positional heuristic is not used as a silent substitute for Rule Set 2."
        ),
        "requires": (
            "Official Azimuth runtime and model files compatible with this Python, "
            "loaded without untrusted pickle, for SpCas9 30-nt context only."
        ),
    }


def deephf_availability() -> dict:
    """DeepHF nao esta embutido.

    Args:
        Nenhum.

    Returns:
        available False, reason, requires, paper.

    Raises:
        Nenhum.
    """
    return {
        "available": False,
        "model": "DeepHF",
        "version": "",
        "paper": DEEPHF_PAPER,
        "nucleases": list(DEEPHF_NUCLEASES),
        "reason": (
            "DeepHF is a trained recurrent network (Wang et al. 2019). Weights and a "
            "deep-learning runtime are not vendored. A house network is not DeepHF."
        ),
        "requires": "Official DeepHF weights plus a pinned inference environment.",
    }


def _na_score(
    *,
    model: str,
    paper: str,
    reason: str,
    nuclease: str,
    context: str,
    sequence_hash: str,
    extra: Optional[Mapping[str, object]] = None,
) -> dict:
    payload = {
        "score": None,
        "status": "UNAVAILABLE",
        "model": model,
        "model_version": "",
        "paper": paper,
        "nuclease": nuclease,
        "input_context": context,
        "sequence_hash": sequence_hash,
        "reason": reason,
        "not_zero": True,
    }
    if extra:
        payload.update(dict(extra))
    record = provenance.analysis_envelope(
        module="CRISPR on-target",
        payload=payload,
        status="UNAVAILABLE",
        algorithm=model,
        parameters={"nuclease": nuclease, "context_length": len(context)},
        source="HelixScope CRISPR on-target",
        input_identifier=sequence_hash,
    )
    record.update(payload)
    record["status"] = "UNAVAILABLE"
    record["score"] = None
    return record


def score_doench_ruleset2(
    context: str,
    *,
    nuclease: str = "SpCas9",
) -> dict:
    """Pede Rule Set 2. Sem modelo real o score e None, nunca 0.

    Args:
        context: Janela de 30 nt (4+20+NGG+3) ou spacer de 20 nt. Nao ha padding.
        nuclease: Apenas SpCas9 e aplicavel ao modelo original.

    Returns:
        Envelope UNAVAILABLE com score None. INVALID context tambem None.

    Raises:
        Nenhum.

    Nota biologica:
        Rule Set 2 / Azimuth 2 foi treinado para SpCas9 NGG com contexto de
        ~30 nt. Nao se aplica a Cas12/Cas13. Ausencia do modelo nao e eficiencia
        zero.
    """
    blob = "".join(ch for ch in str(context or "").upper() if not ch.isspace())
    digest = provenance.sequence_digest(blob)
    nuc = str(nuclease or "").strip()
    info = ruleset2_availability()
    extra = {
        "required_context_nt": RS2_CONTEXT_NT,
        "provided_nt": len(blob),
    }
    if nuc != RS2_NUCLEASE:
        return _na_score(
            model=info["model"],
            paper=info["paper"],
            reason=(
                "Doench Rule Set 2 / Azimuth 2 is an SpCas9 on-target model. "
                f"It is not applied to {nuc or 'an unspecified nuclease'}. Score is N/A."
            ),
            nuclease=nuc,
            context=blob,
            sequence_hash=digest,
            extra=extra,
        )
    if len(blob) not in {RS2_SPACER_NT, RS2_CONTEXT_NT}:
        return _na_score(
            model=info["model"],
            paper=info["paper"],
            reason=(
                "Rule Set 2 expects a 20-nt spacer or the published 30-nt context. "
                "Missing flanking bases are not invented. Score is N/A."
            ),
            nuclease=nuc,
            context=blob,
            sequence_hash=digest,
            extra=extra,
        )
    if any(base not in "ACGT" for base in blob):
        return _na_score(
            model=info["model"],
            paper=info["paper"],
            reason="Non-ACGT context is not scored by Rule Set 2 here. Score is N/A.",
            nuclease=nuc,
            context=blob,
            sequence_hash=digest,
            extra=extra,
        )
    return _na_score(
        model=info["model"],
        paper=info["paper"],
        reason=str(info["reason"]),
        nuclease=nuc,
        context=blob,
        sequence_hash=digest,
        extra=extra,
    )


def score_deephf(
    spacer_pam: str,
    *,
    nuclease: str = "SpCas9",
) -> dict:
    """Pede DeepHF. Sem pesos o score e None, nunca 0.

    Args:
        spacer_pam: Spacer+PAM quando o modelo existir. Sem padding.
        nuclease: Variante Cas9; Cas12 nao e DeepHF SpCas9.

    Returns:
        Envelope UNAVAILABLE com score None.

    Raises:
        Nenhum.

    Nota biologica:
        DeepHF cobre variantes Cas9 publicadas no paper de 2019, nao Cas12a.
        Rede ausente nao e probabilidade zero.
    """
    blob = "".join(ch for ch in str(spacer_pam or "").upper() if not ch.isspace())
    digest = provenance.sequence_digest(blob)
    nuc = str(nuclease or "").strip()
    info = deephf_availability()
    extra = {"provided_nt": len(blob)}
    if nuc not in DEEPHF_NUCLEASES:
        return _na_score(
            model=info["model"],
            paper=info["paper"],
            reason=(
                "DeepHF is not applied to nucleases outside the published Cas9 "
                f"variants. Requested {nuc or 'unspecified'}. Score is N/A."
            ),
            nuclease=nuc,
            context=blob,
            sequence_hash=digest,
            extra=extra,
        )
    if not blob:
        return _na_score(
            model=info["model"],
            paper=info["paper"],
            reason="Empty spacer/PAM is not scored. Score is N/A, not 0.",
            nuclease=nuc,
            context=blob,
            sequence_hash=digest,
            extra=extra,
        )
    return _na_score(
        model=info["model"],
        paper=info["paper"],
        reason=str(info["reason"]),
        nuclease=nuc,
        context=blob,
        sequence_hash=digest,
        extra=extra,
    )
