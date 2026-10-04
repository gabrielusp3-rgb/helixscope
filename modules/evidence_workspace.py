"""Itens de evidencia com proveniencia. Sem score de confianca inventado.

Cada afirmacao visivel no Evidence Workspace aponta para fonte, identificador,
estado (RETRIEVED/COMPUTED/...), timestamp e estado de mapeamento. Conflitos
entre fontes sao listados; tres fontes iguais nao vencem uma fonte divergente
por voto de maioria.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import json
from typing import Any, Mapping, Optional, Sequence

from . import provenance

EVIDENCE_STATUSES: tuple[str, ...] = provenance.EVIDENCE_STATUSES


class EvidenceError(RuntimeError):
    """Item de evidencia malformado.

    Attributes:
        category: INVALID_INPUT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category)


def evidence_item(
    *,
    field: str,
    value: Any,
    source: str,
    identifier: str = "",
    evidence_status: str,
    retrieved_at_utc: str = "",
    version: str = "",
    mapping_status: str = "AVAILABLE",
    note: str = "",
) -> dict:
    """Constroi um item de evidencia auditavel.

    Args:
        field: Nome do campo (ex. most_severe_consequence, rmsd_global).
        value: Valor; None permanece None (N/A), nunca vira 0 silencioso.
        source: Servico ou modulo.
        identifier: Accession, PDB ID, ticket, etc.
        evidence_status: Um de provenance.EVIDENCE_STATUSES.
        retrieved_at_utc: Timestamp ISO.
        version: Versao do software ou do servico.
        mapping_status: AVAILABLE, UNMAPPED, SKIPPED, NOT_FOUND.
        note: Texto livre curto.

    Returns:
        Dict serializavel.

    Raises:
        EvidenceError: INVALID_INPUT se o estado for desconhecido.
    """
    status = str(evidence_status or "").strip().upper()
    if status not in EVIDENCE_STATUSES:
        raise EvidenceError(
            f"Unknown evidence status {evidence_status!r}. Allowed: {EVIDENCE_STATUSES}."
        )
    mapping = str(mapping_status or "AVAILABLE").strip().upper()
    return {
        "field": str(field or "").strip(),
        "value": value,
        "source": str(source or "").strip(),
        "identifier": str(identifier or "").strip(),
        "evidence_status": status,
        "retrieved_at_utc": str(retrieved_at_utc or provenance.utc_now()),
        "version": str(version or provenance.HELIXSCOPE_VERSION),
        "mapping_status": mapping,
        "note": str(note or ""),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def collect_from_alignment(record: Mapping[str, Any]) -> list[dict]:
    """Itens a partir de um alinhamento estrutural.

    Args:
        record: parse_alignment_payload.

    Returns:
        Lista de evidence_item.

    Raises:
        Nenhum.
    """
    stamp = str(record.get("retrieved_at_utc") or provenance.utc_now())
    engine = str(record.get("engine") or "RCSB Alignment API")
    ident = str((record.get("ticket") or record.get("identity_hash") or "")[:40])
    items = [
        evidence_item(
            field="alignment_method",
            value=record.get("method_display") or record.get("method"),
            source=engine,
            identifier=ident,
            evidence_status="RETRIEVED",
            retrieved_at_utc=stamp,
        ),
        evidence_item(
            field="n_aligned_residue_pairs",
            value=record.get("n_aligned_residue_pairs"),
            source=engine,
            identifier=ident,
            evidence_status="RETRIEVED",
            retrieved_at_utc=stamp,
        ),
        evidence_item(
            field="rmsd_global_angstrom",
            value=record.get("rmsd_global_angstrom"),
            source=engine,
            identifier=ident,
            evidence_status="RETRIEVED",
            retrieved_at_utc=stamp,
            note=str(record.get("rmsd_note") or ""),
        ),
        evidence_item(
            field="rmsd_block0_angstrom",
            value=record.get("rmsd_block0_angstrom"),
            source=engine,
            identifier=ident,
            evidence_status="RETRIEVED",
            retrieved_at_utc=stamp,
        ),
        evidence_item(
            field="aln_coverage_percent",
            value=record.get("aln_coverage_percent"),
            source=engine,
            identifier=ident,
            evidence_status="RETRIEVED",
            retrieved_at_utc=stamp,
            note=str(record.get("coverage_note") or ""),
        ),
    ]
    for row in record.get("tm_scores") or []:
        items.append(
            evidence_item(
                field=str(row.get("type") or "TM-score"),
                value=row.get("value"),
                source=engine,
                identifier=ident,
                evidence_status="RETRIEVED",
                retrieved_at_utc=stamp,
                note=str(record.get("tm_score_note") or ""),
            )
        )
    return items


def collect_from_variant_result(result: Mapping[str, Any], *, label: str = "") -> list[dict]:
    """Itens a partir de explore_variant.

    Args:
        result: Dict do Variant Explorer.
        label: Prefixo opcional (Variant A / B).

    Returns:
        Lista de evidence_item.

    Raises:
        Nenhum.
    """
    prefix = f"{label} " if label else ""
    variant = result.get("variant") or {}
    layers = result.get("layers") or {}
    stamp = str(result.get("retrieved_at_utc") or provenance.utc_now())
    items: list[dict] = [
        evidence_item(
            field=f"{prefix}assembly",
            value=variant.get("assembly"),
            source="variant_core",
            identifier=str(variant.get("identity_hash") or "")[:32],
            evidence_status="COMPUTED",
            retrieved_at_utc=stamp,
        ),
        evidence_item(
            field=f"{prefix}location",
            value=(
                f"{variant.get('contig')}:{variant.get('position_1based')} "
                f"{variant.get('ref')}>{variant.get('alt')}"
            ),
            source="variant_core",
            identifier=str(variant.get("identity_hash") or "")[:32],
            evidence_status="COMPUTED",
            retrieved_at_utc=stamp,
        ),
    ]
    consequences = layers.get("consequences") or {}
    if consequences.get("status") == "AVAILABLE":
        data = consequences.get("data") or {}
        items.append(
            evidence_item(
                field=f"{prefix}most_severe_consequence",
                value=data.get("most_severe_consequence"),
                source="Ensembl VEP",
                identifier=str(data.get("allele_string") or ""),
                evidence_status="RETRIEVED",
                retrieved_at_utc=stamp,
                note="Sequence Ontology term preserved verbatim. All transcripts kept.",
            )
        )
        items.append(
            evidence_item(
                field=f"{prefix}transcript_count",
                value=data.get("transcript_count"),
                source="Ensembl VEP",
                identifier="",
                evidence_status="RETRIEVED",
                retrieved_at_utc=stamp,
            )
        )
    protein = layers.get("protein_mapping") or {}
    items.append(
        evidence_item(
            field=f"{prefix}protein_mapping",
            value=(protein.get("data") or {}).get("count") if protein.get("status") == "AVAILABLE" else None,
            source="Ensembl VEP",
            identifier="",
            evidence_status="RETRIEVED" if protein.get("status") == "AVAILABLE" else "UNMAPPED",
            mapping_status=str(protein.get("status") or "UNMAPPED"),
            retrieved_at_utc=stamp,
            note=str(protein.get("reason") or ""),
        )
    )
    clinical = layers.get("clinical_evidence") or {}
    clin_data = clinical.get("data") or {}
    items.append(
        evidence_item(
            field=f"{prefix}clinvar_status",
            value=clin_data.get("status") or clinical.get("status"),
            source="NCBI ClinVar",
            identifier="",
            evidence_status="RETRIEVED" if clinical.get("status") in {"AVAILABLE", "NOT_FOUND"} else "UNAVAILABLE",
            mapping_status=str(clinical.get("status") or ""),
            retrieved_at_utc=stamp,
            note="ClinVar is retrieved submissions, not a HelixScope diagnosis.",
        )
    )
    return items


def source_conflicts(items: Sequence[Mapping[str, Any]]) -> list[dict]:
    """Campos com o mesmo nome base e valores diferentes entre fontes.

    Args:
        items: Lista de evidence_item.

    Returns:
        Lista de conflitos. Nunca escolhe um vencedor.

    Raises:
        Nenhum.

    Nota:
        O nome do campo e comparado apos remover o prefixo 'Variant A ' / 'Variant B '.
        Assembly diferente entre variantes nao e conflito de fonte: e identidade.
    """
    grouped: dict[str, list[dict]] = {}
    for item in items:
        if not isinstance(item, Mapping):
            continue
        field = str(item.get("field") or "")
        if field.startswith("Variant A "):
            key = field[len("Variant A ") :]
        elif field.startswith("Variant B "):
            key = field[len("Variant B ") :]
        else:
            key = field
        if key in {"location", "assembly"}:
            continue
        grouped.setdefault(key, []).append(dict(item))
    conflicts: list[dict] = []
    for field, rows in grouped.items():
        values = {json.dumps(row.get("value"), default=str) for row in rows}
        sources = {str(row.get("source")) for row in rows}
        if len(values) > 1 and len(sources) > 1:
            conflicts.append(
                {
                    "field": field,
                    "items": rows,
                    "policy": (
                        "Sources disagree. HelixScope shows all values and does "
                        "not take a majority vote."
                    ),
                }
            )
    return conflicts


def export_pack(
    *,
    kind: str,
    inputs: Mapping[str, Any],
    items: Sequence[Mapping[str, Any]],
    extra: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Pacote reproduzivel da comparacao.

    Args:
        kind: structures, variants, guides, evolution, evidence.
        inputs: Identidades de entrada.
        items: Evidence items.
        extra: Metricas/hashes adicionais.

    Returns:
        Dict JSON-safe.

    Raises:
        Nenhum.
    """
    pack = {
        "kind": str(kind or ""),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "exported_at_utc": provenance.utc_now(),
        "inputs": dict(inputs),
        "evidence": [dict(item) for item in items],
        "conflicts": source_conflicts(items),
        "confidence_score": None,
        "confidence_note": (
            "HelixScope does not emit an evidence-confidence percentage. "
            "Absence of a score is N/A, not 100."
        ),
    }
    if extra:
        pack["extra"] = dict(extra)
    return pack
