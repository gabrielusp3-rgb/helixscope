"""Orquestracao da cadeia VARIANT -> TRANSCRIPT -> PROTEIN -> STRUCTURE -> MSA.

Cada camada e independente e falha isoladamente: se o ClinVar estiver fora do ar,
o mapeamento estrutural continua; se o InterPro falhar, as consequencias do VEP
permanecem. Nenhuma falha de uma fonte apaga as outras e nenhuma ausencia de
dados e convertida em resultado negativo.

Cadeia implementada nesta fase:
    VARIANT (identidade + assembly obrigatoria)
      -> GENE / TRANSCRIPT / CDS / CODON / PROTEIN  (Ensembl VEP, retrieved)
      -> RESIDUE PROPERTIES                          (computado, sem claim de estabilidade)
      -> DOMAIN                                      (InterPro, retrieved)
      -> PROTEIN REFERENCE                           (UniProtKB, retrieved)
      -> CLINICAL EVIDENCE                           (ClinVar, retrieved, nunca inferido)
      -> MSA COLUMN / CONSERVATION                   (computado sobre MSA real)

O que este modulo NAO faz, por decisao explicita:
    - nao prediz patogenicidade nem efeito funcional;
    - nao afirma estabilizacao ou desestabilizacao de proteina;
    - nao fabrica coordenadas de proteina mutante;
    - nao converte conservacao em patogenicidade;
    - nao escolhe um transcrito canonico e esconde os restantes;
    - nao resolve divergencia entre fontes: mostra a divergencia.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import hashlib
from typing import Any, Callable, Mapping, Optional, Sequence

from . import (
    clinvar_evidence,
    ensembl_vep,
    msa as msa_module,
    protein_analysis,
    protein_domains,
    provenance,
    variant_core,
)

LAYERS: tuple[str, ...] = (
    "variant_identity",
    "consequences",
    "protein_mapping",
    "residue_properties",
    "domains",
    "protein_reference",
    "clinical_evidence",
    "evolutionary_context",
)
"""Ordem canonica das camadas do Variant Explorer."""

LAYER_STATUSES: tuple[str, ...] = (
    "AVAILABLE",
    "UNAVAILABLE",
    "UNMAPPED",
    "NOT_FOUND",
    "ERROR",
    "SKIPPED",
)
"""Estados possiveis de uma camada. UNMAPPED e diferente de NOT_FOUND e de ERROR."""

CHARGE_CLASSES: Mapping[str, str] = {
    "D": "negative",
    "E": "negative",
    "K": "positive",
    "R": "positive",
    "H": "positive_at_low_pH",
}
"""Classe de carga da cadeia lateral a pH fisiologico.

Histidina e marcada a parte porque o seu pKa (~6.0) fica proximo do pH
fisiologico: ela esta parcialmente protonada e o comportamento depende do
microambiente.
"""

RESIDUE_VOLUME_A3: Mapping[str, float] = {
    "G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1,
    "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0,
    "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7,
    "K": 168.6, "R": 173.4, "F": 189.9, "Y": 193.6, "W": 227.8,
}
"""Volume de residuo em angstrom cubico, Zamyatnin (1972), Prog Biophys Mol Biol 24:107."""


class VariantExplorerError(RuntimeError):
    """Falha na orquestracao do Variant Explorer.

    Attributes:
        category: INVALID_INPUT ou RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT").strip().upper()


def layer_result(
    name: str,
    status: str,
    *,
    data: Any = None,
    reason: str = "",
    source: str = "",
    evidence_status: str = "",
) -> dict:
    """Embrulha o resultado de uma camada com estado explicito.

    Args:
        name: Nome da camada, um valor de LAYERS.
        status: Um valor de LAYER_STATUSES.
        data: Conteudo da camada quando disponivel.
        reason: Explicacao textual, obrigatoria quando o status nao e AVAILABLE.
        source: Fonte dos dados.
        evidence_status: RETRIEVED, COMPUTED, PREDICTED, EXPERIMENTAL, MAPPED.

    Returns:
        Dict layer, status, data, reason, source, evidence_status.

    Raises:
        VariantExplorerError: INVALID_INPUT para nome de camada ou status
            desconhecido.

    Nota biologica:
        Distinguir UNMAPPED de NOT_FOUND e de ERROR e cientificamente relevante:
        "a variante nao toca proteina", "nao existe submissao clinica" e "o
        servico falhou" sao conclusoes completamente diferentes.
    """
    if name not in LAYERS:
        raise VariantExplorerError(f"Unknown Variant Explorer layer: {name}.")
    if status not in LAYER_STATUSES:
        raise VariantExplorerError(f"Unknown layer status: {status}.")
    if status != "AVAILABLE" and not str(reason or "").strip():
        raise VariantExplorerError(
            f"Layer {name} reported {status} without a reason. A non-available "
            "layer must always explain itself."
        )
    return {
        "layer": name,
        "status": status,
        "data": data,
        "reason": str(reason or ""),
        "source": str(source or ""),
        "evidence_status": str(evidence_status or ""),
    }


def _run_layer(
    name: str,
    call: Callable[[], dict],
    *,
    source: str,
    evidence_status: str,
    error_types: tuple[type[Exception], ...],
) -> dict:
    """Executa uma camada isolando a sua falha (uso interno)."""
    try:
        return call()
    except error_types as exc:
        category = str(getattr(exc, "category", "ERROR") or "ERROR")
        status = "NOT_FOUND" if category == "NOT_FOUND" else "ERROR"
        return layer_result(
            name,
            status,
            reason=(
                f"{source} layer failed with {category}: {exc}. Other layers were "
                "not affected; this is not evidence about the variant."
            ),
            source=source,
            evidence_status=evidence_status,
        )


def residue_property_change(reference_aa: str, variant_aa: str) -> dict:
    """Compara propriedades objetivas de dois aminoacidos. COMPUTED.

    Args:
        reference_aa: Aminoacido de referencia em codigo de uma letra.
        variant_aa: Aminoacido variante em codigo de uma letra.

    Returns:
        Dict com hidropatia Kyte-Doolittle de cada residuo e a diferenca, classe
        de carga, categoria de cadeia lateral, volume e mudancas booleanas.

    Raises:
        VariantExplorerError: INVALID_INPUT para simbolos fora dos 20 padrao.

    Nota biologica:
        Estes valores sao propriedades tabeladas dos aminoacidos, nao predicoes.
        Uma troca de carga ou de volume descreve a diferenca fisico-quimica; NAO
        diz se a proteina fica estabilizada ou desestabilizada. Essa pergunta
        exige um metodo validado (por exemplo FoldX, Rosetta ddG ou dados
        experimentais) que o HelixScope nao executa.
    """
    ref = str(reference_aa or "").strip().upper()
    alt = str(variant_aa or "").strip().upper()
    for symbol, label in ((ref, "reference"), (alt, "variant")):
        if symbol not in protein_analysis.KYTE_DOOLITTLE:
            raise VariantExplorerError(
                f"The {label} amino acid '{symbol}' is not one of the 20 standard "
                "residues; property comparison is refused rather than guessed."
            )
    ref_hydropathy = protein_analysis.KYTE_DOOLITTLE[ref]
    alt_hydropathy = protein_analysis.KYTE_DOOLITTLE[alt]
    ref_charge = CHARGE_CLASSES.get(ref, "neutral")
    alt_charge = CHARGE_CLASSES.get(alt, "neutral")
    ref_category = _side_chain_category(ref)
    alt_category = _side_chain_category(alt)
    ref_volume = RESIDUE_VOLUME_A3[ref]
    alt_volume = RESIDUE_VOLUME_A3[alt]
    return {
        "status": "COMPUTED",
        "evidence_status": "COMPUTED",
        "reference_amino_acid": ref,
        "variant_amino_acid": alt,
        "is_substitution": ref != alt,
        "hydropathy_reference": ref_hydropathy,
        "hydropathy_variant": alt_hydropathy,
        "hydropathy_difference": round(alt_hydropathy - ref_hydropathy, 3),
        "hydropathy_scale": "Kyte-Doolittle (1982), J Mol Biol 157:105",
        "charge_class_reference": ref_charge,
        "charge_class_variant": alt_charge,
        "charge_class_changed": ref_charge != alt_charge,
        "side_chain_category_reference": ref_category,
        "side_chain_category_variant": alt_category,
        "side_chain_category_changed": ref_category != alt_category,
        "volume_reference_a3": ref_volume,
        "volume_variant_a3": alt_volume,
        "volume_difference_a3": round(alt_volume - ref_volume, 1),
        "volume_scale": "Zamyatnin (1972), Prog Biophys Mol Biol 24:107",
        "proline_involved": "P" in {ref, alt},
        "glycine_involved": "G" in {ref, alt},
        "cysteine_involved": "C" in {ref, alt},
        "stability_claim": None,
        "method": (
            "Tabulated amino-acid properties compared directly. No model was "
            "run."
        ),
        "disclaimer": (
            "These are physico-chemical differences between two amino acids. "
            "HelixScope does not state whether the substitution is "
            "destabilizing, stabilizing, tolerated or damaging: that requires a "
            "validated stability or effect model, which is not integrated."
        ),
    }


def _side_chain_category(symbol: str) -> str:
    for category, members in protein_analysis.AMINO_ACID_CATEGORIES.items():
        if symbol in members:
            return category
    return "unclassified"


def variant_cache_key(
    *,
    variant: Mapping[str, Any],
    layer: str,
    tool: str,
    tool_version: str,
    parameters: Optional[Mapping[str, Any]] = None,
) -> str:
    """Constroi a chave de cache de uma camada do Variant Explorer.

    Args:
        variant: Objeto de `variant_core.build_variant`.
        layer: Nome da camada.
        tool: Ferramenta/serviço que produziu o resultado.
        tool_version: Versao/release da ferramenta ou fonte.
        parameters: Parametros que alteram o resultado.

    Returns:
        SHA-256 hexadecimal da chave.

    Raises:
        VariantExplorerError: INVALID_INPUT se a variante nao tiver identidade.

    Nota biologica:
        A chave inclui assembly, coordenada, alelos, camada, ferramenta e versao.
        Trocar a assembly, o alelo ou a release do VEP produz outra chave, o que
        impede servir a anotacao de uma variante diferente a partir do cache.
    """
    identity = str(variant.get("identity_hash") or "")
    if not identity:
        identifier = str(variant.get("identifier") or "")
        if not identifier:
            raise VariantExplorerError(
                "Cannot build a cache key for a variant without identity."
            )
        identity = hashlib.sha256(identifier.encode("utf-8")).hexdigest()
    parts = [
        identity,
        str(variant.get("assembly") or ""),
        str(layer or ""),
        str(tool or ""),
        str(tool_version or ""),
        provenance.HELIXSCOPE_VERSION,
    ]
    for key in sorted((parameters or {}).keys()):
        parts.append(f"{key}={(parameters or {})[key]}")
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()


def select_protein_consequences(consequences: Mapping[str, Any]) -> list[dict]:
    """Devolve os transcritos com posicao proteica resolvida. Preserva todos.

    Args:
        consequences: Resultado de `ensembl_vep.parse_vep_response`.

    Returns:
        Lista de transcritos com protein_start definido, ordenada por canonical
        primeiro e depois por transcript_id, sem descartar nenhum.

    Raises:
        Nenhum.

    Nota biologica:
        Ordenar por canonical e conveniencia de leitura. Todos os transcritos
        com proteina permanecem na lista: um mesmo nucleotido pode ser missense
        num transcrito e sinonimo noutro, e essa ambiguidade e biologica, nao
        ruido.
    """
    rows = [
        row
        for row in (consequences.get("transcript_consequences") or [])
        if isinstance(row, Mapping) and row.get("protein_start") is not None
    ]
    return sorted(
        rows,
        key=lambda row: (0 if row.get("canonical") else 1, str(row.get("transcript_id") or "")),
    )


def uniprot_accessions_from_consequence(consequence: Mapping[str, Any]) -> list[str]:
    """Extrai accessions UniProt validos de um transcrito anotado pelo VEP.

    Args:
        consequence: Uma linha de transcript_consequences.

    Returns:
        Lista de accessions BASE sem duplicados, em ordem de preferencia
        SwissProt, isoforma e por fim TrEMBL.

    Raises:
        Nenhum.

    Nota biologica:
        SwissProt (revisto manualmente) e preferido a TrEMBL (automatico). O
        sufixo de versao (`P38398.282`) identifica a versao da entrada e o
        sufixo de isoforma (`P38398-1`) identifica uma isoforma da MESMA entrada;
        ambos sao reduzidos ao accession base porque e esse que as APIs do
        UniProt e do InterPro indexam. A isoforma exata continua visivel no campo
        `uniprot_isoform` da consequencia do VEP e nao e perdida.
    """
    found: list[str] = []
    for key in ("swissprot", "uniprot_isoform", "trembl"):
        for raw in consequence.get(key) or []:
            text = str(raw or "").strip().upper()
            if not text:
                continue
            base = text.split(".", 1)[0].split("-", 1)[0]
            try:
                accession = protein_domains.validate_uniprot_accession(base)
            except protein_domains.DomainError:
                continue
            if accession not in found:
                found.append(accession)
    return found


def map_residue_to_msa_column(
    *,
    aligned_rows: Sequence[str],
    row_index: int,
    residue_1based: int,
) -> dict:
    """Traduz uma posicao de residuo numa coluna do alinhamento.

    Args:
        aligned_rows: Linhas alinhadas do MSA, todas do mesmo comprimento.
        row_index: Indice da linha que corresponde a proteina da variante.
        residue_1based: Posicao do residuo na sequencia SEM gaps.

    Returns:
        Dict column_0based, column_1based, residue_in_alignment e status.

    Raises:
        VariantExplorerError: INVALID_INPUT para indice/posicao fora do intervalo.

    Nota biologica:
        A posicao de um residuo e contada na sequencia sem gaps; a coluna e
        contada no alinhamento com gaps. Confundir as duas desloca a leitura de
        conservacao para outro residuo.
    """
    if not aligned_rows:
        raise VariantExplorerError("Alignment has no rows.")
    try:
        index = int(row_index)
        residue = int(residue_1based)
    except (TypeError, ValueError) as exc:
        raise VariantExplorerError("row_index and residue_1based must be integers.") from exc
    if index < 0 or index >= len(aligned_rows):
        raise VariantExplorerError(
            f"row_index {index} is outside the alignment ({len(aligned_rows)} rows)."
        )
    if residue < 1:
        raise VariantExplorerError("residue_1based must be >= 1.")
    row = str(aligned_rows[index])
    seen = 0
    for column, symbol in enumerate(row):
        if symbol in msa_module.GAP_CHARS:
            continue
        seen += 1
        if seen == residue:
            return {
                "status": "AVAILABLE",
                "column_0based": column,
                "column_1based": column + 1,
                "residue_in_alignment": symbol,
                "row_index": index,
                "residue_1based": residue,
                "ungapped_length": sum(
                    1 for char in row if char not in msa_module.GAP_CHARS
                ),
            }
    return {
        "status": "UNMAPPED",
        "column_0based": None,
        "column_1based": None,
        "residue_in_alignment": "",
        "row_index": index,
        "residue_1based": residue,
        "ungapped_length": seen,
        "reason": (
            f"Residue {residue} is past the end of this aligned sequence, which "
            f"has {seen} non-gap residues. The alignment row is probably a "
            "different isoform or a partial sequence."
        ),
    }


def evolutionary_context_at_residue(
    *,
    aligned_rows: Sequence[str],
    labels: Sequence[str],
    row_index: int,
    residue_1based: int,
    reference_amino_acid: str = "",
    variant_amino_acid: str = "",
    high_conservation_threshold: float = 0.9,
) -> dict:
    """Distribuicao de residuos e conservacao na coluna da variante. COMPUTED.

    Args:
        aligned_rows: Linhas alinhadas do MSA real.
        labels: Identificadores das linhas, na mesma ordem.
        row_index: Linha da proteina da variante.
        residue_1based: Posicao do residuo nessa proteina.
        reference_amino_acid: Aminoacido de referencia, para verificar coerencia.
        variant_amino_acid: Aminoacido variante, para ver se ja ocorre na coluna.
        high_conservation_threshold: Limiar declarado de alta conservacao.

    Returns:
        Dict com coluna, conservacao Shannon, fracao de gaps, contagem por
        residuo, presenca do alelo variante noutras especies e avisos.

    Raises:
        VariantExplorerError: INVALID_INPUT quando labels e linhas divergem.

    Nota biologica:
        Uma posicao altamente conservada indica pressao seletiva, o que NAO
        equivale a patogenicidade: existem posicoes conservadas toleradas e
        posicoes variaveis criticas. Encontrar o aminoacido variante noutra
        especie tambem nao prova que a troca e benigna no humano, porque o
        contexto de sequencia inteiro difere.
    """
    if len(labels) != len(aligned_rows):
        raise VariantExplorerError(
            "Alignment labels and rows must have the same length."
        )
    mapping = map_residue_to_msa_column(
        aligned_rows=aligned_rows, row_index=row_index, residue_1based=residue_1based
    )
    if mapping["status"] != "AVAILABLE":
        return {
            "status": "UNMAPPED",
            "mapping": mapping,
            "reason": str(mapping.get("reason") or "Residue could not be mapped."),
        }
    column = int(mapping["column_0based"])
    conservation = msa_module.conservation_shannon(list(aligned_rows), "protein")
    scores = conservation["scores"]
    gap_fractions = conservation["gap_fractions"]
    counts: dict[str, int] = {}
    per_label: list[dict] = []
    for index, row in enumerate(aligned_rows):
        symbol = str(row)[column] if column < len(str(row)) else "-"
        counts[symbol] = counts.get(symbol, 0) + 1
        per_label.append({"label": str(labels[index]), "symbol": symbol})
    total = sum(counts.values())
    non_gap = sum(count for symbol, count in counts.items() if symbol not in msa_module.GAP_CHARS)
    distribution = [
        {
            "symbol": symbol,
            "count": count,
            "fraction_of_all_rows": round(count / total, 6) if total else None,
            "fraction_of_non_gap": (
                round(count / non_gap, 6)
                if non_gap and symbol not in msa_module.GAP_CHARS
                else None
            ),
        }
        for symbol, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]
    score = scores[column] if column < len(scores) else float("nan")
    observed = str(mapping["residue_in_alignment"]).upper()
    reference = str(reference_amino_acid or "").strip().upper()
    variant = str(variant_amino_acid or "").strip().upper()
    warnings: list[str] = []
    if reference and observed and reference != observed:
        warnings.append(
            f"The aligned sequence carries {observed} at this column while the "
            f"variant record declares reference {reference}. The alignment row "
            "may be a different isoform, species or numbering; the conservation "
            "reading is reported but should not be trusted as this variant's "
            "column until the discrepancy is resolved."
        )
    variant_elsewhere = (
        [item["label"] for item in per_label if item["symbol"].upper() == variant]
        if variant
        else []
    )
    is_high = bool(score == score and score >= float(high_conservation_threshold))
    return {
        "status": "AVAILABLE",
        "evidence_status": "COMPUTED",
        "mapping": mapping,
        "column_1based": int(mapping["column_1based"]),
        "conservation_score": None if score != score else round(float(score), 6),
        "conservation_is_nan": score != score,
        "conservation_method": conservation["method"],
        "conservation_threshold": float(high_conservation_threshold),
        "highly_conserved_position": is_high,
        "conservation_wording": (
            "highly conserved position among the aligned sequences"
            if is_high
            else "not highly conserved among the aligned sequences"
        ),
        "gap_fraction": (
            None
            if column >= len(gap_fractions) or gap_fractions[column] != gap_fractions[column]
            else round(float(gap_fractions[column]), 6)
        ),
        "gap_treatment": conservation["gap_treatment"],
        "residue_distribution": distribution,
        "rows_at_column": per_label,
        "row_count": len(aligned_rows),
        "non_gap_count": non_gap,
        "variant_allele_present_in_labels": variant_elsewhere,
        "variant_allele_present_elsewhere": bool(variant_elsewhere),
        "warnings": warnings,
        "disclaimer": (
            "Conservation is computed over the supplied alignment only. A "
            "conserved position is not evidence of pathogenicity, and observing "
            "the variant amino acid in another species is not evidence that the "
            "substitution is tolerated in humans."
        ),
        "retrieved_at_utc": provenance.utc_now(),
    }


def cross_source_gene_agreement(
    *,
    vep: Optional[Mapping[str, Any]] = None,
    clinvar: Optional[Mapping[str, Any]] = None,
    uniprot: Optional[Mapping[str, Any]] = None,
) -> dict:
    """Compara o simbolo de gene entre fontes e expoe divergencia.

    Args:
        vep: Resultado do Ensembl VEP, se disponivel.
        clinvar: Resultado do ClinVar, se disponivel.
        uniprot: Resultado do UniProtKB, se disponivel.

    Returns:
        Dict por fonte com os simbolos observados, `agree`, `symbols` e nota.

    Raises:
        Nenhum.

    Nota biologica:
        Fontes podem discordar legitimamente: o VEP reporta todos os genes
        sobrepostos a variante, enquanto o ClinVar reporta o gene do registo
        submetido. A divergencia e mostrada, nunca resolvida em silencio por
        prioridade nao documentada.
    """
    observed: dict[str, list[str]] = {}
    if vep:
        symbols = sorted(
            {
                str(row.get("gene_symbol") or "").strip()
                for row in (vep.get("transcript_consequences") or [])
                if str(row.get("gene_symbol") or "").strip()
            }
        )
        if symbols:
            observed["Ensembl VEP"] = symbols
    if clinvar:
        symbols = sorted(
            {
                str(gene.get("symbol") or "").strip()
                for record in (clinvar.get("records") or [])
                for gene in (record.get("genes") or [])
                if str(gene.get("symbol") or "").strip()
            }
        )
        if symbols:
            observed["NCBI ClinVar"] = symbols
    if uniprot:
        symbols = sorted({str(name).strip() for name in (uniprot.get("gene_names") or []) if str(name).strip()})
        if symbols:
            observed["UniProtKB"] = symbols
    all_symbols = sorted({symbol for symbols in observed.values() for symbol in symbols})
    shared = (
        sorted(set.intersection(*[set(values) for values in observed.values()]))
        if len(observed) > 1
        else all_symbols
    )
    agree = bool(shared) and len(observed) > 1
    return {
        "sources": observed,
        "symbols": all_symbols,
        "shared_symbols": shared,
        "agree": agree if len(observed) > 1 else None,
        "comparable": len(observed) > 1,
        "note": (
            (
                f"Sources agree on {', '.join(shared)}."
                if agree
                else "Sources report different gene symbols and HelixScope does "
                "not pick a winner. Ensembl VEP lists every gene overlapping the "
                "coordinate, while ClinVar reports the gene of the submitted "
                "record; both can be correct for different reasons."
            )
            if len(observed) > 1
            else "Fewer than two sources returned a gene symbol, so there is "
            "nothing to compare."
        ),
        "policy": (
            "HelixScope has no silent source priority. When sources disagree, "
            "every value is shown with its origin."
        ),
    }


def field_provenance(
    field: str,
    *,
    value: Any,
    source: str,
    evidence_status: str,
    method: str = "",
    version: str = "",
    retrieved_at_utc: str = "",
) -> dict:
    """Anexa origem a um campo individual do Variant Explorer.

    Args:
        field: Nome do campo apresentado.
        value: Valor exibido.
        source: Base de dados ou modulo de origem.
        evidence_status: RETRIEVED, COMPUTED, MAPPED, PREDICTED ou EXPERIMENTAL.
        method: Metodo/algoritmo por tras do valor.
        version: Versao da ferramenta ou release da fonte.
        retrieved_at_utc: Momento da recuperacao.

    Returns:
        Dict com o valor e a sua linhagem.

    Raises:
        VariantExplorerError: INVALID_INPUT para evidence_status desconhecido.

    Nota biologica:
        A pergunta "de onde veio este campo" precisa de resposta ao nivel do
        campo, e nao apenas ao nivel do painel: numa mesma linha podem conviver
        um valor recuperado do VEP, um calculado localmente e um recuperado do
        ClinVar.
    """
    allowed = {"RETRIEVED", "COMPUTED", "MAPPED", "PREDICTED", "EXPERIMENTAL", "UNAVAILABLE"}
    status = str(evidence_status or "").strip().upper()
    if status not in allowed:
        raise VariantExplorerError(
            f"Unknown evidence status '{evidence_status}'. Allowed: "
            f"{', '.join(sorted(allowed))}."
        )
    return {
        "field": str(field or ""),
        "value": value,
        "source": str(source or ""),
        "evidence_status": status,
        "method": str(method or ""),
        "version": str(version or ""),
        "retrieved_at_utc": str(retrieved_at_utc or ""),
    }


def explore_variant(
    *,
    variant: Mapping[str, Any],
    email: str = "",
    enable_vep: bool = True,
    enable_clinvar: bool = True,
    enable_domains: bool = True,
    vep_urlopen_fn=None,
    clinvar_urlopen_fn=None,
    domain_urlopen_fn=None,
) -> dict:
    """Executa a cadeia de camadas para uma variante. FAZ REDE quando ativado.

    Args:
        variant: Objeto de `variant_core.build_variant`.
        email: E-mail de contacto exigido pelo NCBI para a camada ClinVar.
        enable_vep: Consulta consequencias no Ensembl VEP.
        enable_clinvar: Consulta evidencia clinica no ClinVar.
        enable_domains: Consulta dominios no InterPro e proteina no UniProtKB.
        vep_urlopen_fn: Injecao de rede para teste da camada VEP.
        clinvar_urlopen_fn: Injecao de rede para teste da camada ClinVar.
        domain_urlopen_fn: Injecao de rede para teste das camadas de proteina.

    Returns:
        Dict com `variant`, `layers` (dict por nome de camada), `cross_source`,
        `network_disclosures` e proveniencia. Cada camada traz o seu estado.

    Raises:
        VariantExplorerError: INVALID_INPUT se a variante nao tiver identidade
            normalizada nem identificador.

    Nota biologica:
        A ordem das camadas segue a cadeia biologica (genoma, transcrito, CDS,
        codao, proteina, residuo, dominio, evidencia clinica). Uma camada a
        jusante so e tentada quando a montante forneceu o que ela precisa; caso
        contrario fica UNMAPPED com a razao, em vez de produzir um valor vazio.
    """
    if not isinstance(variant, Mapping):
        raise VariantExplorerError("variant must be a mapping.")
    status = str(variant.get("status") or "")
    if status not in {"NORMALIZED", "IDENTIFIER_ONLY"}:
        raise VariantExplorerError(
            "Variant must come from variant_core.build_variant."
        )
    layers: dict[str, dict] = {}
    layers["variant_identity"] = layer_result(
        "variant_identity",
        "AVAILABLE",
        data=dict(variant),
        source="HelixScope variant_core",
        evidence_status="COMPUTED",
    )
    assembly = str(variant.get("assembly") or "")
    contig = str(variant.get("contig") or "")
    position = variant.get("position_1based")
    ref = str(variant.get("ref") or "")
    alt = str(variant.get("alt") or "")

    consequences: Optional[dict] = None
    if not enable_vep:
        layers["consequences"] = layer_result(
            "consequences",
            "SKIPPED",
            reason=(
                "Ensembl VEP was not queried. Molecular consequences are "
                "unknown, not absent."
            ),
            source="Ensembl VEP REST API",
            evidence_status="UNAVAILABLE",
        )
    elif status == "IDENTIFIER_ONLY":
        identifier = str(variant.get("identifier") or "")
        outcome = _run_layer(
            "consequences",
            lambda: layer_result(
                "consequences",
                "AVAILABLE",
                data=ensembl_vep.annotate_identifier(
                    identifier=identifier,
                    assembly=assembly,
                    urlopen_fn=vep_urlopen_fn,
                ),
                source="Ensembl VEP REST API",
                evidence_status="RETRIEVED",
            ),
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            error_types=(ensembl_vep.VepError,),
        )
        layers["consequences"] = outcome
        consequences = outcome["data"] if outcome["status"] == "AVAILABLE" else None
    else:
        ensembl_contig = variant_core.ensembl_contig_name(contig)
        outcome = _run_layer(
            "consequences",
            lambda: layer_result(
                "consequences",
                "AVAILABLE",
                data=ensembl_vep.annotate_region(
                    contig=ensembl_contig,
                    position_1based=int(position),
                    ref=ref,
                    alt=alt,
                    assembly=assembly,
                    urlopen_fn=vep_urlopen_fn,
                ),
                source="Ensembl VEP REST API",
                evidence_status="RETRIEVED",
            ),
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            error_types=(ensembl_vep.VepError,),
        )
        layers["consequences"] = outcome
        consequences = outcome["data"] if outcome["status"] == "AVAILABLE" else None

    protein_rows = select_protein_consequences(consequences) if consequences else []
    if consequences is None:
        layers["protein_mapping"] = layer_result(
            "protein_mapping",
            "UNMAPPED",
            reason=(
                "No consequence layer is available, so no transcript could be "
                "mapped to a protein position."
            ),
            source="Ensembl VEP REST API",
            evidence_status="UNAVAILABLE",
        )
    elif not protein_rows:
        layers["protein_mapping"] = layer_result(
            "protein_mapping",
            "UNMAPPED",
            reason=(
                "Ensembl VEP returned consequences but none of them places this "
                "variant inside a coding sequence, so there is no protein "
                "position. This is a genuine biological outcome, not a failure."
            ),
            source="Ensembl VEP REST API",
            evidence_status="RETRIEVED",
        )
    else:
        layers["protein_mapping"] = layer_result(
            "protein_mapping",
            "AVAILABLE",
            data={
                "transcripts": protein_rows,
                "count": len(protein_rows),
                "canonical_count": sum(1 for row in protein_rows if row.get("canonical")),
                "note": (
                    f"{len(protein_rows)} transcript(s) place this variant in a "
                    "coding sequence. All are kept; none is hidden."
                ),
            },
            source="Ensembl VEP REST API",
            evidence_status="RETRIEVED",
        )

    primary = protein_rows[0] if protein_rows else None
    if not primary:
        layers["residue_properties"] = layer_result(
            "residue_properties",
            "UNMAPPED",
            reason="No protein position is available, so no residue comparison was made.",
            source="HelixScope variant_explorer",
            evidence_status="UNAVAILABLE",
        )
    else:
        reference_aa = str(primary.get("reference_amino_acid") or "")
        variant_aa = str(primary.get("variant_amino_acid") or "")
        if not reference_aa or not variant_aa:
            layers["residue_properties"] = layer_result(
                "residue_properties",
                "UNMAPPED",
                reason=(
                    "The transcript has a protein position but Ensembl VEP did "
                    "not report both amino acids, which happens for synonymous "
                    "and frameshift consequences."
                ),
                source="Ensembl VEP REST API",
                evidence_status="RETRIEVED",
            )
        else:
            try:
                comparison = residue_property_change(reference_aa, variant_aa)
                comparison["transcript_id"] = str(primary.get("transcript_id") or "")
                comparison["protein_id"] = str(primary.get("protein_id") or "")
                comparison["protein_position"] = primary.get("protein_start")
                layers["residue_properties"] = layer_result(
                    "residue_properties",
                    "AVAILABLE",
                    data=comparison,
                    source="HelixScope variant_explorer (tabulated scales)",
                    evidence_status="COMPUTED",
                )
            except VariantExplorerError as exc:
                layers["residue_properties"] = layer_result(
                    "residue_properties",
                    "UNMAPPED",
                    reason=str(exc),
                    source="HelixScope variant_explorer",
                    evidence_status="UNAVAILABLE",
                )

    accessions = uniprot_accessions_from_consequence(primary) if primary else []
    uniprot_data: Optional[dict] = None
    domain_data: Optional[dict] = None
    if not enable_domains:
        reason = (
            "InterPro and UniProt were not queried. Domain context is unknown, "
            "not absent."
        )
        layers["domains"] = layer_result(
            "domains", "SKIPPED", reason=reason, source="InterPro (EMBL-EBI)",
            evidence_status="UNAVAILABLE",
        )
        layers["protein_reference"] = layer_result(
            "protein_reference", "SKIPPED", reason=reason, source="UniProtKB",
            evidence_status="UNAVAILABLE",
        )
    elif not accessions:
        reason = (
            "No UniProt accession was reported for this transcript, so protein "
            "annotation could not be looked up."
        )
        layers["domains"] = layer_result(
            "domains", "UNMAPPED", reason=reason, source="InterPro (EMBL-EBI)",
            evidence_status="UNAVAILABLE",
        )
        layers["protein_reference"] = layer_result(
            "protein_reference", "UNMAPPED", reason=reason, source="UniProtKB",
            evidence_status="UNAVAILABLE",
        )
    else:
        accession = accessions[0]
        protein_outcome = _run_layer(
            "protein_reference",
            lambda: layer_result(
                "protein_reference",
                "AVAILABLE",
                data=protein_domains.fetch_uniprot_protein(
                    accession=accession, urlopen_fn=domain_urlopen_fn
                ),
                source="UniProtKB",
                evidence_status="RETRIEVED",
            ),
            source="UniProtKB",
            evidence_status="RETRIEVED",
            error_types=(protein_domains.DomainError,),
        )
        layers["protein_reference"] = protein_outcome
        uniprot_data = protein_outcome["data"] if protein_outcome["status"] == "AVAILABLE" else None
        residue = primary.get("protein_start") if primary else None
        domain_outcome = _run_layer(
            "domains",
            lambda: layer_result(
                "domains",
                "AVAILABLE",
                data=_domain_payload(
                    accession=accession,
                    residue=residue,
                    urlopen_fn=domain_urlopen_fn,
                ),
                source="InterPro (EMBL-EBI)",
                evidence_status="RETRIEVED",
            ),
            source="InterPro",
            evidence_status="RETRIEVED",
            error_types=(protein_domains.DomainError,),
        )
        layers["domains"] = domain_outcome
        domain_data = domain_outcome["data"] if domain_outcome["status"] == "AVAILABLE" else None

    clinvar_data: Optional[dict] = None
    if not enable_clinvar:
        layers["clinical_evidence"] = layer_result(
            "clinical_evidence",
            "SKIPPED",
            reason=(
                "ClinVar was not queried. No clinical submission was retrieved; "
                "this is not evidence that the variant is benign."
            ),
            source="NCBI ClinVar",
            evidence_status="UNAVAILABLE",
        )
    elif status != "NORMALIZED":
        layers["clinical_evidence"] = layer_result(
            "clinical_evidence",
            "UNMAPPED",
            reason=(
                "ClinVar lookup needs a resolved coordinate to verify that a "
                "retrieved record really belongs to this variant."
            ),
            source="NCBI ClinVar",
            evidence_status="UNAVAILABLE",
        )
    else:
        chrom = variant_core.ensembl_contig_name(contig)
        outcome = _run_layer(
            "clinical_evidence",
            lambda: _clinvar_layer(
                assembly=assembly,
                contig=chrom,
                position_1based=int(position),
                ref=ref,
                alt=alt,
                email=email,
                urlopen_fn=clinvar_urlopen_fn,
            ),
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
            error_types=(clinvar_evidence.ClinVarError,),
        )
        layers["clinical_evidence"] = outcome
        clinvar_data = outcome["data"] if outcome["status"] == "AVAILABLE" else None

    layers["evolutionary_context"] = layer_result(
        "evolutionary_context",
        "SKIPPED",
        reason=(
            "Evolutionary context needs an alignment. Load an MSA whose rows "
            "include this protein and call evolutionary_context_at_residue with "
            "the mapped residue position."
        ),
        source="HelixScope msa",
        evidence_status="UNAVAILABLE",
    )

    return {
        "variant": dict(variant),
        "layers": layers,
        "layer_order": list(LAYERS),
        "cross_source": cross_source_gene_agreement(
            vep=consequences, clinvar=clinvar_data, uniprot=uniprot_data
        ),
        "network_disclosures": {
            "ensembl_vep": ensembl_vep.network_disclosure(assembly),
            "ncbi_clinvar": clinvar_evidence.network_disclosure(),
            "interpro_uniprot": protein_domains.network_disclosure(),
        },
        "status_summary": {
            name: layers.get(name, {}).get("status", "SKIPPED") for name in LAYERS
        },
        "no_cascade_failure": (
            "Every layer is executed independently. A failure in one source is "
            "recorded on that layer only and never removes results from another."
        ),
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def _domain_payload(*, accession: str, residue: Any, urlopen_fn) -> dict:
    """Junta entradas InterPro e cobertura do residuo (uso interno)."""
    entries = protein_domains.fetch_protein_domains(
        accession=accession, urlopen_fn=urlopen_fn
    )
    payload = dict(entries)
    if residue is None:
        payload["coverage"] = None
        payload["coverage_note"] = (
            "No protein residue position was available, so domain coverage was "
            "not evaluated."
        )
        return payload
    payload["coverage"] = protein_domains.domains_covering_residue(
        entries, residue_1based=int(residue)
    )
    payload["coverage_note"] = (
        "Coverage states whether the residue falls inside an annotated range. "
        "It does not demonstrate a functional effect."
    )
    return payload


def _clinvar_layer(
    *,
    assembly: str,
    contig: str,
    position_1based: int,
    ref: str,
    alt: str,
    email: str,
    urlopen_fn,
) -> dict:
    """Consulta ClinVar por coordenada e rotula o resultado (uso interno)."""
    if not str(email or "").strip():
        return layer_result(
            "clinical_evidence",
            "SKIPPED",
            reason=(
                "NCBI requires a contact e-mail for E-utilities. ClinVar was not "
                "queried; no clinical submission was retrieved and that is not "
                "evidence that the variant is benign."
            ),
            source="NCBI ClinVar",
            evidence_status="UNAVAILABLE",
        )
    term = clinvar_evidence.build_search_term(
        assembly=assembly, contig=contig, position_1based=position_1based
    )
    payload = clinvar_evidence.fetch_clinvar_records(
        term=term,
        email=email,
        retmax=10,
        urlopen_fn=urlopen_fn,
        assembly=assembly,
        contig=contig,
        position_1based=position_1based,
        ref=ref,
        alt=alt,
    )
    if payload.get("status") == "NOT_FOUND":
        return layer_result(
            "clinical_evidence",
            "NOT_FOUND",
            data=payload,
            reason=str(payload.get("message") or "No ClinVar record retrieved."),
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
        )
    return layer_result(
        "clinical_evidence",
        "AVAILABLE",
        data=payload,
        source="NCBI ClinVar",
        evidence_status="RETRIEVED",
    )
