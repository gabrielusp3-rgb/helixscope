"""HelixScope: aplicacao Streamlit de analise de bioinformatica.

Orquestrador da interface. Esta primeira metade implementa a barra lateral de
navegacao e as abas de analise de DNA e de RNA; as demais abas (proteina,
alinhamento, NCBI e CRISPR) sao adicionadas na etapa seguinte.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import time
from typing import Any, Mapping, Optional

import pandas as pd
import streamlit as st
from Bio.Seq import Seq

from helixscope_core.compat import (
    alignment,
    blast_search,
    clinvar_evidence,
    comparative,
    complex_structure,
    crispr,
    crispr_assemblies,
    crispr_casoffinder,
    crispr_cfd,
    crispr_offtarget,
    crispr_reference,
    crispr_specificity,
    crispr_structure_catalog,
    dna_analysis,
    dna_3d,
    engine_validation,
    ensembl_vep,
    evidence_workspace,
    genome_annotation,
    genome_coordinates,
    genome_download,
    genome_fasta,
    genome_jobs,
    genome_region,
    genome_store,
    opencl_runtime,
    molecule_selection,
    motif_search,
    msa,
    ncbi_fetch,
    nucleic_geometry,
    phylogeny,
    provenance,
    protein_analysis,
    protein_domains,
    protein_structure,
    rcsb_alignment,
    region_nav,
    resource_admission,
    rna_analysis,
    rna_folding,
    rna_3d,
    scale_profile,
    scientific_checks,
    structure_alignment,
    taxonomy,
    tool_detection,
    thermodynamics,
    tool_registry,
    usalign,
    variant_core,
    variant_explorer,
)
from modules import runtime
from ui import charts
from ui import components as ui_components
from ui import structure_viewer
from ui import workspace as helix_workspace
from ui.styles import inject_custom_css
from ui.navigation import persist_orphaned_widget_keys, render_sidebar_navigation
from ui import shell as helix_shell
from ui.hand_control import render_hand_control

_RUNTIME_STATUS = runtime.ensure_runtime(expected_version=provenance.HELIXSCOPE_VERSION)

# Streamlit on Windows re-runs app.py in the same interpreter. Mass
# importlib.reload of every scientific module on each rerun is not the
# architecture: modules.runtime reloads a module only if a required
# attribute is missing, then asks for a process restart if it is still stale.
# Do not reload protein_structure / BLAST / MSA caches here.

badge = ui_components.badge
badge_row = ui_components.badge_row
glass_card = ui_components.glass_card
html_escape = ui_components.html_escape
inspector_panel = ui_components.inspector_panel
meta_grid = ui_components.meta_grid
render_metric_grid = helix_workspace.render_metric_grid
safe_download_filename = ui_components.safe_download_filename
sequence_display = ui_components.sequence_display
status_badge = ui_components.status_badge

MAX_UPLOAD_BYTES: int = 1_048_576  # FASTA upload ceiling (1 MiB), aligned with Streamlit server.maxUploadSize.

NCBI_MIN_INTERVAL_S: float = 1.0

st.set_page_config(
    page_title="HelixScope",
    layout="wide",
    initial_sidebar_state="expanded",
)
inject_custom_css()
helix_shell.inject_shell_css()
persist_orphaned_widget_keys()
if _RUNTIME_STATUS.get("needs_restart"):
    st.error(str(_RUNTIME_STATUS.get("restart_reason") or "Restart Streamlit to load current HelixScope APIs."))
    st.caption(
        "Stale module APIs: "
        + ", ".join(str(item) for item in (_RUNTIME_STATUS.get("missing") or [])[:12])
    )


def _section_toggle(label: str, *, key: str, expanded: bool = False) -> bool:
    """Toggle de secao recolhivel sem icones Material (uso interno)."""
    return st.toggle(label, value=expanded, key=key)


def _sequence_hashes_match(stored: object, current: object) -> bool:
    """Compara hashes; fallback se o modulo provenance em memoria estiver obsoleto."""
    matcher = getattr(provenance, "hashes_match", None)
    if callable(matcher):
        return bool(matcher(stored, current))
    left = str(stored or "")
    right = str(current or "")
    return bool(left) and left == right


def _warn_if_analyzed_input_changed(widget_text: str, analyzed: str, *, action: str) -> None:
    """Avisa quando o texto visivel ja nao e a sequencia analisada."""
    if not str(widget_text or "").strip() or not str(analyzed or "").strip():
        return
    try:
        live = _resolve_sequence(widget_text, None)
    except ValueError:
        return
    live_n = _normalized_seq(live)
    analyzed_n = _normalized_seq(analyzed)
    if live_n and analyzed_n and not _sequence_hashes_match(
        provenance.sequence_digest(live_n),
        provenance.sequence_digest(analyzed_n),
    ):
        st.warning(
            f"The text area changed after {action}. Results below belong to the "
            f"last {action}, not the current text."
        )


def _hide_stale_result(message: str) -> None:
    st.info(message)


def _validate(seq: str) -> dict:
    """Valida e detecta o tipo de uma sequencia (uso interno)."""
    return dna_analysis.validate_sequence(seq)


def _validate_for(seq: str, molecule: str) -> dict:
    """Valida a sequencia contra um unico tipo molecular (uso interno).

    Recarrega dna_analysis se o processo Streamlit ainda tiver um modulo antigo
    em memoria (comum no Windows quando so app.py e recarregado).
    """
    fn = runtime.require_attr(dna_analysis, "validate_for_molecule")
    if fn is None:
        raise RuntimeError(
            "Sequence validation is unavailable. Stop the app and run: "
            "python -m streamlit run app.py --server.port 8501"
        )
    return fn(seq, molecule)


@st.cache_data(show_spinner=False)
def _composition(seq: str) -> dict:
    """Composicao de nucleotideos, com cache (uso interno)."""
    return dna_analysis.nucleotide_composition(seq)


@st.cache_data(show_spinner=False)
def _gc_content(seq: str) -> float:
    """Conteudo GC, com cache (uso interno)."""
    return dna_analysis.gc_content(seq)


@st.cache_data(show_spinner=False)
def _melting_temperature(seq: str) -> float:
    """Temperatura de melting, com cache (uso interno)."""
    return dna_analysis.melting_temperature(seq)


@st.cache_data(show_spinner=False)
def _melting_temperature_report(seq: str) -> dict:
    """Relatorio de Tm com rotulo de metodo, com cache (uso interno)."""
    return dna_analysis.melting_temperature_report(seq)


@st.cache_data(show_spinner=False)
def _santalucia_tm(seq: str, na_m: float, oligo_nm: float, mg_m: float) -> dict:
    """Tm nearest-neighbor SantaLucia, com cache (uso interno)."""
    fn = getattr(thermodynamics, "santalucia_nearest_neighbor_tm", None)
    if fn is None:
        return {
            "status": "UNAVAILABLE",
            "value_c": None,
            "reason": "SantaLucia engine is not loaded in this process. Restart the app.",
        }
    return fn(seq, na_m=na_m, oligo_nm=oligo_nm, mg_m=mg_m)


@st.cache_data(show_spinner=False)
def _kmer_summary_dna(seq: str, k: int) -> dict:
    """Resumo de k-mers de DNA, com cache (uso interno)."""
    return dna_analysis.kmer_summary(seq, k)


@st.cache_data(show_spinner=False)
def _windowed_profiles(seq: str, window: int, step: int) -> list:
    """Perfis GC/AT/skew/entropia na mesma grade, com cache (uso interno)."""
    return dna_analysis.windowed_profiles(seq, window=window, step=step)


@st.cache_data(show_spinner=False)
def _molecular_weight(seq: str, mol_type: str) -> float:
    """Massa molecular do acido nucleico, com cache (uso interno)."""
    return dna_analysis.molecular_weight(seq, mol_type)


@st.cache_data(show_spinner=False)
def _gc_skew(seq: str, window: int) -> list:
    """GC skew por janela, com cache (uso interno)."""
    return dna_analysis.gc_skew(seq, window)


@st.cache_data(show_spinner=False)
def _restriction_sites(seq: str) -> dict:
    """Sitios de restricao, com cache (uso interno)."""
    return dna_analysis.find_restriction_sites(seq)


@st.cache_data(show_spinner=False)
def _orfs(seq: str, min_length: int) -> list:
    """ORFs nos seis frames, com cache (uso interno)."""
    return dna_analysis.find_orfs(seq, min_length)


@st.cache_data(show_spinner=False)
def _orfs_non_overlapping(orfs: list, min_length: int) -> list:
    """ORFs sem sobreposicao redundante, com cache (uso interno)."""
    return dna_analysis.select_non_overlapping_orfs(orfs, min_length=min_length)


@st.cache_data(show_spinner=False)
def _cpg_islands(seq: str) -> list:
    """Ilhas CpG, com cache (uso interno)."""
    return dna_analysis.cpg_islands(seq)


@st.cache_data(show_spinner=False)
def _at_content(seq: str) -> float:
    """Conteudo AT, com cache (uso interno)."""
    return dna_analysis.at_content(seq)


@st.cache_data(show_spinner=False)
def _dinucleotides(seq: str) -> dict:
    """Frequencias e razao observado/esperado de dinucleotideos (uso interno)."""
    return dna_analysis.dinucleotide_frequencies(seq)


@st.cache_data(show_spinner=False)
def _gc_sliding_window(seq: str, window: int, step: int) -> list:
    """Perfil de GC em janela deslizante, com cache (uso interno)."""
    return dna_analysis.gc_sliding_window(seq, window=window, step=step)


@st.cache_data(show_spinner=False)
def _cumulative_skew_summary(seq: str) -> dict:
    """Skew cumulativo de Lobry, com cache. O grafico nao substitui a soma."""
    return dna_analysis.cumulative_skew_summary(seq)


@st.cache_data(show_spinner=False)
def _shannon_entropy(seq: str) -> float:
    """Entropia de Shannon de DNA, com cache (uso interno)."""
    return dna_analysis.shannon_entropy(seq)


@st.cache_data(show_spinner=False)
def _at_skew(seq: str) -> float:
    """AT skew global de DNA, com cache (uso interno)."""
    return dna_analysis.at_skew(seq)


@st.cache_data(show_spinner=False)
def _at_skew_windows(seq: str, window: int) -> list:
    """AT skew por janela, com cache (uso interno)."""
    return dna_analysis.at_skew_windows(seq, window)


@st.cache_data(show_spinner=False)
def _kmer_counts_dna(seq: str, k: int) -> dict:
    """Contagem de k-mers de DNA, com cache (uso interno)."""
    return dna_analysis.kmer_counts(seq, k)


@st.cache_data(show_spinner=False)
def _entropy_profile(seq: str, window: int, step: int) -> list:
    """Perfil de entropia de DNA, com cache (uso interno)."""
    return dna_analysis.entropy_sliding_window(seq, window=window, step=step)


@st.cache_data(show_spinner=False)
def _rna_shannon_entropy(seq: str) -> float:
    """Entropia de Shannon de RNA, com cache (uso interno)."""
    return rna_analysis.shannon_entropy(seq)


@st.cache_data(show_spinner=False)
def _kmer_counts_rna(seq: str, k: int) -> dict:
    """Contagem de k-mers de RNA, com cache (uso interno)."""
    return rna_analysis.kmer_counts(seq, k)


@st.cache_data(show_spinner=False, ttl=1800)
def _annotation_bundle(accession: str, db: str, _email: str) -> dict:
    """Recupera anotacao oficial e extrai as regioes CDS, com cache (uso interno).

    O e-mail de contato nao entra na chave de cache.
    """
    record = ncbi_fetch.fetch_by_accession(
        accession, _email, db, api_key=ncbi_fetch.api_key_from_environment()
    )
    return {
        "accession": str(record["accession"]),
        "organism": str(record["organism"]),
        "description": str(record["description"]),
        "regions": ncbi_fetch.extract_cds_regions(record),
    }


@st.cache_data(show_spinner=False)
def _resolve_coding_codons(
    sequence: str,
    cds_regions: list,
    source: str,
    min_cds_length: int,
    drop_overlapping: bool,
    allow_whole_sequence: bool,
) -> dict:
    """Resolve regioes codificadoras e extrai codons, com cache (uso interno)."""
    return rna_analysis.resolve_coding_codons(
        sequence,
        cds_regions=cds_regions or None,
        source=source,
        min_cds_length=min_cds_length,
        drop_overlapping=drop_overlapping,
        allow_whole_sequence=allow_whole_sequence,
    )


@st.cache_data(show_spinner=False)
def _codon_usage_table(seq: str) -> pd.DataFrame:
    """Tabela de uso de codons, com cache (uso interno)."""
    return rna_analysis.codon_usage_table(seq)


@st.cache_data(show_spinner=False)
def _codon_adaptation_index(seq: str, organism: str) -> float:
    """Codon Adaptation Index, com cache (uso interno)."""
    return rna_analysis.codon_adaptation_index(seq, organism)


@st.cache_data(show_spinner=False)
def _codon_adaptation_index_report(seq: str, organism: str) -> dict:
    """CAI com tabela de referencia e metodo, com cache (uso interno)."""
    return rna_analysis.codon_adaptation_index_report(seq, organism)


@st.cache_data(show_spinner=False)
def _effective_number_of_codons_report(seq: str) -> dict:
    """ENC com status de aplicabilidade, com cache (uso interno)."""
    return rna_analysis.effective_number_of_codons_report(seq)


@st.cache_data(show_spinner=False)
def _rscu_table(seq: str) -> pd.DataFrame:
    """Tabela de RSCU por codon, com cache (uso interno)."""
    return rna_analysis.relative_synonymous_codon_usage(seq)


@st.cache_data(show_spinner=False)
def _effective_number_of_codons(seq: str) -> float:
    """Numero Efetivo de Codons (ENC), com cache (uso interno)."""
    return rna_analysis.effective_number_of_codons(seq)


@st.cache_data(show_spinner=False)
def _translate(seq: str) -> str:
    """Traducao de mRNA em proteina, com cache (uso interno)."""
    return protein_analysis.translate(seq)


@st.cache_data(show_spinner=False)
def _amino_acid_composition(seq: str) -> dict:
    """Composicao de aminoacidos, com cache (uso interno)."""
    return protein_analysis.amino_acid_composition(seq)


@st.cache_data(show_spinner=False)
def _gravy_index(seq: str) -> float:
    """GRAVY (hidropatia media global), com cache (uso interno)."""
    return protein_analysis.gravy_index(seq)


@st.cache_data(show_spinner=False)
def _net_charge(seq: str, ph: float) -> float:
    """Carga liquida em um pH definido, com cache (uso interno)."""
    return protein_analysis.net_charge(seq, ph)


@st.cache_data(show_spinner=False)
def _aliphatic_index(seq: str) -> float:
    """Indice alifatico de Ikai, com cache (uso interno)."""
    return protein_analysis.aliphatic_index(seq)


@st.cache_data(show_spinner=False)
def _amino_acid_categories(seq: str) -> dict:
    """Categorias de cadeia lateral dos residuos, com cache (uso interno)."""
    return protein_analysis.amino_acid_categories(seq)


@st.cache_data(show_spinner=False)
def _protein_mw(seq: str) -> float:
    """Massa molecular da proteina em kDa, com cache (uso interno)."""
    return protein_analysis.molecular_weight_protein(seq)


@st.cache_data(show_spinner=False)
def _physicochemical_report(seq: str, ph: float = 7.0) -> dict:
    """Pacote de propriedades teoricas da proteina, com cache (uso interno)."""
    return protein_analysis.physicochemical_report(seq, ph=ph)


@st.cache_data(show_spinner=False)
def _isoelectric_point(seq: str) -> float:
    """Ponto isoeletrico, com cache (uso interno)."""
    return protein_analysis.isoelectric_point(seq)


@st.cache_data(show_spinner=False)
def _instability_index(seq: str) -> float:
    """Indice de instabilidade, com cache (uso interno)."""
    return protein_analysis.instability_index(seq)


@st.cache_data(show_spinner=False)
def _hydrophobicity_profile(seq: str, window: int) -> list:
    """Perfil de hidrofobicidade Kyte-Doolittle, com cache (uso interno)."""
    return protein_analysis.hydrophobicity_profile(seq, window)


@st.cache_data(show_spinner=False)
def _aromaticity(seq: str) -> float:
    """Aromaticidade ProtParam, com cache (uso interno)."""
    return protein_analysis.aromaticity(seq)


@st.cache_data(show_spinner=False)
def _charge_profile(seq: str, window: int) -> list:
    """Perfil de carga formal, com cache (uso interno)."""
    return protein_analysis.charge_count_profile(seq, window)


@st.cache_data(show_spinner=False)
def _pairwise(seq1: str, seq2: str, mode: str) -> dict:
    """Alinhamento global ou local, com cache (uso interno)."""
    if mode == "local":
        return alignment.pairwise_local(seq1, seq2)
    return alignment.pairwise_global(seq1, seq2)


@st.cache_data(show_spinner=False)
def _dotplot_matrix(seq1: str, seq2: str):
    """Matriz de dotplot binaria, com cache (uso interno)."""
    return alignment.dotplot_matrix(seq1, seq2)


@st.cache_data(show_spinner=False, ttl=1800)
def _fetch_by_accession(accession: str, db: str, _email: str) -> dict:
    """Recupera um registro do NCBI por acesso, com cache (uso interno).

    O e-mail de contato nao entra na chave de cache.
    """
    return ncbi_fetch.fetch_by_accession(
        accession, _email, db, api_key=ncbi_fetch.api_key_from_environment()
    )


@st.cache_data(show_spinner=False, ttl=1800)
def _search_ncbi(term: str, db: str, _email: str, retstart: int = 0) -> dict:
    """Busca Entrez esearch+esummary, com cache (uso interno)."""
    return ncbi_fetch.search_records(
        term,
        _email,
        db,
        api_key=ncbi_fetch.api_key_from_environment(),
        retstart=retstart,
    )


@st.cache_data(show_spinner=False)
def _design_guides(sequence: str, cas_system: str) -> list:
    """Encontra e classifica guias CRISPR, com cache (uso interno)."""
    return crispr.rank_guides(crispr.find_guides(sequence, cas_system, guide_length=20))


@st.cache_data(show_spinner=False)
def _evaluate_guides(
    sequence: str, cas_system: str, max_mismatches: int, run_off_target: bool
) -> list:
    """Encontra e avalia guias CRISPR com metricas por guia, com cache (interno)."""
    found = crispr.find_guides(sequence, cas_system, guide_length=20)
    evaluated = crispr.evaluate_guides(
        found,
        sequence,
        cas_system,
        max_mismatches=max_mismatches,
        run_off_target=run_off_target,
    )
    for guide in evaluated:
        if not scientific_checks.crispr_guide_in_sequence(sequence, guide):
            raise ValueError(
                "Analysis failed: a CRISPR candidate does not match the input "
                "sequence. No fabricated guide was kept."
            )
    return evaluated


@st.cache_data(show_spinner=False)
def _off_targets(
    guide_sequence: str, sequence: str, max_mismatches: int, pam_sequence: str = ""
) -> list:
    """Varre off-targets de um guia na sequencia colada, com cache (uso interno)."""
    return crispr.find_off_targets(
        guide_sequence,
        sequence,
        max_mismatches,
        cas_system="SpCas9",
        require_pam=True,
        on_target_pam=pam_sequence,
    )


@st.cache_data(show_spinner=False)
def _reference_offtarget_search(
    guide_sequence: str,
    pam_sequence: str,
    strand: str,
    position: int,
    fasta_text: str,
    source: str,
    organism: str,
    assembly: str,
    version: str,
    accession: str,
    upload_filename: str,
    max_mismatches: int,
    include_perfect: bool,
    target_hash: str,
    engine: str = "pam_index",
) -> dict:
    """Busca off-targets na referencia fornecida. A chave inclui assembly, hashes e motor."""
    reference = crispr_reference.load_reference_from_text(
        fasta_text,
        source=source,
        organism=organism,
        assembly=assembly,
        version=version,
        accession=accession,
        upload_filename=upload_filename,
    )
    guide = {
        "guide_sequence": guide_sequence,
        "pam_sequence": pam_sequence,
        "strand": strand,
        "position": position,
    }
    if engine == "cas_offinder":
        return crispr_casoffinder.search_with_cas_offinder(
            guide,
            reference,
            max_mismatches=max_mismatches,
            include_perfect=include_perfect,
            target_hash=target_hash,
        )
    return crispr_offtarget.search_off_targets_in_reference(
        guide,
        reference,
        max_mismatches=max_mismatches,
        include_perfect=include_perfect,
        target_hash=target_hash,
    )


@st.cache_data(show_spinner=False)
def _primer_pair(guide_sequence: str, sequence: str) -> dict:
    """Desenha o par de primers de validacao, com cache (uso interno)."""
    return crispr.design_primer_pair(guide_sequence, sequence)


def _reject_nonportable_upload_name(name: str) -> None:
    """Refuse an upload name that is not a single portable file name.

    Args:
        name: The name reported by the browser. It is never used as a path.

    Raises:
        ValueError: The name contains a separator, a parent segment, or is empty.
    """
    text = str(name or "")
    if (
        not text
        or text in {".", ".."}
        or ".." in text
        or "/" in text
        or "\\" in text
        or ":" in text
        or "\x00" in text
    ):
        raise ValueError(
            "Upload rejected. Use a simple file name without slashes or '..'. "
            "The name is not used as a path, and the file was not analyzed."
        )


def _read_uploaded_text(uploaded_file) -> str:
    """Le o FASTA enviado, recusando nome inseguro, vazio ou acima de MAX_UPLOAD_BYTES."""
    _reject_nonportable_upload_name(str(getattr(uploaded_file, "name", "") or ""))
    data = uploaded_file.getvalue()
    if len(data) > MAX_UPLOAD_BYTES:
        raise ValueError(
            "Uploaded file exceeds 1 MiB. Use a smaller FASTA or paste a "
            "region of interest."
        )
    if len(data) == 0:
        raise ValueError("Uploaded file is empty. Nothing was analyzed.")
    return data.decode("utf-8", errors="replace")


def _enforce_raw_input_size(text: str) -> str:
    """Recusa texto colado acima do limite tecnico de caracteres."""
    payload = text or ""
    if len(payload) > dna_analysis.MAX_RAW_INPUT_CHARS:
        raise ValueError(
            f"Input exceeds the technical limit of "
            f"{dna_analysis.MAX_RAW_INPUT_CHARS:,} characters. Paste a gene, "
            "transcript or contig region, not a chromosome-scale sequence."
        )
    return payload


def _ncbi_allow_request() -> bool:
    """Aplica o intervalo minimo entre consultas NCBI na sessao atual."""
    now = time.monotonic()
    last = float(st.session_state.get("_ncbi_last_monotonic", 0.0))
    if last and (now - last) < NCBI_MIN_INTERVAL_S:
        remaining = NCBI_MIN_INTERVAL_S - (now - last)
        st.warning(
            f"Wait {remaining:.1f}s before the next NCBI request (rate policy)."
        )
        return False
    st.session_state["_ncbi_last_monotonic"] = now
    return True


def _resolve_sequence(text_value: str, uploaded_file) -> str:
    """Resolve a sequencia de entrada a partir de texto ou arquivo FASTA.

    Args:
        text_value: Conteudo colado pelo usuario na area de texto.
        uploaded_file: Arquivo carregado (.fasta/.fa) ou None.

    Returns:
        Sequencia bruta como string. FASTA de um registro devolve os residuos.
        FASTA de varios registros e recusado.

    Raises:
        ValueError: Se o arquivo ou o texto excederem os limites tecnicos, ou se
            o FASTA contiver mais de um registro.

    Nota:
        Prioriza o arquivo sobre o texto quando ambos estao presentes.
    """
    if uploaded_file is not None:
        content = _read_uploaded_text(uploaded_file)
        parsed = dna_analysis.parse_sequence_payload(content)
        return str(parsed["sequence"])
    parsed = dna_analysis.parse_sequence_payload(text_value or "")
    return str(parsed["sequence"])


def _resolve_raw_text(text_value: str, uploaded_file) -> str:
    """Retorna o texto de entrada preservando cabecalhos FASTA (uso interno).

    Args:
        text_value: Conteudo colado pelo usuario na area de texto.
        uploaded_file: Arquivo carregado (.fasta/.fa) ou None.

    Returns:
        Texto bruto, sem remover a linha de descricao do FASTA.

    Raises:
        ValueError: Se o arquivo ou o texto excederem os limites tecnicos.

    Nota:
        O cabecalho e necessario para detectar numeros de acesso, que dao acesso a
        anotacao oficial de regioes codificadoras.
    """
    if uploaded_file is not None:
        content = _read_uploaded_text(uploaded_file)
    else:
        content = _enforce_raw_input_size(text_value or "")
    if (content or "").strip().startswith(">"):
        dna_analysis.parse_sequence_payload(content)
    return content


def _normalized_seq(seq: str) -> str:
    """Remove espacos e padroniza maiusculas para comparar entradas (interno).

    Args:
        seq: Sequencia bruta, possivelmente com espacos ou letras minusculas.

    Returns:
        String continua em maiusculas. String vazia quando seq e vazia.

    Raises:
        Nenhum.
    """
    return "".join((seq or "").split()).upper()


def _format_percent(value: float) -> str:
    """Formata um percentual; NaN e indisponivel, nao zero (interno).

    Args:
        value: Percentual calculado, ou NaN quando o denominador e zero.

    Returns:
        Texto com duas casas ou "Unavailable".

    Raises:
        Nenhum.
    """
    if isinstance(value, float) and math.isnan(value):
        return "Unavailable"
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return "N/A"


def _store_workspace(
    molecule: str,
    sequence: str,
    source: str,
    identifier: str = "",
    accession: str = "",
    version: str = "",
    organism: str = "",
) -> None:
    """Guarda uma unica copia da sequencia ativa e o hash de proveniencia."""
    key = str(molecule or "").strip().lower()
    st.session_state[f"workspace_{key}"] = sequence
    st.session_state[f"workspace_{key}_meta"] = provenance.workspace_snapshot(
        sequence=sequence,
        molecule=molecule,
        source=source,
        identifier=identifier,
        accession=accession,
        version=version,
        organism=organism,
    )


def _export_frame(rows: list) -> pd.DataFrame:
    """DataFrame de exportacao com NaN/None escritos como N/A, nao celula vazia."""
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return frame.map(provenance.csv_cell)


def _render_coding_transparency(resolution: dict, notes: list) -> None:
    """Exibe metodo, contagens e alertas da resolucao de CDS (uso interno).

    Args:
        resolution: Dicionario retornado por rna_analysis.resolve_coding_codons().
        notes: Avisos sobre mecanismos especiais de traducao.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.

    Nota:
        A declaracao do metodo, do numero de regioes e do total de codons e
        obrigatoria para que o usuario saiba se olha anotacao oficial ou predicao.
    """
    render_metric_grid(
        [
            ("Coding Regions", f"{resolution['n_cds']:,}", ""),
            ("Codons Counted", f"{resolution['n_codons']:,}", ""),
            ("Coding Coverage", f"{resolution['coverage_pct']:.2f}", "%"),
            ("Unfiltered Count", f"{resolution['naive_codon_count']:,}", ""),
        ]
    )

    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>Unfiltered Count is the "
        "sequence length divided by three, shown only as a reference for what a "
        "filterless count would produce. It is not a valid codon total.</span>",
        unsafe_allow_html=True,
    )

    for message in resolution["warnings"]:
        st.warning(message)
    for note in notes:
        st.info(note)


def _render_cds_resolver(rna_seq: str) -> Optional[dict]:
    """Renderiza os controles de resolucao de CDS e devolve o resultado (interno).

    Args:
        rna_seq: Sequencia de mRNA ja validada e normalizada.

    Returns:
        Dicionario produzido por rna_analysis.resolve_coding_codons(), ou None
        quando nenhuma regiao codificadora pode ser resolvida com os criterios
        escolhidos.

    Raises:
        Nenhum.

    Nota:
        Nao existe recurso silencioso a contagem da sequencia inteira: essa leitura
        bruta so ocorre quando o usuario a seleciona explicitamente.
    """
    raw_text = str(st.session_state.get("rna_raw_text", ""))
    accession = rna_analysis.detect_accession(raw_text)

    with st.container(border=True):
        st.markdown(
            '<p class="module-title" style="font-size:1rem;">Coding Region '
            "Resolution</p>"
            '<p class="module-desc">Codon counts always come from resolved coding '
            "regions. Official annotation takes precedence over ORF prediction, and "
            "the whole sequence is never counted unless explicitly selected.</p>",
            unsafe_allow_html=True,
        )

        src_col, len_col = st.columns([2, 1])
        with src_col:
            source = st.selectbox(
                "CDS source",
                options=["auto", "annotation", "prediction", "whole"],
                key="rna_cds_source",
            )
        with len_col:
            min_cds = st.number_input(
                "Minimum CDS length (nt)",
                min_value=30,
                max_value=3000,
                value=dna_analysis.MIN_CDS_LENGTH_NT,
                step=30,
                key="rna_min_cds",
            )

        email_col, opt_col = st.columns([2, 1])
        with email_col:
            email = st.text_input(
                "NCBI email (required to fetch official annotation)",
                key="rna_ncbi_email",
                placeholder="your.email@university.edu",
            )
        with opt_col:
            drop_overlapping = st.checkbox(
                "Discard nested ORFs", value=True, key="rna_drop_overlapping"
            )

        cds_regions: list = []
        organism = ""
        if source in {"auto", "annotation"}:
            if not accession:
                st.info(
                    "No GenBank or RefSeq accession detected in the input. Keep the "
                    "FASTA header or paste the accession to enable official "
                    "annotation."
                )
            elif not email.strip():
                st.info(
                    f"Accession {accession} detected. Provide an NCBI email above "
                    "to fetch the official CDS annotation."
                )
            else:
                identifier_kind = ncbi_fetch.classify_identifier(accession)["kind"]
                if identifier_kind in {
                    "protein",
                    "gene_id",
                    "assembly",
                    "organism_name",
                }:
                    st.warning(
                        f"Detected identifier {accession} is {identifier_kind}, not a "
                        "nucleotide accession. Official CDS annotation is fetched "
                        "only from nucleotide records."
                    )
                else:
                    try:
                        with st.spinner(f"Fetching annotation for {accession}..."):
                            bundle = _annotation_bundle(
                                accession, "nucleotide", email.strip()
                            )
                        cds_regions = bundle["regions"]
                        organism = bundle["organism"]
                        if cds_regions:
                            st.success(
                                f"Official annotation for {bundle['accession']} "
                                f"({organism or 'unknown organism'}): "
                                f"{len(cds_regions):,} CDS features."
                            )
                        else:
                            st.warning(
                                f"Record {bundle['accession']} has no annotated CDS "
                                "features."
                            )
                    except (ValueError, RuntimeError, ncbi_fetch.NCBIQueryError) as exc:
                        st.warning(f"NCBI annotation unavailable: {exc}")

        if source == "annotation" and not cds_regions:
            st.error(
                "Official annotation was required but no CDS could be retrieved. "
                "Switch to ORF prediction or provide a valid accession and email."
            )
            return None

        try:
            resolution = _resolve_coding_codons(
                rna_seq,
                cds_regions,
                source,
                int(min_cds),
                bool(drop_overlapping),
                source == "whole",
            )
        except ValueError as exc:
            st.error(str(exc))
            return None

        product_names = [
            str(region.get("product") or region.get("gene") or "")
            for region in resolution["regions"]
        ]
        notes = rna_analysis.special_genome_notes(accession, organism, product_names)
        _render_coding_transparency(resolution, notes)

    return resolution


def _render_alerts(alerts: list) -> None:
    """Exibe os destaques de valores atipicos encontrados (uso interno).

    Args:
        alerts: Frases explicando cada valor atipico e seu significado biologico.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.

    Nota:
        Quando nenhum valor se destaca, a ausencia de alerta tambem e informacao e
        e declarada explicitamente.
    """
    with st.container(border=True):
        st.markdown(
            '<p class="module-title" style="font-size:1rem;">Atypical Values</p>',
            unsafe_allow_html=True,
        )
        if not alerts:
            st.markdown(
                "<span style='color:var(--hs-text-secondary);'>No metric fell outside its typical "
                "range for the thresholds used.</span>",
                unsafe_allow_html=True,
            )
            return
        for alert in alerts:
            st.warning(alert)


def _render_summary_table(rows: list, key: str, file_name: str) -> None:
    """Renderiza a tabela-resumo dos indices calculados (uso interno).

    Args:
        rows: Lista de dicts com as chaves "Metric", "Value" e "Unit".
        key: Chave unica do botao de download.
        file_name: Nome do arquivo CSV oferecido para download.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.

    Nota:
        Concentrar todos os indices em uma tabela unica permite comparar e exportar
        o resultado sem percorrer cada secao.
    """
    if not rows:
        return
    with st.container(border=True):
        st.markdown(
            '<p class="module-title" style="font-size:1rem;">Summary of Indices</p>',
            unsafe_allow_html=True,
        )
        display_rows = [
            {
                "Metric": str(row.get("Metric") or ""),
                "Value": str(provenance.csv_cell(row.get("Value"))),
                "Unit": str(row.get("Unit") or ""),
            }
            for row in rows
        ]
        frame = pd.DataFrame(display_rows, columns=["Metric", "Value", "Unit"])
        st.dataframe(frame, width="stretch", hide_index=True)
        st.download_button(
            "Download summary (CSV)",
            data=_export_frame(rows).to_csv(index=False).encode("utf-8"),
            file_name=file_name,
            mime="text/csv",
            key=key,
        )


def _nucleotide_alerts(
    gc: float,
    dinucleotides: dict,
    enc: Optional[float] = None,
    cai: Optional[float] = None,
    cai_organism: str = "",
    resolution: Optional[dict] = None,
) -> list:
    """Detecta valores atipicos em uma analise de DNA ou RNA (uso interno).

    Args:
        gc: Conteudo GC em porcentagem.
        dinucleotides: Saida de dna_analysis.dinucleotide_frequencies().
        enc: Numero Efetivo de Codons, quando disponivel.
        cai: Codon Adaptation Index, quando disponivel.
        cai_organism: Organismo de referencia usado no CAI.
        resolution: Resultado de rna_analysis.resolve_coding_codons(), quando
            disponivel.

    Returns:
        Lista de frases descrevendo cada desvio e seu significado biologico.

    Raises:
        Nenhum.

    Nota:
        Os limiares sao convencoes de leitura, nao testes estatisticos; servem para
        chamar atencao, nao para concluir.
    """
    alerts: list = []

    if isinstance(gc, float) and math.isnan(gc):
        alerts.append(
            "GC is unavailable: the sequence has no canonical A, C, G or T bases. "
            "This is not a GC of 0%."
        )
    elif gc < 30.0:
        alerts.append(
            f"GC de {gc:.2f}% e baixo. Genomas AT-ricos ocorrem em organismos como "
            "Plasmodium e muitos bacteriofagos; na pratica, reduzem a temperatura "
            "de melting e dificultam o desenho de primers especificos."
        )
    elif gc > 65.0:
        alerts.append(
            f"GC de {gc:.2f}% e alto, padrao de Actinobacteria e de muitos genomas "
            "termofilicos. GC elevado estabiliza estrutura secundaria e costuma "
            "exigir aditivos como DMSO em PCR."
        )

    cpg = dinucleotides.get("CG") if dinucleotides else None
    if isinstance(cpg, dict):
        try:
            odds = float(cpg.get("observed_expected"))
        except (TypeError, ValueError):
            odds = float("nan")
        if not math.isnan(odds) and 0.0 < odds < 0.4:
            alerts.append(
                f"CpG observado/esperado de {odds:.2f} indica forte deplecao do "
                "dinucleotideo CG em relacao ao produto das frequencias de C e G. "
                "Em genomas de vertebrados isso costuma refletir metilacao e "
                "desaminacao C para T. O mesmo numero tambem aparece em sequencias "
                "sinteticas ou muito repetitivas que simplesmente nao colocam C ao "
                "lado de G; nao prova origem de vertebrado."
            )
        elif not math.isnan(odds) and odds >= 0.9:
            alerts.append(
                f"CpG observado/esperado de {odds:.2f} esta proximo do esperado "
                "por azar. Isso aponta ausencia de deplecao por metilacao, "
                "compativel com ilha CpG hipometilada de promotor, com genoma "
                "bacteriano ou de invertebrado, ou com genoma viral."
            )

    if enc is not None:
        if enc < 35.0:
            alerts.append(
                f"ENC of {enc:.2f} indicates strong codon usage bias in the "
                "resolved coding codons. This is a compositional statistic "
                "(Wright 1990), not proof of high expression."
            )
        elif enc > 55.0:
            alerts.append(
                f"ENC of {enc:.2f} indicates nearly uniform synonymous codon use "
                "in this region. It is not a statement about transcription."
            )

    if cai is not None:
        reference = cai_organism or "selected reference table"
        if cai < 0.5:
            alerts.append(
                f"CAI of {cai:.4f} versus the {reference} table is low. CAI is "
                "relative to that embedded codon-usage table, not an experimental "
                "expression measurement."
            )
        elif cai > 0.8:
            alerts.append(
                f"CAI of {cai:.4f} versus the {reference} table is high. That "
                "describes similarity to the selected reference codon usage, not "
                "validated expression in a host."
            )

    if resolution is not None:
        coverage = float(resolution.get("coverage_pct", 0.0))
        method = str(resolution.get("method", ""))
        if method == rna_analysis.CDS_SOURCE_PREDICTION and coverage > 80.0:
            alerts.append(
                f"Cobertura codificante de {coverage:.2f}% obtida por predicao. "
                "Genomas compactos como os virais chegam a esse patamar, mas em "
                "sequencias eucarioticas isso normalmente indica ORFs previstas "
                "por azar em regiao nao codificante."
            )

    return alerts


def _protein_alerts(
    gravy: float,
    instability: float,
    pi: float,
    aliphatic: float,
    categories: dict,
    composition: dict,
) -> list:
    """Detecta valores atipicos em uma analise de proteina (uso interno).

    Args:
        gravy: Valor de GRAVY.
        instability: Indice de instabilidade.
        pi: Ponto isoeletrico.
        aliphatic: Indice alifatico.
        categories: Saida de protein_analysis.amino_acid_categories().
        composition: Saida de protein_analysis.amino_acid_composition().

    Returns:
        Lista de frases descrevendo cada desvio e seu significado biologico.

    Raises:
        Nenhum.

    Nota:
        Os limiares seguem convencoes de leitura da literatura de ProtParam e nao
        substituem validacao experimental.
    """
    alerts: list = []

    if gravy > 0.5:
        alerts.append(
            f"GRAVY of {gravy:.3f} (Kyte-Doolittle mean) is positive: the chain "
            "is hydrophobic on this scale. That is not a membrane-protein "
            "assignment and is not a solubility assay."
        )
    elif gravy < -1.0:
        alerts.append(
            f"GRAVY of {gravy:.3f} (Kyte-Doolittle mean) is strongly negative: "
            "the chain is hydrophilic on this scale. That is not a disordered-"
            "region or ribosomal-protein assignment."
        )

    if instability >= 40.0:
        alerts.append(
            f"Instability index of {instability:.2f} is the Guruprasad "
            "computational index (threshold 40). It is not experimental evidence "
            "of instability in vitro or in the cell."
        )

    if pi < 5.0:
        alerts.append(
            f"Theoretical pI of {pi:.2f} is acidic on the ProtParam/Bjellqvist "
            "scale. It is not an experimentally measured isoelectric point."
        )
    elif pi > 9.0:
        alerts.append(
            f"Theoretical pI of {pi:.2f} is basic on the ProtParam/Bjellqvist "
            "scale. It is not evidence of DNA/RNA binding."
        )

    if aliphatic > 100.0:
        alerts.append(
            f"Aliphatic index of {aliphatic:.2f} (Ikai 1980) is high. It is a "
            "composition statistic, not a thermostability measurement."
        )

    charged = categories.get("charged_total", {})
    if isinstance(charged, dict) and float(charged.get("frequency", 0.0)) > 45.0:
        alerts.append(
            f"{charged.get('frequency', 0.0):.2f}% dos residuos sao carregados. "
            "Fracoes altas de residuos carregados, combinadas com baixa "
            "hidrofobicidade, sao caracteristicas de proteinas intrinsecamente "
            "desordenadas, que nao adotam estrutura estavel isoladamente."
        )

    cysteine = composition.get("C", {})
    if isinstance(cysteine, dict) and float(cysteine.get("frequency", 0.0)) > 5.0:
        alerts.append(
            f"Cisteina representa {cysteine.get('frequency', 0.0):.2f}% dos "
            "residuos, bem acima da media de cerca de 2%. Enriquecimento em "
            "cisteina sugere multiplas pontes dissulfeto e, portanto, proteina "
            "secretada ou extracelular."
        )

    return alerts


def _type_badge(info: dict) -> str:
    """Monta os badges de tipo detectado e validade (uso interno)."""
    type_variant = {"DNA": "dna", "RNA": "rna", "PROTEIN": "prot"}.get(
        info["type"], "warn"
    )
    valid_variant = "dna" if info["is_valid"] else "err"
    valid_text = "Valid" if info["is_valid"] else "Invalid"
    return (
        f'{badge("Type: " + info["type"], type_variant)} '
        f'{badge(valid_text, valid_variant)}'
    )


def _nucleic_highlight_indices(sequence: str, molecule: str) -> list[int]:
    """Motif/CRISPR spans on this sequence. Empty if hashes do not match."""
    digest = provenance.sequence_digest(sequence)
    indices: list[int] = []
    motif = st.session_state.get("motif_result")
    if (
        isinstance(motif, dict)
        and str(motif.get("molecule") or "") == molecule
        and str(motif.get("sequence_hash") or "") == digest
    ):
        for hit in list(motif.get("hits") or []):
            start = int(hit.get("start") or 0)
            end = int(hit.get("end") or 0)
            if 0 <= start < end <= len(sequence):
                indices.extend(range(start, end))
    selection = st.session_state.get("molecule_selection")
    if isinstance(selection, dict) and molecule_selection.selection_applies_to(
        selection,
        sequence_hash=digest,
        molecule=molecule,
    ):
        indices.extend(molecule_selection.query_indices_of(selection))
    feature_source = st.session_state.get("ncbi_feature_source") or {}
    spans = st.session_state.get("ncbi_feature_spans")
    if (
        molecule in {"DNA", "RNA"}
        and isinstance(spans, list)
        and str(feature_source.get("sequence_hash") or "") == digest
    ):
        for span in spans:
            if not isinstance(span, dict):
                continue
            try:
                start = int(span.get("start"))
                end = int(span.get("end"))
            except (TypeError, ValueError):
                continue
            if 0 <= start < end <= len(sequence):
                indices.extend(range(start, end))
    packed = st.session_state.get("crispr_result")
    if (
        molecule == "DNA"
        and isinstance(packed, dict)
        and str(packed.get("sequence_hash") or "") == digest
    ):
        for guide in list(packed.get("guides") or []):
            start = guide.get("start_0based")
            if start is None and guide.get("position"):
                start = int(guide["position"]) - 1
            length = len(str(guide.get("guide_sequence") or ""))
            if start is not None and length:
                end = int(start) + int(length)
                if 0 <= int(start) < end <= len(sequence):
                    indices.extend(range(int(start), end))
            pam_start = guide.get("pam_start_0based")
            pam_seq = str(guide.get("pam_sequence") or "")
            if pam_start is not None and pam_seq:
                pam_end = int(pam_start) + len(pam_seq)
                if 0 <= int(pam_start) < pam_end <= len(sequence):
                    indices.extend(range(int(pam_start), pam_end))
    return sorted(set(indices))


def _lod_select(key: str, default: str = region_nav.LOD_HIGH) -> str:
    """Selector de LOD visual. Nao altera coordenadas cientificas."""
    options = [region_nav.LOD_LOW, region_nav.LOD_MEDIUM, region_nav.LOD_HIGH]
    index = options.index(default) if default in options else 2
    return str(
        st.selectbox(
            "Visual level of detail",
            options,
            index=index,
            format_func=lambda item: {
                region_nav.LOD_LOW: "Low (stride 4, shorter helix window)",
                region_nav.LOD_MEDIUM: "Medium (stride 2)",
                region_nav.LOD_HIGH: "High (all mapped points)",
            }[item],
            key=key,
            help=(
                "Visual only. Omitted residues keep deposited or illustrative "
                "coordinates in the scientific envelope. Points are never interpolated."
            ),
        )
    )


def _consume_3d_pending_position(widget_key: str, length: int) -> None:
    """Aplica pick 3D ao widget de posicao no rerun seguinte. Uso interno."""
    pending_key = f"{widget_key}_from_3d"
    pending = st.session_state.pop(pending_key, None)
    if isinstance(pending, int) and length > 0 and 0 <= pending < length:
        st.session_state[widget_key] = pending


def _apply_plotly_pick_to_sequence(
    event: object,
    *,
    sequence: str,
    live_hash: str,
    scene_hash: str,
    selected: int,
    position_key: str,
    range_start_key: str = "",
    range_end_key: str = "",
) -> None:
    """3D click -> posicao da sequencia se hashes coincidirem. Uso interno."""
    parsed_fn = runtime.require_attr(structure_viewer, "parse_plotly_selection")
    indices_fn = runtime.require_attr(structure_viewer, "parse_plotly_selection_indices")
    if parsed_fn is None or indices_fn is None:
        st.error(
            "3D selection API is missing in this process (stale ui.structure_viewer). "
            "Stop Streamlit and run: python -m streamlit run app.py"
        )
        return
    parsed = parsed_fn(event)
    indices = list(indices_fn(event))
    if not parsed and not indices:
        return
    if live_hash and scene_hash and str(live_hash) != str(scene_hash):
        st.markdown(status_badge("STALE"), unsafe_allow_html=True)
        st.warning("stale/incompatible structure")
        return
    if range_start_key and range_end_key and len(indices) >= 2 and sequence:
        start = max(0, min(indices))
        end = min(len(sequence), max(indices) + 1)
        st.caption(
            f"3D multi-point pick → helix range [{start}, {end}). "
            f"sequence_hash {str(live_hash or '')[:16]}. Range is a view window, "
            "not a biological feature annotation."
        )
        if (
            int(st.session_state.get(range_start_key) or -1) != start
            or int(st.session_state.get(range_end_key) or -1) != end
        ):
            st.session_state[range_start_key] = start
            st.session_state[range_end_key] = end
            st.rerun()
        return
    idx = parsed.get("query_index_0based") if parsed else None
    if idx is None:
        st.info("No mapped sequence position for this 3D point.")
        return
    index = int(idx)
    if not sequence or index < 0 or index >= len(sequence):
        st.info("No mapped sequence position for this 3D point.")
        return
    st.caption(
        f"3D pick → position {index} identity {sequence[index]} "
        f"chain {parsed.get('chain_id') or 'N/A'} "
        f"atom {parsed.get('atom_name') or 'N/A'} "
        f"sequence_hash {str(live_hash or '')[:16]}. "
        "Selection marker is visual, not a scientific score."
    )
    if index != int(selected):
        st.session_state[f"{position_key}_from_3d"] = index
        st.rerun()


def _structure_plotly_chart(
    figure: Any,
    *,
    key: str,
    on_select: str = "ignore",
    height: Optional[int] = None,
) -> Any:
    """Renderiza Scatter3d com altura fixa e sem o tema Streamlit do Plotly.

    Streamlit 1.58 aplica theme='streamlit' e height='content' por omissao.
    O tema substitui cores por placeholders 2D; Scatter3d/gl3d no Firefox
    fica em branco. height='content' numa coluna tambem pode colapsar o canvas.
    """
    kwargs: dict[str, Any] = {
        "width": "stretch",
        "height": int(height) if height is not None else 640,
        "theme": None,
        "key": key,
        "config": {
            "scrollZoom": True,
            "displaylogo": False,
            "responsive": True,
        },
    }
    if on_select != "ignore":
        kwargs["on_select"] = on_select
        kwargs["selection_mode"] = "points"
    with st.container(key=f"{key}_host"):
        return st.plotly_chart(figure, **kwargs)


def _sequence_viewer_panel(
    seq: str,
    seq_type: str,
    *,
    key_prefix: str,
    highlight_indices: Optional[list[int]] = None,
) -> None:
    """Viewer HTML com salto e janela. Metricas usam a sequencia completa."""
    length = len(seq)
    cap = int(scale_profile.MAX_SEQUENCE_VIEWER_RESIDUES)
    jump_col, start_col, end_col, nav_col = st.columns(4)
    with jump_col:
        jump = int(
            st.number_input(
                "Jump (0-based)",
                min_value=0,
                max_value=max(0, length - 1),
                value=0,
                step=1,
                key=f"{key_prefix}_jump",
            )
        )
    with start_col:
        start = int(
            st.number_input(
                "Start (0-based)",
                min_value=0,
                max_value=max(0, length - 1),
                value=0,
                step=1,
                key=f"{key_prefix}_vstart",
            )
        )
    with end_col:
        end = int(
            st.number_input(
                "End (exclusive)",
                min_value=1,
                max_value=max(1, length),
                value=min(length, cap),
                step=1,
                key=f"{key_prefix}_vend",
            )
        )
    with nav_col:
        if st.button("Center on jump", key=f"{key_prefix}_center"):
            half = cap // 2
            st.session_state[f"{key_prefix}_vstart"] = max(0, int(jump) - half)
            st.session_state[f"{key_prefix}_vend"] = min(
                length, max(0, int(jump) - half) + cap
            )
            st.rerun()
    try:
        start, end = region_nav.clamp_region(length, start, end)
    except region_nav.RegionError as exc:
        st.info(str(exc))
        return
    st.caption(
        f"Viewer window [{start}, {end}) of {length} residues. HTML display is "
        f"capped at {cap:,} residues per window. Metrics use the full sequence."
    )
    st.markdown(
        sequence_display(
            seq,
            seq_type,
            highlight_indices=highlight_indices,
            view_start=start,
            view_end=end,
        ),
        unsafe_allow_html=True,
    )


def _render_scene_inspector(scene: dict, sequence: str, selected: int) -> None:
    """Inspector seccionado a partir do scene model. Sem interpretacao extra."""
    residue = ""
    if sequence and 0 <= int(selected) < len(sequence):
        residue = sequence[int(selected)]
    sections = {
        "IDENTITY": [
            ("Structure ID", scene.get("structure_id")),
            ("Kind", scene.get("kind_label") or scene.get("kind")),
            ("Selected residue", residue or "N/A"),
            ("Selected index (0-based)", selected),
        ],
        "SOURCE": [
            ("Source", scene.get("source")),
            ("Status", scene.get("status")),
        ],
        "STRUCTURE": [
            ("Chain", scene.get("chain_id")),
            ("Model", scene.get("model_number")),
            ("Representation", scene.get("representation")),
            ("Visual LOD", scene.get("visual_lod")),
            ("Visual stride", scene.get("visual_stride")),
            ("Partial view", scene.get("partial_view")),
            ("Backbone points drawn", scene.get("n_backbone_points")),
            ("Atoms rendered", scene.get("n_atoms_rendered")),
        ],
        "MAPPING": [
            ("Mapping status", scene.get("mapping_status")),
            ("Coverage", scene.get("coverage_query_with_coordinates")),
            ("View scope", scene.get("view_label") or scene.get("view_scope")),
        ],
        "ANALYSIS": [
            ("Color mode", scene.get("color_mode")),
            ("LOD note", scene.get("lod_note") or "N/A"),
        ],
        "PROVENANCE": [
            ("sequence_hash", scene.get("sequence_hash")),
            ("structure_hash", scene.get("structure_hash")),
        ],
    }
    st.markdown(inspector_panel(sections), unsafe_allow_html=True)


def _render_validated_scene_3d(
    sequence: str,
    envelope: dict,
    *,
    key_prefix: str,
    selected: int = 0,
    highlight_indices: Optional[list[int]] = None,
    position_key: str = "",
    lod: Optional[str] = None,
    range_start_key: str = "",
    range_end_key: str = "",
) -> None:
    """Renderer 3D lazy a partir de um envelope ja validado."""
    if not st.toggle("View 3D Structure", value=False, key=f"{key_prefix}_view_3d"):
        st.caption(
            "3D is not loaded until View 3D Structure is enabled. "
            "Metadata and mapping remain available as text."
        )
        return
    live_hash = provenance.sequence_digest(sequence) if sequence else ""
    stored_hash = str(envelope.get("sequence_hash") or "")
    if live_hash and stored_hash and stored_hash != live_hash:
        st.markdown(status_badge("STALE"), unsafe_allow_html=True)
        st.warning("stale/incompatible structure")
        return
    try:
        from modules import structure_scene
        from ui import structure_viewer
    except Exception:
        st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
        st.error("3D unavailable in this browser/environment.")
        return
    kind = str(envelope.get("kind") or "")
    molecule = str(envelope.get("molecule_type") or "").strip().upper()
    color_mode = "residue_type" if molecule in {"DNA", "RNA"} else "chain"
    statuses = [kind.upper() or "UNAVAILABLE"]
    if envelope.get("partial_view"):
        statuses.append("PARTIAL_VIEW")
    st.markdown(badge_row(statuses), unsafe_allow_html=True)
    st.markdown(
        meta_grid(
            [
                ("Kind", envelope.get("kind_label") or protein_structure.structure_kind_label(kind)),
                ("Structure ID", envelope.get("structure_id")),
                ("Source", envelope.get("source")),
                (
                    "View",
                    (
                        f"{envelope.get('view_start')}–{envelope.get('view_end')}"
                        if envelope.get("view_start") is not None
                        else envelope.get("view_scope")
                    ),
                ),
                ("sequence_hash", (envelope.get("sequence_hash") or "")[:16]),
                ("method", (envelope.get("metadata") or {}).get("method")),
                ("resolution_A", (envelope.get("metadata") or {}).get("resolution_angstrom")),
                ("molecule", envelope.get("molecule_type")),
            ]
        ),
        unsafe_allow_html=True,
    )
    disclaimer = str(
        envelope.get("disclaimer") or (envelope.get("metadata") or {}).get("disclaimer") or ""
    )
    if disclaimer:
        st.caption(disclaimer)
    if kind == "illustrative":
        st.caption(
            "ILLUSTRATIVE helical geometry. Sequence order and base identity are "
            "real. XYZ coordinates are not an experimental structure and not a "
            "predicted fold."
        )
    elif kind == "experimental":
        st.caption("Experimental structure. This is not an exact molecular structure claim.")
    elif kind == "predicted":
        st.caption("Predicted Structure. This is not an experimental measurement.")
    if lod is None:
        lod = _lod_select(f"{key_prefix}_lod")
    else:
        st.caption(
            f"Visual LOD {lod}. Omitted residues keep deposited or illustrative "
            "coordinates in the scientific envelope. Points are never interpolated."
        )
    try:
        scene = structure_scene.build_scene(
            envelope,
            expected_sequence_hash=stored_hash or live_hash,
            representation="backbone",
            color_mode=color_mode,
            selected_query_index=int(selected) if sequence else None,
            highlight_query_indices=highlight_indices or [],
            lod=lod,
        )
    except Exception as exc:
        st.error("3D viewport could not be built from the structure envelope.")
        st.caption(str(exc))
        return
    scene_status = str(scene.get("status") or "")
    if scene_status != "READY":
        st.markdown(status_badge(scene_status), unsafe_allow_html=True)
        st.error(str(scene.get("message") or "Coordinates Unavailable"))
        if scene.get("webgl_note"):
            st.caption(str(scene.get("webgl_note")))
        return
    try:
        figure = structure_viewer.figure_from_scene(scene)
    except Exception as exc:
        if "insufficient" in str(exc).lower():
            st.markdown(status_badge("INSUFFICIENT_DATA"), unsafe_allow_html=True)
            st.error("3D viewport has no finite coordinates to draw.")
        else:
            st.error("3D viewport could not be built from the scene model.")
        st.caption(str(exc))
        return
    if int(scene.get("n_backbone_points") or 0) <= 0 and int(scene.get("n_atoms_rendered") or 0) <= 0:
        st.markdown(status_badge("INSUFFICIENT_DATA"), unsafe_allow_html=True)
        st.error("3D viewport has no finite coordinates to draw.")
        return
    n_points = int(scene.get("n_backbone_points") or 0)
    n_pairs = int(scene.get("n_pair_traces") or 0)
    st.caption(
        f"3D viewport (Plotly Scatter3d). {n_points} backbone points"
        + (f", {n_pairs} base pairs" if n_pairs else "")
        + f", kind={scene.get('kind')}. "
        "Identities in hover come from the mapped sequence. "
        "Not Mol*, not a molecular surface, not DSSP cartoon, and not molecular dynamics. "
        "Rotate, pan and zoom inside the canvas. "
        "Click a mapped backbone point to move the sequence cursor when mapping exists."
    )
    expanded = st.toggle(
        "Expanded viewport",
        value=False,
        key=f"{key_prefix}_expanded_3d",
        help="Taller canvas. Does not change coordinates.",
    )
    from ui.tokens import STRUCTURE_VIEWPORT_HEIGHT, STRUCTURE_VIEWPORT_HEIGHT_EXPANDED

    chart_height = STRUCTURE_VIEWPORT_HEIGHT_EXPANDED if expanded else STRUCTURE_VIEWPORT_HEIGHT
    try:
        figure.update_layout(height=chart_height)
    except Exception:
        pass
    nonce = int(st.session_state.get(f"{key_prefix}_cam_nonce") or 0)
    plot_key = f"{key_prefix}_plotly" if nonce == 0 else f"{key_prefix}_plotly_{nonce}"
    reset_col, hand_col = st.columns([1, 3])
    with reset_col:
        if st.button("Reset View", key=f"{key_prefix}_reset_view"):
            st.session_state[f"{key_prefix}_cam_nonce"] = nonce + 1
            st.rerun()
    with hand_col:
        render_hand_control(
            plot_key=plot_key,
            structure_id=str(
                scene.get("structure_id") or scene.get("content_hash") or key_prefix
            ),
            illustrative=str(scene.get("kind") or "").lower() == "illustrative",
            control_key=str(key_prefix),
        )
    st.markdown('<div class="hs-viewport">', unsafe_allow_html=True)
    try:
        if expanded:
            event = _structure_plotly_chart(
                figure,
                key=plot_key,
                on_select="rerun" if position_key else "ignore",
                height=chart_height,
            )
            _render_scene_inspector(scene, sequence, int(selected) if sequence else 0)
        else:
            main_col, insp_col = helix_shell.split_workspace()
            with main_col:
                event = _structure_plotly_chart(
                    figure,
                    key=plot_key,
                    on_select="rerun" if position_key else "ignore",
                    height=chart_height,
                )
            with insp_col:
                _render_scene_inspector(scene, sequence, int(selected) if sequence else 0)
    except Exception as exc:
        st.error("3D unavailable in this browser/environment.")
        st.caption(str(exc) or str(scene.get("webgl_note") or ""))
        return
    st.markdown("</div>", unsafe_allow_html=True)
    if position_key:
        _apply_plotly_pick_to_sequence(
            event,
            sequence=sequence,
            live_hash=live_hash,
            scene_hash=str(scene.get("sequence_hash") or stored_hash or ""),
            selected=int(selected) if sequence else 0,
            position_key=position_key,
            range_start_key=range_start_key,
            range_end_key=range_end_key,
        )
    if scene.get("partial_view"):
        st.markdown(status_badge("PARTIAL_VIEW"), unsafe_allow_html=True)
        st.caption("Partial structure view. This is not the complete polymer or deposited file.")
    report = structure_scene.structure_3d_report(scene, envelope)
    kind_tag = str(scene.get("kind") or "scene")
    st.download_button(
        "3D scene summary (JSON)",
        data=json.dumps(report, indent=2).encode("utf-8"),
        file_name=f"{safe_download_filename(f'helixscope_{key_prefix}_{kind_tag}_3d')}.json",
        mime="application/json",
        key=f"{key_prefix}_scene_json",
    )


def _render_polymer_structure_search(seq: str, molecule: str, key_prefix: str) -> Optional[dict]:
    """Busca RCSB para DNA/RNA sob demanda. Sem AlphaFold. Sem fetch ao abrir."""
    st.caption(
        f"Experimental {molecule} 3D requires a deposited RCSB structure. "
        "Opening this section does not search. AlphaFold DB is a protein resource."
    )
    pdb_id = st.text_input(
        "PDB ID (optional)",
        key=f"{key_prefix}_pdb",
        placeholder="1BNA",
    )
    seq_search = st.checkbox(
        f"Search RCSB PDB by this {molecule} sequence (MMseqs2 Search API)",
        key=f"{key_prefix}_seq_search",
    )
    if st.button(f"Search {molecule} structure sources", key=f"{key_prefix}_search_btn"):
        try:
            with st.spinner(f"Retrieving {molecule} structure records..."):
                found = protein_structure.search_structure_sources(
                    seq,
                    pdb_id=pdb_id,
                    search_pdb_by_sequence=seq_search,
                    lookup_uniprot=False,
                    lookup_alphafold=False,
                    molecule=molecule,
                )
        except (protein_structure.StructureError, ValueError) as exc:
            _render_structure_error(exc)
        else:
            st.session_state[f"{key_prefix}_hits"] = found
            st.session_state.pop(f"{key_prefix}_result", None)
    found = st.session_state.get(f"{key_prefix}_hits")
    if isinstance(found, dict) and str(found.get("sequence_hash") or "") != provenance.sequence_digest(seq):
        found = None
    if not isinstance(found, dict):
        loaded = st.session_state.get(f"{key_prefix}_result")
        if isinstance(loaded, dict) and str(loaded.get("sequence_hash") or "") != provenance.sequence_digest(seq):
            return None
        return loaded if isinstance(loaded, dict) else None
    st.markdown(status_badge(str(found.get("status") or "UNAVAILABLE")), unsafe_allow_html=True)
    st.caption(str(found.get("message") or ""))
    hits = list(found.get("hits") or [])
    if hits:
        rows = [
            {
                "Source": item.get("source"),
                "Kind": item.get("kind"),
                "ID": item.get("structure_id"),
                "Method": item.get("method") or "N/A",
                "Match": item.get("match_basis"),
            }
            for item in hits
        ]
        st.dataframe(_export_frame(rows), width="stretch", hide_index=True)
        labels = [
            f"{item.get('source')} | {item.get('kind')} | {item.get('structure_id')}"
            for item in hits
        ]
        pick = st.selectbox(
            "Select a retrieved record (HelixScope does not rank a best structure)",
            list(range(len(labels))),
            format_func=lambda i: labels[i],
            key=f"{key_prefix}_pick",
        )
        selected = hits[int(pick)]
        chain_id = st.text_input(
            "Chain ID (optional; empty uses the first matching polymer)",
            key=f"{key_prefix}_chain",
        )
        if st.button("Load selected structure", key=f"{key_prefix}_load"):
            try:
                with st.spinner("Retrieving and validating structure coordinates..."):
                    structure_text = ""
                    if selected.get("bundled_fixture"):
                        reader = getattr(protein_structure, "read_bundled_mmcif", None)
                        if callable(reader):
                            structure_text = reader(str(selected.get("bundled_fixture")))
                    loaded = protein_structure.load_protein_structure(
                        seq,
                        source=str(selected.get("source") or ""),
                        structure_id=str(selected.get("structure_id") or ""),
                        chain_id=str(chain_id or "").strip(),
                        metadata=selected,
                        structure_text=structure_text,
                        cache=st.session_state.setdefault(f"{key_prefix}_cache", {}),
                        molecule=molecule,
                    )
            except (protein_structure.StructureError, ValueError) as exc:
                _render_structure_error(exc)
            else:
                st.session_state[f"{key_prefix}_result"] = loaded
    loaded = st.session_state.get(f"{key_prefix}_result")
    if isinstance(loaded, dict) and str(loaded.get("sequence_hash") or "") != provenance.sequence_digest(seq):
        return None
    return loaded if isinstance(loaded, dict) else None


def _render_dna_3d_section(dna_seq: str) -> None:
    """DNA experimental (RCSB) e helice B-DNA ilustrativa, rotulos distintos."""
    st.caption(
        "A DNA sequence is not a 3D structure. Experimental coordinates come from "
        "a deposited file. Without a matching deposition the default path is an "
        "ILLUSTRATIVE canonical B-DNA helix (rise 3.38 A, 36 deg/bp, 10 bp/turn "
        "fibre parameters). It is not the cellular conformation. "
        f"Analysis of {len(dna_seq):,} nt is independent of the 3D window "
        f"(max {region_nav.MAX_HELIX_NT_HIGH} consecutive nt at high LOD). "
        "PREDICTED DNA 3D is UNAVAILABLE unless a deposited predicted envelope "
        "is attached; HelixScope does not fold DNA."
    )
    highlights = _nucleic_highlight_indices(dna_seq, "DNA")
    _consume_3d_pending_position("dna_3d_pos", len(dna_seq))
    selected = int(
        st.number_input(
            "Inspect DNA position (0-based)",
            min_value=0,
            max_value=max(0, len(dna_seq) - 1),
            value=0,
            step=1,
            key="dna_3d_pos",
        )
    )
    st.caption(f"Selected base {selected}: {dna_seq[selected] if dna_seq else 'N/A'}")
    if highlights:
        st.caption(
            f"{len(highlights)} positions highlighted from motif and/or CRISPR "
            "spans on this same sequence hash."
        )
    lod = _lod_select("dna_3d_lod")
    cap = region_nav.helix_limit_for_lod(lod)
    window_hash_key = "dna_3d_window_hash"
    digest = provenance.sequence_digest(dna_seq)
    if st.session_state.get(window_hash_key) != digest:
        st.session_state["dna_3d_rstart"] = 0
        st.session_state["dna_3d_rend"] = min(len(dna_seq), cap)
        st.session_state[window_hash_key] = digest
    region_cols = st.columns(3)
    with region_cols[0]:
        region_start = int(
            st.number_input(
                "Helix start (0-based)",
                min_value=0,
                max_value=max(0, len(dna_seq) - 1),
                step=1,
                key="dna_3d_rstart",
            )
        )
    with region_cols[1]:
        region_end = int(
            st.number_input(
                "Helix end (exclusive)",
                min_value=1,
                max_value=max(1, len(dna_seq)),
                step=1,
                key="dna_3d_rend",
            )
        )
    with region_cols[2]:
        center_on_selected = st.checkbox(
            "Center helix on selected base",
            value=False,
            key="dna_3d_center",
        )
    prefer_mode = str(
        st.selectbox(
            "DNA 3D path",
            ["auto", "illustrative", "experimental", "predicted"],
            index=0,
            key="dna_3d_prefer",
            help=(
                "auto: deposited experimental if the sequence hash matches, else "
                "ILLUSTRATIVE. predicted stays UNAVAILABLE without a deposited "
                "predicted envelope."
            ),
        )
    )
    st.caption(
        f"Atomic-detail helix draws a consecutive window of at most {cap} nt "
        f"(LOD {lod}). Omitted flanks are not interpolated. PARTIAL VIEW is shown "
        "when the window is not the full sequence. 50k nt is analysis/navigation "
        "scale, not simultaneous atomic detail."
    )
    window_nt = max(0, int(region_end) - int(region_start))
    if window_nt < 12:
        st.caption(
            f"The illustrative helix window is only {window_nt} nt "
            f"[{region_start}, {region_end}). Increase Helix end (exclusive) "
            f"up to {cap} nt at this LOD to draw more of the sequence."
        )
    deposited = _render_polymer_structure_search(dna_seq, "DNA", "dna_struct")
    prefer = None if prefer_mode == "auto" else prefer_mode
    representation = dna_3d.representation_for_sequence(
        dna_seq,
        deposited=deposited,
        prefer=prefer,
        region_start=region_start,
        region_end=region_end,
        lod=lod,
        center=selected if center_on_selected else None,
        highlight_indices=highlights,
    )
    st.session_state["dna_3d_representation"] = representation
    st.markdown(
        status_badge(str(representation.get("status") or "UNAVAILABLE")),
        unsafe_allow_html=True,
    )
    if representation.get("partial_view"):
        st.markdown(status_badge("PARTIAL_VIEW"), unsafe_allow_html=True)
        st.caption(
            "PARTIAL VIEW: contiguous window only. Omitted sequence is not drawn "
            "and not interpolated."
        )
    disclaimer = str(representation.get("disclaimer") or representation.get("reason") or "")
    if disclaimer:
        st.caption(disclaimer)
    mapped = dna_3d.position_to_3d(representation, selected)
    st.caption(
        f"Sequence position {selected} -> 3D mapping status {mapped.get('status')}."
    )
    envelope = representation.get("envelope")
    if isinstance(envelope, dict):
        _render_validated_scene_3d(
            dna_seq,
            envelope,
            key_prefix="dna_illust",
            selected=selected,
            highlight_indices=highlights,
            position_key="dna_3d_pos",
            lod=lod,
            range_start_key="dna_3d_rstart",
            range_end_key="dna_3d_rend",
        )
    elif str(representation.get("kind") or "") == dna_3d.KIND_UNAVAILABLE:
        st.info(str(representation.get("reason") or "No DNA 3D representation."))


def _render_rna_3d_section(rna_seq: str, fold_result: Optional[dict] = None) -> None:
    """RNA experimental (RCSB) e helice A-RNA ilustrativa. MFE nao vira 3D."""
    st.caption(
        "ViennaRNA MFE is secondary structure. It does not produce 3D coordinates. "
        "Dot-bracket is not converted into a physical fold. Experimental RNA 3D "
        "requires a deposited structure. The default A-RNA helix is ILLUSTRATIVE "
        "and linear; it is not FARFAR, RNAComposer or SimRNA."
    )
    highlights = _nucleic_highlight_indices(rna_seq, "RNA")
    _consume_3d_pending_position("rna_3d_pos", len(rna_seq))
    selected = int(
        st.number_input(
            "Inspect RNA 3D position (0-based)",
            min_value=0,
            max_value=max(0, len(rna_seq) - 1),
            value=0,
            step=1,
            key="rna_3d_pos",
        )
    )
    lod = _lod_select("rna_3d_lod")
    cap = region_nav.helix_limit_for_lod(lod)
    rna_window_key = "rna_3d_window_hash"
    rna_digest = provenance.sequence_digest(rna_seq)
    if st.session_state.get(rna_window_key) != rna_digest:
        st.session_state["rna_3d_rstart"] = 0
        st.session_state["rna_3d_rend"] = min(len(rna_seq), cap)
        st.session_state[rna_window_key] = rna_digest
    rna_region = st.columns(3)
    with rna_region[0]:
        rna_start = int(
            st.number_input(
                "RNA helix start (0-based)",
                min_value=0,
                max_value=max(0, len(rna_seq) - 1),
                step=1,
                key="rna_3d_rstart",
            )
        )
    with rna_region[1]:
        rna_end = int(
            st.number_input(
                "RNA helix end (exclusive)",
                min_value=1,
                max_value=max(1, len(rna_seq)),
                step=1,
                key="rna_3d_rend",
            )
        )
    with rna_region[2]:
        rna_center = st.checkbox(
            "Center helix on selected base",
            value=False,
            key="rna_3d_center",
        )
    prefer_mode = str(
        st.selectbox(
            "RNA 3D path",
            ["auto", "illustrative", "experimental", "predicted"],
            index=0,
            key="rna_3d_prefer",
            help=(
                "auto: deposited experimental if the hash matches, else "
                "ILLUSTRATIVE A-RNA. MFE never becomes experimental 3D."
            ),
        )
    )
    st.caption(
        f"Illustrative A-RNA draws at most {cap} consecutive nt. MFE/dot-bracket "
        "do not produce these coordinates."
    )
    deposited = _render_polymer_structure_search(rna_seq, "RNA", "rna_struct")
    prefer = None if prefer_mode == "auto" else prefer_mode
    representation = rna_3d.representation_for_sequence(
        rna_seq,
        deposited=deposited,
        prefer=prefer,
        region_start=rna_start,
        region_end=rna_end,
        lod=lod,
        center=selected if rna_center else None,
        highlight_indices=highlights,
    )
    st.session_state["rna_3d_representation"] = representation
    st.markdown(
        status_badge(str(representation.get("status") or "UNAVAILABLE")),
        unsafe_allow_html=True,
    )
    if representation.get("partial_view"):
        st.markdown(status_badge("PARTIAL_VIEW"), unsafe_allow_html=True)
        st.caption(
            "PARTIAL VIEW: contiguous window only. Omitted sequence is not drawn "
            "and not interpolated."
        )
    disclaimer = str(representation.get("disclaimer") or representation.get("reason") or "")
    if disclaimer:
        st.caption(disclaimer)
    fold = fold_result if isinstance(fold_result, dict) else st.session_state.get("rna_fold_result")
    pairing = rna_3d.pairing_sync(representation, fold if isinstance(fold, dict) else None)
    st.caption(str(pairing.get("note") or ""))
    st.caption(
        f"Geometry kind={pairing.get('geometry_kind')}; pairing status="
        f"{pairing.get('pairing_status')} (MFE is not experimental 3D)."
    )
    pair_choices: list[tuple[int, int]] = []
    for item in list(pairing.get("illustrative_helix_pairs") or []):
        if isinstance(item, dict):
            left = item.get("plus_index_0based")
            if left is not None:
                pair_choices.append((int(left), int(left)))
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            pair_choices.append((int(item[0]), int(item[1])))
    if str(pairing.get("pairing_status") or "") == "PREDICTED":
        for item in list(pairing.get("mfe_pairs") or []):
            if isinstance(item, dict):
                left = item.get("i")
                right = item.get("j")
                if left is not None and right is not None:
                    pair_choices.append((int(left), int(right)))
    unique_pairs: list[tuple[int, int]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for pair in pair_choices:
        key = (min(pair), max(pair))
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        unique_pairs.append(pair)
    if unique_pairs:
        labels = [
            f"{left}-{right} ({pairing.get('pairing_status')})"
            for left, right in unique_pairs
        ]
        pair_index = int(
            st.selectbox(
                "Secondary / helix pair (sequence indices)",
                list(range(len(labels))),
                format_func=lambda i: labels[i],
                key="rna_pair_pick",
            )
        )
        if st.checkbox("Highlight selected pair in 3D", value=False, key="rna_pair_highlight"):
            left, right = unique_pairs[pair_index]
            highlights = list(highlights) + [left, right]
            st.caption(
                "Pair highlight uses sequence indices. PREDICTED MFE pairs do not "
                "relabel ILLUSTRATIVE XYZ as experimental coordinates."
            )
    envelope = representation.get("envelope")
    if isinstance(envelope, dict):
        _render_validated_scene_3d(
            rna_seq,
            envelope,
            key_prefix="rna_illust",
            selected=selected,
            highlight_indices=highlights,
            position_key="rna_3d_pos",
            lod=lod,
        )
    elif str(representation.get("kind") or "") == rna_3d.KIND_UNAVAILABLE:
        st.info(str(representation.get("reason") or "No RNA 3D representation."))


def _render_crispr_complex_section(dna_seq: str, guides: list) -> None:
    """Fetch RCSB Cas-guide-DNA complexes. No synthetic Cas9 assembly."""
    complexes = crispr_structure_catalog.list_experimental_complexes()
    labels = [f"{item['pdb_id']} | {item['title']}" for item in complexes]
    pick = st.selectbox(
        "Deposited Cas-guide-DNA complex (RCSB PDB)",
        list(range(len(labels))),
        format_func=lambda i: labels[i],
        key="crispr_complex_pick",
    )
    chosen = complexes[int(pick)]
    st.caption(str(chosen.get("notes") or ""))
    guide_seq = ""
    pam = "NGG"
    if guides:
        guide_seq = str(guides[0].get("guide_sequence") or "")
        pam = str(guides[0].get("pam_sequence") or "NGG")
    if st.button("Fetch deposited Cas complex from RCSB", key="crispr_complex_fetch"):
        try:
            with st.spinner("Retrieving deposited mmCIF from RCSB Files API..."):
                structure_text = ""
                fixture_name = f"{str(chosen['pdb_id']).upper()}.cif"
                reader = getattr(protein_structure, "read_bundled_mmcif", None)
                if callable(reader):
                    try:
                        structure_text = reader(fixture_name)
                    except protein_structure.StructureError:
                        structure_text = ""
                parsed = protein_structure.load_deposited_macromolecule(
                    source="RCSB PDB",
                    structure_id=str(chosen["pdb_id"]),
                    structure_text=structure_text,
                    cache=st.session_state.setdefault("crispr_complex_cache", {}),
                )
                annotated = complex_structure.annotate_parsed_complex(
                    parsed,
                    structure_id=str(chosen["pdb_id"]),
                    source="RCSB PDB",
                    metadata=parsed.get("metadata") if isinstance(parsed.get("metadata"), dict) else parsed,
                    guide_sequence=guide_seq,
                    pam=pam,
                )
        except (protein_structure.StructureError, complex_structure.ComplexError, ValueError) as exc:
            _render_structure_error(exc)
        else:
            annotated["bound_guide_hash"] = (
                provenance.sequence_digest(guide_seq) if guide_seq else ""
            )
            annotated["bound_target_hash"] = (
                provenance.sequence_digest(str(dna_seq or "")) if dna_seq else ""
            )
            st.session_state["crispr_complex"] = annotated
    annotated = st.session_state.get("crispr_complex")
    if not isinstance(annotated, dict):
        return
    live_guide_hash = provenance.sequence_digest(guide_seq) if guide_seq else ""
    live_target_hash = provenance.sequence_digest(str(dna_seq or "")) if dna_seq else ""
    stored_guide = str(annotated.get("bound_guide_hash") or "")
    stored_target = str(annotated.get("bound_target_hash") or "")
    if (stored_guide and live_guide_hash and stored_guide != live_guide_hash) or (
        stored_target and live_target_hash and stored_target != live_target_hash
    ):
        st.markdown(status_badge("STALE"), unsafe_allow_html=True)
        st.warning(
            "stale/incompatible structure: the deposited complex was fetched for a "
            "different guide or target hash. Fetch again after Design Guides."
        )
        return
    st.markdown(status_badge(str(annotated.get("status") or "RETRIEVED")), unsafe_allow_html=True)
    st.caption(str(annotated.get("disclaimer") or ""))
    chain_rows = [
        {
            "Chain": item.get("chain_id"),
            "Role": item.get("role"),
            "Molecule": item.get("molecule"),
            "Length": item.get("length"),
            "Description": item.get("description") or "",
        }
        for item in list(annotated.get("chains") or [])
    ]
    st.dataframe(_export_frame(chain_rows), width="stretch", hide_index=True)
    legend = helix_workspace.chain_role_legend_html(list(annotated.get("chains") or []))
    if legend:
        st.markdown(legend, unsafe_allow_html=True)
        st.caption(
            "Chain roles are those already classified from deposited polymer "
            "descriptions (for example 4UN3). Unknown is not inferred. "
            "This is not a synthetic Cas-guide-DNA assembly."
        )
    guide_map = annotated.get("guide_mapping") or {}
    pam_map = annotated.get("pam_mapping") or {}
    st.caption(
        f"Guide mapping: {guide_map.get('status')} "
        f"{guide_map.get('reason') or guide_map.get('method') or ''}. "
        f"PAM mapping: {pam_map.get('status')} "
        f"{pam_map.get('reason') or pam_map.get('method') or ''}."
    )
    try:
        envelope = complex_structure.deposited_to_structure_3d_input(
            annotated,
            query_sequence=str(dna_seq or ""),
            query_molecule="DNA",
        )
    except complex_structure.ComplexError as exc:
        st.error(str(exc))
        return
    highlights: list[int] = []
    highlights.extend(
        complex_structure.highlight_indices_for_chain_span(
            envelope.get("mappings") or [],
            str(guide_map.get("chain_id") or ""),
            guide_map.get("start_0based"),
            guide_map.get("end_0based"),
        )
    )
    highlights.extend(
        complex_structure.highlight_indices_for_chain_span(
            envelope.get("mappings") or [],
            str(pam_map.get("chain_id") or ""),
            pam_map.get("start_0based"),
            pam_map.get("end_0based"),
        )
    )
    query = str(dna_seq or "")
    _consume_3d_pending_position("crispr_complex_pos", len(query))
    selected_complex = 0
    if query:
        selected_complex = int(
            st.number_input(
                "Inspect CRISPR-complex query position (0-based)",
                min_value=0,
                max_value=max(0, len(query) - 1),
                value=0,
                step=1,
                key="crispr_complex_pos",
            )
        )
    _render_validated_scene_3d(
        query,
        envelope,
        key_prefix="crispr_complex",
        selected=selected_complex,
        highlight_indices=highlights,
        position_key="crispr_complex_pos" if query else "",
    )


def render_dna_analysis() -> None:
    """Renderiza a aba de analise de DNA.

    Coleta a sequencia (texto ou FASTA), valida, e ao clicar em Analyze exibe
    metricas, composicao, GC skew, sitios de restricao, ORFs, ilhas CpG e o
    visualizador da sequencia.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    in_col, side_col = st.columns([2, 1])
    with in_col:
        text_value = st.text_area("Paste a DNA sequence", height=200, key="dna_text")
    with side_col:
        uploaded = st.file_uploader(
            "Or upload a FASTA file", type=["fasta", "fa"], key="dna_file"
        )
        st.caption(
            "A rejected upload is not analyzed. Use a simple file name such as sample.fasta."
        )
        analyze = st.button("Analyze", key="dna_analyze")

    if analyze:
        try:
            resolved = _resolve_sequence(text_value, uploaded)
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state["dna_input"] = resolved
        st.session_state.pop("dna_3d_pos_from_3d", None)
        st.session_state.pop("dna_3d_representation", None)
        st.session_state.pop("dna_struct_result", None)
        st.session_state.pop("dna_struct_hits", None)
        if not resolved.strip():
            st.error("Provide a DNA sequence or a FASTA file to analyze.")
            return

    seq = st.session_state.get("dna_input", "")
    if not seq:
        return
    _warn_if_analyzed_input_changed(text_value, seq, action="Analyze")

    info = _validate_for(seq, "DNA")
    st.markdown(_type_badge(info), unsafe_allow_html=True)
    st.markdown(status_badge("COMPUTED"), unsafe_allow_html=True)

    if not info["is_valid"]:
        st.error(info["rejection_reason"] or "Provide a valid DNA sequence (A, T, C, G only).")
        return

    clean = info["sequence"]
    _store_workspace("DNA", clean, source="user sequence", identifier="")
    length = info["length"]
    gc = _gc_content(clean)
    st.markdown(
        helix_workspace.result_header_html(
            "DNA analysis",
            status="COMPUTED",
            lines=[f"{length:,} bp", f"GC {_format_percent(gc)}%"],
        ),
        unsafe_allow_html=True,
    )
    tm_report = _melting_temperature_report(clean)
    if tm_report["status"] == "COMPUTED" and tm_report["value_c"] is not None:
        tm_value = f"{tm_report['value_c']:.2f}"
        tm_method_label = str(tm_report["method"])
    else:
        tm_value = "Unavailable"
        tm_method_label = "Unavailable for this sequence/method"
        
    mw = _molecular_weight(clean, "DNA")
    if isinstance(mw, float) and math.isnan(mw):
        mw_str = "Unavailable"
        mw_unit = ""
    elif mw >= 1e6:
        mw_str = f"{mw / 1e6:.2f}"
        mw_unit = "MDa"
    elif mw >= 1e3:
        mw_str = f"{mw / 1e3:.2f}"
        mw_unit = "kDa"
    else:
        mw_str = f"{mw:.2f}"
        mw_unit = "Da"

    at = _at_content(clean)
    entropy = _shannon_entropy(clean)
    skew_at = _at_skew(clean)
    skew_gc = dna_analysis.gc_skew_global(clean)
    skew_at_label = "Unavailable" if math.isnan(skew_at) else f"{skew_at:.4f}"
    skew_gc_label = "Unavailable" if math.isnan(skew_gc) else f"{skew_gc:.4f}"
    helix_workspace.render_explanation(
        "dna",
        {
            "status": "COMPUTED",
            "length": length,
            "gc_percent": gc,
            "at_percent": at,
            "method": "nucleotide composition",
            "source": "HelixScope DNA analysis",
        },
    )
    helix_workspace.render_explanation(
        "dna_tm",
        {
            "status": tm_report.get("status"),
            "value_c": tm_report.get("value_c"),
            "method": tm_report.get("method"),
            "source": "HelixScope DNA analysis",
        },
    )

    entropy_label = (
        "Unavailable"
        if isinstance(entropy, float) and math.isnan(entropy)
        else f"{entropy:.4f}"
    )
    render_metric_grid(
        [
            ("Length", str(length), "bp"),
            (
                "GC Content",
                _format_percent(gc),
                "%" if _format_percent(gc) != "Unavailable" else "",
            ),
            (
                "AT Content",
                _format_percent(at),
                "%" if _format_percent(at) != "Unavailable" else "",
            ),
            ("Tm", tm_value, "C" if "Unavailable" not in tm_value else ""),
            ("MW", mw_str, mw_unit),
            (
                "Shannon entropy",
                entropy_label,
                "bits" if entropy_label != "Unavailable" else "",
            ),
            ("AT skew", skew_at_label, ""),
            ("GC skew", skew_gc_label, ""),
        ]
    )
    cumulative = _cumulative_skew_summary(clean)
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Tm method: "
        f"{html_escape(tm_method_label)}. Wallace for fewer than 14 canonical bp "
        "and salt-adjusted (0.05 M Na+) from 14 to 200 bp. That value is not "
        "nearest-neighbor. Sequences longer than 200 bp or without canonical "
        "bases are N/A, not 0. GC skew is (G-C)/(G+C); AT "
        "skew is (A-T)/(A+T); a zero denominator is N/A, not 0. "
        f"INPUT_LENGTH {cumulative['input_length']:,}. "
        f"ANALYZED_LENGTH {cumulative['analyzed_length']:,}. "
        f"VISUALIZED_LENGTH {cumulative['visualized_length']:,} "
        f"({cumulative['visualization']}, step {cumulative['plot_step']}). "
        f"Lobry cumulative GC {cumulative['final_gc']:+d}, "
        f"AT {cumulative['final_at']:+d}. "
        f"{html_escape(str(cumulative['method']))}</span>",
        unsafe_allow_html=True,
    )
    tm_na_m = st.number_input(
        "SantaLucia Na+ (M)",
        min_value=0.001,
        max_value=1.5,
        value=0.05,
        step=0.01,
        format="%.3f",
        key="dna_santalucia_na",
    )
    tm_oligo_nm = st.number_input(
        "SantaLucia total oligo (nM)",
        min_value=1.0,
        max_value=1.0e6,
        value=250.0,
        step=10.0,
        key="dna_santalucia_oligo",
    )
    nn_tm = _santalucia_tm(clean, float(tm_na_m), float(tm_oligo_nm), 0.0)
    nn_value = (
        f"{nn_tm['value_c']:.2f}"
        if nn_tm.get("status") == "COMPUTED" and nn_tm.get("value_c") is not None
        else "N/A"
    )
    render_metric_grid(
        [
            (
                "Tm SantaLucia NN",
                nn_value,
                "C" if nn_value != "N/A" else "",
            )
        ],
        columns=4,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(str(nn_tm.get('method') or 'SantaLucia nearest-neighbor'))}. "
        f"{html_escape(str(nn_tm.get('citation') or ''))} "
        f"Na+={tm_na_m} M; oligo={tm_oligo_nm} nM; Mg2+=0 (Mg correction is not mixed in). "
        f"{html_escape(str(nn_tm.get('reason') or ''))}"
        "</span>",
        unsafe_allow_html=True,
    )

    if _section_toggle("Nucleotide Composition", key="dna_nuc_comp"):
        with st.container(border=True):
            composition = _composition(clean)
            st.plotly_chart(
                charts.nucleotide_bar_chart(composition), width="stretch"
            )
            comp_rows = [
                {"Symbol": symbol, "Count": data["count"], "Frequency_pct": data["frequency"]}
                for symbol, data in composition.items()
                if data["count"] > 0
            ]
            st.dataframe(pd.DataFrame(comp_rows), width="stretch")
            if not dna_analysis.composition_is_consistent(clean):
                st.warning(
                    "Composition counts do not sum to sequence length. "
                    "The result is marked inconsistent and was not adjusted."
                )

    if _section_toggle("Dinucleotide Frequencies and CpG Odds Ratio", key="dna_dinuc"):
        with st.container(border=True):
            dinucleotides = _dinucleotides(clean)
            st.plotly_chart(
                charts.dinucleotide_heatmap(dinucleotides), width="stretch"
            )
            dinuc_rows = [
                {
                    "Dinucleotide": pair,
                    "Count": data["count"],
                    "Frequency_pct": data["frequency"],
                    "Observed_Expected": data["observed_expected"],
                }
                for pair, data in dinucleotides.items()
            ]
            dinuc_df = pd.DataFrame(dinuc_rows).sort_values(
                by="Observed_Expected", ascending=False
            )
            st.dataframe(dinuc_df, width="stretch", hide_index=True)
            st.download_button(
                "Download dinucleotides (CSV)",
                data=_export_frame(dinuc_rows).to_csv(index=False).encode("utf-8"),
                file_name="dna_dinucleotides.csv",
                mime="text/csv",
                key="dna_dinuc_download",
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Observed/Expected near "
                "1.0 means adjacent bases behave independently. Values far from 1.0 "
                "indicate selection or mutational pressure on that dinucleotide. "
                "When the expected frequency is zero the ratio is N/A, not 0."
                "</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("GC Skew", key="dna_gc_skew"):
        with st.container(border=True):
            plan = scale_profile.analysis_plan(length, molecule="DNA")
            window = int(plan["suggested_window"]) if length >= 100 else length
            st.caption(
                f"Non-overlapping windows of {window} bp. Suggested for this length "
                f"by the scale profile (analysis target {scale_profile.TARGET_DNA_NT:,} nt). "
                "This changes window size, not the (G-C)/(G+C) formula."
            )
            skew = _gc_skew(clean, window)
            st.plotly_chart(charts.gc_skew_plot(skew, window), width="stretch")

    if _section_toggle("GC Content in Sliding Window", key="dna_gc_window"):
        with st.container(border=True):
            plan = scale_profile.analysis_plan(length, molecule="DNA")
            max_window = max(10, min(2000, length))
            default_window = min(max_window, int(plan["suggested_window"]))
            default_step = min(max_window, max(1, int(plan["suggested_step"])))
            win_col, step_col = st.columns(2)
            with win_col:
                gc_window = st.number_input(
                    "Window size (bp)",
                    min_value=10,
                    max_value=max_window,
                    value=default_window,
                    step=10,
                    key="dna_gc_window_size",
                )
            with step_col:
                gc_step = st.number_input(
                    "Step (bp)",
                    min_value=1,
                    max_value=max_window,
                    value=default_step,
                    step=5,
                    key="dna_gc_window_step",
                )
            st.caption(
                f"Display/memory cap is {scale_profile.MAX_WINDOW_PROFILE_ROWS:,} windows. "
                "If the grid is too dense, increase the step. The GC formula is unchanged."
            )
            try:
                profile = _gc_sliding_window(clean, int(gc_window), int(gc_step))
                st.plotly_chart(
                    charts.gc_sliding_window_plot(profile), width="stretch"
                )
                windows = _windowed_profiles(clean, int(gc_window), int(gc_step))
                if windows:
                    window_df = _export_frame(windows)
                    st.dataframe(window_df, width="stretch", hide_index=True)
                    st.download_button(
                        "Download windowed profiles (CSV)",
                        data=window_df.to_csv(index=False).encode("utf-8"),
                        file_name="dna_windowed_profiles.csv",
                        mime="text/csv",
                        key="dna_windowed_download",
                    )
                    st.markdown(
                        "<span style='color:var(--hs-text-secondary);font-size:12px;'>Coordinates "
                        "are 0-based, end exclusive. Alphabet: ACGT. window_size "
                        f"and step_size are {int(gc_window)} and {int(gc_step)}. "
                        "Undefined metrics are N/A, not 0.</span>",
                        unsafe_allow_html=True,
                    )
            except ValueError as exc:
                st.info(str(exc))

    if _section_toggle("Entropy, k-mers and AT skew", key="dna_stats"):
        with st.container(border=True):
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Shannon entropy uses "
                "canonical A,C,G,T only (max 2 bits). k-mers ignore windows that "
                "contain non-ACGT symbols. AT skew is (A-T)/(A+T).</span>",
                unsafe_allow_html=True,
            )
            max_window = max(10, min(2000, length))
            plan = scale_profile.analysis_plan(length, molecule="DNA")
            default_window = min(max_window, int(plan["suggested_window"]))
            k_col, win_col, step_col = st.columns(3)
            with k_col:
                kmer_k = st.number_input(
                    "k-mer length",
                    min_value=1,
                    max_value=5,
                    value=3,
                    key="dna_kmer_k",
                )
            with win_col:
                ent_window = st.number_input(
                    "Entropy window (bp)",
                    min_value=10,
                    max_value=max_window,
                    value=default_window,
                    key="dna_entropy_window",
                )
            with step_col:
                ent_step = st.number_input(
                    "Entropy step (bp)",
                    min_value=1,
                    max_value=max_window,
                    value=min(max_window, max(1, int(plan["suggested_step"]))),
                    key="dna_entropy_step",
                )
            kmers = _kmer_counts_dna(clean, int(kmer_k))
            summary = _kmer_summary_dna(clean, int(kmer_k))
            if kmers:
                st.plotly_chart(charts.kmer_bar_chart(kmers), width="stretch")
                diversity = provenance.csv_cell(summary.get("diversity"))
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);font-size:12px;'>k={html_escape(str(summary['k']))}; "
                    f"valid windows={summary['valid_windows']:,}; unique="
                    f"{summary['unique_kmers']:,} / {summary['possible_kmers']:,} "
                    f"possible; diversity={html_escape(str(diversity))}; "
                    f"status={html_escape(str(summary['status']))}. k is capped at "
                    f"{dna_analysis.MAX_KMER_K} to avoid memory explosion.</span>",
                    unsafe_allow_html=True,
                )
                if summary.get("top"):
                    st.dataframe(
                        _export_frame(summary["top"]),
                        width="stretch",
                        hide_index=True,
                    )
            else:
                st.info("Insufficient data for k-mer counts.")
            try:
                ent_profile = _entropy_profile(clean, int(ent_window), int(ent_step))
            except ValueError as exc:
                st.info(str(exc))
                ent_profile = []
            if ent_profile:
                st.plotly_chart(charts.entropy_profile_plot(ent_profile), width="stretch")
            else:
                st.info("Insufficient data for an entropy profile at this window size.")
            plan_skew = scale_profile.analysis_plan(length, molecule="DNA")
            at_window = int(plan_skew["suggested_window"]) if length >= 100 else length
            try:
                at_profile = _at_skew_windows(clean, at_window)
                st.plotly_chart(
                    charts.gc_skew_plot(at_profile, at_window, y_title="AT Skew"),
                    width="stretch",
                )
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);font-size:12px;'>The plot above is "
                    "AT skew versus position, using the same renderer as GC skew. "
                    "Y values are (A-T)/(A+T), not GC skew.</span>",
                    unsafe_allow_html=True,
                )
            except ValueError as exc:
                st.info(str(exc))

    if _section_toggle("Sequence map (predicted ORFs)", key="dna_map"):
        with st.container(border=True):
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Predicted ORFs, not "
                "confirmed genes. Coordinates are mapped onto the input strand. "
                "Method: ATG start, standard stops, non-overlapping selection, "
                f"minimum {dna_analysis.MIN_CDS_LENGTH_NT} nt.</span>",
                unsafe_allow_html=True,
            )
            map_orfs = _orfs_non_overlapping(
                _orfs(clean, dna_analysis.MIN_CDS_LENGTH_NT),
                dna_analysis.MIN_CDS_LENGTH_NT,
            )
            features = []
            for orf in map_orfs:
                try:
                    start, end = dna_analysis.orf_coordinates_on_input(orf, length)
                except ValueError:
                    continue
                features.append(
                    {
                        "start": start,
                        "end": end,
                        "label": f"Predicted ORF frame {orf['frame']}",
                    }
                )
            if features:
                st.plotly_chart(
                    charts.sequence_feature_map(
                        length, [{"name": "Predicted ORFs", "features": features}]
                    ),
                    width="stretch",
                )
            else:
                st.info(
                    "Insufficient data: no predicted ORFs at the default minimum "
                    "length."
                )

    if _section_toggle("Restriction Sites", key="dna_restriction"):
        with st.container(border=True):
            sites = _restriction_sites(clean)
            site_rows = [
                {
                    "Enzyme": enzyme,
                    "Pattern": data["pattern"],
                    "Sites": ", ".join(str(p) for p in data["sites"]) or "-",
                    "Count": data["count"],
                }
                for enzyme, data in sites.items()
            ]
            st.dataframe(pd.DataFrame(site_rows), width="stretch")
            hit_rows = [
                {
                    "enzyme": hit["enzyme"],
                    "site": hit["site"],
                    "position": hit["position"],
                    "strand": hit["strand"],
                }
                for data in sites.values()
                for hit in data.get("hits") or []
            ]
            if hit_rows:
                st.dataframe(
                    pd.DataFrame(hit_rows), width="stretch", hide_index=True
                )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Positions are 1-based "
                "starts of the recognition site on the input strand, not the "
                "cleavage cut reported by Bio.Restriction.search. Enzymes are the "
                "six names in RESTRICTION_ENZYME_NAMES (EcoRI, BamHI, HindIII, "
                "NcoI, NdeI, XhoI). Palindromes are reported on strand + of the "
                "input sequence.</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("Open Reading Frames", key="dna_orf"):
        with st.container(border=True):
            ctrl_col, opt_col = st.columns([1, 2])
            with ctrl_col:
                min_orf = st.number_input(
                    "Minimum ORF length (nt)",
                    min_value=30,
                    max_value=3000,
                    value=dna_analysis.MIN_CDS_LENGTH_NT,
                    step=30,
                    key="dna_orf_min_length",
                )
            with opt_col:
                drop_nested = st.checkbox(
                    "Discard nested and overlapping ORFs, keeping one per locus",
                    value=True,
                    key="dna_orf_drop_nested",
                )

            raw_orfs = _orfs(clean, int(min_orf))
            orfs = (
                _orfs_non_overlapping(raw_orfs, int(min_orf))
                if drop_nested
                else raw_orfs
            )

            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Method: ORF prediction "
                f"with ATG start, standard stop codons and minimum length of "
                f"{int(min_orf)} nt. Status of every row is PREDICTED, not gene. "
                "Coordinates start/end and input_start/input_end are 0-based "
                "half-open [start, end). "
                f"Raw candidates: {len(raw_orfs):,}. "
                f"Reported: {len(orfs):,}.</span>",
                unsafe_allow_html=True,
            )
            n_codons = sum(int(orf.get("length_bp") or 0) for orf in orfs) // 3
            helix_workspace.render_explanation(
                "dna_orf",
                {
                    "status": "PREDICTED",
                    "n_orfs": len(orfs),
                    "n_codons": n_codons,
                    "method": "ATG start, standard stop, non-overlapping filter",
                    "source": "HelixScope DNA analysis",
                    "method_note": (
                        f"Minimum length {int(min_orf)} nt. "
                        f"Raw candidates {len(raw_orfs)}; reported {len(orfs)}."
                    ),
                },
            )

            if orfs:
                orf_df = pd.DataFrame(orfs)
                st.dataframe(orf_df, width="stretch")
                st.download_button(
                    "Download ORFs (CSV)",
                    data=orf_df.to_csv(index=False).encode("utf-8"),
                    file_name="orfs.csv",
                    mime="text/csv",
                    key="dna_orf_download",
                )
                st.plotly_chart(
                    charts.orf_length_histogram(
                        [int(orf["length_bp"]) for orf in orfs]
                    ),
                    width="stretch",
                )
            else:
                st.info(
                    f"No ORFs found with the minimum length of {int(min_orf)} nt. "
                    "Lower the threshold or confirm that the input contains a "
                    "coding region."
                )

    if _section_toggle("CpG Islands", key="dna_cpg"):
        with st.container(border=True):
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Prediction using "
                "Gardiner-Garden and Frommer (1987) thresholds (GC at least 50% "
                "and CpG O/E at least 0.6) with a 200 bp window and a 50 bp step. "
                "The original paper scanned with step 1.</span>",
                unsafe_allow_html=True,
            )
            islands = _cpg_islands(clean)
            cpg_positions = dna_analysis.cpg_dinucleotide_positions(clean)
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Direct CpG count on "
                f"the input strand: {len(cpg_positions):,} dinucleotides. "
                "Listed island coordinates are 0-based and must fall inside the "
                "sequence.</span>",
                unsafe_allow_html=True,
            )
            if islands:
                cpg_df = pd.DataFrame(islands).rename(
                    columns={
                        "start": "Start",
                        "end": "End",
                        "gc_percent": "GC%",
                        "cpg_oe": "CpG O/E",
                        "cpg_count": "CpG count",
                        "length_bp": "Length_bp",
                    }
                )
                st.dataframe(cpg_df, width="stretch")
            else:
                st.info("No CpG islands detected.")

    if _section_toggle("Sequence Viewer", key="dna_seq_view"):
        with st.container(border=True):
            meta = st.session_state.get("workspace_dna_meta") or {}
            digest = html_escape(str(meta.get("input_hash") or ""))
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Workspace DNA hash "
                f"{digest[:16]}... length {length}. Viewer and metrics use this "
                "same sequence.</span>",
                unsafe_allow_html=True,
            )
            _sequence_viewer_panel(clean, "DNA", key_prefix="dna_view")

    if _section_toggle("DNA 3D", key="dna_3d"):
        with st.container(border=True):
            _render_dna_3d_section(clean)

    dinucleotides = _dinucleotides(clean)
    cpg_odds = dinucleotides.get("CG", {}).get("observed_expected", float("nan"))
    summary_rows = [
        {"Metric": "Length", "Value": length, "Unit": "bp"},
        {"Metric": "GC content", "Value": gc, "Unit": "%"},
        {"Metric": "AT content", "Value": at, "Unit": "%"},
        {"Metric": "Melting temperature", "Value": tm_value, "Unit": "C"},
        {"Metric": "Tm method", "Value": tm_method_label, "Unit": ""},
        {"Metric": "Molecular weight", "Value": f"{mw_str} {mw_unit}".strip(), "Unit": ""},
        {"Metric": "CpG observed/expected", "Value": cpg_odds, "Unit": "ratio"},
        {"Metric": "Shannon entropy", "Value": entropy, "Unit": "bits"},
        {"Metric": "AT skew", "Value": skew_at, "Unit": ""},
        {"Metric": "GC skew", "Value": skew_gc, "Unit": ""},
        {
            "Metric": "DNA 3D coordinates",
            "Value": (
                "Valid DNA has an ILLUSTRATIVE helix path by default; experimental "
                "3D only from a deposited RCSB structure whose sequence hash matches. "
                "Sequence is not coordinates. PREDICTED DNA 3D is UNAVAILABLE without "
                "a deposited predicted envelope."
            ),
            "Unit": "",
        },
    ]
    _render_summary_table(summary_rows, "dna_summary_download", "dna_summary.csv")
    report = provenance.analysis_envelope(
        module="DNA",
        payload={
            "length_bp": length,
            "gc_percent": gc,
            "at_percent": at,
            "shannon_entropy_bits": entropy,
            "at_skew": skew_at,
            "gc_skew": skew_gc,
            "tm_value_c": tm_report.get("value_c"),
            "tm_method": tm_method_label,
            "tm_status": tm_report.get("status"),
            "tm_santalucia_value_c": nn_tm.get("value_c"),
            "tm_santalucia_status": nn_tm.get("status"),
            "tm_santalucia_parameter_set": nn_tm.get("parameter_set"),
            "molecular_weight": None if mw_str == "N/A" else f"{mw_str} {mw_unit}".strip(),
            "cpg_observed_expected": cpg_odds,
            "composition_consistent": dna_analysis.composition_is_consistent(clean),
            "methods": {
                "gc": "canonical A,C,G,T denominator; NaN if no canonical bases",
                "entropy": "Shannon on ACGT",
                "orfs": "ATG start, standard stops, predicted not confirmed genes",
                "tm": "Wallace <14 bp; salt-adjusted 14-200 bp; unavailable >200 bp",
                "tm_santalucia": "SantaLucia & Hicks 2004 NN + SantaLucia 1998 Na entropy; not Wallace",
            },
        },
        status="COMPUTED",
        algorithm="HelixScope DNA composition, entropy, ORF prediction",
        parameters={"length_bp": length},
        sequence=clean,
        source="user sequence",
    )
    st.download_button(
        "Download analysis report (JSON)",
        data=json.dumps(report, indent=2).encode("utf-8"),
        file_name="helixscope_dna_report.json",
        mime="application/json",
        key="dna_report_json",
    )
    _render_alerts(_nucleotide_alerts(gc, dinucleotides))


_FOLDING_ERROR_LABELS: dict[str, str] = {
    "TIMEOUT": "Timeout",
    "TOOL_NOT_INSTALLED": "Tool unavailable",
    "TOOL_FAILED": "Tool failure",
    "JOB_FAILED": "Tool failure",
    "INVALID_INPUT": "Invalid input",
    "PARSING_ERROR": "Parsing error",
    "RESOURCE_LIMIT": "Resource limit",
    "NOT_FOUND": "Not found",
    "NETWORK_ERROR": "Network error",
    "NO_STRUCTURE": "No validated structure found",
    "SERVICE_UNAVAILABLE": "Service unavailable",
    "RATE_LIMITED": "Rate limited",
}


def _render_folding_error(exc: Exception) -> None:
    """Exibe falha de folding classificada, nunca como estrutura desenhada."""
    if isinstance(exc, rna_folding.FoldingError):
        label = _FOLDING_ERROR_LABELS.get(exc.category, exc.category)
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(f"RNA structure prediction failed. {label}. {exc}")
        return
    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
    st.error(f"RNA structure prediction failed. {exc}")


def _render_rna_structure_section(rna_seq: str) -> Optional[dict]:
    """Folding sob demanda. Nao dispara ao abrir a aba nem ao analisar composicao.

    Args:
        rna_seq: RNA ja validada nesta aba.

    Returns:
        Envelope PREDICTED se houver resultado na sessao para este hash; senao None.
    """
    availability = rna_folding.tool_availability()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(availability['reason'])} "
        "A predicted structure is not experimental, not functional annotation, "
        "and not a 3D model."
        "</span>",
        unsafe_allow_html=True,
    )
    if not availability["available"]:
        st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
        st.info("RNA folding unavailable in this environment.")
        for tool in availability.get("recommended_tools") or []:
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);'>- {html_escape(tool)}</span>",
                unsafe_allow_html=True,
            )
        return None

    backends = rna_folding.available_backend_ids()
    backend = st.selectbox(
        "ViennaRNA backend (only tools that can actually run)",
        backends,
        key="rna_fold_backend",
    )
    fold_region = st.checkbox(
        "Fold a selected subsequence (not native local structure)",
        key="rna_fold_region",
    )
    region_start = 0
    region_end = len(rna_seq)
    if fold_region:
        region_start = int(
            st.number_input(
                "Region start (0-based)",
                min_value=0,
                max_value=max(0, len(rna_seq) - 1),
                value=0,
                step=1,
                key="rna_fold_start",
            )
        )
        region_end = int(
            st.number_input(
                "Region end (exclusive, 0-based)",
                min_value=region_start + 1,
                max_value=len(rna_seq),
                value=len(rna_seq),
                step=1,
                key="rna_fold_end",
            )
        )
        st.caption(
            "Folding performed on selected subsequence. This is not the native "
            "local structure of that region inside the full RNA."
        )

    meta = st.session_state.get("workspace_rna_meta") or {}
    seq_hash = provenance.sequence_digest(rna_seq)
    source = "user sequence"
    accession = ""
    version = ""
    organism = ""
    identifier = ""
    retrieved_at = ""
    database = ""
    if str(meta.get("input_hash") or "") == seq_hash:
        source = str(meta.get("source") or source)
        accession = str(meta.get("accession") or "")
        version = str(meta.get("version") or "")
        organism = str(meta.get("organism") or "")
        identifier = str(meta.get("identifier") or "")
        retrieved_at = str(meta.get("timestamp") or "")

    folded = rna_seq[region_start:region_end] if fold_region else rna_seq
    tool_version = ""
    for item in availability.get("backends") or []:
        if str(item.get("id")) == str(backend):
            tool_version = str(item.get("version") or "")
            break
    parameters = rna_folding.default_fold_parameters(
        "selected subsequence" if fold_region else "full sequence"
    )
    key = rna_folding.cache_key(
        sequence_hash=provenance.sequence_digest(folded),
        backend=str(backend),
        tool_version=tool_version,
        parameters=parameters,
    )
    admission = rna_folding.folding_admission(len(rna_seq))
    st.caption(
        f"INPUT_LENGTH {len(rna_seq):,}. "
        f"Global fold: {admission['global_fold']}. "
        f"Local long-RNA: {admission['local_long_rna']}. "
        f"{admission['reason']}"
    )
    if st.button("Local pair probabilities (ViennaRNA pfl_fold)", key="rna_local_run"):
        try:
            with st.spinner("Computing local pair probabilities..."):
                local = rna_folding.local_long_rna_analysis(rna_seq)
        except (rna_folding.FoldingError, ValueError) as exc:
            _render_folding_error(exc)
        else:
            st.session_state["rna_local_result"] = local
    local_result = st.session_state.get("rna_local_result")
    if isinstance(local_result, dict) and str(local_result.get("sequence_hash") or "") == seq_hash:
        st.markdown(status_badge("PREDICTED"), unsafe_allow_html=True)
        st.caption(
            f"Mode {local_result.get('mode')}. "
            f"INPUT_LENGTH {local_result.get('input_length'):,}. "
            f"ANALYZED_LENGTH {local_result.get('analyzed_length'):,}. "
            f"VISUALIZED_LENGTH {local_result.get('visualized_length'):,} "
            f"(pair preview, not the full catalogue). "
            f"Window {local_result.get('window_size')}. "
            f"Max pair span {local_result.get('max_pair_span')}. "
            f"Engine {local_result.get('engine')} {local_result.get('engine_version')}. "
            f"{local_result.get('disclaimer')}"
        )
        render_metric_grid(
            [
                ("Local pairs at cutoff", str(local_result.get("pair_count")), ""),
                (
                    "Max pair probability",
                    "Unavailable"
                    if local_result.get("max_pair_probability") is None
                    else f"{float(local_result['max_pair_probability']):.4f}",
                    "",
                ),
                (
                    "Mean unpaired probability",
                    "Unavailable"
                    if local_result.get("unpaired_probability_mean") is None
                    else f"{float(local_result['unpaired_probability_mean']):.4f}",
                    local_result.get("unpaired_probability_status") or "",
                ),
            ]
        )
    if st.button("Predict secondary structure (ViennaRNA MFE)", key="rna_fold_run"):
        cache = st.session_state.setdefault("rna_fold_cache", {})
        cached = cache.get(key)
        if isinstance(cached, dict):
            st.session_state["rna_fold_result"] = rna_folding.mark_cached_result(cached)
            st.session_state["rna_fold_cache_key"] = key
        else:
            try:
                with st.spinner("Computing RNA secondary structure..."):
                    result = rna_folding.fold_rna(
                        rna_seq,
                        source=source,
                        identifier=identifier,
                        accession=accession,
                        version=version,
                        organism=organism,
                        database=database,
                        retrieved_at=retrieved_at,
                        region_start=region_start if fold_region else None,
                        region_end=region_end if fold_region else None,
                        full_sequence=rna_seq if fold_region else "",
                        backend=str(backend),
                    )
            except (rna_folding.FoldingError, ValueError) as exc:
                _render_folding_error(exc)
            else:
                live_key = rna_folding.cache_key(
                    sequence_hash=str(result.get("sequence_hash") or ""),
                    backend=str(result.get("backend") or backend),
                    tool_version=str(result.get("tool_version") or ""),
                    parameters=dict(result.get("parameters") or parameters),
                )
                cache[live_key] = result
                st.session_state["rna_fold_result"] = result
                st.session_state["rna_fold_cache_key"] = live_key

    result = st.session_state.get("rna_fold_result")
    if not isinstance(result, dict):
        return None
    folded_hash = provenance.sequence_digest(folded)
    if str(result.get("sequence_hash")) not in {folded_hash, seq_hash} and str(
        result.get("full_sequence_hash") or ""
    ) != seq_hash:
        return None

    st.markdown(status_badge("PREDICTED"), unsafe_allow_html=True)
    if str(result.get("cache_status") or "") == "cached":
        st.markdown(status_badge("CACHED"), unsafe_allow_html=True)
    helix_workspace.render_explanation("rna_fold", result)
    version_text = result.get("tool_version") or "unavailable"
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"tool {html_escape(result.get('tool'))} | version {html_escape(version_text)} | "
        f"method {html_escape(result.get('method'))} | "
        f"sequence source {html_escape(result.get('sequence_source'))} | "
        f"structure source {html_escape(result.get('structure_source'))}"
        "</span>",
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(result.get('disclaimer'))}</span>",
        unsafe_allow_html=True,
    )
    mfe = result.get("mfe_kcal_mol")
    render_metric_grid(
        [
            (
                "Minimum Free Energy",
                f"{float(mfe):.2f}" if scientific_checks.mfe_kcal_mol_is_valid(mfe) else "Unavailable",
                str(result.get("mfe_unit") or "kcal/mol")
                if scientific_checks.mfe_kcal_mol_is_valid(mfe)
                else "",
            )
        ],
        columns=4,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Method: "
        f"{html_escape(result.get('method'))}. Temperature: "
        f"{html_escape((result.get('parameters') or {}).get('temperature_c'))} C "
        f"({html_escape((result.get('parameters') or {}).get('temperature_source'))})."
        "</span>",
        unsafe_allow_html=True,
    )
    st.code(
        str(result.get("sequence") or "") + "\n" + str(result.get("dot_bracket") or ""),
        language="text",
    )
    selected = int(
        st.number_input(
            "Inspect RNA position (0-based)",
            min_value=0,
            max_value=max(0, len(str(result.get("sequence") or "")) - 1),
            value=0,
            step=1,
            key="rna_fold_pos",
        )
    )
    pair = rna_folding.pair_at_position(result, selected)
    residue = str(result.get("sequence") or "")[selected]
    if pair is None:
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Position {selected}: "
            f"{html_escape(residue)} unpaired in the predicted MFE structure.</span>",
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Position {selected}: "
            f"{html_escape(pair.get('base'))} paired with "
            f"{html_escape(pair.get('paired_position_0based'))}:"
            f"{html_escape(pair.get('paired_base'))}.</span>",
            unsafe_allow_html=True,
        )
    try:
        paths = rna_folding.arc_paths(len(str(result.get("sequence") or "")), list(result.get("base_pairs") or []))
        st.plotly_chart(
            charts.rna_arc_diagram(str(result.get("sequence") or ""), paths, selected),
            width="stretch",
        )
        layout = rna_folding.circular_layout(str(result.get("sequence") or ""), list(result.get("base_pairs") or []))
        st.plotly_chart(charts.rna_circular_pairs(layout, selected), width="stretch")
    except (ValueError, rna_folding.FoldingError) as exc:
        st.info(str(exc))
    pair_rows = [
        {
            "Position_0based": item.get("position_0based"),
            "Base": item.get("base"),
            "Paired_0based": item.get("paired_position_0based"),
            "Paired_base": item.get("paired_base"),
        }
        for item in list(result.get("base_pairs") or [])
    ]
    if pair_rows:
        st.dataframe(_export_frame(pair_rows), width="stretch", hide_index=True)
    else:
        st.info("No base pairs in this predicted MFE structure.")
    element_rows = [
        {
            "Type": item.get("type"),
            "Start_0based": item.get("start"),
            "End_exclusive": item.get("end"),
            "Label": item.get("label"),
        }
        for item in list(result.get("elements") or [])
    ]
    if element_rows:
        st.caption(
            "Structural element labels are parsed from the predicted dot-bracket. "
            "They are not genes, CDS or regulatory annotations."
        )
        st.dataframe(_export_frame(element_rows), width="stretch", hide_index=True)
    export_cols = st.columns(3)
    with export_cols[0]:
        st.download_button(
            "Dot-bracket (DBN)",
            data=rna_folding.export_dot_bracket(result).encode("utf-8"),
            file_name="helixscope_rna_structure.dbn",
            mime="text/plain",
            key="rna_dbn",
        )
    with export_cols[1]:
        st.download_button(
            "Base pairs (CSV)",
            data=rna_folding.export_pairs_csv(result).encode("utf-8"),
            file_name="helixscope_rna_pairs.csv",
            mime="text/csv",
            key="rna_pairs_csv",
        )
    with export_cols[2]:
        st.download_button(
            "Structure JSON",
            data=json.dumps(rna_folding.export_result_bundle(result), indent=2).encode("utf-8"),
            file_name="helixscope_rna_structure.json",
            mime="application/json",
            key="rna_fold_json",
        )
    return result


def render_rna_analysis() -> None:
    """Renderiza a aba de analise de RNA.

    Usa somente a sequencia colada ou enviada nesta aba. Nao reutiliza o DNA da
    aba DNA e nao transcreve DNA. Aceita apenas RNA no alfabeto A, U, C e G.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    in_col, side_col = st.columns([2, 1])
    with in_col:
        text_value = st.text_area(
            "Paste an RNA sequence", height=200, key="rna_text"
        )
    with side_col:
        uploaded = st.file_uploader(
            "Or upload a FASTA file", type=["fasta", "fa"], key="rna_file"
        )
        st.caption(
            "A rejected upload is not analyzed. Use a simple file name such as sample.fasta."
        )
        analyze = st.button("Analyze", key="rna_analyze")

    if analyze:
        try:
            resolved = _resolve_sequence(text_value, uploaded)
            raw_text = _resolve_raw_text(text_value, uploaded)
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state["rna_input"] = resolved
        st.session_state.pop("rna_3d_pos_from_3d", None)
        st.session_state.pop("rna_3d_representation", None)
        st.session_state.pop("rna_struct_result", None)
        st.session_state.pop("rna_struct_hits", None)
        st.session_state["rna_raw_text"] = raw_text
        if not resolved.strip():
            st.error("Provide an RNA sequence or a FASTA file to analyze.")
            return

    raw = st.session_state.get("rna_input", "")
    if not raw:
        return
    _warn_if_analyzed_input_changed(text_value, raw, action="Analyze")

    info = _validate_for(raw, "RNA")
    st.markdown(_type_badge(info), unsafe_allow_html=True)

    dna_pasted = _normalized_seq(str(st.session_state.get("dna_input", "")))
    rna_pasted = _normalized_seq(raw)
    if dna_pasted and rna_pasted and dna_pasted == rna_pasted:
        st.warning(
            "The string in this tab is identical to the DNA tab input. "
            "The RNA tab does not import DNA-tab results and does not transcribe DNA."
        )

    if not info["is_valid"]:
        st.error(info["rejection_reason"] or "Provide a valid RNA sequence (A, U, C, G only).")
        return

    rna_seq = info["sequence"]
    _store_workspace("RNA", rna_seq, source="user sequence", identifier="")
    dna_view = rna_seq.replace("U", "T")
    gc = _gc_content(dna_view)
    st.markdown(badge("Analysis molecule: RNA as pasted", "rna"), unsafe_allow_html=True)
    st.markdown(status_badge("COMPUTED"), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.result_header_html(
            "RNA analysis",
            status="COMPUTED",
            lines=[f"{len(rna_seq):,} nt", f"GC {_format_percent(gc)}%"],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation(
        "rna",
        {
            "status": "COMPUTED",
            "length": len(rna_seq),
            "gc_percent": gc,
            "method": "RNA composition",
            "source": "HelixScope RNA analysis",
        },
    )
    at = _at_content(dna_view)
    rna_entropy = _rna_shannon_entropy(rna_seq)

    rna_entropy_label = (
        "N/A"
        if isinstance(rna_entropy, float) and math.isnan(rna_entropy)
        else f"{rna_entropy:.4f}"
    )
    render_metric_grid(
        [
            ("Length", str(len(rna_seq)), "nt"),
            (
                "GC Content",
                _format_percent(gc),
                "%" if _format_percent(gc) != "Unavailable" else "",
            ),
            (
                "AU Content",
                _format_percent(at),
                "%" if _format_percent(at) != "Unavailable" else "",
            ),
            (
                "Shannon entropy",
                rna_entropy_label,
                "bits" if rna_entropy_label != "N/A" else "",
            ),
        ]
    )
    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>GC and AU are computed on "
        "this tab's sequence after U is treated as T in the formula. Identical "
        "GC% to the DNA tab means the same base composition after T/U "
        "substitution, not that DNA results were copied.</span>",
        unsafe_allow_html=True,
    )

    resolution = _render_cds_resolver(rna_seq)
    fold_result = None
    if _section_toggle("Secondary Structure", key="rna_structure"):
        with st.container(border=True):
            fold_result = _render_rna_structure_section(rna_seq)
    if _section_toggle("RNA 3D", key="rna_3d"):
        with st.container(border=True):
            _render_rna_3d_section(rna_seq, fold_result=fold_result)
    if resolution is None:
        st.info(
            "Codon-usage sections require a resolved coding region. "
            "Predicted RNA secondary structure does not."
        )
        return
    coding_codons = str(resolution["codons"])

    dinucleotides = _dinucleotides(dna_view)

    if _section_toggle("Nucleotide Composition", key="rna_nuc_comp"):
        with st.container(border=True):
            composition = _composition(rna_seq)
            st.plotly_chart(
                charts.nucleotide_bar_chart(composition), width="stretch"
            )
            comp_rows = [
                {
                    "Symbol": symbol,
                    "Count": data["count"],
                    "Frequency_pct": data["frequency"],
                }
                for symbol, data in composition.items()
                if data["count"] > 0
            ]
            st.dataframe(pd.DataFrame(comp_rows), width="stretch", hide_index=True)

    if _section_toggle(
        "Dinucleotide Frequencies and CpG Odds Ratio", key="rna_dinuc"
    ):
        with st.container(border=True):
            st.plotly_chart(
                charts.dinucleotide_heatmap(dinucleotides), width="stretch"
            )
            dinuc_rows = [
                {
                    "Dinucleotide": pair,
                    "Count": data["count"],
                    "Frequency_pct": data["frequency"],
                    "Observed_Expected": data["observed_expected"],
                }
                for pair, data in dinucleotides.items()
            ]
            st.dataframe(
                pd.DataFrame(dinuc_rows).sort_values(
                    by="Observed_Expected", ascending=False
                ),
                width="stretch",
                hide_index=True,
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Counted on the DNA "
                "equivalent of the sequence, so U appears as T.</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("GC Content in Sliding Window", key="rna_gc_window"):
        with st.container(border=True):
            plan = scale_profile.analysis_plan(len(rna_seq), molecule="RNA")
            max_window = max(10, min(2000, len(rna_seq)))
            default_window = min(max_window, int(plan["suggested_window"]))
            default_step = min(max_window, max(1, int(plan["suggested_step"])))
            win_col, step_col = st.columns(2)
            with win_col:
                gc_window = st.number_input(
                    "Window size (nt)",
                    min_value=10,
                    max_value=max_window,
                    value=default_window,
                    step=10,
                    key="rna_gc_window_size",
                )
            with step_col:
                gc_step = st.number_input(
                    "Step (nt)",
                    min_value=1,
                    max_value=max_window,
                    value=default_step,
                    step=5,
                    key="rna_gc_window_step",
                )
            st.caption(
                f"Display/memory cap is {scale_profile.MAX_WINDOW_PROFILE_ROWS:,} windows. "
                "ViennaRNA MFE remains limited to "
                f"{rna_folding.MAX_FOLD_NT} nt; that is a folding limit, not an analysis limit."
            )
            try:
                profile = _gc_sliding_window(dna_view, int(gc_window), int(gc_step))
                st.plotly_chart(
                    charts.gc_sliding_window_plot(profile), width="stretch"
                )
            except ValueError as exc:
                st.info(str(exc))

    if _section_toggle("Codon Usage Table", key="rna_codon_usage"):
        with st.container(border=True):
            try:
                codon_df = _codon_usage_table(coding_codons)
                st.plotly_chart(
                    charts.codon_usage_bar_chart(codon_df, "Absolute_pct"),
                    width="stretch",
                )
                st.plotly_chart(
                    charts.codon_usage_heatmap(codon_df, "Synonymous_pct"),
                    width="stretch",
                )
                st.dataframe(codon_df, width="stretch", hide_index=True)
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);font-size:12px;'>Synonymous_pct is "
                    "the share of the codon within its amino acid family. "
                    "Absolute_pct is its share of all codons counted. Both come "
                    f"from the {resolution['n_codons']:,} codons of "
                    f"{resolution['n_cds']:,} coding regions, never from the raw "
                    "sequence.</span>",
                    unsafe_allow_html=True,
                )
                st.download_button(
                    "Download codon usage (CSV)",
                    data=codon_df.to_csv(index=False).encode("utf-8"),
                    file_name="codon_usage.csv",
                    mime="text/csv",
                    key="rna_codon_download",
                )
            except ValueError as exc:
                st.info(str(exc))

    if _section_toggle("Relative Synonymous Codon Usage (RSCU)", key="rna_rscu"):
        with st.container(border=True):
            try:
                rscu_df = _rscu_table(coding_codons)
                st.plotly_chart(charts.rscu_bar_chart(rscu_df), width="stretch")
                st.dataframe(rscu_df, width="stretch", hide_index=True)
                st.download_button(
                    "Download RSCU (CSV)",
                    data=rscu_df.to_csv(index=False).encode("utf-8"),
                    file_name="rscu.csv",
                    mime="text/csv",
                    key="rna_rscu_download",
                )
            except ValueError as exc:
                st.info(str(exc))

    enc_value: Optional[float] = None
    cai_value: Optional[float] = None
    organism = "human"
    if _section_toggle("Codon Bias Indices (ENC and CAI)", key="rna_bias"):
        with st.container(border=True):
            organism = st.selectbox(
                "CAI reference organism",
                ["human", "ecoli", "yeast"],
                key="rna_organism",
            )
            st.markdown(
                badge(f"CAI reference: {organism} codon usage table", "rna"),
                unsafe_allow_html=True,
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>The reference table is "
                "an approximate Kazusa-style codon usage table embedded in the "
                "project, expressed as frequency per thousand codons. CAI is not a "
                "universal property of the sequence: it depends on this table. No "
                "organism is inferred from the sequence.</span>",
                unsafe_allow_html=True,
            )

            enc_report = _effective_number_of_codons_report(coding_codons)
            if enc_report["status"] == "COMPUTED":
                enc_value = float(enc_report["value"])
                enc_shown = f"{enc_value:.2f}"
            else:
                enc_shown = "N/A"
                st.info(enc_report.get("reason") or "ENC is not applicable.")
            cai_shown = "N/A"
            try:
                cai_report = _codon_adaptation_index_report(coding_codons, organism)
                if cai_report["status"] == "HEURISTIC" and not (
                    isinstance(cai_report["value"], float)
                    and math.isnan(cai_report["value"])
                ):
                    cai_value = float(cai_report["value"])
                    cai_shown = f"{cai_value:.4f}"
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);font-size:12px;'>CAI table: "
                    f"{html_escape(str(cai_report['table_source']))} Method: "
                    f"{html_escape(str(cai_report['method']))} Informative "
                    f"codons: {cai_report['n_informative']}.</span>",
                    unsafe_allow_html=True,
                )
            except ValueError as exc:
                st.info(str(exc))
            render_metric_grid(
                [
                    ("ENC (Nc)", enc_shown, ""),
                    ("CAI", cai_shown, ""),
                ]
            )

            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>ENC is Wright (1990) "
                "Nc, applicable only when synonymous families allow F. Range 20 to "
                "61 when computed. CAI is a heuristic inspired by Sharp and Li "
                "(1987) against the selected table, not an experimental expression "
                f"measurement. Codons: {resolution['n_codons']:,} from "
                f"{resolution['method_label']}.</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("Translated Protein Composition", key="rna_translated"):
        with st.container(border=True):
            try:
                translated = _translate(coding_codons)
            except (ValueError, RuntimeError) as exc:
                translated = ""
                st.info(str(exc))
            if translated:
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);font-size:12px;'>Translation of the "
                    "concatenated coding regions, stopping at the first stop codon "
                    "of the concatenation. It summarises amino acid usage, not the "
                    "sequence of any individual protein.</span>",
                    unsafe_allow_html=True,
                )
                aa_composition = _amino_acid_composition(translated)
                st.plotly_chart(
                    charts.amino_acid_bar_chart(aa_composition), width="stretch"
                )
                aa_rows = [
                    {
                        "Amino acid": symbol,
                        "Count": data["count"],
                        "Frequency_pct": data["frequency"],
                    }
                    for symbol, data in aa_composition.items()
                    if data["count"] > 0
                ]
                st.dataframe(
                    pd.DataFrame(aa_rows), width="stretch", hide_index=True
                )

    if _section_toggle("Resolved Coding Regions", key="rna_regions"):
        with st.container(border=True):
            region_rows = [
                {
                    "Label": region.get("label")
                    or region.get("product")
                    or region.get("gene")
                    or "CDS",
                    "Gene": region.get("gene", ""),
                    "Product": region.get("product", ""),
                    "Strand": region.get("strand", "+"),
                    "Start_0based": region.get("start", region["segments"][0][0]),
                    "End_exclusive": region.get("end", region["segments"][-1][1]),
                    "Length_nt": region.get(
                        "length_nt",
                        sum(int(s[1]) - int(s[0]) for s in region["segments"]),
                    ),
                    "Segments": len(region["segments"]),
                }
                for region in resolution["regions"]
            ]
            region_df = pd.DataFrame(region_rows)
            st.dataframe(region_df, width="stretch", hide_index=True)
            st.caption("Start_0based / End_exclusive are 0-based half-open coordinates on the pasted RNA.")
            st.plotly_chart(
                charts.orf_length_histogram(region_df["Length_nt"].tolist()),
                width="stretch",
            )
            st.download_button(
                "Download coding regions (CSV)",
                data=region_df.to_csv(index=False).encode("utf-8"),
                file_name="coding_regions.csv",
                mime="text/csv",
                key="rna_regions_download",
            )

    if _section_toggle("Sequence Viewer", key="rna_seq_view"):
        with st.container(border=True):
            _sequence_viewer_panel(rna_seq, "RNA", key_prefix="rna_view")

    if _section_toggle("RNA k-mers", key="rna_kmers"):
        with st.container(border=True):
            rna_k = st.number_input(
                "k-mer length",
                min_value=1,
                max_value=5,
                value=3,
                key="rna_kmer_k",
            )
            rna_kmers = _kmer_counts_rna(rna_seq, int(rna_k))
            if rna_kmers:
                st.plotly_chart(charts.kmer_bar_chart(rna_kmers), width="stretch")
            else:
                st.info("Insufficient data for RNA k-mer counts.")

    if fold_result:
        structure_summary = "PREDICTED (ViennaRNA MFE)"
        structure_report = {
            "status": fold_result.get("status"),
            "tool": fold_result.get("tool"),
            "mfe_kcal_mol": fold_result.get("mfe_kcal_mol"),
            "sequence_hash": fold_result.get("sequence_hash"),
        }
    elif rna_folding.tool_availability()["available"]:
        structure_summary = "Not computed (click Predict in Secondary Structure)"
        structure_report = "NOT_COMPUTED"
    else:
        structure_summary = "UNAVAILABLE"
        structure_report = "UNAVAILABLE"

    summary_rows = [
        {"Metric": "Length", "Value": len(rna_seq), "Unit": "nt"},
        {"Metric": "GC content", "Value": gc, "Unit": "%"},
        {"Metric": "AU content", "Value": at, "Unit": "%"},
        {"Metric": "Shannon entropy", "Value": rna_entropy, "Unit": "bits"},
        {"Metric": "RNA secondary structure", "Value": structure_summary, "Unit": ""},
        {
            "Metric": "RNA 3D coordinates",
            "Value": "UNAVAILABLE (ViennaRNA MFE is not a 3D model)",
            "Unit": "",
        },
        {
            "Metric": "CpG observed/expected",
            "Value": dinucleotides.get("CG", {}).get("observed_expected", float("nan")),
            "Unit": "ratio",
        },
        {"Metric": "CDS resolution method", "Value": resolution["method_label"], "Unit": ""},
        {"Metric": "Coding regions", "Value": resolution["n_cds"], "Unit": "count"},
        {"Metric": "Codons counted", "Value": resolution["n_codons"], "Unit": "count"},
        {"Metric": "Coding coverage", "Value": resolution["coverage_pct"], "Unit": "%"},
    ]
    if enc_value is not None:
        summary_rows.append({"Metric": "ENC (Nc)", "Value": enc_value, "Unit": "codons"})
    if cai_value is not None:
        summary_rows.append(
            {"Metric": f"CAI ({organism})", "Value": cai_value, "Unit": "0 to 1"}
        )
    _render_summary_table(summary_rows, "rna_summary_download", "rna_summary.csv")
    rna_report = provenance.analysis_envelope(
        module="RNA",
        payload={
            "length_nt": len(rna_seq),
            "gc_percent": gc,
            "au_percent": at,
            "shannon_entropy_bits": rna_entropy,
            "cds_method": resolution["method_label"],
            "n_cds": resolution["n_cds"],
            "n_codons": resolution["n_codons"],
            "coding_coverage_pct": resolution["coverage_pct"],
            "enc": enc_value,
            "cai": cai_value,
            "cai_organism": organism,
            "rna_secondary_structure": structure_report,
        },
        status="COMPUTED",
        algorithm="HelixScope RNA composition and codon usage",
        parameters={
            "cds_source": resolution.get("source"),
            "cai_organism": organism,
        },
        sequence=rna_seq,
        source="user sequence",
    )
    st.download_button(
        "Download analysis report (JSON)",
        data=json.dumps(rna_report, indent=2).encode("utf-8"),
        file_name="helixscope_rna_report.json",
        mime="application/json",
        key="rna_report_json",
    )
    _render_alerts(
        _nucleotide_alerts(
            gc, dinucleotides, enc_value, cai_value, organism, resolution
        )
    )


def _na_text(value: object) -> str:
    """Mostra N/A para ausente; nao converte ausente em 0 (interno)."""
    if value is None or value == "":
        return "N/A"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "N/A"
    return str(value)


def _msa_conservation_overlay_for_protein(protein_seq: str) -> Optional[dict]:
    """Overlay de conservacao MSA apenas se o hash da proteina estiver no MSA."""
    result = st.session_state.get("msa_result")
    if not isinstance(result, dict):
        return None
    digest = provenance.sequence_digest(protein_seq)
    hashes = [str(item) for item in list(result.get("input_hashes") or [])]
    if digest not in hashes:
        return None
    rows = list(result.get("rows") or [])
    maps = list(result.get("coordinate_maps") or [])
    scores = list((result.get("conservation") or {}).get("scores") or [])
    row_index = None
    for index, row in enumerate(rows):
        if str(row.get("hash") or "") == digest:
            row_index = index
            break
    if row_index is None and hashes:
        try:
            row_index = hashes.index(digest)
        except ValueError:
            return None
    if row_index is None or row_index >= len(maps):
        return None
    overlay: dict = {}
    for column, position in enumerate(maps[row_index]):
        if position is None or column >= len(scores):
            continue
        overlay[int(position)] = scores[column]
    return overlay or None


def _selection_fingerprint(selection: object) -> str:
    """Identifica uma MoleculeSelection sem recalcular ciencia."""
    if not isinstance(selection, dict):
        return ""
    return molecule_selection.selection_identity_hash(selection)


def _session_nucleic_structure_envelope() -> Optional[dict]:
    """Envelope DNA/RNA depositado ou ilustrativo ja validado nesta sessao."""
    for key in ("dna_struct_result", "rna_struct_result"):
        item = st.session_state.get(key)
        if isinstance(item, dict) and item.get("sequence_hash"):
            return item
    for key in (
        "dna_3d_representation",
        "rna_3d_representation",
        "dna_illustrative",
        "rna_illustrative",
    ):
        item = st.session_state.get(key)
        if not isinstance(item, dict):
            continue
        envelope = item.get("envelope")
        if isinstance(envelope, dict) and envelope.get("sequence_hash"):
            return envelope
        if item.get("sequence_hash"):
            return item
    return None


def _store_molecule_selection(selection: dict) -> None:
    """Guarda a selecao partilhada e forca uma aplicacao unica ao inspector."""
    st.session_state["molecule_selection"] = selection
    st.session_state.pop("molecule_selection_applied", None)


def _active_molecule_selection(
    protein_seq: str,
    structure_result: Optional[dict] = None,
) -> Optional[dict]:
    """Devolve a selecao se for da mesma proteina (e estrutura, se conhecida)."""
    from modules import molecule_selection

    selection = st.session_state.get("molecule_selection")
    structure_id = None
    structure_hash = None
    if isinstance(structure_result, dict):
        structure_id = structure_result.get("structure_id")
        structure_hash = structure_result.get("content_hash")
    if not molecule_selection.selection_applies_to(
        selection,
        sequence_hash=provenance.sequence_digest(protein_seq),
        molecule="PROTEIN",
        structure_id=structure_id,
        structure_hash=structure_hash,
    ):
        return None
    return dict(selection)


def _apply_pending_molecule_selection(protein_seq: str) -> None:
    """Aplica uma selecao nova ao indice de inspecao uma vez, antes do widget."""
    from modules import molecule_selection

    selection = st.session_state.get("molecule_selection")
    if not molecule_selection.selection_applies_to(
        selection,
        sequence_hash=provenance.sequence_digest(protein_seq),
        molecule="PROTEIN",
    ):
        return
    fingerprint = _selection_fingerprint(selection)
    if st.session_state.get("molecule_selection_applied") == fingerprint:
        return
    primary = molecule_selection.primary_index(selection)
    if primary is not None and 0 <= primary < len(protein_seq):
        st.session_state["prot_struct_pos"] = int(primary)
    st.session_state["molecule_selection_applied"] = fingerprint


def _protein_highlight_indices(
    protein_seq: str,
    selected: int,
    structure_result: Optional[dict],
) -> tuple[int, list[int], Optional[dict]]:
    """Sincroniza inspecao e highlights. Nao inventa residuos."""
    from modules import molecule_selection

    selection = _active_molecule_selection(protein_seq, None)
    if selection is None:
        return selected, [selected], None
    status = str(selection.get("status") or "")
    if status == "UNMAPPED":
        return selected, [], selection
    if status not in {"READY", ""}:
        return selected, [selected], selection
    indices = [
        index
        for index in molecule_selection.query_indices_of(selection)
        if 0 <= index < len(protein_seq)
    ]
    if indices and selected not in indices:
        rewritten = molecule_selection.selection_from_sequence_position(
            molecule="PROTEIN",
            sequence=protein_seq,
            position=selected,
        )
        if isinstance(structure_result, dict):
            rewritten["structure_id"] = structure_result.get("structure_id")
            rewritten["structure_hash"] = structure_result.get("content_hash")
        _store_molecule_selection(rewritten)
        st.session_state["molecule_selection_applied"] = _selection_fingerprint(rewritten)
        return selected, [selected], rewritten
    structure_selection = _active_molecule_selection(protein_seq, structure_result)
    return selected, (indices or [selected]), structure_selection or selection


def _atom_record_for_residue(result: dict, residue: Mapping[str, Any], atom_name: str = "CA") -> Optional[dict]:
    """Localiza um ATOM depositado para o residuo mapeado. Nao inventa coordenadas."""
    wanted_chain = str(residue.get("chain_id") or "")
    wanted_label = residue.get("label_seq_id")
    wanted_auth = residue.get("auth_seq_id")
    wanted_model = residue.get("model")
    wanted_ins = str(residue.get("insertion_code") or "")
    for atom in list(result.get("atoms") or []):
        if str(atom.get("group") or "ATOM") != "ATOM":
            continue
        if str(atom.get("atom_name") or "") != str(atom_name):
            continue
        chain = str(atom.get("auth_asym_id") or atom.get("label_asym_id") or "")
        if wanted_chain and chain != wanted_chain:
            continue
        if wanted_model is not None and int(atom.get("model") or 1) != int(wanted_model):
            continue
        if wanted_label is not None and atom.get("label_seq_id") != wanted_label:
            continue
        if wanted_auth is not None and atom.get("auth_seq_id") != wanted_auth:
            continue
        if str(atom.get("insertion_code") or "") != wanted_ins:
            continue
        return dict(atom)
    return None


def _render_structure_inspector(
    protein_seq: str,
    result: dict,
    selected: int,
    selection: Optional[dict],
) -> None:
    """Painel de dados complementar ao viewport. Nao e tooltip."""
    from modules import molecule_selection

    st.caption("Structure inspector (data panel; values are not hover-only)")
    residue = protein_structure.query_position_to_residue(result, selected)
    source = str(result.get("source") or "")
    kind = str(result.get("kind") or "")
    mapping_status = str(result.get("mapping_status") or "UNCERTAIN")
    selection_source = str((selection or {}).get("source") or "sequence")
    feature_label = str((selection or {}).get("feature_label") or (selection or {}).get("label") or "")
    atom_name = str((selection or {}).get("atom_name") or "CA")
    atom = None
    if residue is not None:
        atom = _atom_record_for_residue(result, residue, atom_name if atom_name else "CA")
        if atom is None:
            atom = _atom_record_for_residue(result, residue, "CA")
    ca = (residue or {}).get("ca") or {}
    if atom:
        xyz_text = f"{atom.get('x')}, {atom.get('y')}, {atom.get('z')}"
        occupancy_text = _na_text(atom.get("occupancy"))
        element_text = _na_text(atom.get("element"))
        atom_label = _na_text(atom.get("atom_name"))
    elif residue and residue.get("has_coordinates") and ca:
        xyz_text = f"{ca.get('x')}, {ca.get('y')}, {ca.get('z')}"
        occupancy_text = "N/A"
        element_text = "N/A"
        atom_label = "CA"
    else:
        xyz_text = "Coordinates unavailable"
        occupancy_text = "N/A"
        element_text = "N/A"
        atom_label = "N/A"
    query_aa = protein_seq[selected] if 0 <= selected < len(protein_seq) else ""
    sections = {
        "IDENTITY": [
            ("Residue", query_aa or _na_text((residue or {}).get("query_residue"))),
            ("Internal position (0-based)", selected),
            ("Chain", _na_text((residue or {}).get("chain_id"))),
            ("label_seq_id", _na_text((residue or {}).get("label_seq_id"))),
            ("auth_seq_id", _na_text((residue or {}).get("auth_seq_id"))),
            ("Insertion", (residue or {}).get("insertion_code") or ""),
        ],
        "SOURCE": [
            ("Source", source),
            ("Kind", kind),
            ("Selection source", selection_source),
            ("Feature", feature_label or "N/A"),
            ("Feature source", str((selection or {}).get("feature_source") or "")),
        ],
        "STRUCTURE": [
            ("Coordinates", xyz_text),
            ("Atom", atom_label),
            ("Element", element_text),
            ("Occupancy", occupancy_text),
            ("Model", _na_text((residue or {}).get("model") or result.get("model_number"))),
        ],
        "MAPPING": [
            ("Mapping", _na_text((residue or {}).get("mapping_row_status") or mapping_status)),
            ("Has coordinates", (residue or {}).get("has_coordinates")),
        ],
        "ANALYSIS": [
            ("Selected query index", selected),
        ],
        "PROVENANCE": [
            ("sequence_hash", result.get("sequence_hash")),
            ("structure_hash", result.get("content_hash") or result.get("structure_hash")),
            ("structure_id", result.get("structure_id")),
        ],
    }
    st.markdown(inspector_panel(sections), unsafe_allow_html=True)
    if residue is None:
        st.info("No mapped structural residue")
    elif not residue.get("has_coordinates"):
        st.info("Coordinates unavailable for this residue")
    if selection and str(selection.get("status") or "") == "UNMAPPED":
        st.info(str(selection.get("message") or molecule_selection.GAP_MESSAGE))
    if mapping_status == "BEST_EFFORT":
        st.caption("Mapping status is BEST_EFFORT. This is not exact sequence mapping.")
    st.caption(
        "Internal position, label_seq_id and auth_seq_id are shown separately. "
        "They are not interchangeable. Lighting and materials are visual, not "
        "confidence, activity or importance."
    )


def _render_msa_protein_3d_link(msa_result: dict, column: int) -> None:
    """MSA column -> selecao proteica. Gap nao cria residuo 3D."""
    from modules import molecule_selection

    protein_seq = str(st.session_state.get("protein_input") or "")
    if not protein_seq:
        st.caption(
            "Analyze a matching protein in the Protein tab to highlight this "
            "MSA column on 3D. A column is not a structure."
        )
        return
    digest = provenance.sequence_digest(protein_seq)
    hashes = [str(item) for item in list(msa_result.get("input_hashes") or [])]
    if digest not in hashes:
        st.caption(
            "This MSA does not include the current Protein tab sequence. "
            "No 3D highlight is created."
        )
        return
    if st.button("Highlight this MSA column on protein 3D", key="msa_3d_highlight"):
        structure = st.session_state.get("prot_struct_result")
        selection = molecule_selection.selection_from_msa_column(
            msa_result,
            int(column),
            sequence=protein_seq,
            molecule="PROTEIN",
            structure_result=structure if isinstance(structure, dict) else None,
        )
        _store_molecule_selection(selection)
        if str(selection.get("status") or "") == "UNMAPPED":
            st.info(str(selection.get("message") or molecule_selection.GAP_MESSAGE))
        elif str(selection.get("status") or "") == "STALE":
            st.warning(str(selection.get("message") or "stale/incompatible structure"))
        else:
            st.caption(
                "MSA column mapped to the analysis sequence. Open Protein Structure "
                "and View 3D Structure to see the highlight. Conservation is among "
                "analyzed sequences, not a functional residue annotation."
            )
            st.rerun()


def _render_structure_error(exc: Exception) -> None:
    """Exibe falha estrutural classificada, nunca como coordenadas desenhadas."""
    if isinstance(exc, protein_structure.StructureError):
        label = _FOLDING_ERROR_LABELS.get(exc.category, exc.category)
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        if exc.category == "NO_STRUCTURE":
            st.info("No validated structure found.")
            st.caption(str(exc))
            return
        st.error(f"Structure retrieval failed. {label}. {exc}")
        st.caption(f"category={exc.category} source=structure_retrieval")
        return
    if isinstance(exc, complex_structure.ComplexError):
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(f"Complex structure failed. {exc.category}. {exc}")
        return
    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
    st.error(f"Structure retrieval failed. {exc}")


def _render_protein_structure_3d(
    protein_seq: str,
    result: dict,
    selected: int,
    highlight_indices: Optional[list[int]] = None,
    selection: Optional[dict] = None,
) -> None:
    """Renderer 3D lazy. Isolado: falha nao derruba a analise de proteina."""
    if not st.toggle("View 3D Structure", value=False, key="prot_view_3d"):
        st.caption(
            "3D is not loaded until View 3D Structure is enabled. "
            "Opening the Protein tab does not start the renderer. "
            "Sequence, mapping table and structure metadata remain available above."
        )
        _render_structure_inspector(protein_seq, result, selected, selection)
        return
    current_hash = provenance.sequence_digest(protein_seq)
    if str(result.get("sequence_hash") or "") != current_hash:
        st.markdown(status_badge("STALE"), unsafe_allow_html=True)
        st.warning("stale/incompatible structure")
        return
    try:
        from modules import structure_scene
        from ui import structure_viewer
    except Exception:
        st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
        st.error("3D unavailable in this browser/environment.")
        st.caption("Sequence, mapping table and structure metadata remain available above.")
        return
    kind = str(result.get("kind") or "")
    statuses = [kind.upper() or "UNAVAILABLE"]
    st.markdown(badge_row(statuses), unsafe_allow_html=True)
    st.caption(protein_structure.structure_kind_label(kind))
    method_text = _na_text(result.get("method"))
    resolution = result.get("resolution_angstrom")
    if resolution is None or (
        isinstance(resolution, float) and (math.isnan(resolution) or math.isinf(resolution))
    ):
        resolution_text = "N/A"
    else:
        try:
            resolution_text = f"{float(resolution):.2f}"
        except (TypeError, ValueError):
            resolution_text = "N/A"
    st.markdown(
        meta_grid(
            [
                ("Structure ID", _na_text(result.get("structure_id"))),
                ("Source", _na_text(result.get("source"))),
                ("Kind", _na_text(kind)),
                ("Method", method_text),
                ("Resolution (A)", resolution_text),
                ("Chain", _na_text((result.get("selected_chain") or {}).get("chain_id"))),
                ("Model", _na_text(result.get("model_number"))),
                ("Coverage", _na_text(result.get("coverage_query_with_coordinates"))),
                ("sequence_hash", _na_text(result.get("sequence_hash"))),
            ]
        ),
        unsafe_allow_html=True,
    )
    if kind == "experimental":
        st.caption("Experimental structure. This is not an exact molecular structure claim.")
    if kind == "predicted":
        st.caption("Predicted Structure. This is not an experimental measurement.")
    if kind == "illustrative":
        st.caption("Illustrative representation. Coordinates are not an experimental measurement.")
    st.caption(
        "This is a static 3D visualization of deposited coordinates, not a "
        "molecular simulation and not molecular dynamics."
    )
    models = list(result.get("models") or [result.get("model_number") or 1])
    model_options = [int(item) for item in models]
    if len(model_options) > 1:
        st.caption(
            "Multiple NMR models are present. The displayed model is a selected "
            "conformer, not the true structure."
        )
        model_filter = int(
            st.selectbox(
                "NMR model to visualize",
                model_options,
                index=model_options.index(int(result.get("model_number") or model_options[0]))
                if int(result.get("model_number") or model_options[0]) in model_options
                else 0,
                key="prot_3d_model",
            )
        )
    else:
        model_filter = int(result.get("model_number") or 1)
    chain_ids = []
    for item in list(result.get("chains") or []):
        cid = str(item.get("chain_id") or "")
        if cid and cid not in chain_ids:
            chain_ids.append(cid)
    default_chain = str((result.get("selected_chain") or {}).get("chain_id") or "")
    if len(chain_ids) > 1:
        chain_filter = str(
            st.selectbox(
                "Chain to visualize",
                chain_ids,
                index=chain_ids.index(default_chain) if default_chain in chain_ids else 0,
                key="prot_3d_chain",
            )
        )
    else:
        chain_filter = default_chain
    visible_chains = None
    if len(chain_ids) > 1:
        visible_chains = st.multiselect(
            "Visible chains",
            chain_ids,
            default=chain_ids,
            key="prot_3d_visible",
        )
        if chain_filter and chain_filter not in visible_chains:
            st.caption(
                "Mapped chain is hidden. Backbone of the mapped polymer is not drawn. "
                "This is a visual choice and does not change the scientific structure."
            )
        if not visible_chains:
            st.info("Coordinates unavailable. Select at least one chain.")
            _render_structure_inspector(protein_seq, result, selected, selection)
            return
    representation = st.selectbox(
        "Representation",
        ["backbone", "ca", "atoms"],
        format_func=lambda item: {
            "backbone": "Backbone (CA trace)",
            "ca": "CA atoms + backbone",
            "atoms": "All ATOM records",
        }[item],
        key="prot_3d_rep",
    )
    st.caption(
        "Cartoon ribbon (DSSP/STRIDE) and molecular surface are UNAVAILABLE. "
        "Backbone is a geometric CA trace of deposited coordinates."
    )
    lod = _lod_select("prot_3d_lod")
    color_options = ["chain", "residue_type", "selection", "hydropathy", "charge"]
    if kind == "predicted":
        color_options.append("plddt")
    if _msa_conservation_overlay_for_protein(protein_seq) is not None:
        color_options.append("conservation")
    color_mode = st.selectbox(
        "Color mode",
        color_options,
        format_func=lambda item: {
            "chain": "Chain",
            "residue_type": "Residue type",
            "selection": "Selection",
            "plddt": "pLDDT (deposited B-factor, predicted only)",
            "hydropathy": "Hydropathy (Kyte-Doolittle)",
            "charge": "Formal side-chain charge (K/R D/E, not net_charge)",
            "conservation": "MSA conservation among analyzed sequences",
        }[item],
        key="prot_3d_color",
    )
    overlay = None
    if color_mode == "hydropathy":
        overlay = {
            index: protein_analysis.KYTE_DOOLITTLE.get(residue)
            for index, residue in enumerate(protein_seq)
            if residue in protein_analysis.KYTE_DOOLITTLE
        }
    elif color_mode == "conservation":
        overlay = _msa_conservation_overlay_for_protein(protein_seq)
    elif color_mode == "charge":
        charges = protein_analysis.per_residue_formal_charge(protein_seq)
        overlay = {index: float(value) for index, value in enumerate(charges)}
    cam_cols = st.columns(3)
    with cam_cols[0]:
        fit = st.button("Fit", key="prot_3d_fit")
    with cam_cols[1]:
        reset = st.button("Reset", key="prot_3d_reset")
    with cam_cols[2]:
        center_sel = st.button("Center selection", key="prot_3d_center")
    if fit or reset:
        st.session_state["prot_3d_cam"] = int(st.session_state.get("prot_3d_cam") or 0) + 1
        st.session_state["prot_3d_center_sel"] = False
    if center_sel:
        st.session_state["prot_3d_center_sel"] = True
    try:
        envelope = protein_structure.structure_3d_input(result)
        scene_highlights = list(highlight_indices) if highlight_indices is not None else [selected]
        selected_for_scene = selected if scene_highlights or selected is not None else None
        if selection and str(selection.get("status") or "") == "UNMAPPED":
            selected_for_scene = None
            scene_highlights = []
        scene = structure_scene.build_scene(
            envelope,
            expected_sequence_hash=current_hash,
            representation=str(representation),
            color_mode=str(color_mode),
            selected_query_index=selected_for_scene,
            highlight_query_indices=scene_highlights,
            chain_filter=chain_filter,
            model_filter=model_filter,
            overlay=overlay,
            visible_chains=visible_chains,
            lod=lod,
        )
    except Exception:
        st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
        st.error("3D unavailable in this browser/environment.")
        st.caption("Sequence, mapping table and structure metadata remain available above.")
        _render_structure_inspector(protein_seq, result, selected, selection)
        return
    status = str(scene.get("status") or "ERROR")
    st.markdown(status_badge(status), unsafe_allow_html=True)
    if scene.get("partial_view"):
        st.markdown(status_badge("PARTIAL_VIEW"), unsafe_allow_html=True)
        st.caption("Partial structure view. This is not the complete deposited file.")
    else:
        st.caption(str(scene.get("view_label") or "Selected chain"))
    if scene.get("mapping_status"):
        st.caption(
            f"Mapping status: {scene.get('mapping_status')}. "
            f"{scene.get('mapping_status_note') or ''}"
        )
        if str(scene.get("mapping_status") or "") == "BEST_EFFORT":
            st.caption("BEST_EFFORT is not exact sequence mapping.")
    st.caption(str(scene.get("representation_note") or ""))
    st.caption(str(scene.get("color_note") or ""))
    st.caption(str(scene.get("webgl_note") or ""))
    if scene.get("lod_note"):
        st.caption(str(scene.get("lod_note")))
    timings = scene.get("timings_ms") or {}
    if timings.get("scene_construction") is not None:
        st.caption(
            f"Scene construction: {timings.get('scene_construction')} ms on this machine. "
            "Not an FPS claim and not a GPU-acceleration claim."
        )
    st.caption(
        "Rotate, pan and zoom use the Plotly viewport. Fit and Reset restore the "
        "camera. Center selection moves the view to the selected residue when "
        "coordinates exist. Previous/Next residue navigate without requiring a 3D click."
    )
    expanded = st.toggle(
        "Expanded viewport",
        value=False,
        key="prot_expanded_3d",
        help="Taller canvas. Does not change coordinates.",
    )
    from ui.tokens import STRUCTURE_VIEWPORT_HEIGHT, STRUCTURE_VIEWPORT_HEIGHT_EXPANDED

    chart_height = STRUCTURE_VIEWPORT_HEIGHT_EXPANDED if expanded else STRUCTURE_VIEWPORT_HEIGHT
    view_col, insp_col = st.columns([3, 2], gap="medium")
    with insp_col:
        _render_structure_inspector(protein_seq, result, selected, selection)
    if status != "READY":
        with view_col:
            if status == "UNAVAILABLE":
                st.info("Coordinates Unavailable")
            elif status == "STALE":
                st.warning("stale/incompatible structure")
            elif status == "RESOURCE_LIMIT":
                st.error(str(scene.get("message") or "RESOURCE_LIMIT"))
            else:
                st.error(str(scene.get("message") or "3D view failed. Structural tables remain available."))
        return
    try:
        fig = structure_viewer.figure_from_scene(
            scene,
            center_selection=bool(st.session_state.get("prot_3d_center_sel")),
        )
        try:
            fig.update_layout(height=chart_height)
        except Exception:
            pass
        chart_key = f"prot_3d_chart_{int(st.session_state.get('prot_3d_cam') or 0)}"
        with view_col:
            render_hand_control(
                plot_key=chart_key,
                structure_id=str(
                    result.get("structure_id") or result.get("content_hash") or "protein"
                ),
                illustrative=False,
                control_key="prot_3d",
            )
            st.markdown('<div class="hs-viewport">', unsafe_allow_html=True)
            event = _structure_plotly_chart(
                fig,
                key=chart_key,
                on_select="rerun",
                height=chart_height,
            )
            st.markdown("</div>", unsafe_allow_html=True)
        parsed = structure_viewer.parse_plotly_selection(event)
        if parsed:
            from modules import molecule_selection

            pick = molecule_selection.selection_from_structure_pick(
                sequence_hash=current_hash,
                structure_id=result.get("structure_id"),
                structure_hash=result.get("content_hash"),
                chain_id=parsed.get("chain_id"),
                model=parsed.get("model"),
                query_index=parsed.get("query_index_0based"),
                atom_name=parsed.get("atom_name"),
                label_seq_id=parsed.get("label_seq_id"),
                auth_seq_id=parsed.get("auth_seq_id"),
            )
            _store_molecule_selection(pick)
            query_index = parsed.get("query_index_0based")
            if query_index is None:
                st.info("No mapped structural residue")
            elif int(query_index) != int(selected):
                st.session_state["prot_3d_pending_pos"] = int(query_index)
                st.rerun()
    except Exception as exc:
        if "insufficient" in str(exc).lower():
            st.markdown(status_badge("INSUFFICIENT_DATA"), unsafe_allow_html=True)
            st.error("3D viewport has no finite coordinates to draw.")
        else:
            st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
            st.error("3D unavailable in this browser/environment.")
        st.caption("Sequence, mapping table and structure metadata remain available above.")
        return
    report = structure_scene.structure_3d_report(scene, envelope)
    st.download_button(
        "3D structure summary (JSON)",
        data=json.dumps(report, indent=2).encode("utf-8"),
        file_name="helixscope_structure_3d.json",
        mime="application/json",
        key="prot_3d_report",
    )


def _render_protein_structure_section(protein_seq: str) -> Optional[dict]:
    """Busca e carrega estruturas reais sob demanda. Nao dispara no Analyze.

    Args:
        protein_seq: Proteina ja validada nesta aba.

    Returns:
        Envelope carregado na sessao, ou None.
    """
    sources = protein_structure.remote_source_status()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(sources['rcsb_pdb']['method'])}. "
        f"{html_escape(sources['alphafold_db']['method'])}. "
        f"{html_escape(sources['uniprot']['method'])}. "
        "A protein sequence does not imply a structure. Predicted models stay "
        "predicted. Experimental structures stay experimental. 3D rendering uses "
        "Plotly Scatter3d on validated coordinates only after View 3D Structure. "
        "Opening this section does not fetch a structure."
        "</span>",
        unsafe_allow_html=True,
    )
    dssp = protein_structure.dssp_availability()
    st.caption(
        f"DSSP/STRIDE: {dssp.get('status')} "
        f"{html_escape(dssp.get('version') or 'version unavailable')}. "
        f"{html_escape(dssp.get('reason'))} HelixScope does not invent "
        "secondary-structure assignments or call a house heuristic STRIDE. "
        "DSSP is geometric assignment of a resolved structure, not sequence prediction."
    )
    pdb_id = st.text_input("PDB ID (optional)", key="prot_struct_pdb", placeholder="1CRN")
    uniprot = st.text_input(
        "UniProt accession (optional)", key="prot_struct_uniprot", placeholder="P01308"
    )
    alphafold_id = st.text_input(
        "AlphaFold identifier or UniProt (optional)",
        key="prot_struct_af_id",
        placeholder="AF-P01308-F1",
    )
    seq_search = st.checkbox(
        "Search RCSB PDB by sequence (MMseqs2 Search API, experimental entries)",
        key="prot_struct_seq_search",
    )
    uniprot_lookup = st.checkbox(
        "Look up UniProt PDB and AlphaFoldDB cross-references",
        key="prot_struct_uniprot_lookup",
    )
    af_lookup = st.checkbox(
        "Look up AlphaFold DB prediction metadata",
        key="prot_struct_af_lookup",
    )
    if st.button("Search structure sources", key="prot_struct_search_btn"):
        try:
            with st.spinner("Retrieving protein structure records..."):
                found = protein_structure.search_structure_sources(
                    protein_seq,
                    pdb_id=pdb_id,
                    uniprot=uniprot,
                    alphafold_id=alphafold_id,
                    search_pdb_by_sequence=seq_search,
                    lookup_uniprot=uniprot_lookup,
                    lookup_alphafold=af_lookup,
                )
        except (protein_structure.StructureError, ValueError) as exc:
            _render_structure_error(exc)
        else:
            st.session_state["prot_struct_hits"] = found
            st.session_state.pop("prot_struct_result", None)

    found = st.session_state.get("prot_struct_hits")
    if isinstance(found, dict) and str(found.get("sequence_hash") or "") != provenance.sequence_digest(
        protein_seq
    ):
        found = None
    if isinstance(found, dict):
        st.markdown(status_badge(str(found.get("status") or "UNAVAILABLE")), unsafe_allow_html=True)
        st.caption(str(found.get("message") or ""))
        st.caption(str(found.get("disclaimer") or ""))
        for item in list(found.get("errors") or []):
            st.caption(
                f"{item.get('source')}: {item.get('category')}: {item.get('message')}"
            )
        hits = list(found.get("hits") or [])
        if hits:
            rows = [
                {
                    "Source": item.get("source"),
                    "Kind": item.get("kind"),
                    "ID": item.get("structure_id"),
                    "Method": item.get("method") or "N/A",
                    "Resolution_A": item.get("resolution_angstrom")
                    if item.get("resolution_angstrom") is not None
                    else "N/A",
                    "RCSB_score": item.get("rcsb_score")
                    if item.get("rcsb_score") is not None
                    else "N/A",
                    "Match": item.get("match_basis"),
                }
                for item in hits
            ]
            st.dataframe(_export_frame(rows), width="stretch", hide_index=True)
            labels = [
                f"{item.get('source')} | {item.get('kind')} | {item.get('structure_id')}"
                for item in hits
            ]
            pick = st.selectbox(
                "Select a retrieved record (HelixScope does not rank a best structure)",
                list(range(len(labels))),
                format_func=lambda i: labels[i],
                key="prot_struct_pick",
            )
            selected = hits[int(pick)]
            chain_id = st.text_input(
                "Chain ID (optional; empty uses the first polymer chain)",
                key="prot_struct_chain",
            )
            model_text = st.text_input(
                "NMR model number (optional; empty uses the first declared model)",
                key="prot_struct_model",
            )
            if st.button("Load selected structure", key="prot_struct_load"):
                model_number = None
                if str(model_text or "").strip():
                    try:
                        model_number = int(str(model_text).strip())
                    except ValueError:
                        _render_structure_error(
                            protein_structure.StructureError(
                                "NMR model number must be an integer.",
                                "INVALID_INPUT",
                            )
                        )
                        model_number = False
                if model_number is not False:
                    try:
                        with st.spinner("Retrieving and validating structure coordinates..."):
                            structure_text = ""
                            if selected.get("bundled_fixture"):
                                reader = getattr(protein_structure, "read_bundled_mmcif", None)
                                if callable(reader):
                                    structure_text = reader(
                                        str(selected.get("bundled_fixture"))
                                    )
                            loaded = protein_structure.load_protein_structure(
                                protein_seq,
                                source=str(selected.get("source") or ""),
                                structure_id=str(selected.get("structure_id") or ""),
                                chain_id=str(chain_id or "").strip(),
                                model_number=model_number,
                                cif_url=str(selected.get("cif_url") or ""),
                                metadata=selected,
                                structure_text=structure_text,
                                cache=st.session_state.setdefault("prot_struct_cache", {}),
                            )
                    except (protein_structure.StructureError, ValueError) as exc:
                        _render_structure_error(exc)
                    else:
                        st.session_state["prot_struct_result"] = loaded
                        st.session_state.pop("prot_view_3d", None)
                        st.session_state.pop("prot_3d_cam", None)
                        st.session_state.pop("prot_3d_visible", None)
                        st.session_state.pop("prot_3d_chain", None)
                        st.session_state.pop("prot_3d_model", None)
                        st.session_state.pop("molecule_selection", None)
                        st.session_state.pop("molecule_selection_applied", None)

    result = st.session_state.get("prot_struct_result")
    if not isinstance(result, dict):
        return None
    if str(result.get("sequence_hash")) != provenance.sequence_digest(protein_seq):
        return None
    st.markdown(status_badge("RETRIEVED"), unsafe_allow_html=True)
    st.markdown(status_badge(str(result.get("kind") or "").upper()), unsafe_allow_html=True)
    if str(result.get("cache_status") or "") == "cached":
        st.markdown(status_badge("CACHED"), unsafe_allow_html=True)
    resolution = result.get("resolution_angstrom")
    if resolution is None or (
        isinstance(resolution, float) and (math.isnan(resolution) or math.isinf(resolution))
    ):
        resolution_text = "N/A"
    else:
        try:
            resolution_text = f"{float(resolution):.2f}"
        except (TypeError, ValueError):
            resolution_text = "N/A"
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"Source: {html_escape(result.get('source'))} | Kind: {html_escape(result.get('kind'))} | "
        f"ID: {html_escape(result.get('structure_id'))} | Method: {html_escape(result.get('method') or 'N/A')} | "
        f"Resolution: {html_escape(resolution_text)} A | Chain: "
        f"{html_escape((result.get('selected_chain') or {}).get('chain_id'))} | Model: "
        f"{html_escape(result.get('model_number'))}"
        "</span>",
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{html_escape(result.get('disclaimer'))}</span>",
        unsafe_allow_html=True,
    )
    st.caption(str(result.get("model_note") or ""))
    helix_workspace.render_explanation("structure", result)
    if result.get("sequences_identical") is False:
        st.warning(
            "The analysis sequence is not identical to the structure polymer sequence. "
            "Mapping used Needleman-Wunsch (BLOSUM62). This is not perfect identity."
        )
    coverage = result.get("coverage_query_with_coordinates")
    if coverage is None or (isinstance(coverage, float) and math.isnan(coverage)):
        coverage_text = "N/A"
    else:
        coverage_text = f"{float(coverage) * 100:.1f}"
    identity = result.get("identity_ungapped_pct")
    if identity is None or (isinstance(identity, float) and math.isnan(identity)):
        identity_text = "N/A"
    elif scientific_checks.percent_is_valid(identity):
        identity_text = f"{float(identity):.2f}"
    else:
        identity_text = "N/A"
    confidence = result.get("confidence_data")
    coverage_items = [
        (
            "Coordinate coverage",
            coverage_text,
            "% of analysis residues with coordinates" if coverage_text != "N/A" else "",
        ),
        (
            "Ungapped identity (mapping alignment)",
            identity_text,
            "%" if identity_text != "N/A" else "",
        ),
    ]
    if isinstance(confidence, dict) and confidence.get("mean") is not None:
        coverage_items.append(
            (
                "Mean B-factor pLDDT (predicted models)",
                f"{float(confidence['mean']):.2f}",
                str(confidence.get("range") or ""),
            )
        )
    render_metric_grid(coverage_items)
    st.caption(
        str(result.get("coverage_formula") or "")
        + ". Missing coordinates mean the residue exists in the sequence but has no atoms."
    )
    if isinstance(confidence, dict) and confidence.get("metric"):
        st.caption(str(confidence.get("metric") or ""))
    mapping_status = str(result.get("mapping_status") or "UNCERTAIN")
    st.markdown(status_badge(mapping_status), unsafe_allow_html=True)
    st.caption(str(result.get("mapping_status_note") or ""))
    if mapping_status in {"PARTIAL", "BEST_EFFORT"}:
        st.warning(
            "This mapping is not exact sequence mapping. "
            + str(result.get("mapping_status_note") or "")
        )
    pending = st.session_state.pop("prot_3d_pending_pos", None)
    if isinstance(pending, int) and 0 <= pending < len(protein_seq):
        st.session_state["prot_struct_pos"] = pending
    _apply_pending_molecule_selection(protein_seq)
    nav_cols = st.columns([3, 1, 1])
    with nav_cols[0]:
        selected = int(
            st.number_input(
                "Inspect analysis residue (0-based)",
                min_value=0,
                max_value=max(0, len(protein_seq) - 1),
                value=0,
                step=1,
                key="prot_struct_pos",
            )
        )
    with nav_cols[1]:
        if st.button("Previous residue", key="prot_res_prev") and selected > 0:
            st.session_state["prot_struct_pos"] = selected - 1
            st.rerun()
    with nav_cols[2]:
        if st.button("Next residue", key="prot_res_next") and selected < len(protein_seq) - 1:
            st.session_state["prot_struct_pos"] = selected + 1
            st.rerun()
    selected, highlight_indices, selection = _protein_highlight_indices(
        protein_seq, selected, result
    )
    _render_protein_structure_3d(
        protein_seq,
        result,
        selected,
        highlight_indices=highlight_indices,
        selection=selection,
    )
    _sequence_viewer_panel(
        protein_seq,
        "PROT",
        key_prefix="prot_view",
        highlight_indices=highlight_indices or None,
    )
    map_rows = [
        {
            "Query_0based": item.get("query_index_0based"),
            "AA": item.get("query_residue"),
            "Chain": item.get("chain_id"),
            "label_seq_id": _na_text(item.get("label_seq_id")),
            "auth_seq_id": _na_text(item.get("auth_seq_id")),
            "Insertion": item.get("insertion_code"),
            "Has_coords": item.get("has_coordinates"),
            "Mapping": item.get("mapping_row_status"),
        }
        for item in list(result.get("residue_mapping") or [])
    ]
    st.dataframe(_export_frame(map_rows), width="stretch", hide_index=True)
    st.caption(str(result.get("alt_location_policy") or ""))
    chains = list(result.get("chains") or [])
    if chains:
        st.caption(
            "Chains are not concatenated. "
            + "; ".join(
                f"{item.get('chain_id')} model {item.get('model')} "
                f"({item.get('length')} aa)"
                for item in chains
            )
        )
    deposited_ss = result.get("deposited_secondary_structure") or {}
    if deposited_ss.get("helices") or deposited_ss.get("sheets"):
        st.caption(
            "Deposited secondary-structure annotation from mmCIF "
            f"{deposited_ss.get('source')} — assignment in the coordinate file, "
            "not a sequence prediction and not DSSP/STRIDE."
        )
        st.caption(
            f"Helix records: {deposited_ss.get('n_helix_records')}; "
            f"sheet records: {deposited_ss.get('n_sheet_records')}."
        )
    if st.button("Compute CA-CA contact map (8 A)", key="prot_contact_map"):
        try:
            contacts = protein_structure.ca_contact_map(result, threshold_angstrom=8.0)
        except protein_structure.StructureError as exc:
            _render_structure_error(exc)
        else:
            contacts["structure_hash"] = result.get("content_hash")
            st.session_state["prot_contacts"] = contacts
    contacts = st.session_state.get("prot_contacts")
    if isinstance(contacts, dict):
        if str(contacts.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_contacts", None)
        else:
            st.caption(
                f"{contacts.get('method')} | threshold "
                f"{contacts.get('threshold_angstrom')} A | "
                f"{contacts.get('n_pairs')} pairs among "
                f"{contacts.get('n_residues_with_ca')} CA atoms. "
                f"{contacts.get('disclaimer')}"
            )
            helix_workspace.render_explanation(
                "contacts",
                {
                    "status": contacts.get("status") or "COMPUTED",
                    "n_contacts": contacts.get("n_pairs"),
                    "method": contacts.get("method") or "CA-CA distance",
                    "source": "HelixScope structure analysis",
                },
            )
            pair_rows = [
                {
                    "i": item.get("i"),
                    "j": item.get("j"),
                    "distance_angstrom": item.get("distance_angstrom"),
                }
                for item in list(contacts.get("pairs") or [])[:80]
            ]
            if pair_rows:
                st.dataframe(_export_frame(pair_rows), width="stretch", hide_index=True)
                if int(contacts.get("n_pairs") or 0) > 80:
                    st.caption("Table truncated to 80 pairs. The count above is complete.")
    map_dssp = getattr(protein_structure, "map_dssp_assignments", None)
    assign_stride = getattr(protein_structure, "assign_secondary_structure_stride", None)
    interchain_fn = getattr(protein_structure, "interchain_contacts", None)
    sasa_fn = getattr(protein_structure, "shrake_rupley_sasa", None)
    raw_text = str(result.get("raw_structure_text") or "")
    if dssp.get("available") and st.button("Assign DSSP (local mkdssp)", key="prot_dssp_run"):
        dssp_result = protein_structure.assign_secondary_structure_dssp(raw_text)
        if map_dssp and dssp_result.get("status") == "EXPERIMENTAL":
            chain = result.get("selected_chain") or {}
            dssp_result["mapping"] = map_dssp(
                list(dssp_result.get("residues") or []),
                list(chain.get("residues") or []),
            )
        dssp_result["structure_hash"] = result.get("content_hash")
        st.session_state["prot_dssp"] = dssp_result
    if st.button(
        "Assign DSSP via PDB-REDO API (remote)",
        key="prot_dssp_remote_run",
        help=(
            "Sends these coordinates to the official PDB-REDO DSSP HTTP API. "
            "Not a local mkdssp binary. Not sequence prediction. "
            "Not LIVE_VALIDATED local."
        ),
    ):
        dssp_result = protein_structure.assign_secondary_structure_dssp_remote(
            raw_text,
            coordinate_kind=str(result.get("kind") or "experimental"),
        )
        mapped_ok = str(dssp_result.get("status") or "") in {"EXPERIMENTAL", "COMPUTED"}
        if map_dssp and mapped_ok:
            chain = result.get("selected_chain") or {}
            dssp_result["mapping"] = map_dssp(
                list(dssp_result.get("residues") or []),
                list(chain.get("residues") or []),
            )
        dssp_result["structure_hash"] = result.get("content_hash")
        st.session_state["prot_dssp"] = dssp_result
    dssp_state = st.session_state.get("prot_dssp")
    if isinstance(dssp_state, dict):
        if str(dssp_state.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_dssp", None)
            st.markdown(status_badge("STALE"), unsafe_allow_html=True)
            st.caption("DSSP assignment belonged to a previous structure.")
        else:
            st.markdown(
                status_badge(str(dssp_state.get("status") or "UNAVAILABLE")),
                unsafe_allow_html=True,
            )
            st.caption(
                f"DSSP {dssp_state.get('tool') or ''} "
                f"{html_escape(str(dssp_state.get('version') or ''))}. "
                f"{html_escape(str(dssp_state.get('reason') or dssp_state.get('note') or ''))} "
                "DSSP is geometric assignment, not sequence prediction, and not a cartoon invented from consecutive residues."
            )
            mapping = dssp_state.get("mapping") or {}
            mapped_rows = list(mapping.get("rows") or [])[:80]
            if mapped_rows:
                st.dataframe(pd.DataFrame(mapped_rows), width="stretch", hide_index=True)
    stride_available = protein_structure.stride_availability().get("available")
    if assign_stride and stride_available and st.button("Assign STRIDE", key="prot_stride_run"):
        stride_result = assign_stride(raw_text)
        stride_result["structure_hash"] = result.get("content_hash")
        st.session_state["prot_stride"] = stride_result
    stride_state = st.session_state.get("prot_stride")
    if isinstance(stride_state, dict):
        if str(stride_state.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_stride", None)
        else:
            st.markdown(
                status_badge(str(stride_state.get("status") or "UNAVAILABLE")),
                unsafe_allow_html=True,
            )
            st.caption(
                "STRIDE is separate from DSSP. Results are not mixed. "
                f"{html_escape(str(stride_state.get('reason') or stride_state.get('note') or ''))}"
            )
    edtsurf_ready = protein_structure.edtsurf_availability().get("available")
    if edtsurf_ready and st.button(
        "Generate EDTSurf molecular surface (MS, probe 1.4 A)",
        key="prot_edtsurf_run",
        help=(
            "Runs the official EDTSurf executable on these coordinates. "
            "This is a triangulated molecular surface, not Shrake-Rupley SASA."
        ),
    ):
        kind_key = str(result.get("kind") or "experimental").strip().lower()
        if kind_key == "illustrative":
            st.caption("Molecular surface is not generated from ILLUSTRATIVE geometry.")
        else:
            surf = protein_structure.generate_molecular_surface(
                kind=kind_key if kind_key in {"experimental", "predicted"} else "experimental",
                structure_hash=str(result.get("content_hash") or ""),
                parsed={"atoms": list(result.get("atoms") or [])},
                surface_code=3,
                probe_radius=1.4,
            )
            surf["structure_hash"] = result.get("content_hash")
            st.session_state["prot_edtsurf"] = surf
    surf_state = st.session_state.get("prot_edtsurf")
    if isinstance(surf_state, dict):
        if str(surf_state.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_edtsurf", None)
        else:
            st.markdown(
                status_badge(str(surf_state.get("status") or "UNAVAILABLE")),
                unsafe_allow_html=True,
            )
            st.caption(
                f"{html_escape(str(surf_state.get('backend') or surf_state.get('tool') or 'surface'))} "
                f"{html_escape(str(surf_state.get('surface_definition') or ''))} "
                f"probe {html_escape(str(surf_state.get('probe_radius_angstrom') or ''))} A. "
                f"{html_escape(str(surf_state.get('reason') or ''))} "
                "This mesh is not Shrake-Rupley SASA."
            )
            st.caption(
                f"vertices {int(surf_state.get('n_vertices') or 0)} · "
                f"faces {int(surf_state.get('n_faces') or 0)}"
            )
    if interchain_fn and st.button("Inter-chain CA contacts (8 A)", key="prot_interchain"):
        parsed_atoms = result if result.get("atoms") else {
            "atoms": list((result.get("selected_chain") or {}).get("atoms") or [])
        }
        try:
            inter = interchain_fn(parsed_atoms, threshold_angstrom=8.0, atom_name="CA")
        except protein_structure.StructureError as exc:
            _render_structure_error(exc)
        else:
            inter["structure_hash"] = result.get("content_hash")
            st.session_state["prot_interchain"] = inter
    inter_state = st.session_state.get("prot_interchain")
    if isinstance(inter_state, dict):
        if str(inter_state.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_interchain", None)
        else:
            st.caption(
                f"{html_escape(str(inter_state.get('method')))} | "
                f"{inter_state.get('n_pairs')} pairs. "
                f"{html_escape(str(inter_state.get('disclaimer') or ''))}"
            )
    if sasa_fn and st.button("Shrake-Rupley SASA (not MSMS surface)", key="prot_sasa"):
        parsed_atoms = result if result.get("atoms") else {
            "atoms": list((result.get("selected_chain") or {}).get("atoms") or [])
        }
        try:
            sasa = sasa_fn(parsed_atoms)
        except protein_structure.StructureError as exc:
            _render_structure_error(exc)
        else:
            sasa["structure_hash"] = result.get("content_hash")
            st.session_state["prot_sasa"] = sasa
    sasa_state = st.session_state.get("prot_sasa")
    if isinstance(sasa_state, dict):
        if str(sasa_state.get("structure_hash") or "") != str(result.get("content_hash") or ""):
            st.session_state.pop("prot_sasa", None)
        else:
            st.markdown(
                status_badge(str(sasa_state.get("status") or "UNAVAILABLE")),
                unsafe_allow_html=True,
            )
            total = sasa_state.get("sasa_angstrom2")
            st.caption(
                "Shrake-Rupley SASA total: "
                f"{helix_workspace.unavailable_text(total)} A^2. "
                f"Probe {sasa_state.get('probe_radius_angstrom')} A. "
                "Solvent-accessible area, not a molecular surface mesh."
            )
            helix_workspace.render_explanation("sasa", sasa_state)
    export_cols = st.columns(3)
    with export_cols[0]:
        st.download_button(
            "Residue mapping (CSV)",
            data=protein_structure.export_mapping_csv(result).encode("utf-8"),
            file_name="helixscope_structure_mapping.csv",
            mime="text/csv",
            key="prot_struct_csv",
        )
    with export_cols[1]:
        st.download_button(
            "Structure JSON",
            data=json.dumps(protein_structure.export_result_bundle(result), indent=2).encode("utf-8"),
            file_name="helixscope_structure.json",
            mime="application/json",
            key="prot_struct_json",
        )
    with export_cols[2]:
        raw = str(result.get("raw_structure_text") or "")
        if raw:
            st.download_button(
                "Source mmCIF/PDB",
                data=raw.encode("utf-8"),
                file_name=ui_components.safe_download_filename(
                    f"{result.get('structure_id') or 'structure'}.cif"
                ),
                mime="chemical/x-cif",
                key="prot_struct_cif",
            )
    return result


def render_protein_analysis() -> None:
    """Renderiza a aba de analise de proteina.

    Aceita apenas sequencias detectadas como proteina. Ao clicar em Analyze
    exibe metricas, um badge de estabilidade, a composicao de aminoacidos, o
    perfil de hidrofobicidade e o visualizador da sequencia.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    in_col, side_col = st.columns([2, 1])
    with in_col:
        text_value = st.text_area("Paste a protein sequence", height=200, key="prot_text")
    with side_col:
        uploaded = st.file_uploader(
            "Or upload a FASTA file", type=["fasta", "fa"], key="prot_file"
        )
        st.caption(
            "A rejected upload is not analyzed. Use a simple file name such as sample.fasta."
        )
        analyze = st.button("Analyze", key="prot_analyze")

    if analyze:
        try:
            resolved = _resolve_sequence(text_value, uploaded)
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state["protein_input"] = resolved
        st.session_state.pop("prot_struct_hits", None)
        st.session_state.pop("prot_struct_result", None)
        st.session_state.pop("prot_view_3d", None)
        st.session_state.pop("prot_3d_cam", None)
        st.session_state.pop("prot_3d_visible", None)
        st.session_state.pop("prot_3d_chain", None)
        st.session_state.pop("prot_3d_model", None)
        st.session_state.pop("molecule_selection", None)
        st.session_state.pop("molecule_selection_applied", None)
        st.session_state.pop("prot_struct_pos", None)
        if not resolved.strip():
            st.error("Provide a protein sequence or a FASTA file to analyze.")
            return

    seq = st.session_state.get("protein_input", "")
    if not seq:
        return
    _warn_if_analyzed_input_changed(text_value, seq, action="Analyze")

    info = _validate_for(seq, "PROTEIN")
    st.markdown(_type_badge(info), unsafe_allow_html=True)

    if not info["is_valid"]:
        st.error(
            info["rejection_reason"]
            or "Provide a valid protein sequence (20 standard amino acids only)."
        )
        return

    clean = info["sequence"]
    _store_workspace("PROTEIN", clean, source="user sequence", identifier="")
    st.markdown(status_badge("COMPUTED"), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.result_header_html(
            "Protein analysis",
            status="COMPUTED",
            lines=[f"{info['length']:,} aa"],
        ),
        unsafe_allow_html=True,
    )
    try:
        phys = _physicochemical_report(clean, 7.0)
        mw = float(phys["molecular_weight_kda"])
        pi = float(phys["isoelectric_point"])
        instability = float(phys["instability_index"])
        gravy = float(phys["gravy"])
        aliphatic = float(phys["aliphatic_index"])
        charge = float(phys["net_charge"])
        aroma = float(phys["aromaticity"])
        categories = _amino_acid_categories(clean)
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))
        return

    helix_workspace.render_explanation(
        "protein",
        {
            "status": "COMPUTED",
            "length": len(clean),
            "molecular_weight_kda": mw,
            "isoelectric_point": pi,
            "gravy": gravy,
            "method": "protein physicochemical indices",
            "source": "HelixScope protein analysis",
        },
    )
    render_metric_grid(
        [
            ("Length", str(len(clean)), "aa"),
            ("MW", f"{mw:.3f}", "kDa"),
            ("Theoretical pI", f"{pi:.2f}", ""),
            ("Instability Index", f"{instability:.2f}", ""),
            ("GRAVY (Kyte-Doolittle)", f"{gravy:.3f}", ""),
            ("Net Charge (pH 7.0)", f"{charge:+.2f}", ""),
            ("Aliphatic Index", f"{aliphatic:.2f}", ""),
            (
                "Charged Residues",
                f"{categories['charged_total']['frequency']:.2f}",
                "%",
            ),
            ("Aromaticity (ProtParam)", f"{aroma:.4f}", ""),
            (
                "Extinction 280 reduced Cys",
                str(phys["extinction_coefficient_280_reduced"]),
                "M^-1 cm^-1",
            ),
            (
                "Extinction 280 cystine",
                str(phys["extinction_coefficient_280_cystine"]),
                "M^-1 cm^-1",
            ),
            (
                "Shannon complexity",
                f"{float(phys['sequence_complexity_bits']):.4f}",
                "bits",
            ),
        ]
    )
    st.caption(
        f"{phys['extinction_method']} {phys['sequence_complexity_method']} "
        "SEG low-complexity is UNAVAILABLE."
    )

    if instability < 40:
        st.markdown(
            badge("Computational instability index < 40 (Guruprasad)", "dna"),
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            badge("Computational instability index >= 40 (Guruprasad)", "warn"),
            unsafe_allow_html=True,
        )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(str(phys['isoelectric_point_method']))}. "
        f"{html_escape(str(phys['gravy_method']))} hydropathy. "
        f"Net charge at pH {phys['net_charge_ph']}: "
        f"{html_escape(str(phys['net_charge_method']))}. "
        f"{html_escape(str(phys['instability_index_method']))}.</span>",
        unsafe_allow_html=True,
    )

    composition = _amino_acid_composition(clean)

    if _section_toggle("Amino Acid Composition", key="prot_aa_comp"):
        with st.container(border=True):
            st.plotly_chart(
                charts.amino_acid_bar_chart(composition), width="stretch"
            )
            st.plotly_chart(charts.amino_acid_radar(composition), width="stretch")
            comp_rows = [
                {"Amino acid": symbol, "Count": data["count"], "Frequency_pct": data["frequency"]}
                for symbol, data in composition.items()
                if data["count"] > 0
            ]
            st.dataframe(pd.DataFrame(comp_rows), width="stretch", hide_index=True)
            if not protein_analysis.composition_is_consistent(clean):
                st.warning(
                    "Amino-acid counts do not sum to protein length. "
                    "The result is marked inconsistent and was not adjusted."
                )

    if _section_toggle("Residue Categories", key="prot_categories"):
        with st.container(border=True):
            st.plotly_chart(
                charts.amino_acid_category_chart(categories), width="stretch"
            )
            category_rows = [
                {
                    "Category": name.replace("_", " "),
                    "Residues": data["residues"],
                    "Count": data["count"],
                    "Frequency_pct": data["frequency"],
                }
                for name, data in categories.items()
            ]
            st.dataframe(
                pd.DataFrame(category_rows), width="stretch", hide_index=True
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>The first five "
                "categories are mutually exclusive and sum to 100%. Charged total "
                "and hydrophobic overlap with them and must not be added to the "
                "rest.</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("Hydrophobicity Profile (Kyte-Doolittle)", key="prot_hydro"):
        with st.container(border=True):
            window = 9 if len(clean) >= 9 else len(clean)
            profile = _hydrophobicity_profile(clean, window)
            st.plotly_chart(
                charts.hydrophobicity_plot(profile, clean), width="stretch"
            )
            records = protein_analysis.hydrophobicity_profile_records(clean, window)
            st.dataframe(_export_frame(records), width="stretch", hide_index=True)
            st.markdown(
                "<span style='color:var(--hs-text-secondary);'>Method: Kyte-Doolittle windowed "
                f"hydropathy (window={window}). GRAVY "
                f"{gravy:.3f} is the chain mean on the same scale. This is "
                "hydropathy, not TMHMM, SignalP or a domain assignment.</span>",
                unsafe_allow_html=True,
            )

    if _section_toggle("Formal charge profile", key="prot_charge"):
        with st.container(border=True):
            charge_window = 9 if len(clean) >= 9 else len(clean)
            charge_profile = _charge_profile(clean, charge_window)
            st.plotly_chart(
                charts.charge_profile_plot(charge_profile, clean), width="stretch"
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Window mean of "
                "(K+R minus D+E) / window. Histidine is omitted because it is "
                "only partially protonated at pH 7. Full-chain net charge at pH 7 "
                "uses ProtParam Henderson-Hasselbalch. This is not a transmembrane "
                "or domain predictor.</span>",
                unsafe_allow_html=True,
            )

    structure_result = None
    if _section_toggle("Protein Structure (PDB / AlphaFold / UniProt)", key="prot_structure"):
        with st.container(border=True):
            structure_result = _render_protein_structure_section(clean)
    if not isinstance(structure_result, dict):
        stored = st.session_state.get("prot_struct_result")
        if (
            isinstance(stored, dict)
            and str(stored.get("sequence_hash") or "") == provenance.sequence_digest(clean)
        ):
            structure_result = stored

    if _section_toggle("Secondary Structure", key="prot_secondary"):
        with st.container(border=True):
            availability = protein_analysis.secondary_structure_availability()
            st.error(
                "Secondary structure prediction is not available in HelixScope."
            )
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);'>{availability['reason']}</span>",
                unsafe_allow_html=True,
            )
            st.markdown(
                "<span style='color:var(--hs-text-secondary);'>Use one of these instead:</span>",
                unsafe_allow_html=True,
            )
            for tool in availability["recommended_tools"]:
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);'>- {html_escape(tool)}</span>",
                    unsafe_allow_html=True,
                )
            features = protein_analysis.feature_prediction_availability()
            st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
            for name, info_feat in features.items():
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
                    f"{html_escape(name)}: {html_escape(info_feat['reason'])}</span>",
                    unsafe_allow_html=True,
                )

    if _section_toggle("Sequence Viewer", key="prot_seq_view"):
        with st.container(border=True):
            selected_struct = st.session_state.get("prot_struct_pos")
            selection = _active_molecule_selection(clean, structure_result)
            highlight = None
            if selection and str(selection.get("status") or "") == "READY":
                from modules import molecule_selection as _mol_sel

                highlight = _mol_sel.query_indices_of(selection)
            elif isinstance(selected_struct, int):
                highlight = [int(selected_struct)]
            st.markdown(
                sequence_display(clean, "PROT", highlight_indices=highlight),
                unsafe_allow_html=True,
            )

    summary_rows = [
        {"Metric": "Length", "Value": len(clean), "Unit": "aa"},
        {"Metric": "Molecular weight", "Value": mw, "Unit": "kDa"},
        {"Metric": "Theoretical pI", "Value": pi, "Unit": "pH"},
        {"Metric": "Net charge at pH 7.0", "Value": charge, "Unit": "charge"},
        {"Metric": "GRAVY (Kyte-Doolittle)", "Value": gravy, "Unit": "hydropathy"},
        {"Metric": "Aliphatic index (Ikai 1980)", "Value": aliphatic, "Unit": ""},
        {"Metric": "Instability index (Guruprasad, computational)", "Value": instability, "Unit": ""},
        {"Metric": "Aromaticity (ProtParam)", "Value": aroma, "Unit": "fraction"},
        {
            "Metric": "Charged residues",
            "Value": categories["charged_total"]["frequency"],
            "Unit": "%",
        },
        {
            "Metric": "Hydrophobic residues (Kyte-Doolittle)",
            "Value": categories["hydrophobic_kyte_doolittle"]["frequency"],
            "Unit": "%",
        },
        {"Metric": "Secondary structure", "Value": "UNAVAILABLE", "Unit": ""},
        {"Metric": "Domains / signal / TM / conservation", "Value": "UNAVAILABLE", "Unit": ""},
        {
            "Metric": "Retrieved protein structure",
            "Value": (
                f"{structure_result.get('status')} {structure_result.get('kind')} "
                f"{structure_result.get('source')} {structure_result.get('structure_id')} "
                f"mapping={structure_result.get('mapping_status')}"
                if isinstance(structure_result, dict)
                else "Not retrieved (open Protein Structure)"
            ),
            "Unit": "",
        },
        {
            "Metric": "3D renderer",
            "Value": (
                "plotly Scatter3d"
                if isinstance(structure_result, dict) and st.session_state.get("prot_view_3d")
                else "Plotly Scatter3d (View 3D Structure, not auto-loaded)"
            ),
            "Unit": "",
        },
        {
            "Metric": "3D selection source",
            "Value": (
                str((_active_molecule_selection(clean, structure_result) or {}).get("source") or "sequence inspect residue")
            ),
            "Unit": "",
        },
    ]
    _render_summary_table(
        summary_rows, "prot_summary_download", "protein_summary.csv"
    )
    protein_report = provenance.analysis_envelope(
        module="PROTEIN",
        payload={
            **phys,
            "composition_consistent": protein_analysis.composition_is_consistent(clean),
            "secondary_structure": "UNAVAILABLE",
            "retrieved_structure": (
                {
                    "status": structure_result.get("status"),
                    "kind": structure_result.get("kind"),
                    "source": structure_result.get("source"),
                    "structure_id": structure_result.get("structure_id"),
                    "chain_id": (structure_result.get("selected_chain") or {}).get("chain_id"),
                    "model_number": structure_result.get("model_number"),
                    "coverage_query_with_coordinates": structure_result.get(
                        "coverage_query_with_coordinates"
                    ),
                    "mapping_status": structure_result.get("mapping_status"),
                }
                if isinstance(structure_result, dict)
                else None
            ),
            "renderer_3d": (
                "plotly"
                if isinstance(structure_result, dict) and st.session_state.get("prot_view_3d")
                else "data-only"
            ),
        },
        status="COMPUTED",
        algorithm="HelixScope protein physicochemical report",
        parameters={"ph": 7.0},
        sequence=clean,
        source="user sequence",
    )
    st.download_button(
        "Download analysis report (JSON)",
        data=json.dumps(protein_report, indent=2).encode("utf-8"),
        file_name="helixscope_protein_report.json",
        mime="application/json",
        key="prot_report_json",
    )
    _render_alerts(
        _protein_alerts(gravy, instability, pi, aliphatic, categories, composition)
    )


def _match_line(aligned_seq1: str, aligned_seq2: str) -> str:
    """Monta a linha de correspondencia de um alinhamento (uso interno).

    Args:
        aligned_seq1: Primeira sequencia alinhada (com gaps).
        aligned_seq2: Segunda sequencia alinhada (com gaps).

    Returns:
        String em que cada posicao recebe "|" para match, "." para mismatch e
        espaco para gap.

    Raises:
        Nenhum.
    """
    symbols = []
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == "-" or b == "-":
            symbols.append(" ")
        elif a == b:
            symbols.append("|")
        else:
            symbols.append(".")
    return "".join(symbols)


def render_alignment() -> None:
    """Renderiza a aba de alinhamento de sequencias.

    Recebe duas sequencias, executa alinhamento global ou local, exibe metricas
    e a visualizacao do alinhamento, e opcionalmente um dot plot.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    has_align = isinstance(st.session_state.get("align_result"), dict)
    with st.expander("Sequences and method", expanded=not has_align):
        cols = st.columns(2)
        with cols[0]:
            seq1 = st.text_area("Sequence 1", height=140, key="align_seq1")
        with cols[1]:
            seq2 = st.text_area("Sequence 2", height=140, key="align_seq2")

        algorithm = st.radio(
            "Algorithm",
            ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"],
            key="align_algorithm",
        )
        translate_nucleic = st.checkbox(
            "Translate nucleic acids (frame +1, stop at first stop) then align as protein",
            value=False,
            key="align_translate",
        )
        show_dotplot = st.checkbox("Show Dot Plot", value=False, key="align_dotplot")
        align_clicked = st.button("Align", key="align_button")
        msa_info = alignment.multiple_alignment_availability()
        blast_info = alignment.blast_search_availability()
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Pairwise only. "
            f"{html_escape(msa_info['reason'])} {html_escape(blast_info['reason'])}</span>",
            unsafe_allow_html=True,
        )

    if align_clicked:
        info1 = _validate(seq1)
        info2 = _validate(seq2)
        if not info1["is_valid"] or not info2["is_valid"]:
            st.session_state.pop("align_result", None)
            reason1 = info1.get("rejection_reason") or "sequence 1 is invalid"
            reason2 = info2.get("rejection_reason") or "sequence 2 is invalid"
            st.error(
                "Provide two valid sequences of the same molecule type "
                "(DNA, RNA or protein). "
                f"Sequence 1: {reason1} Sequence 2: {reason2}"
            )
            return
        if info1["type"] != info2["type"]:
            st.session_state.pop("align_result", None)
            st.error(
                f"Cannot align {info1['type']} with {info2['type']}. "
                "Both sequences must be the same molecule type."
            )
            return
        seq_a = info1["sequence"]
        seq_b = info2["sequence"]
        translation_note = ""
        if translate_nucleic:
            if info1["type"] not in {"DNA", "RNA"}:
                st.session_state.pop("align_result", None)
                st.error(
                    "Translation before alignment applies only to DNA or RNA pairs."
                )
                return
            try:
                if info1["type"] == "DNA":
                    seq_a = str(Seq(seq_a).translate(table=1, to_stop=True))
                    seq_b = str(Seq(seq_b).translate(table=1, to_stop=True))
                else:
                    seq_a = protein_analysis.translate(seq_a)
                    seq_b = protein_analysis.translate(seq_b)
            except (ValueError, RuntimeError) as exc:
                st.session_state.pop("align_result", None)
                st.error(str(exc))
                return
            if not seq_a or not seq_b:
                st.session_state.pop("align_result", None)
                st.error(
                    "Translation produced an empty protein. Check that each "
                    "sequence has a complete codon before the first stop."
                )
                return
            translation_note = (
                "Translated in frame +1 until the first stop codon. "
                "Not a predicted CDS or gene model."
            )
        mode = "local" if algorithm.startswith("Local") else "global"
        try:
            st.session_state["align_result"] = {
                "result": _pairwise(seq_a, seq_b, mode),
                "seq1": seq_a,
                "seq2": seq_b,
                "seq1_hash": provenance.sequence_digest(seq_a),
                "seq2_hash": provenance.sequence_digest(seq_b),
                "translate": bool(translate_nucleic),
                "translation_note": translation_note,
            }
        except ValueError as exc:
            st.session_state.pop("align_result", None)
            st.error(str(exc))
            return

    packed = st.session_state.get("align_result")
    if not packed:
        return
    current_1 = _validate(seq1)
    current_2 = _validate(seq2)
    live_a = str(current_1.get("sequence") or "")
    live_b = str(current_2.get("sequence") or "")
    if packed.get("translate") and current_1.get("is_valid") and current_2.get("is_valid"):
        try:
            if current_1.get("type") == "DNA":
                live_a = str(Seq(live_a).translate(table=1, to_stop=True))
                live_b = str(Seq(live_b).translate(table=1, to_stop=True))
            elif current_1.get("type") == "RNA":
                live_a = protein_analysis.translate(live_a)
                live_b = protein_analysis.translate(live_b)
        except (ValueError, RuntimeError):
            live_a, live_b = "", ""
    if not _sequence_hashes_match(
        packed.get("seq1_hash") or provenance.sequence_digest(str(packed.get("seq1") or "")),
        provenance.sequence_digest(live_a),
    ) or not _sequence_hashes_match(
        packed.get("seq2_hash") or provenance.sequence_digest(str(packed.get("seq2") or "")),
        provenance.sequence_digest(live_b),
    ):
        _hide_stale_result(
            "Stored pairwise alignment belongs to different sequences. Click Align to refresh."
        )
        return

    result = packed["result"]
    info1_seq = packed["seq1"]
    info2_seq = packed["seq2"]
    translation_note = str(packed.get("translation_note") or "")
    if result.get("aligned_seq1") or result.get("aligned_seq2"):
        if not scientific_checks.alignment_rows_are_consistent(
            str(result.get("aligned_seq1") or ""),
            str(result.get("aligned_seq2") or ""),
        ):
            st.error(
                "Analysis failed: aligned sequences have different lengths. "
                "The alignment was not displayed as a scientific result."
            )
            return

    st.markdown(status_badge(str(result.get("status") or "COMPUTED")), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.result_header_html(
            "Pairwise alignment",
            status=str(result.get("status") or "COMPUTED"),
            lines=[
                str(result.get("method") or result.get("algorithm") or ""),
                f"length {result.get('length')} · identity {result.get('identity_pct')}",
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation("alignment", result)
    ungapped = result.get("ungapped_identity_pct")
    ungapped_label = (
        "N/A"
        if ungapped is None or (isinstance(ungapped, float) and math.isnan(ungapped))
        else f"{float(ungapped):.2f}"
    )
    cov1 = result.get("coverage_seq1_pct")
    cov2 = result.get("coverage_seq2_pct")
    cov1_label = (
        "N/A"
        if cov1 is None or (isinstance(cov1, float) and math.isnan(cov1))
        else f"{float(cov1):.2f}"
    )
    cov2_label = (
        "N/A"
        if cov2 is None or (isinstance(cov2, float) and math.isnan(cov2))
        else f"{float(cov2):.2f}"
    )
    render_metric_grid(
        [
            ("Score", f"{result['score']:.2f}", ""),
            ("Identity (incl. gap columns)", f"{result['identity_pct']:.2f}", "%"),
            ("Gaps", str(result["gaps"]), ""),
            ("Aligned Length", str(result["length"]), ""),
            ("Matches", str(result.get("matches", "-")), ""),
            ("Mismatches", str(result.get("mismatches", "-")), ""),
            (
                "Similarity",
                f"{float(result.get('similarity_pct', result['identity_pct'])):.2f}",
                "%",
            ),
            (
                "Ungapped identity",
                ungapped_label,
                "%" if ungapped_label != "N/A" else "",
            ),
            (
                "Coverage seq1",
                cov1_label,
                "%" if cov1_label != "N/A" else "",
            ),
            (
                "Coverage seq2",
                cov2_label,
                "%" if cov2_label != "N/A" else "",
            ),
        ]
    )
    st.caption(
        "Pairwise identity includes gap columns. Ungapped identity uses matches / "
        "(matches + mismatches). Coverage is non-gap residues / original length. "
        "These are not BLAST HSP identity or MSA column identity."
    )
    scoring = html_escape(str(result.get("scoring") or "identity match=2 mismatch=-1"))
    molecule = html_escape(str(result.get("molecule") or "unknown"))
    method = html_escape(str(result.get("method") or "pairwise"))
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{method} via Biopython "
        f"PairwiseAligner. Molecule: {molecule}. Scoring: {scoring}. Gap open "
        "-5.0, gap extend -0.5. Identity is matches divided by aligned length "
        "including gap columns. Ungapped identity is matches / (matches + "
        "mismatches) and is a separate metric. Coverage seq1 = non-gap residues "
        "in row 1 / original length of sequence 1 (same for seq2). Similarity "
        "equals identity for DNA/RNA and BLOSUM positives for protein.</span>",
        unsafe_allow_html=True,
    )
    if translation_note:
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{html_escape(translation_note)}</span>",
            unsafe_allow_html=True,
        )

    if _section_toggle("Alignment View", key="align_view"):
        with st.container(border=True):
            match = _match_line(result["aligned_seq1"], result["aligned_seq2"])
            block = "\n".join(
                [result["aligned_seq1"], match, result["aligned_seq2"]]
            )
            st.markdown(
                f"<pre style=\"font-family:'JetBrains Mono',monospace;font-size:13px;"
                "color:var(--hs-text);background:rgba(255,255,255,0.02);"
                "border:1px solid rgba(255,255,255,0.06);padding:14px;border-radius:12px;"
                f"overflow:auto;\">{html_escape(block)}</pre>",
                unsafe_allow_html=True,
            )
            if result.get("aligned_seq1") and result.get("aligned_seq2"):
                try:
                    export_text = alignment.format_pairwise_export(
                        result, "seq1", "seq2"
                    )
                    st.download_button(
                        "Download alignment (text)",
                        data=export_text.encode("utf-8"),
                        file_name="helixscope_alignment.txt",
                        mime="text/plain",
                        key="align_export_txt",
                    )
                except ValueError as exc:
                    st.info(str(exc))
                try:
                    classes = alignment.classify_alignment_columns(
                        str(result["aligned_seq1"]),
                        str(result["aligned_seq2"]),
                    )
                    st.plotly_chart(
                        charts.alignment_column_map(classes),
                        width="stretch",
                    )
                except ValueError as exc:
                    st.info(str(exc))
    align_report = provenance.analysis_envelope(
        module="ALIGNMENT",
        payload={
            "score": result.get("score"),
            "identity_pct": result.get("identity_pct"),
            "ungapped_identity_pct": result.get("ungapped_identity_pct"),
            "coverage_seq1_pct": result.get("coverage_seq1_pct"),
            "coverage_seq2_pct": result.get("coverage_seq2_pct"),
            "matches": result.get("matches"),
            "mismatches": result.get("mismatches"),
            "gaps": result.get("gaps"),
            "alignment_length": result.get("length"),
            "method": result.get("method"),
            "scoring": result.get("scoring"),
            "molecule": result.get("molecule"),
            "status": result.get("status"),
        },
        status=str(result.get("status") or "COMPUTED"),
        algorithm=str(result.get("method") or "pairwise"),
        parameters=result.get("parameters") or {},
        sequence=str(info1_seq),
        source="user sequence",
    )
    st.download_button(
        "Download alignment report (JSON)",
        data=json.dumps(align_report, indent=2).encode("utf-8"),
        file_name="helixscope_alignment_report.json",
        mime="application/json",
        key="align_report_json",
    )

    if show_dotplot:
        if _section_toggle("Dot Plot", key="align_dotplot_view"):
            with st.container(border=True):
                if len(info1_seq) > alignment.MAX_DOTPLOT_LENGTH or len(
                    info2_seq
                ) > alignment.MAX_DOTPLOT_LENGTH:
                    st.warning(
                        f"Dot plot is limited to {alignment.MAX_DOTPLOT_LENGTH:,} "
                        "residues per sequence."
                    )
                else:
                    matrix = _dotplot_matrix(info1_seq, info2_seq)
                    st.plotly_chart(
                        charts.alignment_dotplot(matrix, "Sequence 1", "Sequence 2"),
                        width="stretch",
                    )


def _format_feature_qualifiers(qualifiers: object) -> str:
    """Concatena qualifiers GenBank com valores reais (uso interno)."""
    if not isinstance(qualifiers, dict) or not qualifiers:
        return "-"
    parts: list[str] = []
    for name, value in sorted(qualifiers.items()):
        if isinstance(value, list):
            text = "; ".join(str(item) for item in value)
        elif value is None:
            text = ""
        else:
            text = str(value)
        if len(text) > 500:
            text = text[:500] + "..."
        parts.append(f"{name}={text}")
    return " | ".join(parts) if parts else "-"


def _render_ncbi_error(exc: Exception) -> None:
    """Exibe um erro NCBI com mensagem especifica por categoria (interno)."""
    category = getattr(exc, "category", "internal")
    failure_kind = getattr(exc, "failure_kind", category)
    message = str(exc)
    if failure_kind == "timeout":
        st.error(message)
        st.info(
            "NCBI request timed out. This is not the same as zero records. "
            "No fabricated record was generated."
        )
        return
    if failure_kind == "rate_limited":
        st.error(message)
        st.info(
            "NCBI rate limited the request. Wait and retry. "
            "No fabricated record was generated."
        )
        return
    if category == "invalid_input":
        st.error(message)
        return
    if category == "not_found":
        st.error(message)
        return
    if category == "wrong_database":
        st.error(message)
        return
    if category == "source_unavailable":
        st.error(message)
        st.info("NCBI unavailable. No fabricated record was generated.")
        return
    if category == "configuration":
        st.error(message)
        return
    st.error(f"NCBI request failed: {message}")


def render_ncbi_fetch() -> None:
    """Renderiza a aba de recuperacao de registros no NCBI.

    Classifica o identificador, consulta o banco Entrez apropriado (nucleotide,
    protein ou gene) e exibe descricao, metricas, features, sequencia quando
    existir, registros relacionados e um botao de download em JSON.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    st.markdown(
        '<p style="color:var(--hs-text-secondary);font-size:13px;margin-bottom:16px;">'
        "NCBI Entrez requires a valid contact email for every request. "
        "Numeric IDs are treated as GeneIDs. Names of organisms are not "
        "accessions.</p>",
        unsafe_allow_html=True,
    )
    email = st.text_input(
        "NCBI email",
        key="ncbi_email",
        placeholder="your.email@university.edu",
    )
    blast_info = ncbi_fetch.homology_search_availability()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{html_escape(blast_info['reason'])} "
        "Open the BLAST tab to submit a real NCBI BLAST job. Pairwise alignment "
        "is not BLAST.</span>",
        unsafe_allow_html=True,
    )

    st.markdown("#### Search")
    search_col, search_db_col = st.columns([2, 1])
    with search_col:
        search_term = st.text_input(
            "NCBI Search query",
            placeholder="BRCA1 AND human[orgn]",
            key="ncbi_search_term",
        )
    with search_db_col:
        search_db = st.selectbox(
            "Search database",
            list(ncbi_fetch.SEARCHABLE_DATABASES),
            key="ncbi_search_db",
        )
    search_clicked = st.button("Search NCBI", key="ncbi_search_button")
    if search_clicked:
        if not email.strip():
            st.error("Enter your NCBI email above before searching.")
        elif not search_term.strip():
            st.error("Provide a search query.")
        elif _ncbi_allow_request():
            try:
                with st.spinner("Searching NCBI..."):
                    st.session_state["ncbi_search"] = _search_ncbi(
                        search_term.strip(), search_db, email.strip(), 0
                    )
            except ncbi_fetch.NCBIQueryError as exc:
                st.session_state.pop("ncbi_search", None)
                _render_ncbi_error(exc)
            except ValueError as exc:
                st.session_state.pop("ncbi_search", None)
                st.error(str(exc))
            except RuntimeError as exc:
                st.session_state.pop("ncbi_search", None)
                _render_ncbi_error(exc)

    search_payload = st.session_state.get("ncbi_search")
    if isinstance(search_payload, dict):
        st.markdown(status_badge("RETRIEVED"), unsafe_allow_html=True)
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Source: "
            f"{html_escape(search_payload.get('source'))} | "
            f"retrieved {html_escape(search_payload.get('retrieved_at_utc'))} | "
            f"NCBI count {html_escape(search_payload.get('count'))} | "
            f"offset {html_escape(search_payload.get('retstart'))} | "
            f"page size {html_escape(search_payload.get('retmax'))} | "
            f"local cache TTL 1800 s</span>",
            unsafe_allow_html=True,
        )
        hits = list(search_payload.get("records") or [])
        if search_payload.get("status") == "no_records_found" or not hits:
            st.info(
                "Search completed with zero records. This is not a request "
                "failure and not BLAST NO_HITS."
            )
        else:
            hit_rows = [
                {
                    "Accession": item.get("accession"),
                    "Title": item.get("title"),
                    "Organism": item.get("organism"),
                    "Length": item.get("length") if item.get("length") is not None else "Not available",
                    "Database": item.get("database"),
                    "Extra": item.get("extra"),
                }
                for item in hits
            ]
            st.dataframe(pd.DataFrame(hit_rows), width="stretch", hide_index=True)
            labels = [
                f"{item.get('accession')} | {str(item.get('title') or '')[:80]}"
                for item in hits
            ]
            chosen = st.selectbox("Selected search hit", labels, key="ncbi_hit_choice")
            if st.button("Use selected accession", key="ncbi_use_hit"):
                index = labels.index(chosen)
                st.session_state["ncbi_accession"] = str(hits[index].get("accession") or "")
                st.rerun()
            total_count = int(search_payload.get("count") or 0)
            retstart = int(search_payload.get("retstart") or 0)
            retmax = int(search_payload.get("retmax") or 15)
            page_cols = st.columns(2)
            query_text = str(search_payload.get("query") or search_term)
            query_db = str(search_payload.get("database") or search_db)
            with page_cols[0]:
                if retstart > 0 and st.button("Previous NCBI page", key="ncbi_prev_page"):
                    if _ncbi_allow_request():
                        new_start = max(0, retstart - retmax)
                        try:
                            st.session_state["ncbi_search"] = _search_ncbi(
                                query_text, query_db, email.strip(), new_start
                            )
                            st.rerun()
                        except ncbi_fetch.NCBIQueryError as exc:
                            _render_ncbi_error(exc)
            with page_cols[1]:
                if (retstart + len(hits) < total_count) and st.button(
                    "Next NCBI page", key="ncbi_next_page"
                ):
                    if _ncbi_allow_request():
                        try:
                            st.session_state["ncbi_search"] = _search_ncbi(
                                query_text, query_db, email.strip(), retstart + retmax
                            )
                            st.rerun()
                        except ncbi_fetch.NCBIQueryError as exc:
                            _render_ncbi_error(exc)

    st.markdown("#### Direct accession fetch")
    acc_col, db_col = st.columns([2, 1])
    with acc_col:
        accession = st.text_input(
            "Accession or GeneID", placeholder="NM_007294.4", key="ncbi_accession"
        )
    with db_col:
        db = st.selectbox("Database", ["nucleotide", "protein", "gene"], key="ncbi_db")

    fetch_clicked = st.button("Fetch", key="ncbi_fetch_button")

    if fetch_clicked:
        if not email.strip():
            st.error("Enter your NCBI email above before fetching records.")
            return
        if not accession.strip():
            st.error("Provide an accession number or a GeneID.")
            return
        if not _ncbi_allow_request():
            return
        identifier_info = ncbi_fetch.classify_identifier(accession.strip())
        try:
            with st.spinner("Fetching from NCBI..."):
                record = _fetch_by_accession(
                    accession.strip(), db, email.strip()
                )
        except ncbi_fetch.NCBIQueryError as exc:
            st.session_state.pop("ncbi_record", None)
            st.markdown(
                badge(
                    f"Detected: {identifier_info['kind']} -> "
                    f"{identifier_info['suggested_db'] or 'none'}",
                    "rna",
                ),
                unsafe_allow_html=True,
            )
            _render_ncbi_error(exc)
            return
        except ValueError as exc:
            st.session_state.pop("ncbi_record", None)
            st.error(str(exc))
            return
        except RuntimeError as exc:
            st.session_state.pop("ncbi_record", None)
            _render_ncbi_error(exc)
            return
        st.session_state["ncbi_record"] = record
        st.session_state["ncbi_identifier"] = identifier_info

    record = st.session_state.get("ncbi_record")
    if not isinstance(record, dict):
        if not email.strip():
            st.info("Enter your NCBI email above before fetching records.")
        else:
            st.info("Fetch a record to inspect features, sequence and downloads.")
        return

    identifier_info = st.session_state.get("ncbi_identifier") or {
        "kind": record.get("identifier_kind", "unknown"),
        "suggested_db": record.get("database"),
    }
    st.markdown(
        badge(
            f"Detected: {identifier_info.get('kind')} -> "
            f"{identifier_info.get('suggested_db') or 'none'}",
            "rna",
        ),
        unsafe_allow_html=True,
    )
    st.markdown(
        helix_workspace.result_header_html(
            str(record.get("accession") or record.get("identifier") or "NCBI record"),
            status="RETRIEVED",
            lines=[
                str(record.get("organism") or "Unknown organism"),
                str(record.get("database") or ""),
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation("ncbi", record)

    resolved_db = str(record.get("database") or db)
    record_kind = str(record.get("record_kind") or "sequence")
    sequence_available = bool(record.get("sequence_available", True))
    gene_meta = record.get("gene") if isinstance(record.get("gene"), dict) else None

    source_line = "Source: NCBI Entrez"
    if record.get("retrieved_at_utc"):
        source_line += f" | retrieved {html_escape(record['retrieved_at_utc'])}"
    source_line += " | local cache TTL 1800 s"
    if record.get("genbank_date"):
        source_line += f" | GenBank date {html_escape(record['genbank_date'])}"
    description = html_escape(record.get("description") or "")
    organism = html_escape(record.get("organism") or "Unknown")
    resolved_db_safe = html_escape(resolved_db)
    record_kind_safe = html_escape(record_kind)
    explorer = ncbi_fetch.record_explorer_fields(record)
    explorer_bits = []
    for key in (
        "accession",
        "version",
        "molecule",
        "topology",
        "length",
    ):
        if key in explorer:
            explorer_bits.append(
                f"{html_escape(key)}: {html_escape(explorer[key])}"
            )
    taxonomy = explorer.get("taxonomy") or []
    taxonomy_line = ""
    if taxonomy:
        taxonomy_line = html_escape("; ".join(str(item) for item in taxonomy[:12]))
    card_content = (
        f"<div style='color:var(--hs-text);margin-bottom:6px;'>{description}</div>"
        f"<div style='color:var(--hs-text-secondary);'>Organism: {organism}</div>"
        f"<div style='color:var(--hs-text-secondary);'>Database: {resolved_db_safe} "
        f"({record_kind_safe})</div>"
        f"<div style='color:var(--hs-text-secondary);'>{source_line}</div>"
    )
    if explorer_bits:
        card_content += (
            "<div style='color:var(--hs-text-secondary);'>"
            + " | ".join(explorer_bits)
            + "</div>"
        )
    if taxonomy_line:
        card_content += (
            f"<div style='color:var(--hs-text-secondary);'>Taxonomy: {taxonomy_line}</div>"
        )
    st.markdown(glass_card(str(record["accession"]), card_content), unsafe_allow_html=True)

    for warning in record.get("warnings") or []:
        st.warning(str(warning))

    if gene_meta:
        render_metric_grid(
            [
                ("GeneID", str(gene_meta.get("id") or record["accession"]), ""),
                ("Symbol", str(gene_meta.get("symbol") or "-"), ""),
                ("Chromosome", str(gene_meta.get("chromosome") or "-"), ""),
                ("Map", str(gene_meta.get("map_location") or "-"), ""),
            ]
        )
        if gene_meta.get("aliases"):
            aliases = ", ".join(str(item) for item in gene_meta["aliases"])
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);'>Aliases: "
                f"{html_escape(aliases)}</span>",
                unsafe_allow_html=True,
            )
        if gene_meta.get("summary"):
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);'>{html_escape(gene_meta['summary'])}</span>",
                unsafe_allow_html=True,
            )
    else:
        reported_length = int(record.get("length") or 0)
        if reported_length > 0:
            length_unit = "aa" if resolved_db == "protein" else "bp"
            length_value = f"{reported_length:,}"
        else:
            length_unit = ""
            length_value = "undefined"
        render_metric_grid(
            [
                ("Length", length_value, length_unit),
                ("Feature Count", str(len(record["features"])), ""),
            ]
        )

    if record.get("related_records"):
        if _section_toggle("Related Records", key="ncbi_related", expanded=True):
            with st.container(border=True):
                related_rows = [
                    {
                        "Database": item.get("database", ""),
                        "ID": item.get("id", ""),
                        "Role": item.get("role", ""),
                        "Detail": item.get("detail", ""),
                    }
                    for item in record["related_records"]
                    if isinstance(item, dict)
                ]
                st.dataframe(pd.DataFrame(related_rows), width="stretch", hide_index=True)

    if explorer.get("references") and _section_toggle(
        "References", key="ncbi_references"
    ):
        with st.container(border=True):
            ref_rows = []
            for item in explorer["references"]:
                if not isinstance(item, dict):
                    continue
                ref_rows.append(
                    {
                        "Title": item.get("title") or "",
                        "Authors": item.get("authors") or "",
                        "Journal": item.get("journal") or "",
                        "PubMed": item.get("pubmed") or "",
                    }
                )
            if ref_rows:
                st.dataframe(pd.DataFrame(ref_rows), width="stretch", hide_index=True)

    if _section_toggle("Features Table", key="ncbi_features"):
        with st.container(border=True):
            annotated = ncbi_fetch.extract_annotated_features(record)
            if not record["features"]:
                if record_kind == "gene":
                    st.info(
                        "Gene records do not include GenBank sequence features. "
                        "Use a related nucleotide accession to inspect CDS."
                    )
                else:
                    st.info("This record has no annotated features.")
            else:
                feature_page_size = 40
                total_features = len(annotated) if annotated else len(record["features"])
                page = st.number_input(
                    "Feature page",
                    min_value=1,
                    max_value=max(1, math.ceil(total_features / feature_page_size)),
                    value=1,
                    step=1,
                    key="ncbi_feature_page",
                )
                if annotated:
                    start = (int(page) - 1) * feature_page_size
                    page_items = annotated[start : start + feature_page_size]
                    feature_rows = [
                        {
                            "Type": item.get("type"),
                            "Status": item.get("status"),
                            "Strand": item.get("strand"),
                            "Join": item.get("is_join"),
                            "Complement": item.get("is_complement"),
                            "Start": item.get("start")
                            if item.get("start") is not None
                            else "UNAVAILABLE",
                            "End": item.get("end")
                            if item.get("end") is not None
                            else "UNAVAILABLE",
                            "Gene": item.get("gene") or "",
                            "Product": item.get("product") or "",
                            "Location": item.get("location"),
                            "Note": item.get("unavailable_reason") or "",
                        }
                        for item in page_items
                    ]
                    st.dataframe(pd.DataFrame(feature_rows), width="stretch")
                    recoverable = [
                        item
                        for item in page_items
                        if item.get("status") == "AVAILABLE" and item.get("segments")
                    ]
                    if recoverable and sequence_available:
                        labels = [
                            f"{item.get('type')} {item.get('start')}-{item.get('end')} "
                            f"{item.get('gene') or item.get('product') or ''}".strip()
                            for item in recoverable
                        ]
                        chosen_feature = st.selectbox(
                            "Feature region for later analysis",
                            labels,
                            key="ncbi_feature_choice",
                        )
                        if st.button(
                            "Store feature coordinates in workspace",
                            key="ncbi_feature_span",
                        ):
                            selected = recoverable[labels.index(chosen_feature)]
                            try:
                                spans = ncbi_fetch.feature_spans_for_analysis(
                                    selected, str(record.get("sequence") or "")
                                )
                            except ValueError as exc:
                                st.error(str(exc))
                            else:
                                st.session_state["ncbi_feature_spans"] = spans
                                st.session_state["ncbi_feature_source"] = {
                                    "accession": record.get("accession"),
                                    "type": selected.get("type"),
                                    "status": "RETRIEVED",
                                    "database": resolved_db,
                                    "sequence_hash": provenance.sequence_digest(
                                        str(record.get("sequence") or "")
                                    )
                                    if record.get("sequence")
                                    else "",
                                }
                                st.success(
                                    f"Stored {len(spans)} retrieved feature span(s). "
                                    "Coordinates are NCBI GenBank, not CRISPR sites."
                                )
                else:
                    start = (int(page) - 1) * feature_page_size
                    page_raw = record["features"][start : start + feature_page_size]
                    feature_rows = [
                        {
                            "Type": feature["type"],
                            "Location": feature["location"],
                            "Qualifiers": _format_feature_qualifiers(
                                feature.get("qualifiers")
                            ),
                        }
                        for feature in page_raw
                    ]
                    st.dataframe(pd.DataFrame(feature_rows), width="stretch")
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);font-size:12px;'>Values come from "
                    "the NCBI GenBank record. Feature Start/End are Biopython "
                    "0-based half-open coordinates on the plus strand, not GenBank "
                    "1-based inclusive numbering and not PDB residue numbers. "
                    "Features whose coordinates cannot be recovered are marked "
                    "UNAVAILABLE; no sequence is invented."
                    "</span>",
                    unsafe_allow_html=True,
                )
                stored_spans = st.session_state.get("ncbi_feature_spans")
                feature_source = st.session_state.get("ncbi_feature_source") or {}
                if stored_spans:
                    st.caption(
                        "Source: NCBI | Feature: "
                        + str(feature_source.get("type") or "")
                        + " | Accession: "
                        + str(feature_source.get("accession") or record.get("accession") or "")
                    )
                    if str(resolved_db or "") != "protein":
                        st.caption(
                            "NCBI feature is on a nucleotide record. It is not mapped "
                            "onto a protein 3D structure."
                        )
                        dna_seq = str(st.session_state.get("dna_input") or "")
                        rna_seq = str(st.session_state.get("rna_input") or "")
                        ncbi_seq = str(record.get("sequence") or "")
                        ncbi_hash = provenance.sequence_digest(ncbi_seq) if ncbi_seq else ""
                        if (
                            dna_seq
                            and ncbi_hash
                            and _sequence_hashes_match(
                                ncbi_hash, provenance.sequence_digest(dna_seq)
                            )
                        ):
                            if st.button(
                                "Highlight this NCBI feature on DNA 3D",
                                key="ncbi_dna_3d_highlight",
                            ):
                                from modules import molecule_selection

                                selection = molecule_selection.selection_from_ncbi_spans(
                                    stored_spans,
                                    sequence=dna_seq,
                                    molecule="DNA",
                                    accession=record.get("accession"),
                                )
                                _store_molecule_selection(selection)
                                st.caption(
                                    "NCBI feature spans mapped to the DNA analysis sequence. "
                                    "Open DNA 3D. Partial mapping is shown as partial."
                                )
                                st.rerun()
                        elif (
                            rna_seq
                            and ncbi_hash
                            and _sequence_hashes_match(
                                ncbi_hash, provenance.sequence_digest(rna_seq)
                            )
                        ):
                            if st.button(
                                "Highlight this NCBI feature on RNA 3D",
                                key="ncbi_rna_3d_highlight",
                            ):
                                from modules import molecule_selection

                                selection = molecule_selection.selection_from_ncbi_spans(
                                    stored_spans,
                                    sequence=rna_seq,
                                    molecule="RNA",
                                    accession=record.get("accession"),
                                )
                                _store_molecule_selection(selection)
                                st.caption(
                                    "NCBI feature spans mapped to the RNA analysis sequence. "
                                    "Open RNA 3D. Partial mapping is shown as partial."
                                )
                                st.rerun()
                        elif dna_seq or rna_seq:
                            st.caption(
                                "NCBI nucleotide sequence does not match the DNA/RNA tab "
                                "sequence. No 3D mapping is invented."
                            )
                        else:
                            st.caption(
                                "Analyze the matching NCBI nucleotide in the DNA or RNA "
                                "tab to highlight this feature on 3D."
                            )
                    else:
                        protein_seq = str(st.session_state.get("protein_input") or "")
                        ncbi_seq = str(record.get("sequence") or "")
                        if (
                            protein_seq
                            and ncbi_seq
                            and _sequence_hashes_match(
                                provenance.sequence_digest(ncbi_seq),
                                provenance.sequence_digest(protein_seq),
                            )
                        ):
                            if st.button(
                                "Highlight this NCBI feature on protein 3D",
                                key="ncbi_3d_highlight",
                            ):
                                from modules import molecule_selection

                                selection = molecule_selection.selection_from_ncbi_spans(
                                    stored_spans,
                                    sequence=protein_seq,
                                    molecule="PROTEIN",
                                    accession=record.get("accession"),
                                )
                                _store_molecule_selection(selection)
                                st.caption(
                                    "NCBI feature spans mapped to the analysis sequence. "
                                    "Open Protein Structure and View 3D Structure. "
                                    "Partial mapping is shown as partial."
                                )
                                st.rerun()
                        elif protein_seq:
                            st.caption(
                                "NCBI protein sequence does not match the Protein tab sequence. "
                                "No 3D mapping is invented."
                            )
                        else:
                            st.caption(
                                "Analyze the matching NCBI protein in the Protein tab to "
                                "highlight this feature on 3D."
                            )

    if _section_toggle("Sequence", key="ncbi_sequence"):
        with st.container(border=True):
            official_length = int(record.get("length") or 0)
            downloaded = str(record.get("sequence") or "")
            if sequence_available and downloaded:
                seq_type = "PROT" if resolved_db == "protein" else "DNA"
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Downloaded from "
                    f"NCBI Entrez: {len(downloaded):,} "
                    f"{'aa' if seq_type == 'PROT' else 'bp'}.</span>",
                    unsafe_allow_html=True,
                )
                st.markdown(
                    sequence_display(downloaded, seq_type), unsafe_allow_html=True
                )
                if st.button("Send sequence to workspace", key="ncbi_to_workspace"):
                    try:
                        payload = ncbi_fetch.sequence_for_workspace(record)
                    except ValueError as exc:
                        st.error(str(exc))
                    else:
                        _store_workspace(
                            payload["molecule"],
                            payload["sequence"],
                            source=payload["source"],
                            identifier=payload["version"] or payload["accession"],
                            accession=payload["accession"],
                            version=payload["version"],
                            organism=payload["organism"],
                        )
                        widget_key = {
                            "DNA": "dna_text",
                            "RNA": "rna_text",
                            "PROTEIN": "prot_text",
                        }.get(payload["molecule"])
                        if widget_key:
                            st.session_state[widget_key] = payload["sequence"]
                        input_key = {
                            "DNA": "dna_input",
                            "RNA": "rna_input",
                            "PROTEIN": "protein_input",
                        }.get(payload["molecule"])
                        if input_key:
                            st.session_state[input_key] = payload["sequence"]
                        st.session_state["blast_query_prefill"] = payload["sequence"]
                        st.success(
                            "Sequence copied to the workspace without recoding. "
                            f"Hash {payload['hash'][:12]}."
                        )
                if st.button("Add this sequence to MSA set", key="ncbi_to_msa"):
                    try:
                        member = msa.member_from_ncbi_record(record)
                    except msa.MsaError as exc:
                        st.error(str(exc))
                    else:
                        collection = st.session_state.setdefault("msa_collection", [])
                        if len(collection) >= msa.MAX_SEQUENCES:
                            st.error(
                                f"MSA set is limited to {msa.MAX_SEQUENCES} sequences."
                            )
                        else:
                            collection.append(member)
                            st.session_state["msa_collection"] = collection
                            st.success(
                                f"Added {member['identifier']} to the MSA set. "
                                "This is not an MSA until Clustal Omega/MAFFT/MUSCLE runs."
                            )
            else:
                unit = "aa" if resolved_db == "protein" else "bp"
                length_text = (
                    f"{official_length:,} {unit}"
                    if official_length > 0
                    else "not deposited in this accession"
                )
                st.info(
                    f"Sequence bases were not downloaded. NCBI reports length "
                    f"{length_text}. Record kind: {record_kind}. HelixScope does "
                    "not invent nucleotides or residues."
                )
                if record_kind == "large_sequence":
                    st.markdown(
                        "<span style='color:var(--hs-text-secondary);font-size:12px;'>To inspect "
                        "bases, use a smaller accession (mRNA, CDS or protein), "
                        "not a chromosome-scale WGS scaffold.</span>",
                        unsafe_allow_html=True,
                    )

    st.download_button(
        "Download record (JSON)",
        data=json.dumps(record, indent=2, default=str).encode("utf-8"),
        file_name=f"{safe_download_filename(record.get('accession'), 'ncbi_record')}.json",
        mime="application/json",
        key="ncbi_download",
    )


def _render_blast_error(exc: Exception) -> None:
    """Exibe falha BLAST classificada, nunca como NO_HITS (interno)."""
    if isinstance(exc, blast_search.BlastError):
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(str(exc))
        st.markdown(
            "<span style='color:var(--hs-text-secondary);font-size:12px;'>This failure is not "
            "NO_HITS. No fabricated homology was inserted.</span>",
            unsafe_allow_html=True,
        )
        return
    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
    st.error(str(exc))


def _blast_alignment_html(hit: dict) -> str:
    """Monta o visualizador de alinhamento BLAST com texto escapado."""
    qseq = html_escape(hit.get("qseq") or "")
    midline = html_escape(hit.get("midline") or "")
    hseq = html_escape(hit.get("hseq") or "")
    q_from = html_escape(hit.get("query_from"))
    q_to = html_escape(hit.get("query_to"))
    h_from = html_escape(hit.get("hit_from"))
    h_to = html_escape(hit.get("hit_to"))
    acc = html_escape(hit.get("accession") or hit.get("hit_id"))
    return (
        "<pre style='white-space:pre-wrap;word-break:break-all;color:var(--hs-text);"
        "font-size:12px;line-height:1.45;'>"
        f"Query  {q_from}-{q_to}\n{qseq}\n{midline}\n{hseq}\n"
        f"Sbjct  {h_from}-{h_to}  {acc}"
        "</pre>"
    )


def render_blast_search() -> None:
    """Submete e acompanha um job NCBI BLAST real (Common URL API).

    Args:
        Nenhum.

    Returns:
        None. Escreve a aba BLAST.

    Raises:
        Nenhum.
    """
    availability = blast_search.ncbi_blast_availability()
    local = blast_search.local_blast_availability()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(availability['reason'])} "
        f"{html_escape(local.get('reason'))} "
        "Sequence similarity is not confirmed biological function. "
        "Pairwise Smith-Waterman in the Alignment tab is not BLAST. "
        "Remote NCBI BLAST is never silently replaced by BLAST+ local. "
        f"BLAST+ executables: "
        f"{'detected' if local.get('executables_detected') else 'not installed'}. "
        f"Local database: "
        f"{'AVAILABLE' if (local.get('database') or {}).get('available') else 'UNAVAILABLE'}."
        "</span>",
        unsafe_allow_html=True,
    )
    if not str(st.session_state.get("blast_email") or "").strip():
        ncbi_email = str(st.session_state.get("ncbi_email") or "")
        if ncbi_email:
            st.session_state["blast_email"] = ncbi_email
    prefill = st.session_state.pop("blast_query_prefill", None)
    if prefill:
        st.session_state["blast_query"] = prefill
    if not str(st.session_state.get("blast_query") or "").strip():
        for key in ("workspace_dna", "workspace_rna", "workspace_protein"):
            seq = str(st.session_state.get(key) or "")
            if seq:
                st.session_state["blast_query"] = seq
                break

    st.markdown("#### Query")
    email = st.text_input(
        "NCBI BLAST email",
        key="blast_email",
        placeholder="your.email@university.edu",
    )
    query = st.text_area("BLAST query", height=160, key="blast_query")
    info = dna_analysis.validate_sequence(query or "")
    if query.strip() and not info["is_valid"]:
        st.warning(info.get("rejection_reason") or "Query is not a valid molecule.")
        programs: list = []
    elif info["is_valid"]:
        programs = blast_search.programs_for_molecule(info["type"])
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Detected molecule: "
            f"{html_escape(info['type'])} ({info['length']} residues). "
            "RNA queries are recoded U to T only for blastn/blastx/tblastx, "
            "and that recoding is recorded in job provenance.</span>",
            unsafe_allow_html=True,
        )
    else:
        programs = []
        st.info("Paste a DNA, RNA or protein query to see compatible BLAST programs.")

    if not programs:
        program = ""
        databases: tuple = ()
    else:
        st.markdown("#### Remote configuration")
        program = st.selectbox("BLAST program", programs, key="blast_program")
        databases = blast_search.databases_for_program(program)
    database = (
        st.selectbox("NCBI BLAST database", list(databases), key="blast_database")
        if databases
        else ""
    )
    expect = st.number_input(
        "Expect threshold (sent to NCBI, not computed here)",
        min_value=1e-20,
        max_value=1000.0,
        value=10.0,
        format="%.4g",
        key="blast_expect",
    )
    hitlist = st.number_input(
        "Hit list size",
        min_value=1,
        max_value=int(blast_search.MAX_HITLIST_SIZE),
        value=int(blast_search.DEFAULT_HITLIST_SIZE),
        step=1,
        key="blast_hitlist",
    )
    blast_backend_options = ["ncbi_remote"]
    if local.get("available"):
        blast_backend_options.append("blastplus_local")
    else:
        st.caption(
            "Local BLAST+ search is UNAVAILABLE until an official BLAST+ binary "
            "and HELIXSCOPE_BLAST_DB prefix with index files are present. "
            "Missing local search does not start a remote job."
        )
    blast_backend = st.radio(
        "BLAST backend",
        blast_backend_options,
        format_func=lambda item: (
            "NCBI remote BLAST (Blast.cgi)"
            if item == "ncbi_remote"
            else f"BLAST+ local ({local.get('version') or 'version not reported'}; declared DB)"
        ),
        key="blast_backend",
    )
    if blast_backend == "blastplus_local" and program:
        if st.button("Run local BLAST+", key="blast_local_run"):
            try:
                prepared = blast_search.prepare_query(query, program)
                db_identity = str((local.get("database") or {}).get("identity") or "")
                key = blast_search.cache_key(
                    query_hash=str(prepared["query_hash"]),
                    program=program,
                    database=db_identity,
                    expect=float(expect),
                    hitlist_size=int(hitlist),
                    backend="blastplus_local",
                    database_identity=db_identity,
                    tool_version=str(local.get("version") or ""),
                )
                cache = st.session_state.setdefault("blast_result_cache", {})
                cached = cache.get(key)
                if isinstance(cached, dict):
                    st.session_state["blast_result"] = blast_search.mark_cached_result(cached)
                    st.session_state.pop("blast_job", None)
                else:
                    result = blast_search.run_local_blast(
                        program=program,
                        query=query,
                        molecule=str(info.get("type") or ""),
                    )
                    st.session_state["blast_result"] = result
                    st.session_state.pop("blast_job", None)
                    cache[key] = result
            except blast_search.BlastError as exc:
                st.markdown(status_badge(exc.category), unsafe_allow_html=True)
                st.error(str(exc))
    st.markdown("#### Execution")
    submit = st.button(
        "Submit NCBI BLAST",
        key="blast_submit",
        disabled=blast_backend != "ncbi_remote",
    )
    if submit:
        if not email.strip():
            st.error("Enter a contact email before submitting BLAST.")
        elif not program:
            st.error("Choose a BLAST program compatible with the query molecule.")
        else:
            existing = st.session_state.get("blast_job")
            if isinstance(existing, dict) and existing.get("status") == "WAITING":
                st.error(
                    "A BLAST job is already waiting in this session. "
                    "Check its status or cancel it before submitting another."
                )
            else:
                try:
                    prepared = blast_search.prepare_query(query, program)
                    key = blast_search.cache_key(
                        query_hash=str(prepared["query_hash"]),
                        program=program,
                        database=database,
                        expect=float(expect),
                        hitlist_size=int(hitlist),
                    )
                    cache = st.session_state.setdefault("blast_result_cache", {})
                    cached = cache.get(key)
                    if isinstance(cached, dict):
                        st.session_state["blast_result"] = blast_search.mark_cached_result(
                            cached
                        )
                        st.session_state.pop("blast_job", None)
                    else:
                        job = blast_search.submit_search(
                            query,
                            email.strip(),
                            program=program,
                            database=database,
                            expect=float(expect),
                            hitlist_size=int(hitlist),
                        )
                        st.session_state["blast_job"] = job
                        st.session_state.pop("blast_result", None)
                        st.session_state["blast_cache_key"] = key
                except (blast_search.BlastError, ValueError) as exc:
                    _render_blast_error(exc)

    job = st.session_state.get("blast_job")
    if isinstance(job, dict):
        st.markdown(
            helix_workspace.job_state_html(
                str(job.get("status") or "WAITING"),
                f"RID {job.get('rid')} · {job.get('program')} · {job.get('database')}",
            ),
            unsafe_allow_html=True,
        )
        st.markdown(status_badge(str(job.get("status") or "WAITING")), unsafe_allow_html=True)
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>RID "
            f"{html_escape(job.get('rid'))} | program "
            f"{html_escape(job.get('program'))} | database "
            f"{html_escape(job.get('database'))} | submitted "
            f"{html_escape(job.get('submitted_at_utc'))} | query hash "
            f"{html_escape(str(job.get('query_hash') or '')[:16])}</span>",
            unsafe_allow_html=True,
        )
        wait = blast_search.seconds_until_next_poll(job)
        if wait > 0:
            st.info(
                f"NCBI policy: wait {wait:.0f}s before polling this RID. "
                "The interface is not blocked; click Check BLAST status later."
            )
        poll_cols = st.columns(2)
        with poll_cols[0]:
            if st.button(
                "Check BLAST status",
                key="blast_poll",
                disabled=wait > 0,
            ):
                try:
                    updated = blast_search.poll_search(job, email.strip())
                    st.session_state["blast_job"] = updated
                    if updated.get("status") == "READY":
                        result = blast_search.retrieve_search(updated, email.strip())
                        st.session_state["blast_result"] = result
                        key = st.session_state.get("blast_cache_key")
                        if key:
                            cache = st.session_state.setdefault("blast_result_cache", {})
                            if len(cache) >= blast_search.MAX_SESSION_CACHE:
                                oldest = next(iter(cache))
                                cache.pop(oldest, None)
                            cache[key] = result
                except (blast_search.BlastError, ValueError) as exc:
                    _render_blast_error(exc)
        with poll_cols[1]:
            if st.button("Cancel local BLAST job", key="blast_cancel"):
                st.session_state.pop("blast_job", None)
                st.rerun()

    result = st.session_state.get("blast_result")
    if isinstance(result, dict):
        current_hash = ""
        if program:
            try:
                current_hash = str(blast_search.prepare_query(query, program).get("query_hash") or "")
            except (blast_search.BlastError, ValueError):
                current_hash = ""
        same_program = str(result.get("program") or "").lower() == str(program or "").lower()
        same_database = str(result.get("database") or "") == str(database or "")
        if (
            not same_program
            or not same_database
            or not _sequence_hashes_match(result.get("query_hash"), current_hash)
        ):
            _hide_stale_result(
                "Stored BLAST result belongs to a different query, program or database. "
                "Submit a new search to refresh."
            )
            result = None
    if not isinstance(result, dict):
        return

    status = str(result.get("status") or "")
    st.markdown(status_badge(status if status else "ERROR"), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.result_header_html(
            f"BLAST {result.get('program') or ''} {result.get('database') or ''}".strip(),
            status=status if status else "ERROR",
            lines=[
                f"RID {result.get('rid') or 'Unavailable'}",
                f"version {result.get('blast_version') or 'not reported'}",
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation("blast", result)
    cache_status = str(result.get("cache_status") or "live")
    if cache_status == "cached":
        st.markdown(status_badge("CACHED"), unsafe_allow_html=True)
    provenance_line = (
        f"Source: {html_escape(result.get('source') or 'NCBI BLAST')} | "
        f"program {html_escape(result.get('program'))} | "
        f"database {html_escape(result.get('database'))} | "
        f"RID {html_escape(result.get('rid'))} | "
        f"BLAST version {html_escape(result.get('blast_version') or 'not reported')} | "
        f"retrieved {html_escape(result.get('retrieved_at_utc'))} | "
        f"cache {html_escape(cache_status)}"
    )
    if result.get("served_from_cache_at_utc"):
        provenance_line += (
            f" | served from cache {html_escape(result.get('served_from_cache_at_utc'))}"
        )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{provenance_line}</span>",
        unsafe_allow_html=True,
    )
    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>BLAST reports sequence "
        "similarity, not a confirmed gene function, domain annotation or "
        "orthology call. Identity % is NCBI Hsp_identity / Hsp_align-len. "
        "Query coverage % uses NCBI 1-based HSP coordinates "
        "(|query_to - query_from| + 1) / query_length. These are not "
        "Needleman-Wunsch identity or MSA column identity.</span>",
        unsafe_allow_html=True,
    )
    st.caption(
        "BLAST is not a structure source. A later PDB or AlphaFold retrieval for "
        "a BLAST subject is a separate NCBI/RCSB/AlphaFold fetch with its own "
        "provenance. There is no BLAST structure."
    )
    for warning in result.get("parse_warnings") or []:
        st.warning(str(warning))

    hits = list(result.get("hits") or [])
    if status == "NO_HITS":
        st.info(
            "NCBI BLAST reported no hits for this query and database. "
            "This is NO_HITS, not a timeout or service failure."
        )
    elif not hits:
        st.error("BLAST result has no hits and is not marked NO_HITS.")
        return
    else:
        query_len = int(result.get("query_length") or 0)
        valid_hits = []
        for hit in hits:
            if not scientific_checks.blast_hsp_is_valid(
                hit, query_len=query_len, hit_len=int(hit.get("hit_length") or 0)
            ):
                st.error(
                    "A BLAST hit failed validation and was not displayed. "
                    "This is ERROR, not NO_HITS."
                )
                continue
            valid_hits.append(hit)
        if not valid_hits:
            st.error("All BLAST hits failed display validation.")
            return
        page = st.number_input(
            "Hit page",
            min_value=1,
            max_value=max(1, math.ceil(len(valid_hits) / blast_search.MAX_DISPLAY_HITS)),
            value=1,
            step=1,
            key="blast_hit_page",
        )
        paged = blast_search.paginate_hits(valid_hits, int(page))
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
            f"{paged['n_hits']} hits | page {paged['page']} of {paged['n_pages']}"
            "</span>",
            unsafe_allow_html=True,
        )
        table_rows = [
            {
                "Accession": item.get("accession"),
                "Organism": item.get("organism") or "",
                "Bit score": item.get("bit_score"),
                "E-value": item.get("evalue"),
                "Identities": item.get("identities"),
                "Alignment length": item.get("alignment_length"),
                "Identity %": item.get("identity_pct"),
                "Query coverage %": item.get("query_coverage_pct"),
                "Description": item.get("description"),
            }
            for item in paged["items"]
        ]
        st.dataframe(pd.DataFrame(table_rows), width="stretch", hide_index=True)
        hit_labels = [
            f"{item.get('accession')} | {str(item.get('description') or '')[:60]}"
            for item in paged["items"]
        ]
        selected_hits = st.multiselect(
            "Select BLAST hits to retrieve full NCBI sequences for the MSA set",
            hit_labels,
            key="blast_msa_hits",
            help="HSP qseq/hseq is not the full subject. NCBI fetch is required. Not an MSA.",
        )
        if st.button("Fetch selected subjects into MSA set", key="blast_to_msa"):
            if status not in {"READY"}:
                st.error(
                    "BLAST subjects can be fetched only from a READY result. "
                    "No sequence was invented from a failed BLAST job."
                )
            elif not email.strip():
                st.error("Enter the NCBI BLAST email before fetching subject records.")
            elif not selected_hits:
                st.error("Select at least one BLAST hit. Hits are not added automatically.")
            elif len(selected_hits) > msa.MAX_BLAST_HITS_TO_FETCH:
                st.error(
                    f"At most {msa.MAX_BLAST_HITS_TO_FETCH} BLAST subjects can be fetched at once."
                )
            else:
                collection = st.session_state.setdefault("msa_collection", [])
                for label in selected_hits:
                    item = paged["items"][hit_labels.index(label)]
                    try:
                        pending = msa.pending_member_from_blast_hit(item, result)
                    except msa.MsaError as exc:
                        st.error(str(exc))
                        continue
                    if not _ncbi_allow_request():
                        break
                    try:
                        record = _fetch_by_accession(
                            str(pending["accession"]),
                            "nucleotide" if str(result.get("program") or "") in {
                                "blastn",
                                "tblastn",
                                "tblastx",
                            } else "protein",
                            email.strip(),
                        )
                    except (ncbi_fetch.NCBIQueryError, ValueError, RuntimeError) as exc:
                        st.error(
                            f"Sequence unavailable for {pending['accession']}: {exc}. "
                            "The HSP was not used as a substitute."
                        )
                        continue
                    try:
                        member = msa.complete_blast_member_with_ncbi(pending, record)
                    except msa.MsaError as exc:
                        st.error(str(exc))
                        continue
                    if len(collection) >= msa.MAX_SEQUENCES:
                        st.error(f"MSA set is limited to {msa.MAX_SEQUENCES} sequences.")
                        break
                    collection.append(member)
                    st.success(
                        f"Retrieved NCBI sequence {member['identifier']} for BLAST hit "
                        f"{pending['accession']}. This is not an MSA."
                    )
                st.session_state["msa_collection"] = collection
        try:
            st.plotly_chart(
                charts.blast_hit_ranking_chart(valid_hits),
                width="stretch",
            )
            st.plotly_chart(charts.blast_evalue_chart(valid_hits), width="stretch")
            st.plotly_chart(charts.blast_identity_chart(valid_hits), width="stretch")
            st.plotly_chart(charts.blast_coverage_chart(valid_hits), width="stretch")
            if query_len > 0:
                st.plotly_chart(
                    charts.blast_alignment_overview(valid_hits, query_len),
                    width="stretch",
                )
        except ValueError as exc:
            st.warning(str(exc))
        labels = [
            f"{item.get('accession')} E={item.get('evalue')}"
            for item in paged["items"]
        ]
        chosen = st.selectbox("Alignment viewer", labels, key="blast_align_choice")
        selected = paged["items"][labels.index(chosen)]
        st.markdown(_blast_alignment_html(selected), unsafe_allow_html=True)

    csv_text = blast_search.export_hits_table(result)
    bundle = blast_search.export_result_bundle(result)
    export_cols = st.columns(2)
    with export_cols[0]:
        st.download_button(
            "Download BLAST hits (CSV)",
            data=csv_text.encode("utf-8"),
            file_name=f"{safe_download_filename(result.get('rid'), 'blast_hits')}.csv",
            mime="text/csv",
            key="blast_csv",
        )
    with export_cols[1]:
        st.download_button(
            "Download BLAST result (JSON)",
            data=json.dumps(bundle, indent=2).encode("utf-8"),
            file_name=f"{safe_download_filename(result.get('rid'), 'blast_result')}.json",
            mime="application/json",
            key="blast_json",
        )


def _store_msa_result(result: dict) -> None:
    """Guarda o MSA e invalida a arvore se a identidade do alinhamento mudou."""
    previous = st.session_state.get("msa_result")
    previous_hash = ""
    if isinstance(previous, dict):
        previous_hash = str(previous.get("alignment_hash") or "")
    new_hash = str(result.get("alignment_hash") or "")
    st.session_state["msa_result"] = result
    if previous_hash and new_hash and previous_hash != new_hash:
        st.session_state.pop("phylo_result", None)
        st.session_state.pop("phylo_selected_hash", None)
        st.session_state.pop("phylo_selected_tree_id", None)


def _render_phylo_error(exc: Exception) -> None:
    """Exibe falha filogenetica classificada."""
    if isinstance(exc, phylogeny.PhylogenyError):
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(str(exc))
        return
    if isinstance(exc, taxonomy.TaxonomyError):
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.warning(str(exc))
        return
    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
    st.error(str(exc))


def _render_evolution_linked_workspace(msa_result: dict, tree: dict) -> None:
    """TREE | MSA | STRUCTURE using session objects. No new phylogeny job.

    Args:
        msa_result: Completed MSA result already in session.
        tree: Completed phylogeny result already in session.

    Returns:
        None.

    Raises:
        Nenhum.

    Nota:
        Tree-to-residue is not direct. Mapping passes through MSA identity.
        Conservation is not pathogenicity or function.
    """
    st.markdown("#### TREE | MSA | STRUCTURE")
    st.caption(
        "Linked view of objects already in this session. Selecting a column does "
        "not rerun IQ-TREE. Highly conserved is not functionally essential. "
        "There is no direct tree-to-residue shortcut."
    )
    helix_workspace.render_explanation(
        "evolution",
        {
            "status": "COMPUTED",
            "conservation_label": "column conservation in this MSA",
            "method": "MSA column conservation",
            "source": "HelixScope linked TREE | MSA | STRUCTURE view",
        },
    )
    tree_col, msa_col, struct_col = helix_shell.three_column_science()
    with tree_col:
        st.caption("Phylogeny")
        st.markdown(status_badge(str(tree.get("status") or "COMPUTED")), unsafe_allow_html=True)
        st.caption(
            f"{tree.get('method_label')} · model {tree.get('model')} · "
            f"hash {(str(tree.get('alignment_hash') or ''))[:12]}"
        )
        leaf = str(st.session_state.get("phylo_selected_tree_id") or "")
        if leaf:
            st.caption(f"Selected leaf {leaf}. MSA row emphasis uses this identity when present.")
        else:
            st.caption("Select a leaf in the tree above to emphasize the matching MSA row.")
    with msa_col:
        st.caption("MSA")
        length = int(msa_result.get("alignment_length") or 0)
        col_idx = int(st.session_state.get("msa_column") or 0)
        if length > 0:
            col_idx = int(
                st.number_input(
                    "Linked MSA column (0-based)",
                    min_value=0,
                    max_value=max(0, length - 1),
                    value=min(col_idx, max(0, length - 1)),
                    step=1,
                    key="evo_linked_column",
                )
            )
        try:
            detail = msa.column_detail(msa_result, col_idx)
        except msa.MsaError as exc:
            st.caption(str(exc))
        else:
            cons = detail.get("conservation")
            cons_text = (
                "Unavailable"
                if isinstance(cons, float) and math.isnan(cons)
                else f"{float(cons):.3f}"
            )
            st.caption(
                f"Column {detail.get('column')} consensus {detail.get('consensus')} "
                f"class {detail.get('variation_class')} conservation {cons_text}. "
                "Highly conserved position, not functionally essential."
            )
    with struct_col:
        st.caption("Structure")
        envelope = st.session_state.get("prot_struct_result")
        if not isinstance(envelope, dict):
            st.caption("No protein structure in this session. Load one in Protein Analysis.")
            return
        kind_label = str(envelope.get("kind_label") or envelope.get("kind") or "UNAVAILABLE")
        st.markdown(status_badge(kind_label), unsafe_allow_html=True)
        st.caption(f"{envelope.get('structure_id') or 'N/A'} · {kind_label}")
        member = st.text_input(
            "MSA member id",
            value=str(st.session_state.get("evo_member_id") or ""),
            key="evo_member_id",
        )
        if member.strip():
            try:
                mapped = comparative.structure_residue_to_msa_column(
                    msa_result=msa_result,
                    member_id=member.strip(),
                    ungapped_index=0,
                )
            except Exception as exc:
                st.caption(str(exc))
            else:
                st.markdown(status_badge(str(mapped.get("status") or "UNMAPPED")), unsafe_allow_html=True)
                st.caption(
                    f"First ungapped residue maps to MSA column "
                    f"{mapped.get('column')} ({mapped.get('status')}). "
                    "UNMAPPED stays UNMAPPED."
                )


def _render_phylogeny_panel(msa_result: dict) -> None:
    """Inferencia filogenetica explicita a partir do MSA completo."""
    st.markdown("### Phylogenetic inference")
    st.caption(
        "Inference from this completed MSA. Not a taxonomic tree and not true "
        "evolutionary history. Distance is computational. Opening this panel "
        "does not fetch taxonomy or compute a tree."
    )
    try:
        phylogeny.validate_msa_for_phylogeny(msa_result)
    except phylogeny.PhylogenyError as exc:
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.warning(str(exc))
        st.caption("Build phylogenetic tree is blocked until the MSA is valid for inference.")
        return
    availability = phylogeny.tool_availability()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(availability['summary_reason'])} Biopython "
        f"{html_escape(availability.get('biopython_version'))}."
        "</span>",
        unsafe_allow_html=True,
    )
    molecule = str(msa_result.get("molecule") or "DNA")
    method_labels = {
        phylogeny.METHOD_NJ: "Neighbor-Joining (Saitou and Nei 1987)",
        phylogeny.METHOD_UPGMA: "UPGMA (ultrametric / molecular-clock assumption)",
    }
    method_ids = [phylogeny.METHOD_NJ, phylogeny.METHOD_UPGMA]
    if availability["fasttree"]["available"]:
        method_labels[phylogeny.METHOD_FASTTREE] = (
            f"FastTree approximate ML ({availability['fasttree'].get('version') or 'version not reported'})"
        )
        method_ids.append(phylogeny.METHOD_FASTTREE)
    else:
        st.caption(
            "FastTree is UNAVAILABLE (executable not installed). HelixScope "
            "does not reimplement FastTree."
        )
    if availability["iqtree"]["available"]:
        method_labels[phylogeny.METHOD_IQTREE] = (
            f"IQ-TREE ML ({availability['iqtree'].get('version') or 'version not reported'})"
        )
        method_ids.append(phylogeny.METHOD_IQTREE)
    else:
        st.caption(
            "IQ-TREE is UNAVAILABLE (executable not installed). HelixScope "
            "does not reimplement IQ-TREE."
        )
    method = st.selectbox(
        "Phylogenetic method",
        method_ids,
        format_func=lambda item: method_labels[item],
        key="phylo_method",
    )
    distance_options = list(
        availability["distance_methods"].get(molecule) or [phylogeny.DISTANCE_P]
    )
    distance_labels = {
        phylogeny.DISTANCE_P: "p-distance (pairwise deletion)",
        phylogeny.DISTANCE_JC69: "Jukes-Cantor 1969",
        phylogeny.DISTANCE_K2P: "Kimura 1980 (K2P)",
        phylogeny.DISTANCE_BLOSUM62: (
            "Biopython BLOSUM62 scoring-matrix distance (not WAG/LG)"
        ),
    }
    if method in {phylogeny.METHOD_FASTTREE, phylogeny.METHOD_IQTREE}:
        st.caption(
            "Maximum-likelihood tools use their own substitution model. "
            "The NJ distance matrix is not applied."
        )
        if method == phylogeny.METHOD_IQTREE and molecule in {"DNA", "RNA"}:
            distance_options = [phylogeny.DISTANCE_JC69, phylogeny.DISTANCE_K2P]
        elif method == phylogeny.METHOD_IQTREE:
            distance_options = [phylogeny.DISTANCE_P]
    distance = st.selectbox(
        "Distance / nucleotide model",
        distance_options,
        format_func=lambda item: distance_labels.get(item, item),
        key="phylo_distance",
        disabled=method == phylogeny.METHOD_FASTTREE,
    )
    if method == phylogeny.METHOD_UPGMA:
        st.info(
            "UPGMA is rooted and assumes a molecular clock / ultrametric "
            "distances. It is not a universal phylogenetic method."
        )
        rooting = phylogeny.ROOTING_UPGMA
        outgroup = ""
    else:
        rooting = st.selectbox(
            "Rooting",
            [
                phylogeny.ROOTING_UNROOTED,
                phylogeny.ROOTING_MIDPOINT,
                phylogeny.ROOTING_OUTGROUP,
            ],
            format_func=lambda item: {
                phylogeny.ROOTING_UNROOTED: "Unrooted (drawing node is not a biological root)",
                phylogeny.ROOTING_MIDPOINT: "Midpoint rooted (not a biological root)",
                phylogeny.ROOTING_OUTGROUP: "Outgroup rooted (user-selected leaf)",
            }[item],
            key="phylo_rooting",
        )
        outgroup_index = 0
        if rooting == phylogeny.ROOTING_OUTGROUP:
            rows = list(msa_result.get("rows") or [])
            outgroup_index = st.selectbox(
                "Outgroup sequence",
                list(range(len(rows))),
                format_func=lambda i: (
                    f"{rows[i].get('identifier')} "
                    f"({str(rows[i].get('hash') or '')[:8]})"
                ),
                key="phylo_outgroup_index",
            )
    bootstrap = st.number_input(
        "Bootstrap replicates (Felsenstein 1985; 0 = support N/A, not 0)",
        min_value=0,
        max_value=phylogeny.MAX_BOOTSTRAP_REPLICATES,
        value=0,
        step=1,
        key="phylo_bootstrap",
        disabled=method in {phylogeny.METHOD_FASTTREE, phylogeny.METHOD_IQTREE},
    )
    iqtree_model_mode = phylogeny.IQTREE_MODEL_USER
    ufboot_replicates = 0
    alrt_replicates = 0
    if method == phylogeny.METHOD_IQTREE:
        iqtree_model_mode = st.selectbox(
            "IQ-TREE substitution model",
            [phylogeny.IQTREE_MODEL_USER, phylogeny.IQTREE_MODEL_MFP],
            format_func=lambda item: (
                "User model (JC / K2P / WAG mapped from the menu above)"
                if item == phylogeny.IQTREE_MODEL_USER
                else "ModelFinder Plus (MFP; BIC by default)"
            ),
            key="phylo_iqtree_model_mode",
        )
        ufboot_on = st.checkbox(
            "UFBoot (ultrafast bootstrap, 1000 replicates; not Felsenstein bootstrap)",
            value=False,
            key="phylo_ufboot",
        )
        alrt_on = st.checkbox(
            "SH-aLRT (1000 replicates; not UFBoot)",
            value=False,
            key="phylo_alrt",
        )
        ufboot_replicates = phylogeny.MIN_UFBOOT_REPLICATES if ufboot_on else 0
        alrt_replicates = phylogeny.MIN_ALRT_REPLICATES if alrt_on else 0
        st.caption(
            "IQ-TREE is the installed binary. FastTree is a different approximate ML method. "
            "ModelFinder, UFBoot and SH-aLRT run only when requested and only if IQ-TREE is available."
        )
    if st.button("Build phylogenetic tree", key="phylo_build"):
        try:
            outgroup_id = ""
            if method != phylogeny.METHOD_UPGMA and rooting == phylogeny.ROOTING_OUTGROUP:
                records, _transforms = phylogeny.unique_leaf_labels(
                    [str(row.get("identifier") or "") for row in list(msa_result.get("rows") or [])],
                    [str(row.get("hash") or "") for row in list(msa_result.get("rows") or [])],
                )
                if 0 <= int(outgroup_index) < len(records):
                    outgroup_id = records[int(outgroup_index)]["tree_id"]
            cache = st.session_state.setdefault("phylo_cache", {})
            tree = phylogeny.infer_phylogeny(
                msa_result,
                method=str(method),
                distance_model=str(distance),
                rooting=str(rooting),
                outgroup_tree_id=outgroup_id,
                bootstrap_replicates=int(bootstrap),
                cache=cache,
                iqtree_model_mode=str(iqtree_model_mode),
                ufboot_replicates=int(ufboot_replicates),
                alrt_replicates=int(alrt_replicates),
            )
            st.session_state["phylo_result"] = tree
            st.session_state["phylo_msa_hash"] = str(msa_result.get("alignment_hash") or "")
        except phylogeny.PhylogenyError as exc:
            _render_phylo_error(exc)
            st.session_state.pop("phylo_result", None)

    tree = st.session_state.get("phylo_result")
    live_hash = str(msa_result.get("alignment_hash") or "")
    if isinstance(tree, dict) and str(tree.get("alignment_hash") or "") != live_hash:
        _hide_stale_result(
            "Stored phylogenetic tree belongs to a previous MSA. "
            "Build the tree again."
        )
        st.session_state.pop("phylo_result", None)
        tree = None
    if not isinstance(tree, dict):
        return
    st.markdown(status_badge(str(tree.get("status") or "COMPUTED")), unsafe_allow_html=True)
    if str(tree.get("cache_status") or "") == "cached":
        st.markdown(status_badge("CACHED"), unsafe_allow_html=True)
    helix_workspace.render_explanation("phylogeny", tree)
    rooting_info = tree.get("rooting") or {}
    support = tree.get("support") or {}
    st.markdown(
        meta_grid(
            [
                ("Method", tree.get("method_label")),
                ("Model", tree.get("model")),
                ("Tool", tree.get("tool")),
                ("Version", tree.get("tool_version") or "not reported"),
                (
                    "Rooting",
                    f"{rooting_info.get('method')} "
                    f"({'rooted' if rooting_info.get('rooted') else 'unrooted'})",
                ),
                (
                    "Support",
                    "N/A"
                    if str(support.get("status")) != "COMPUTED"
                    else (support.get("label") or "support"),
                ),
                (
                    "Support methods",
                    ", ".join(
                        str(item.get("method") or "")
                        for item in list(support.get("methods") or [])
                    )
                    or "N/A",
                ),
                ("Selected model", tree.get("selected_model") or tree.get("model") or "N/A"),
                ("Elapsed (ms)", f"{float(tree.get('elapsed_ms') or 0):.1f}"),
            ]
        ),
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(tree.get('disclaimer'))} "
        f"{html_escape(rooting_info.get('note'))}"
        "</span>",
        unsafe_allow_html=True,
    )
    if int(tree.get("negative_branch_lengths") or 0) > 0:
        st.warning(
            "Neighbor-Joining produced one or more negative branch lengths. "
            "They are algorithm output and were not replaced by 1.0."
        )
    leaf_ids = [str(item.get("tree_id") or "") for item in list(tree.get("leaves") or [])]
    workspace_seq = str(st.session_state.get(f"workspace_{molecule.lower()}") or "")
    leaf_query = str(st.text_input("Search leaf id / organism / accession", key="phylo_leaf_search") or "").strip().lower()
    if leaf_query:
        filtered_ids = []
        for item in list(tree.get("leaves") or []):
            blob = " ".join(
                [
                    str(item.get("tree_id") or ""),
                    str(item.get("original_id") or ""),
                    str(item.get("accession") or ""),
                    str(item.get("organism_display") or ""),
                ]
            ).lower()
            if leaf_query in blob:
                filtered_ids.append(str(item.get("tree_id") or ""))
        if not filtered_ids:
            st.info("No leaf matches that search. The tree is unchanged.")
            filtered_ids = leaf_ids
    else:
        filtered_ids = leaf_ids
    selected_id = st.selectbox(
        "Select leaf (highlights the MSA sequence with the same hash)",
        filtered_ids,
        key="phylo_leaf_select",
    )
    selected_leaf = phylogeny.leaf_by_tree_id(tree, str(selected_id))
    if selected_leaf:
        st.session_state["phylo_selected_hash"] = selected_leaf.get("sequence_hash")
        try:
            from modules import molecule_selection

            selection = molecule_selection.selection_from_tree_leaf(
                tree,
                str(selected_id),
                sequence=workspace_seq,
                molecule=molecule,
            )
            _store_molecule_selection(selection)
        except Exception as exc:
            st.caption(str(exc))
        tax = selected_leaf.get("taxonomy") or {}
        st.markdown(
            meta_grid(
                [
                    ("Leaf", selected_leaf.get("original_id")),
                    ("Hash", str(selected_leaf.get("sequence_hash") or "")[:16]),
                    ("Accession", selected_leaf.get("accession") or "N/A"),
                    ("Organism", selected_leaf.get("organism_display") or "unknown / user-provided"),
                    ("Taxon ID", (tax or {}).get("taxon_id") or "N/A"),
                    ("Rank", (tax or {}).get("rank") or "N/A"),
                    ("Taxonomy source", (tax or {}).get("source") or "N/A"),
                ]
            ),
            unsafe_allow_html=True,
        )
        availability = molecule_selection.structure_availability_for_hash(
            str(selected_leaf.get("sequence_hash") or ""),
            protein_result=st.session_state.get("prot_struct_result"),
            nucleic_result=_session_nucleic_structure_envelope(),
            complex_result=st.session_state.get("crispr_complex"),
        )
        if availability["status"] == "AVAILABLE":
            st.caption(
                f"Validated structure in this session: {availability.get('kind_label')} "
                f"{availability.get('structure_id') or ''} "
                f"({availability.get('source')}). Open the matching analysis tab "
                "to view 3D. HelixScope did not fetch a PDB from this leaf click."
            )
        else:
            st.caption("No validated structure available.")
        leaf_seq = molecule_selection.ungapped_sequence_for_hash(
            msa_result, str(selected_leaf.get("sequence_hash") or "")
        )
        if leaf_seq and molecule == "PROTEIN":
            if st.button("Send this leaf sequence to Protein analysis", key="phylo_to_protein"):
                _store_workspace(
                    "PROTEIN",
                    leaf_seq,
                    source="phylogeny leaf",
                    identifier=str(selected_leaf.get("original_id") or ""),
                    accession=str(selected_leaf.get("accession") or ""),
                    organism=str(selected_leaf.get("organism_display") or ""),
                )
                st.caption("Stored in the Protein workspace. Open the Protein tab to analyze.")
        elif leaf_seq and molecule == "DNA":
            if st.button("Send this leaf sequence to DNA analysis", key="phylo_to_dna"):
                _store_workspace(
                    "DNA",
                    leaf_seq,
                    source="phylogeny leaf",
                    identifier=str(selected_leaf.get("original_id") or ""),
                    accession=str(selected_leaf.get("accession") or ""),
                    organism=str(selected_leaf.get("organism_display") or ""),
                )
                st.caption("Stored in the DNA workspace. Open the DNA tab to analyze.")
        elif leaf_seq and molecule == "RNA":
            if st.button("Send this leaf sequence to RNA analysis", key="phylo_to_rna"):
                _store_workspace(
                    "RNA",
                    leaf_seq,
                    source="phylogeny leaf",
                    identifier=str(selected_leaf.get("original_id") or ""),
                    accession=str(selected_leaf.get("accession") or ""),
                    organism=str(selected_leaf.get("organism_display") or ""),
                )
                st.caption("Stored in the RNA workspace. Open the RNA tab to analyze.")
    layout_mode = st.selectbox(
        "Tree layout (visual only; topology unchanged)",
        ["rectangular", "radial"],
        key="phylo_layout",
    )
    color_rank = st.selectbox(
        "Color leaves by NCBI rank (visual only; requires retrieved taxonomy)",
        ["none", "species", "genus", "family", "order", "class", "phylum"],
        key="phylo_color_rank",
    )
    color_groups = None
    if color_rank != "none":
        grouping = taxonomy.color_groups_by_rank(tree, str(color_rank))
        color_groups = grouping.get("groups") or {}
        st.caption(str(grouping.get("disclaimer") or ""))
        legend = grouping.get("legend") or {}
        if legend:
            st.caption("Legend: " + ", ".join(f"{html_escape(k)} ({v})" for k, v in legend.items()))
    layouts = tree.get("layouts") or {}
    layout = layouts.get(layout_mode)
    if not layout:
        layout = phylogeny.layout_phylogram(tree, mode=str(layout_mode))
    try:
        figure = charts.phylogenetic_tree_figure(
            layout,
            selected_tree_id=str(selected_id or ""),
            color_groups=color_groups,
            show_branch_lengths=bool(tree.get("branch_lengths_present")),
            show_support=str(support.get("status")) == "COMPUTED",
            rooted_label=(
                "rooted" if rooting_info.get("rooted") else "unrooted"
            ),
        )
        st.plotly_chart(figure, width="stretch", key="phylo_chart")
    except ValueError as exc:
        st.info(str(exc))
    _render_evolution_linked_workspace(msa_result, tree)
    if str(support.get("status")) != "COMPUTED":
        st.caption("Bootstrap support: N/A (no replicates). Absence is not 0.")
    elif support.get("truncated"):
        st.warning(
            "Bootstrap replicate count was truncated to the resource limit. "
            "Support is not from a complete requested set."
        )
    taxonomy_email = str(
        st.session_state.get("msa_email")
        or st.session_state.get("ncbi_email")
        or st.session_state.get("blast_email")
        or ""
    )
    st.caption(
        "Retrieve NCBI Taxonomy only for organisms already declared on MSA "
        "members. This does not recompute the tree. User-provided sequences "
        "without an organism stay unknown / user-provided."
    )
    if st.button("Retrieve NCBI taxonomy for declared organisms", key="phylo_taxonomy"):
        try:
            tax_cache = st.session_state.setdefault("taxonomy_cache", {})
            updated = taxonomy.attach_ncbi_taxonomy(
                tree,
                email=taxonomy_email,
                cache=tax_cache,
            )
            st.session_state["phylo_result"] = updated
            layer = updated.get("taxonomy_layer") or {}
            if int(layer.get("n_failed") or 0) > 0:
                st.warning(
                    "Some taxonomy lookups failed. The phylogenetic tree is unchanged. "
                    + "; ".join(
                        str(item.get("error") or "")
                        for item in list(layer.get("failures") or [])[:4]
                    )
                )
            tree = updated
        except taxonomy.TaxonomyError as exc:
            _render_phylo_error(exc)
    leaf_table = []
    for leaf in list(tree.get("leaves") or []):
        tax = leaf.get("taxonomy") or {}
        leaf_table.append(
            {
                "Leaf": leaf.get("tree_id"),
                "Original": leaf.get("original_id"),
                "Hash": str(leaf.get("sequence_hash") or "")[:12],
                "Accession": leaf.get("accession") or "",
                "Organism": leaf.get("organism_display") or "unknown / user-provided",
                "Taxon ID": (tax or {}).get("taxon_id") or "N/A",
                "Rank": (tax or {}).get("rank") or "N/A",
                "Source": (tax or {}).get("source") or leaf.get("metadata_source") or "N/A",
            }
        )
    st.dataframe(pd.DataFrame(leaf_table), width="stretch", hide_index=True)
    export_tree = st.columns(3)
    report = phylogeny.scientific_report(tree)
    with export_tree[0]:
        st.download_button(
            "Newick",
            data=str(tree.get("newick") or "").encode("utf-8"),
            file_name=f"{safe_download_filename('phylogeny', 'helixscope_tree')}.nwk",
            mime="text/plain",
            key="phylo_nwk",
        )
    with export_tree[1]:
        st.download_button(
            "Tree JSON report",
            data=json.dumps(report, indent=2).encode("utf-8"),
            file_name=f"{safe_download_filename('phylogeny', 'helixscope_tree')}.json",
            mime="application/json",
            key="phylo_json",
        )
    with export_tree[2]:
        st.download_button(
            "Leaf metadata CSV",
            data=phylogeny.export_leaf_csv(tree).encode("utf-8"),
            file_name=f"{safe_download_filename('phylogeny', 'helixscope_tree_leaves')}.csv",
            mime="text/csv",
            key="phylo_csv",
        )


def _render_msa_error(exc: Exception) -> None:
    """Exibe falha MSA classificada, nunca como UNAVAILABLE de sucesso."""
    if isinstance(exc, msa.MsaError):
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(str(exc))
        return
    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
    st.error(str(exc))


def render_msa_analysis() -> None:
    """Conjunto comparativo, MSA real e analise de colunas.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    availability = msa.tool_availability()
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(availability['summary_reason'])} "
        "BLAST hits are not an MSA. Consensus is not a biological sequence. "
        "Conservation is not function. Phylogenetic inference is a separate "
        "explicit action on a completed MSA, not this alignment step. "
        "IQ-TREE and FastTree stay UNAVAILABLE unless those official "
        "executables are installed; HelixScope does not reimplement them."
        "</span>",
        unsafe_allow_html=True,
    )
    st.caption(
        "EMBL-EBI jobs are submitted and polled from this session (minimum "
        f"{int(msa.MIN_POLL_INTERVAL_S)}s between status checks). Local "
        f"clustalo/MAFFT/MUSCLE, when installed, run in-process with a "
        f"{int(msa.LOCAL_TIMEOUT_S)}s timeout. There is no distributed worker pool."
    )
    st.markdown("#### Dataset setup")
    backends = msa.available_backend_ids()
    collection = st.session_state.setdefault("msa_collection", [])
    if not isinstance(collection, list):
        collection = []
        st.session_state["msa_collection"] = collection

    fasta_text = st.text_area(
        "Add sequences (multi-FASTA is kept as multiple records)",
        height=140,
        key="msa_fasta_text",
    )
    if st.button("Add FASTA records to MSA set", key="msa_add_fasta"):
        try:
            records = msa.parse_unaligned_fasta(fasta_text)
            if not records:
                raise msa.MsaError("No FASTA records to add.", "INVALID_INPUT")
            added = 0
            for record in records:
                if len(collection) >= msa.MAX_SEQUENCES:
                    st.error(f"MSA set is limited to {msa.MAX_SEQUENCES} sequences.")
                    break
                member = msa.build_collection_member(
                    sequence=str(record.get("sequence") or ""),
                    identifier=str(record.get("identifier") or ""),
                    source="user input",
                )
                collection.append(member)
                added += 1
            st.session_state["msa_collection"] = collection
            if added:
                st.success(f"Added {added} sequence(s). None were silently dropped.")
        except (msa.MsaError, ValueError) as exc:
            _render_msa_error(exc)

    if st.button("Load textarea as completed prealigned MSA", key="msa_load_prealigned"):
        try:
            aligned_result = msa.import_prealigned_fasta(fasta_text)
            _store_msa_result(aligned_result)
            st.session_state["msa_collection"] = list(aligned_result.get("input_members") or [])
            st.success(
                "Loaded user-declared prealigned FASTA as a completed MSA. "
                "Columns were not inferred by Clustal/MAFFT/MUSCLE."
            )
            st.rerun()
        except (msa.MsaError, ValueError) as exc:
            _render_msa_error(exc)

    ws_options = []
    for molecule, key in (
        ("DNA", "workspace_dna"),
        ("RNA", "workspace_rna"),
        ("PROTEIN", "workspace_protein"),
    ):
        seq = str(st.session_state.get(key) or "")
        if seq:
            ws_options.append(molecule)
    if ws_options:
        chosen_ws = st.multiselect(
            "Add from analysis workspace",
            ws_options,
            key="msa_from_workspace",
        )
        if st.button("Add workspace sequences", key="msa_add_workspace"):
            for molecule in chosen_ws:
                seq = str(st.session_state.get(f"workspace_{molecule.lower()}") or "")
                meta = st.session_state.get(f"workspace_{molecule.lower()}_meta") or {}
                try:
                    member = msa.build_collection_member(
                        sequence=seq,
                        identifier=str(meta.get("identifier") or molecule),
                        source="workspace",
                        molecule=molecule,
                        accession=str(meta.get("accession") or ""),
                        version=str(meta.get("version") or ""),
                        organism=str(meta.get("organism") or ""),
                    )
                except msa.MsaError as exc:
                    _render_msa_error(exc)
                    continue
                if len(collection) >= msa.MAX_SEQUENCES:
                    st.error(f"MSA set is limited to {msa.MAX_SEQUENCES} sequences.")
                    break
                collection.append(member)
            st.session_state["msa_collection"] = collection

    if collection:
        table_rows = []
        for index, item in enumerate(collection):
            table_rows.append(
                {
                    "Order": index,
                    "Identifier": item.get("identifier"),
                    "Source": item.get("source"),
                    "Accession": item.get("accession") or "",
                    "Version": item.get("version") or "",
                    "Organism": item.get("organism") or "",
                    "Molecule": item.get("molecule"),
                    "Length": item.get("length"),
                    "Hash": str(item.get("hash") or "")[:12],
                    "BLAST RID": item.get("blast_rid") or "",
                }
            )
        st.dataframe(pd.DataFrame(table_rows), width="stretch", hide_index=True)
        try:
            preview = msa.validate_collection(collection)
            groups = preview.get("identical_groups") or []
        except msa.MsaError:
            groups = msa.identical_sequence_groups(collection)
        if groups:
            for group in groups:
                st.info(
                    "Identical sequences (same hash), not silently removed: "
                    + ", ".join(str(label) for label in group.get("identifiers") or [])
                )
        move_index = st.number_input(
            "Sequence order index",
            min_value=0,
            max_value=max(0, len(collection) - 1),
            value=0,
            step=1,
            key="msa_order_index",
        )
        cols = st.columns(3)
        with cols[0]:
            if st.button("Move up", key="msa_move_up") and int(move_index) > 0:
                i = int(move_index)
                collection[i - 1], collection[i] = collection[i], collection[i - 1]
                st.session_state["msa_collection"] = collection
                st.rerun()
        with cols[1]:
            if st.button("Move down", key="msa_move_down") and int(move_index) < len(collection) - 1:
                i = int(move_index)
                collection[i + 1], collection[i] = collection[i], collection[i + 1]
                st.session_state["msa_collection"] = collection
                st.rerun()
        with cols[2]:
            if st.button("Remove at index", key="msa_remove"):
                i = int(move_index)
                if 0 <= i < len(collection):
                    collection.pop(i)
                    st.session_state["msa_collection"] = collection
                    st.rerun()
        if st.button("Clear MSA set", key="msa_clear"):
            st.session_state["msa_collection"] = []
            st.session_state.pop("msa_job", None)
            st.session_state.pop("msa_result", None)
            st.session_state.pop("phylo_result", None)
            st.rerun()
    else:
        st.info("Add at least two sequences of the same molecule type.")

    if not str(st.session_state.get("msa_email") or "").strip():
        shared = str(
            st.session_state.get("blast_email")
            or st.session_state.get("ncbi_email")
            or ""
        )
        if shared:
            st.session_state["msa_email"] = shared
    email = st.text_input(
        "Contact email (required for EMBL-EBI Clustal Omega)",
        key="msa_email",
        placeholder="your.email@university.edu",
    )
    st.markdown("#### Engine configuration")
    backend = st.selectbox(
        "MSA backend (only tools that can actually run)",
        backends,
        key="msa_backend",
    )
    submit = st.button("Run multiple sequence alignment", key="msa_submit")
    if submit:
        existing = st.session_state.get("msa_job")
        if isinstance(existing, dict) and existing.get("status") in {"QUEUED", "RUNNING"}:
            st.error("An MSA job is already running in this session. Poll or cancel it.")
        else:
            try:
                validated = msa.validate_collection(collection)
                key = msa.cache_key(
                    hashes=validated["hashes"],
                    backend=str(backend),
                    tool_version="",
                    parameters={"backend": backend, "outfmt": "fa", "order": "input"},
                )
                cache = st.session_state.setdefault("msa_result_cache", {})
                cached = cache.get(key)
                if isinstance(cached, dict):
                    _store_msa_result(msa.mark_cached_result(cached))
                    st.session_state.pop("msa_job", None)
                else:
                    job = msa.submit_msa(
                        collection, email.strip(), backend=str(backend)
                    )
                    st.session_state["msa_job"] = job
                    st.session_state["msa_cache_key"] = key
                    if job.get("status") == "COMPLETED" and job.get("result"):
                        _store_msa_result(job["result"])
                        cache[key] = job["result"]
                    else:
                        st.session_state.pop("msa_result", None)
                        st.session_state.pop("phylo_result", None)
            except (msa.MsaError, ValueError) as exc:
                _render_msa_error(exc)

    job = st.session_state.get("msa_job")
    if isinstance(job, dict):
        st.markdown(
            helix_workspace.job_state_html(
                str(job.get("status") or "QUEUED"),
                f"{job.get('tool')} · {job.get('backend')} · {job.get('job_id') or job.get('id') or ''}",
            ),
            unsafe_allow_html=True,
        )
        st.markdown(status_badge(str(job.get("status") or "QUEUED")), unsafe_allow_html=True)
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>tool "
            f"{html_escape(job.get('tool'))} | backend "
            f"{html_escape(job.get('backend'))} | job "
            f"{html_escape(job.get('job_id'))} | order "
            f"{html_escape(', '.join(job.get('identifiers') or []))}</span>",
            unsafe_allow_html=True,
        )
        if str(job.get("backend")) == "ebi_clustalo" and job.get("status") in {
            "QUEUED",
            "RUNNING",
        }:
            wait = msa.seconds_until_next_poll(job)
            if wait > 0:
                st.info(
                    f"Wait {wait:.0f}s before polling this EMBL-EBI job. "
                    "The interface is not blocked."
                )
            poll_cols = st.columns(2)
            with poll_cols[0]:
                if st.button("Check MSA status", key="msa_poll", disabled=wait > 0):
                    try:
                        updated = msa.poll_msa(job, email.strip())
                        st.session_state["msa_job"] = updated
                        if updated.get("status") == "COMPLETED":
                            result = msa.retrieve_msa(updated, email.strip())
                            _store_msa_result(result)
                            key = st.session_state.get("msa_cache_key")
                            if key:
                                cache = st.session_state.setdefault("msa_result_cache", {})
                                cache[key] = result
                    except (msa.MsaError, ValueError) as exc:
                        _render_msa_error(exc)
            with poll_cols[1]:
                if st.button("Stop tracking this job", key="msa_cancel"):
                    job["status"] = "CANCELLED"
                    st.session_state["msa_job"] = job
                    st.info(
                        "This session will not retrieve the result. "
                        "A remote EMBL-EBI job is not killed from here."
                    )
                    st.rerun()

    result = st.session_state.get("msa_result")
    if isinstance(result, dict):
        current_hashes = [str(item.get("hash") or "") for item in collection]
        stored_hashes = [str(item) for item in list(result.get("input_hashes") or [])]
        if current_hashes != stored_hashes:
            _hide_stale_result(
                "Stored MSA belongs to a different sequence collection. "
                "Run MSA again to refresh."
            )
            st.session_state.pop("phylo_result", None)
            result = None
    if not isinstance(result, dict):
        st.caption("No completed MSA is loaded.")
        st.caption(
            "Run MSA or load a valid pre-aligned FASTA before inferring a "
            "phylogenetic tree. Phylogeny engines are not missing; a tree is "
            "not invented from an empty alignment."
        )
        return
    st.markdown(status_badge(str(result.get("status") or "COMPLETED")), unsafe_allow_html=True)
    if str(result.get("cache_status") or "") == "cached":
        st.markdown(status_badge("CACHED"), unsafe_allow_html=True)
    n_rows = int(result.get("n_sequences") or len(result.get("rows") or []))
    length = int(result.get("alignment_length") or 0)
    st.markdown(
        helix_workspace.result_header_html(
            "Multiple sequence alignment",
            status=str(result.get("status") or "COMPLETED"),
            lines=[
                f"{n_rows} sequences · {length} columns",
                f"{result.get('tool') or 'tool'} {result.get('tool_version') or 'version not reported'}",
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation("msa", result)
    version_text = result.get("tool_version") or "not reported by the tool"
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"tool {html_escape(result.get('tool'))} | version {html_escape(version_text)} | "
        f"method {html_escape(result.get('method'))} | "
        f"length {html_escape(result.get('alignment_length'))} | "
        f"retrieved {html_escape(result.get('retrieved_at_utc'))}"
        "</span>",
        unsafe_allow_html=True,
    )
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(result.get('disclaimer'))}</span>",
        unsafe_allow_html=True,
    )
    row_start = st.number_input(
        "Viewer row start",
        min_value=0,
        max_value=max(0, n_rows - 1),
        value=0,
        step=1,
        key="msa_row_start",
    )
    col_start = st.number_input(
        "Viewer column start",
        min_value=0,
        max_value=max(0, length - 1),
        value=0,
        step=msa.MAX_DISPLAY_COLUMNS,
        key="msa_col_start",
    )
    try:
        window = msa.viewer_window(
            result, row_start=int(row_start), col_start=int(col_start)
        )
    except msa.MsaError as exc:
        _render_msa_error(exc)
        return
    lines = [
        f"columns {window['col_start']}-{window['col_end'] - 1}  "
        f"rows {window['row_start']}-{window['row_end'] - 1}"
    ]
    for row in window["rows"]:
        ident = str(row.get("identifier") or "")[:16].ljust(16)
        mark = ""
        selected_hash = str(st.session_state.get("phylo_selected_hash") or "")
        if selected_hash and str(row.get("hash") or "") == selected_hash:
            mark = "> "
        lines.append(f"{mark}{ident} {row.get('aligned_slice')}")
    lines.append("consensus".ljust(16) + " " + str(window.get("consensus_slice") or ""))
    st.markdown(
        f"<div class='hs-msa-viewer'><pre style='white-space:pre;overflow:auto;color:var(--hs-text);font-size:12px;margin:0;'>"
        + html_escape("\n".join(lines))
        + "</pre></div>",
        unsafe_allow_html=True,
    )
    selected_col = st.number_input(
        "Inspect alignment column (0-based)",
        min_value=0,
        max_value=max(0, length - 1),
        value=int(window["col_start"]),
        step=1,
        key="msa_column",
    )
    try:
        detail = msa.column_detail(result, int(selected_col))
    except msa.MsaError as exc:
        _render_msa_error(exc)
    else:
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Column {detail['column']} | "
            f"consensus {html_escape(detail['consensus'])} | "
            f"class {html_escape(detail['variation_class'])} | "
            f"gaps {detail['n_gaps']} | "
            f"original 0-based positions {html_escape(detail['original_positions_0based'])}"
            "</span>",
            unsafe_allow_html=True,
        )
        residue_rows = [
            {
                "Identifier": ident,
                "Residue": residue,
                "Original_0based": "gap" if pos is None else str(pos),
            }
            for ident, residue, pos in zip(
                detail["identifiers"],
                detail["residues"],
                detail["original_positions_0based"],
            )
        ]
        st.dataframe(pd.DataFrame(residue_rows), width="stretch", hide_index=True)
        cons_score = detail.get("conservation")
        if not (
            isinstance(cons_score, float) and math.isnan(cons_score)
        ) and not scientific_checks.conservation_score_is_valid(float(cons_score)):
            st.error("Conservation score failed validation and was not displayed.")
        else:
            cons_shown = (
                "N/A"
                if isinstance(cons_score, float) and math.isnan(cons_score)
                else f"{float(cons_score):.3f}"
            )
            render_metric_grid(
                [("Column conservation", cons_shown, "")],
                columns=4,
            )
        _render_msa_protein_3d_link(result, int(selected_col))
    consensus = result.get("consensus") or {}
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(consensus.get('disclaimer'))} "
        f"{html_escape(consensus.get('method'))}</span>",
        unsafe_allow_html=True,
    )
    conservation = result.get("conservation") or {}
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
        f"{html_escape(conservation.get('disclaimer'))} "
        f"{html_escape(conservation.get('method'))}</span>",
        unsafe_allow_html=True,
    )
    try:
        cols = list(range(window["col_start"], window["col_end"]))
        scores = list(window.get("conservation_slice") or [])
        cons_chars = list(window.get("consensus_slice") or "")
        st.plotly_chart(
            charts.msa_conservation_chart(cols, scores, cons_chars),
            width="stretch",
        )
    except ValueError:
        st.info("No conservation scores in this window to plot.")
    identity = result.get("identity_matrix") or {}
    matrix = identity.get("matrix") or []
    labels = list(result.get("input_order") or [])
    if matrix:
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>"
            f"{html_escape(identity.get('method'))}. "
            f"Denominator: {html_escape(identity.get('denominator'))}. "
            "Values are ungapped pairwise fractions in [0, 1], not percent "
            "and not Needleman-Wunsch.</span>",
            unsafe_allow_html=True,
        )
        ident_rows = []
        for i, row in enumerate(matrix):
            rec = {"Sequence": labels[i] if i < len(labels) else i}
            for j, value in enumerate(row):
                name = labels[j] if j < len(labels) else str(j)
                if isinstance(value, float) and math.isnan(value):
                    rec[name] = "N/A"
                else:
                    rec[name] = f"{float(value):.4f}"
            ident_rows.append(rec)
        st.dataframe(pd.DataFrame(ident_rows), width="stretch", hide_index=True)
    _render_phylogeny_panel(result)
    export_cols = st.columns(4)
    with export_cols[0]:
        st.download_button(
            "Aligned FASTA",
            data=msa.export_aligned_fasta(result).encode("utf-8"),
            file_name="helixscope_msa.fa",
            mime="text/plain",
            key="msa_fa",
        )
    with export_cols[1]:
        st.download_button(
            "CLUSTAL-like",
            data=msa.export_clustal_like(result).encode("utf-8"),
            file_name="helixscope_msa.aln",
            mime="text/plain",
            key="msa_clu",
        )
    with export_cols[2]:
        st.download_button(
            "Column CSV",
            data=msa.export_column_csv(result).encode("utf-8"),
            file_name="helixscope_msa_columns.csv",
            mime="text/csv",
            key="msa_csv",
        )
    with export_cols[3]:
        st.download_button(
            "MSA JSON",
            data=json.dumps(msa.export_result_bundle(result), indent=2).encode("utf-8"),
            file_name="helixscope_msa.json",
            mime="application/json",
            key="msa_json",
        )


def _risk_variant(level: str) -> str:
    """Mapeia um nivel de risco de off-target para a variante de badge (interno)."""
    return {"Low": "dna", "Moderate": "warn", "High": "err"}.get(level, "warn")


def _int_or_default(value: Any, default: int) -> int:
    """Inteiro com default so para None/vazio. 0 permanece 0."""
    if value is None or value == "":
        return int(default)
    return int(value)


def _guide_quality_flags(guide: dict) -> dict:
    """Resume os sinais de qualidade de um guia enriquecido (uso interno).

    Args:
        guide: Guia produzido por crispr.evaluate_guides.

    Returns:
        Dicionario com "valid_gc" (bool), "poly_t" (bool, presenca de sinal),
        "hairpin_risk" (str) e "high_quality" (bool).

    Raises:
        Nenhum.

    Nota:
        Um guia de alta qualidade combina GC na faixa, ausencia de poli-T, baixo
        risco de grampo e, quando disponivel, especificidade alta.
    """
    gc_raw = guide.get("gc_content")
    if gc_raw is None or (isinstance(gc_raw, float) and math.isnan(gc_raw)):
        valid_gc = False
    else:
        gc = float(gc_raw)
        valid_gc = 40.0 <= gc <= 80.0
    poly_t = bool(guide.get("poly_t", {}).get("has_signal", False))
    hairpin = str(guide.get("self_complementarity", {}).get("risk", "low"))
    proxy = guide.get("specificity_proxy")
    high_quality = (
        valid_gc
        and not poly_t
        and hairpin != "high"
        and float(guide.get("doench_score", 0.0)) >= 0.5
        and (proxy is None or float(proxy) >= 0.8)
    )
    return {
        "valid_gc": valid_gc,
        "poly_t": poly_t,
        "hairpin_risk": hairpin,
        "high_quality": high_quality,
    }


def _render_crispr_model_transparency() -> None:
    """Declara na interface o que e calculado e o que permanece indisponivel."""
    availability = crispr.model_availability()
    with st.container(border=True):
        st.markdown(
            '<p class="module-title" style="font-size:1rem;">What can and cannot '
            "be computed here</p>"
            '<p class="module-desc">Rule Set 2, Azimuth, DeepHF remain unavailable. '
            "Genome-wide MIT/CFD aggregates exist only after a Cas-OFFinder job "
            "status COMPLETED_FULL_REFERENCE on a checksum-verified READY public "
            "assembly. Deposited Cas-guide-DNA complexes (RCSB 4UN3, 4OO8) can be "
            "fetched on an explicit action; HelixScope does not assemble "
            "Cas9+guide+DNA procedurally. Per-site CFD (Doench 2016 matrices) and "
            "MIT/CFD aggregates over a COMPLETED provided-reference search are "
            "computed and labelled. HelixScope will not invent missing numbers. "
            "A user FASTA, workspace DNA or NCBI nucleotide record already in this "
            "session can be searched explicitly; that search is not genome-wide. "
            "Cas-OFFinder runs only if the real binary is installed. Opening this "
            "tab does not download a genome. This deploy has no authentication, "
            "per-user quotas or public/cloud workers.</p>",
            unsafe_allow_html=True,
        )
        for name, info in availability.items():
            color = "#34D399" if info.get("available") else "#F87171"
            state = "AVAILABLE" if info.get("available") else "UNAVAILABLE"
            st.markdown(
                f"<p style='margin:0 0 8px 0;'><strong style='color:{color};'>"
                f"{html_escape(name)}</strong> "
                f"<span style='color:var(--hs-text-secondary);font-size:12px;'>{state}</span>"
                f"<br><span style='color:var(--hs-text-secondary);font-size:13px;'>"
                f"{html_escape(info['reason'])} Requires: "
                f"{html_escape(info['requires'])}.</span></p>",
                unsafe_allow_html=True,
            )


def render_motif_search() -> None:
    """Renderiza a aba de busca de motivos em DNA, RNA ou proteina.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes na pagina.

    Raises:
        Nenhum.
    """
    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>Exact or IUPAC pattern "
        "search on the sequence you provide. Hits are not PWM scores, Pfam "
        "domains or experimentally validated sites.</span>",
        unsafe_allow_html=True,
    )
    molecule = st.selectbox(
        "Molecule",
        ["DNA", "RNA", "PROTEIN"],
        key="motif_molecule",
    )
    workspace_key = {
        "DNA": "workspace_dna",
        "RNA": "workspace_rna",
        "PROTEIN": "workspace_protein",
    }[molecule]
    workspace_seq = str(st.session_state.get(workspace_key) or "")
    widget_key = f"motif_text_{molecule}"
    if not str(st.session_state.get(widget_key) or "").strip() and workspace_seq:
        st.session_state[widget_key] = workspace_seq
    text_value = st.text_area(
        "Sequence",
        height=160,
        key=widget_key,
    )
    pattern = st.text_input(
        "Motif pattern",
        placeholder="GAATTC, NGG  (comma-separated; max 20)",
        key="motif_pattern",
        help="DNA/RNA: IUPAC. Protein: exact letters plus X. Multiple patterns separated by comma or newline.",
    )
    strand_filter = "+"
    if molecule != "PROTEIN":
        strand_choice = st.selectbox(
            "Strand filter",
            ["both", "+", "-"],
            key="motif_strand_filter",
        )
        strand_filter = "" if strand_choice == "both" else strand_choice
    search_clicked = st.button("Search motif", key="motif_search_button")
    if search_clicked:
        try:
            resolved = _resolve_sequence(text_value, None)
        except ValueError as exc:
            st.error(str(exc))
            return
        info = _validate_for(resolved, molecule)
        if not info["is_valid"]:
            st.error(info["rejection_reason"] or "Invalid sequence for this molecule.")
            return
        seq = info["sequence"]
        try:
            patterns = motif_search.parse_motif_patterns(pattern)
            hits = motif_search.find_motifs(seq, patterns, molecule)
            if molecule != "PROTEIN" and strand_filter:
                hits = motif_search.filter_motif_hits(hits, strand=strand_filter)
            for hit in hits:
                start = int(hit["start"])
                end = int(hit["end"])
                if str(hit.get("strand") or "+") == "+":
                    if not scientific_checks.motif_hit_matches_sequence(
                        seq, start, end, str(hit.get("match") or "")
                    ):
                        raise ValueError(
                            "Analysis failed: a motif hit does not match the "
                            "sequence at the reported coordinates."
                        )
                elif not scientific_checks.span_is_valid(start, end, len(seq), base=0):
                    raise ValueError(
                        "Analysis failed: a motif hit falls outside the sequence."
                    )
        except ValueError as exc:
            message = str(exc)
            if "exceeded" in message.lower() or "limited" in message.lower():
                st.error("Analysis limit reached")
                st.info(message)
                return
            st.error(message)
            return
        st.session_state["motif_result"] = {
            "molecule": molecule,
            "pattern": pattern,
            "patterns": patterns,
            "sequence_length": len(seq),
            "sequence_hash": provenance.sequence_digest(seq),
            "hits": hits,
            "status": "COMPUTED",
        }

    packed = st.session_state.get("motif_result")
    current_motif_hash = ""
    try:
        live_motif = _resolve_sequence(text_value, None)
        live_info = _validate_for(live_motif, molecule)
        if live_info.get("is_valid"):
            current_motif_hash = provenance.sequence_digest(str(live_info["sequence"]))
    except (ValueError, RuntimeError):
        current_motif_hash = ""
    if isinstance(packed, dict) and (
        packed.get("molecule") != molecule
        or str(packed.get("pattern") or "") != str(pattern or "")
        or not _sequence_hashes_match(packed.get("sequence_hash"), current_motif_hash)
    ):
        _hide_stale_result(
            "Stored motif hits belong to a different sequence, molecule or pattern. "
            "Search again to refresh."
        )
        packed = None
    if not isinstance(packed, dict):
        return
    st.markdown(status_badge("COMPUTED"), unsafe_allow_html=True)
    hits = list(packed.get("hits") or [])
    summary = motif_search.summarize_motif_hits(hits)
    st.markdown(
        helix_workspace.result_header_html(
            "Motif search",
            status="COMPUTED",
            lines=[
                f"{summary['unique_loci']} unique loci",
                (
                    f"{summary['forward_hits']} forward, "
                    f"{summary['reverse_complement_hits']} reverse-complement"
                ),
                f"{packed.get('molecule')} · {packed.get('sequence_length')} residues",
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation(
        "motif",
        {
            "status": "COMPUTED",
            "n_hits": summary["literal_hits"],
            "unique_loci": summary["unique_loci"],
            "forward_hits": summary["forward_hits"],
            "reverse_complement_hits": summary["reverse_complement_hits"],
            "pattern": packed.get("pattern"),
            "method": "IUPAC motif scan",
            "source": "HelixScope motif search",
        },
    )
    render_metric_grid(
        [
            ("Unique loci", str(summary["unique_loci"]), ""),
            ("Forward", str(summary["forward_hits"]), ""),
            ("Reverse complement", str(summary["reverse_complement_hits"]), ""),
        ],
        columns=3,
    )
    if not hits:
        st.info("No motif matches in the provided sequence.")
        return
    rows = [
        {
            "Pattern": item.get("pattern", packed.get("pattern", "")),
            "Start_0based": item["start"],
            "End_exclusive": item["end"],
            "Strand": item["strand"] if str(packed.get("molecule") or "") != "PROTEIN" else "n/a",
            "Match": item["match"],
            "Method": "IUPAC/exact pattern",
        }
        for item in hits
    ]
    table = pd.DataFrame(rows)
    st.dataframe(table, width="stretch", hide_index=True)
    if str(packed.get("molecule") or "") == "PROTEIN":
        protein_seq = str(st.session_state.get("protein_input") or "")
        if protein_seq and _sequence_hashes_match(
            packed.get("sequence_hash"), provenance.sequence_digest(protein_seq)
        ):
            labels = [
                f"{item.get('start')}-{item.get('end')} {item.get('match')}"
                for item in hits
            ]
            chosen = st.selectbox("Motif hit for protein 3D", labels, key="motif_3d_hit")
            if st.button("Highlight this motif on protein 3D", key="motif_3d_highlight"):
                from modules import molecule_selection

                hit = hits[labels.index(chosen)]
                selection = molecule_selection.selection_from_motif_hit(
                    hit, molecule="PROTEIN", sequence=protein_seq
                )
                _store_molecule_selection(selection)
                st.caption(
                    "Motif span mapped to the analysis sequence. Open Protein Structure "
                    "and View 3D Structure. Partial mapping is shown as partial; "
                    "unmapped residues are not invented."
                )
                st.rerun()
        elif protein_seq:
            st.caption(
                "Motif hits belong to a different protein sequence than the Protein tab. "
                "No 3D highlight is created."
            )
        else:
            st.caption(
                "Analyze the same protein in the Protein tab to highlight this motif on 3D."
            )
    else:
        st.caption(
            "This motif search is not protein. DNA/RNA motif hits are not mapped "
            "onto a protein structure."
        )
    st.download_button(
        "Download motif hits (CSV)",
        data=_export_frame(rows).to_csv(index=False).encode("utf-8"),
        file_name="helixscope_motif_hits.csv",
        mime="text/csv",
        key="motif_hits_download",
    )
    motif_report = provenance.analysis_envelope(
        module="MOTIF",
        payload={
            "molecule": packed.get("molecule"),
            "patterns": packed.get("patterns") or [packed.get("pattern")],
            "n_hits": len(hits),
            "hits": hits,
        },
        status="COMPUTED",
        algorithm="IUPAC/exact motif search",
        parameters={
            "molecule": packed.get("molecule"),
            "max_hits": motif_search.MAX_MOTIF_HITS,
            "max_patterns": motif_search.MAX_MOTIF_PATTERNS,
        },
        source="user sequence",
        input_identifier=str(packed.get("sequence_hash") or ""),
    )
    st.download_button(
        "Download motif report (JSON)",
        data=json.dumps(motif_report, indent=2).encode("utf-8"),
        file_name="helixscope_motif_report.json",
        mime="application/json",
        key="motif_report_json",
    )


def _render_crispr_genome_reference(guides: list, sequence: str) -> None:
    """Painel de referencia genomica verificada e busca Cas-OFFinder.

    Args:
        guides: Guias SpCas9 ja desenhados nesta sessao.
        sequence: DNA de desenho atual.

    Returns:
        None. Escreve a secao na pagina.

    Raises:
        Nenhum.
    """
    if not _section_toggle(
        "Reference Genome",
        key="crispr_genome_panel",
        expanded=True,
    ):
        return
    with st.container(border=True):
        st.markdown(
            "<span style='color:var(--hs-text-secondary);font-size:12px;'>A filename such as "
            "hg38.fa is not GRCh38. READY means the local genomic.fna.gz MD5 "
            "matched NCBI md5checksums.txt, the FASTA was indexed, and the "
            "index sidecar records the FASTA SHA-256. Opening this tab does "
            "not download a genome. Genome-wide is only COMPLETED_FULL_REFERENCE "
            "on a public READY assembly. The HelixScope TEST REFERENCE is never "
            "a public genome. Local workers are not cloud workers. Annotation "
            "is UNKNOWN unless a matching GFF/GTF for the same assembly is "
            "loaded. Feature overlap is not a functional or pathogenic claim.</span>",
            unsafe_allow_html=True,
        )
        rows = genome_store.list_local_references()
        public_rows = [r for r in rows if r.get("kind") != "synthetic_test_reference"]
        test_rows = [r for r in rows if r.get("kind") == "synthetic_test_reference"]
        table = []
        for item in public_rows:
            table.append(
                {
                    "Assembly": item.get("assembly"),
                    "Accession": item.get("accession"),
                    "Organism": item.get("organism"),
                    "Status": item.get("local_status"),
                    "Verified": item.get("ready"),
                    "Size": item.get("estimated_size_label"),
                    "Official MD5": item.get("official_compressed_md5"),
                    "SHA-256": (item.get("file_sha256") or "")[:16] or "N/A",
                    "Source": item.get("source"),
                }
            )
        if table:
            st.dataframe(pd.DataFrame(table), width="stretch", hide_index=True)
        labels = [
            f"{item.get('assembly')} [{item.get('local_status')}]"
            for item in public_rows
        ]
        selected_label = st.selectbox(
            "Public assembly",
            labels or ["(catalog empty)"],
            key="crispr_gw_assembly_label",
        )
        selected = public_rows[labels.index(selected_label)] if labels else None
        if selected is not None:
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Organism: "
                f"{html_escape(selected.get('organism'))}. "
                f"Assembly: {html_escape(selected.get('assembly'))}. "
                f"Accession: {html_escape(selected.get('accession'))}. "
                f"Release: {html_escape(selected.get('release'))}. "
                f"Level: {html_escape(selected.get('assembly_level'))}. "
                f"Official compressed MD5 (NCBI md5checksums.txt "
                f"{html_escape(crispr_assemblies.CATALOG_SNAPSHOT)}): "
                f"{html_escape(selected.get('official_compressed_md5'))}. "
                f"Local status: {html_escape(selected.get('local_status'))}. "
                f"Checksum status: {html_escape(selected.get('checksum_status'))}. "
                f"Estimated size: {html_escape(selected.get('estimated_size_label'))}. "
                f"Download URL: {html_escape(selected.get('download_url'))}. "
                f"Destination: data/references/{html_escape(selected.get('id'))}/ "
                f"(or HELIXSCOPE_REFERENCE_DIR). "
                f"License: {html_escape(selected.get('license_note'))}</span>",
                unsafe_allow_html=True,
            )
            disk = genome_store.disk_status(
                need_bytes=int(selected.get("estimated_compressed_bytes") or 0) * 3
            )
            st.caption(
                f"Free disk: {int(disk['free_bytes']) / (1024 ** 3):.1f} GiB. "
                f"{'Enough for this download estimate.' if disk.get('ok') else disk.get('reason')}"
            )
            compressed = int(selected.get("estimated_compressed_bytes") or 0)
            if compressed:
                preview = resource_admission.admit_reference_search(
                    assembly_id=str(selected.get("id") or ""),
                    fasta_bytes=compressed * 3,
                    engine_available=True,
                )
                st.caption(
                    f"Cas-OFFinder memory gate (estimate from compressed size x3, "
                    f"not a measured FASTA): {preview.get('decision')}. "
                    f"{preview.get('reason')}"
                )
            ident = str(selected.get("id") or "")
            if selected.get("local_status") != genome_store.STATUS_READY:
                if st.button(
                    "Download / configure this assembly",
                    key="crispr_gw_download",
                    help="Explicit user action. Not started by opening CRISPR.",
                ):
                    try:
                        spawned = genome_download.spawn_download_assembly(ident)
                    except genome_download.GenomeDownloadError as exc:
                        st.error(str(exc))
                    else:
                        st.info(
                            f"Download worker pid {spawned.get('pid')} started. "
                            f"Status will move DOWNLOADING to VERIFYING to READY "
                            f"only after NCBI MD5 matches. Estimated "
                            f"{spawned.get('estimated_size_label')}."
                        )
                        st.rerun()
            else:
                st.success(
                    f"READY (verified): {selected.get('assembly')} "
                    f"{selected.get('accession')} SHA-256 "
                    f"{(selected.get('file_sha256') or '')[:16]}..."
                )
        test_item = test_rows[0] if test_rows else None
        if test_item is not None:
            st.caption(
                "TEST REFERENCE (synthetic; never labelled as GRCh38/T2T/GRCm39 "
                f"or genome-wide). Status: {test_item.get('local_status')}."
            )
            st.markdown(status_badge("TEST_ONLY"), unsafe_allow_html=True)
            if test_item.get("local_status") != genome_store.STATUS_READY:
                if st.button(
                    "Install HelixScope TEST REFERENCE",
                    key="crispr_gw_install_test",
                ):
                    try:
                        genome_store.install_test_reference()
                    except genome_store.GenomeStoreError as exc:
                        st.error(str(exc))
                    else:
                        st.rerun()
            else:
                st.caption(
                    "TEST REFERENCE is READY (fixture SHA-256 after copy+index). "
                    "Not a public assembly."
                )

        ready_choices = [
            item
            for item in rows
            if item.get("ready")
        ]
        if not ready_choices:
            st.info(
                "No READY reference. Install the TEST REFERENCE for a small "
                "engine test, or download a public assembly explicitly. "
                "Missing Cas-OFFinder is UNAVAILABLE, not a silent fallback."
            )
            return
        if not guides:
            st.info("Design SpCas9 guides first. Genome search needs a 20 nt spacer.")
            return
        ready_labels = [
            (
                f"{item.get('assembly')} "
                f"{'(TEST REFERENCE, not genome-wide)' if item.get('not_a_public_assembly') else '(public READY)'}"
            )
            for item in ready_choices
        ]
        ready_label = st.selectbox(
            "READY reference to search",
            ready_labels,
            key="crispr_gw_ready_label",
        )
        ready_item = ready_choices[ready_labels.index(ready_label)]
        is_test = bool(ready_item.get("not_a_public_assembly"))
        scope_options = [
            "CURRENT REGION (pasted DNA / local scan only)",
            "PROVIDED REFERENCE (user FASTA / workspace / NCBI record)",
            (
                "TEST REFERENCE (complete fixture; not genome-wide)"
                if is_test
                else "FULL GENOME (READY public assembly; genome-wide only if the job completes)"
            ),
        ]
        st.radio(
            "Search scope (informational; this panel runs only the READY reference)",
            scope_options,
            index=2,
            key="crispr_gw_scope_radio",
        )
        explicit = st.checkbox(
            (
                "Search the TEST REFERENCE with Cas-OFFinder (not genome-wide)"
                if is_test
                else (
                    "Search the complete READY assembly with Cas-OFFinder "
                    "(genome-wide only if COMPLETED_FULL_REFERENCE)"
                )
            ),
            value=False,
            key="crispr_gw_explicit",
        )
        guide_labels = [
            f"#{g.get('rank')} {g.get('guide_sequence')} {g.get('pam_sequence')} "
            f"{g.get('strand')} @{g.get('position')}"
            for g in guides[: min(len(guides), crispr.MAX_REFERENCE_GUIDES)]
        ]
        selected_guide_label = st.selectbox(
            "Guide (SpCas9 20 nt + NGG)",
            guide_labels,
            key="crispr_gw_guide",
        )
        selected_guide = guides[guide_labels.index(selected_guide_label)]
        mm = st.slider(
            "Cas-OFFinder mismatch threshold",
            min_value=0,
            max_value=4,
            value=3,
            key="crispr_gw_mm",
        )
        include_nag = st.checkbox(
            "Also search NAG (in addition to NGG)",
            value=True,
            key="crispr_gw_nag",
        )
        dna_bulge = st.number_input(
            "DNA bulge (Cas-OFFinder; 0 = Hamming spacer, no DNA bulge)",
            min_value=0,
            max_value=2,
            value=0,
            key="crispr_gw_dna_bulge",
        )
        rna_bulge = st.number_input(
            "RNA bulge (Cas-OFFinder; 0 = no RNA bulge)",
            min_value=0,
            max_value=2,
            value=0,
            key="crispr_gw_rna_bulge",
        )
        cof = crispr_casoffinder.detect_cas_offinder()
        live_cas = engine_validation.live_record("Cas-OFFinder")
        live_status = str((live_cas or {}).get("status") or "")
        if live_status == engine_validation.STATUS_LIVE_VALIDATED:
            engine_state = "LIVE_VALIDATED"
        elif cof.get("available"):
            engine_state = "DETECTED"
        else:
            engine_state = "NOT_INSTALLED"
        ocl = opencl_runtime.diagnose()
        st.caption(
            f"Engine: Cas-OFFinder. Status: {engine_state}. "
            f"Detected binary: {bool(cof.get('available'))}. "
            f"Version: {cof.get('version') or 'N/A'}. "
            f"Stable policy: {crispr_casoffinder.CAS_OFFINDER_STABLE_VERSION} "
            f"(Cas-OFFinder 3 not used). Device: "
            f"{(ocl.get('device_selection') or {}).get('device') or ocl.get('preferred_device')}. "
            f"GPU is not required. OpenCL platforms: "
            f"{(ocl.get('platforms') or {}).get('n_platforms')}. "
            f"{'OpenCL ready.' if ocl.get('opencl_ready_for_cas_offinder') else ocl.get('blocker') or 'OpenCL not ready.'} "
            f"Bulges requested: DNA {int(dna_bulge)} RNA {int(rna_bulge)} "
            f"(native bulges UNAVAILABLE on 2.4.1). Threads: 1."
        )
        if not cof.get("available"):
            st.warning(str(cof.get("reason") or "Cas-OFFinder is not installed."))
        elif engine_state != "LIVE_VALIDATED":
            st.warning(
                "Cas-OFFinder requires a working OpenCL runtime. DETECTED is not "
                "LIVE_VALIDATED. FULL GENOME stays locked until a TEST REFERENCE "
                "run completes with re-extracted hits."
            )
        ident = str(ready_item.get("id") or "")
        ready_rec = genome_store.ready_record(ident)
        gate = genome_jobs.full_genome_search_preflight(ident)
        admission = gate.get("resource_admission") or {}
        if admission:
            st.caption(
                f"Resource admission for this READY FASTA: "
                f"{admission.get('decision') or gate.get('resource_status')}. "
                f"{admission.get('reason') or ''}"
            )
        if is_test:
            full_genome_ok = False
            run_allowed = bool(explicit and cof.get("available") and ready_rec)
        else:
            full_genome_ok = bool(gate.get("allowed"))
            run_allowed = bool(explicit and full_genome_ok and ready_rec)
        if not is_test and not full_genome_ok:
            st.info(str(gate.get("reason") or "FULL GENOME is locked."))
        current_identity = ""
        if ready_rec is not None:
            current_identity = genome_jobs.cache_identity(
                guide_hash=provenance.sequence_digest(
                    str(selected_guide.get("guide_sequence") or "")
                ),
                assembly_id=ident,
                reference_sha256=str(ready_rec["manifest"].get("file_sha256") or ""),
                engine="cas_offinder",
                engine_version=str(cof.get("version") or ""),
                max_mismatches=int(mm),
                dna_bulge=int(dna_bulge),
                rna_bulge=int(rna_bulge),
                include_nag=bool(include_nag),
            )
        col_run, col_cancel, col_refresh = st.columns(3)
        with col_run:
            run_clicked = st.button(
                "Start Cas-OFFinder job",
                key="crispr_gw_run",
                disabled=not run_allowed,
            )
        with col_cancel:
            cancel_clicked = st.button("Cancel job", key="crispr_gw_cancel")
        with col_refresh:
            refresh_clicked = st.button("Refresh job status", key="crispr_gw_refresh")
        if run_clicked:
            try:
                status = genome_jobs.submit_casoffinder_job(
                    guide=selected_guide,
                    assembly_id=ident,
                    max_mismatches=int(mm),
                    dna_bulge=int(dna_bulge),
                    rna_bulge=int(rna_bulge),
                    include_nag=bool(include_nag),
                    threads=1,
                    spawn=True,
                )
            except genome_jobs.GenomeJobError as exc:
                st.error(str(exc))
            else:
                st.session_state["crispr_gw_job_id"] = status.get("job_id")
                st.rerun()
        job_id = str(st.session_state.get("crispr_gw_job_id") or "")
        if cancel_clicked and job_id:
            genome_jobs.request_cancel(job_id)
            st.rerun()
        if refresh_clicked:
            st.rerun()
        if not job_id:
            return
        status = genome_jobs.load_status(job_id) or {}
        job = genome_jobs.load_job(job_id) or {}
        result = genome_jobs.load_result(job_id)
        status_name = str(status.get("status") or "")
        stale = bool(
            current_identity
            and job.get("cache_identity")
            and current_identity != job.get("cache_identity")
        )
        st.markdown(
            f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Job "
            f"{html_escape(job_id[:12])}... status "
            f"<strong>{html_escape(status_name)}</strong>. "
            f"{html_escape(status.get('reason') or job.get('progress_note') or '')} "
            f"Progress: {html_escape(status.get('stage') or job.get('stage') or job.get('progress') or status.get('progress') or 'indeterminate')} "
            f"(Cas-OFFinder does not report a percent). "
            f"Cache hit: {html_escape(status.get('cache_hit'))}. "
            f"Stale vs current widgets: {html_escape(stale)}.</span>",
            unsafe_allow_html=True,
        )
        if stale:
            st.warning(
                "Guide, assembly, engine version or parameters changed. "
                "The displayed job is not the current query. Start a new job."
            )
        if status_name == genome_jobs.JOB_RUNNING:
            st.caption("Worker is running. Status is indeterminate until the process exits.")
            time.sleep(1.2)
            st.rerun()
        if result is None:
            return
        genome_wide_flag = bool(result.get("genome_wide"))
        if genome_wide_flag and status_name != genome_jobs.JOB_COMPLETED_FULL_REFERENCE:
            genome_wide_flag = False
        scores = result.get("guide_specificity") or crispr_specificity.scores_from_search(result)
        mit = (scores.get("mit_sguide") or {}) if isinstance(scores, dict) else {}
        cfd = (scores.get("cfd_sguide") or {}) if isinstance(scores, dict) else {}
        assembly_name = result.get("verified_assembly") or result.get("assembly_declared") or ident
        mit_label = (
            f"MIT guide specificity — {assembly_name} genome-wide"
            if mit.get("genome_wide")
            else "MIT guide specificity — this search scope (not genome-wide)"
        )
        cfd_label = (
            f"CFD guide specificity — {assembly_name} genome-wide"
            if cfd.get("genome_wide")
            else "CFD guide specificity — this search scope (not genome-wide)"
        )
        st.caption(
            f"{mit_label}: {mit.get('score') if mit.get('score') is not None else 'N/A'}. "
            f"{cfd_label}: {cfd.get('score') if cfd.get('score') is not None else 'N/A'}. "
            f"Status {status_name}. Truncated={result.get('truncated')}. "
            f"Complete={html_escape((result.get('search_completeness') or {}).get('complete'))}. "
            f"Additional exact sites: {mit.get('n_additional_exact')}. "
            "Partial/timeout/cancel/resource-limit scores stay N/A."
        )
        st.markdown(
            "<p class='module-desc'>Provenance</p>",
            unsafe_allow_html=True,
        )
        st.caption(
            f"Assembly: {html_escape(result.get('verified_assembly') or result.get('assembly_declared'))}. "
            f"Accession: {html_escape(result.get('accession'))}. "
            f"Scope: {html_escape(result.get('scope_label'))}. "
            f"Reference SHA-256: {html_escape((result.get('engine_input_manifest') or {}).get('observed_sha256') or result.get('reference_identity_hash'))}. "
            f"Engine: Cas-OFFinder {(result.get('tool') or {}).get('version')}. "
            f"Device: {(result.get('tool') or {}).get('device') or 'C'}. "
            f"Guide: {html_escape(result.get('guide_sequence'))}. "
            f"Mismatches: {html_escape(result.get('max_mismatches'))}. "
            f"Completion: {html_escape(status_name)}. "
            f"Genome-wide flag: {html_escape(result.get('genome_wide'))}. "
            f"Timestamp: {html_escape(result.get('timestamp'))}. "
            f"Software: {html_escape(result.get('software_version'))}."
        )
        hits = list(result.get("hits") or [])
        annotation_pack = genome_annotation.load_ready_annotation(ident)
        if annotation_pack is None:
            hits = genome_annotation.annotate_hits(
                hits, None, assembly_id=ident, annotation_sha256=""
            )
        else:
            hits = genome_annotation.annotate_hits(
                hits,
                annotation_pack["index"],
                assembly_id=ident,
                annotation_sha256=str(annotation_pack.get("file_sha256") or ""),
            )
        if status_name not in {
            genome_jobs.JOB_COMPLETED_FULL_REFERENCE,
            genome_jobs.JOB_COMPLETED_TEST_REFERENCE,
            crispr_offtarget.JOB_COMPLETED,
        }:
            st.info("Hits are not a complete catalogue for this status.")
            if status_name == genome_jobs.JOB_RESOURCE_LIMIT:
                st.warning(
                    "Hit cap reached. Result is not genome-wide. "
                    "Aggregate specificity is N/A."
                )
            return
        roles = sorted({str(h.get("site_role") or "") for h in hits} | {""})
        pams = sorted({str(h.get("PAM") or "") for h in hits} | {""})
        filter_mm = st.selectbox(
            "Filter mismatches",
            ["all"] + [str(i) for i in range(0, 5)],
            key="crispr_gw_filter_mm",
        )
        filter_role = st.selectbox("Filter role", ["all"] + [r for r in roles if r], key="crispr_gw_filter_role")
        filter_pam = st.selectbox("Filter PAM", ["all"] + [p for p in pams if p], key="crispr_gw_filter_pam")
        filter_feature = st.selectbox(
            "Filter genomic context",
            ["all", "coding", "noncoding", "intergenic", "exact", "off-target"],
            key="crispr_gw_filter_feature",
        )
        sort_by = st.selectbox(
            "Sort",
            ["mismatches", "cfd", "hsu", "coordinate", "feature"],
            key="crispr_gw_sort",
        )
        filtered = []
        for hit in hits:
            if filter_mm != "all" and str(hit.get("mismatches")) != filter_mm:
                continue
            if filter_role != "all" and str(hit.get("site_role") or "") != filter_role:
                continue
            if filter_pam != "all" and str(hit.get("PAM") or "") != filter_pam:
                continue
            ctx = hit.get("genomic_context") or {}
            label = str(hit.get("feature_annotation") or ctx.get("primary_label") or "UNKNOWN")
            if filter_feature == "coding" and not ctx.get("coding"):
                continue
            if filter_feature == "noncoding" and not ctx.get("noncoding"):
                continue
            if filter_feature == "intergenic" and not ctx.get("intergenic") and label != "intergenic":
                continue
            if filter_feature == "exact" and str(hit.get("site_role") or "") not in {
                "on_target",
                "additional_exact",
            }:
                continue
            if filter_feature == "off-target" and str(hit.get("site_role") or "") == "on_target":
                continue
            filtered.append(hit)
        reverse = sort_by in {"cfd", "hsu"}
        def _sort_key(hit: dict):
            if sort_by == "cfd":
                raw = hit.get("cfd_score")
                return float(raw) if isinstance(raw, (int, float)) else -1.0
            if sort_by == "hsu":
                raw = hit.get("score")
                return float(raw) if isinstance(raw, (int, float)) else -1.0
            if sort_by == "coordinate":
                return (str(hit.get("chromosome_or_contig") or ""), int(hit.get("start_0based") or 0))
            if sort_by == "feature":
                return str(hit.get("feature_annotation") or "")
            return int(hit.get("mismatches") or 0)
        filtered.sort(key=_sort_key, reverse=reverse)
        page_size = 25
        n_pages = max(1, (len(filtered) + page_size - 1) // page_size)
        page = int(
            st.number_input(
                "Results page",
                min_value=1,
                max_value=n_pages,
                value=1,
                key="crispr_gw_page",
            )
        )
        start = (page - 1) * page_size
        page_hits = filtered[start : start + page_size]
        st.caption(
            f"{len(filtered)} hit(s) after filters; showing {len(page_hits)} on page "
            f"{page}/{n_pages}. Not rendering the full list at once. "
            "Sorting uses existing metrics; it is not a risk prediction. "
            f"Feature annotation from matching GFF/GTF only; otherwise UNKNOWN."
        )
        table_rows = []
        for hit in page_hits:
            disp = genome_coordinates.internal_to_display(
                int(hit.get("start_0based") or 0),
                int(hit.get("end_0based") or 0),
            )
            ctx = hit.get("genomic_context") or {}
            genes = ctx.get("genes") or []
            gene_txt = ", ".join(
                str(g.get("symbol") or g.get("gene_id") or "")
                for g in genes[:4]
            )
            table_rows.append(
                {
                    "contig": hit.get("chromosome_or_contig"),
                    "start_1based": disp.get("start_1based"),
                    "end_1based": disp.get("end_1based"),
                    "strand": hit.get("strand"),
                    "sequence": hit.get("target_sequence"),
                    "PAM": hit.get("PAM"),
                    "mismatches": hit.get("mismatches"),
                    "CFD": hit.get("cfd_score") if hit.get("cfd_score") is not None else "N/A",
                    "Hsu": hit.get("score") if hit.get("score") is not None else "N/A",
                    "role": hit.get("site_role"),
                    "context": hit.get("feature_annotation") or "UNKNOWN",
                    "gene": gene_txt or "N/A",
                }
            )
        if table_rows:
            st.dataframe(pd.DataFrame(table_rows), width="stretch", hide_index=True)
        if not page_hits:
            return
        detail_labels = [
            f"{h.get('chromosome_or_contig')} {h.get('start_0based')}-{h.get('end_0based')} {h.get('strand')}"
            for h in page_hits
        ]
        detail_label = st.selectbox("Hit detail", detail_labels, key="crispr_gw_detail")
        detail = page_hits[detail_labels.index(detail_label)]
        mm_pos = detail.get("mismatch_positions") or []
        ctx = detail.get("genomic_context") or {}
        gene_bits = []
        for gene in list(ctx.get("genes") or [])[:8]:
            gene_bits.append(
                f"{html_escape(gene.get('symbol') or '')} "
                f"({html_escape(gene.get('gene_id') or '')}, "
                f"{html_escape(gene.get('strand') or '')})"
            )
        st.caption(
            f"Alignment vs spacer {selected_guide.get('guide_sequence')}: "
            f"target {detail.get('target_sequence')} PAM {detail.get('PAM')} "
            f"mismatch positions {mm_pos}. "
            f"Context: {html_escape(detail.get('feature_annotation') or 'UNKNOWN')}. "
            f"Genes: {'; '.join(gene_bits) if gene_bits else 'none'}. "
            f"Transcripts listed (no principal-transcript pick): "
            f"{html_escape(', '.join(str(t) for t in (ctx.get('transcripts') or [])[:8]) or 'none')}. "
            f"{html_escape(ctx.get('effect_note') or 'No functional effect is claimed.')}"
        )
        if ready_rec is not None:
            contig = str(detail.get("chromosome_or_contig") or "")
            fai = genome_fasta.load_fai(ready_rec["fai_path"])
            row = (fai.get("contigs") or {}).get(contig)
            if row:
                start0 = max(0, int(detail.get("start_0based") or 0) - 40)
                end0 = min(int(row["length"]), int(detail.get("end_0based") or 0) + 40)
                try:
                    window = genome_fasta.fetch_sequence(
                        ready_rec["fasta_path"],
                        row,
                        start0,
                        end0,
                        strand="+",
                    )
                except genome_fasta.GenomeFastaError as exc:
                    st.caption(str(exc))
                else:
                    features = []
                    if annotation_pack is not None:
                        features = genome_annotation.overlapping_features(
                            annotation_pack["index"],
                            contig,
                            start0,
                            end0,
                        )
                    try:
                        view = genome_region.local_region_view(
                            contig=contig,
                            contig_length=int(row["length"]),
                            hit=detail,
                            features=features,
                            flank_nt=40,
                            sequence=window,
                        )
                    except genome_region.RegionViewError as exc:
                        st.caption(str(exc))
                    else:
                        axis = view.get("axis") or {}
                        disp_axis = axis.get("display") or {}
                        st.caption(
                            f"Local region {html_escape(contig)} "
                            f"display {disp_axis.get('start_1based')}-{disp_axis.get('end_1based')} "
                            f"(internal {axis.get('start_0based')}-{axis.get('end_0based')}, "
                            "0-based half-open). Not a full-genome track."
                        )
                        st.code(
                            f"{contig}:{start0}-{end0} (internal 0-based half-open, + strand)\n{window}",
                            language="text",
                        )
        st.caption(
            "A genome hit is not a 3D complex. Fetch 4UN3/4OO8 only via the "
            "structure panel when a deposited complex exists. No phylogenetic "
            "claim is made from these coordinates."
        )


def _render_crispr_provided_reference(
    result: dict,
    guides: list,
    sequence: str,
    live_hash: str,
    max_mismatches: int,
) -> None:
    """Painel de busca off-target em referencia fornecida (acao explicita).

    Args:
        result: Pacote crispr_result da sessao.
        guides: Guias avaliados.
        sequence: DNA de desenho atual.
        live_hash: Hash da sequencia colada.

    Returns:
        None. Escreve a secao na pagina.

    Raises:
        Nenhum.
    """
    if _section_toggle(
        "Off-target search in a provided reference",
        key="crispr_reference_panel",
        expanded=True,
    ):
        with st.container(border=True):
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>This search uses "
                "DNA you already have: an uploaded FASTA, workspace DNA, or an "
                "NCBI nucleotide record fetched in this session. HelixScope does "
                "not download hg38/mm39, does not accept an arbitrary genome URL "
                "or filesystem path, and does not call this search genome-wide. "
                "The upload filename is not the assembly. Click Search; opening "
                "this tab does not start a search. This session has no login, "
                "per-user quotas or distributed workers.</span>",
                unsafe_allow_html=True,
            )
            source_choice = st.radio(
                "Reference source",
                [
                    "Upload FASTA",
                    "Workspace DNA",
                    "NCBI nucleotide in this session",
                ],
                key="crispr_ref_source",
            )
            org = st.text_input(
                "Organism (user-declared unless NCBI supplies it)",
                key="crispr_ref_organism",
            )
            assembly = st.text_input(
                "Assembly (user-declared; not inferred from filename)",
                key="crispr_ref_assembly",
            )
            version = st.text_input(
                "Assembly version / release (user-declared)",
                key="crispr_ref_version",
            )
            catalog = crispr_assemblies.lookup_declaration(
                assembly, declared_organism=org
            )
            if catalog.get("catalog_match"):
                match = catalog["catalog_match"]
                st.caption(
                    "Catalog snapshot "
                    f"{catalog.get('catalog_snapshot')}: '{assembly or org}' "
                    f"matches {match.get('assembly')} "
                    f"({match.get('accession')}, {match.get('source')}). "
                    f"Status: {catalog.get('assembly_identity_status')}. "
                    "This does not verify that the FASTA content is that genome. "
                    "Checksum was not retrieved. Nothing was downloaded."
                )
            elif assembly:
                st.caption(
                    "Declared assembly did not match the HelixScope catalog "
                    "snapshot. It remains a user declaration, not a verified "
                    "public assembly."
                )
            uploaded = None
            fasta_text = ""
            source_label = "user FASTA upload"
            accession = ""
            filename = ""
            if source_choice == "Upload FASTA":
                uploaded = st.file_uploader(
                    "Reference FASTA (.fa, .fasta, optional .gz)",
                    type=["fa", "fasta", "fna", "gz"],
                    key="crispr_ref_upload",
                )
                st.caption(
                    "A rejected upload is not analyzed. Use a simple file name such as sample.fasta."
                )
            guide_labels = [
                f"#{g.get('rank')} {g.get('guide_sequence')} {g.get('pam_sequence')} "
                f"{g.get('strand')} @{g.get('position')}"
                for g in guides[: min(len(guides), crispr.MAX_REFERENCE_GUIDES)]
            ]
            if len(guides) > crispr.MAX_REFERENCE_GUIDES:
                st.info(
                    f"Reference search is limited to the top "
                    f"{crispr.MAX_REFERENCE_GUIDES} guides by current ranking "
                    "(resource cap, not a biological rule)."
                )
            selected_label = st.selectbox(
                "Guide to search",
                guide_labels,
                key="crispr_ref_guide",
            )
            selected = guides[guide_labels.index(selected_label)]
            try:
                span = crispr.guide_target_span(selected, sequence)
                st.caption(
                    f"Protospacer span on pasted DNA (internal 0-based sense): "
                    f"{span['start']}-{span['end']}, strand {span['strand']}. "
                    f"Structure mapping: {span['structure_mapping']}. "
                    "This is not a Cas9 3D complex."
                )
            except crispr.CrisprError as exc:
                st.caption(str(exc))
            include_perfect = st.checkbox(
                "Include 0-mismatch sites (perfect matches with NGG/NAG)",
                value=True,
                key="crispr_ref_include_perfect",
            )
            cof = crispr_casoffinder.detect_cas_offinder()
            engine_choices = ["Built-in PAM index"]
            if cof.get("available"):
                engine_choices.append("Cas-OFFinder")
            engine_label = st.radio(
                "Search engine",
                engine_choices,
                key="crispr_ref_engine",
                help=(
                    "Cas-OFFinder appears only when the real binary is detected. "
                    "HelixScope does not reimplement it and does not switch "
                    "silently to the PAM index."
                ),
            )
            engine_id = (
                "cas_offinder"
                if str(engine_label).startswith("Cas-OFFinder")
                else "pam_index"
            )
            if engine_id == "cas_offinder":
                if cof.get("available"):
                    st.caption(
                        "Cas-OFFinder detected: "
                        f"{html_escape(cof.get('version') or 'version not reported')} "
                        f"({html_escape(cof.get('report_path') or 'executable')}). "
                        "Device C (CPU). Bulges not requested. Search is the "
                        "provided FASTA, not genome-wide."
                    )
                else:
                    st.warning(
                        "Cas-OFFinder is not installed. A search with this engine "
                        "will be UNAVAILABLE. The built-in PAM index is not used "
                        "as a silent substitute."
                    )
            search_clicked = st.button(
                "Search selected reference",
                key="crispr_ref_search_button",
            )
            if search_clicked:
                try:
                    if source_choice == "Upload FASTA":
                        if uploaded is None:
                            raise crispr_reference.ReferenceError(
                                "Upload a FASTA before searching. No genome was assumed.",
                                "INVALID_INPUT",
                            )
                        filename = str(getattr(uploaded, "name", "") or "")
                        _reject_nonportable_upload_name(filename)
                        raw_reference = uploaded.getvalue()
                        if not raw_reference:
                            raise crispr_reference.ReferenceError(
                                "Uploaded file is empty. Nothing was analyzed.",
                                "INVALID_INPUT",
                            )
                        fasta_text = crispr_reference.decode_reference_payload(
                            raw_reference
                        )
                        source_label = "user FASTA upload"
                    elif source_choice == "Workspace DNA":
                        workspace = str(st.session_state.get("workspace_dna") or "")
                        meta = st.session_state.get("workspace_dna_meta") or {}
                        loaded = crispr_reference.reference_from_workspace_dna(
                            workspace,
                            identifier=str(meta.get("identifier") or "workspace_dna"),
                            organism=org or str(meta.get("organism") or ""),
                            accession=str(meta.get("accession") or ""),
                            version=version or str(meta.get("version") or ""),
                        )
                        fasta_text = (
                            f">{loaded['contigs'][0]['identifier']}\n"
                            f"{loaded['contigs'][0]['sequence']}\n"
                        )
                        source_label = "workspace DNA"
                        accession = str(meta.get("accession") or "")
                        if not org:
                            org = str(meta.get("organism") or "")
                    else:
                        record = st.session_state.get("ncbi_record")
                        loaded = crispr_reference.reference_from_ncbi_record(
                            record if isinstance(record, dict) else {}
                        )
                        fasta_text = (
                            f">{loaded['contigs'][0]['identifier']}\n"
                            f"{loaded['contigs'][0]['sequence']}\n"
                        )
                        source_label = "NCBI Entrez"
                        accession = str(loaded.get("accession") or "")
                        if not org:
                            org = str(loaded.get("organism_declared") or "")
                        if not version:
                            version = str(loaded.get("version_declared") or "")
                    with st.spinner(
                        "Running Cas-OFFinder on the provided FASTA..."
                        if engine_id == "cas_offinder"
                        else "Indexing PAM sites and searching the provided reference..."
                    ):
                        packed = _reference_offtarget_search(
                            str(selected.get("guide_sequence") or ""),
                            str(selected.get("pam_sequence") or ""),
                            str(selected.get("strand") or "+"),
                            int(selected.get("position") or 0),
                            fasta_text,
                            source_label,
                            org,
                            assembly,
                            version,
                            accession,
                            filename,
                            int(max_mismatches),
                            bool(include_perfect),
                            live_hash,
                            engine_id,
                        )
                    packed["reference_meta"] = {
                        "source": source_label,
                        "organism_declared": org,
                        "assembly_declared": assembly,
                        "version_declared": version,
                        "accession": accession,
                        "upload_filename": filename,
                        "fasta_text": fasta_text,
                    }
                    packed["guide_snapshot"] = {
                        "guide_sequence": selected.get("guide_sequence"),
                        "pam_sequence": selected.get("pam_sequence"),
                        "strand": selected.get("strand"),
                        "position": selected.get("position"),
                        "efficiency_score": selected.get("doench_score"),
                        "efficiency_method": selected.get("efficiency_method"),
                    }
                    packed["target_hash"] = live_hash
                    st.session_state["crispr_reference_search"] = packed
                except crispr_reference.ReferenceError as exc:
                    st.session_state["crispr_reference_search"] = {
                        "status": exc.category,
                        "reason": str(exc),
                        "verified_hit_count": None,
                        "hits": [],
                        "genome_wide": False,
                        "target_hash": live_hash,
                        "truncated": False,
                    }
                    st.markdown(
                        status_badge(exc.category), unsafe_allow_html=True
                    )
                    st.error(str(exc))
                    return
                except (ValueError, TypeError) as exc:
                    st.session_state["crispr_reference_search"] = {
                        "status": "ERROR",
                        "reason": str(exc),
                        "verified_hit_count": None,
                        "hits": [],
                        "genome_wide": False,
                        "target_hash": live_hash,
                    }
                    st.markdown(status_badge("ERROR"), unsafe_allow_html=True)
                    st.error(str(exc))
                    return

            packed = st.session_state.get("crispr_reference_search")
            if not isinstance(packed, dict):
                st.markdown(status_badge("NOT_RUN"), unsafe_allow_html=True)
                st.info(
                    "Reference search has not been run. Hit counts are not shown "
                    "as zero."
                )
                return
            stored_target = str(packed.get("target_hash") or "")
            if stored_target and not _sequence_hashes_match(stored_target, live_hash):
                _hide_stale_result(
                    "Stored reference off-target results belong to a different "
                    "target DNA. Search again."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            if packed.get("max_mismatches") not in {None, max_mismatches} and int(
                packed.get("max_mismatches") or -1
            ) != int(max_mismatches):
                _hide_stale_result(
                    "Mismatch parameter changed. Previous reference hits are not reused."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            if packed.get("include_perfect") not in {None, include_perfect} and bool(
                packed.get("include_perfect")
            ) != bool(include_perfect):
                _hide_stale_result(
                    "Perfect-site inclusion changed. Previous reference hits are not reused."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            expected_algorithm = (
                crispr_casoffinder.ALGORITHM_CAS_OFFINDER
                if engine_id == "cas_offinder"
                else crispr.ALGORITHM_SPCAS9_PAM_INDEX
            )
            if packed.get("algorithm") and packed.get("algorithm") != expected_algorithm:
                _hide_stale_result(
                    "Search engine changed. Previous reference hits are not reused."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            selected_guide_seq = str(selected.get("guide_sequence") or "")
            if (
                packed.get("guide_sequence")
                and packed.get("guide_sequence") != selected_guide_seq
            ):
                _hide_stale_result(
                    "Stored reference hits belong to a different guide. Search again."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            meta = packed.get("reference_meta") or {}
            if assembly and str(meta.get("assembly_declared") or "") not in {
                assembly,
                str(packed.get("assembly_declared") or ""),
            }:
                _hide_stale_result(
                    "Declared assembly changed. Previous off-target results are "
                    "not reused."
                )
                st.session_state.pop("crispr_reference_search", None)
                return
            status = str(packed.get("status") or "NOT_RUN")
            st.markdown(status_badge(status), unsafe_allow_html=True)
            if packed.get("reason"):
                st.caption(str(packed["reason"]))
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>"
                f"Scope: {html_escape(packed.get('scope_label') or packed.get('scope') or 'provided_reference')}. "
                f"Algorithm: {html_escape(packed.get('algorithm') or '')}. "
                f"Organism declared: {html_escape(packed.get('organism_declared') or 'none')}. "
                f"Assembly declared: {html_escape(packed.get('assembly_declared') or 'none')}. "
                f"Verified assembly: {html_escape(packed.get('verified_assembly') or 'none')}. "
                f"Assembly identity: {html_escape(packed.get('assembly_identity_status') or 'n/a')}. "
                f"Version declared: {html_escape(packed.get('version_declared') or 'none')}. "
                f"Reference hash: {html_escape((packed.get('reference_identity_hash') or '')[:16])}. "
                f"Index hash: {html_escape((packed.get('index_hash') or '')[:16])}. "
                f"Genome-wide: {html_escape(packed.get('genome_wide'))}. "
                f"Workers: in-process only."
                "</span>",
                unsafe_allow_html=True,
            )
            snap = packed.get("guide_snapshot") or packed
            pam_text = str(snap.get("pam_sequence") or packed.get("pam_sequence") or "")
            strand_text = str(snap.get("strand") or packed.get("strand") or "")
            verified = packed.get("verified_hit_count")
            render_metric_grid(
                [
                    (
                        "Guide",
                        str(snap.get("guide_sequence") or packed.get("guide_sequence") or "")[:12]
                        or "Unavailable",
                        "",
                    ),
                    ("PAM", pam_text or "Unavailable", ""),
                    ("Strand", strand_text or "Unavailable", ""),
                    (
                        "Verified hits",
                        "not reported" if verified is None else str(verified),
                        "",
                    ),
                ]
            )
            if packed.get("truncated"):
                st.warning(
                    "Result is truncated. It is not a complete off-target catalogue."
                )
            disclosure = packed.get("algorithm_disclosure") or {}
            if disclosure:
                st.caption(
                    "Mismatch model: "
                    f"{disclosure.get('mismatch_model')} | PAM model: "
                    f"{disclosure.get('pam_model')} | Seed: "
                    f"{disclosure.get('seed_model')} | Bulges: "
                    f"{disclosure.get('bulge_support')} | Score: "
                    f"{disclosure.get('score_model')}"
                )
            spec = packed.get("guide_specificity") or {}
            mit = spec.get("mit_sguide") or {}
            cfd_g = spec.get("cfd_sguide") or {}
            spec_items: list[tuple[str, str, str]] = []
            if mit.get("score") is not None:
                spec_items.append(
                    ("MIT Sguide (search scope)", f"{mit.get('score'):.2f}", "")
                )
            if cfd_g.get("score") is not None:
                spec_items.append(
                    ("CFD specificity (search scope)", f"{cfd_g.get('score'):.2f}", "")
                )
            if spec_items:
                render_metric_grid(spec_items)
            if mit.get("score") is None:
                st.caption(
                    "MIT Sguide: not computed "
                    f"({mit.get('status') or packed.get('status')}). "
                    f"{mit.get('method') or 'CRISPOR calcMitGuideScore'}. "
                    "Scale 0-100. This is not 0 and not genome-wide MIT."
                )
            else:
                st.caption(
                    f"Method: {mit.get('method')} | Model: {mit.get('model')} | "
                    f"Version: {mit.get('version')} | Scale: {mit.get('scale')} | "
                    f"CRISPOR integer: {mit.get('score_crispor_integer')} | "
                    f"Sites included: {mit.get('n_sites_included')} | "
                    f"Intended excluded: {mit.get('n_intended_excluded')} | "
                    f"Additional exact: {mit.get('n_additional_exact')} | "
                    f"{mit.get('reason')}"
                )
            if cfd_g.get("score") is None:
                st.caption(
                    "CFD guide specificity: not computed "
                    f"({cfd_g.get('status') or packed.get('status')}). "
                    "Per-site CFD remains 0-1 when hits exist. This is not 0."
                )
            else:
                st.caption(
                    f"Method: {cfd_g.get('method')} | Model: {cfd_g.get('model')} | "
                    f"Version: {cfd_g.get('version')} | Scale: {cfd_g.get('scale')} | "
                    f"Per-hit CFD scale: {cfd_g.get('per_hit_cfd_scale')} | "
                    f"{cfd_g.get('reason')}"
                )
            if status == "COMPLETED" and packed.get("verified_hit_count") == 0:
                st.info(
                    "0 verified hits under selected method/reference. "
                    "This is not a safety claim."
                )
            hits = list(packed.get("hits") or [])
            show_hits = status in {"COMPLETED", "RESOURCE_LIMIT"} and hits
            if show_hits:
                if status == "RESOURCE_LIMIT":
                    st.warning(
                        "These hits are a truncated prefix. verified_hit_count "
                        "is not reported as a complete total."
                    )
                summary = crispr_offtarget.mismatch_summary_from_hits(
                    hits, _int_or_default(packed.get("max_mismatches"), 3)
                )
                if sum(summary.values()) > 0:
                    st.plotly_chart(
                        charts.off_target_mismatch_histogram(summary),
                        width="stretch",
                    )
                contig_ids = {str(h.get("chromosome_or_contig") or "") for h in hits}
                if len(contig_ids) > 1:
                    st.plotly_chart(
                        charts.off_target_contig_histogram(hits),
                        width="stretch",
                    )
                pam_counts = {}
                for hit in hits:
                    pam_name = str(hit.get("pam_class") or hit.get("PAM") or "")
                    if pam_name:
                        pam_counts[pam_name] = pam_counts.get(pam_name, 0) + 1
                if pam_counts:
                    st.plotly_chart(
                        charts.off_target_pam_histogram(pam_counts),
                        width="stretch",
                    )
                hsu_values = [
                    float(hit["score"])
                    for hit in hits
                    if hit.get("score") is not None
                ]
                if hsu_values:
                    st.plotly_chart(
                        charts.off_target_score_histogram(
                            hsu_values,
                            title="Hsu 2013 single-hit scores (0-100) of verified sites",
                            x_title="Hsu single-hit score",
                        ),
                        width="stretch",
                    )
                page_size = 25
                n_pages = max(1, (len(hits) + page_size - 1) // page_size)
                page = st.number_input(
                    "Hit page",
                    min_value=1,
                    max_value=n_pages,
                    value=1,
                    key="crispr_ref_hit_page",
                )
                start = (int(page) - 1) * page_size
                page_hits = hits[start : start + page_size]
                rows = []
                for hit in page_hits:
                    rows.append(
                        {
                            "Contig": hit.get("chromosome_or_contig"),
                            "Position": hit.get("position"),
                            "Strand": hit.get("strand"),
                            "Target": hit.get("target_sequence"),
                            "PAM": hit.get("PAM"),
                            "PAM class": hit.get("pam_class"),
                            "Mismatches": hit.get("mismatches"),
                            "Mismatch positions": ",".join(
                                str(p) for p in hit.get("mismatch_positions") or []
                            ),
                            "Site role": hit.get("site_role"),
                            "Hsu hit score": hit.get("score"),
                            "Hsu method": hit.get("score_method"),
                            "Hsu model": hit.get("score_model"),
                            "Hsu scale": "0-100",
                            "CFD": hit.get("cfd_score"),
                            "CFD method": hit.get("cfd_method"),
                            "CFD scale": hit.get("cfd_scale"),
                        }
                    )
                st.caption(
                    f"Showing hits {start + 1}-{start + len(page_hits)} of "
                    f"{len(hits)} verified hits on this page window. "
                    "Not a genome map."
                )
                st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
                if status == "COMPLETED":
                    try:
                        export_rows = crispr.export_offtarget_rows(packed)
                    except crispr.CrisprError as exc:
                        st.info(str(exc))
                    else:
                        st.download_button(
                            "Download verified hits (CSV)",
                            data=_export_frame(export_rows).to_csv(index=False).encode("utf-8"),
                            file_name=(
                                f"{safe_download_filename(result.get('name'), 'crispr')}"
                                "_offtargets.csv"
                            ),
                            mime="text/csv",
                            key="crispr_ref_hits_csv",
                        )
                else:
                    st.info(
                        "Truncated hits are not exported as a complete catalogue."
                    )
            elif status in {"TIMEOUT", "ERROR", "INVALID_INPUT", "INCOMPLETE"}:
                st.info(
                    "No hit catalogue is shown as complete for this status. "
                    "Counts are not converted to zero hits."
                )
            record = st.session_state.get("ncbi_record")
            if (
                status == "COMPLETED"
                and hits
                and isinstance(record, dict)
                and packed.get("scope") == crispr_reference.SEARCH_SCOPE_NCBI_RECORD
            ):
                features = ncbi_fetch.extract_annotated_features(record)
                contig_len = int((record.get("length") or 0) or len(str(record.get("sequence") or "")))
                contexts = []
                for hit in hits[:25]:
                    overlap = crispr_offtarget.overlapping_ncbi_features(
                        hit, features, contig_length=contig_len
                    )
                    for item in overlap:
                        contexts.append(
                            {
                                "Contig": hit.get("chromosome_or_contig"),
                                "Hit position": hit.get("position"),
                                "Feature": item.get("type"),
                                "Gene": item.get("gene") or "-",
                                "Overlap nt": item.get("overlap_nt"),
                                "Note": item.get("claim"),
                            }
                        )
                if contexts:
                    st.caption(
                        "NCBI feature overlap on this record. Overlap is not "
                        "gene disruption."
                    )
                    st.dataframe(
                        pd.DataFrame(contexts), width="stretch", hide_index=True
                    )


def render_crispr() -> None:
    """Renderiza a aba de desenho de guias CRISPR.

    Coleta a sequencia alvo e parametros, encontra e avalia guias com metricas por
    guia, declara explicitamente quais scores reais nao estao disponiveis, ordena
    por qualidade combinada, exibe tabelas, graficos, adaptacoes de editores de
    base e prime, analise de off-target local e desenho de primers de validacao.

    Args:
        Nenhum.

    Returns:
        None. Escreve componentes diretamente na pagina.

    Raises:
        Nenhum.
    """
    in_col, opt_col = st.columns([2, 1])
    with in_col:
        workspace_dna = str(st.session_state.get("workspace_dna") or "")
        draft_seq = str(st.session_state.get("crispr_draft_seq") or "")
        draft_name = str(st.session_state.get("crispr_draft_name") or "")
        if not str(st.session_state.get("crispr_text") or "").strip():
            if draft_seq:
                st.session_state["crispr_text"] = draft_seq
            elif workspace_dna:
                st.session_state["crispr_text"] = workspace_dna
        if not str(st.session_state.get("crispr_name") or "").strip() and draft_name:
            st.session_state["crispr_name"] = draft_name
        sequence_text = st.text_area(
            "Target DNA (SpCas9)",
            height=220,
            key="crispr_text",
            placeholder="Paste target DNA (A, T, C, G)",
            help=(
                "DNA with A, T, C, G. SpCas9 needs NGG and at least 23 bp. "
                "This field starts empty. Sequence is not inserted automatically."
            ),
        )
    with opt_col:
        sequence_name = st.text_input(
            "Sequence Name", placeholder="Optional name", key="crispr_name"
        )
        st.markdown(badge("SpCas9  NGG  20 nt spacer", "dna"), unsafe_allow_html=True)
        max_mismatches = st.slider(
            "Max Off-target Mismatches", min_value=1, max_value=4, value=3, key="crispr_mm"
        )
        run_off_target = st.checkbox(
            "Run Off-target Analysis (local scan)",
            value=True,
            key="crispr_offtarget",
        )
    st.session_state["crispr_draft_seq"] = str(sequence_text or "")
    st.session_state["crispr_draft_name"] = str(sequence_name or "")
    st.markdown(
        "<span style='color:var(--hs-text-secondary);'>This tab is SpCas9 only (the densest and most "
        "used DNA nuclease here). Other Cas enzymes remain in the Python module "
        "but are not offered in the UI. Off-target scan is local, requires an "
        f"adjacent NGG or NAG, and is refused above {crispr.MAX_OFF_TARGET_SCAN_NT:,} "
        "nt. It is not genome-wide. Reference FASTA search is a separate explicit "
        "action. The Reference Genome panel downloads only after an explicit "
        "button; opening this tab does not start a genome download.</span>",
        unsafe_allow_html=True,
    )
    pdb_ids = ", ".join(
        item["pdb_id"] for item in crispr_structure_catalog.list_experimental_complexes()
    )
    st.caption(
        "HelixScope does not fabricate a Cas protein or a DNA-Cas complex from a "
        "guide. Catalog PDB identifiers "
        f"({pdb_ids}) can be fetched from the RCSB Files API by explicit action "
        "below. Guide and PAM highlights require an exact subsequence match on "
        "the deposited polymer. A protospacer span on the pasted DNA is a "
        "sequence interval, not a structure. Off-target genomic coordinates are "
        "not 3D."
    )
    packed = st.session_state.get("crispr_result") or {}
    _render_crispr_complex_section(
        str(packed.get("sequence") or sequence_text),
        list(packed.get("guides") or []),
    )

    if st.button("Design Guides", key="crispr_button"):
        cas_system = "SpCas9"
        st.session_state.pop("crispr_complex", None)
        info = _validate_for(sequence_text, "DNA")
        if not info["is_valid"]:
            st.session_state.pop("crispr_result", None)
            st.error(
                info["rejection_reason"]
                or "Provide a valid DNA sequence (A, T, C, G) for SpCas9."
            )
            diagnosis = crispr.diagnose_spcas9_target(sequence_text)
            for message in diagnosis["messages"]:
                st.warning(message)
            return
        diagnosis = crispr.diagnose_spcas9_target(info["sequence"])
        if len(info["sequence"]) > crispr.MAX_OFF_TARGET_SCAN_NT:
            run_off_target = False
            st.info(
                "Local off-target scan skipped: sequence longer than "
                f"{crispr.MAX_OFF_TARGET_SCAN_NT:,} nt."
            )
        try:
            with st.spinner("Scanning SpCas9 NGG sites and evaluating guides..."):
                guides_ranked = _evaluate_guides(
                    info["sequence"], cas_system, max_mismatches, run_off_target
                )
        except ValueError as exc:
            st.error(str(exc))
            return
        st.session_state["crispr_result"] = {
            "sequence": info["sequence"],
            "sequence_hash": provenance.sequence_digest(info["sequence"]),
            "name": sequence_name or "Target",
            "cas_system": cas_system,
            "max_mismatches": max_mismatches,
            "run_off_target": run_off_target,
            "guides": guides_ranked,
            "diagnosis": diagnosis,
        }

    result = st.session_state.get("crispr_result")
    if result and result.get("cas_system") != "SpCas9":
        st.session_state.pop("crispr_result", None)
        st.session_state.pop("crispr_reference_search", None)
        result = None
    live_crispr = _validate_for(sequence_text, "DNA")
    live_hash = (
        provenance.sequence_digest(str(live_crispr["sequence"]))
        if live_crispr.get("is_valid")
        else ""
    )
    stored_hash = str((result or {}).get("sequence_hash") or "")
    if not stored_hash and result:
        stored_hash = provenance.sequence_digest(str(result.get("sequence") or ""))
    if result and not _sequence_hashes_match(stored_hash, live_hash):
        _hide_stale_result(
            "Stored CRISPR guides belong to a different DNA sequence. "
            "Click Design Guides to refresh."
        )
        st.session_state.pop("crispr_reference_search", None)
        result = None
    if not result:
        _render_crispr_model_transparency()
        _render_crispr_genome_reference([], "")
        st.info(
            "Paste DNA with at least one NGG (20 nt spacer + PAM). Click Design "
            "Guides. Target DNA and sequence name start empty. The published "
            "EMX1 oligo remains in tests and documentation only; it is not "
            "inserted automatically. The Reference Genome panel above does not "
            "download assemblies until you click Download."
        )
        return

    guides = result["guides"]
    cas_system = result["cas_system"]
    info = crispr.system_info(cas_system)
    sequence = result["sequence"]
    diagnosis = result.get("diagnosis") or crispr.diagnose_spcas9_target(sequence)
    st.markdown(status_badge("PREDICTED"), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.result_header_html(
            "SpCas9 guide design",
            status="PREDICTED",
            lines=[
                f"{len(guides)} guides",
                f"PAM {info.get('pam_display') or 'NGG'} · local off-target "
                f"{'on' if result.get('run_off_target') else 'off'}",
            ],
        ),
        unsafe_allow_html=True,
    )
    helix_workspace.render_explanation(
        "crispr",
        {
            "status": "PREDICTED",
            "n_guides": len(guides),
            "pam": info.get("pam_display") or "NGG",
            "reference_scope": result.get("reference_id") or result.get("reference_scope") or "",
            "test_reference": bool(result.get("test_reference") or result.get("not_a_public_assembly")),
            "method": "SpCas9 guide scan",
            "source": "HelixScope CRISPR",
        },
    )
    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>Computational prediction. "
        "Requires experimental validation. Not experimentally confirmed.</span>",
        unsafe_allow_html=True,
    )

    render_metric_grid(
        [
            ("Length", str(diagnosis.get("length", len(sequence))), "bp"),
            ("NGG both strands", str(diagnosis.get("ngg_total", 0)), ""),
            ("Complete sites", str(diagnosis.get("guides_possible", 0)), ""),
            ("NAG (not designed)", str(diagnosis.get("nag_total", 0)), ""),
        ]
    )
    st.markdown(
        "<span style='color:var(--hs-text-secondary);font-size:12px;'>NGG counts are motif counts "
        "on both strands. Complete sites are NGG with 20 nt upstream. NAG is a "
        "known reduced-activity PAM (Hsu et al., 2013) and is not used to design "
        "guides here.</span>",
        unsafe_allow_html=True,
    )

    if not guides:
        _render_crispr_model_transparency()
        st.warning(
            "No SpCas9 guide RNAs found: every complete site needs 20 nt + NGG "
            "on the same strand."
        )
        for message in diagnosis.get("messages") or []:
            st.info(message)
        return

    best = guides[0]
    valid_guides = [g for g in guides if _guide_quality_flags(g)["valid_gc"] and not _guide_quality_flags(g)["poly_t"]]
    high_quality = [g for g in guides if _guide_quality_flags(g)["high_quality"]]

    render_metric_grid(
        [
            ("Guides Found", str(len(guides)), ""),
            ("Valid Guides", str(len(valid_guides)), ""),
            ("Highest heuristic efficiency", f"{best['doench_score']:.4f}", ""),
        ]
    )
    st.markdown(badge(info["label"], "dna"), unsafe_allow_html=True)

    st.markdown(
        f"<span style='color:var(--hs-text-secondary);'>System: {html_escape(info['label'])} - PAM "
        f"{html_escape(info['pam_display'])} - target "
        f"{html_escape(info['target_molecule'])}. {html_escape(info['notes'])}"
        "</span>",
        unsafe_allow_html=True,
    )

    _render_crispr_model_transparency()

    ready = best["guide_sequence"]
    pam_note = f" + PAM {best['pam_sequence']}" if best.get("pam_sequence") else ""
    st.markdown(
        badge(
            f"Highest-scoring spacer under the selected heuristic: {ready}{pam_note}",
            "prot",
        ),
        unsafe_allow_html=True,
    )
    order = best.get("order") or crispr.order_ready_oligos(best)
    st.markdown(
        f"<span style='color:var(--hs-text-secondary);'>RNA alphabet (T to U, spacer only, no "
        f"tracrRNA scaffold): {html_escape(order.get('spacer_rna', ''))}</span>",
        unsafe_allow_html=True,
    )
    mh = best.get("microhomology") or {}
    if mh.get("length"):
        st.info(
            f"Observed microhomology at the estimated cut: {mh['length']} bp "
            f"({mh.get('motif')}). This is a sequence observation, not an indel "
            "or frameshift prediction."
        )
    seed = best.get("seed") or {}
    if seed.get("count", 0) > 1:
        st.warning(
            f"The 12 nt seed {seed.get('seed')} occurs {seed['count']} times in "
            "the pasted sequence (both strands). That is not a genome-wide count."
        )

    if info["editor"] == "base_editor":
        window = crispr.base_editing_window(best, cas_system)
        st.info(
            f"Base editor {window['editor']}: edit window at protospacer positions "
            f"{window['window'][0]}-{window['window'][1]}. {window['product']}."
        )
    elif info["editor"] == "prime_editor":
        pe_notes = crispr.prime_editing_notes(best, cas_system)
        st.info(
            f"Prime editor: nick at position {pe_notes['nick_position']} of the "
            f"pasted sequence. {pe_notes['reason']} Recommended design tools: "
            f"{', '.join(pe_notes['recommended_tools'])}."
        )
    elif info["target_molecule"] == "RNA":
        st.info(
            "Cas13 targets single-stranded RNA and shows collateral cleavage after "
            "activation, so there is no site-specific DNA cut and no cut-site map."
        )

    display_count = min(len(guides), 30)
    guide_filter = str(
        st.text_input("Filter guides by sequence or PAM", key="crispr_guide_filter") or ""
    ).strip().upper()
    visible_guides = guides
    if guide_filter:
        visible_guides = [
            item
            for item in guides
            if guide_filter in str(item.get("guide_sequence") or "").upper()
            or guide_filter in str(item.get("pam_sequence") or "").upper()
            or guide_filter in str(item.get("strand") or "").upper()
        ]
        if not visible_guides:
            st.info("No guide matches that filter. Ranked list is unchanged.")
            visible_guides = guides
    display_count = min(len(visible_guides), 30)
    if _section_toggle(
        f"Ranked Guides (showing {display_count} of {len(visible_guides)} / {len(guides)} scanned)",
        key="crispr_ranked",
        expanded=True,
    ):
        with st.container(border=True):
            if high_quality:
                st.markdown(
                    badge(f"{len(high_quality)} high-quality guide(s)", "dna"),
                    unsafe_allow_html=True,
                )
            rows = []
            for guide in visible_guides[:display_count]:
                flags = _guide_quality_flags(guide)
                proxy = guide.get("specificity_proxy")
                composite = guide.get("composite_score")
                risk = guide.get("risk")
                rows.append(
                    {
                        "Rank": guide.get("rank"),
                        "Guide": guide.get("guide_sequence"),
                        "PAM": guide.get("pam_sequence") or "-",
                        "Pos": guide.get("position"),
                        "Strand": guide.get("strand"),
                        "GC%": guide.get("gc_content"),
                        "Efficiency heuristic": guide.get("doench_score"),
                        "Efficiency method": guide.get(
                            "efficiency_method", crispr.HEURISTIC_EFFICIENCY_METHOD
                        ),
                        "Doench RS2": (
                            "N/A"
                            if guide.get("ruleset2_score") is None
                            else guide.get("ruleset2_score")
                        ),
                        "DeepHF": (
                            "N/A"
                            if guide.get("deephf_score") is None
                            else guide.get("deephf_score")
                        ),
                        "Poly-T": "yes" if flags["poly_t"] else "no",
                        "Hairpin": flags["hairpin_risk"],
                        "Seed copies": (guide.get("seed") or {}).get("count", "-"),
                        "MH bp": (guide.get("microhomology") or {}).get("length", "-"),
                        "RE sites": ", ".join(guide.get("restriction_sites") or []) or "-",
                        "Specificity": "-" if proxy is None else f"{proxy:.4f}",
                        "Specificity method": (
                            "-"
                            if proxy is None
                            else guide.get("specificity_method") or "-"
                        ),
                        "Local class": risk["level"] if risk else "-",
                        "Composite": "-" if composite is None else f"{composite:.4f}",
                        "High quality": "yes" if flags["high_quality"] else "no",
                    }
                )
            st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
            st.markdown(
                "<span style='color:var(--hs-text-secondary);font-size:12px;'>Efficiency heuristic is a "
                "positional score inspired by Rule Set 2, not the published "
                "model. Doench RS2 and DeepHF are N/A unless those official models "
                "are loaded; absence is never shown as 0. Specificity is a local proxy from the pasted sequence, not "
                "CFD or MIT. Seed copies count exact 12 nt PAM-proximal matches on "
                "both strands of the pasted DNA (1 means unique here). MH bp is "
                "observed microhomology at the cut, not an indel predictor. RE "
                "sites are exact matches of EcoRI/BamHI/HindIII/NcoI/NdeI/XhoI in "
                "the 23 nt target. Seed copies are computed for the top "
                f"{crispr.MAX_GUIDES_DETAILED} guides by heuristic efficiency; "
                "blank seed in the CSV means not calculated, not zero copies. "
                "Composite ranks by both only when the off-target "
                "scan is on; otherwise ranking uses efficiency alone.</span>",
                unsafe_allow_html=True,
            )
            report = crispr.generate_report(guides, result["name"])
            st.download_button(
                "Download guide table (CSV)",
                data=report.to_csv().encode("utf-8"),
                file_name=(
                    f"{safe_download_filename(result['name'], 'crispr')}"
                    "_crispr_guides.csv"
                ),
                mime="text/csv",
                key="crispr_download",
            )
            sci_report = crispr.scientific_crispr_report(
                target_sequence=sequence,
                guide=best,
                search=st.session_state.get("crispr_reference_search"),
            )
            st.download_button(
                "Download scientific report (JSON)",
                data=json.dumps(sci_report, indent=2).encode("utf-8"),
                file_name=(
                    f"{safe_download_filename(result['name'], 'crispr')}"
                    "_crispr_report.json"
                ),
                mime="application/json",
                key="crispr_scientific_report",
            )

    if _section_toggle("Guide Comparison and Ranking Charts", key="crispr_charts"):
        with st.container(border=True):
            st.plotly_chart(
                charts.guide_comparison_chart(guides[: min(len(guides), 8)]),
                width="stretch",
            )
            if result["run_off_target"]:
                try:
                    st.plotly_chart(
                        charts.guide_ranking_scatter(guides), width="stretch"
                    )
                except ValueError as exc:
                    st.info(str(exc))
            else:
                st.info(
                    "Enable the off-target scan to plot efficiency against local "
                    "specificity."
                )

    if info["target_molecule"] != "RNA" and _section_toggle(
        "Cut Site and PAM Location", key="crispr_cutmap"
    ):
        with st.container(border=True):
            st.plotly_chart(
                charts.cut_site_map(len(sequence), best, cas_system),
                width="stretch",
            )
            if best.get("cut_site") is not None:
                st.markdown(
                    f"<span style='color:var(--hs-text-secondary);'>Estimated cut site at position "
                    f"{best['cut_site']} of the pasted sequence. Position relative "
                    "to a gene requires an annotation not provided here.</span>",
                    unsafe_allow_html=True,
                )

    if _section_toggle("Guide Sequence Viewer", key="crispr_guide_view"):
        with st.container(border=True):
            viewer_seq = best["guide_sequence"] + best.get("pam_sequence", "")
            st.markdown(
                sequence_display(viewer_seq.replace("U", "T"), "DNA"),
                unsafe_allow_html=True,
            )
            hairpin = best.get("self_complementarity", {})
            poly = best.get("poly_t", {})
            render_metric_grid(
                [
                    ("GC Content", f"{best['gc_content']:.2f}", "%"),
                    ("Heuristic efficiency", f"{best['doench_score']:.4f}", ""),
                    ("Hairpin Stem", str(hairpin.get("max_stem", 0)), "bp"),
                    ("Max poly-T", str(poly.get("max_run", 0)), ""),
                ]
            )
            if float(best["gc_content"]) < 40.0 or float(best["gc_content"]) > 80.0:
                st.warning(
                    f"GC of {best['gc_content']:.2f}% is outside the 40-80% range "
                    "usually recommended for reliable activity."
                )
            if poly.get("has_signal"):
                st.warning(
                    f"Poly-T run of {poly.get('max_run')} may terminate Pol III "
                    "transcription and truncate the guide."
                )
            if hairpin.get("risk") == "high":
                st.warning(
                    f"Self-complementary stem of {hairpin.get('max_stem')} bp may "
                    "form a hairpin that competes with sgRNA scaffold folding."
                )

    if info["editor"] == "base_editor" and _section_toggle(
        "Base Editing Window", key="crispr_be_window"
    ):
        with st.container(border=True):
            window = crispr.base_editing_window(best, cas_system)
            st.markdown(
                f"<span style='color:var(--hs-text-secondary);'>Editor {window['editor']}, window at "
                f"protospacer positions {window['window'][0]}-{window['window'][1]}. "
                f"{window['product']}.</span>",
                unsafe_allow_html=True,
            )
            if not window["editable"]:
                st.warning(
                    "The best guide has no editable base in the window; a different "
                    "guide is needed to make this edit."
                )

    if result["run_off_target"] and info["target_molecule"] != "RNA":
        if _section_toggle("Off-target Analysis (local)", key="crispr_offtarget_view", expanded=True):
            with st.container(border=True):
                summary = best.get("off_target_summary")
                risk = best.get("risk")
                if risk:
                    st.markdown(
                        badge(
                            f"Local class (pasted sequence): {risk['level']}",
                            _risk_variant(risk["level"]),
                        ),
                        unsafe_allow_html=True,
                    )
                    st.markdown(
                        f"<span style='color:var(--hs-text-secondary);font-size:12px;'>Criterion: "
                        f"{html_escape(risk['criterion'])} This is not a safety "
                        "claim and is not genome-wide.</span>",
                        unsafe_allow_html=True,
                    )
                if summary and sum(summary.values()) > 0:
                    st.plotly_chart(
                        charts.off_target_mismatch_histogram(summary),
                        width="stretch",
                    )
                with st.spinner("Listing local off-target sites..."):
                    off_targets = _off_targets(
                        best["guide_sequence"],
                        sequence,
                        result["max_mismatches"],
                        str(best.get("pam_sequence") or ""),
                    )
                if not off_targets:
                    st.info(
                        "0 verified hits under selected method/reference "
                        "(pasted sequence, adjacent NGG or NAG, mismatches at "
                        "or below the slider). This is not a safety claim and "
                        "not 'no risk'."
                    )
                else:
                    off_rows = [
                        {
                            "Position": off["position"],
                            "Strand": off["strand"],
                            "Off-target Sequence": off["off_target_sequence"],
                            "PAM": off.get("pam_sequence") or "-",
                            "PAM class": off.get("pam_class") or "-",
                            "Mismatches": off["mismatches"],
                            "Mismatch positions": ",".join(
                                str(p) for p in off.get("mismatch_positions") or []
                            ),
                            "Local weighted score": off["risk_score"],
                            "Local score method": off.get("risk_score_method")
                            or "local weighted mismatch heuristic",
                            "Hsu hit score": (
                                "-"
                                if off.get("hsu_hit_score") is None
                                else off.get("hsu_hit_score")
                            ),
                            "Hsu method": off.get("hsu_hit_score_method") or "-",
                        }
                        for off in off_targets
                    ]
                    st.dataframe(
                        pd.DataFrame(off_rows), width="stretch", hide_index=True
                    )
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);'>Scope: pasted sequence only. "
                    "A window is counted only if it is followed by NGG or NAG "
                    "on that strand. This is not genome-wide. BLAST similarity "
                    "is not CRISPR specificity. Hsu single-hit is not the MIT "
                    "Specificity Score.</span>",
                    unsafe_allow_html=True,
                )

    _render_crispr_genome_reference(guides, sequence)
    _render_crispr_provided_reference(
        result, guides, sequence, live_hash, _int_or_default(result.get("max_mismatches"), 3)
    )

    if info["target_molecule"] != "RNA" and _section_toggle(
        "Primer Design for Validation", key="crispr_primers"
    ):
        with st.container(border=True):
            try:
                primers = _primer_pair(best["guide_sequence"], sequence)
            except ValueError as exc:
                st.warning(str(exc))
            else:
                render_metric_grid(
                    [
                        ("Forward Tm", f"{primers['forward_tm']:.2f}", "C"),
                        ("Reverse Tm", f"{primers['reverse_tm']:.2f}", "C"),
                        ("Amplicon Size", str(primers["amplicon_size"]), "bp"),
                    ]
                )
                st.markdown(sequence_display(primers["forward_primer"], "DNA"), unsafe_allow_html=True)
                st.markdown(sequence_display(primers["reverse_primer"], "DNA"), unsafe_allow_html=True)
                st.markdown(
                    "<span style='color:var(--hs-text-secondary);'>These primers are designed for "
                    "Sanger sequencing validation of CRISPR editing efficiency via "
                    "TIDE or ICE analysis.</span>",
                    unsafe_allow_html=True,
                )


def render_variant_explorer() -> None:
    """Renderiza a aba Variant Explorer.

    Recebe a identidade de uma variante (assembly obrigatoria), normaliza-a
    localmente e, apenas quando o utilizador pede explicitamente, consulta o
    Ensembl VEP, o NCBI ClinVar e o InterPro/UniProtKB. Cada camada mostra o seu
    proprio estado; a falha de uma nao apaga as outras.

    Args:
        Nenhum.

    Returns:
        None. Escreve diretamente na interface Streamlit.

    Raises:
        Nenhum. Erros de entrada e de servico sao apresentados como estado.

    Nota biologica:
        HelixScope nao prediz efeito de variante. As consequencias moleculares
        sao termos Sequence Ontology recuperados do Ensembl VEP, as afirmacoes
        clinicas sao submissoes recuperadas do ClinVar e os dominios sao
        entradas recuperadas do InterPro.
    """
    left, right = st.columns([3, 2])
    with left:
        text = st.text_input(
            "Variant",
            key="variant_input",
            placeholder="17 43093557 C G",
            help=(
                "Supported: VCF-like '17 43093557 C G', SPDI "
                "'NC_000017.11:43093556:C:G', HGVS g. substitution "
                "'NC_000017.11:g.43093557C>G', rsID, or ClinVar VCV accession. "
                "HGVS c./p. is not parsed locally; it is sent to the Ensembl VEP "
                "HGVS endpoint."
            ),
        )
    with right:
        assembly = st.selectbox(
            "Assembly (required)",
            options=["GRCh38.p14", "GRCh37"],
            key="variant_assembly",
            help=(
                "Selects both the coordinate system and the Ensembl REST host. "
                "GRCh38 uses rest.ensembl.org and GRCh37 uses "
                "grch37.rest.ensembl.org."
            ),
        )

    with st.expander("Remote services and what leaves this machine", expanded=False):
        disclosures = {
            "Ensembl VEP": ensembl_vep.network_disclosure(assembly),
            "NCBI ClinVar": clinvar_evidence.network_disclosure(),
            "InterPro / UniProtKB": protein_domains.network_disclosure(),
        }
        for name, disclosure in disclosures.items():
            host = disclosure.get("host") or ", ".join(
                str(service.get("host")) for service in disclosure.get("services", [])
            )
            sent = "; ".join(str(item) for item in disclosure.get("data_sent", []))
            st.markdown(
                ui_components.glass_card(
                    f"{name} ({host})",
                    f"<b>Sent:</b> {ui_components.html_escape(sent)}<br>"
                    f"<b>Rate limit policy:</b> "
                    f"{ui_components.html_escape(str(disclosure.get('rate_limit_policy') or ''))}",
                    "default",
                ),
                unsafe_allow_html=True,
            )
        st.caption(
            "Nothing is sent until you press Annotate. Disabling a source below "
            "means its answer is UNKNOWN, not absent."
        )

    col_a, col_b, col_c, col_d = st.columns(4)
    with col_a:
        use_vep = st.toggle("Ensembl VEP", value=True, key="variant_use_vep")
    with col_b:
        use_clinvar = st.toggle("NCBI ClinVar", value=True, key="variant_use_clinvar")
    with col_c:
        use_domains = st.toggle("InterPro / UniProt", value=True, key="variant_use_domains")
    with col_d:
        email = st.text_input(
            "Contact e-mail (NCBI policy)",
            key="variant_email",
            help=(
                "NCBI requires a contact address for E-utilities. Without it the "
                "ClinVar layer is skipped, which is not evidence of absence."
            ),
        )

    if not st.button("Annotate variant", key="variant_button", type="primary"):
        stored = st.session_state.get("variant_result")
        if isinstance(stored, Mapping):
            _render_variant_result(stored)
        return

    try:
        variant = variant_core.build_variant(
            text=text, assembly=assembly, source=variant_core.SOURCE_USER
        )
    except variant_core.VariantInputError as exc:
        st.markdown(ui_components.status_badge(exc.category), unsafe_allow_html=True)
        st.error(str(exc))
        return

    with st.spinner("Retrieving annotation from official services"):
        result = variant_explorer.explore_variant(
            variant=variant,
            email=email,
            enable_vep=bool(use_vep),
            enable_clinvar=bool(use_clinvar),
            enable_domains=bool(use_domains),
        )

    st.session_state["variant_result"] = result
    _render_variant_result(result)


def _render_variant_result(result: Mapping[str, Any]) -> None:
    """Apresenta o resultado do Variant Explorer camada a camada (uso interno).

    Args:
        result: Dict devolvido por variant_explorer.explore_variant.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    variant = result["variant"]
    layers = result["layers"]

    loc = ""
    if variant.get("contig") and variant.get("position_1based") is not None:
        loc = (
            f"{variant.get('contig')}:{variant.get('position_1based')} "
            f"{variant.get('ref') or ''}>{variant.get('alt') or ''}"
        )
    cons = layers.get("consequences") or {}
    cons_data = cons.get("data") if isinstance(cons.get("data"), Mapping) else {}
    protein = layers.get("protein_mapping") or {}
    protein_data = protein.get("data") if isinstance(protein.get("data"), Mapping) else {}
    transcripts = list(protein_data.get("transcripts") or [])
    gene = ""
    hgvsp = ""
    if transcripts:
        gene = str(transcripts[0].get("gene_symbol") or "")
        hgvsp = str(transcripts[0].get("hgvsp") or transcripts[0].get("amino_acids") or "")
    journey_status_genome = "COMPUTED" if loc else "UNMAPPED"
    journey_status_tx = "AVAILABLE" if cons.get("status") == "AVAILABLE" else str(cons.get("status") or "UNAVAILABLE")
    journey_status_prot = (
        "AVAILABLE" if protein.get("status") == "AVAILABLE" else str(protein.get("status") or "UNMAPPED")
    )
    st.markdown(
        helix_workspace.variant_journey_html(
            (
                ("Genome", f"{variant.get('assembly') or ''} {loc}".strip(), journey_status_genome),
                ("Gene", gene or (cons_data.get("most_severe_consequence") or ""), journey_status_tx),
                ("Transcript", str(cons_data.get("transcript_count") or "Unavailable"), journey_status_tx),
                ("Protein", hgvsp, journey_status_prot),
            )
        ),
        unsafe_allow_html=True,
    )
    clinvar_layer = layers.get("clinvar") or {}
    clinvar_data = clinvar_layer.get("data") if isinstance(clinvar_layer.get("data"), Mapping) else {}
    helix_workspace.render_explanation(
        "variant",
        {
            "status": str(result.get("status") or cons.get("status") or "RETRIEVED"),
            "hgvs": loc,
            "assembly": variant.get("assembly"),
            "most_severe_consequence": cons_data.get("most_severe_consequence"),
            "clinvar_significance": clinvar_data.get("clinical_significance")
            or clinvar_data.get("description")
            or clinvar_layer.get("status"),
            "method": "Ensembl VEP / ClinVar retrieval",
            "source": "Variant Explorer",
        },
    )

    st.markdown("#### Variant identity")
    st.markdown(
        ui_components.badge_row(
            [str(layers[name]["status"]) for name in result["layer_order"]]
        ),
        unsafe_allow_html=True,
    )
    st.caption(
        " | ".join(
            f"{name.replace('_', ' ')}: {layers[name]['status']}"
            for name in result["layer_order"]
        )
    )
    identity_rows = [
        ("Assembly", variant["assembly"] or "N/A"),
        ("Assembly accession", variant["accession"] or "N/A"),
        ("Input format", variant["format"]),
        ("Contig", variant["contig"] or "N/A"),
        (
            "Position (1-based)",
            "N/A" if variant["position_1based"] is None else str(variant["position_1based"]),
        ),
        ("Reference allele", variant["ref"] or "(empty, pure insertion)"),
        ("Alternate allele", variant["alt"] or "(empty, pure deletion)"),
        ("Class", variant["variant_kind"] or "N/A"),
        ("Normalization", variant["normalization"]),
        ("REF validation", variant["ref_validation"]),
        ("Identity hash", (variant["identity_hash"] or "N/A")[:32]),
        ("Effect computed by HelixScope", "None. HelixScope does not predict effect."),
    ]
    st.markdown(ui_components.meta_grid(identity_rows), unsafe_allow_html=True)
    st.caption(variant.get("normalization_method") or variant.get("ref_validation_reason") or "")
    if variant["ref_validation"] == variant_core.REF_NOT_CHECKED:
        st.info(variant["ref_validation_reason"])

    consequences = layers["consequences"]
    st.markdown("#### Molecular consequences (Ensembl VEP)")
    if consequences["status"] != "AVAILABLE":
        st.markdown(ui_components.status_badge(consequences["status"]), unsafe_allow_html=True)
        st.warning(consequences["reason"])
    else:
        data = consequences["data"]
        st.markdown(
            ui_components.meta_grid(
                [
                    ("Source", data["source"]),
                    ("Assembly reported by VEP", data["assembly_name"] or "N/A"),
                    ("Allele string", data["allele_string"] or "N/A"),
                    ("Most severe consequence", data["most_severe_consequence"] or "N/A"),
                    ("Transcripts returned", str(data["transcript_count"])),
                    ("Vocabulary", data["consequence_vocabulary"]),
                ]
            ),
            unsafe_allow_html=True,
        )
        st.caption(data["transcript_note"])
        rows = data["transcript_consequences"]
        if rows:
            frame = pd.DataFrame(
                [
                    {
                        "transcript": row["transcript_id"],
                        "gene": row["gene_symbol"],
                        "biotype": row["biotype"],
                        "canonical": row["canonical"],
                        "consequence_terms (Sequence Ontology)": ", ".join(
                            row["consequence_terms"]
                        ),
                        "impact": row["impact"],
                        "protein_position": row["protein_start"],
                        "amino_acids": row["amino_acids"],
                        "codons": row["codons"],
                        "HGVSc": row["hgvsc"],
                        "HGVSp": row["hgvsp"],
                    }
                    for row in rows
                ]
            )
            st.dataframe(frame, width="stretch", hide_index=True)
            tx_ids = [str(row.get("transcript_id") or "") for row in rows if row.get("transcript_id")]
            if tx_ids:
                st.selectbox(
                    "Selected transcript",
                    tx_ids,
                    key="variant_selected_transcript",
                    help="Dependent protein panels use this transcript. HelixScope does not rank a canonical choice.",
                )
            st.caption(
                "Every transcript returned by VEP is listed. Consequence terms "
                "are Sequence Ontology terms shown exactly as the service "
                "returned them."
            )
        if data["colocated_variants"]:
            st.caption(
                "Colocated variants reported by Ensembl: "
                + ", ".join(
                    f"{item['id']} ({', '.join(item['clin_sig']) or 'no clinical significance reported'})"
                    for item in data["colocated_variants"]
                )
            )

    st.markdown("#### Protein and residue context")
    protein = layers["protein_mapping"]
    if protein["status"] != "AVAILABLE":
        st.markdown(ui_components.status_badge(protein["status"]), unsafe_allow_html=True)
        st.info(protein["reason"])
    else:
        choices = list(protein["data"]["transcripts"])
        selected_id = str(st.session_state.get("variant_selected_transcript") or "")
        primary = choices[0]
        for item in choices:
            if str(item.get("transcript_id") or "") == selected_id:
                primary = item
                break
        st.caption("Selected transcript (explicit; not a hidden canonical ranking).")
        st.markdown(
            ui_components.meta_grid(
                [
                    ("Transcript shown", primary["transcript_id"]),
                    ("Canonical", "yes" if primary["canonical"] else "no"),
                    ("Protein", primary["protein_id"] or "N/A"),
                    ("Protein position", str(primary["protein_start"])),
                    ("Codon change", primary["codons"] or "N/A"),
                    ("Amino acid change", primary["amino_acids"] or "N/A"),
                    ("Transcripts with a protein position", str(protein["data"]["count"])),
                ]
            ),
            unsafe_allow_html=True,
        )
        st.caption(
            "This is the first transcript in the table above, not a hidden "
            "canonical choice. The full list stays visible."
        )
        residue = layers["residue_properties"]
        if residue["status"] == "AVAILABLE":
            change = residue["data"]
            st.markdown(ui_components.status_badge("COMPUTED"), unsafe_allow_html=True)
            st.markdown(
                ui_components.meta_grid(
                    [
                        (
                            "Hydropathy (Kyte-Doolittle)",
                            f"{change['hydropathy_reference']} -> "
                            f"{change['hydropathy_variant']} "
                            f"(delta {change['hydropathy_difference']})",
                        ),
                        (
                            "Charge class",
                            f"{change['charge_class_reference']} -> "
                            f"{change['charge_class_variant']}",
                        ),
                        (
                            "Side-chain category",
                            f"{change['side_chain_category_reference']} -> "
                            f"{change['side_chain_category_variant']}",
                        ),
                        (
                            "Residue volume (A^3)",
                            f"{change['volume_reference_a3']} -> "
                            f"{change['volume_variant_a3']} "
                            f"(delta {change['volume_difference_a3']})",
                        ),
                        ("Stability claim", "None. No stability model is integrated."),
                    ]
                ),
                unsafe_allow_html=True,
            )
            st.caption(change["disclaimer"])
        else:
            st.info(residue["reason"])

    st.markdown("#### Domain context (InterPro) and protein reference (UniProtKB)")
    reference = layers["protein_reference"]
    if reference["status"] == "AVAILABLE":
        entry = reference["data"]
        st.markdown(
            ui_components.meta_grid(
                [
                    ("UniProt accession", entry["accession"]),
                    ("Entry name", entry["entry_name"]),
                    ("Protein name", entry["protein_name"] or "N/A"),
                    ("Organism", entry["organism"] or "N/A"),
                    ("Sequence length", str(entry["sequence_length"] or "N/A")),
                    ("Reviewed (Swiss-Prot)", "yes" if entry["reviewed"] else "no"),
                ]
            ),
            unsafe_allow_html=True,
        )
    else:
        st.markdown(ui_components.status_badge(reference["status"]), unsafe_allow_html=True)
        st.info(reference["reason"])

    domains = layers["domains"]
    if domains["status"] == "AVAILABLE":
        payload = domains["data"]
        coverage = payload.get("coverage") or {}
        covering = coverage.get("covering") or []
        if covering:
            frame = pd.DataFrame(
                [
                    {
                        "entry": item["entry_accession"],
                        "name": item["name"],
                        "type": item["type"],
                        "source_database": item["source_database"],
                        "range(s) covering residue": ", ".join(
                            f"{fragment['start_1based']}-{fragment['end_1based']}"
                            for fragment in item["matching_fragments"]
                        ),
                        "member databases": ", ".join(
                            f"{member['source_database']}:{member['accession']}"
                            for member in item["member_databases"]
                        ),
                    }
                    for item in covering
                ]
            )
            st.dataframe(frame, width="stretch", hide_index=True)
        st.caption(
            str(coverage.get("note") or payload.get("message") or "")
            + " "
            + str(payload.get("scientific_limit") or "")
        )
    else:
        st.markdown(ui_components.status_badge(domains["status"]), unsafe_allow_html=True)
        st.info(domains["reason"])

    st.markdown("#### Clinical evidence retrieved from ClinVar")
    clinical = layers["clinical_evidence"]
    if clinical["status"] not in {"AVAILABLE", "NOT_FOUND"} or not clinical["data"]:
        st.markdown(ui_components.status_badge(clinical["status"]), unsafe_allow_html=True)
        st.info(clinical["reason"])
    else:
        data = clinical["data"]
        if data["status"] == "NOT_FOUND":
            st.markdown(ui_components.status_badge("NOT_FOUND"), unsafe_allow_html=True)
            st.info(data["message"])
        else:
            st.caption(data["identity_note"])
            records = data["records"]
            frame = pd.DataFrame(
                [
                    {
                        "accession": record["accession_version"],
                        "title": record["title"],
                        "germline classification (submitters)": ", ".join(
                            record["germline_descriptions"]
                        )
                        or "not provided",
                        "review status": ", ".join(
                            block["review_status"] for block in record["classifications"]
                        ),
                        "conflict": record["has_conflict"],
                        "position matches query": record.get("identity_match"),
                        "allele verified": record.get("allele_match"),
                        "submitted records (SCV)": record["submission_count"],
                    }
                    for record in records
                ]
            )
            st.dataframe(frame, width="stretch", hide_index=True)
            st.caption(
                "ClinVar classification is the submitted text in the table. "
                "HelixScope does not color benign as green or pathogenic as red, "
                "and does not treat color as a diagnosis."
            )
            if data.get("allele_unverified_count"):
                st.warning(
                    "ClinVar's summary does not publish reference/alternate "
                    "alleles for these locations, so allele identity is "
                    "UNVERIFIED. Different variants at the same base can carry "
                    "different classifications, so do not read a record as this "
                    "variant's classification on position alone."
                )
            if data["has_conflict"]:
                st.warning(
                    "At least one record has conflicting submitted "
                    "classifications. HelixScope shows the conflict and does not "
                    "resolve it."
                )
        st.caption(clinvar_evidence.CLINICAL_DISCLAIMER)

    st.markdown("#### Cross-source evidence")
    cross = result["cross_source"]
    if cross["comparable"]:
        st.markdown(
            ui_components.meta_grid(
                [(source, ", ".join(symbols)) for source, symbols in cross["sources"].items()]
            ),
            unsafe_allow_html=True,
        )
    st.caption(f"{cross['note']} {cross['policy']}")
    st.caption(result["no_cascade_failure"])


def _render_evidence_table(items: list) -> None:
    """Tabela de evidencias com proveniencia (uso interno)."""
    if not items:
        st.caption("No evidence items in this comparison.")
        return
    rows = []
    for item in items:
        if not isinstance(item, Mapping):
            continue
        rows.append(
            {
                "Field": str(item.get("field") or ""),
                "Value": provenance.csv_cell(item.get("value")),
                "Source": str(item.get("source") or ""),
                "Status": str(item.get("evidence_status") or ""),
                "Mapping": str(item.get("mapping_status") or ""),
                "Identifier": str(item.get("identifier") or ""),
                "Retrieved": str(item.get("retrieved_at_utc") or ""),
                "Version": str(item.get("version") or ""),
            }
        )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)


def render_compare_workspace() -> None:
    """Aba Compare: estruturas, variantes, proteinas, guias, evolucao, evidencia.

    Um unico tab com modos. Nao e um redesign visual. RMSD/TM-score so aparecem
    quando a RCSB Alignment API (ou parser) os fornece. Superposicao usa uma
    copia transformada; coordenadas originais nao sao sobrescritas.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum. Erros de servico aparecem como estado.

    Nota biologica:
        Semelhanca estrutural nao e homologia nem funcao. Conservacao de MSA
        nao e essencialidade. Variantes nao sao ordenadas por patogenicidade.
        MIT/CFD so comparam no mesmo ambito. Sem GRCh38 COMPLETED_FULL_REFERENCE
        nao existe comparacao genome-wide.
    """
    na_status = usalign.nucleic_structural_comparison_availability()
    st.caption(
        f"DNA/RNA structural comparison is {na_status.get('status')}: "
        f"{na_status.get('reason')}. Block RMSD and global RMSD stay distinct "
        "fields when the RCSB Alignment API returns both."
    )
    mode = st.radio(
        "Compare",
        ["Structures", "Variants", "Proteins", "Guides", "Evolution", "Evidence"],
        horizontal=True,
        key="compare_mode",
    )
    if mode == "Structures":
        _render_compare_structures()
    elif mode == "Variants":
        _render_compare_variants()
    elif mode == "Proteins":
        _render_compare_proteins()
    elif mode == "Guides":
        _render_compare_guides()
    elif mode == "Evolution":
        _render_compare_evolution()
    else:
        _render_compare_evidence()


def _render_compare_structures() -> None:
    """Modo estruturas: RCSB Alignment API + superposicao (uso interno)."""
    methods = rcsb_alignment.supported_methods()
    labels = [str(row["display_name"]) for row in methods]
    disclosure = rcsb_alignment.network_disclosure()
    with st.expander("Remote service: what leaves this machine", expanded=True):
        st.caption(
            f"Host {disclosure['host']}. Sent: "
            + "; ".join(str(item) for item in disclosure["data_sent"])
            + ". Not sent: "
            + "; ".join(str(item) for item in disclosure["not_sent"])
            + f" {disclosure['rate_limit_policy']}"
        )
        st.caption(str(disclosure.get("engine_location") or ""))
    cols = st.columns(4)
    with cols[0]:
        ref_id = st.text_input(
            "Reference structure",
            value="8HSK",
            key="cmp_ref_id",
            help="PDB ID (experimental) or AF_ / AF-P01308-F1 (predicted CSM).",
        )
    with cols[1]:
        ref_chain = st.text_input("Reference chain", value="A", key="cmp_ref_chain")
    with cols[2]:
        tgt_id = st.text_input(
            "Target structure",
            value="8HSF",
            key="cmp_tgt_id",
            help="PDB, AF_ CSM, or AF-UniProt-F1. Predicted stays PREDICTED.",
        )
    with cols[3]:
        tgt_chain = st.text_input("Target chain", value="A", key="cmp_tgt_chain")
    us_info = usalign.detect_usalign()
    backend_options = ["RCSB Alignment API"]
    if us_info.get("available"):
        backend_options.append("US-align")
    backend = st.radio(
        "Alignment backend",
        backend_options,
        horizontal=True,
        key="cmp_backend",
        help=(
            "RCSB Alignment API is remote. US-align is a local executable. "
            "They are never merged into one RMSD or TM-score."
        ),
    )
    method_label = st.selectbox("Alignment method", labels, key="cmp_method")
    api_name = next(row["api_name"] for row in methods if row["display_name"] == method_label)
    kind = rcsb_alignment.method_kind(api_name)
    st.caption(
        f"API method {api_name} ({kind}). Rigid methods apply one transform to a "
        "coordinate copy. Flexible multi-block results keep the table; the 3D "
        "overlay is PARTIAL (block 0) and is not a single rigid superposition."
    )
    fetch_coords = st.checkbox(
        "Download coordinates for superposition (allowlisted RCSB/AlphaFold files)",
        value=True,
        key="cmp_fetch_coords",
    )
    us_mol = "prot"
    if backend == "US-align":
        us_mol = st.selectbox(
            "US-align molecule type",
            ["prot", "RNA", "auto"],
            index=0,
            key="cmp_usalign_mol",
            help="RNA includes DNA per upstream US-align help. auto may mix polymers.",
        )
        st.caption(
            f"LOCAL US-align {us_info.get('version') or 'version not reported'}. "
            "This does not call the RCSB Alignment API. 3D overlay remains "
            "UNAVAILABLE on this backend because US-align rotates structure 1 "
            "onto structure 2, which is not the RCSB target-transform path."
        )
    if not st.button("Align structures", key="cmp_align_btn", type="primary"):
        stored = st.session_state.get("compare_structures")
        if isinstance(stored, dict):
            _show_structure_comparison(stored)
        return
    try:
        ref_parsed = None
        tgt_parsed = None
        need_coords = bool(fetch_coords) or backend == "US-align"
        if need_coords:
            with st.spinner("Retrieving allowlisted coordinates"):
                ref_load = comparative.load_parsed_structure_for_entry(ref_id)
                tgt_load = comparative.load_parsed_structure_for_entry(tgt_id)
            if ref_load.get("status") == "AVAILABLE":
                ref_parsed = ref_load["parsed"]
            else:
                st.caption(
                    f"Reference coordinates {ref_load.get('status')}: {ref_load.get('reason')}"
                )
            if tgt_load.get("status") == "AVAILABLE":
                tgt_parsed = tgt_load["parsed"]
            else:
                st.caption(
                    f"Target coordinates {tgt_load.get('status')}: {tgt_load.get('reason')}"
                )
        if backend == "US-align":
            if ref_parsed is None or tgt_parsed is None:
                raise comparative.ComparativeError(
                    "US-align requires loaded coordinates for both structures.",
                    "UNAVAILABLE",
                )
            with st.spinner("Running local US-align"):
                result = comparative.compare_structures_usalign(
                    reference_parsed=ref_parsed,
                    target_parsed=tgt_parsed,
                    reference_entry=ref_id,
                    target_entry=tgt_id,
                    reference_chain=ref_chain,
                    target_chain=tgt_chain,
                    mol=str(us_mol),
                    record_validation=True,
                )
        else:
            with st.spinner("Submitting pairwise alignment to RCSB"):
                result = comparative.compare_structures_remote(
                    reference_entry=ref_id,
                    target_entry=tgt_id,
                    reference_chain=ref_chain,
                    target_chain=tgt_chain,
                    method=api_name,
                    reference_parsed=ref_parsed,
                    target_parsed=tgt_parsed,
                )
            engine_validation.record_remote_validation(
                "RCSB Alignment API",
                ok=True,
                version="api/v1",
                details={
                    "reference": ref_id,
                    "target": tgt_id,
                    "method": api_name,
                    "n_pairs": (result.get("alignment") or {}).get("n_aligned_residue_pairs"),
                    "rmsd_global": (result.get("alignment") or {}).get("rmsd_global_angstrom"),
                },
            )
        st.session_state["compare_structures"] = result
        st.session_state["compare_export"] = result.get("export")
        _show_structure_comparison(result)
    except (
        rcsb_alignment.AlignmentApiError,
        structure_alignment.StructureAlignmentError,
        comparative.ComparativeError,
        usalign.USAlignError,
    ) as exc:
        category = getattr(exc, "category", "ERROR")
        st.markdown(status_badge(category), unsafe_allow_html=True)
        st.error(str(exc))


def _show_structure_comparison(result: Mapping[str, Any]) -> None:
    """Mostra metricas, mapeamento e overlay (uso interno)."""
    alignment_rec = result.get("alignment") or {}
    st.markdown(status_badge(str(alignment_rec.get("status") or "RETRIEVED")), unsafe_allow_html=True)
    st.markdown(
        helix_workspace.compare_hero_html(alignment_rec),
        unsafe_allow_html=True,
    )
    n_pairs = alignment_rec.get("n_aligned_residue_pairs")
    rmsd_b = alignment_rec.get("rmsd_block0_angstrom")
    atoms_fitted = alignment_rec.get("atoms_fitted")
    st.markdown(
        f"Aligned residue pairs: {helix_workspace.unavailable_text(n_pairs)}. "
        f"Block RMSD: {helix_workspace.unavailable_text(None if rmsd_b is None else f'{rmsd_b} A')}. "
        f"Atoms fitted: {helix_workspace.unavailable_text(atoms_fitted)}."
    )
    helix_workspace.render_explanation("compare", alignment_rec)
    helix_workspace.render_explanation("compare_rmsd_block", alignment_rec)
    st.caption(
        str(alignment_rec.get("rmsd_note") or "")
        + " Block RMSD and global RMSD are distinct. TM-score is not a similarity percent."
    )
    ref = alignment_rec.get("reference") or {}
    tgt = alignment_rec.get("target") or {}
    st.markdown(
        badge_row(
            [
                str(ref.get("kind_label") or "EXPERIMENTAL"),
                str(tgt.get("kind_label") or "EXPERIMENTAL"),
                str(alignment_rec.get("method_kind") or ""),
                str(alignment_rec.get("superposition_visual") or ""),
            ]
        ),
        unsafe_allow_html=True,
    )
    st.markdown(
        meta_grid(structure_alignment.comparison_summary_rows(alignment_rec)),
        unsafe_allow_html=True,
    )
    st.caption(str(alignment_rec.get("rmsd_note") or ""))
    st.caption(str(alignment_rec.get("tm_score_note") or ""))
    st.caption(str(alignment_rec.get("coverage_note") or ""))
    st.caption(str(alignment_rec.get("atoms_fitted_note") or ""))
    st.caption(str(alignment_rec.get("disclaimer") or ""))
    check = result.get("rmsd_agreement")
    if isinstance(check, Mapping):
        if check.get("comparable") and check.get("agree"):
            st.caption(
                f"Independent CA RMSD after transform: "
                f"{check.get('local_rmsd_angstrom'):.3f} A vs API "
                f"{check.get('api_rmsd_angstrom')} A (delta "
                f"{check.get('delta_angstrom'):.3f} A; "
                f"tolerance {check.get('abs_tolerance_angstrom')} A abs or "
                f"{100 * float(check.get('rel_tolerance') or 0):.0f}% rel)."
            )
        elif check.get("comparable"):
            st.warning(
                "Independent RMSD disagrees with the API beyond documented "
                f"tolerance (local {check.get('local_rmsd_angstrom')}, API "
                f"{check.get('api_rmsd_angstrom')}). Superposition is still a "
                "transformed copy; the disagreement is shown, not hidden."
            )
        else:
            st.caption(str(check.get("reason") or "Independent RMSD not comparable."))
    bundle = result.get("superposition")
    if isinstance(bundle, Mapping) and bundle.get("visual_status") == "UNAVAILABLE":
        st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
        st.caption(str(bundle.get("note") or "Superposition visual unavailable."))
    elif isinstance(bundle, Mapping):
        if bundle.get("visual_status") == "PARTIAL":
            st.markdown(status_badge("PARTIAL"), unsafe_allow_html=True)
            st.caption(str(bundle.get("note") or ""))
        try:
            view_mode = st.radio(
                "Viewport",
                ["Overlay", "A", "B"],
                horizontal=True,
                key="cmp_view_mode",
                help="Visibility only. Overlay shows reference plus the transformed target copy.",
            )
            mode_map = {"Overlay": "overlay", "A": "reference", "B": "target"}
            figure = structure_viewer.figure_from_superposition(
                bundle,
                view_mode=mode_map.get(str(view_mode), "overlay"),
            )
        except ValueError:
            st.markdown(status_badge("UNAVAILABLE"), unsafe_allow_html=True)
            st.caption("Not enough finite CA coordinates to draw the overlay.")
        else:
            render_hand_control(
                plot_key="cmp_superposition",
                structure_id=str(
                    alignment_rec.get("alignment_id")
                    or alignment_rec.get("job_id")
                    or "compare"
                ),
                illustrative=False,
                control_key="cmp",
            )
            event = _structure_plotly_chart(
                figure,
                key="cmp_superposition",
                on_select="rerun",
                height=640,
            )
            picked = structure_viewer.parse_superposition_selection(event)
            if picked and picked.get("label_seq_id") is not None:
                lookup = structure_alignment.lookup_correspondence(
                    alignment_rec,
                    side=str(picked.get("side") or "reference"),
                    asym_id=str(picked.get("asym_id") or ""),
                    label_seq_id=int(picked["label_seq_id"]),
                )
                st.markdown(status_badge(str(lookup.get("status"))), unsafe_allow_html=True)
                st.caption(
                    f"Clicked {lookup.get('clicked')} -> "
                    f"{lookup.get('corresponding') or 'UNMAPPED'}. "
                    f"{lookup.get('reason') or ''}"
                )
    else:
        st.caption(
            "Superposition requires downloaded coordinates. Alignment metrics "
            "above are from the remote API and do not need local atoms."
        )
    look_cols = st.columns(3)
    with look_cols[0]:
        side = st.selectbox("Lookup side", ["reference", "target"], key="cmp_lookup_side")
    with look_cols[1]:
        chain = st.text_input("Lookup chain", value=str(ref.get("asym_id") or "A"), key="cmp_lookup_chain")
    with look_cols[2]:
        seq = st.number_input("label_seq_id (1-based)", min_value=1, value=2, step=1, key="cmp_lookup_seq")
    looked = structure_alignment.lookup_correspondence(
        alignment_rec,
        side=str(side),
        asym_id=str(chain),
        label_seq_id=int(seq),
    )
    st.markdown(status_badge(str(looked.get("status"))), unsafe_allow_html=True)
    st.caption(
        f"{looked.get('clicked')} corresponds to {looked.get('corresponding') or 'UNMAPPED'}. "
        "Nearest-residue mapping is not used."
    )
    export = result.get("export") or {}
    st.download_button(
        "Download comparison JSON",
        data=json.dumps(export, indent=2, default=str),
        file_name=safe_download_filename(
            f"helixscope_structure_compare_{ref.get('entry_id')}_{tgt.get('entry_id')}.json"
        ),
        mime="application/json",
        key="cmp_struct_dl",
    )


def _render_compare_variants() -> None:
    """Modo variantes lado a lado (uso interno)."""
    st.caption(
        "Two Variant Explorer results, same assembly recommended. HelixScope "
        "does not rank which variant is worse. Transcript ambiguity is preserved."
    )
    left, right = st.columns(2)
    with left:
        text_a = st.text_input("Variant A", value="17 43093557 C G", key="cmp_var_a")
    with right:
        text_b = st.text_input("Variant B", value="17 43094500 T C", key="cmp_var_b")
    assembly = st.selectbox(
        "Assembly (required)",
        options=["GRCh38.p14", "GRCh37"],
        key="cmp_var_assembly",
    )
    email = st.text_input("Contact e-mail (NCBI ClinVar policy)", key="cmp_var_email")
    if not st.button("Compare variants", key="cmp_var_btn", type="primary"):
        stored = st.session_state.get("compare_variants")
        if isinstance(stored, dict):
            _show_variant_comparison(stored)
        return
    try:
        va = variant_core.build_variant(text=text_a, assembly=assembly, source=variant_core.SOURCE_USER)
        vb = variant_core.build_variant(text=text_b, assembly=assembly, source=variant_core.SOURCE_USER)
    except variant_core.VariantInputError as exc:
        st.markdown(status_badge(exc.category), unsafe_allow_html=True)
        st.error(str(exc))
        return
    with st.spinner("Retrieving both variants from official services"):
        ra = variant_explorer.explore_variant(variant=va, email=email)
        rb = variant_explorer.explore_variant(variant=vb, email=email)
        compared = comparative.compare_variants(ra, rb)
    st.session_state["compare_variants"] = compared
    st.session_state["compare_export"] = compared.get("export")
    _show_variant_comparison(compared)


def _show_variant_comparison(compared: Mapping[str, Any]) -> None:
    """Cartões A/B sem ranking (uso interno)."""
    st.caption(str(compared.get("ranking_note") or ""))
    mutant = compared.get("mutant_structure") or {}
    st.markdown(status_badge(str(mutant.get("status") or "UNAVAILABLE")), unsafe_allow_html=True)
    st.caption(str(mutant.get("reason") or ""))
    side = compared.get("side_by_side") or {}
    col_a, col_b = st.columns(2)
    for column, key, title in ((col_a, "a", "Variant A"), (col_b, "b", "Variant B")):
        card = side.get(key) or {}
        with column:
            st.markdown(f"**{title}**")
            st.markdown(
                meta_grid(
                    [
                        ("Assembly", card.get("assembly")),
                        ("Location", f"{card.get('contig')}:{card.get('position_1based')} {card.get('ref')}>{card.get('alt')}"),
                        ("Most severe consequence", card.get("most_severe_consequence")),
                        ("Transcripts", card.get("transcript_count")),
                        ("HGVS p.", card.get("hgvsp")),
                        ("Amino acids", card.get("amino_acids")),
                        ("Protein position", card.get("protein_position")),
                        ("ClinVar", card.get("clinvar_status")),
                        ("Domains", card.get("domains_status")),
                    ]
                ),
                unsafe_allow_html=True,
            )
    conflicts = compared.get("conflicts") or []
    if conflicts:
        st.markdown(status_badge("WARNING"), unsafe_allow_html=True)
        st.caption("Sources disagree. All values are shown; HelixScope does not take a majority vote.")
        for row in conflicts:
            st.caption(f"{row.get('field')}: {row.get('policy')}")
    export = compared.get("export") or {}
    st.download_button(
        "Download variant comparison JSON",
        data=json.dumps(export, indent=2, default=str),
        file_name=safe_download_filename("helixscope_variant_compare.json"),
        mime="application/json",
        key="cmp_var_dl",
    )


def _render_compare_proteins() -> None:
    """Identidade de sequencia + InterPro real (uso interno)."""
    col_a, col_b = st.columns(2)
    with col_a:
        seq_a = st.text_area("Protein A", height=120, key="cmp_prot_a")
        acc_a = st.text_input("UniProt A (optional InterPro)", key="cmp_prot_acc_a")
    with col_b:
        seq_b = st.text_area("Protein B", height=120, key="cmp_prot_b")
        acc_b = st.text_input("UniProt B (optional InterPro)", key="cmp_prot_acc_b")
    if not st.button("Compare proteins", key="cmp_prot_btn"):
        return
    domains_a = None
    domains_b = None
    try:
        if acc_a.strip():
            domains_a = protein_domains.fetch_protein_domains(accession=acc_a.strip())
        if acc_b.strip():
            domains_b = protein_domains.fetch_protein_domains(accession=acc_b.strip())
        result = comparative.compare_proteins(
            sequence_a=seq_a,
            sequence_b=seq_b,
            identifier_a=acc_a.strip() or "protein_a",
            identifier_b=acc_b.strip() or "protein_b",
            domains_a=domains_a,
            domains_b=domains_b,
            structure_comparison=st.session_state.get("compare_structures"),
        )
    except (comparative.ComparativeError, Exception) as exc:
        st.error(str(exc))
        return
    ident = result.get("identity") or {}
    st.markdown(
        meta_grid(
            [
                ("Method", ident.get("method")),
                ("Identity % (gapped)", ident.get("identity_pct")),
                ("Ungapped identity %", ident.get("ungapped_identity_pct")),
                ("Coverage A %", ident.get("coverage_seq1_pct")),
                ("Coverage B %", ident.get("coverage_seq2_pct")),
                ("Hash A", ident.get("hash_a")),
                ("Hash B", ident.get("hash_b")),
            ]
        ),
        unsafe_allow_html=True,
    )
    st.caption(str(ident.get("identity_definition") or ""))
    st.caption(str(result.get("disclaimer") or ""))
    st.caption(str((result.get("domains") or {}).get("note") or ""))
    st.session_state["compare_export"] = result.get("export")


def _render_compare_guides() -> None:
    """Comparacao de guias CRISPR com gate de ambito (uso interno)."""
    col_a, col_b = st.columns(2)
    with col_a:
        g_a = st.text_input("Guide A (20 nt proto-spacer)", key="cmp_guide_a")
        pam_a = st.text_input("PAM A", value="NGG", key="cmp_pam_a")
    with col_b:
        g_b = st.text_input("Guide B (20 nt proto-spacer)", key="cmp_guide_b")
        pam_b = st.text_input("PAM B", value="NGG", key="cmp_pam_b")
    search_a = st.session_state.get("crispr_offtarget_result_a") or st.session_state.get("crispr_search_result")
    search_b = st.session_state.get("crispr_offtarget_result_b")
    if not st.button("Compare guides", key="cmp_guide_btn"):
        return
    try:
        result = comparative.compare_guides(
            {"guide_sequence": g_a, "pam": pam_a},
            {"guide_sequence": g_b, "pam": pam_b},
            search_a=search_a if isinstance(search_a, Mapping) else None,
            search_b=search_b if isinstance(search_b, Mapping) else None,
        )
    except comparative.ComparativeError as exc:
        st.error(str(exc))
        return
    if not result.get("specificity_comparable"):
        st.markdown(status_badge("SCOPE_MISMATCH"), unsafe_allow_html=True)
    st.caption(str(result.get("specificity_reason") or ""))
    st.caption(str(result.get("genome_wide_note") or ""))
    side = result.get("side_by_side") or {}
    cols = st.columns(2)
    for column, key, title in ((cols[0], "a", "Guide A"), (cols[1], "b", "Guide B")):
        card = side.get(key) or {}
        with column:
            st.markdown(f"**{title}**")
            flags = []
            if card.get("test_only"):
                flags.append("TEST-ONLY")
            if card.get("genome_wide"):
                flags.append("GENOME_WIDE")
            if flags:
                st.markdown(badge_row(flags), unsafe_allow_html=True)
            st.markdown(
                meta_grid(
                    [
                        ("Sequence", card.get("guide_sequence")),
                        ("PAM", card.get("pam")),
                        ("System", card.get("system")),
                        ("Hits", card.get("n_hits")),
                        ("MIT Sguide", card.get("mit_sguide")),
                        ("MIT genome-wide", card.get("mit_genome_wide")),
                        ("CFD Sguide", card.get("cfd_sguide")),
                        ("CFD genome-wide", card.get("cfd_genome_wide")),
                        ("Scope", card.get("search_scope")),
                    ]
                ),
                unsafe_allow_html=True,
            )
    st.session_state["compare_export"] = result.get("export")


def _render_compare_evolution() -> None:
    """MSA -> residuo / padrao associado a grupo (uso interno)."""
    st.caption(
        "Paste a user-declared prealigned FASTA. Conservation is Shannon among "
        "analyzed sequences, not functional essentiality. Group-associated "
        "residue patterns are not adaptive mutations and not dN/dS."
    )
    fasta = st.text_area("Prealigned FASTA", height=180, key="cmp_evo_fasta")
    group = st.text_input("Group identifiers (comma-separated)", key="cmp_evo_group")
    column = st.number_input("MSA column (0-based)", min_value=0, value=0, step=1, key="cmp_evo_col")
    member = st.text_input("Member id for structure/MSA mapping", key="cmp_evo_member")
    if not st.button("Inspect evolutionary columns", key="cmp_evo_btn"):
        return
    try:
        msa_result = msa.import_prealigned_fasta(fasta)
        detail = msa.column_detail(msa_result, int(column))
        ids = [item.strip() for item in str(group).split(",") if item.strip()]
        associated = (
            comparative.group_associated_positions(msa_result=msa_result, group_ids=ids)
            if ids
            else None
        )
        mapped_col = None
        if member.strip():
            mapped_col = comparative.structure_residue_to_msa_column(
                msa_result=msa_result,
                member_id=member.strip(),
                ungapped_index=0,
            )
    except Exception as exc:
        st.error(str(exc))
        return
    st.caption(
        f"Column {detail.get('column')} conservation={detail.get('conservation')} "
        f"class={detail.get('variation_class')} gaps={detail.get('n_gaps')}. "
        "Highly conserved is not functionally essential."
    )
    if associated:
        st.caption(str(associated.get("statement") or ""))
        st.caption(f"{associated.get('n_positions')} group-associated positions.")
    if mapped_col:
        st.caption(f"First ungapped residue of {member} maps to MSA column {mapped_col.get('column')} ({mapped_col.get('status')}).")


def _render_compare_evidence() -> None:
    """Painel WHERE DID THIS COME FROM? (uso interno)."""
    export = st.session_state.get("compare_export")
    struct = st.session_state.get("compare_structures")
    variants = st.session_state.get("compare_variants")
    pack = export
    if not isinstance(pack, Mapping) and isinstance(struct, Mapping):
        pack = struct.get("export")
    if not isinstance(pack, Mapping) and isinstance(variants, Mapping):
        pack = variants.get("export")
    if not isinstance(pack, Mapping):
        st.caption("Run a comparison first. Evidence items are not invented.")
        return
    st.markdown(status_badge("RETRIEVED"), unsafe_allow_html=True)
    st.caption(str(pack.get("confidence_note") or ""))
    st.caption(f"Exported {pack.get('exported_at_utc')} version {pack.get('software_version')}.")
    conflicts = pack.get("conflicts") or []
    helix_workspace.render_explanation(
        "evidence",
        {
            "status": "RETRIEVED",
            "n_conflicts": len(conflicts) if isinstance(conflicts, list) else 0,
            "method": "evidence workspace",
            "source": "HelixScope",
        },
    )
    if conflicts:
        st.markdown(status_badge("PARTIAL"), unsafe_allow_html=True)
        st.markdown("**CONFLICTING EVIDENCE**")
        st.caption("Source conflicts (no majority vote). HelixScope does not decide truth.")
        for row in conflicts:
            st.caption(f"{row.get('field')}: {row.get('policy')}")
    for item in list(pack.get("evidence") or []):
        if isinstance(item, Mapping):
            st.markdown(helix_workspace.evidence_item_html(item), unsafe_allow_html=True)
    _render_evidence_table(list(pack.get("evidence") or []))
    st.download_button(
        "Download evidence pack",
        data=json.dumps(pack, indent=2, default=str),
        file_name=safe_download_filename("helixscope_evidence_pack.json"),
        mime="application/json",
        key="cmp_ev_dl",
    )


def render_phylogeny_workspace() -> None:
    """Phylogeny module: same panel as MSA, only when a completed MSA exists.

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    result = st.session_state.get("msa_result")
    collection = st.session_state.get("msa_collection")
    if not isinstance(collection, list):
        collection = []
    if isinstance(result, dict):
        current_hashes = [str(item.get("hash") or "") for item in collection]
        stored_hashes = [str(item) for item in list(result.get("input_hashes") or [])]
        if current_hashes != stored_hashes:
            _hide_stale_result(
                "Stored MSA belongs to a different sequence collection. "
                "Run MSA again to refresh."
            )
            result = None
    if not isinstance(result, dict):
        st.caption("No completed MSA is loaded.")
        st.caption(
            "Run MSA or load a valid pre-aligned FASTA before inferring a "
            "phylogenetic tree. Phylogeny engines are not missing; a tree is "
            "not invented from an empty alignment."
        )
        return
    _render_phylogeny_panel(result)


def _dispatch_active_module(module_id: str) -> None:
    """Render one scientific module. Does not recompute stored results.

    Args:
        module_id: Semantic id from ui.navigation.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    mapping = {
        "overview": helix_shell.render_overview,
        "settings": helix_shell.render_settings,
        "dna": render_dna_analysis,
        "rna": render_rna_analysis,
        "protein": render_protein_analysis,
        "alignment": render_alignment,
        "motif": render_motif_search,
        "crispr": render_crispr,
        "variant": render_variant_explorer,
        "ncbi": render_ncbi_fetch,
        "blast": render_blast_search,
        "references": helix_shell.render_references,
        "msa": render_msa_analysis,
        "phylogeny": render_phylogeny_workspace,
        "structure_3d": helix_shell.render_structure_hub,
        "compare": render_compare_workspace,
    }
    renderer = mapping.get(module_id)
    if renderer is None:
        st.error(f"Unknown module: {module_id}")
        return
    helix_shell.render_module_page(module_id, renderer)


_ACTIVE_MODULE = render_sidebar_navigation()
helix_shell.render_topbar(_ACTIVE_MODULE)
_dispatch_active_module(_ACTIVE_MODULE)
