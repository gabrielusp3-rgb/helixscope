"""Camada deterministica de interpretacao cientifica para nao especialistas.

Consome resultados ja calculados. Nao recomputa GC, Tm, RMSD, TM-score,
conservacao, CFD, MIT, identidade, consequencias de variante nem estruturas.
Nenhuma funcao importa Streamlit. Nao ha LLM nem texto generativo em runtime.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional

from . import provenance

EXPERIMENTAL_CLAIM_MARKERS: tuple[str, ...] = (
    "experimental coordinates",
    "experimentally determined coordinates",
    "measured experimentally",
)
"""Frases que uma explicacao PREDICTED/ILLUSTRATIVE nunca pode conter."""

DIAGNOSIS_MARKERS: tuple[str, ...] = (
    "this mutation causes disease",
    "this residue is essential",
    "this motif is functional",
    "this contact proves interaction",
    "this tree proves ancestry",
    "this guide is safe",
    "this sequence is healthy",
    "this variant is pathogenic",
)
"""Afirmacoes diagnosticas proibidas, independentemente do resultado."""


@dataclass(frozen=True)
class Explanation:
    """Interpretacao de duas camadas sobre um resultado ja existente.

    Attributes:
        kind: Identificador do tipo de resultado.
        title: Titulo visivel, com o valor real quando aplicavel.
        status: Estado de evidencia copiado do resultado (nunca inventado).
        plain_meaning: Camada 1, uma a tres frases.
        why_this_result: Porque este numero ou figura apareceu.
        interpretation: Como ler o valor sem extrapolar.
        limitations: O que o resultado nao prova.
        method: Metodo ou ferramenta ja declarados no resultado.
        source: Proveniencia ja declarada no resultado.
    """

    kind: str
    title: str
    status: str
    plain_meaning: str
    why_this_result: str
    interpretation: str
    limitations: str
    method: str
    source: str


def display_number(
    value: object,
    *,
    digits: int = 2,
    suffix: str = "",
) -> str:
    """Formata um numero cientifico. None/NaN nunca viram 0.

    Args:
        value: Quantidade ja calculada, ou ausente.
        digits: Casas decimais quando o valor e numerico finito.
        suffix: Unidade opcional (%, A, kcal/mol).

    Returns:
        Texto numerico ou "Unavailable".

    Raises:
        Nenhum.
    """
    if value is None or value == "":
        return "Unavailable"
    if isinstance(value, bool):
        return "Unavailable"
    try:
        number = float(value)
    except (TypeError, ValueError):
        text = str(value).strip()
        if not text or text.upper() in {"N/A", "NONE", "NULL", "NAN"}:
            return "Unavailable"
        return text
    if not math.isfinite(number):
        return "Unavailable"
    rendered = f"{number:.{int(digits)}f}"
    return f"{rendered}{suffix}" if suffix else rendered


def display_text(value: object, *, empty: str = "Unavailable") -> str:
    """Texto de campo ausente. String vazia nao vira sucesso.

    Args:
        value: Campo do resultado.
        empty: Rotulo para ausencia.

    Returns:
        Texto ou empty.

    Raises:
        Nenhum.
    """
    if value is None:
        return empty
    text = str(value).strip()
    if not text or text.upper() in {"N/A", "NONE", "NULL", "NAN"}:
        return empty
    return text


def _status_of(values: Mapping[str, Any], default: str = "COMPUTED") -> str:
    raw = values.get("status") or values.get("evidence_status") or default
    token = str(raw or default).strip().upper().replace(" ", "_")
    if token in provenance.EVIDENCE_STATUSES or token in {
        "LIVE_VALIDATED",
        "REMOTE_VALIDATED",
        "TEST_ONLY",
        "HEURISTIC",
        "NOT_INSTALLED",
        "DETECTED",
        "COMPLETED",
        "NO_HITS",
        "READY",
        "AVAILABLE",
    }:
        return token
    return str(raw or default).strip() or default


def _limitations_for_status(status: str, extra: str) -> str:
    token = str(status or "").upper()
    prefix = ""
    if token == "PREDICTED":
        prefix = (
            "This is a computational prediction, not an experimentally determined "
            "structure or measurement. "
        )
    elif token == "ILLUSTRATIVE":
        prefix = (
            "This geometry is illustrative. It is not experimental coordinates "
            "and must not be treated as a deposited structure. "
        )
    elif token == "HEURISTIC":
        prefix = (
            "This value is a local heuristic, not a published trained model. "
        )
    elif token == "PARTIAL":
        prefix = "This result is partial: some requested fields were not obtained. "
    elif token in {"UNAVAILABLE", "NOT_INSTALLED", "ERROR"}:
        prefix = (
            "This computation is unavailable. Absence is not a numeric zero "
            "and not a successful result. "
        )
    elif token == "RESOURCE_LIMIT":
        prefix = (
            "Hardware or resource admission blocked this computation. "
            "That is not a scientific negative finding. "
        )
    elif token == "EXPERIMENTAL":
        prefix = (
            "Coordinates or records retrieved from an experiment are still a "
            "model of that experiment, not every conformation the molecule can adopt. "
        )
    body = str(extra or "").strip()
    return (prefix + body).strip()


def _method_source(values: Mapping[str, Any], *, method: str = "", source: str = "") -> tuple[str, str]:
    used_method = display_text(values.get("method") or values.get("algorithm") or method, empty="")
    used_source = display_text(
        values.get("source") or values.get("engine") or values.get("tool") or source,
        empty="",
    )
    return used_method or method, used_source or source


def _unavailable_explanation(kind: str, values: Mapping[str, Any], title: str) -> Explanation:
    status = _status_of(values, "UNAVAILABLE")
    reason = display_text(values.get("reason") or values.get("limitations") or "")
    method, source = _method_source(values)
    return Explanation(
        kind=kind,
        title=title,
        status=status,
        plain_meaning=(
            "This result is not available in the current session. "
            "HelixScope did not invent a substitute value."
        ),
        why_this_result=reason or "The underlying method did not return a usable result.",
        interpretation="Treat missing fields as unknown, not as zero or as a negative finding.",
        limitations=_limitations_for_status(
            status,
            "Do not read an empty panel as evidence that the molecule lacks the property.",
        ),
        method=method or "not run",
        source=source or "HelixScope",
    )


def _dna(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    length = display_number(values.get("length"), digits=0)
    gc = display_number(values.get("gc_percent") if values.get("gc_percent") is not None else values.get("gc_content"), digits=2, suffix="%")
    at = display_number(values.get("at_percent") if values.get("at_percent") is not None else values.get("at_content"), digits=2, suffix="%")
    if status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("dna", values, "DNA ANALYSIS")
    method, source = _method_source(values, method="composition from the analyzed sequence", source="HelixScope")
    return Explanation(
        kind="dna",
        title=f"DNA ANALYSIS — {length} bp, GC {gc}",
        status=status,
        plain_meaning=(
            f"The analyzed DNA string contains {length} base pairs. "
            f"{gc} of the canonical nucleotides are guanine (G) or cytosine (C); "
            f"AT content is {at}."
        ),
        why_this_result=(
            "These percentages are counted from the sequence you submitted after "
            "normalization. Ambiguous bases are excluded from the GC/AT denominator."
        ),
        interpretation=(
            "GC content describes composition and can influence duplex stability. "
            "Length is the number of residues in the analyzed string, not genome size "
            "unless you pasted a complete chromosome."
        ),
        limitations=_limitations_for_status(
            status,
            "GC content alone does not identify the organism, determine gene function, "
            "or indicate that a sequence is biologically good or bad.",
        ),
        method=method,
        source=source,
    )


def _dna_tm(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    tm = display_number(values.get("value_c"), digits=2, suffix=" C")
    method, source = _method_source(values, method="melting-temperature regime", source="HelixScope")
    if tm == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return Explanation(
            kind="dna_tm",
            title="MELTING TEMPERATURE — Unavailable",
            status=status if status != "COMPUTED" else "UNAVAILABLE",
            plain_meaning=(
                "A melting temperature was not computed for this sequence and method. "
                "That is not a Tm of 0 C."
            ),
            why_this_result=display_text(
                values.get("reason") or values.get("method") or "The selected Tm regime does not apply."
            ),
            interpretation="Use a regime that matches oligo length and salt, or leave Tm unavailable.",
            limitations=_limitations_for_status(
                "UNAVAILABLE",
                "HelixScope does not substitute Wallace or SantaLucia values when the regime is out of range.",
            ),
            method=method,
            source=source,
        )
    return Explanation(
        kind="dna_tm",
        title=f"MELTING TEMPERATURE — {tm}",
        status=status,
        plain_meaning=(
            f"Under the declared method, the estimated duplex melting temperature is {tm}."
        ),
        why_this_result=(
            f"Method: {method}. The value depends on oligo length, composition and, "
            "when applicable, salt and strand concentration."
        ),
        interpretation=(
            "Tm is an estimate for a nucleic-acid duplex under the stated conditions, "
            "not a PCR cycling protocol by itself."
        ),
        limitations=_limitations_for_status(
            status,
            "This Tm is not a guarantee of primer performance and is not nearest-neighbor "
            "unless the method name says SantaLucia.",
        ),
        method=method,
        source=source,
    )


def _dna_orf(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    n_orfs = display_number(values.get("n_orfs"), digits=0)
    n_codons = display_number(values.get("n_codons"), digits=0)
    method, source = _method_source(values, method="ORF/CDS selection", source="HelixScope")
    return Explanation(
        kind="dna_orf",
        title=f"OPEN READING FRAMES — {n_orfs}",
        status=status,
        plain_meaning=(
            f"HelixScope reports {n_orfs} selected ORF/CDS regions and {n_codons} codons "
            "from those regions, not from slicing the whole sequence into triplets."
        ),
        why_this_result=display_text(
            values.get("method_note")
            or "ORFs require a start codon, a stop codon, a minimum length, and overlap filtering."
        ),
        interpretation=(
            "Codon totals must correspond to real CDS/ORF spans. Nested ORFs that share "
            "nucleotides are not all counted as independent coding sequences."
        ),
        limitations=_limitations_for_status(
            status,
            "Predicted ORFs are not experimental gene annotations. Absence of an ORF "
            "is not proof that the region is non-coding in vivo.",
        ),
        method=method,
        source=source,
    )


def _rna(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    length = display_number(values.get("length"), digits=0)
    gc_raw = values.get("gc_percent")
    if gc_raw is None:
        gc_raw = values.get("gc_content")
    gc = display_number(gc_raw, digits=2, suffix="%")
    method, source = _method_source(values, method="RNA composition", source="HelixScope")
    return Explanation(
        kind="rna",
        title=f"RNA ANALYSIS — {length} nt, GC {gc}",
        status=status,
        plain_meaning=(
            f"The analyzed RNA string contains {length} nucleotides. "
            f"GC content among canonical bases is {gc}."
        ),
        why_this_result="Counts come from the RNA alphabet after validation (A, C, G, U).",
        interpretation=(
            "The RNA tab analyzes the pasted RNA. It does not silently transcribe DNA "
            "from another tab."
        ),
        limitations=_limitations_for_status(
            status,
            "Composition does not assign function, and this panel is not an RNA 3D fold.",
        ),
        method=method,
        source=source,
    )


def _rna_fold(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "PREDICTED")
    mfe = display_number(values.get("mfe_kcal_mol"), digits=2, suffix=" kcal/mol")
    backend = display_text(values.get("backend") or values.get("tool"))
    method, source = _method_source(
        values,
        method="ViennaRNA MFE",
        source="ViennaRNA",
    )
    if status in {"UNAVAILABLE", "TOOL_NOT_INSTALLED", "ERROR"} or mfe == "Unavailable":
        return _unavailable_explanation("rna_fold", values, "RNA SECONDARY STRUCTURE")
    return Explanation(
        kind="rna_fold",
        title=f"PREDICTED RNA SECONDARY STRUCTURE — MFE {mfe}",
        status="PREDICTED",
        plain_meaning=(
            f"ViennaRNA predicted a minimum-free-energy secondary structure with "
            f"reported MFE {mfe}. Dot-bracket is a 2D pairing diagram, not a 3D fold."
        ),
        why_this_result=(
            f"Backend {backend}. The algorithm is thermodynamic MFE (Zuker-style DP) "
            "with the library energy parameters, not an experimental structure."
        ),
        interpretation=(
            "A more negative MFE means a more stable predicted fold under the model, "
            "not that this conformation was observed in a cell. Unpaired positions "
            "are unpaired in this MFE, not proven single-stranded in vivo."
        ),
        limitations=_limitations_for_status(
            "PREDICTED",
            "This is not an experimental RNA structure and is not converted into RNA 3D "
            "coordinates. Pseudoknots are not represented in this MFE alphabet.",
        ),
        method=method,
        source=source,
    )


def _protein(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    length = display_number(values.get("length"), digits=0)
    mw = display_number(values.get("molecular_weight_kda"), digits=3, suffix=" kDa")
    pi = display_number(values.get("isoelectric_point"), digits=2)
    gravy = display_number(values.get("gravy"), digits=3)
    method, source = _method_source(values, method="protein physicochemical indices", source="HelixScope")
    return Explanation(
        kind="protein",
        title=f"PROTEIN ANALYSIS — {length} aa, MW {mw}, pI {pi}",
        status=status,
        plain_meaning=(
            f"The analyzed protein has {length} amino acids, approximate molecular "
            f"weight {mw}, theoretical isoelectric point {pi}, and GRAVY {gravy}."
        ),
        why_this_result=(
            "These indices are computed from the amino-acid string using the declared "
            "scales (pI, GRAVY, instability, aliphatic index, charge at the stated pH)."
        ),
        interpretation=(
            "pI is theoretical for the isolated polypeptide. GRAVY is a hydropathy "
            "average, not a localization prediction."
        ),
        limitations=_limitations_for_status(
            status,
            "Sequence indices do not establish fold, activity, or pathogenicity.",
        ),
        method=method,
        source=source,
    )


def _alignment(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    identity = display_number(values.get("identity_pct"), digits=2, suffix="%")
    score = display_number(values.get("score"), digits=2)
    length = display_number(values.get("length"), digits=0)
    algorithm = display_text(values.get("algorithm") or values.get("method"))
    method, source = _method_source(values, method=algorithm, source="HelixScope pairwise alignment")
    return Explanation(
        kind="alignment",
        title=f"PAIRWISE ALIGNMENT — identity {identity}, score {score}",
        status=status,
        plain_meaning=(
            f"{algorithm} aligned two sequences to length {length} columns with "
            f"identity {identity} (gap columns included in the stated denominator) "
            f"and score {score}."
        ),
        why_this_result="The score uses the selected substitution matrix and gap penalties.",
        interpretation=(
            "Identity is a percentage of aligned columns under the declared formula. "
            "Global (Needleman-Wunsch) and local (Smith-Waterman) alignments answer "
            "different questions and are not interchangeable."
        ),
        limitations=_limitations_for_status(
            status,
            "Alignment identity is not phylogenetic proof, not homology by itself, "
            "and not a functional annotation.",
        ),
        method=method,
        source=source,
    )


def _motif(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    n_hits = display_number(values.get("n_hits"), digits=0)
    pattern = display_text(values.get("pattern"))
    method, source = _method_source(values, method="IUPAC motif scan", source="HelixScope")
    return Explanation(
        kind="motif",
        title=f"MOTIF SEARCH — {n_hits} hits",
        status=status,
        plain_meaning=(
            f"The pattern {pattern} occurs {n_hits} time(s) in the analyzed sequence "
            "on the reported strand(s)."
        ),
        why_this_result="Matches are string/IUPAC occurrences with 0-based half-open coordinates internally.",
        interpretation="A motif hit is a sequence occurrence, not a demonstrated binding event.",
        limitations=_limitations_for_status(
            status,
            "Motif occurrence alone does not prove biological function.",
        ),
        method=method,
        source=source,
    )


def _ncbi(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "RETRIEVED")
    accession = display_text(values.get("accession") or values.get("identifier"))
    organism = display_text(values.get("organism"))
    method, source = _method_source(values, method="NCBI Entrez", source="NCBI")
    return Explanation(
        kind="ncbi",
        title=f"NCBI RECORD — {accession}",
        status=status,
        plain_meaning=(
            f"This record was retrieved from NCBI. Accession/version shown: {accession}. "
            f"Organism field: {organism}."
        ),
        why_this_result="Identifiers and feature tables come from the Entrez record, not from local invention.",
        interpretation="The organism string is what NCBI stored for this accession, not a taxonomic reclassification by HelixScope.",
        limitations=_limitations_for_status(
            status,
            "A retrieved record can be outdated or misannotated at the source. HelixScope does not correct NCBI.",
        ),
        method=method,
        source=source,
    )


def _blast(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    backend = display_text(values.get("backend") or values.get("source"), empty="BLAST")
    n_hits = display_number(values.get("n_hits"), digits=0)
    program = display_text(values.get("program"))
    local = bool(values.get("not_remote_ncbi"))
    method, source = _method_source(values, method=program, source=backend)
    if status in {"NO_HITS"}:
        meaning = (
            f"{backend} ({program}) reported no hits above the search settings. "
            "NO_HITS is not a parse failure and not a fabricated alignment."
        )
    else:
        meaning = (
            f"{backend} ({program}) returned {n_hits} hit(s). "
            + (
                "This is local BLAST+ against a declared database, not NCBI nt/nr."
                if local
                else "This is remote NCBI BLAST, not a local BLAST+ search."
            )
        )
    return Explanation(
        kind="blast",
        title=f"BLAST — {program}, {n_hits} hits",
        status=status,
        plain_meaning=meaning,
        why_this_result=(
            "Bit scores, E-values and HSP coordinates come from the BLAST engine. "
            "Local and remote backends remain distinct."
        ),
        interpretation=(
            "E-value is a database-size-dependent statistic, not a probability of homology "
            "in the biological sense by itself. Identity is per HSP as reported."
        ),
        limitations=_limitations_for_status(
            status,
            "A hit is not a gene name assignment. A tiny local fixture database is not nt or nr.",
        ),
        method=method,
        source=source,
    )


def _msa(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    n_seq = display_number(values.get("n_sequences"), digits=0)
    length = display_number(values.get("alignment_length"), digits=0)
    tool = display_text(values.get("tool"))
    method, source = _method_source(values, method=display_text(values.get("method")), source=tool)
    return Explanation(
        kind="msa",
        title=f"MULTIPLE SEQUENCE ALIGNMENT — {n_seq} sequences, {length} columns",
        status=status,
        plain_meaning=(
            f"This MSA contains {n_seq} aligned sequences and {length} columns. "
            f"Tool: {tool}."
        ),
        why_this_result="Every row is padded to the same column count. Consensus and conservation are derived from these columns.",
        interpretation=(
            "Conservation in this alignment applies to the sequences included here. "
            "A conserved column is not automatically a catalytic residue."
        ),
        limitations=_limitations_for_status(
            status,
            "MSA conservation does not prove function or pathogenicity. Prealigned input is not a Clustal run unless the tool name says so.",
        ),
        method=method,
        source=source,
    )


def _phylogeny(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    tool = display_text(values.get("tool") or values.get("method"))
    n_leaves = display_number(values.get("n_leaves"), digits=0)
    method, source = _method_source(values, method=display_text(values.get("algorithm") or tool), source=tool)
    extra_limit = ""
    if "fasttree" in tool.lower() or "fasttree" in method.lower():
        extra_limit = (
            "FastTree is approximate maximum-likelihood. It is not the same inference as IQ-TREE."
        )
    elif "iq-tree" in tool.lower() or "iqtree" in method.lower():
        extra_limit = (
            "UFBoot and SH-aLRT are support statistics on this alignment, not taxonomic proof."
        )
    elif "neighbor" in method.lower() or method.upper() == "NJ":
        extra_limit = "Neighbor-Joining is a distance method, not a substitution-model ML search."
    elif "upgma" in method.lower():
        extra_limit = "UPGMA assumes a molecular clock and can distort branch lengths when that is false."
    return Explanation(
        kind="phylogeny",
        title=f"PHYLOGENETIC TREE — {tool}, {n_leaves} leaves",
        status=status,
        plain_meaning=(
            f"A tree with {n_leaves} leaves was inferred by {tool}. Branch lengths, "
            "when present, are in the units of that method."
        ),
        why_this_result="The tree is computed from the current MSA (or distances derived from it), not from a taxonomy database.",
        interpretation=(
            "Topology is a hypothesis for this alignment. Bootstrap or UFBoot values "
            "summarize resampling support, not certainty of ancestry."
        ),
        limitations=_limitations_for_status(
            status,
            extra_limit
            + " A tree does not prove species ancestry or that sequences are orthologs.",
        ),
        method=method,
        source=source,
    )


def _crispr(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "PREDICTED")
    n_guides = display_number(values.get("n_guides"), digits=0)
    pam = display_text(values.get("pam") or "NGG")
    scope = display_text(values.get("reference_scope") or values.get("reference_id"))
    method, source = _method_source(values, method="SpCas9 guide scan", source="HelixScope CRISPR")
    test_ref = bool(values.get("test_reference") or "TEST" in scope.upper())
    extra = (
        "Off-target search on a TEST REFERENCE is not a genome-wide result."
        if test_ref
        else "Off-target scope is only as wide as the declared reference."
    )
    return Explanation(
        kind="crispr",
        title=f"CRISPR GUIDES — {n_guides} candidates, PAM {pam}",
        status="PREDICTED",
        plain_meaning=(
            f"HelixScope listed {n_guides} candidate guide(s) for PAM {pam}. "
            f"Reference scope: {scope}."
        ),
        why_this_result=(
            "Guides are computed from the pasted DNA. On-target and specificity numbers "
            "are only those HelixScope can calculate offline; trained genome-wide scores "
            "are declared unavailable when models are absent."
        ),
        interpretation=(
            "A high local heuristic score is not a safety certificate. MIT/CFD labels "
            "must match the actual method; a local proxy is not CFD or MIT on GRCh38."
        ),
        limitations=_limitations_for_status("PREDICTED", extra + " This guide is not experimentally validated."),
        method=method,
        source=source,
    )


def _variant(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    hgvs = display_text(values.get("hgvs") or values.get("variant_id"))
    assembly = display_text(values.get("assembly"))
    consequence = display_text(values.get("most_severe_consequence"))
    clinvar = display_text(values.get("clinvar_significance") or values.get("clinvar"))
    method, source = _method_source(values, method="VEP / ClinVar retrieval", source="Ensembl VEP / NCBI ClinVar")
    return Explanation(
        kind="variant",
        title=f"VARIANT — {hgvs} ({assembly})",
        status=status,
        plain_meaning=(
            f"The genomic variant {hgvs} is interpreted on assembly {assembly}. "
            f"VEP annotates consequence as {consequence}."
        ),
        why_this_result=(
            "Transcript choice, protein change and database records are retrieved or "
            "mapped; they are not invented coordinates."
        ),
        interpretation=(
            f"ClinVar reports {clinvar}. That is a submitted classification in ClinVar, "
            "not a HelixScope diagnosis."
        ),
        limitations=_limitations_for_status(
            status,
            "ClinVar does not make this variant pathogenic by itself. Conservation or "
            "a protein change does not prove disease. Mapping ambiguity remains UNMAPPED "
            "when transcripts disagree.",
        ),
        method=method,
        source=source,
    )


def _structure(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    kind_label = display_text(values.get("kind") or values.get("kind_label") or status).lower()
    structure_id = display_text(values.get("structure_id") or values.get("entry_id"))
    source_name = display_text(values.get("source") or values.get("structure_source"))
    method, source = _method_source(values, method=display_text(values.get("method")), source=source_name)
    if "illustrative" in kind_label or status == "ILLUSTRATIVE":
        return Explanation(
            kind="structure",
            title=f"ILLUSTRATIVE STRUCTURE — {structure_id}",
            status="ILLUSTRATIVE",
            plain_meaning=(
                "The displayed nucleic-acid geometry is an illustrative model generated "
                "from sequence, not coordinates measured in an experiment."
            ),
            why_this_result="HelixScope builds a regular helix for visualization when no deposited structure is loaded.",
            interpretation="Use it to inspect sequence register, not as a crystallographic model.",
            limitations=_limitations_for_status(
                "ILLUSTRATIVE",
                "This is not a PDB deposition and not a molecular-dynamics trajectory.",
            ),
            method=method or "illustrative nucleic-acid geometry",
            source=source,
        )
    if "predict" in kind_label or status == "PREDICTED" or "alphafold" in source_name.lower():
        return Explanation(
            kind="structure",
            title=f"PREDICTED STRUCTURE — {structure_id}",
            status="PREDICTED",
            plain_meaning=(
                f"Structure {structure_id} was produced by a computational prediction "
                f"model ({source_name}), not by a direct experimental structure determination."
            ),
            why_this_result="AlphaFold or equivalent predicted coordinates were retrieved from the stated database.",
            interpretation=(
                "Prediction confidence, when shown (for example pLDDT), applies to this "
                "model. Low-confidence regions should not be treated as experimental fact."
            ),
            limitations=_limitations_for_status(
                "PREDICTED",
                "A predicted model is not an experimental structure and does not show every conformation.",
            ),
            method=method or "computational structure prediction",
            source=source,
        )
    return Explanation(
        kind="structure",
        title=f"EXPERIMENTAL STRUCTURE — {structure_id}",
        status="EXPERIMENTAL" if status in {"EXPERIMENTAL", "RETRIEVED", "COMPUTED"} else status,
        plain_meaning=(
            f"The displayed coordinates for {structure_id} were retrieved from an "
            f"experimentally determined structure ({source_name})."
        ),
        why_this_result="Atoms come from the PDB/mmCIF file, not from a generated helix.",
        interpretation=(
            "This is a structural model derived from an experiment (crystal, NMR, or EM "
            "as deposited). It is not a molecular-dynamics simulation."
        ),
        limitations=_limitations_for_status(
            "EXPERIMENTAL",
            "Deposited coordinates do not show every conformation the molecule can adopt.",
        ),
        method=method or "PDB/mmCIF coordinates",
        source=source,
    )


def _compare(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "RETRIEVED")
    engine = display_text(values.get("engine") or values.get("source"))
    rmsd_g = display_number(values.get("rmsd_global_angstrom"), digits=2, suffix=" A")
    rmsd_b = display_number(values.get("rmsd_block0_angstrom"), digits=2, suffix=" A")
    tm = "Unavailable"
    tm_scores = values.get("tm_scores")
    if isinstance(tm_scores, (list, tuple)) and tm_scores:
        first = tm_scores[0] if isinstance(tm_scores[0], Mapping) else {}
        tm = display_number(first.get("value") if isinstance(first, Mapping) else tm_scores[0], digits=2)
    elif values.get("tm_score") is not None:
        tm = display_number(values.get("tm_score"), digits=2)
    n_pairs = display_number(values.get("n_aligned_residue_pairs"), digits=0)
    method, source = _method_source(values, method=display_text(values.get("method_display") or values.get("method")), source=engine)
    location = display_text(values.get("engine_location"), empty="unspecified")
    return Explanation(
        kind="compare",
        title=f"STRUCTURAL COMPARISON — global RMSD {rmsd_g}",
        status=status,
        plain_meaning=(
            f"Engine {engine} ({location}) aligned {n_pairs} residue pairs. "
            f"Global RMSD is {rmsd_g}; block-0 RMSD is {rmsd_b}; TM-score is {tm}."
        ),
        why_this_result=(
            "RMSD and TM-score are the engine-reported values for the atoms it fitted. "
            "Block RMSD and global RMSD are kept separate and are not averaged."
        ),
        interpretation=(
            "A lower RMSD means closer geometric agreement for the compared atoms after "
            "superposition. RMSD depends strongly on coverage and which atoms were aligned. "
            "TM-score is a length-normalized fold similarity, not a percent identity."
        ),
        limitations=_limitations_for_status(
            status,
            "RMSD or TM-score does not by itself establish equivalent biological function. "
            "Unmapped residues are unmapped; HelixScope does not assign the nearest residue. "
            "US-align and the RCSB Alignment API are different backends.",
        ),
        method=method,
        source=source,
    )


def _compare_rmsd_block(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "RETRIEVED")
    rmsd_b = display_number(values.get("rmsd_block0_angstrom"), digits=2, suffix=" A")
    method, source = _method_source(values, method="block RMSD", source=display_text(values.get("engine")))
    return Explanation(
        kind="compare_rmsd_block",
        title=f"BLOCK RMSD — {rmsd_b}",
        status=status,
        plain_meaning=(
            f"Block RMSD is {rmsd_b}. It applies to the aligned atoms of block 0, "
            "not to the global summary RMSD."
        ),
        why_this_result="Flexible or multi-block methods may report per-block RMSD distinct from the global value.",
        interpretation="Do not substitute block RMSD for global RMSD or the reverse.",
        limitations=_limitations_for_status(
            status,
            "A single block overlay is PARTIAL when additional blocks exist. "
            "This is not a TM-score.",
        ),
        method=method,
        source=source,
    )


def _evidence(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    n_conflicts = display_number(values.get("n_conflicts"), digits=0)
    method, source = _method_source(values, method="evidence workspace", source="HelixScope")
    return Explanation(
        kind="evidence",
        title=f"EVIDENCE — {n_conflicts} documented conflict(s)",
        status=status,
        plain_meaning=(
            "Each row is a retrieved or computed statement with its source. "
            f"HelixScope recorded {n_conflicts} conflict(s) between sources when they disagree."
        ),
        why_this_result="Conflicts are listed; they are not resolved by majority vote.",
        interpretation="When sources disagree, both remain visible. There is no confidence_score invented to pick a winner.",
        limitations=_limitations_for_status(
            status,
            "Evidence assembly is not a diagnosis and not a probability of truth.",
        ),
        method=method,
        source=source,
    )


def _evolution(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    conservation = display_text(values.get("conservation_label") or values.get("conservation"))
    method, source = _method_source(values, method="MSA column conservation", source="HelixScope")
    return Explanation(
        kind="evolution",
        title=f"CONSERVATION — {conservation}",
        status=status,
        plain_meaning=(
            f"Within this MSA, the selected column is described as {conservation}."
        ),
        why_this_result="Conservation is computed from the aligned sequences currently loaded, optionally viewed next to a tree or structure when those objects exist.",
        interpretation=(
            "High conservation can suggest evolutionary constraint in this dataset. "
            "It does not prove the residue is essential or that mutation is pathogenic."
        ),
        limitations=_limitations_for_status(
            status,
            "Conservation is not function. High conservation does not prove that the "
            "residue is essential or that mutation is pathogenic. A tree linked to "
            "this MSA is not independent proof of ancestry.",
        ),
        method=method,
        source=source,
    )


def _sasa(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    sasa = display_number(values.get("sasa_angstrom2"), digits=1, suffix=" A^2")
    method, source = _method_source(values, method="Shrake-Rupley", source="HelixScope")
    if sasa == "Unavailable":
        return _unavailable_explanation("sasa", values, "SOLVENT-ACCESSIBLE SURFACE AREA")
    return Explanation(
        kind="sasa",
        title=f"SASA — {sasa}",
        status=status,
        plain_meaning=(
            f"Shrake-Rupley solvent-accessible surface area for the selected atoms is {sasa}."
        ),
        why_this_result="The calculation uses deposited or predicted coordinates already loaded, not a sequence hydropathy scale.",
        interpretation="SASA is a geometric accessibility estimate, not a binding affinity.",
        limitations=_limitations_for_status(
            status,
            "SASA does not prove a residue is functional or that two chains interact.",
        ),
        method=method,
        source=source,
    )


def _contacts(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    n_contacts = display_number(values.get("n_contacts"), digits=0)
    method, source = _method_source(values, method="distance contacts", source="HelixScope")
    return Explanation(
        kind="contacts",
        title=f"ATOMIC CONTACTS — {n_contacts}",
        status=status,
        plain_meaning=(
            f"HelixScope counted {n_contacts} interatomic contact(s) under the declared distance cutoff."
        ),
        why_this_result="Contacts are geometric distances in the loaded coordinates.",
        interpretation="A contact is proximity in this model, not proof of a biological interaction.",
        limitations=_limitations_for_status(
            status,
            "This contact does not prove interaction in the cell.",
        ),
        method=method,
        source=source,
    )


def _entropy(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    bits = display_number(values.get("shannon_entropy") or values.get("entropy"), digits=3, suffix=" bits")
    if bits == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("entropy", values, "SHANNON ENTROPY — Unavailable")
    method, source = _method_source(values, method="Shannon entropy over the analyzed alphabet", source="HelixScope")
    return Explanation(
        kind="entropy",
        title=f"SHANNON ENTROPY — {bits}",
        status=status,
        plain_meaning=(
            f"The residue composition of this string has Shannon entropy {bits}. "
            "Higher entropy means the letters are more evenly used."
        ),
        why_this_result="Entropy is counted from the normalized sequence, not from a genome annotation.",
        interpretation=(
            "Low entropy can mean a biased alphabet (for example many A/T). "
            "It does not mean the sequence is unstructured in a cell."
        ),
        limitations=_limitations_for_status(
            status,
            "Entropy is a composition statistic. It is not complexity of an organism and not evidence of function.",
        ),
        method=method,
        source=source,
    )


def _gravy(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    gravy = display_number(values.get("gravy"), digits=3)
    if gravy == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("gravy", values, "GRAVY — Unavailable")
    method, source = _method_source(values, method="Kyte-Doolittle GRAVY", source="HelixScope")
    return Explanation(
        kind="gravy",
        title=f"GRAVY — {gravy}",
        status=status,
        plain_meaning=(
            f"The Grand Average of Hydropathy (GRAVY) of this protein string is {gravy}. "
            "Positive values are more hydrophobic on the Kyte-Doolittle scale."
        ),
        why_this_result="Each residue contributes its Kyte-Doolittle hydropathy; the mean is GRAVY.",
        interpretation="GRAVY is a sequence average. It is not a measured solubility or membrane-insertion assay.",
        limitations=_limitations_for_status(
            status,
            "GRAVY does not identify a protein, prove localization, or replace experimental hydrophobicity data.",
        ),
        method=method,
        source=source,
    )


def _pi(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    pi = display_number(values.get("isoelectric_point") or values.get("pI"), digits=2)
    if pi == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("pi", values, "ISOELECTRIC POINT — Unavailable")
    method, source = _method_source(
        values,
        method="theoretical pI (ProtParam / Bjellqvist pKa)",
        source="HelixScope",
    )
    return Explanation(
        kind="pi",
        title=f"ISOELECTRIC POINT — {pi}",
        status=status,
        plain_meaning=(
            f"The theoretical isoelectric point of this standard-residue string is {pi}. "
            "That is the pH where the computed net charge is zero under the declared pKa set."
        ),
        why_this_result="pI is solved from residue pKa values. It is not an IEF measurement.",
        interpretation="A theoretical pI can differ from an experimental pI because of modifications, folding, and the pKa model.",
        limitations=_limitations_for_status(
            status,
            "This pI is theoretical. It is not a clinical marker and not an experimental isoelectric focusing result.",
        ),
        method=method,
        source=source,
    )


def _instability(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    index = display_number(values.get("instability_index"), digits=2)
    if index == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("instability", values, "INSTABILITY INDEX — Unavailable")
    method, source = _method_source(values, method="Guruprasad instability index via ProtParam", source="HelixScope")
    return Explanation(
        kind="instability",
        title=f"INSTABILITY INDEX — {index}",
        status=status,
        plain_meaning=(
            f"The Guruprasad instability index of this protein string is {index}. "
            "Values above 40 are historically called unstable in that dipeptide statistic."
        ),
        why_this_result="The index is a weighted dipeptide statistic from the sequence, not a measured half-life.",
        interpretation="Do not read this number as in-cell protein stability or degradation rate.",
        limitations=_limitations_for_status(
            status,
            "The instability index is a heuristic sequence statistic. It is not an experimental stability assay.",
        ),
        method=method,
        source=source,
    )


def _rscu(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    method, source = _method_source(values, method="relative synonymous codon usage", source="HelixScope")
    n = display_number(values.get("n_codons") or values.get("codon_count"), digits=0)
    return Explanation(
        kind="rscu",
        title="RSCU",
        status=status,
        plain_meaning=(
            "Relative synonymous codon usage (RSCU) compares how often each codon is used "
            "versus an equal-usage expectation within its amino-acid family."
        ),
        why_this_result=(
            f"RSCU was computed from the coding-codon table derived for this RNA string "
            f"({n} codons counted when that count is available)."
        ),
        interpretation=(
            "RSCU near 1.0 means that synonym is used about as often as a uniform model. "
            "It is not proof of translational selection in an organism."
        ),
        limitations=_limitations_for_status(
            status,
            "RSCU requires a coding-frame codon table. Arbitrary non-coding RNA is not a CDS. "
            "RSCU does not identify the species.",
        ),
        method=method,
        source=source,
    )


def _enc(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    enc = display_number(values.get("enc") or values.get("value"), digits=2)
    method, source = _method_source(values, method="effective number of codons (ENC)", source="HelixScope")
    if enc == "Unavailable" and status not in {"UNAVAILABLE", "ERROR"}:
        status = str(values.get("status") or status)
    return Explanation(
        kind="enc",
        title=f"ENC — {enc}",
        status=status,
        plain_meaning=(
            f"The effective number of codons (ENC) for this codon table is {enc}. "
            "ENC near 61 means little synonymous bias; lower values mean stronger codon preference in this string."
        ),
        why_this_result="ENC is computed from synonymous-codon frequencies in the analyzed codon set.",
        interpretation="ENC is a sequence statistic of codon evenness. It is not an organism-level codon-adaptation proof.",
        limitations=_limitations_for_status(
            status,
            "ENC is only meaningful for a coding codon set. It does not identify the gene or the species.",
        ),
        method=method,
        source=source,
    )


def _cai(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    cai = display_number(values.get("cai") or values.get("value"), digits=3)
    method, source = _method_source(values, method="codon adaptation index against the declared reference table", source="HelixScope")
    return Explanation(
        kind="cai",
        title=f"CAI — {cai}",
        status=status,
        plain_meaning=(
            f"The codon adaptation index (CAI) for this codon table is {cai}. "
            "It scores how close synonymous codon use is to a reference highly expressed set."
        ),
        why_this_result="CAI uses the reference table declared in the result, not an inferred organism.",
        interpretation="CAI is relative to that reference table. It is not proof of high expression in a cell.",
        limitations=_limitations_for_status(
            status,
            "CAI is not a measured expression level. The reference table must match the intended comparison.",
        ),
        method=method,
        source=source,
    )


def _identity_metric(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    identity = display_number(values.get("identity_pct"), digits=2, suffix="%")
    method, source = _method_source(values, method="pairwise identity over the alignment columns", source="HelixScope")
    if identity == "Unavailable" or status in {"UNAVAILABLE", "ERROR"}:
        return _unavailable_explanation("identity", values, "ALIGNMENT IDENTITY — Unavailable")
    return Explanation(
        kind="identity",
        title=f"ALIGNMENT IDENTITY — {identity}",
        status=status,
        plain_meaning=f"Identical residues occupy {identity} of the columns in this pairwise alignment.",
        why_this_result="Identity is counted on the aligned strings returned by the declared method (Needleman-Wunsch or Smith-Waterman).",
        interpretation="Identity is not BLAST bit score and not evolutionary distance.",
        limitations=_limitations_for_status(
            status,
            "This identity is for this pair and this alignment. It does not prove homology beyond the computed alignment.",
        ),
        method=method,
        source=source,
    )


def _evalue(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values)
    evalue = display_text(values.get("evalue") or values.get("E-value"))
    method, source = _method_source(values, method="BLAST expect value from the search engine", source="BLAST")
    return Explanation(
        kind="evalue",
        title=f"E-VALUE — {evalue}",
        status=status,
        plain_meaning=(
            f"The BLAST expect value (E-value) for this hit is {evalue}. "
            "It estimates how many alignments of this score would be expected by chance in a search of this database size."
        ),
        why_this_result="The E-value comes from the BLAST engine that produced this hit, not from HelixScope rewriting the statistic.",
        interpretation="A smaller E-value means a less chance-like score under the BLAST model. It is not a probability that the sequences share a biological function.",
        limitations=_limitations_for_status(
            status,
            "E-value depends on database size and search parameters. It is not a p-value of function.",
        ),
        method=method,
        source=source,
    )


def _mit(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "HEURISTIC")
    mit = display_number(values.get("mit_score") or values.get("mit"), digits=2)
    method, source = _method_source(values, method="local MIT-like specificity proxy as implemented", source="HelixScope")
    return Explanation(
        kind="mit",
        title=f"MIT SPECIFICITY PROXY — {mit}",
        status=status,
        plain_meaning=(
            f"The MIT-like specificity proxy for this guide is {mit}. "
            "It is a local sequence statistic, not a genome-wide MIT specificity score from a full reference."
        ),
        why_this_result="HelixScope reports the MIT-related number that the CRISPR module actually computed, or Unavailable.",
        interpretation="Do not treat this as a published genome-wide MIT score unless the result says a full reference was scanned.",
        limitations=_limitations_for_status(
            status,
            "A sequence-local or test-reference MIT proxy is not CFD, not Azimuth, and not a guarantee that the guide is safe.",
        ),
        method=method,
        source=source,
    )


def _cfd(values: Mapping[str, Any]) -> Explanation:
    status = _status_of(values, "HEURISTIC")
    cfd = display_number(values.get("cfd_score") or values.get("cfd"), digits=3)
    method, source = _method_source(values, method="local CFD-like specificity proxy as implemented", source="HelixScope")
    return Explanation(
        kind="cfd",
        title=f"CFD SPECIFICITY PROXY — {cfd}",
        status=status,
        plain_meaning=(
            f"The CFD-like specificity proxy for this guide is {cfd}. "
            "It is not a genome-wide Cutting Frequency Determination score unless a real genome scan produced it."
        ),
        why_this_result="The CRISPR module returns CFD-related numbers only when that implementation actually ran.",
        interpretation="Unavailable means the published CFD model or a genome-wide scan was not applied. That is not a CFD of 0.",
        limitations=_limitations_for_status(
            status,
            "Do not conclude that a guide is clinically safe from this proxy. Local scan is not genome-wide.",
        ),
        method=method,
        source=source,
    )


BUILDERS: dict[str, Callable[[Mapping[str, Any]], Explanation]] = {
    "dna": _dna,
    "dna_tm": _dna_tm,
    "dna_orf": _dna_orf,
    "rna": _rna,
    "rna_fold": _rna_fold,
    "protein": _protein,
    "alignment": _alignment,
    "motif": _motif,
    "ncbi": _ncbi,
    "blast": _blast,
    "msa": _msa,
    "phylogeny": _phylogeny,
    "crispr": _crispr,
    "variant": _variant,
    "structure": _structure,
    "compare": _compare,
    "compare_rmsd_block": _compare_rmsd_block,
    "evidence": _evidence,
    "evolution": _evolution,
    "sasa": _sasa,
    "contacts": _contacts,
    "entropy": _entropy,
    "gravy": _gravy,
    "pi": _pi,
    "instability": _instability,
    "rscu": _rscu,
    "enc": _enc,
    "cai": _cai,
    "identity": _identity_metric,
    "evalue": _evalue,
    "mit": _mit,
    "cfd": _cfd,
}


def explain(kind: str, values: Optional[Mapping[str, Any]] = None) -> Explanation:
    """Constroi a interpretacao deterministica de um resultado existente.

    Args:
        kind: Chave em BUILDERS (dna, rna_fold, compare, ...).
        values: Campos ja calculados. None e tratado como resultado indisponivel.

    Returns:
        Explanation imutavel. Valores ausentes aparecem como Unavailable.

    Raises:
        ValueError: Se kind nao for suportado.

    Nota biologica:
        O texto nunca e mais confiante que o status. PREDICTED nao e
        EXPERIMENTAL. ClinVar nao e diagnostico. Conservacao nao e funcao.
    """
    key = str(kind or "").strip().lower()
    payload: Mapping[str, Any] = values if isinstance(values, Mapping) else {}
    builder = BUILDERS.get(key)
    if builder is None:
        raise ValueError(f"Explanation kind '{kind}' is not defined.")
    record = builder(payload)
    _assert_status_honesty(record)
    return record


def explain_metrics(domain: str, values: Optional[Mapping[str, Any]] = None) -> dict[str, dict]:
    """Build per-metric two-layer explanations from an already computed result.

    Args:
        domain: dna, rna, protein, alignment, motif, blast, crispr, variant, ncbi, msa, phylogeny.
        values: Domain result fields. Missing keys are omitted, never filled with 0.

    Returns:
        Mapping of metric id to explanation dict (plain_meaning + method/limitations).

    Raises:
        Nenhum. Unknown domains return an empty mapping.

    Nota biologica:
        Each metric keeps the method that produced it. Unavailable stays unavailable.
    """
    payload: Mapping[str, Any] = values if isinstance(values, Mapping) else {}
    plans = DOMAIN_METRIC_PLANS.get(str(domain or "").strip().lower(), ())
    out: dict[str, dict] = {}
    for metric_id, kind, extractor in plans:
        extracted = extractor(payload)
        if extracted is None:
            continue
        if not isinstance(extracted, Mapping):
            continue
        record = explain(kind, extracted)
        out[metric_id] = explanation_as_dict(record)
    return out


def _nested(values: Mapping[str, Any], key: str) -> Mapping[str, Any] | None:
    item = values.get(key)
    return item if isinstance(item, Mapping) else None


def _if_present(values: Mapping[str, Any], *keys: str) -> Mapping[str, Any] | None:
    if not keys:
        return values if values else None
    if any(key in values and values.get(key) is not None for key in keys):
        return values
    return None


def _dna_tm_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    report = _nested(values, "tm_report")
    return report if report else None


def _dna_orf_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    orfs = values.get("orfs")
    if not isinstance(orfs, list) and "orfs_raw_count" not in values:
        return None
    n = len(orfs) if isinstance(orfs, list) else values.get("orfs_raw_count")
    return {
        "status": "PREDICTED",
        "n_orfs": n,
        "n_codons": values.get("n_codons"),
        "method": "ATG start, in-frame stop, minimum length, optional non-overlap filter",
    }


def _fold_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    fold = _nested(values, "fold")
    return fold if fold else None


def _rscu_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    if values.get("rscu") is None:
        return None
    cai = _nested(values, "cai") or {}
    enc = _nested(values, "enc") or {}
    return {
        "status": str(cai.get("status") or enc.get("status") or values.get("codon_metrics_status") or "COMPUTED"),
        "n_codons": cai.get("n_codons") or enc.get("n_codons"),
        "method": "RSCU from the codon table of this RNA string",
    }


def _enc_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    enc = _nested(values, "enc")
    if not enc:
        return None
    return {
        "status": enc.get("status") or "COMPUTED",
        "enc": enc.get("enc") or enc.get("value") or enc.get("ENC"),
        "method": enc.get("method") or "effective number of codons",
    }


def _cai_values(values: Mapping[str, Any]) -> Mapping[str, Any] | None:
    cai = _nested(values, "cai")
    if not cai:
        return None
    return {
        "status": cai.get("status") or "COMPUTED",
        "cai": cai.get("cai") or cai.get("value") or cai.get("CAI"),
        "method": cai.get("method") or "codon adaptation index",
    }


DOMAIN_METRIC_PLANS: dict[str, tuple[tuple[str, str, Callable[[Mapping[str, Any]], Mapping[str, Any] | None]], ...]] = {
    "dna": (
        ("gc_content", "dna", _if_present),
        ("tm", "dna_tm", _dna_tm_values),
        ("orfs", "dna_orf", _dna_orf_values),
        ("shannon_entropy", "entropy", lambda values: _if_present(values, "shannon_entropy")),
    ),
    "rna": (
        ("gc_content", "rna", _if_present),
        ("shannon_entropy", "entropy", lambda values: _if_present(values, "shannon_entropy")),
        ("rscu", "rscu", _rscu_values),
        ("enc", "enc", _enc_values),
        ("cai", "cai", _cai_values),
        ("fold", "rna_fold", _fold_values),
    ),
    "protein": (
        ("physicochemical", "protein", _if_present),
        ("gravy", "gravy", lambda values: _if_present(values, "gravy")),
        ("isoelectric_point", "pi", lambda values: _if_present(values, "isoelectric_point")),
        ("instability_index", "instability", lambda values: _if_present(values, "instability_index")),
    ),
    "alignment": (("identity_pct", "identity", lambda values: _if_present(values, "identity_pct")),),
    "motif": (("hits", "motif", _if_present),),
    "blast": (("evalue", "evalue", lambda values: _if_present(values, "evalue") or values),),
    "crispr": (
        ("guides", "crispr", _if_present),
        ("mit", "mit", lambda values: _if_present(values, "mit_score", "mit")),
        ("cfd", "cfd", lambda values: _if_present(values, "cfd_score", "cfd")),
    ),
    "variant": (("identity", "variant", _if_present),),
    "ncbi": (("record", "ncbi", _if_present),),
    "msa": (("alignment", "msa", _if_present),),
    "phylogeny": (("tree", "phylogeny", _if_present),),
}


def _assert_status_honesty(record: Explanation) -> None:
    blob = " ".join(
        [
            record.title,
            record.plain_meaning,
            record.why_this_result,
            record.interpretation,
            record.limitations,
        ]
    ).lower()
    for banned in DIAGNOSIS_MARKERS:
        if banned in blob:
            raise ValueError(f"Explanation contains a prohibited diagnostic claim: {banned}")
    token = record.status.upper()
    if token in {"PREDICTED", "ILLUSTRATIVE", "HEURISTIC"}:
        for marker in EXPERIMENTAL_CLAIM_MARKERS:
            if marker in blob and "not" not in blob:
                raise ValueError("Non-experimental explanation claimed experimental coordinates.")
        if token == "PREDICTED" and "experimental structure" in blob and "not" not in blob:
            raise ValueError("PREDICTED explanation claimed an experimental structure.")
    if token in {"UNAVAILABLE", "NOT_INSTALLED", "RESOURCE_LIMIT"}:
        if "successfully computed" in blob or "the computation succeeded" in blob:
            raise ValueError("Unavailable explanation implied success.")


def explanation_as_dict(record: Explanation) -> dict:
    """Dict JSON-safe para testes e export.

    Args:
        record: Explanation.

    Returns:
        Dict com os campos publicos.

    Raises:
        Nenhum.
    """
    return {
        "kind": record.kind,
        "title": record.title,
        "status": record.status,
        "plain_meaning": record.plain_meaning,
        "why_this_result": record.why_this_result,
        "interpretation": record.interpretation,
        "limitations": record.limitations,
        "method": record.method,
        "source": record.source,
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
