"""Phase 19 information architecture: module registry and sidebar routing.

UI/navigation metadata only. Scientific results must not be stored here.
Frontend code in this module does not recompute GC, RMSD, alignments, or scores.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import streamlit as st

ACTIVE_MODULE_KEY: str = "helix_active_module"
MODULE_HISTORY_KEY: str = "helix_module_history"
MODULE_FORWARD_KEY: str = "helix_module_forward"
DEFAULT_MODULE_ID: str = "overview"
NAV_BUTTON_PREFIX: str = "nav_mod_"
MAX_MODULE_HISTORY: int = 24

GROUP_WORKSPACE: str = "Workspace"
GROUP_ANALYSIS: str = "Analysis"
GROUP_DISCOVERY: str = "Discovery"
GROUP_DATA: str = "Data"
GROUP_EVOLUTION: str = "Evolution"
GROUP_STRUCTURE: str = "Structure"

GROUP_ORDER: tuple[str, ...] = (
    GROUP_WORKSPACE,
    GROUP_ANALYSIS,
    GROUP_DISCOVERY,
    GROUP_DATA,
    GROUP_EVOLUTION,
    GROUP_STRUCTURE,
)


@dataclass(frozen=True)
class ModuleSpec:
    """Navigation metadata for one HelixScope module.

    Attributes:
        id: Semantic identifier used in session state. Never a visual tab index.
        label: Sidebar label. Scientific names are preserved (not marketing).
        group: Sidebar section label from GROUP_ORDER.
        title: Contextual header title.
        summary: One-line description shown in the compact header.
        method: How the calculation or retrieval works.
        scope: What the result means and does not mean.
        privacy: What data may leave the machine, if anything.
        provenance: Where results come from when a run exists.
        privacy_expanded: When True, the privacy expander starts open.
        renderer: Dispatch key resolved in app.py. Not a scientific object.
    """

    id: str
    label: str
    group: str
    title: str
    summary: str
    method: str
    scope: str
    privacy: str = ""
    provenance: str = ""
    privacy_expanded: bool = False
    renderer: str = ""


MODULES: tuple[ModuleSpec, ...] = (
    ModuleSpec(
        id="overview",
        label="Overview",
        group=GROUP_WORKSPACE,
        title="Workspace overview",
        summary="Active scientific objects in this session and engine status. Not a marketing page.",
        method="This page reads Streamlit session state and the optional-tool registry snapshot. It does not recompute GC, alignments, RMSD, or variant consequences.",
        scope="Presence of an object means it was loaded or computed in this session. Absence is not a scientific negative result.",
        privacy="Overview does not send sequences or identifiers off this machine.",
        provenance="Engine rows come from modules.tool_registry.collect_tool_snapshot (PATH/env detection, no disk scan).",
        renderer="overview",
    ),
    ModuleSpec(
        id="settings",
        label="Scientific Engines",
        group=GROUP_WORKSPACE,
        title="Scientific engines",
        summary="Local and remote engines that have a real validation on this machine. Unvalidated tools are omitted.",
        method="Status is assembled from existing detectors (ViennaRNA, FastTree, DSSP, IQ-TREE, Cas-OFFinder, BLAST+, US-align, RCSB Alignment API). HelixScope does not reimplement missing binaries.",
        scope="DETECTED is not LIVE_VALIDATED. REMOTE_VALIDATED is not a local engine. TEST REFERENCE is not a public genome. REFERENCE READY is not genome-wide completed.",
        privacy="This panel does not submit jobs. Individual modules disclose network use when the user runs them.",
        provenance="Snapshot identity follows modules.tool_registry and engine_validation records already stored on this machine.",
        renderer="settings",
    ),
    ModuleSpec(
        id="dna",
        label="DNA",
        group=GROUP_ANALYSIS,
        title="DNA analysis",
        summary="Independent of the RNA module.",
        method="Validate DNA, compute GC, Tm, restriction sites, ORFs and CpG islands. Displayed numbers come from modules.dna_analysis, not from this layout.",
        scope="This module does not transcribe DNA into RNA. DNA is rejected, not transcribed, in the RNA module. Illustrative 3D is not an experimental structure.",
        privacy="DNA analysis is local. Optional RCSB/AlphaFold retrieval in related structure sections is on-demand and disclosed there.",
        provenance="COMPUTED metrics carry sequence_hash and software_version from modules.provenance.",
        renderer="dna",
    ),
    ModuleSpec(
        id="rna",
        label="RNA",
        group=GROUP_ANALYSIS,
        title="RNA analysis",
        summary="Independent of the DNA module. Analyzes RNA (A, U, C, G).",
        method="DNA is rejected, not transcribed. Secondary structure is ViennaRNA MFE when installed, otherwise UNAVAILABLE. HelixScope does not invent a fold.",
        scope="PREDICTED fold is not EXPERIMENTAL structure. Codon-usage sections require a resolved coding region.",
        privacy="Folding, when available, runs locally. Structure retrieval is on-demand and disclosed in the structure section.",
        provenance="Fold records include tool version when ViennaRNA is present; otherwise status is UNAVAILABLE.",
        renderer="rna",
    ),
    ModuleSpec(
        id="protein",
        label="Protein",
        group=GROUP_ANALYSIS,
        title="Protein analysis",
        summary="Translate and profile amino acid composition, pI, hydrophobicity and stability.",
        method="Experimental PDB and predicted AlphaFold records are retrieved on demand. Local secondary-structure assignment stays UNAVAILABLE unless DSSP/STRIDE is installed.",
        scope="PREDICTED != EXPERIMENTAL. Loading a sequence does not imply a structure. SASA is a local geometric calculation when run, not a wet-lab measurement.",
        privacy="RCSB and AlphaFold requests leave this machine only when the user searches or loads a structure.",
        provenance="Structure kind and source are stored on the result object. The UI must not relabel PREDICTED as EXPERIMENTAL.",
        renderer="protein",
    ),
    ModuleSpec(
        id="alignment",
        label="Alignment",
        group=GROUP_ANALYSIS,
        title="Sequence alignment",
        summary="Global and local pairwise alignment with optional dot plot visualization.",
        method="Global is Needleman-Wunsch. Local is Smith-Waterman. Optional translation uses frame +1 and stops at the first stop, then aligns as protein.",
        scope="PAIRWISE ALIGNMENT is not MSA. SMITH-WATERMAN is not BLAST. Pairwise identity is not a phylogenetic tree.",
        privacy="Pairwise alignment is local. It does not submit sequences to NCBI.",
        provenance="Alignment metrics come from modules.alignment result objects stored in session state.",
        renderer="alignment",
    ),
    ModuleSpec(
        id="motif",
        label="Motif Search",
        group=GROUP_DISCOVERY,
        title="Motif search",
        summary="Exact and IUPAC motif search on DNA, RNA or protein.",
        method="Exact or IUPAC pattern search on the sequence you provide. Hits are not PWM scores, Pfam domains or experimentally validated sites.",
        scope="Motif hits are pattern matches, not annotated domains. A hit is not a proven binding site.",
        privacy="Motif search is local. Nothing is sent to motif databases.",
        provenance="Hit tables come from modules.motif_search for the sequence and patterns in this session.",
        renderer="motif",
    ),
    ModuleSpec(
        id="crispr",
        label="CRISPR Design",
        group=GROUP_DISCOVERY,
        title="CRISPR design",
        summary="Scan Cas sites on the provided sequence. Efficiency is a labeled heuristic, not published Rule Set 2.",
        method="SpCas9 (NGG), SaCas9, Cas12a/Cas12b, Cas13, base editors and prime editor use the metrics adapted in modules.crispr. On-target Doench Rule Set 2 / DeepHF and genome-wide CFD/MIT are declared unavailable offline via crispr.model_availability().",
        scope="Off-target in a provided FASTA is explicit search, not genome-wide. TEST REFERENCE is not a public genome. REFERENCE READY is not GENOME-WIDE COMPLETED. DETECTED is not LIVE_VALIDATED.",
        privacy="Cas-OFFinder runs locally when installed. Genome FASTA stays on this machine. No genome-wide MIT/CFD is claimed without COMPLETED_FULL_REFERENCE.",
        provenance="Guide tables and reference-search packs keep status, reference identity, and software_version from the CRISPR modules.",
        renderer="crispr",
    ),
    ModuleSpec(
        id="variant",
        label="Variant Explorer",
        group=GROUP_DISCOVERY,
        title="Variant Explorer",
        summary="Variant identity with mandatory assembly. HelixScope predicts nothing here.",
        method="Consequences come from Ensembl VEP, clinical claims are ClinVar submissions retrieved verbatim, domains come from InterPro. Genome coordinates depend on the selected assembly.",
        scope="HelixScope records variant identity and retrieves annotation from official services. It does not predict pathogenicity, does not compute a variant effect and does not build mutant structures. A genomic coordinate without an assembly is refused, because chr17:43093557 C>G means different bases in GRCh37 and GRCh38. CLINVAR EVIDENCE is not a HelixScope diagnosis. CONSERVATION is not PATHOGENICITY.",
        privacy="Nothing is sent until you press Annotate. Disabling a source means its answer is UNKNOWN, not absent. NCBI requires a contact address for E-utilities.",
        provenance="Each layer (identity, VEP, ClinVar, InterPro/UniProt) keeps its own status on the variant result object.",
        privacy_expanded=True,
        renderer="variant",
    ),
    ModuleSpec(
        id="ncbi",
        label="NCBI Fetch",
        group=GROUP_DATA,
        title="NCBI Fetch",
        summary="Search and retrieve real Entrez records. Homology search is NCBI BLAST, not this module.",
        method="NCBI Entrez requires a valid contact email for every request. Numeric IDs are treated as GeneIDs. Names of organisms are not accessions.",
        scope="Retrieved records are RETRIEVED, not computed. Search with zero records is not BLAST NO_HITS. Pairwise alignment is not BLAST.",
        privacy="Search queries, accessions, and the contact email leave this machine for NCBI Entrez when you search or fetch.",
        provenance="Each payload stores source, retrieved_at_utc, and NCBI count/offset metadata.",
        privacy_expanded=True,
        renderer="ncbi",
    ),
    ModuleSpec(
        id="blast",
        label="NCBI BLAST",
        group=GROUP_DATA,
        title="NCBI BLAST",
        summary="Real NCBI BLAST Common URL API. No local scores, no fake E-values.",
        method="Remote NCBI BLAST is never silently replaced by BLAST+ local. Sequence similarity is not confirmed biological function.",
        scope="SMITH-WATERMAN is not BLAST. PAIRWISE ALIGNMENT is not BLAST. BLAST hits are not an MSA.",
        privacy="The query sequence, program, database, and contact email are sent to NCBI when you submit. BLAST+ local, when used, stays on this machine and requires an explicit local database.",
        provenance="Job identity, RID, database, and cache keys follow modules.blast_search, including software_version.",
        privacy_expanded=True,
        renderer="blast",
    ),
    ModuleSpec(
        id="references",
        label="References",
        group=GROUP_DATA,
        title="Reference genomes",
        summary="Local catalog status for assemblies HelixScope already knows how to install. Not a download page that bypasses resource admission.",
        method="Rows come from modules.genome_store.list_local_references plus resource-admission memory facts. HelixScope does not scan the disk for arbitrary FASTA files here.",
        scope="TEST REFERENCE is not a public genome. NOT_INSTALLED is not READY. READY is not COMPLETED_FULL_REFERENCE. GRCh38.p14 is GCF_000001405.40. RESOURCE_LIMIT is a scientific/machine fact. No genome-wide MIT or CFD is claimed.",
        privacy="This page does not download genomes. Install and Cas-OFFinder search stay in CRISPR Design / genome jobs.",
        provenance="Manifest software_version and file_sha256 are shown when a local manifest exists.",
        renderer="references",
    ),
    ModuleSpec(
        id="msa",
        label="MSA",
        group=GROUP_EVOLUTION,
        title="Multiple sequence alignment",
        summary="Real Clustal Omega (EMBL-EBI) or local MAFFT/MUSCLE/clustalo when installed.",
        method="BLAST hits are not an MSA. Consensus is not a biological sequence. Conservation is not function. Phylogenetic inference is a separate explicit action on a completed MSA, not this alignment step. IQ-TREE and FastTree stay UNAVAILABLE unless those official executables are installed; HelixScope does not reimplement them.",
        scope="PAIRWISE ALIGNMENT is not MSA. MSA is not phylogeny. TAXONOMY is not phylogeny.",
        privacy="EMBL-EBI jobs are submitted and polled from this session when that backend is chosen. Local clustalo/MAFFT/MUSCLE, when installed, run in-process.",
        provenance="MSA result objects store tool, tool_version, method, alignment_hash, and input hashes.",
        privacy_expanded=True,
        renderer="msa",
    ),
    ModuleSpec(
        id="phylogeny",
        label="Phylogeny",
        group=GROUP_EVOLUTION,
        title="Phylogeny",
        summary="Explicit tree inference on a completed MSA. Not automatic. Not a taxonomic authority tree.",
        method="Neighbor-Joining and UPGMA are local distance methods. IQ-TREE ML, ModelFinder, UFBoot and SH-aLRT run only when the official IQ-TREE executable is installed. FastTree stays NOT_INSTALLED unless that binary is present.",
        scope="MSA is not phylogeny. TAXONOMY is not phylogeny. A drawing root is not a confirmed biological root unless the method says otherwise. Branch lengths and support values are not rewritten by this layout.",
        privacy="Local NJ/UPGMA stay on this machine. IQ-TREE runs locally when installed. No sequences are sent to a tree web service from this page.",
        provenance="Tree objects store method, model, alignment_hash, n_leaves, and software_version. Status remains COMPUTED or UNAVAILABLE/ERROR as produced by modules.phylogeny.",
        renderer="phylogeny",
    ),
    ModuleSpec(
        id="structure_3d",
        label="3D Viewer",
        group=GROUP_STRUCTURE,
        title="3D viewer",
        summary="Session inventory of loaded structure representations. Rendering engines remain the validated Plotly/WebGL viewers in DNA, RNA, Protein, CRISPR and Compare.",
        method="This page lists structure objects already stored in session state. It does not fetch PDB files, does not superpose, and does not invent coordinates.",
        scope="ILLUSTRATIVE is not EXPERIMENTAL. PREDICTED is not EXPERIMENTAL. STRUCTURAL CONTACT is not proven function. A listed object is not automatically the active Compare pair.",
        privacy="No new remote fetch is triggered by opening this page.",
        provenance="Each listed object keeps the kind/status already stored on it (EXPERIMENTAL, PREDICTED, ILLUSTRATIVE, STALE, UNAVAILABLE).",
        renderer="structure_3d",
    ),
    ModuleSpec(
        id="compare",
        label="Structure Compare",
        group=GROUP_STRUCTURE,
        title="Compare",
        summary="Comparative science: structure alignment via the RCSB Alignment API, plus variant, protein, guide, evolution and evidence modes. Similarity is not causality.",
        method="Protein structure alignment is the RCSB PDB Structure Alignment API (remote). Modes Structures, Variants, Proteins, Guides, Evolution and Evidence stay inside this workspace. HelixScope does not invent RMSD, TM-score, mutant coordinates, pathogenicity ranks, dN/dS, or genome-wide CRISPR scores. Block RMSD and global RMSD are distinct fields.",
        scope="STRUCTURAL SIMILARITY is not causality. STRUCTURAL CONTACT is not proven function. CONSERVATION is not PATHOGENICITY. MIT/CFD only compare in the same declared scope. Without GRCh38 COMPLETED_FULL_REFERENCE there is no genome-wide comparison. REMOTE_VALIDATED is not a local engine.",
        privacy="Structure IDs and alignment method are sent to alignment.rcsb.org when you align. Coordinate download for superposition is explicit. Original coordinates are not overwritten; superposition uses a transformed copy.",
        provenance="Alignment records store method, pair identity, RMSD fields, TM-score, coverage, and software_version from modules.rcsb_alignment / structure_alignment.",
        privacy_expanded=True,
        renderer="compare",
    ),
)

MODULE_BY_ID: dict[str, ModuleSpec] = {item.id: item for item in MODULES}


def module_ids() -> tuple[str, ...]:
    """Return registry ids in sidebar order.

    Returns:
        Semantic module identifiers.

    Raises:
        Nenhum.
    """
    return tuple(item.id for item in MODULES)


def get_module(module_id: str) -> ModuleSpec:
    """Return one module spec.

    Args:
        module_id: Semantic id.

    Returns:
        Frozen ModuleSpec.

    Raises:
        KeyError: Unknown id.
    """
    return MODULE_BY_ID[module_id]


def grouped_modules() -> list[tuple[str, tuple[ModuleSpec, ...]]]:
    """Group registry entries for sidebar rendering.

    Returns:
        Pairs (group label, modules in that group). Empty groups are omitted.

    Raises:
        Nenhum.
    """
    rows: list[tuple[str, tuple[ModuleSpec, ...]]] = []
    for group in GROUP_ORDER:
        items = tuple(item for item in MODULES if item.group == group)
        if items:
            rows.append((group, items))
    return rows


SCIENTIFIC_SESSION_KEYS: frozenset[str] = frozenset(
    {
        ACTIVE_MODULE_KEY,
        "dna_input",
        "rna_input",
        "rna_raw_text",
        "align_result",
        "msa_result",
        "msa_collection",
        "msa_job",
        "msa_cache_key",
        "msa_column",
        "phylo_result",
        "blast_result",
        "blast_job",
        "crispr_result",
        "crispr_complex",
        "crispr_reference_search",
        "crispr_draft_seq",
        "crispr_draft_name",
        "compare_structures",
        "compare_variants",
        "compare_export",
        "prot_struct_result",
        "prot_struct_hits",
        "dna_3d_representation",
        "dna_struct_result",
        "rna_3d_representation",
        "molecule_selection",
        "workspace_dna",
        "workspace_rna",
        "workspace_protein",
        "protein_input",
        "ncbi_search",
        "ncbi_record",
    }
)


def persist_orphaned_widget_keys(
    extra_skip_prefixes: Iterable[str] = (),
) -> None:
    """Re-save scientific result objects only.

    Streamlit 1.58 forbids assigning button and file_uploader keys through
    session_state. Those widget types are never copied here. Result objects
    such as dna_input and msa_result are not widgets; re-saving them documents
    the session boundary and does not recompute science.

    Args:
        extra_skip_prefixes: Unused; kept for call-site compatibility.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    del extra_skip_prefixes
    for key in SCIENTIFIC_SESSION_KEYS:
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]


def active_module_id() -> str:
    """Return the semantic active module, defaulting to overview.

    Returns:
        A registry id.

    Raises:
        Nenhum.
    """
    current = str(st.session_state.get(ACTIVE_MODULE_KEY) or "")
    if current not in MODULE_BY_ID:
        current = DEFAULT_MODULE_ID
        st.session_state[ACTIVE_MODULE_KEY] = current
    return current


def _known_modules() -> set[str]:
    return set(MODULE_BY_ID)


def set_active_module(module_id: str) -> None:
    """Abre um modulo e guarda o anterior para Voltar.

    Args:
        module_id: Semantic id.

    Returns:
        None.

    Raises:
        Nenhum. Unknown ids are ignored.
    """
    from modules.session_history import visit_module

    if module_id not in MODULE_BY_ID:
        return
    current = active_module_id()
    back, forward, target = visit_module(
        list(st.session_state.get(MODULE_HISTORY_KEY) or []),
        list(st.session_state.get(MODULE_FORWARD_KEY) or []),
        current,
        module_id,
        known=_known_modules(),
    )
    st.session_state[MODULE_HISTORY_KEY] = back
    st.session_state[MODULE_FORWARD_KEY] = forward
    st.session_state[ACTIVE_MODULE_KEY] = target


def can_return_to_previous_module() -> bool:
    """True when Voltar tem destino.

    Returns:
        Se a pilha de volta nao esta vazia.

    Raises:
        Nenhum.
    """
    return bool(st.session_state.get(MODULE_HISTORY_KEY) or [])


def can_advance_to_next_module() -> bool:
    """True when Avancar tem destino.

    Returns:
        Se a pilha de avanco nao esta vazia.

    Raises:
        Nenhum.
    """
    return bool(st.session_state.get(MODULE_FORWARD_KEY) or [])


def return_to_previous_module() -> str:
    """Restaura o modulo visitado imediatamente antes.

    Returns:
        O modulo agora ativo.

    Raises:
        Nenhum.
    """
    from modules.session_history import go_back

    back, forward, target = go_back(
        list(st.session_state.get(MODULE_HISTORY_KEY) or []),
        list(st.session_state.get(MODULE_FORWARD_KEY) or []),
        active_module_id(),
        known=_known_modules(),
    )
    st.session_state[MODULE_HISTORY_KEY] = back
    st.session_state[MODULE_FORWARD_KEY] = forward
    st.session_state[ACTIVE_MODULE_KEY] = target
    return target


def advance_to_next_module() -> str:
    """Restaura o modulo abandonado pelo ultimo Voltar.

    Returns:
        O modulo agora ativo.

    Raises:
        Nenhum.
    """
    from modules.session_history import go_forward

    back, forward, target = go_forward(
        list(st.session_state.get(MODULE_HISTORY_KEY) or []),
        list(st.session_state.get(MODULE_FORWARD_KEY) or []),
        active_module_id(),
        known=_known_modules(),
    )
    st.session_state[MODULE_HISTORY_KEY] = back
    st.session_state[MODULE_FORWARD_KEY] = forward
    st.session_state[ACTIVE_MODULE_KEY] = target
    return target


def render_back_control() -> None:
    """Setas de voltar e avancar. Permanecem visiveis mesmo sem destino.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    back_col, forward_col = st.columns(2, gap="small")
    with back_col:
        if st.button(
            "← Voltar",
            key="helix_back",
            disabled=not can_return_to_previous_module(),
            help="Return to the previous module",
        ):
            return_to_previous_module()
            st.rerun()
    with forward_col:
        if st.button(
            "→ Avançar",
            key="helix_forward",
            disabled=not can_advance_to_next_module(),
            help="Return to the module left by Voltar",
        ):
            advance_to_next_module()
            st.rerun()


def render_sidebar_navigation() -> str:
    """Render grouped sidebar navigation and return the active module id.

    Args:
        Nenhum.

    Returns:
        Active semantic module id.

    Raises:
        Nenhum.
    """
    current = active_module_id()
    with st.sidebar:
        st.markdown(
            '<p class="helix-nav-brand">HelixScope</p>'
            '<p class="helix-nav-tag">Scientific workstation</p>',
            unsafe_allow_html=True,
        )
        for group, items in grouped_modules():
            st.caption(group)
            for spec in items:
                is_active = spec.id == current
                clicked = st.button(
                    spec.label,
                    key=f"{NAV_BUTTON_PREFIX}{spec.id}",
                    type="primary" if is_active else "secondary",
                    width="stretch",
                )
                if clicked and not is_active:
                    set_active_module(spec.id)
                    st.rerun()
    return active_module_id()


def workspace_object_catalog() -> tuple[tuple[str, str, str], ...]:
    """Session keys the overview may list. Values are not computed here.

    Returns:
        Tuples of (session_key, label, module_id).

    Raises:
        Nenhum.
    """
    return (
        ("dna_input", "DNA sequence", "dna"),
        ("rna_input", "RNA sequence", "rna"),
        ("align_result", "Pairwise alignment", "alignment"),
        ("msa_result", "MSA", "msa"),
        ("phylo_result", "Phylogeny", "phylogeny"),
        ("blast_result", "BLAST result", "blast"),
        ("crispr_result", "CRISPR guides", "crispr"),
        ("variant_result", "Variant", "variant"),
        ("prot_struct_result", "Protein structure", "protein"),
        ("dna_3d_representation", "DNA 3D representation", "dna"),
        ("rna_3d_representation", "RNA 3D representation", "rna"),
        ("crispr_complex", "CRISPR complex", "crispr"),
        ("molecule_selection", "Molecule selection", "structure_3d"),
    )
