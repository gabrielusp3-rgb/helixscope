"""HelixScope scientific core: Streamlit-free public API.

Implementation remains in ``modules/`` (one science). This package is the
stable facade FastAPI and workers must import. Product version stays
``0.24.3-19`` in ``modules.provenance.HELIXSCOPE_VERSION``.
"""

from __future__ import annotations

from helixscope_core.alignment import pairwise_align
from helixscope_core.capabilities import get_capabilities
from helixscope_core.crispr import design_guides
from helixscope_core.dna import analyze_dna
from helixscope_core.evidence import evidence_conflict_pack
from helixscope_core.explain import explain_metrics, explain_result
from helixscope_core.motif import search_motifs
from helixscope_core.phylogeny import infer_tree
from helixscope_core.protein import analyze_protein, retrieve_uniprot_entry
from helixscope_core.references import list_references
from helixscope_core.rna import analyze_rna, transcribe_dna
from helixscope_core.variants import identify_variant

CORE_PACKAGE_VERSION: str = "0.24.3.19"
"""PEP 440 package metadata. Not a product release bump."""

__all__ = (
    "CORE_PACKAGE_VERSION",
    "analyze_dna",
    "analyze_rna",
    "analyze_protein",
    "retrieve_uniprot_entry",
    "transcribe_dna",
    "pairwise_align",
    "search_motifs",
    "explain_metrics",
    "explain_result",
    "infer_tree",
    "design_guides",
    "identify_variant",
    "evidence_conflict_pack",
    "get_capabilities",
    "list_references",
)
