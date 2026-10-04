"""Facade identity: one science implementation, re-exported not copied."""

from __future__ import annotations

from helixscope_core import alignment as core_alignment
from helixscope_core import crispr as core_crispr
from helixscope_core import dna as core_dna
from helixscope_core import explain as core_explain
from helixscope_core import motif as core_motif
from helixscope_core import msa as core_msa
from helixscope_core import phylogeny as core_phylo
from helixscope_core import protein as core_protein
from helixscope_core import rna as core_rna
from helixscope_core import variants as core_variants
from modules import (
    alignment,
    crispr,
    dna_analysis,
    explain,
    motif_search,
    msa,
    phylogeny,
    protein_analysis,
    rna_analysis,
    variant_core,
)


def test_reexported_callables_are_the_same_objects() -> None:
    assert core_dna.gc_content is dna_analysis.gc_content
    assert core_dna.find_orfs is dna_analysis.find_orfs
    assert core_rna.transcribe is rna_analysis.transcribe
    assert core_protein.physicochemical_report is protein_analysis.physicochemical_report
    assert core_alignment.pairwise_global is alignment.pairwise_global
    assert core_alignment.pairwise_local is alignment.pairwise_local
    assert core_motif.find_motif is motif_search.find_motif
    assert core_explain.explain is explain.explain
    assert core_msa.import_prealigned_fasta is msa.import_prealigned_fasta
    assert core_msa.parse_unaligned_fasta is msa.parse_unaligned_fasta
    assert core_phylo.infer_phylogeny is phylogeny.infer_phylogeny
    assert core_crispr.find_guides is crispr.find_guides
    assert core_crispr.model_availability is crispr.model_availability
    assert core_variants.build_variant is variant_core.build_variant
