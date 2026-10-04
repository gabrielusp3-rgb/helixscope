"""Domain exceptions already used by modules/. One taxonomy, no new hierarchy."""

from __future__ import annotations

from modules.blast_search import BlastError
from modules.cds_mapping import CdsMappingError
from modules.clinvar_evidence import ClinVarError
from modules.comparative import ComparativeError
from modules.complex_structure import ComplexError
from modules.crispr import CrisprError
from modules.ensembl_vep import VepError
from modules.evidence_workspace import EvidenceError
from modules.genome_coordinates import CoordinateError
from modules.genome_download import GenomeDownloadError
from modules.genome_jobs import GenomeJobError
from modules.genome_store import GenomeStoreError
from modules.msa import MsaError
from modules.ncbi_fetch import NCBIQueryError
from modules.opencl_runtime import OpenCLError
from modules.phylogeny import PhylogenyError
from modules.protein_domains import DomainError
from modules.protein_structure import StructureError
from modules.rcsb_alignment import AlignmentApiError
from modules.rna_folding import FoldingError
from modules.taxonomy import TaxonomyError
from modules.thermodynamics import ThermodynamicsError
from modules.usalign import USAlignError
from modules.variant_core import VariantInputError
from modules.variant_explorer import VariantExplorerError

__all__ = (
    "AlignmentApiError",
    "BlastError",
    "CdsMappingError",
    "ClinVarError",
    "ComparativeError",
    "ComplexError",
    "CoordinateError",
    "CrisprError",
    "DomainError",
    "EvidenceError",
    "FoldingError",
    "GenomeDownloadError",
    "GenomeJobError",
    "GenomeStoreError",
    "MsaError",
    "NCBIQueryError",
    "OpenCLError",
    "PhylogenyError",
    "StructureError",
    "TaxonomyError",
    "ThermodynamicsError",
    "USAlignError",
    "VariantExplorerError",
    "VariantInputError",
    "VepError",
)
