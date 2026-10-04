"""Shared transport schemas. Result payloads stay extra-allow dicts after to_transport."""

from __future__ import annotations

from typing import Any, Optional

from pydantic import BaseModel, ConfigDict, Field

from helixscope_api import API_CONTRACT_VERSION, PRODUCT_VERSION


class ApiEnvelope(BaseModel):
    """Consistent success wrapper. Domain meaning stays inside result."""

    model_config = ConfigDict(extra="forbid")

    request_id: str
    api_version: str = API_CONTRACT_VERSION
    product_version: str = PRODUCT_VERSION
    scientific_status: Optional[str] = None
    execution_status: Optional[str] = None
    input_hash: Optional[str] = None
    warnings: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)
    result: dict[str, Any]


class ErrorDetail(BaseModel):
    code: str
    message: str
    details: dict[str, Any] = Field(default_factory=dict)
    request_id: str


class ErrorEnvelope(BaseModel):
    error: ErrorDetail


class SequenceRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str = Field(
        ...,
        min_length=1,
        max_length=250_000,
        examples=["ATGCATGCATGC"],
    )


class DnaAnalyzeRequest(SequenceRequest):
    min_orf_length: int = Field(150, ge=3, le=10_000)
    include_orfs: bool = False
    include_restriction: bool = False
    include_profiles: bool = False
    include_explanation: bool = True
    include_cpg: bool = False
    include_kmers: bool = False
    include_santalucia: bool = False
    drop_overlapping_orfs: bool = True
    profile_window: int = Field(100, ge=10, le=2000)
    profile_step: int = Field(50, ge=1, le=2000)
    kmer_k: int = Field(3, ge=1, le=5)
    kmer_top_n: int = Field(30, ge=1, le=256)
    santalucia_na_m: float = Field(0.05, ge=0.001, le=1.5)
    santalucia_oligo_nm: float = Field(250.0, ge=1.0, le=1.0e6)


class DnaHelix3dRequest(SequenceRequest):
    prefer: str = Field("illustrative", pattern="^(illustrative|experimental|predicted|auto)$")
    lod: str = Field("high", pattern="^(low|medium|high)$")
    region_start: int = Field(0, ge=0)
    region_end: Optional[int] = Field(default=None, ge=0)
    inspect_position: int = Field(0, ge=0)
    center_on_selected: bool = False


class RnaAnalyzeRequest(SequenceRequest):
    fold: bool = False
    include_codon_metrics: bool = True
    include_explanation: bool = True
    include_profiles: bool = False
    include_kmers: bool = False
    backend: str = "auto"
    profile_window: int = Field(100, ge=10, le=2000)
    profile_step: int = Field(50, ge=1, le=2000)
    kmer_k: int = Field(3, ge=1, le=5)


class RnaHelix3dRequest(SequenceRequest):
    prefer: str = Field("illustrative", pattern="^(illustrative|experimental|predicted|auto)$")
    lod: str = Field("high", pattern="^(low|medium|high)$")
    region_start: int = Field(0, ge=0)
    region_end: Optional[int] = Field(default=None, ge=0)
    inspect_position: int = Field(0, ge=0)
    center_on_selected: bool = False


class TranscribeRequest(SequenceRequest):
    pass


class ProteinAnalyzeRequest(SequenceRequest):
    ph: float = Field(7.0, ge=0.0, le=14.0)
    include_explanation: bool = True


class ProteinUniprotRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    accession: str = Field(..., min_length=1, max_length=32)
    include_domains: bool = False
    include_explanation: bool = True


class PairwiseAlignRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    seq1: str = Field(..., min_length=1, max_length=8_000)
    seq2: str = Field(..., min_length=1, max_length=8_000)
    mode: str = Field("global", pattern="^(global|local)$")
    include_explanation: bool = True
    include_dotplot: bool = False
    translate_nucleic: bool = False


class MotifSearchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str = Field(..., min_length=1, max_length=250_000)
    pattern: str = Field(..., min_length=1, max_length=500)
    molecule: str = Field("DNA", pattern="^(DNA|RNA|PROTEIN)$")
    include_explanation: bool = True


class ExplainRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: str = Field(..., min_length=1, max_length=80)
    values: dict[str, Any] = Field(default_factory=dict)


class NcbiFetchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    accession: str = Field(..., min_length=1, max_length=64)
    email: str = Field(..., min_length=3, max_length=254)
    db: str = Field("nucleotide", max_length=32)
    api_key: Optional[str] = Field(default=None, max_length=128)


class NcbiSearchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term: str = Field(..., min_length=1, max_length=200)
    email: str = Field(..., min_length=3, max_length=254)
    db: str = Field("nucleotide", max_length=32)
    retmax: int = Field(15, ge=1, le=20)
    api_key: Optional[str] = Field(default=None, max_length=128)


class CrisprGuidesRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str = Field(..., min_length=1, max_length=250_000)
    cas_system: str = Field("SpCas9", max_length=80)
    max_mismatches: int = Field(3, ge=0, le=8)
    run_off_target: bool = False
    include_explanation: bool = True


class VariantIdentifyRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    text: str = Field(..., min_length=1, max_length=500)
    assembly: str = Field(..., min_length=1, max_length=80)


class VariantExploreRequest(VariantIdentifyRequest):
    email: str = Field("", max_length=254)
    enable_vep: bool = False
    enable_clinvar: bool = False
    enable_domains: bool = False


class EvidenceItemIn(BaseModel):
    model_config = ConfigDict(extra="forbid")

    field: str
    value: Any = None
    source: str
    evidence_status: str = "RETRIEVED"
    retrieved_at_utc: str = ""
    identifier: str = ""
    mapping_status: str = ""
    note: str = ""


class EvidencePackRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: str = "variants"
    items: list[EvidenceItemIn]
    inputs: dict[str, Any] = Field(default_factory=dict)


class MsaPrealignedRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fasta: str = Field(..., min_length=1, max_length=2_000_000)
    format: str = Field("auto", pattern="^(auto|fasta|clustal|stockholm|phylip)$")


class PhylogenyInferRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fasta: str = Field(..., min_length=1, max_length=2_000_000)
    method: str = Field(..., min_length=1, max_length=40)
    distance_model: str = "p_distance"
    alignment_hash: Optional[str] = Field(default=None, max_length=128)


class PhylogenyTaxonomyRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    tree: dict[str, Any]
    email: str = Field(..., min_length=3, max_length=254)
    declared_organisms: dict[str, str] = Field(default_factory=dict)


class StructureFixtureRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field(..., pattern="^(1CRN|1BNA|1RNA)$")


class StructureMappedRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field(..., pattern="^(1CRN|1BNA|1RNA)$")
    query_sequence: str = Field(..., min_length=1, max_length=10_000)
    molecule: str = Field("PROTEIN", pattern="^(PROTEIN|DNA|RNA)$")


class StructureMappingRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequences_identical: bool
    polymer_index_mode: bool
    observed_coordinate_residues: int = Field(..., ge=0)
    polymer_length: int = Field(..., ge=0)


class StructureCoordinateFileRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field(..., pattern="^(1CRN|1BNA|1RNA)$")


class StructureComplexRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field(..., pattern="^(4UN3|4OO8)$")


class CompareParseRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    payload: dict[str, Any]


class CompareRemoteJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reference_entry: str = Field(..., min_length=4, max_length=80, pattern="^[0-9A-Za-z_-]+$")
    target_entry: str = Field(..., min_length=4, max_length=80, pattern="^[0-9A-Za-z_-]+$")
    reference_chain: str = Field("A", max_length=8)
    target_chain: str = Field("A", max_length=8)
    method: str = Field("tm-align", max_length=40)
    fetch_coordinates: bool = True


class CompareUSalignJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reference_entry: str = Field(..., pattern="^(1CRN|1BNA|1RNA)$")
    target_entry: str = Field(..., pattern="^(1CRN|1BNA|1RNA)$")
    reference_chain: str = Field("A", max_length=8)
    target_chain: str = Field("A", max_length=8)


class BlastLocalJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query: str = Field(..., min_length=1, max_length=8_000)
    program: str = Field("blastn", max_length=16)
    molecule: str = Field("DNA", pattern="^(DNA|RNA|PROTEIN)$")


class BlastRemoteJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query: str = Field(..., min_length=1, max_length=8_000)
    email: str = Field(..., min_length=3, max_length=254)
    program: str = Field("blastn", max_length=16)


class MsaEngineJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fasta: str = Field(..., min_length=1, max_length=2_000_000)
    email: str = Field(..., min_length=3, max_length=254)
    backend: str = Field("ebi_clustalo", max_length=40)


class CasOffinderJobRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    guide_sequence: str = Field(..., min_length=10, max_length=40)
    assembly_id: str = Field(..., min_length=1, max_length=80)
    pam: str = Field("NGG", max_length=16)
    cas_system: str = Field("SpCas9", min_length=1, max_length=80)


class JobAccepted(BaseModel):
    request_id: str
    job_id: str
    kind: str
    execution_status: str
    poll_url: str


class RnaTranslateRequest(SequenceRequest):
    include_explanation: bool = True


class ProteinStructureSearchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str = Field(..., min_length=1, max_length=10_000)
    molecule: str = Field("PROTEIN", pattern="^(PROTEIN|DNA|RNA)$")
    pdb_id: str = Field("", max_length=32)
    uniprot: str = Field("", max_length=32)
    alphafold_id: str = Field("", max_length=80)
    search_pdb_by_sequence: bool = False
    lookup_uniprot: bool = False
    lookup_alphafold: bool = False
    identity_cutoff: float = Field(0.9, ge=0.0, le=1.0)


class ProteinStructureLoadRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field(..., min_length=4, max_length=80, pattern="^[0-9A-Za-z_-]+$")
    sequence: str = Field("", max_length=10_000)
    molecule: str = Field("PROTEIN", pattern="^(PROTEIN|DNA|RNA)$")
    chain: str = Field("", max_length=8)
    bundled: bool = False


class ProteinStructureAnalyzeRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str = Field("1CRN", pattern="^(1CRN|1BNA|1RNA)$")
    run_sasa: bool = True
    run_dssp: bool = False
    run_dssp_remote: bool = False
    run_stride: bool = False
    run_edtsurf: bool = False


class ProteinAlphafoldRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    uniprot: str = Field(..., min_length=1, max_length=32)


class CompareVariantsRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    text_a: str = Field(..., min_length=1, max_length=500)
    assembly_a: str = Field(..., min_length=1, max_length=80)
    text_b: str = Field(..., min_length=1, max_length=500)
    assembly_b: str = Field(..., min_length=1, max_length=80)


class CompareProteinsRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence_a: str = Field(..., min_length=1, max_length=8_000)
    sequence_b: str = Field(..., min_length=1, max_length=8_000)
    identifier_a: str = Field("protein_a", max_length=80)
    identifier_b: str = Field("protein_b", max_length=80)


class CompareGuidesRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    guide_a: dict[str, Any]
    guide_b: dict[str, Any]


class CompareEvolutionRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    fasta: str = Field(..., min_length=1, max_length=2_000_000)
    group_ids: list[str] = Field(default_factory=list)
    column: int = Field(0, ge=0)
    member_id: str = Field("", max_length=80)


class CompareEvidenceRequest(CompareVariantsRequest):
    """Evidence mode reuses two identified variant identities."""


class ReferenceAdmissionRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    assembly_id: str = Field(..., min_length=1, max_length=80)


class ReferenceInstallRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    assembly_id: str = Field(..., min_length=1, max_length=80)
    confirm: bool = False
