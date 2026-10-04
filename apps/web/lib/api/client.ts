import type { components } from "@/lib/api/generated/schema";
import { classifyHttpError, HelixApiError } from "@/lib/api/errors";
import { asRecord } from "@/lib/api/numeric";

export type ApiEnvelope = components["schemas"]["ApiEnvelope"];
export type JobAccepted = components["schemas"]["JobAccepted"];
export type DnaAnalyzeRequest = components["schemas"]["DnaAnalyzeRequest"];
export type RnaAnalyzeRequest = components["schemas"]["RnaAnalyzeRequest"];
export type TranscribeRequest = components["schemas"]["TranscribeRequest"];
export type ProteinAnalyzeRequest = components["schemas"]["ProteinAnalyzeRequest"];
export type ProteinUniprotRequest = components["schemas"]["ProteinUniprotRequest"];
export type PairwiseAlignRequest = components["schemas"]["PairwiseAlignRequest"];
export type MotifSearchRequest = components["schemas"]["MotifSearchRequest"];
export type CrisprGuidesRequest = components["schemas"]["CrisprGuidesRequest"];
export type NcbiFetchRequest = components["schemas"]["NcbiFetchRequest"];
export type NcbiSearchRequest = components["schemas"]["NcbiSearchRequest"];
export type MsaPrealignedRequest = components["schemas"]["MsaPrealignedRequest"];
export type MsaEngineJobRequest = components["schemas"]["MsaEngineJobRequest"];
export type PhylogenyInferRequest = components["schemas"]["PhylogenyInferRequest"];
export type StructureFixtureRequest = components["schemas"]["StructureFixtureRequest"];
export type StructureMappedRequest = components["schemas"]["StructureMappedRequest"];
export type StructureMappingRequest = components["schemas"]["StructureMappingRequest"];
export type CompareParseRequest = components["schemas"]["CompareParseRequest"];
export type CompareRemoteJobRequest = components["schemas"]["CompareRemoteJobRequest"];
export type CompareUSalignJobRequest = components["schemas"]["CompareUSalignJobRequest"];
export type VariantIdentifyRequest = components["schemas"]["VariantIdentifyRequest"];
export type VariantExploreRequest = components["schemas"]["VariantExploreRequest"];
export type EvidencePackRequest = components["schemas"]["EvidencePackRequest"];
export type ExplainRequest = components["schemas"]["ExplainRequest"];
export type BlastLocalJobRequest = components["schemas"]["BlastLocalJobRequest"];
export type BlastRemoteJobRequest = components["schemas"]["BlastRemoteJobRequest"];
export type CasOffinderJobRequest = components["schemas"]["CasOffinderJobRequest"];

export type JsonMap = Record<string, unknown>;

export function apiBaseUrl(): string {
  const raw = process.env.NEXT_PUBLIC_HELIX_API_BASE_URL;
  if (typeof raw === "string" && raw.trim() !== "") {
    return raw.trim().replace(/\/$/, "");
  }
  // NEXT_PUBLIC_VERCEL_ENV is injected by Vercel into the browser bundle. process.env.VERCEL is not.
  if (process.env.NEXT_PUBLIC_VERCEL_ENV) {
    return "";
  }
  if (typeof window !== "undefined") {
    const host = window.location.hostname;
    if (host.endsWith(".vercel.app") || host.endsWith(".vercel.sh")) {
      return "";
    }
  }
  return "http://127.0.0.1:8000";
}

function newRequestId(): string {
  if (typeof crypto !== "undefined" && "randomUUID" in crypto) {
    return crypto.randomUUID();
  }
  return `hs-${Date.now()}-${Math.random().toString(16).slice(2)}`;
}

async function parseBody(response: Response): Promise<unknown> {
  const text = await response.text();
  if (!text) return {};
  try {
    return JSON.parse(text) as unknown;
  } catch {
    throw new HelixApiError({
      kind: "API_DOMAIN",
      code: "PARSE_ERROR",
      message: "The API returned a non-JSON body.",
      status: response.status,
    });
  }
}

export async function helixRequest<T>(
  path: string,
  init: {
    method?: string;
    body?: unknown;
    signal?: AbortSignal;
    requestId?: string;
    headers?: Record<string, string>;
  } = {},
): Promise<{ data: T; requestId: string; status: number }> {
  const requestId = init.requestId || newRequestId();
  const headers: Record<string, string> = {
    Accept: "application/json",
    "X-Request-ID": requestId,
    ...init.headers,
  };
  let body: BodyInit | undefined;
  if (init.body instanceof FormData) {
    body = init.body;
  } else if (init.body !== undefined) {
    headers["Content-Type"] = "application/json";
    body = JSON.stringify(init.body);
  }
  let response: Response;
  try {
    response = await fetch(`${apiBaseUrl()}${path}`, {
      method: init.method ?? "GET",
      headers,
      body,
      signal: init.signal,
    });
  } catch (error) {
    if (error instanceof DOMException && error.name === "AbortError") throw error;
    throw new HelixApiError({
      kind: "NETWORK",
      code: "NETWORK",
      message: "Cannot reach HelixScope API v1. Confirm FastAPI is running.",
      requestId,
    });
  }
  const payload = await parseBody(response);
  const returnedId =
    (asRecord(payload).request_id as string | undefined) ||
    (asRecord(asRecord(payload).error).request_id as string | undefined) ||
    requestId;
  if (!response.ok) {
    const err = asRecord(asRecord(payload).error);
    const code = String(err.code || `HTTP_${response.status}`);
    throw new HelixApiError({
      kind: classifyHttpError(response.status, code),
      code,
      message: String(err.message || response.statusText || "Request failed"),
      status: response.status,
      requestId: String(err.request_id || returnedId),
      details: asRecord(err.details),
    });
  }
  return { data: payload as T, requestId: returnedId, status: response.status };
}

export const helixApi = {
  healthLive: (signal?: AbortSignal) =>
    helixRequest<JsonMap>("/health/live", { signal }),
  healthReady: (signal?: AbortSignal) =>
    helixRequest<JsonMap>("/health/ready", { signal }),
  systemIdentity: (signal?: AbortSignal) =>
    helixRequest<JsonMap>("/api/v1/system/identity", { signal }),
  capabilities: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/system/capabilities", { signal }),
  references: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/references", { signal }),
  dnaAnalyze: (body: DnaAnalyzeRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/dna/analyze", { method: "POST", body, signal }),
  dnaHelix3d: (
    body: {
      sequence: string;
      prefer: string;
      lod: string;
      region_start: number;
      region_end: number | null;
      inspect_position: number;
      center_on_selected: boolean;
    },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/dna/helix-3d", { method: "POST", body, signal }),
  rnaAnalyze: (body: RnaAnalyzeRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/rna/analyze", { method: "POST", body, signal }),
  rnaHelix3d: (
    body: {
      sequence: string;
      prefer: string;
      lod: string;
      region_start: number;
      region_end: number | null;
      inspect_position: number;
      center_on_selected: boolean;
    },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/rna/helix-3d", { method: "POST", body, signal }),
  rnaTranscribe: (body: TranscribeRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/rna/transcribe", { method: "POST", body, signal }),
  rnaTranslateCoding: (body: { sequence: string; include_explanation?: boolean }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/rna/translate-coding", { method: "POST", body, signal }),
  proteinAnalyze: (body: ProteinAnalyzeRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/protein/analyze", { method: "POST", body, signal }),
  proteinUniprot: (body: ProteinUniprotRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/protein/uniprot", { method: "POST", body, signal }),
  proteinStructureCapabilities: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/protein/structure/capabilities", { signal }),
  proteinStructureSearch: (
    body: {
      sequence: string;
      molecule?: string;
      pdb_id?: string;
      uniprot?: string;
      alphafold_id?: string;
      search_pdb_by_sequence?: boolean;
      lookup_uniprot?: boolean;
      lookup_alphafold?: boolean;
      identity_cutoff?: number;
    },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/protein/structure/search", { method: "POST", body, signal }),
  proteinStructureLoad: (
    body: { structure_id: string; sequence?: string; molecule?: string; chain?: string; bundled?: boolean },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/protein/structure/load", { method: "POST", body, signal }),
  proteinStructureAnalyze: (
    body: {
      structure_id: "1CRN" | "1BNA" | "1RNA";
      run_sasa?: boolean;
      run_dssp?: boolean;
      run_dssp_remote?: boolean;
      run_stride?: boolean;
      run_edtsurf?: boolean;
    },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/protein/structure/analyze", { method: "POST", body, signal }),
  proteinAlphafold: (body: { uniprot: string }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/protein/structure/alphafold", { method: "POST", body, signal }),
  alignmentPairwise: (body: PairwiseAlignRequest & { translate_nucleic?: boolean }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/alignment/pairwise", { method: "POST", body, signal }),
  motifSearch: (body: MotifSearchRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/motif/search", { method: "POST", body, signal }),
  explain: (body: ExplainRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/explain", { method: "POST", body, signal }),
  ncbiFetch: (body: NcbiFetchRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/ncbi/fetch", { method: "POST", body, signal }),
  ncbiSearch: (body: NcbiSearchRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/ncbi/search", { method: "POST", body, signal }),
  crisprAvailability: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/crispr/availability", { signal }),
  crisprSystems: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/crispr/systems", { signal }),
  crisprGuides: (body: CrisprGuidesRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/crispr/guides", { method: "POST", body, signal }),
  casOffinderAvailability: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/crispr/casoffinder/availability", { signal }),
  casOffinderSubmit: (body: CasOffinderJobRequest & { cas_system?: string }, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/crispr/casoffinder/jobs", { method: "POST", body, signal }),
  variantIdentify: (body: VariantIdentifyRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/variants/identify", { method: "POST", body, signal }),
  variantExplore: (body: VariantExploreRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/variants/explore", { method: "POST", body, signal }),
  variantExploreSubmit: (body: VariantExploreRequest, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/variants/jobs", { method: "POST", body, signal }),
  evidencePack: (body: EvidencePackRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/evidence/pack", { method: "POST", body, signal }),
  structureCatalog: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/catalog/experimental-complexes", { signal }),
  structureExperimentalComplex: (body: { structure_id: "4UN3" | "4OO8" }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/experimental-complex", { method: "POST", body, signal }),
  structureFixture: (body: StructureFixtureRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/fixture", { method: "POST", body, signal }),
  structureFixtureMapped: (body: StructureMappedRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/fixture/mapped", { method: "POST", body, signal }),
  structureClassifyMapping: (body: StructureMappingRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/mapping/classify", { method: "POST", body, signal }),
  structureCoordinateFile: (body: { structure_id: "1CRN" | "1BNA" | "1RNA" }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/structures/coordinate-file", { method: "POST", body, signal }),
  compareParse: (body: CompareParseRequest, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/compare/parse", { method: "POST", body, signal }),
  compareRemoteSubmit: (
    body: CompareRemoteJobRequest & { reference_chain?: string; target_chain?: string; method?: string },
    signal?: AbortSignal,
  ) => helixRequest<JobAccepted>("/api/v1/compare/jobs/remote", { method: "POST", body, signal }),
  compareUSalignSubmit: (body: CompareUSalignJobRequest, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/compare/jobs/usalign", { method: "POST", body, signal }),
  compareModes: (signal?: AbortSignal) => helixRequest<ApiEnvelope>("/api/v1/compare/modes", { signal }),
  compareMethods: (signal?: AbortSignal) => helixRequest<ApiEnvelope>("/api/v1/compare/methods", { signal }),
  compareVariants: (
    body: { text_a: string; assembly_a: string; text_b: string; assembly_b: string },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/compare/variants", { method: "POST", body, signal }),
  compareProteins: (
    body: { sequence_a: string; sequence_b: string; identifier_a?: string; identifier_b?: string },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/compare/proteins", { method: "POST", body, signal }),
  compareGuides: (
    body: { guide_a: Record<string, unknown>; guide_b: Record<string, unknown> },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/compare/guides", { method: "POST", body, signal }),
  compareEvolution: (
    body: { fasta: string; group_ids?: string[]; column?: number; member_id?: string },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/compare/evolution", { method: "POST", body, signal }),
  compareEvidence: (
    body: { text_a: string; assembly_a: string; text_b: string; assembly_b: string },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/compare/evidence", { method: "POST", body, signal }),
  msaAvailability: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/msa/availability", { signal }),
  msaPrealigned: (body: MsaPrealignedRequest & { format?: string }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/msa/prealigned", { method: "POST", body, signal }),
  msaEngineSubmit: (body: MsaEngineJobRequest, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/msa/jobs", { method: "POST", body, signal }),
  phylogenyInfer: (body: PhylogenyInferRequest & { alignment_hash?: string | null }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/phylogeny/infer", { method: "POST", body, signal }),
  phylogenySubmit: (body: PhylogenyInferRequest & { alignment_hash?: string | null }, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/phylogeny/jobs", { method: "POST", body, signal }),
  phylogenyTaxonomy: (
    body: { tree: Record<string, unknown>; email: string; declared_organisms?: Record<string, string> },
    signal?: AbortSignal,
  ) => helixRequest<ApiEnvelope>("/api/v1/phylogeny/taxonomy", { method: "POST", body, signal }),
  referencesAdmission: (body: { assembly_id: string }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/references/admission", { method: "POST", body, signal }),
  referencesInstall: (body: { assembly_id: string; confirm: boolean }, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/references/install", { method: "POST", body, signal }),
  blastAvailability: (signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>("/api/v1/blast/availability", { signal }),
  blastLocalSubmit: (body: BlastLocalJobRequest, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/blast/jobs/local", { method: "POST", body, signal }),
  blastRemoteSubmit: (body: BlastRemoteJobRequest, signal?: AbortSignal) =>
    helixRequest<JobAccepted>("/api/v1/blast/jobs/remote", { method: "POST", body, signal }),
  jobGet: (jobId: string, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>(`/api/v1/jobs/${encodeURIComponent(jobId)}`, { signal }),
  jobCancel: (jobId: string, signal?: AbortSignal) =>
    helixRequest<ApiEnvelope>(`/api/v1/jobs/${encodeURIComponent(jobId)}`, {
      method: "DELETE",
      signal,
    }),
};

export async function uploadPrealignedFasta(file: File, signal?: AbortSignal) {
  const form = new FormData();
  form.append("file", file);
  return helixRequest<ApiEnvelope>("/api/v1/msa/prealigned/upload", {
    method: "POST",
    body: form,
    signal,
  });
}
