"use client";

import {
  createContext,
  useCallback,
  useContext,
  useMemo,
  useState,
  type ReactNode,
} from "react";
import type { ApiEnvelope } from "@/lib/api/client";

export type TransferKind =
  | "dna-to-crispr"
  | "dna-to-blast"
  | "protein-to-structure"
  | "msa-to-phylogeny"
  | "tree-leaf-to-msa"
  | "blast-to-msa"
  | "ncbi-to-dna"
  | "ncbi-to-rna"
  | "ncbi-to-protein"
  | "ncbi-to-alignment"
  | "ncbi-to-blast"
  | "ncbi-to-msa"
  | "uniprot-to-protein"
  | "rna-to-protein"
  | "crispr-to-casoffinder"
  | "structure-to-compare"
  | "variant-to-evidence";

export type TransferRecord = {
  kind: TransferKind;
  at: string;
  note: string;
  payload: Record<string, string>;
};

export type PresentationMode = "beginner" | "expert";

type Drafts = {
  dna: {
    sequence: string;
    includeOrfs: boolean;
    includeProfiles: boolean;
    includeRestriction: boolean;
    includeCpg: boolean;
    includeKmers: boolean;
    includeSantalucia: boolean;
    dropNested: boolean;
    minOrfLength: number;
    profileWindow: number;
    profileStep: number;
    kmerK: number;
    santaluciaNa: string;
    santaluciaOligo: string;
  };
  rna: { sequence: string; fold: boolean; codon: boolean; includeProfiles: boolean; includeKmers: boolean };
  protein: { sequence: string; ph: string; uniprotAccession: string; includeDomains: boolean };
  alignment: { seq1: string; seq2: string; mode: "global" | "local"; translateNucleic: boolean };
  motif: { sequence: string; pattern: string; molecule: "DNA" | "RNA" | "PROTEIN" };
  crispr: { sequence: string; sequenceName: string; cas: string; offTarget: boolean };
  ncbi: { accession: string; email: string; term: string; db: string };
  blast: { query: string; email: string; program: string; remote: boolean };
  msa: { fasta: string; email: string; backend: string; format: "auto" | "fasta" | "clustal" | "stockholm" | "phylip" };
  phylogeny: {
    fasta: string;
    method: string;
    distance: string;
    alignmentHash: string;
    taxonomyEmail: string;
    declaredOrganisms: string;
  };
  variant: {
    text: string;
    assembly: string;
    vep: boolean;
    clinvar: boolean;
    domains: boolean;
    email: string;
    evidenceJson: string;
  };
  compare: {
    json: string;
    reference: string;
    target: string;
    referenceChain: string;
    targetChain: string;
    method: string;
    mode: "structures" | "variants" | "proteins" | "guides" | "evolution" | "evidence";
    textA: string;
    textB: string;
    assemblyA: string;
    assemblyB: string;
    proteinA: string;
    proteinB: string;
    guideA: string;
    guideB: string;
    evolutionFasta: string;
  };
  viewer: { structureId: "1CRN" | "1BNA" | "1RNA" };
};

type Results = Partial<{
  dna: ApiEnvelope;
  rna: ApiEnvelope;
  protein: ApiEnvelope;
  proteinUniprot: ApiEnvelope;
  proteinStructure: ApiEnvelope;
  rnaTranslate: ApiEnvelope;
  alignment: ApiEnvelope;
  motif: ApiEnvelope;
  crispr: ApiEnvelope;
  ncbi: ApiEnvelope;
  blast: ApiEnvelope;
  msa: ApiEnvelope;
  phylogeny: ApiEnvelope;
  variant: ApiEnvelope;
  evidence: ApiEnvelope;
  compare: ApiEnvelope;
  viewer: ApiEnvelope;
  engines: ApiEnvelope;
  references: ApiEnvelope;
  referencesAdmission: ApiEnvelope;
  blastAvailability: ApiEnvelope;
  msaAvailability: ApiEnvelope;
  crisprAvailability: ApiEnvelope;
}>;

const emptyDrafts = (): Drafts => ({
  dna: {
    sequence: "",
    includeOrfs: true,
    includeProfiles: true,
    includeRestriction: true,
    includeCpg: true,
    includeKmers: true,
    includeSantalucia: true,
    dropNested: true,
    minOrfLength: 150,
    profileWindow: 100,
    profileStep: 50,
    kmerK: 3,
    santaluciaNa: "0.05",
    santaluciaOligo: "250",
  },
  rna: { sequence: "", fold: false, codon: true, includeProfiles: true, includeKmers: true },
  protein: { sequence: "", ph: "7.0", uniprotAccession: "", includeDomains: false },
  alignment: { seq1: "", seq2: "", mode: "global", translateNucleic: false },
  motif: { sequence: "", pattern: "", molecule: "DNA" },
  crispr: { sequence: "", sequenceName: "", cas: "SpCas9", offTarget: false },
  ncbi: { accession: "", email: "", term: "", db: "nucleotide" },
  blast: { query: "", email: "", program: "blastn", remote: false },
  msa: { fasta: "", email: "", backend: "ebi_clustalo", format: "auto" },
  phylogeny: {
    fasta: "",
    method: "neighbor_joining",
    distance: "p_distance",
    alignmentHash: "",
    taxonomyEmail: "",
    declaredOrganisms: "",
  },
  variant: {
    text: "",
    assembly: "GRCh38.p14",
    vep: false,
    clinvar: false,
    domains: false,
    email: "",
    evidenceJson: "",
  },
  compare: {
    json: "",
    reference: "8HSK",
    target: "8HSF",
    referenceChain: "A",
    targetChain: "A",
    method: "tm-align",
    mode: "structures",
    textA: "",
    textB: "",
    assemblyA: "GRCh38.p14",
    assemblyB: "GRCh38.p14",
    proteinA: "",
    proteinB: "",
    guideA: "",
    guideB: "",
    evolutionFasta: "",
  },
  viewer: { structureId: "1CRN" },
});

type WorkspaceValue = {
  drafts: Drafts;
  results: Results;
  transfers: TransferRecord[];
  selectedMsaId: string | null;
  selectedTreeLeaf: string | null;
  presentationMode: PresentationMode;
  setPresentationMode: (mode: PresentationMode) => void;
  patchDraft: <K extends keyof Drafts>(key: K, patch: Partial<Drafts[K]>) => void;
  setResult: <K extends keyof Results>(key: K, value: Results[K] | undefined) => void;
  recordTransfer: (transfer: Omit<TransferRecord, "at">) => void;
  setSelectedMsaId: (id: string | null) => void;
  setSelectedTreeLeaf: (id: string | null) => void;
};

const WorkspaceContext = createContext<WorkspaceValue | null>(null);

export function WorkstationProvider({ children }: { children: ReactNode }) {
  const [drafts, setDrafts] = useState<Drafts>(emptyDrafts);
  const [results, setResults] = useState<Results>({});
  const [transfers, setTransfers] = useState<TransferRecord[]>([]);
  const [selectedMsaId, setSelectedMsaId] = useState<string | null>(null);
  const [selectedTreeLeaf, setSelectedTreeLeaf] = useState<string | null>(null);
  const [presentationMode, setPresentationMode] = useState<PresentationMode>("beginner");

  const patchDraft = useCallback(<K extends keyof Drafts>(key: K, patch: Partial<Drafts[K]>) => {
    setDrafts((prev) => ({ ...prev, [key]: { ...prev[key], ...patch } }));
  }, []);

  const setResult = useCallback(<K extends keyof Results>(key: K, value: Results[K] | undefined) => {
    setResults((prev) => {
      if (value === undefined) {
        const next = { ...prev };
        delete next[key];
        return next;
      }
      return { ...prev, [key]: value };
    });
  }, []);

  const recordTransfer = useCallback((transfer: Omit<TransferRecord, "at">) => {
    setTransfers((prev) => [
      { ...transfer, at: new Date().toISOString() },
      ...prev.slice(0, 11),
    ]);
  }, []);

  const value = useMemo(
    () => ({
      drafts,
      results,
      transfers,
      selectedMsaId,
      selectedTreeLeaf,
      presentationMode,
      setPresentationMode,
      patchDraft,
      setResult,
      recordTransfer,
      setSelectedMsaId,
      setSelectedTreeLeaf,
    }),
    [
      drafts,
      results,
      transfers,
      selectedMsaId,
      selectedTreeLeaf,
      presentationMode,
      patchDraft,
      setResult,
      recordTransfer,
    ],
  );

  return <WorkspaceContext.Provider value={value}>{children}</WorkspaceContext.Provider>;
}

export function useWorkspace(): WorkspaceValue {
  const ctx = useContext(WorkspaceContext);
  if (!ctx) throw new Error("useWorkspace requires WorkstationProvider");
  return ctx;
}
