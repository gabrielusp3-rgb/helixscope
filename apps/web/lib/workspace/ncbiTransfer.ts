import { asRecord } from "@/lib/api/numeric";

export type NcbiTarget = "dna" | "rna" | "protein" | "alignment" | "blast" | "msa";

export function ncbiRecord(result: Record<string, unknown>): Record<string, unknown> {
  return asRecord(result.record ?? result);
}

export function ncbiRetrievedSequence(result: Record<string, unknown>): string {
  const record = ncbiRecord(result);
  return String(record.sequence ?? result.sequence ?? "").replace(/\s+/g, "");
}

export function ncbiAccession(result: Record<string, unknown>): string {
  const record = ncbiRecord(result);
  return String(record.accession ?? record.Accession ?? record.accession_version ?? result.accession ?? "");
}

export function ncbiIsPubmed(result: Record<string, unknown>, db: string): boolean {
  const record = ncbiRecord(result);
  return db === "pubmed" || String(record.record_kind || "") === "pubmed";
}

export function suggestedNcbiTarget(result: Record<string, unknown>, db: string): NcbiTarget | null {
  if (ncbiIsPubmed(result, db)) return null;
  const sequence = ncbiRetrievedSequence(result);
  if (!sequence) return null;
  if (db === "protein") return "protein";
  const record = ncbiRecord(result);
  const mol = String(record.molecule || record.mol_type || "").toUpperCase();
  if (mol.includes("PROTEIN") || mol.includes("AA")) return "protein";
  if (/T/.test(sequence)) return "dna";
  if (/U/.test(sequence)) return "rna";
  if (mol.includes("RNA") && !mol.includes("DNA")) return "rna";
  return "dna";
}
