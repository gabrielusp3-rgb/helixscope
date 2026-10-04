import { asRecordRows, chartNumeric } from "@/lib/api/numeric";

/** Pivot API codon_usage records onto AminoAcid x Codon. Does not recompute percentages. */
export function codonSynonymousMatrix(codonUsage: unknown): {
  x: string[];
  y: string[];
  z: Array<Array<number | null>>;
} {
  const rows = asRecordRows(codonUsage, "Codon");
  const x = Array.from(new Set(rows.map((row) => String(row.Codon ?? row.codon))));
  const y = Array.from(new Set(rows.map((row) => String(row.AminoAcid ?? row.amino_acid))));
  const lookup = new Map<string, number | null>();
  for (const row of rows) {
    const codon = String(row.Codon ?? row.codon);
    const amino = String(row.AminoAcid ?? row.amino_acid);
    lookup.set(`${amino}\0${codon}`, chartNumeric(row.Synonymous_pct));
  }
  const z = y.map((amino) => x.map((codon) => lookup.get(`${amino}\0${codon}`) ?? null));
  return { x, y, z };
}
