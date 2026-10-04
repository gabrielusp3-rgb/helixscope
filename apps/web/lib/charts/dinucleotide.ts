import { asRecord, chartNumeric } from "@/lib/api/numeric";

/** Reshape Core dinucleotide dict into a 4x4 O/E matrix. Does not recompute O/E. */
export function dinucleotideOeMatrix(dinucleotides: unknown): {
  x: string[];
  y: string[];
  z: Array<Array<number | null>>;
} {
  const rec = asRecord(dinucleotides);
  const bases = ["A", "C", "G", "T"];
  const z = bases.map((first) =>
    bases.map((second) => {
      const cell = asRecord(rec[`${first}${second}`]);
      return chartNumeric(cell.observed_expected);
    }),
  );
  return { x: bases, y: bases, z };
}
