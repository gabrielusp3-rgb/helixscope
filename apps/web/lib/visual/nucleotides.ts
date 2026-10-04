/** Shared nucleotide/residue presentation colors. Not biological computation. */

export const NUCLEOTIDE_COLORS: Record<string, string> = {
  A: "#3dd68c",
  C: "#5b9dff",
  G: "#f0c14a",
  T: "#ff6b8a",
  U: "#c084fc",
  N: "#6b88b0",
  "-": "#243044",
  ".": "#243044",
};

export const AMINO_COLORS: Record<string, string> = {
  D: "#f87171",
  E: "#f87171",
  K: "#60a5fa",
  R: "#60a5fa",
  H: "#818cf8",
  F: "#fbbf24",
  Y: "#fbbf24",
  W: "#fbbf24",
  A: "#94a3b8",
  G: "#94a3b8",
  I: "#94a3b8",
  L: "#94a3b8",
  M: "#94a3b8",
  V: "#94a3b8",
  S: "#34d399",
  T: "#34d399",
  N: "#34d399",
  Q: "#34d399",
  C: "#f59e0b",
  P: "#c084fc",
};

export function residueColor(symbol: string, molecule?: string): string {
  const ch = String(symbol || "").toUpperCase();
  if (molecule === "PROTEIN") return AMINO_COLORS[ch] || "#6b88b0";
  return NUCLEOTIDE_COLORS[ch] || "#6b88b0";
}

export const PLOT_PAPER = "#01040B";
export const PLOT_INNER = "#020711";
export const PLOT_GRID = "rgba(145, 180, 228, 0.10)";
export const PLOT_AXIS = "rgba(175, 203, 241, 0.28)";
export const PLOT_FONT = "#AFCBF1";
