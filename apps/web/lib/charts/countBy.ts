/** Presentation grouping of already-computed API rows. Does not invent hits. */
export function countByLabel(
  rows: Array<Record<string, unknown>>,
  key: string,
): { labels: string[]; values: number[] } {
  const counts = new Map<string, number>();
  for (const row of rows) {
    const name = String(row[key] ?? "").trim() || "unnamed";
    counts.set(name, (counts.get(name) ?? 0) + 1);
  }
  const labels = Array.from(counts.keys()).sort();
  return { labels, values: labels.map((label) => counts.get(label) ?? 0) };
}
