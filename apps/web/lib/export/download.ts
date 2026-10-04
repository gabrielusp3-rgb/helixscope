import { formatScientificValue } from "@/lib/api/numeric";

/** Client-side file download from already-fetched API payloads. Not a scientific source. */

export function downloadText(filename: string, text: string, mime = "text/plain;charset=utf-8"): void {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.rel = "noopener";
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
}

export function downloadJson(filename: string, value: unknown): void {
  downloadText(filename, `${JSON.stringify(value, null, 2)}\n`, "application/json;charset=utf-8");
}

export function fastaRecord(identifier: string, sequence: string): string {
  const id = identifier.replace(/\s+/g, "_") || "sequence";
  const wrapped = sequence.replace(/(.{80})/g, "$1\n").replace(/\n$/u, "");
  return `>${id}\n${wrapped}\n`;
}

export function downloadCsv(filename: string, rows: Array<Record<string, unknown>>): void {
  if (rows.length === 0) {
    downloadText(filename, "", "text/csv;charset=utf-8");
    return;
  }
  const columns = Array.from(new Set(rows.flatMap((row) => Object.keys(row))));
  const escape = (value: unknown): string =>
    `"${formatScientificValue(value, "").replace(/"/g, '""')}"`;
  const lines = [
    columns.join(","),
    ...rows.map((row) => columns.map((column) => escape(row[column])).join(",")),
  ];
  downloadText(filename, `${lines.join("\n")}\n`, "text/csv;charset=utf-8");
}
