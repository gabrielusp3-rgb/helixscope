import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, resolve } from "node:path";
import { describe, expect, it } from "vitest";

const ROOT = resolve(__dirname, "../..");
const SKIP = new Set(["node_modules", ".next", "playwright-report", "test-results", "public", "scripts"]);

function files(dir: string, acc: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (SKIP.has(name)) continue;
    const full = join(dir, name);
    const st = statSync(full);
    if (st.isDirectory()) files(full, acc);
    else if (/\.(ts|tsx|js|mjs)$/.test(name) && !full.includes("generated") && !full.includes("tests")) acc.push(full);
  }
  return acc;
}

describe("no client-side biology", () => {
  it("only lib/api/client.ts calls fetch", () => {
    const hits: string[] = [];
    for (const file of files(ROOT)) {
      if (file.replace(/\\/g, "/").endsWith("lib/api/client.ts")) continue;
      const text = readFileSync(file, "utf8");
      if (/\bfetch\s*\(/.test(text)) hits.push(file);
    }
    expect(hits).toEqual([]);
  });

  it("does not compute GC, Tm, RMSD, identity, or codon counts in the workstation", () => {
    const forbidden = [
      /function\s+computeGc/i,
      /gc_content\s*=\s*\(/,
      /identity_pct\s*=\s*\(/,
      /rmsd\s*=\s*Math/,
      /length\s*\/\s*3\s*\*\s*100/,
    ];
    const hits: string[] = [];
    for (const file of files(ROOT)) {
      const text = readFileSync(file, "utf8");
      if (forbidden.some((re) => re.test(text))) hits.push(file);
    }
    expect(hits).toEqual([]);
  });
});
