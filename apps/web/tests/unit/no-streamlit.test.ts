import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, resolve } from "node:path";
import { describe, expect, it } from "vitest";

const ROOT = resolve(__dirname, "../..");
const SKIP = new Set(["node_modules", ".next", "playwright-report", "test-results", "public"]);

function files(dir: string, acc: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (SKIP.has(name)) continue;
    const full = join(dir, name);
    const st = statSync(full);
    if (st.isDirectory()) files(full, acc);
    else if (/\.(ts|tsx|js|mjs|css)$/.test(name) && !full.includes("generated")) acc.push(full);
  }
  return acc;
}

describe("Prompt 4 prohibitions", () => {
  it("has no Streamlit dependency or st-key DOM lookup in the Next.js app", () => {
    const hits: string[] = [];
    for (const file of files(ROOT)) {
      if (file.endsWith("no-streamlit.test.ts")) continue;
      const text = readFileSync(file, "utf8");
      if (/\bstreamlit\b|\.st-key-|st\.session_state/i.test(text)) hits.push(file);
    }
    expect(hits).toEqual([]);
  });
});
