/**
 * Generate TypeScript types from the committed FastAPI OpenAPI artifact.
 * Do not hand-edit lib/api/generated/.
 */
import { execFileSync } from "node:child_process";
import { createHash } from "node:crypto";
import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(here, "../../..");
const openapiPath = resolve(repoRoot, "services/api/openapi.json");
const outDir = resolve(here, "../lib/api/generated");
const outFile = resolve(outDir, "schema.d.ts");
const hashFile = resolve(outDir, "openapi.sha256");

mkdirSync(outDir, { recursive: true });
const bytes = readFileSync(openapiPath);
const hash = createHash("sha256").update(bytes).digest("hex");

execFileSync(
  process.execPath,
  [
    resolve(here, "../node_modules/openapi-typescript/bin/cli.js"),
    openapiPath,
    "-o",
    outFile,
  ],
  { stdio: "inherit" },
);

writeFileSync(hashFile, `${hash}\n`, "utf8");
writeFileSync(
  resolve(outDir, "contract.ts"),
  [
    "/** Generated API contract identity. Do not edit by hand. */",
    `export const EXPECTED_OPENAPI_SHA256 = ${JSON.stringify(hash)};`,
    "",
  ].join("\n"),
  "utf8",
);
writeFileSync(
  resolve(outDir, "README.md"),
  [
    "Generated API types. Do not edit by hand.",
    "",
    "Source: services/api/openapi.json",
    "Regenerate: npm run api:types",
    "Drift check: npm run api:types:check",
    "",
  ].join("\n"),
  "utf8",
);
console.log(`openapi sha256 ${hash}`);
