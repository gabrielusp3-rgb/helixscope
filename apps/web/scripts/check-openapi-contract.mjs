/**
 * Fail when services/api/openapi.json changed but TypeScript types were not regenerated.
 */
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(here, "../../..");
const openapiPath = resolve(repoRoot, "services/api/openapi.json");
const hashFile = resolve(here, "../lib/api/generated/openapi.sha256");
const contractFile = resolve(here, "../lib/api/generated/contract.ts");
const schemaFile = resolve(here, "../lib/api/generated/schema.d.ts");

const current = createHash("sha256").update(readFileSync(openapiPath)).digest("hex").trim();
let recorded = "";
try {
  recorded = readFileSync(hashFile, "utf8").trim();
} catch {
  console.error("Missing lib/api/generated/openapi.sha256. Run npm run api:types.");
  process.exit(1);
}
try {
  readFileSync(schemaFile);
} catch {
  console.error("Missing lib/api/generated/schema.d.ts. Run npm run api:types.");
  process.exit(1);
}
if (current !== recorded) {
  console.error("OpenAPI contract drift detected.");
  console.error(`  committed openapi.json sha256: ${current}`);
  console.error(`  generated types sha256:       ${recorded}`);
  console.error("Run: npm run api:types");
  process.exit(1);
}
let contract = "";
try {
  contract = readFileSync(contractFile, "utf8");
} catch {
  console.error("Missing lib/api/generated/contract.ts. Run npm run api:types.");
  process.exit(1);
}
if (!contract.includes(current)) {
  console.error("Generated contract.ts does not contain the committed OpenAPI sha256.");
  process.exit(1);
}
console.log("OpenAPI TypeScript contract is current.");
