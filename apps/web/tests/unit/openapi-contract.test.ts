import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

describe("OpenAPI TypeScript contract", () => {
  it("matches committed services/api/openapi.json", () => {
    const root = resolve(__dirname, "../../../../");
    const openapi = readFileSync(resolve(root, "services/api/openapi.json"));
    const recorded = readFileSync(
      resolve(__dirname, "../../lib/api/generated/openapi.sha256"),
      "utf8",
    ).trim();
    const current = createHash("sha256").update(openapi).digest("hex");
    expect(current).toBe(recorded);
    const contract = readFileSync(
      resolve(__dirname, "../../lib/api/generated/contract.ts"),
      "utf8",
    );
    expect(contract).toContain(current);
  });
});
