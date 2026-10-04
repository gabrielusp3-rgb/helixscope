import { describe, expect, it } from "vitest";
import { evaluateApiContract } from "@/lib/api/contractEvaluate";

const HASH_A = "a".repeat(64);
const HASH_B = "b".repeat(64);

describe("API contract identity", () => {
  it("unlocks science only when hashes match", () => {
    expect(evaluateApiContract(HASH_A, HASH_A)).toEqual({ status: "ok", scienceLocked: false });
  });

  it("locks science on OpenAPI hash mismatch", () => {
    expect(evaluateApiContract(HASH_A, HASH_B)).toEqual({ status: "mismatch", scienceLocked: true });
  });

  it("locks science when the running hash is not a sha256 hex digest", () => {
    expect(evaluateApiContract(HASH_A, "not-a-hash")).toEqual({
      status: "mismatch",
      scienceLocked: true,
    });
  });
});
