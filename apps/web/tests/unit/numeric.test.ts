import { describe, expect, it } from "vitest";
import { chartNumeric, formatScientificValue, isTransportSentinel } from "@/lib/api/numeric";

describe("scientific display", () => {
  it("never renders transport NaN as 0", () => {
    expect(formatScientificValue({ value: null, value_state: "NAN" })).toBe("Undefined");
    expect(chartNumeric({ value: null, value_state: "NAN" })).toBeNull();
    expect(formatScientificValue(43093556)).toBe("43093556");
    expect(formatScientificValue(50000)).toBe("50000");
    expect(formatScientificValue(0.5)).toBe("0.5");
    expect(formatScientificValue(null)).toBe("N/A");
    expect(isTransportSentinel({ value: null, value_state: "NAN" })).toBe(true);
  });
});
