import { describe, expect, it } from "vitest";
import { compareDisplayResult, overlayCoordinates } from "@/lib/compare/overlay";

describe("compare overlay coordinates", () => {
  it("does not invent XYZ from residue pairs", () => {
    const overlay = overlayCoordinates({
      residue_pairs: [{ ref: 1, tgt: 2 }],
      n_aligned_residue_pairs: 1,
    });
    expect(overlay).toEqual([]);
  });

  it("uses Core superposition atom copies when present", () => {
    const overlay = overlayCoordinates({
      superposition: {
        reference_atoms: [{ x: 1, y: 2, z: 3, name: "CA" }],
        transformed_target_atoms: [{ x: 4, y: 5, z: 6, name: "CA" }],
      },
    });
    expect(overlay).toEqual([
      { x: 1, y: 2, z: 3, name: "CA" },
      { x: 4, y: 5, z: 6, name: "CA" },
    ]);
  });

  it("flattens nested job alignment without dropping superposition", () => {
    const display = compareDisplayResult({
      mode: "structures",
      alignment: { rmsd_global_angstrom: 1.2, n_aligned_residue_pairs: 9 },
      superposition: { transformed_target_atoms: [{ x: 0, y: 0, z: 0 }] },
    });
    expect(display.rmsd_global_angstrom).toBe(1.2);
    expect(overlayCoordinates(display)).toHaveLength(1);
  });
});
