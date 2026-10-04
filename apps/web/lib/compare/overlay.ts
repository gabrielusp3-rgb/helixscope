import { asList, asRecord } from "@/lib/api/numeric";

function xyzAtoms(value: unknown): Array<Record<string, unknown>> {
  return asList(value)
    .map((item) => asRecord(item))
    .filter((atom) => [atom.x, atom.y, atom.z].every((coord) => typeof coord === "number" && Number.isFinite(coord)));
}

export function compareDisplayResult(result: Record<string, unknown>): Record<string, unknown> {
  const alignment = asRecord(result.alignment);
  if (Object.keys(alignment).length === 0) return result;
  return {
    ...alignment,
    superposition: result.superposition ?? alignment.superposition,
    coordinate_fetch: result.coordinate_fetch,
    mode: result.mode,
    rmsd_agreement: result.rmsd_agreement,
    evidence: result.evidence,
    export: result.export,
    network: result.network,
    engine_location: result.engine_location,
  };
}

export function overlayCoordinates(result: Record<string, unknown>): Array<Record<string, unknown>> {
  const superposition = asRecord(result.superposition);
  const candidates = [
    result.overlay_atoms,
    result.transformed_coordinates,
    result.cartesian_overlay,
    superposition.coordinates,
    superposition.atoms,
    superposition.transformed_target_atoms,
    superposition.reference_atoms,
  ];
  const overlay: Array<Record<string, unknown>> = [];
  const reference = xyzAtoms(superposition.reference_atoms);
  const transformed = xyzAtoms(superposition.transformed_target_atoms);
  if (reference.length > 0 || transformed.length > 0) {
    return [...reference, ...transformed];
  }
  for (const candidate of candidates) {
    const xyz = xyzAtoms(candidate);
    if (xyz.length > 0) overlay.push(...xyz);
  }
  return overlay;
}
