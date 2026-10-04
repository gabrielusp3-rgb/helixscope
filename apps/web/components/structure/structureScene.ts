export type StructureScene = {
  coordinates?: unknown;
  kind?: unknown;
  source?: unknown;
  structure_id?: unknown;
  mapping_status?: unknown;
  n_atoms?: unknown;
  structure_hash?: unknown;
  status?: unknown;
};

export const BUNDLED_MOLSTAR_IDS = new Set(["1CRN", "1BNA", "1RNA"]);

export function evidenceBanner(scene: StructureScene | Record<string, unknown> | null): string {
  const rec = scene && typeof scene === "object" ? (scene as Record<string, unknown>) : {};
  const kind = String(rec.kind || rec.status || "").toLowerCase();
  const id = String(rec.structure_id || "");
  const mapping = String(rec.mapping_status || "");
  if (kind.includes("illustrative") || String(rec.status).includes("ILLUSTRATIVE")) {
    return `ILLUSTRATIVE · sequence-derived ${id || "helix"}`;
  }
  if (kind.includes("predicted") || String(rec.status).includes("PREDICTED")) {
    return `PREDICTED · ${id || "AlphaFold"}`;
  }
  if (id === "4UN3" || mapping === "UNCERTAIN") {
    return `EXPERIMENTAL STRUCTURE / UNCERTAIN GUIDE MAPPING · ${id || "complex"}`;
  }
  if (id) return `EXPERIMENTAL · PDB ${id}`;
  return "COORDINATES FROM API";
}
