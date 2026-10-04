export const STATUS_ORDER = [
  "EXPERIMENTAL",
  "PREDICTED",
  "ILLUSTRATIVE",
  "COMPUTED",
  "RETRIEVED",
  "HEURISTIC",
  "PARTIAL",
  "UNCERTAIN",
  "UNAVAILABLE",
  "RESOURCE_LIMIT",
  "LIVE_VALIDATED",
  "REMOTE_VALIDATED",
  "TEST_ONLY",
  "NOT_INSTALLED",
  "ERROR",
] as const;

export type ScientificStatus = (typeof STATUS_ORDER)[number] | string;

const ALIASES: Record<string, string> = {
  VALIDATED: "EXPERIMENTAL",
  WARNING: "PARTIAL",
  NO_HITS: "PARTIAL",
  WAITING: "PARTIAL",
  READY: "RETRIEVED",
  TIMEOUT: "ERROR",
  RATE_LIMITED: "RESOURCE_LIMIT",
  FAILED: "ERROR",
  CACHED: "RETRIEVED",
  QUEUED: "PARTIAL",
  RUNNING: "PARTIAL",
  COMPLETED: "COMPUTED",
  live_validated: "LIVE_VALIDATED",
  remote_validated: "REMOTE_VALIDATED",
  not_installed: "NOT_INSTALLED",
  test_only: "TEST_ONLY",
};

export function normalizeStatus(raw: unknown): string {
  const text = String(raw ?? "").trim();
  if (!text) return "UNAVAILABLE";
  if (ALIASES[text]) return ALIASES[text];
  const upper = text.replace(/ /g, "_").toUpperCase();
  if (ALIASES[upper]) return ALIASES[upper];
  return upper;
}

export function statusClass(raw: unknown): string {
  const status = normalizeStatus(raw);
  const known = [
    "EXPERIMENTAL",
    "PREDICTED",
    "ILLUSTRATIVE",
    "COMPUTED",
    "RETRIEVED",
    "HEURISTIC",
    "PARTIAL",
    "UNCERTAIN",
    "UNAVAILABLE",
    "RESOURCE_LIMIT",
    "LIVE_VALIDATED",
    "REMOTE_VALIDATED",
    "TEST_ONLY",
    "NOT_INSTALLED",
    "ERROR",
  ];
  return known.includes(status) ? status : "UNAVAILABLE";
}
