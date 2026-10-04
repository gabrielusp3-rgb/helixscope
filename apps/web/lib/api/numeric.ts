/** Transport numeric display. Never converts NaN/unavailable into 0. */

export type ValueState = "NAN" | "INF" | "-INF";

export type TransportNumber =
  | number
  | string
  | boolean
  | null
  | undefined
  | { value: number | null; value_state?: ValueState | string };

export function isTransportSentinel(
  value: unknown,
): value is { value: number | null; value_state?: string } {
  return Boolean(
    value &&
      typeof value === "object" &&
      !Array.isArray(value) &&
      "value_state" in value,
  );
}

export function formatScientificValue(value: unknown, fallback = "N/A"): string {
  if (value === null || value === undefined) {
    return fallback;
  }
  if (isTransportSentinel(value)) {
    const state = String(value.value_state || "").toUpperCase();
    if (state === "NAN") return "Undefined";
    if (state === "INF") return "+Infinity";
    if (state === "-INF") return "-Infinity";
    if (value.value === null || value.value === undefined) return fallback;
    return formatFiniteNumber(value.value);
  }
  if (typeof value === "number") {
    if (Number.isNaN(value)) return "Undefined";
    if (!Number.isFinite(value)) return value > 0 ? "+Infinity" : "-Infinity";
    return formatFiniteNumber(value);
  }
  if (typeof value === "boolean") return value ? "true" : "false";
  if (typeof value === "string") return value;
  return fallback;
}

export function formatFiniteNumber(value: number): string {
  if (Object.is(value, -0) || value === 0) return "0";
  if (Number.isInteger(value)) return String(value);
  const abs = Math.abs(value);
  if (abs >= 1000 || abs < 0.001) return value.toPrecision(6);
  return String(Number(value.toPrecision(8)));
}

/** Chart Y values: NaN/unavailable become null gaps, never 0. */
export function chartNumeric(value: unknown): number | null {
  if (typeof value === "number") {
    return Number.isFinite(value) ? value : null;
  }
  if (isTransportSentinel(value)) {
    if (value.value_state) return null;
    return typeof value.value === "number" && Number.isFinite(value.value)
      ? value.value
      : null;
  }
  return null;
}

export function asRecord(value: unknown): Record<string, unknown> {
  if (value && typeof value === "object" && !Array.isArray(value)) {
    return value as Record<string, unknown>;
  }
  return {};
}

export function asList(value: unknown): unknown[] {
  return Array.isArray(value) ? value : [];
}

export function asRecordRows(value: unknown, keyName = "id"): Array<Record<string, unknown>> {
  if (Array.isArray(value)) {
    return value.map((row) => asRecord(row));
  }
  const rec = asRecord(value);
  if (Array.isArray(rec.records)) {
    return rec.records.map((row) => asRecord(row));
  }
  return Object.entries(rec).map(([key, spec]) => {
    if (spec && typeof spec === "object" && !Array.isArray(spec)) {
      return { [keyName]: key, ...asRecord(spec) };
    }
    return { [keyName]: key, value: spec };
  });
}
