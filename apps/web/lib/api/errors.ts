export type ErrorKind =
  | "NETWORK"
  | "API_DOMAIN"
  | "SCIENTIFIC_UNAVAILABLE"
  | "JOB_FAILURE"
  | "TIMEOUT"
  | "ABORTED";

export class HelixApiError extends Error {
  readonly kind: ErrorKind;
  readonly code: string;
  readonly status: number;
  readonly requestId: string;
  readonly details: Record<string, unknown>;

  constructor(init: {
    kind: ErrorKind;
    code: string;
    message: string;
    status?: number;
    requestId?: string;
    details?: Record<string, unknown>;
  }) {
    super(init.message);
    this.name = "HelixApiError";
    this.kind = init.kind;
    this.code = init.code;
    this.status = init.status ?? 0;
    this.requestId = init.requestId ?? "";
    this.details = init.details ?? {};
  }
}

const DOMAIN_CODES = new Set([
  "INVALID_INPUT",
  "NO_RECORD",
  "JOB_NOT_FOUND",
  "PARSE_ERROR",
  "RESOURCE_LIMIT",
  "CANCEL_NOT_SUPPORTED",
  "ENGINE_NOT_INSTALLED",
  "TIMEOUT",
  "RATE_LIMITED",
  "REMOTE_UNAVAILABLE",
  "ENGINE_FAILED",
  "INTERNAL",
  "INVALID_FORMAT",
  "NO_HITS",
  "METHOD_NOT_APPLICABLE",
  "MAPPING_UNCERTAIN",
]);

export function classifyHttpError(status: number, code: string): ErrorKind {
  const normalized = code.toUpperCase();
  if (normalized === "TIMEOUT" || status === 504) return "TIMEOUT";
  if (normalized === "ENGINE_NOT_INSTALLED" || normalized === "REMOTE_UNAVAILABLE") {
    return "SCIENTIFIC_UNAVAILABLE";
  }
  if (DOMAIN_CODES.has(normalized) || status >= 400) return "API_DOMAIN";
  return "API_DOMAIN";
}

export function userFacingError(error: unknown): {
  title: string;
  message: string;
  code: string;
  kind: ErrorKind;
  details?: Record<string, unknown>;
} {
  if (error instanceof HelixApiError) {
    return {
      title: titleFor(error),
      message: error.message,
      code: error.code,
      kind: error.kind,
      details: error.details,
    };
  }
  if (error instanceof DOMException && error.name === "AbortError") {
    return {
      title: "Request cancelled",
      message: "The in-flight request was aborted because a newer request started.",
      code: "ABORTED",
      kind: "ABORTED",
    };
  }
  return {
    title: "Backend unavailable",
    message: error instanceof Error ? error.message : "The HelixScope API could not be reached.",
    code: "NETWORK",
    kind: "NETWORK",
  };
}

function titleFor(error: HelixApiError): string {
  switch (error.code) {
    case "RESOURCE_LIMIT":
      return "Resource limit";
    case "ENGINE_NOT_INSTALLED":
      return "Engine not installed";
    case "RATE_LIMITED":
      return "Rate limited";
    case "REMOTE_UNAVAILABLE":
      return "Remote service unavailable";
    case "INVALID_INPUT":
      return "Invalid input";
    case "CANCEL_NOT_SUPPORTED":
      return "Cancellation not supported";
    case "JOB_NOT_FOUND":
      return "Job not found";
    case "NO_RECORD":
      return "No record";
    case "INVALID_FORMAT":
      return "Unrecognized format";
    case "NO_HITS":
      return "No hits";
    case "METHOD_NOT_APPLICABLE":
      return "Method not applicable";
    case "MAPPING_UNCERTAIN":
      return "Mapping uncertain";
    case "TIMEOUT":
      return "Timed out";
    default:
      if (error.kind === "NETWORK") return "Backend unavailable";
      if (error.kind === "SCIENTIFIC_UNAVAILABLE") return "Scientific result unavailable";
      if (error.kind === "JOB_FAILURE") return "Job failed";
      return "API error";
  }
}
