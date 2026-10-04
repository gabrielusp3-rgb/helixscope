import { helixApi, type ApiEnvelope, type JobAccepted } from "@/lib/api/client";
import { HelixApiError } from "@/lib/api/errors";
import { asRecord } from "@/lib/api/numeric";

export const TERMINAL_JOB_STATES = new Set([
  "COMPLETED",
  "FAILED",
  "CANCELLED",
  "TIMEOUT",
  "RESOURCE_LIMIT",
]);

export function sleep(ms: number, signal?: AbortSignal): Promise<void> {
  return new Promise((resolve, reject) => {
    if (signal?.aborted) {
      reject(new DOMException("Aborted", "AbortError"));
      return;
    }
    const timer = setTimeout(resolve, ms);
    signal?.addEventListener(
      "abort",
      () => {
        clearTimeout(timer);
        reject(new DOMException("Aborted", "AbortError"));
      },
      { once: true },
    );
  });
}

export async function pollJob(
  accepted: JobAccepted,
  options: {
    signal?: AbortSignal;
    expectedHash?: string | null;
    intervalMs?: number;
    maxMs?: number;
  } = {},
): Promise<ApiEnvelope> {
  const jobId = accepted.job_id;
  const started = Date.now();
  const maxMs = options.maxMs ?? 180_000;
  let interval = options.intervalMs ?? 750;
  let last: ApiEnvelope | null = null;
  while (Date.now() - started < maxMs) {
    const { data } = await helixApi.jobGet(jobId, options.signal);
    last = data;
    const snap = asRecord(data.result);
    const execution = String(data.execution_status || snap.execution_status || "");
    if (options.expectedHash && data.input_hash && options.expectedHash !== data.input_hash) {
      throw new HelixApiError({
        kind: "JOB_FAILURE",
        code: "STALE_RESULT",
        message: "Job result input_hash does not match the submitted request.",
        details: { expected: options.expectedHash, got: data.input_hash },
      });
    }
    if (TERMINAL_JOB_STATES.has(execution)) {
      if (execution === "FAILED" || execution === "TIMEOUT" || execution === "RESOURCE_LIMIT") {
        const err = asRecord(snap.error);
        throw new HelixApiError({
          kind: execution === "RESOURCE_LIMIT" ? "SCIENTIFIC_UNAVAILABLE" : "JOB_FAILURE",
          code: String(err.code || execution),
          message: String(err.message || `Job ended with ${execution}.`),
          details: snap,
        });
      }
      return data;
    }
    await sleep(interval, options.signal);
    interval = Math.min(2000, Math.round(interval * 1.15));
  }
  throw new HelixApiError({
    kind: "TIMEOUT",
    code: "TIMEOUT",
    message: "Job polling stopped without a terminal state.",
    details: { job_id: jobId, last: last ?? undefined },
  });
}

export function jobScientificResult(envelope: ApiEnvelope): Record<string, unknown> {
  const snap = asRecord(envelope.result);
  const nested = asRecord(snap.result);
  if (Object.keys(nested).length > 0) return nested;
  return snap;
}
