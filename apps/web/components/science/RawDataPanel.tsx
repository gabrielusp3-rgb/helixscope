"use client";

import { asRecord, formatScientificValue } from "@/lib/api/numeric";
import { downloadJson, downloadText } from "@/lib/export/download";

const BLOCKED = ["password", "api_key", "secret", "token", "authorization"];

function redact(value: unknown): unknown {
  if (Array.isArray(value)) return value.map(redact);
  if (!value || typeof value !== "object") return value;
  const rec = asRecord(value);
  const out: Record<string, unknown> = {};
  for (const [key, item] of Object.entries(rec)) {
    const low = key.toLowerCase();
    if (BLOCKED.some((part) => low.includes(part))) continue;
    if (low.includes("path") && typeof item === "string" && (item.includes("\\") || item.startsWith("/"))) {
      continue;
    }
    out[key] = redact(item);
  }
  return out;
}

export function RawDataPanel({
  result,
  extra,
  fasta,
  fastaName,
}: {
  result: unknown;
  extra?: Record<string, unknown>;
  fasta?: string;
  fastaName?: string;
}) {
  const rec = asRecord(result);
  const payload = redact({ ...rec, ...(extra || {}) }) as Record<string, unknown>;
  const hash = String(payload.sequence_hash || payload.input_hash || rec.sequence_hash || "");
  const sequence = String(payload.sequence || rec.sequence || "");
  return (
    <section className="hs-panel science-surface" data-testid="raw-data-panel">
      <h2>Advanced / raw data</h2>
      <p className="hs-topbar-meta">
        Exact Core/API fields. Secrets, API keys and private filesystem paths are omitted.
      </p>
      <dl className="hs-raw-meta">
        <div>
          <dt>Input hash</dt>
          <dd className="mono">{hash || "Unavailable"}</dd>
        </div>
        <div>
          <dt>Normalized length</dt>
          <dd className="mono">{sequence ? String(sequence.length) : formatScientificValue(payload.length)}</dd>
        </div>
        <div>
          <dt>Status</dt>
          <dd className="mono">{formatScientificValue(payload.status)}</dd>
        </div>
      </dl>
      {sequence ? (
        <details className="disclosure">
          <summary>Normalized sequence</summary>
          <pre className="hs-seq-viewer mono">{sequence.replace(/(.{80})/g, "$1\n").trim()}</pre>
        </details>
      ) : null}
      <details className="disclosure">
        <summary>Structured JSON</summary>
        <pre className="mono hs-raw-json">{JSON.stringify(payload, null, 2)}</pre>
      </details>
      <div className="hs-actions">
        <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-raw.json", payload)}>
          Download JSON
        </button>
        {fasta ? (
          <button
            type="button"
            className="btn-secondary"
            onClick={() => downloadText(fastaName || "helixscope.fasta", fasta)}
          >
            Download FASTA
          </button>
        ) : null}
      </div>
    </section>
  );
}
