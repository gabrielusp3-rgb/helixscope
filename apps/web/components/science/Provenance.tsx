"use client";

import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import type { ApiEnvelope } from "@/lib/api/client";
import { ScientificStatus } from "@/components/science/SciencePrimitives";

export function ProvenancePanel({
  envelope,
  extra,
}: {
  envelope?: ApiEnvelope | null;
  extra?: unknown;
}) {
  const result = asRecord(envelope?.result);
  const provenance = asRecord(extra ?? result.provenance);
  const limitations = asList(envelope?.limitations).concat(asList(result.limitations));
  const warnings = asList(envelope?.warnings).concat(asList(result.warnings));
  return (
    <section className="hs-panel science-surface" data-testid="provenance-panel">
      <h2>Provenance</h2>
      <div className="chip-row">
        <ScientificStatus status={envelope?.scientific_status ?? provenance.status} />
        {envelope?.execution_status ? (
          <span className="hs-topbar-meta">execution {envelope.execution_status}</span>
        ) : null}
      </div>
      <dl className="hs-table-wrap" style={{ padding: 10 }}>
        {row("request_id", envelope?.request_id)}
        {row("api_version", envelope?.api_version)}
        {row("product_version", envelope?.product_version)}
        {row("input_hash", envelope?.input_hash ?? provenance.input_hash)}
        {row("source", provenance.source)}
        {row("algorithm", provenance.algorithm)}
        {row("software", provenance.software)}
        {row("computed_timestamp", provenance.computed_timestamp)}
        {row("retrieved_at_utc", provenance.retrieved_at_utc ?? provenance.retrieval_timestamp)}
      </dl>
      {warnings.length > 0 ? (
        <div>
          <h3>Warnings</h3>
          {warnings.map((item, index) => (
            <p key={index}>{String(item)}</p>
          ))}
        </div>
      ) : null}
      {limitations.length > 0 ? (
        <div>
          <h3>Limitations</h3>
          {limitations.map((item, index) => (
            <p key={index}>{String(item)}</p>
          ))}
        </div>
      ) : null}
    </section>
  );
}

function row(label: string, value: unknown) {
  if (value === undefined || value === null || value === "") return null;
  return (
    <div style={{ display: "grid", gridTemplateColumns: "180px 1fr", gap: 8, marginBottom: 6 }}>
      <dt className="hs-topbar-meta">{label}</dt>
      <dd className="mono" style={{ margin: 0, padding: "4px 8px" }}>
        {formatScientificValue(value)}
      </dd>
    </div>
  );
}

export function InterpretationPanel({ explanation }: { explanation: unknown }) {
  const rec = asRecord(explanation);
  const meaning = String(rec.plain_meaning || rec.interpretation || "");
  if (!meaning && !rec.method) {
    return null;
  }
  return (
    <section className="hs-panel science-surface" data-testid="interpretation-panel">
      <h2 data-testid="what-this-means">What this means</h2>
      <p>{meaning}</p>
      {rec.why_this_result ? <p>{String(rec.why_this_result)}</p> : null}
      <p className="hs-topbar-meta">
        Do not treat this summary as a diagnosis, organism assignment, or experimental measurement
        unless the status field says so.
      </p>
      <details className="disclosure" data-testid="method-limitations">
        <summary>Scientific method & limitations</summary>
        <p>
          <strong>Method. </strong>
          {String(rec.method || "See provenance.")}
        </p>
        <p>
          <strong>Assumptions and range. </strong>
          {String(rec.interpretation || meaning)}
        </p>
        <p>
          <strong>Limitations. </strong>
          {String(rec.limitations || "See provenance limitations.")}
        </p>
        <p>
          <strong>Status. </strong>
          {formatScientificValue(rec.status)}
        </p>
        <p>
          <strong>Source. </strong>
          {String(rec.source || "HelixScope API")}
          {rec.software_version ? ` (${String(rec.software_version)})` : ""}
        </p>
      </details>
    </section>
  );
}

export function MetricExplanations({ explanations }: { explanations: unknown }) {
  const rec = asRecord(explanations);
  const entries = Object.entries(rec);
  if (entries.length === 0) return null;
  return (
    <section className="hs-panel science-surface" data-testid="metric-explanations">
      <h2>Metric explanations</h2>
      <p className="hs-topbar-meta">
        Each card is derived from the API result and the method that produced that metric. Missing
        values stay unavailable; they are not replaced with zero.
      </p>
      {entries.map(([id, item]) => {
        const row = asRecord(item);
        return (
          <details key={id} className="disclosure">
            <summary>
              {String(row.title || id)} · {formatScientificValue(row.status)}
            </summary>
            <p>{String(row.plain_meaning || "")}</p>
            {row.why_this_result ? <p>{String(row.why_this_result)}</p> : null}
            <p className="hs-topbar-meta">
              Method: {String(row.method || "See provenance.")} Source: {String(row.source || "")}
            </p>
            <p>{String(row.limitations || "")}</p>
          </details>
        );
      })}
    </section>
  );
}

export function MethodDisclosure({ text }: { text: unknown }) {
  if (!text) return null;
  return (
    <details className="disclosure science-surface hs-panel">
      <summary>Method disclosure</summary>
      <p>{String(text)}</p>
    </details>
  );
}
