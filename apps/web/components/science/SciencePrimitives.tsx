"use client";

import type { CSSProperties, ReactNode } from "react";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { userFacingError } from "@/lib/api/errors";
import { statusClass } from "@/lib/status";
import { TERMINAL_JOB_STATES } from "@/lib/api/jobs";

export type Metric = {
  id: string;
  label: string;
  value: unknown;
  unit?: string;
};

export function ScientificStatus({ status }: { status: unknown }) {
  const label = statusClass(status);
  const cls = label.toLowerCase().replace(/_/g, "-");
  const style = styleFor(label);
  return (
    <span className="hs-status" data-status={label} style={style} title={label}>
      <span className={`status-mark status-${cls}`} aria-hidden="true" />
      {label}
    </span>
  );
}

function styleFor(label: string): CSSProperties {
  const map: Record<string, [string, string, string]> = {
    EXPERIMENTAL: ["var(--status-experimental-fg)", "var(--status-experimental-bg)", "var(--status-experimental-border)"],
    PREDICTED: ["var(--status-predicted-fg)", "var(--status-predicted-bg)", "var(--status-predicted-border)"],
    ILLUSTRATIVE: ["var(--status-illustrative-fg)", "var(--status-illustrative-bg)", "var(--status-illustrative-border)"],
    COMPUTED: ["var(--status-computed-fg)", "var(--status-computed-bg)", "var(--status-computed-border)"],
    RETRIEVED: ["var(--status-retrieved-fg)", "var(--status-retrieved-bg)", "var(--status-retrieved-border)"],
    HEURISTIC: ["var(--status-heuristic-fg)", "var(--status-heuristic-bg)", "var(--status-heuristic-border)"],
    PARTIAL: ["var(--status-partial-fg)", "var(--status-partial-bg)", "var(--status-partial-border)"],
    UNCERTAIN: ["var(--status-uncertain-fg)", "var(--status-uncertain-bg)", "var(--status-uncertain-border)"],
    UNAVAILABLE: ["var(--status-unavailable-fg)", "var(--status-unavailable-bg)", "var(--status-unavailable-border)"],
    NOT_INSTALLED: ["var(--status-unavailable-fg)", "var(--status-unavailable-bg)", "var(--status-unavailable-border)"],
    RESOURCE_LIMIT: ["var(--status-resource-fg)", "var(--status-resource-bg)", "var(--status-resource-border)"],
    ERROR: ["var(--status-error-fg)", "var(--status-error-bg)", "var(--status-error-border)"],
    LIVE_VALIDATED: ["var(--status-live-fg)", "var(--status-live-bg)", "var(--status-live-border)"],
    REMOTE_VALIDATED: ["var(--status-remote-fg)", "var(--status-remote-bg)", "var(--status-remote-border)"],
    TEST_ONLY: ["var(--status-test-fg)", "var(--status-test-bg)", "var(--status-test-border)"],
  };
  const [color, background, borderColor] = map[label] ?? map.UNAVAILABLE;
  return { color, background, borderColor };
}

export function ScientificMetric({ label, value, unit, testId }: Metric & { testId?: string }) {
  return (
    <div className="metric-card" data-testid={testId}>
      <div className="label">{label}</div>
      <div className="value-row">
        <span className="value" data-metric-value>
          {formatScientificValue(value)}
        </span>
        {unit ? <span className="unit">{unit}</span> : null}
      </div>
    </div>
  );
}

export function MetricGrid({ metrics }: { metrics: Metric[] }) {
  return (
    <div className="metric-grid" data-testid="metric-grid">
      {metrics.map((metric) => (
        <ScientificMetric key={metric.id} {...metric} testId={`metric-${metric.id}`} />
      ))}
    </div>
  );
}

export function EmptyState({ title, children }: { title: string; children?: ReactNode }) {
  return (
    <div className="empty-state science-surface hs-panel" data-testid="empty-state">
      <h3>{title}</h3>
      <div>{children}</div>
    </div>
  );
}

export function ScientificError({ error }: { error: unknown }) {
  const info = userFacingError(error);
  return (
    <div className="error-box" role="alert" data-testid="scientific-error" data-error-code={info.code}>
      <strong>{info.title}</strong>
      <div>{info.message}</div>
      <div className="hs-topbar-meta">
        {info.kind} / {info.code}
      </div>
      {info.details && Object.keys(info.details).length > 0 ? (
        <details className="disclosure">
          <summary>Technical details</summary>
          <pre className="mono">{JSON.stringify(info.details, null, 2)}</pre>
        </details>
      ) : null}
    </div>
  );
}

export function ScientificWarning({ children }: { children: ReactNode }) {
  return (
    <div className="warning" role="status">
      {children}
    </div>
  );
}

export function JobState({
  execution,
  jobId,
  onCancel,
  cancelSupported,
}: {
  execution: string;
  jobId?: string;
  onCancel?: () => void;
  cancelSupported?: boolean;
}) {
  const terminal = TERMINAL_JOB_STATES.has(execution);
  return (
    <div className="hs-panel science-surface" data-testid="job-state" data-execution={execution}>
      <div className="chip-row">
        <ScientificStatus status={execution === "COMPLETED" ? "COMPUTED" : execution} />
        <span className="hs-topbar-meta">Execution {execution || "UNKNOWN"}</span>
        {jobId ? <span className="mono" style={{ padding: "2px 6px" }}>{jobId}</span> : null}
      </div>
      <p className="hs-lede" style={{ marginBottom: 0 }}>
        Job progress percentages are not displayed. The API does not invent completion percent.
      </p>
      {onCancel && cancelSupported && !terminal ? (
        <button type="button" className="btn-secondary" onClick={onCancel}>
          Cancel queued job
        </button>
      ) : null}
    </div>
  );
}

export function SequenceInput({
  id,
  label,
  value,
  onChange,
  rows = 8,
  hint,
  disabled,
}: {
  id: string;
  label: string;
  value: string;
  onChange: (value: string) => void;
  rows?: number;
  hint?: string;
  disabled?: boolean;
}) {
  return (
    <div className="field">
      <label htmlFor={id}>{label}</label>
      {hint ? <p className="hs-topbar-meta">{hint}</p> : null}
      <textarea
        id={id}
        rows={rows}
        value={value}
        onChange={(event) => onChange(event.target.value)}
        spellCheck={false}
        autoComplete="off"
        disabled={disabled}
      />
    </div>
  );
}

export function InvalidMoleculeResult({
  expected,
  result,
}: {
  expected: string;
  result: Record<string, unknown>;
}) {
  const validation = asRecord(result.validation);
  const invalid = asList(validation.invalid_chars).map((item) => String(item));
  const reason = String(
    validation.rejection_reason || result.rejection_reason || "The sequence is not valid for this module.",
  );
  const detected = String(validation.type || validation.detected || "");
  return (
    <div className="error-box" role="alert" data-testid="invalid-input">
      <strong>Invalid input</strong>
      <p>Expected {expected}.</p>
      {detected ? <p>Detected type: {detected}.</p> : null}
      <p>{reason}</p>
      {invalid.length > 0 ? <p>Unsupported characters: {invalid.join(", ")}.</p> : null}
    </div>
  );
}

export function ResultSection({
  title,
  children,
}: {
  title: string;
  children: ReactNode;
}) {
  return (
    <section className="hs-panel science-surface">
      <h2>{title}</h2>
      {children}
    </section>
  );
}

export function ExternalDisclosure({ service, children }: { service: string; children?: ReactNode }) {
  return (
    <ScientificWarning>
      <div data-testid="external-disclosure">
        External transmission: this operation may send relevant sequence or identifier data to {service}.
        HelixScope does not keep accounts in this workstation. {children}
      </div>
    </ScientificWarning>
  );
}
