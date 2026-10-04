"use client";

import type { ReactNode } from "react";

export function SequenceViewer({
  sequence,
  hash,
  testId = "sequence-viewer",
}: {
  sequence: string;
  hash?: string;
  testId?: string;
}) {
  const wrapped = sequence.replace(/(.{80})/g, "$1\n").trim();
  return (
    <section className="hs-panel science-surface" data-testid={testId}>
      <h2>Sequence viewer</h2>
      <p className="hs-topbar-meta">
        Length {sequence.length}. Hash {hash ? String(hash).slice(0, 16) : "from API"}. This is the
        analyzed residue string, not a FASTA header.
      </p>
      <pre className="hs-seq-viewer mono">{wrapped || "(empty)"}</pre>
    </section>
  );
}

export function FastaFileField({
  id,
  label = "Or upload a FASTA file",
  disabled,
  onText,
}: {
  id: string;
  label?: string;
  disabled?: boolean;
  onText: (text: string) => void;
}) {
  return (
    <div className="field">
      <label htmlFor={id}>{label}</label>
      <input
        id={id}
        type="file"
        accept=".fa,.fasta,.txt,.fna"
        disabled={disabled}
        onChange={(event) => {
          const file = event.target.files?.[0];
          if (!file) return;
          void file.text().then(onText);
        }}
      />
    </div>
  );
}

export function CompactInput({
  compact,
  summary,
  children,
}: {
  compact: boolean;
  summary: string;
  children: ReactNode;
}) {
  if (!compact) {
    return <section className="hs-panel glass-medium">{children}</section>;
  }
  return (
    <details className="hs-panel glass-medium hs-input-compact">
      <summary>{summary}</summary>
      {children}
    </details>
  );
}
