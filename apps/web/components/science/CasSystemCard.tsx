"use client";

import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { ScientificStatus } from "@/components/science/SciencePrimitives";

export function CasSystemCard({ system }: { system: Record<string, unknown> | undefined }) {
  const rec = asRecord(system);
  if (!rec.canonical_key) return null;
  const pfs = asRecord(rec.pfs);
  const unsupported = asList(rec.unsupported_metrics).map((row) => asRecord(row));
  const limitations = asList(rec.scientific_limitations);
  const mode = String(rec.mode || rec.editor || "nuclease");
  return (
    <section className="hs-panel science-surface" data-testid="cas-system-card">
      <h2>Selected Cas system</h2>
      <p className="hs-lede">
        Semantics come from GET /api/v1/crispr/systems for {String(rec.canonical_key)}. This is not a
        generic SpCas9 form.
      </p>
      <dl className="hs-figure-meta">
        <div>
          <dt>Canonical key</dt>
          <dd data-testid="cas-canonical">{String(rec.canonical_key)}</dd>
        </div>
        <div>
          <dt>Target molecule</dt>
          <dd data-testid="cas-target-molecule">{String(rec.target_molecule || "DNA")}</dd>
        </div>
        <div>
          <dt>Mode</dt>
          <dd data-testid="cas-mode">{mode}</dd>
        </div>
        <div>
          <dt>Guide length</dt>
          <dd>{formatScientificValue(rec.guide_length)} nt</dd>
        </div>
        <div>
          <dt>{String(rec.target_molecule) === "RNA" ? "PFS" : "PAM"}</dt>
          <dd data-testid="cas-pam">
            {String(rec.pam_display || rec.pam || pfs.note || "see registry")}
          </dd>
        </div>
        <div>
          <dt>Cut / edit</dt>
          <dd>{formatScientificValue(rec.cut_type ?? rec.edit_from)}</dd>
        </div>
      </dl>
      {String(rec.editor) === "base_editor" ? (
        <p data-testid="cas-base-editor">
          Base editor {formatScientificValue(rec.edit_from)} to {formatScientificValue(rec.edit_to)}. This
          is a nick/edit system, not a double-strand-break nuclease.
        </p>
      ) : null}
      {String(rec.editor) === "prime_editor" ? (
        <p data-testid="cas-prime-editor">
          Prime Editor is PARTIAL in Core. Nick/PBS/RTT notes may exist; HelixScope does not invent a
          pegRNA.
        </p>
      ) : null}
      {String(rec.target_molecule) === "RNA" ? (
        <p data-testid="cas-rna-note">
          RNA-targeting nuclease. DNA cut-site metrics are METHOD_NOT_APPLICABLE. Provide RNA (A/U/C/G).
        </p>
      ) : null}
      <p className="hs-topbar-meta">
        Heuristic efficiency <ScientificStatus status={rec.heuristic_efficiency_status || "HEURISTIC"} />.
        Published on-target models remain UNAVAILABLE unless the registry says otherwise.
      </p>
      {limitations.length > 0 ? (
        <ul>
          {limitations.map((item, index) => (
            <li key={index}>{String(item)}</li>
          ))}
        </ul>
      ) : null}
      {unsupported.length > 0 ? (
        <p className="hs-topbar-meta">
          Unsupported:{" "}
          {unsupported
            .map((row) => `${String(row.metric || "")} (${String(row.status || "")})`)
            .join("; ")}
        </p>
      ) : null}
    </section>
  );
}

export function offTargetScopeLabel(scope: unknown, notGenomeWide: unknown): string {
  const text = String(scope || "").toUpperCase();
  if (text.includes("COMPLETED_FULL") || text.includes("FULL READY") || text.includes("GENOME_WIDE_READY")) {
    return "FULL READY REFERENCE";
  }
  if (text.includes("TEST") || text.includes("HELIXSCOPE_TEST")) {
    return "TEST REFERENCE";
  }
  if (notGenomeWide === true || text.includes("PROVIDED") || text.includes("LOCAL") || text.includes("NOT_SCANNED") || !text) {
    return "SEQUENCE LOCAL";
  }
  return String(scope);
}
