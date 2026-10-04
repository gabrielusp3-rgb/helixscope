"use client";

import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { jobScientificResult, pollJob } from "@/lib/api/jobs";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, ExternalDisclosure, MetricGrid, ScientificError, ScientificWarning, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { RawDataPanel } from "@/components/science/RawDataPanel";

export function VariantWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.variant;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [jobExec, setJobExec] = useState("");
  const envelope = results.variant;
  const evidence = results.evidence;
  const result = asRecord(envelope?.result);
  const identity = identityRecord(result);
  const pack = asRecord(evidence?.result);

  async function identify() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.variantIdentify({ text: draft.text, assembly: draft.assembly }, gate.signal);
      if (gate.isCurrent()) setResult("variant", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function explore() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.variantExplore(
        {
          text: draft.text,
          assembly: draft.assembly,
          email: draft.email,
          enable_vep: draft.vep,
          enable_clinvar: draft.clinvar,
          enable_domains: draft.domains,
        },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("variant", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function exploreJob() {
    const gate = next();
    setError(null);
    try {
      const { data: accepted } = await helixApi.variantExploreSubmit(
        {
          text: draft.text,
          assembly: draft.assembly,
          email: draft.email,
          enable_vep: draft.vep,
          enable_clinvar: draft.clinvar,
          enable_domains: draft.domains,
        },
        gate.signal,
      );
      setJobExec(accepted.execution_status);
      const done = await pollJob(accepted, { signal: gate.signal });
      if (gate.isCurrent()) {
        setJobExec(String(done.execution_status || ""));
        setResult("variant", { ...done, result: jobScientificResult(done) });
      }
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function packEvidence() {
    const gate = next();
    setError(null);
    try {
      const parsed = JSON.parse(draft.evidenceJson || "[]") as unknown;
      const items = asList(parsed).map((row) => {
        const rec = asRecord(row);
        return {
          field: String(rec.field || ""),
          value: rec.value,
          source: String(rec.source || ""),
          evidence_status: String(rec.evidence_status || "RETRIEVED"),
          retrieved_at_utc: String(rec.retrieved_at_utc || ""),
          identifier: String(rec.identifier || ""),
          mapping_status: String(rec.mapping_status || ""),
          note: String(rec.note || ""),
        };
      });
      const { data } = await helixApi.evidencePack({ kind: "variants", items }, gate.signal);
      if (gate.isCurrent()) setResult("evidence", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Discovery</p>
      <h1>Variant Explorer</h1>
      <ScientificWarning>
        ClinVar reports classifications from ClinVar. HelixScope does not issue a diagnosis.
      </ScientificWarning>
      {(draft.vep || draft.clinvar || draft.domains) && (
        <ExternalDisclosure service="Ensembl VEP, NCBI ClinVar, InterPro and/or UniProt" />
      )}
      <CompactInput compact={Boolean(envelope)} summary="Edit variant input">
        <div className="field">
          <label htmlFor="vtext">Variant</label>
          <input id="vtext" value={draft.text} onChange={(e) => patchDraft("variant", { text: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="asm">Assembly</label>
          <input id="asm" value={draft.assembly} onChange={(e) => patchDraft("variant", { assembly: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="vemail">Email (remote layers only)</label>
          <input id="vemail" value={draft.email} onChange={(e) => patchDraft("variant", { email: e.target.value })} />
        </div>
        <label><input type="checkbox" checked={draft.vep} onChange={(e) => patchDraft("variant", { vep: e.target.checked })} /> VEP</label>
        <label><input type="checkbox" checked={draft.clinvar} onChange={(e) => patchDraft("variant", { clinvar: e.target.checked })} /> ClinVar</label>
        <label><input type="checkbox" checked={draft.domains} onChange={(e) => patchDraft("variant", { domains: e.target.checked })} /> Domains</label>
        <div className="hs-actions">
          <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void identify()}>Identify</button>
          <button type="button" className="btn-secondary" disabled={scienceLocked} onClick={() => void explore()}>Explore</button>
          <button type="button" className="btn-secondary" onClick={() => void exploreJob()}>Explore as job</button>
        </div>
        {jobExec ? <p data-testid="variant-job">{jobExec}</p> : null}
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No variant investigation yet" /> : (
        <>
          <AnalysisNav
            items={[
              { id: "variant-identity", label: "Identity" },
              { id: "variant-transcripts", label: "Transcripts" },
              { id: "variant-vep", label: "VEP", available: Boolean(asRecord(result.vep).status || asList(result.vep).length) },
              { id: "variant-clinvar", label: "ClinVar", available: Boolean(asRecord(result.clinvar).status || result.clinvar_significance) },
              { id: "variant-conflicts", label: "Conflicts" },
              { id: "variant-meaning", label: "Interpretation" },
              { id: "variant-raw", label: "Raw data" },
            ]}
          />
          <p className="kicker">Variant to assembly to gene to transcript to protein consequence to external evidence — only mapped API fields.</p>
          <ResultAnchor id="variant-identity">
          <MetricGrid
            metrics={[
              { id: "hash", label: "Identity hash", value: identity.identity_hash },
              { id: "pos", label: "Position 0-based", value: identity.position_0based },
              { id: "ref", label: "REF", value: identity.ref ?? identity.REF },
              { id: "alt", label: "ALT", value: identity.alt ?? identity.ALT },
              { id: "chr", label: "Region", value: identity.chromosome ?? identity.region ?? identity.contig },
              { id: "asm", label: "Assembly", value: identity.assembly ?? draft.assembly },
              { id: "effect", label: "Effect status", value: identity.effect_status },
              { id: "gene", label: "Gene", value: identity.gene_symbol ?? identity.gene ?? result.gene_symbol },
              { id: "tx", label: "Transcript", value: identity.transcript_id ?? result.transcript_id },
              { id: "prot", label: "Protein consequence", value: identity.protein_consequence ?? result.protein_consequence },
            ]}
          />
          <p data-testid="identity-hash" className="mono">{String(identity.identity_hash || "")}</p>
          <p data-testid="variant-position" className="mono">
            position_0based {String(identity.position_0based ?? "")}
          </p>
          <p data-testid="no-diagnosis">ClinVar reports {formatScientificValue(result.clinvar_significance ?? asRecord(result.clinvar).clinical_significance)}. HelixScope does not say pathogenic.</p>
          </ResultAnchor>
          <ResultAnchor id="variant-transcripts">
          <DataTable rows={asList(result.transcripts).map((row) => asRecord(row))} caption="Transcripts" />
          </ResultAnchor>
          {asRecord(result.vep).status || asList(result.vep).length ? (
            <ResultAnchor id="variant-vep">
              <section className="hs-panel science-surface">
                <h2>Ensembl VEP</h2>
                <p className="hs-lede">VEP annotations remain a retrieved layer. They are not a HelixScope diagnosis.</p>
                <MetricGrid
                  metrics={[
                    { id: "vep-status", label: "VEP status", value: asRecord(result.vep).status },
                    { id: "vep-csq", label: "Consequence", value: asRecord(result.vep).most_severe_consequence ?? asRecord(result.vep).consequence },
                  ]}
                />
                <DataTable rows={asList(asRecord(result.vep).transcript_consequences).map((row) => asRecord(row))} caption="VEP transcripts" emptyLabel="0 VEP transcript rows." />
              </section>
            </ResultAnchor>
          ) : null}
          {asRecord(result.clinvar).status || result.clinvar_significance ? (
            <ResultAnchor id="variant-clinvar">
              <section className="hs-panel science-surface">
                <h2>ClinVar</h2>
                <p className="hs-lede">ClinVar reports classifications from ClinVar. HelixScope does not issue a diagnosis.</p>
                <MetricGrid
                  metrics={[
                    { id: "clinvar-status", label: "ClinVar status", value: asRecord(result.clinvar).status },
                    { id: "clinvar-sig", label: "Clinical significance", value: result.clinvar_significance ?? asRecord(result.clinvar).clinical_significance },
                  ]}
                />
              </section>
            </ResultAnchor>
          ) : null}
          <ResultAnchor id="variant-conflicts">
          <DataTable rows={asList(result.conflicts).map((row) => asRecord(row))} caption="Conflicts" />
          </ResultAnchor>
          <ResultAnchor id="variant-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="variant-raw">
          <RawDataPanel result={result} extra={identity} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} extra={result} />
        </>
      )}
      <section className="hs-panel glass-medium">
        <h2>Evidence pack</h2>
        <p className="hs-lede">Pack retrieved items. Confidence remains null. This is not a diagnosis score.</p>
        <SequenceInput
          id="evidence-json"
          label="Evidence items JSON"
          value={draft.evidenceJson}
          onChange={(evidenceJson) => patchDraft("variant", { evidenceJson })}
          rows={8}
        />
          <button type="button" className="btn-primary" onClick={() => void packEvidence()}>Pack evidence</button>
          <button
            type="button"
            className="btn-secondary"
            data-testid="variant-to-compare"
            onClick={() => {
              patchDraft("compare", {
                mode: "evidence",
                textA: draft.text,
                assemblyA: draft.assembly,
              });
              recordTransfer({
                kind: "variant-to-evidence",
                note: "Variant identity transferred to Compare evidence mode",
                payload: { text: draft.text, assembly: draft.assembly },
              });
            }}
          >
            Send variant to Compare evidence
          </button>
      </section>
      {evidence ? (
        <section className="hs-panel science-surface">
          <MetricGrid
            metrics={[
              { id: "confidence", label: "Confidence score", value: pack.confidence_score },
              { id: "kind", label: "Kind", value: pack.kind },
            ]}
          />
          <DataTable rows={asList(pack.conflicts).map((row) => asRecord(row))} caption="Conflicts" testId="evidence-conflicts" />
          <ProvenancePanel envelope={evidence} extra={pack} />
        </section>
      ) : null}
    </div>
  );
}

function identityRecord(result: Record<string, unknown>): Record<string, unknown> {
  if (result.identity_hash || result.position_0based !== undefined) return result;
  const nested = asRecord(result.variant);
  if (nested.identity_hash || nested.position_0based !== undefined) return nested;
  const layer = asRecord(asRecord(asRecord(result.layers).variant_identity).data);
  if (layer.identity_hash || layer.position_0based !== undefined) return layer;
  return result;
}
