"use client";

import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, MetricGrid, ScientificError, ScientificWarning, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { ScienceFeatureMap } from "@/components/charts/ScienceFeatureMap";

export function MotifWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult } = useWorkspace();
  const draft = drafts.motif;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const envelope = results.motif;
  const result = asRecord(envelope?.result);
  const hits = asList(result.hits).map((row) => asRecord(Array.isArray(row) ? {} : row));
  const tracks = asList(result.sequence_feature_tracks);
  const computed = Boolean(envelope) && !error;

  async function run() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.motifSearch(
        { sequence: draft.sequence, pattern: draft.pattern, molecule: draft.molecule, include_explanation: true },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("motif", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Discovery</p>
      <h1>Motif Search</h1>
      <ScientificWarning>Motif occurrence is not proven biological function.</ScientificWarning>
      <CompactInput compact={computed} summary="Edit motif input">
        <div className="field">
          <label htmlFor="mol">Molecule</label>
          <select id="mol" value={draft.molecule} onChange={(e) => patchDraft("motif", { molecule: e.target.value as "DNA" | "RNA" | "PROTEIN" })}>
            <option>DNA</option>
            <option>RNA</option>
            <option>PROTEIN</option>
          </select>
        </div>
        <SequenceInput id="motif-seq" label="Sequence" value={draft.sequence} onChange={(sequence) => patchDraft("motif", { sequence })} />
        <div className="field">
          <label htmlFor="pattern">IUPAC motif</label>
          <input id="pattern" value={draft.pattern} onChange={(e) => patchDraft("motif", { pattern: e.target.value })} />
        </div>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void run()}>Search</button>
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No motif search yet" /> : (
        <>
          <AnalysisNav
            items={[
              { id: "motif-summary", label: "Overview" },
              { id: "motif-hits-section", label: "Hits" },
              { id: "motif-meaning-section", label: "Interpretation" },
              { id: "motif-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="motif-summary">
          <MetricGrid
            metrics={[
              { id: "n", label: "Hit count", value: hits.length || result.n_hits },
              { id: "pattern", label: "Pattern", value: result.pattern ?? draft.pattern },
              { id: "molecule", label: "Molecule", value: result.molecule ?? draft.molecule },
              { id: "strand", label: "Strand", value: result.strand ?? result.search_strand },
              { id: "hits-status", label: "Hits status", value: result.hits_status },
            ]}
          />
          <p data-testid="motif-meaning">{String(result.hit_meaning || "")}</p>
          <p data-testid="motif-pattern">pattern {String(result.pattern || draft.pattern)}</p>
          </ResultAnchor>
          <ResultAnchor id="motif-hits-section">
          {tracks.length > 0 ? (
            <InfographicFrame
              title="Motif hit track"
              dataSource="HelixScope Core motif hits"
              method="IUPAC pattern search on the submitted molecule string"
              units="coordinates"
              status={result.status}
              shows="Where the pattern matched the submitted sequence"
              doesNotShow="Biological function of the motif"
            >
              <ScienceFeatureMap
                length={Number(result.length || String(result.sequence || draft.sequence).length) || 0}
                tracks={tracks}
                testId="motif-feature-map"
              />
            </InfographicFrame>
          ) : null}
          <DataTable rows={hits} testId="motif-hits" emptyLabel="0 hits." />
          </ResultAnchor>
          <ResultAnchor id="motif-meaning-section">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="motif-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} />
        </>
      )}
    </div>
  );
}
