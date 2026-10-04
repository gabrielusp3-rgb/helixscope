"use client";

import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, MetricGrid, ScientificError, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { ScienceHeatmap } from "@/components/charts/ScienceHeatmap";
import { ScienceColumnMap } from "@/components/charts/ScienceColumnMap";
import { downloadText } from "@/lib/export/download";

export function AlignmentWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult } = useWorkspace();
  const draft = drafts.alignment;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const envelope = results.alignment;
  const result = asRecord(envelope?.result);
  const computed = Boolean(envelope) && !error;

  async function run() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.alignmentPairwise(
        { seq1: draft.seq1, seq2: draft.seq2, mode: draft.mode, include_explanation: true, include_dotplot: true, translate_nucleic: draft.translateNucleic },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("alignment", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Analysis</p>
      <h1>Alignment</h1>
      <CompactInput compact={computed} summary="Edit alignment input">
        <SequenceInput id="seq1" label="Sequence A" value={draft.seq1} onChange={(seq1) => patchDraft("alignment", { seq1 })} />
        <SequenceInput id="seq2" label="Sequence B" value={draft.seq2} onChange={(seq2) => patchDraft("alignment", { seq2 })} />
        <div className="field">
          <label htmlFor="mode">Method</label>
          <select id="mode" value={draft.mode} onChange={(e) => patchDraft("alignment", { mode: e.target.value as "global" | "local" })}>
            <option value="global">global (Needleman-Wunsch)</option>
            <option value="local">local (Smith-Waterman)</option>
          </select>
        </div>
        <label>
          <input
            type="checkbox"
            data-testid="align-translate"
            checked={draft.translateNucleic}
            onChange={(e) => patchDraft("alignment", { translateNucleic: e.target.checked })}
          />{" "}
          Translate nucleic acid before protein alignment (frame +1, stop at first stop)
        </label>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void run()}>Run</button>
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No alignment yet" /> : (
        <>
          <AnalysisNav
            items={[
              { id: "aln-summary", label: "Overview" },
              { id: "aln-map", label: "Column map", available: asList(result.column_classes).length > 0 },
              { id: "aln-dotplot", label: "Dot plot", available: Boolean(result.dotplot) },
              { id: "aln-meaning", label: "Interpretation" },
              { id: "aln-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="aln-summary">
          <MetricGrid
            metrics={[
              { id: "method", label: "Method", value: result.method },
              { id: "id", label: "Identity", value: result.identity_pct, unit: "%" },
              { id: "score", label: "Score", value: result.score },
              { id: "gaps", label: "Gaps", value: result.gaps },
              { id: "len", label: "Length", value: result.length ?? result.alignment_length },
            ]}
          />
          <pre className="mono" data-testid="aligned-strings">
            {String(result.aligned_seq1 || "")}
            {"\n"}
            {String(result.aligned_seq2 || "")}
          </pre>
          </ResultAnchor>
          {asList(result.column_classes).length > 0 ? (
            <ResultAnchor id="aln-map">
            <section className="hs-panel science-surface">
              <h2>Column class map</h2>
              <InfographicFrame
                title="Match / mismatch / gap"
                dataSource="HelixScope Core classify_alignment_columns"
                method={String(result.method || "pairwise alignment")}
                units="column class"
                status={result.status}
                shows="Which aligned columns are match, mismatch or gap"
                doesNotShow="Homology, phylogeny, or structural similarity"
              >
                <ScienceColumnMap
                  classes={asList(result.column_classes)}
                  title="Match / mismatch / gap (API classify_alignment_columns)"
                  testId="alignment-column-map"
                />
              </InfographicFrame>
            </section>
            </ResultAnchor>
          ) : null}
          {result.dotplot ? (
            <ResultAnchor id="aln-dotplot">
            <section className="hs-panel science-surface">
              <h2>Dot plot</h2>
              <InfographicFrame
                title="Pairwise dot plot"
                dataSource="HelixScope Core alignment dotplot"
                method="Identity matrix of the compared residues from Core"
                units="match indicator"
                status={result.status}
                shows="Where the two strings share identical residues"
                doesNotShow="A substitute identity percent or E-value"
              >
                <ScienceHeatmap matrix={result.dotplot} title="Pairwise dot plot (API matrix)" testId="alignment-dotplot" />
              </InfographicFrame>
            </section>
            </ResultAnchor>
          ) : null}
          <button
            type="button"
            className="btn-secondary"
            onClick={() =>
              downloadText(
                "helixscope-alignment.txt",
                `>${String(result.method || "alignment")}_1\n${String(result.aligned_seq1 || "")}\n>${String(result.method || "alignment")}_2\n${String(result.aligned_seq2 || "")}\n`,
              )
            }
          >
            Download aligned sequences
          </button>
          <ResultAnchor id="aln-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="aln-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} />
        </>
      )}
    </div>
  );
}
