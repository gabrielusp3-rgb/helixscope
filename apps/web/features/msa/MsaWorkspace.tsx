"use client";

import { useRouter } from "next/navigation";
import { useEffect, useState } from "react";
import { helixApi, uploadPrealignedFasta } from "@/lib/api/client";
import { jobScientificResult, pollJob } from "@/lib/api/jobs";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, ExternalDisclosure, JobState, ScientificError, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, ProvenancePanel } from "@/components/science/Provenance";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { MsaViewer } from "@/components/evolution/MsaViewer";
import { downloadText, fastaRecord } from "@/lib/export/download";
import { ScienceChart } from "@/components/charts/ScienceChart";

export function MsaWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.msa;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [execution, setExecution] = useState("");
  const envelope = results.msa;
  const result = asRecord(envelope?.result);
  const availability = asRecord(results.msaAvailability?.result);

  useEffect(() => {
    void helixApi.msaAvailability().then(({ data }) => setResult("msaAvailability", data)).catch(setError);
  }, [setResult]);

  async function importPrealigned() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.msaPrealigned({ fasta: draft.fasta, format: draft.format }, gate.signal);
      if (gate.isCurrent()) setResult("msa", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function runEngine() {
    const gate = next();
    setError(null);
    try {
      const { data: accepted } = await helixApi.msaEngineSubmit(
        { fasta: draft.fasta, email: draft.email, backend: draft.backend },
        gate.signal,
      );
      setExecution(accepted.execution_status);
      const done = await pollJob(accepted, { signal: gate.signal });
      if (gate.isCurrent()) {
        setExecution(String(done.execution_status || ""));
        setResult("msa", { ...done, result: jobScientificResult(done) });
      }
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Evolution</p>
      <h1>MSA</h1>
      <p data-testid="msa-availability">
        MSA backends {formatScientificValue(availability.status ?? asRecord(availability.ebi_clustalo).available)}
      </p>
      <CompactInput compact={Boolean(envelope)} summary="Edit MSA input">
        <SequenceInput id="msa-fasta" label="Alignment text" value={draft.fasta} onChange={(fasta) => patchDraft("msa", { fasta })} rows={12} />
        <div className="field">
          <label htmlFor="msa-format">Input format</label>
          <select id="msa-format" data-testid="msa-format" value={draft.format} onChange={(e) => patchDraft("msa", { format: e.target.value as typeof draft.format })}>
            <option value="auto">auto (FASTA / Clustal / Stockholm / PHYLIP)</option>
            <option value="fasta">FASTA (aligned)</option>
            <option value="clustal">Clustal</option>
            <option value="stockholm">Stockholm</option>
            <option value="phylip">PHYLIP</option>
          </select>
        </div>
        <div className="field">
          <label htmlFor="msa-file">Bounded upload</label>
          <input
            id="msa-file"
            type="file"
            accept=".fa,.fasta,.txt"
            onChange={(event) => {
              const file = event.target.files?.[0];
              if (!file) return;
              void uploadPrealignedFasta(file).then(({ data }) => setResult("msa", data)).catch(setError);
            }}
          />
        </div>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void importPrealigned()}>Import prealigned FASTA</button>
        <ExternalDisclosure service="EBI Clustal Omega or a local MSA binary">Engine jobs are optional.</ExternalDisclosure>
        <div className="field">
          <label htmlFor="email-msa">Email (EBI)</label>
          <input id="email-msa" value={draft.email} onChange={(e) => patchDraft("msa", { email: e.target.value })} />
        </div>
        <button type="button" className="btn-secondary" onClick={() => void runEngine()}>Submit engine job</button>
        {execution ? <JobState execution={execution} /> : null}
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No MSA loaded" /> : (
        <>
          <AnalysisNav
            items={[
              { id: "msa-summary", label: "Overview" },
              { id: "msa-conservation", label: "Conservation", available: asList(asRecord(result.conservation).scores).length > 0 },
              { id: "msa-viewer", label: "Viewer" },
              { id: "msa-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="msa-summary">
          <p data-testid="msa-hash">alignment_hash {formatScientificValue(result.alignment_hash)}</p>
          <p>n_sequences {formatScientificValue(result.n_sequences)} · status {String(result.status)}</p>
          </ResultAnchor>
          {asList(asRecord(result.conservation).scores).length > 0 ? (
            <ResultAnchor id="msa-conservation">
            <section className="hs-panel science-surface">
              <h2>Conservation</h2>
              <InfographicFrame
                title="Column conservation"
                dataSource="HelixScope Core MSA conservation.scores"
                method="Core conservation scores; not recomputed in the browser"
                units="conservation score"
                status={result.status}
                shows="Per-column conservation of the loaded alignment"
                doesNotShow="Functional sites or phylogenetic support"
              >
                <ScienceChart
                  x={asList(asRecord(result.conservation).scores).map((_, index) => index)}
                  y={asList(asRecord(result.conservation).scores)}
                  xTitle="column"
                  yTitle="conservation"
                  title="MSA conservation (API scores)"
                  testId="msa-conservation"
                />
              </InfographicFrame>
            </section>
            </ResultAnchor>
          ) : null}
          <ResultAnchor id="msa-viewer">
          <MsaViewer result={result} />
          </ResultAnchor>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              const aligned = asList(result.rows)
                .map((row) => asRecord(row))
                .map((row) => fastaRecord(String(row.identifier || "seq"), String(row.aligned || row.sequence || "")))
                .join("");
              const fasta = aligned || String(result.raw_alignment || draft.fasta);
              patchDraft("phylogeny", { fasta, alignmentHash: String(result.alignment_hash || "") });
              recordTransfer({
                kind: "msa-to-phylogeny",
                note: "Typed MSA identity transferred for tree inference",
                payload: {
                  alignment_hash: String(result.alignment_hash || ""),
                  molecule: String(result.molecule || ""),
                  engine: String(result.tool || result.method || result.engine || ""),
                },
              });
              router.push("/evolution/phylogeny");
            }}
          >
            Send MSA to Phylogeny
          </button>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              const fasta = asList(result.rows)
                .map((row) => asRecord(row))
                .map((row) => fastaRecord(String(row.identifier || "seq"), String(row.aligned || "")))
                .join("");
              downloadText("helixscope-msa.fasta", fasta);
            }}
          >
            Download aligned FASTA
          </button>
          <ResultAnchor id="msa-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <InterpretationPanel explanation={result.explanation} />
          <ProvenancePanel envelope={envelope} extra={result} />
        </>
      )}
    </div>
  );
}
