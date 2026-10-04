"use client";

import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { helixApi } from "@/lib/api/client";
import { jobScientificResult, pollJob } from "@/lib/api/jobs";
import { asList, asRecord, chartNumeric, formatScientificValue } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, ExternalDisclosure, JobState, ScientificError, ScientificWarning, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { downloadCsv, downloadJson, fastaRecord } from "@/lib/export/download";
import { ScienceBarChart } from "@/components/charts/ScienceBarChart";
import { ScienceFeatureMap } from "@/components/charts/ScienceFeatureMap";

const PCT_RANGE: [number, number] = [0, 100];

export function BlastWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.blast;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [jobId, setJobId] = useState("");
  const [execution, setExecution] = useState("");
  const [selected, setSelected] = useState<string[]>([]);
  const [transferNote, setTransferNote] = useState("");
  const envelope = results.blast;
  const availability = asRecord(results.blastAvailability?.result);
  const result = envelope ? jobScientificResult(envelope) : {};
  const hits = asList(asRecord(result).hits ?? asRecord(result).alignments).map((row) => asRecord(row));
  const labels = hits.map((hit, index) => String(hit.accession ?? hit.hit_id ?? hit.subject ?? `hit_${index + 1}`));
  const bitScores = hits.map((hit) => hit.bit_score);
  const evalues = hits.map((hit) => hit.evalue);
  const identities = hits.map((hit) => hit.identity_pct);
  const coverage = hits.map((hit) => hit.query_coverage_pct);
  const evalueNumbers = evalues.map((value) => chartNumeric(value));
  const logE = evalueNumbers.length > 0 && evalueNumbers.every((value) => value !== null && value > 0);
  const queryLength = Number(asRecord(result).query_length ?? 0);
  const spanTracks = [
    {
      name: "HSPs",
      features: hits.map((hit, index) => {
        const from = chartNumeric(hit.query_from);
        const to = chartNumeric(hit.query_to);
        const low = from !== null && to !== null ? Math.min(from, to) : null;
        const high = from !== null && to !== null ? Math.max(from, to) : null;
        return {
          start: low,
          end: high !== null && low !== null && high === low ? high + 1 : high,
          label: labels[index],
        };
      }),
    },
  ];

  useEffect(() => {
    void helixApi.blastAvailability().then(({ data }) => setResult("blastAvailability", data)).catch(setError);
  }, [setResult]);

  async function run() {
    const gate = next();
    setError(null);
    setSelected([]);
    setTransferNote("");
    try {
      const accepted = draft.remote
        ? (await helixApi.blastRemoteSubmit({ query: draft.query, email: draft.email, program: draft.program }, gate.signal)).data
        : (await helixApi.blastLocalSubmit({ query: draft.query, program: draft.program, molecule: "DNA" }, gate.signal)).data;
      setJobId(accepted.job_id);
      setExecution(accepted.execution_status);
      const done = await pollJob(accepted, { signal: gate.signal });
      if (gate.isCurrent()) {
        setExecution(String(done.execution_status || ""));
        setResult("blast", done);
      }
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function fetchHitsIntoMsa() {
    if (selected.length > 0 && !(draft.email || drafts.ncbi.email)) {
      setError(
        new Error("NCBI fetch needs an email. BLAST HSP alignments are not used as sequences."),
      );
      return;
    }
    setError(null);
    setTransferNote("");
    const parts = [fastaRecord("blast_query", draft.query)];
    const missing: string[] = [];
    for (const accession of selected) {
      try {
        const { data } = await helixApi.ncbiFetch({
          accession,
          email: draft.email || drafts.ncbi.email,
          db: drafts.ncbi.db || "nucleotide",
        });
        const record = asRecord(asRecord(data.result).record ?? data.result);
        const sequence = String(record.sequence || asRecord(data.result).sequence || "");
        if (!sequence) {
          missing.push(accession);
          continue;
        }
        parts.push(fastaRecord(String(record.accession || accession), sequence));
      } catch {
        missing.push(accession);
      }
    }
    patchDraft("msa", { fasta: parts.join(""), email: draft.email || drafts.ncbi.email });
    recordTransfer({
      kind: "blast-to-msa",
      note: "BLAST query plus NCBI records. HSP qseq/hseq was not used as a sequence.",
      payload: { accessions: selected.join(","), missing: missing.join(",") },
    });
    setTransferNote(
      missing.length
        ? `MSA draft received the query plus fetched records. NCBI unavailable for: ${missing.join(", ")}. HSP alignments are not full subjects.`
        : "MSA draft received the BLAST query plus NCBI subject records. HSP alignments were not treated as sequences.",
    );
    router.push("/evolution/msa");
  }

  return (
    <div className="hs-page">
      <p className="kicker">Data</p>
      <h1>BLAST</h1>
      <ScientificWarning>
        NCBI remote BLAST and local BLAST+ are different engines. A tiny local database is not nt/nr.
        BLAST hits are not MSA members until the NCBI record is retrieved by accession.
      </ScientificWarning>
      <p data-testid="blast-backend">
        {draft.remote ? "NCBI BLAST" : "LOCAL BLAST+"}
      </p>
      <p data-testid="blast-not-nt">
        local_is_not_nt_nr {formatScientificValue(availability.local_is_not_nt_nr ?? true)}
      </p>
      <p data-testid="blast-local-status">
        Local BLAST {formatScientificValue(asRecord(availability.local).available)} ·{" "}
        {String(asRecord(availability.local).reason || asRecord(availability.local).status || "")}
      </p>
      {draft.remote ? <ExternalDisclosure service="NCBI BLAST" /> : null}
      <CompactInput compact={Boolean(envelope)} summary="Edit BLAST query">
        <label>
          <input type="checkbox" checked={draft.remote} onChange={(e) => patchDraft("blast", { remote: e.target.checked })} /> NCBI remote BLAST
        </label>
        <SequenceInput id="blast-q" label="Query" value={draft.query} onChange={(query) => patchDraft("blast", { query })} />
        <div className="field">
          <label htmlFor="blast-email">Email (NCBI remote BLAST and BLAST-to-MSA fetch)</label>
          <input id="blast-email" value={draft.email} onChange={(e) => patchDraft("blast", { email: e.target.value })} />
        </div>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void run()}>Submit job</button>
      </CompactInput>
      {jobId ? <JobState execution={execution} jobId={jobId} /> : null}
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No BLAST job completed in this session" /> : (
        <>
          <AnalysisNav
            items={[
              { id: "blast-charts", label: "Hit charts", available: hits.length > 0 },
              { id: "blast-hits-section", label: "Hits" },
              { id: "blast-msa", label: "BLAST to MSA", available: hits.length > 0 },
              { id: "blast-raw", label: "Raw data" },
            ]}
          />
          {hits.length > 0 ? (
            <ResultAnchor id="blast-charts">
            <section className="hs-panel science-surface">
              <h2>Hit charts (API hit fields)</h2>
              <InfographicFrame
                title="Bit score ranking"
                dataSource="BLAST hit table (local BLAST+ or NCBI)"
                method="Reported bit scores; not recomputed in the browser"
                units="bit score"
                status={asRecord(result).status}
                shows="Relative ranking of returned hits by bit score"
                doesNotShow="That a tiny local database is nt/nr"
              >
                <ScienceBarChart
                  labels={labels}
                  values={bitScores}
                  orientation="h"
                  xTitle="Bit score (NCBI)"
                  yTitle="Subject accession"
                  title="BLAST hit ranking (NCBI bit score)"
                  testId="blast-hit-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="E-value"
                dataSource="BLAST hit table"
                method="Reported E-values from the engine XML/API"
                units="E-value"
                status={asRecord(result).status}
                shows="Reported chance of a hit this strong in a database of this size"
                doesNotShow="A substitute for identity. Zero may be XML underflow, not a computed zero."
              >
                <ScienceBarChart
                  labels={labels}
                  values={evalues}
                  logY={logE}
                  xTitle="Subject accession"
                  yTitle={logE ? "E-value (NCBI, log scale)" : "E-value (NCBI; 0 is XML underflow, not a substitute)"}
                  title="BLAST E-value (NCBI)"
                  testId="blast-evalue-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Identity percent"
                dataSource="BLAST HSP identity counts"
                method="identity_pct from the engine result"
                units="percent"
                status={asRecord(result).status}
                shows="Percent identity of each reported HSP"
                doesNotShow="The complete subject sequence"
              >
                <ScienceBarChart
                  labels={labels}
                  values={identities}
                  yRange={PCT_RANGE}
                  xTitle="Subject accession"
                  yTitle="Identity %"
                  title="BLAST identity percent (NCBI identity counts)"
                  testId="blast-identity-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Query coverage"
                dataSource="BLAST HSP query coordinates"
                method="query_coverage_pct from the engine result"
                units="percent"
                status={asRecord(result).status}
                shows="How much of the query is covered by each HSP"
                doesNotShow="Full-subject coverage or taxonomy"
              >
                <ScienceBarChart
                  labels={labels}
                  values={coverage}
                  yRange={PCT_RANGE}
                  xTitle="Subject accession"
                  yTitle="Query coverage %"
                  title="BLAST query coverage (NCBI HSP coordinates)"
                  testId="blast-coverage-chart"
                />
              </InfographicFrame>
              {queryLength > 0 ? (
                <InfographicFrame
                  title="HSP query intervals"
                  dataSource="BLAST HSP query_from / query_to"
                  method="1-based HSP coordinates from the engine"
                  units="query coordinates"
                  status={asRecord(result).status}
                  shows="Where each HSP sits on the query"
                  doesNotShow="The full subject sequence; HSP aligned strings are not MSA members"
                >
                  <ScienceFeatureMap
                    length={queryLength}
                    tracks={spanTracks}
                    title="BLAST HSP query intervals (NCBI 1-based HSP coordinates)"
                    testId="blast-hsp-overview"
                  />
                </InfographicFrame>
              ) : null}
            </section>
            </ResultAnchor>
          ) : null}
          <ResultAnchor id="blast-hits-section">
          <DataTable rows={hits} testId="blast-hits" />
          </ResultAnchor>
          {hits.length > 0 ? (
            <button type="button" className="btn-secondary" onClick={() => downloadCsv("blast_hits.csv", hits)}>
              Download hits (CSV)
            </button>
          ) : null}
          {hits.length > 0 ? (
            <ResultAnchor id="blast-msa">
            <section className="hs-panel science-surface">
              <h2>Fetch selected subjects into MSA set</h2>
              <ScientificWarning>
                BLAST HSP (qseq/hseq) is not the full subject sequence. Selected accessions are retrieved
                from NCBI. Hits without a retrieved record stay out of the MSA draft.
              </ScientificWarning>
              <ExternalDisclosure service="NCBI Entrez" />
              {labels.map((label) => (
                <label key={label}>
                  <input
                    type="checkbox"
                    checked={selected.includes(label)}
                    onChange={(e) =>
                      setSelected((prev) =>
                        e.target.checked ? [...prev, label] : prev.filter((item) => item !== label),
                      )
                    }
                  />{" "}
                  {label}
                </label>
              ))}
              <button
                type="button"
                className="btn-secondary"
                data-testid="blast-to-msa"
                disabled={scienceLocked || (!draft.query && selected.length === 0)}
                onClick={() => void fetchHitsIntoMsa()}
              >
                Fetch selected subjects into MSA set
              </button>
              {transferNote ? <p className="hs-topbar-meta">{transferNote}</p> : null}
            </section>
            </ResultAnchor>
          ) : null}
          <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-blast.json", result)}>
            Download BLAST JSON
          </button>
          <InterpretationPanel explanation={asRecord(result).explanation} />
          <ResultAnchor id="blast-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} />
        </>
      )}
    </div>
  );
}
