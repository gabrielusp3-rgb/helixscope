"use client";

import dynamic from "next/dynamic";
import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { jobScientificResult, pollJob } from "@/lib/api/jobs";
import { asList, asRecord } from "@/lib/api/numeric";
import { compareDisplayResult, overlayCoordinates } from "@/lib/compare/overlay";
import { downloadJson } from "@/lib/export/download";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, ExternalDisclosure, JobState, MetricGrid, ScientificError, ScientificWarning, SequenceInput } from "@/components/science/SciencePrimitives";
import { ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";

const StructureViewportPlotly = dynamic(
  () => import("@/components/structure/StructureViewportPlotly").then((m) => m.StructureViewportPlotly),
  { ssr: false },
);

type CompareMode = "structures" | "variants" | "proteins" | "guides" | "evolution" | "evidence";

export function CompareWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult } = useWorkspace();
  const draft = drafts.compare;
  const [error, setError] = useState<unknown>(null);
  const [execution, setExecution] = useState("");
  const [methods, setMethods] = useState<Array<Record<string, unknown>>>([]);
  const envelope = results.compare;
  const raw = asRecord(envelope?.result);
  const result = compareDisplayResult(raw);
  const tmScores = asList(result.tm_scores).map((row) => asRecord(row));
  const coverage = result.aln_coverage_percent;
  const overlay = overlayCoordinates(result);
  const fetchMeta = asRecord(result.coordinate_fetch);
  const [jobNote, setJobNote] = useState("");
  const mode = draft.mode;

  useEffect(() => {
    void helixApi
      .compareMethods()
      .then(({ data }) => setMethods(asList(asRecord(data.result).methods).map((row) => asRecord(row))))
      .catch(() => null);
  }, []);

  async function parse() {
    setError(null);
    setExecution("");
    try {
      const payload = JSON.parse(draft.json) as Record<string, unknown>;
      const { data } = await helixApi.compareParse({ payload });
      setResult("compare", data);
    } catch (err) {
      if (err instanceof SyntaxError) {
        setError({ message: "Alignment JSON is not valid JSON.", name: "SyntaxError" });
      } else {
        setError(err);
      }
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Structure</p>
      <h1>Structure Compare</h1>
      <ScientificWarning>
        Block RMSD, global RMSD, TM-score, coverage and aligned residue pairs are distinct API fields. They are never merged.
        Form values are sent unchanged. US-align accepts only bundled 1CRN, 1BNA, or 1RNA.
      </ScientificWarning>
      <div className="field">
        <label htmlFor="cmp-mode">Mode</label>
        <select
          id="cmp-mode"
          data-testid="compare-mode"
          value={mode}
          onChange={(e) => patchDraft("compare", { mode: e.target.value as CompareMode })}
        >
          <option value="structures">Structures</option>
          <option value="variants">Variants</option>
          <option value="proteins">Proteins</option>
          <option value="guides">Guides</option>
          <option value="evolution">Evolution</option>
          <option value="evidence">Evidence</option>
        </select>
      </div>
      {mode === "structures" ? (
      <section className="hs-panel glass-medium">
        <div className="field">
          <label htmlFor="ref">Reference</label>
          <input id="ref" data-testid="compare-reference" value={draft.reference} onChange={(e) => patchDraft("compare", { reference: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="tgt">Target</label>
          <input id="tgt" data-testid="compare-target" value={draft.target} onChange={(e) => patchDraft("compare", { target: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="ref-chain">Reference chain</label>
          <input id="ref-chain" data-testid="compare-ref-chain" value={draft.referenceChain} onChange={(e) => patchDraft("compare", { referenceChain: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="tgt-chain">Target chain</label>
          <input id="tgt-chain" data-testid="compare-tgt-chain" value={draft.targetChain} onChange={(e) => patchDraft("compare", { targetChain: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="cmp-method">RCSB method</label>
          <select id="cmp-method" data-testid="compare-method" value={draft.method} onChange={(e) => patchDraft("compare", { method: e.target.value })}>
            {(methods.length ? methods : [{ api_name: "tm-align", display_name: "TM-align" }]).map((row) => (
              <option key={String(row.api_name)} value={String(row.api_name)}>
                {String(row.display_name || row.api_name)}
              </option>
            ))}
          </select>
        </div>
        <SequenceInput id="cmp-json" label="Alignment JSON payload (API parse)" value={draft.json} onChange={(json) => patchDraft("compare", { json })} rows={10} />
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void parse()}>Parse alignment</button>
        <ExternalDisclosure service="RCSB Alignment API">Remote jobs are optional and explicit.</ExternalDisclosure>
        <button
          type="button"
          className="btn-secondary"
          data-testid="compare-rcsb-submit"
          onClick={async () => {
            setError(null);
            try {
              const { data } = await helixApi.compareRemoteSubmit({
                reference_entry: draft.reference,
                target_entry: draft.target,
                reference_chain: draft.referenceChain,
                target_chain: draft.targetChain,
                method: draft.method,
                fetch_coordinates: true,
              });
              setJobNote(`RCSB job ${data.job_id} ${data.execution_status}`);
              setExecution(data.execution_status);
              const done = await pollJob(data);
              setExecution(String(done.execution_status || ""));
              setResult("compare", { ...done, result: jobScientificResult(done) });
            } catch (err) {
              setError(err);
            }
          }}
        >
          Submit RCSB remote job
        </button>
        <button
          type="button"
          className="btn-secondary"
          data-testid="compare-usalign-submit"
          onClick={async () => {
            setError(null);
            try {
              const { data } = await helixApi.compareUSalignSubmit({
                reference_entry: draft.reference as "1CRN" | "1BNA" | "1RNA",
                target_entry: draft.target as "1CRN" | "1BNA" | "1RNA",
                reference_chain: draft.referenceChain,
                target_chain: draft.targetChain,
              });
              setJobNote(`US-align job ${data.job_id} ${data.execution_status}`);
              setExecution(data.execution_status);
              const done = await pollJob(data);
              setExecution(String(done.execution_status || ""));
              setResult("compare", { ...done, result: jobScientificResult(done) });
            } catch (err) {
              setError(err);
            }
          }}
        >
          Submit US-align job
        </button>
        {jobNote ? <p data-testid="compare-job">{jobNote}</p> : null}
        {execution ? <JobState execution={execution} /> : null}
      </section>
      ) : null}
      {mode === "variants" || mode === "evidence" ? (
        <section className="hs-panel glass-medium">
          <div className="field">
            <label htmlFor="var-a">Variant A</label>
            <input id="var-a" data-testid="compare-text-a" value={draft.textA} onChange={(e) => patchDraft("compare", { textA: e.target.value })} />
          </div>
          <div className="field">
            <label htmlFor="asm-a">Assembly A</label>
            <input id="asm-a" value={draft.assemblyA} onChange={(e) => patchDraft("compare", { assemblyA: e.target.value })} />
          </div>
          <div className="field">
            <label htmlFor="var-b">Variant B</label>
            <input id="var-b" data-testid="compare-text-b" value={draft.textB} onChange={(e) => patchDraft("compare", { textB: e.target.value })} />
          </div>
          <div className="field">
            <label htmlFor="asm-b">Assembly B</label>
            <input id="asm-b" value={draft.assemblyB} onChange={(e) => patchDraft("compare", { assemblyB: e.target.value })} />
          </div>
          <button
            type="button"
            className="btn-primary"
            data-testid="compare-variants-run"
            disabled={scienceLocked}
            onClick={async () => {
              setError(null);
              try {
                const body = {
                  text_a: draft.textA,
                  assembly_a: draft.assemblyA,
                  text_b: draft.textB,
                  assembly_b: draft.assemblyB,
                };
                const { data } = mode === "evidence" ? await helixApi.compareEvidence(body) : await helixApi.compareVariants(body);
                setResult("compare", data);
              } catch (err) {
                setError(err);
              }
            }}
          >
            {mode === "evidence" ? "Compare evidence" : "Compare variants"}
          </button>
        </section>
      ) : null}
      {mode === "proteins" ? (
        <section className="hs-panel glass-medium">
          <SequenceInput id="prot-a" label="Protein A" value={draft.proteinA} onChange={(proteinA) => patchDraft("compare", { proteinA })} />
          <SequenceInput id="prot-b" label="Protein B" value={draft.proteinB} onChange={(proteinB) => patchDraft("compare", { proteinB })} />
          <button
            type="button"
            className="btn-primary"
            data-testid="compare-proteins-run"
            disabled={scienceLocked}
            onClick={async () => {
              setError(null);
              try {
                const { data } = await helixApi.compareProteins({
                  sequence_a: draft.proteinA,
                  sequence_b: draft.proteinB,
                  identifier_a: "protein_a",
                  identifier_b: "protein_b",
                });
                setResult("compare", data);
              } catch (err) {
                setError(err);
              }
            }}
          >
            Compare proteins
          </button>
        </section>
      ) : null}
      {mode === "guides" ? (
        <section className="hs-panel glass-medium">
          <SequenceInput id="guide-a" label="Guide A sequence" value={draft.guideA} onChange={(guideA) => patchDraft("compare", { guideA })} />
          <SequenceInput id="guide-b" label="Guide B sequence" value={draft.guideB} onChange={(guideB) => patchDraft("compare", { guideB })} />
          <button
            type="button"
            className="btn-primary"
            data-testid="compare-guides-run"
            disabled={scienceLocked}
            onClick={async () => {
              setError(null);
              try {
                const { data } = await helixApi.compareGuides({
                  guide_a: { guide_sequence: draft.guideA },
                  guide_b: { guide_sequence: draft.guideB },
                });
                setResult("compare", data);
              } catch (err) {
                setError(err);
              }
            }}
          >
            Compare guides
          </button>
        </section>
      ) : null}
      {mode === "evolution" ? (
        <section className="hs-panel glass-medium">
          <SequenceInput id="evo-fasta" label="Prealigned FASTA" value={draft.evolutionFasta} onChange={(evolutionFasta) => patchDraft("compare", { evolutionFasta })} rows={10} />
          <button
            type="button"
            className="btn-primary"
            data-testid="compare-evolution-run"
            disabled={scienceLocked}
            onClick={async () => {
              setError(null);
              try {
                const { data } = await helixApi.compareEvolution({ fasta: draft.evolutionFasta, column: 0 });
                setResult("compare", data);
              } catch (err) {
                setError(err);
              }
            }}
          >
            Inspect MSA evolution
          </button>
        </section>
      ) : null}
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No compare result" /> : (
        <>
          {mode === "structures" ? (
          <>
          <div data-testid="compare-hero">
            <MetricGrid
              metrics={[
                { id: "block", label: "Block RMSD", value: result.rmsd_block0_angstrom, unit: "A" },
                { id: "global", label: "Global RMSD", value: result.rmsd_global_angstrom, unit: "A" },
                { id: "tm", label: "TM-score", value: tmScores[0]?.value ?? result.tm_score },
                { id: "cov", label: "Coverage", value: Array.isArray(coverage) ? coverage.join(" / ") : coverage, unit: "%" },
                { id: "pairs", label: "Aligned residue pairs", value: result.n_aligned_residue_pairs },
              ]}
            />
          </div>
          <p>
            {String(asRecord(result.reference).entry_id)} vs {String(asRecord(result.target).entry_id)} · visual{" "}
            {String(result.superposition_visual)} · original coordinates preserved{" "}
            {String(result.original_coordinates_preserved)}
          </p>
          <p data-testid="compare-overlay-status">
            Overlay {String(fetchMeta.overlay_status || (overlay.length > 0 ? "AVAILABLE" : "UNAVAILABLE"))}.{" "}
            {String(fetchMeta.note || "")}
          </p>
          {overlay.length > 0 ? (
            <StructureViewportPlotly
              scene={{
                coordinates: overlay,
                kind: "experimental",
                status: "EXPERIMENTAL",
                structure_id: `${String(asRecord(result.reference).entry_id || draft.reference)}-vs-${String(asRecord(result.target).entry_id || draft.target)}`,
                source: "compare superposition copy",
                mapping_status: result.mapping_status,
              }}
            />
          ) : (
            <p data-testid="compare-overlay-empty">
              Cartesian overlay is empty. Residue pairs are not converted into invented coordinates.
              Overlay XYZ appears only when Core superposition_bundle returns real atom copies.
            </p>
          )}
          <DataTable rows={asList(result.residue_pairs).map((row) => asRecord(row))} />
          </>
          ) : (
            <p data-testid="compare-mode-result">mode {String(raw.mode || mode)} · ranking {String(raw.ranking ?? "none")}</p>
          )}
          <button
            type="button"
            className="btn-secondary"
            onClick={() => downloadJson("helixscope-compare.json", result)}
          >
            Download compare JSON
          </button>
          <ProvenancePanel envelope={envelope} extra={result} />
        </>
      )}
    </div>
  );
}
