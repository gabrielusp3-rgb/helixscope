"use client";

import dynamic from "next/dynamic";
import { useRouter } from "next/navigation";
import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, asRecordRows, chartNumeric, formatScientificValue } from "@/lib/api/numeric";
import { dinucleotideOeMatrix } from "@/lib/charts/dinucleotide";
import { codonSynonymousMatrix } from "@/lib/charts/codonUsage";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import {
  EmptyState,
  InvalidMoleculeResult,
  MetricGrid,
  ScientificError,
  ScientificStatus,
  ScientificWarning,
  SequenceInput,
} from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput, FastaFileField, SequenceViewer } from "@/components/science/ResultChrome";
import { AnalysisNav, AnonymousSequenceNotice, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { ScienceArcDiagram } from "@/components/charts/ScienceArcDiagram";
import { ScienceCircularPairs } from "@/components/charts/ScienceCircularPairs";
import { ScienceBarChart } from "@/components/charts/ScienceBarChart";
import { ScienceChart } from "@/components/charts/ScienceChart";
import { ScienceHeatmap } from "@/components/charts/ScienceHeatmap";
import { downloadCsv, downloadJson, downloadText, fastaRecord } from "@/lib/export/download";
import type { Orbit } from "@/lib/gesture/math";

const StructureViewportPlotly = dynamic(
  () => import("@/components/structure/StructureViewportPlotly").then((m) => m.StructureViewportPlotly),
  { ssr: false },
);

const DEFAULT_ORBIT: Orbit = { yaw: 0.7, pitch: 0.35, radius: 2.4 };

export function RnaWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.rna;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [loading, setLoading] = useState(false);
  const [helix, setHelix] = useState<Record<string, unknown> | null>(null);
  const [orbit, setOrbit] = useState<Orbit>(DEFAULT_ORBIT);
  const [helixLod, setHelixLod] = useState("high");
  const [helixStart, setHelixStart] = useState(0);
  const [helixEnd, setHelixEnd] = useState(400);
  const [helixInspect, setHelixInspect] = useState(0);
  const [helixCenter, setHelixCenter] = useState(false);
  const [helixPrefer, setHelixPrefer] = useState("illustrative");
  const envelope = results.rna;
  const translateEnvelope = results.rnaTranslate;
  const result = asRecord(envelope?.result);
  const translated = asRecord(translateEnvelope?.result);
  const cai = asRecord(result.cai);
  const enc = asRecord(result.enc);
  const fold = asRecord(result.fold ?? result.vienna ?? result.structure);
  const rscuRows = asRecordRows(result.rscu, "codon");
  const codonRows = asRecordRows(result.codon_usage, "Codon");
  const compositionRows = asRecordRows(result.composition, "symbol").filter((row) => Number(row.count) > 0);
  const dinucRows = asRecordRows(result.dinucleotides, "dinucleotide");
  const dinucHeat = dinucleotideOeMatrix(result.dinucleotides);
  const profiles = asList(result.windowed_profiles).map((row) => asRecord(row));
  const kmers = asRecord(result.kmer_summary);
  const kmerTop = asList(kmers.top).map((row) => asRecord(row));
  const status = String(result.status || "");
  const invalid = status === "ERROR";
  const computed = status === "COMPUTED";
  const busy = loading || scienceLocked;
  const codonHeat = codonSynonymousMatrix(result.codon_usage);
  const seqLen = Number(result.length) || 0;
  const helixViewer = asRecord(helix?.viewer);
  const x = profiles.map((row, index) => chartNumeric(row.midpoint ?? row.center ?? row.position) ?? index);
  const compositionBar = compositionRows.filter((row) => ["A", "U", "T", "C", "G"].includes(String(row.symbol)));

  async function run() {
    const gate = next();
    setLoading(true);
    setError(null);
    setResult("rna", undefined);
    setHelix(null);
    try {
      const { data } = await helixApi.rnaAnalyze(
        {
          sequence: draft.sequence,
          fold: draft.fold,
          include_codon_metrics: draft.codon,
          include_profiles: draft.includeProfiles,
          include_kmers: draft.includeKmers,
          include_explanation: true,
          backend: "auto",
          profile_window: 100,
          profile_step: 50,
          kmer_k: 3,
        },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("rna", data);
      const nextLen = Number(asRecord(data.result).length) || 0;
      setHelixStart(0);
      setHelixEnd(nextLen || 1);
      setHelixInspect(0);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    } finally {
      if (gate.isCurrent()) setLoading(false);
    }
  }

  return (
    <div className="hs-page">
      <p className="kicker">Analysis</p>
      <h1>RNA</h1>
      <p className="hs-lede">
        Codon metrics, RSCU, ENC, CAI, dinucleotide O/E and optional ViennaRNA folding come from the API
        only. MFE is PREDICTED secondary structure, never experimental 3D. Arbitrary RNA is not assumed
        to be coding.
      </p>
      {draft.codon ? (
        <ScientificWarning>
          Codon usage, RSCU, ENC and CAI are coding-sequence metrics. They are not proof that this RNA is
          a CDS. Use Translate coding composition for an explicit coding-context check.
        </ScientificWarning>
      ) : null}
      <CompactInput compact={computed} summary="Edit RNA input">
        <SequenceInput
          id="rna-seq"
          label="RNA sequence"
          hint="Canonical RNA (A, U, C, G) or single-record FASTA. DNA with T is rejected here; use Transcribe DNA to RNA first."
          value={draft.sequence}
          onChange={(sequence) => patchDraft("rna", { sequence })}
          disabled={busy}
        />
        <FastaFileField id="rna-file" disabled={busy} onText={(sequence) => patchDraft("rna", { sequence })} />
        <label>
          <input type="checkbox" checked={draft.fold} onChange={(e) => patchDraft("rna", { fold: e.target.checked })} /> Fold (ViennaRNA if installed)
        </label>
        <label>
          <input type="checkbox" checked={draft.codon} onChange={(e) => patchDraft("rna", { codon: e.target.checked })} /> Codon metrics
        </label>
        <label>
          <input type="checkbox" checked={draft.includeProfiles} onChange={(e) => patchDraft("rna", { includeProfiles: e.target.checked })} /> Windowed profiles
        </label>
        <label>
          <input type="checkbox" checked={draft.includeKmers} onChange={(e) => patchDraft("rna", { includeKmers: e.target.checked })} /> k-mers
        </label>
        <div className="hs-actions">
          <button type="button" className="btn-primary" disabled={busy} onClick={() => void run()}>Analyze</button>
          <button
            type="button"
            className="btn-secondary"
            data-testid="rna-translate"
            disabled={busy}
            onClick={async () => {
              const gate = next();
              setError(null);
              try {
                const { data } = await helixApi.rnaTranslateCoding({ sequence: draft.sequence }, gate.signal);
                if (gate.isCurrent()) setResult("rnaTranslate", data);
              } catch (err) {
                if (gate.isCurrent()) setError(err);
              }
            }}
          >
            Translate coding composition
          </button>
          <button type="button" className="btn-ghost" onClick={() => patchDraft("rna", { sequence: "" })}>Clear</button>
          <button
            type="button"
            className="btn-secondary"
            onClick={async () => {
              const gate = next();
              setError(null);
              try {
                const { data } = await helixApi.rnaTranscribe({ sequence: draft.sequence }, gate.signal);
                const rna = String(asRecord(data.result).rna || "");
                if (gate.isCurrent()) patchDraft("rna", { sequence: rna });
              } catch (err) {
                if (gate.isCurrent()) setError(err);
              }
            }}
          >
            Transcribe DNA to RNA
          </button>
          {computed ? (
            <>
              <button
                type="button"
                className="btn-secondary"
                onClick={() =>
                  downloadText(
                    "helixscope-rna.fasta",
                    fastaRecord(String(result.sequence_hash || "rna"), String(result.sequence || draft.sequence)),
                  )
                }
              >
                Download FASTA
              </button>
              <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-rna.json", result)}>
                Download JSON
              </button>
            </>
          ) : null}
        </div>
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {translateEnvelope ? (
        <section className="hs-panel science-surface" data-testid="rna-translate-result">
          <h2>Translated protein composition (coding context)</h2>
          <p>
            status {String(translated.status)} · context {String(translated.coding_context || translated.method || "")} ·
            frame {String(translated.frame ?? "")}
          </p>
          <p className="mono">{String(translated.translation || "")}</p>
          <p className="hs-lede">Arbitrary RNA is not assumed to be a gene. Only resolved CDS/ORFs are translated.</p>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              const protein = String(translated.translation || "");
              if (!protein) return;
              patchDraft("protein", { sequence: protein });
              recordTransfer({
                kind: "rna-to-protein",
                note: "Coding-context translation transferred to Protein analysis",
                payload: { coding_context: String(translated.coding_context || "") },
              });
              router.push("/analysis/protein");
            }}
          >
            Send translation to Protein
          </button>
        </section>
      ) : null}
      {invalid ? <InvalidMoleculeResult expected="canonical RNA (A, U, C, G), optionally as single-record FASTA" result={result} /> : null}
      {!envelope && !loading && !error ? <EmptyState title="No RNA result yet" /> : null}
      {computed ? (
        <>
          <AnalysisNav
            items={[
              { id: "rna-summary", label: "Overview" },
              { id: "rna-composition", label: "Composition", available: compositionBar.length > 0 },
              { id: "rna-codon", label: "Codon", available: codonRows.length > 0 },
              { id: "rna-rscu", label: "RSCU", available: rscuRows.length > 0 },
              { id: "rna-fold", label: "Folding", available: Boolean(draft.fold || fold.status || fold.dot_bracket) },
              { id: "rna-meaning", label: "Interpretation" },
              { id: "rna-raw", label: "Raw data" },
              { id: "rna-provenance", label: "Provenance" },
            ]}
          />
          <AnonymousSequenceNotice />
          <ResultAnchor id="rna-summary">
          <MetricGrid
            metrics={[
              { id: "len", label: "Length", value: result.length },
              { id: "gc", label: "GC", value: result.gc_content, unit: "%" },
              { id: "au", label: "AU", value: result.au_content, unit: "%" },
              { id: "entropy", label: "Entropy", value: result.shannon_entropy },
              { id: "cai", label: "CAI", value: cai.value },
              { id: "enc", label: "ENC", value: enc.value },
              { id: "cai_status", label: "CAI status", value: cai.status },
              { id: "enc_status", label: "ENC status", value: enc.status },
            ]}
          />
          <p className="hs-topbar-meta">
            Codon metrics, RSCU, ENC and CAI are codon-table statistics on the submitted RNA string.
            They are not proof that the paste is a real coding sequence. {String(result.identity_note || "")}
          </p>
          </ResultAnchor>
          {compositionBar.length > 0 ? (
            <ResultAnchor id="rna-composition">
            <section className="hs-panel science-surface">
              <h2>Composition</h2>
              <InfographicFrame
                title="RNA nucleotide composition"
                dataSource="HelixScope Core RNA composition"
                method="Direct A/U/C/G counts on the normalized RNA string"
                units="count"
                status={result.status}
                shows="Nucleotide counts of the analyzed RNA"
                doesNotShow="Transcript identity or experimental expression"
              >
                <ScienceBarChart
                  labels={compositionBar.map((row) => row.symbol)}
                  values={compositionBar.map((row) => row.count)}
                  title="RNA nucleotide counts (API composition)"
                  testId="rna-composition-bar"
                />
              </InfographicFrame>
              <DataTable rows={compositionRows} caption="Nucleotide composition" testId="rna-composition" />
            </section>
            </ResultAnchor>
          ) : null}
          {dinucRows.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Dinucleotide O/E</h2>
              <p className="hs-topbar-meta">
                Core dinucleotide_frequencies treats U as T, same as the legacy workstation. Axes are A,C,G,T.
              </p>
              <InfographicFrame
                title="RNA dinucleotide O/E"
                dataSource="HelixScope Core dinucleotide_frequencies"
                method="U counted as T; observed/expected pair ratios"
                units="O/E"
                status={result.status}
                shows="Adjacent-base over/under-representation"
                doesNotShow="Secondary-structure pairing from this heatmap alone"
              >
                <ScienceHeatmap
                  z={dinucHeat.z}
                  x={dinucHeat.x}
                  y={dinucHeat.y}
                  zmid={1}
                  title="RNA dinucleotide observed/expected (API; U counted as T)"
                  testId="rna-dinuc-heatmap"
                />
              </InfographicFrame>
              <DataTable rows={dinucRows} caption="Dinucleotides" />
            </section>
          ) : null}
          {profiles.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Windowed GC</h2>
              <ScienceChart
                x={x}
                y={profiles.map((row) => row.gc_percent)}
                xTitle="position"
                yTitle="GC %"
                title="Windowed GC on RNA (API arrays; U treated as T)"
                testId="rna-gc-chart"
              />
            </section>
          ) : null}
          {draft.fold || fold.status || fold.mfe_kcal_mol !== undefined || fold.dot_bracket ? (
            <ResultAnchor id="rna-fold">
            <section className="hs-panel science-surface" data-testid="rna-fold">
              <h2>Fold</h2>
              <ScientificWarning>RNA 2D prediction is PREDICTED, never experimental structure.</ScientificWarning>
              <p data-testid="fold-status">
                Fold status {formatScientificValue(fold.status)} <ScientificStatus status="PREDICTED" />
              </p>
              <p>MFE {formatScientificValue(fold.mfe_kcal_mol ?? result.mfe_kcal_mol)} kcal/mol</p>
              <pre className="mono">{String(fold.dot_bracket ?? fold.structure ?? result.dot_bracket ?? "")}</pre>
              {asList(fold.base_pairs).length > 0 ? (
                <InfographicFrame
                  title="MFE base-pair arcs"
                  dataSource="ViennaRNA Python path via HelixScope Core"
                  method="Minimum free energy secondary structure; layout from Core pair list"
                  units="kcal/mol; pair indices"
                  status="PREDICTED"
                  shows="Predicted paired positions from the returned dot-bracket"
                  doesNotShow="RNA 3D coordinates or experimental structure"
                >
                  <ScienceArcDiagram
                    length={Number(result.length) || 0}
                    pairs={fold.base_pairs}
                    testId="rna-fold-arcs"
                  />
                </InfographicFrame>
              ) : null}
              {asList(asRecord(fold.circular_layout).x).length > 0 ? (
                <ScienceCircularPairs layout={fold.circular_layout} testId="rna-fold-circular" />
              ) : null}
            </section>
            </ResultAnchor>
          ) : null}
          {codonRows.length > 0 ? (
            <ResultAnchor id="rna-codon">
            <section className="hs-panel science-surface">
              <h2>Codon usage</h2>
              <InfographicFrame
                title="Codon Absolute_pct"
                dataSource="HelixScope Core codon_usage"
                method="Frame +1 codon table on the submitted RNA; not a validated CDS unless labeled"
                units="percent of codon table"
                status={result.status}
                shows="Relative abundance of each codon in the codon table of this string"
                doesNotShow="That this RNA is a real coding sequence or an organism codon preference"
              >
                <ScienceBarChart
                  labels={codonRows.map((row) => row.Codon ?? row.codon)}
                  values={codonRows.map((row) => row.Absolute_pct ?? row.Count)}
                  title="Codon Absolute_pct (API table)"
                  testId="rna-codon-bar"
                />
              </InfographicFrame>
              {codonHeat.x.length > 0 ? (
                <InfographicFrame
                  title="Synonymous codon percents"
                  dataSource="HelixScope Core codon_usage"
                  method="Synonymous_pct from the API codon table"
                  units="percent within amino acid"
                  status={result.status}
                  shows="How synonymous codons are distributed in this table"
                  doesNotShow="Selection, tRNA abundance, or species identity"
                >
                  <ScienceHeatmap
                    z={codonHeat.z}
                    x={codonHeat.x}
                    y={codonHeat.y}
                    title="Codon Synonymous_pct (API table; not recomputed)"
                    testId="rna-codon-heatmap"
                  />
                </InfographicFrame>
              ) : null}
              <DataTable rows={codonRows} caption="Codon usage" testId="rna-codon-usage" />
              <button type="button" className="btn-secondary" onClick={() => downloadCsv("codon_usage.csv", codonRows)}>
                Download codon usage (CSV)
              </button>
            </section>
            </ResultAnchor>
          ) : null}
          {rscuRows.length > 0 ? (
            <ResultAnchor id="rna-rscu">
            <section className="hs-panel science-surface">
              <h2>RSCU</h2>
              <InfographicFrame
                title="RSCU"
                dataSource="HelixScope Core rscu"
                method="Relative synonymous codon usage on the codon table"
                units="RSCU"
                status={result.status}
                shows="Which synonymous codons are used more or less than equal usage (1.0)"
                doesNotShow="Proof of a biological CDS or translational efficiency"
              >
                <ScienceBarChart
                  labels={rscuRows.map((row) => row.codon ?? row.Codon)}
                  values={rscuRows.map((row) => row.RSCU ?? row.rscu ?? row.value)}
                  title="RSCU (API table)"
                  testId="rna-rscu-bar"
                />
              </InfographicFrame>
              <DataTable rows={rscuRows} caption="RSCU" emptyLabel="0 codon rows." />
            </section>
            </ResultAnchor>
          ) : null}
          {kmerTop.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>k-mers</h2>
              <ScienceBarChart
                labels={kmerTop.map((row) => row.kmer)}
                values={kmerTop.map((row) => row.count)}
                title="RNA k-mers (API ranking; U treated as T)"
                testId="rna-kmer-bar"
              />
            </section>
          ) : null}
          <SequenceViewer sequence={String(result.sequence || "")} hash={String(result.sequence_hash || "")} testId="rna-sequence-viewer" />
          <section className="hs-panel science-surface hs-viz-hero">
            <h2>RNA 3D</h2>
            <ScientificWarning>
              Illustrative A-RNA is not experimental. ViennaRNA MFE does not yield 3D coordinates. 1RNA is
              a deposited experimental RNA structure.
            </ScientificWarning>
            <div className="hs-numeric-row">
              <div className="field">
                <label htmlFor="rna-3d-lod">Visual LOD</label>
                <select id="rna-3d-lod" value={helixLod} onChange={(e) => setHelixLod(e.target.value)}>
                  <option value="high">high</option>
                  <option value="medium">medium</option>
                  <option value="low">low</option>
                </select>
              </div>
              <div className="field">
                <label htmlFor="rna-3d-start">RNA helix start (0-based)</label>
                <input
                  id="rna-3d-start"
                  type="number"
                  min={0}
                  max={Math.max(0, seqLen - 1)}
                  value={helixStart}
                  onChange={(e) => setHelixStart(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="rna-3d-end">RNA helix end (exclusive)</label>
                <input
                  id="rna-3d-end"
                  type="number"
                  min={1}
                  max={Math.max(1, seqLen)}
                  value={helixEnd}
                  onChange={(e) => setHelixEnd(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="rna-3d-pos">Inspect RNA 3D position (0-based)</label>
                <input
                  id="rna-3d-pos"
                  type="number"
                  min={0}
                  max={Math.max(0, seqLen - 1)}
                  value={helixInspect}
                  onChange={(e) => setHelixInspect(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="rna-3d-prefer">RNA 3D path</label>
                <select id="rna-3d-prefer" value={helixPrefer} onChange={(e) => setHelixPrefer(e.target.value)}>
                  <option value="auto">auto</option>
                  <option value="illustrative">illustrative</option>
                  <option value="experimental">experimental</option>
                  <option value="predicted">predicted</option>
                </select>
              </div>
            </div>
            <label>
              <input type="checkbox" checked={helixCenter} onChange={(e) => setHelixCenter(e.target.checked)} /> Center
              helix on selected base
            </label>
            <div className="hs-actions">
              <button
                type="button"
                className="btn-secondary"
                onClick={async () => {
                  setError(null);
                  try {
                    const { data } = await helixApi.rnaHelix3d({
                      sequence: String(result.sequence || draft.sequence),
                      prefer: helixPrefer,
                      lod: helixLod,
                      region_start: helixStart,
                      region_end: helixEnd,
                      inspect_position: helixInspect,
                      center_on_selected: helixCenter,
                    });
                    setHelix(asRecord(data.result));
                  } catch (err) {
                    setError(err);
                  }
                }}
              >
                Load illustrative A-RNA helix
              </button>
              <button
                type="button"
                className="btn-secondary"
                onClick={() => {
                  patchDraft("viewer", { structureId: "1RNA" });
                  router.push("/structure/viewer");
                }}
              >
                Open experimental 1RNA
              </button>
            </div>
            {helix ? (
              <>
                <p data-testid="rna-3d-kind">
                  kind {String(helix.kind || "")} · status {String(helix.status || "")}
                  {helix.partial_view ? " · PARTIAL_VIEW" : ""}
                </p>
                <p data-testid="rna-3d-mapping">
                  Sequence position {String(helix.inspect_position)} base {String(helix.inspect_base ?? "N/A")} mapping{" "}
                  {String(asRecord(helix.position_mapping).status || "")}
                </p>
                <p className="hs-topbar-meta">{String(helix.disclaimer || helix.kind_label || helix.reason || "")}</p>
                {helixViewer.coordinates || helixViewer.atoms ? (
                  <StructureViewportPlotly
                    scene={{ ...helixViewer, kind: helix.kind, status: helix.status, structure_id: "illustrative-arna" }}
                    orbit={orbit}
                    onOrbit={setOrbit}
                  />
                ) : null}
              </>
            ) : null}
          </section>
          <ResultAnchor id="rna-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="rna-raw">
          <RawDataPanel
            result={result}
            fasta={fastaRecord(String(result.sequence_hash || "rna"), String(result.sequence || draft.sequence))}
            fastaName="helixscope-rna.fasta"
          />
          </ResultAnchor>
          <ResultAnchor id="rna-provenance">
          <ProvenancePanel envelope={envelope} />
          </ResultAnchor>
        </>
      ) : null}
    </div>
  );
}
