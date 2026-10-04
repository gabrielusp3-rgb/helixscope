"use client";

import dynamic from "next/dynamic";
import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { pollJob } from "@/lib/api/jobs";
import { asList, asRecord, asRecordRows, chartNumeric, formatScientificValue } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import {
  EmptyState,
  ExternalDisclosure,
  JobState,
  MetricGrid,
  ScientificError,
  ScientificStatus,
  ScientificWarning,
  SequenceInput,
} from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput, FastaFileField, SequenceViewer } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { CasSystemCard, offTargetScopeLabel } from "@/components/science/CasSystemCard";
import { CrisprLocusViewer } from "@/components/tracks/CrisprLocusViewer";
import { ScienceBarChart } from "@/components/charts/ScienceBarChart";
import { ScienceFeatureMap } from "@/components/charts/ScienceFeatureMap";
import { ScienceHistogram } from "@/components/charts/ScienceHistogram";
import { ScienceScatter } from "@/components/charts/ScienceScatter";
import { downloadCsv, downloadJson } from "@/lib/export/download";
import { countByLabel } from "@/lib/charts/countBy";
import type { Orbit } from "@/lib/gesture/math";

const StructureViewportPlotly = dynamic(
  () => import("@/components/structure/StructureViewportPlotly").then((m) => m.StructureViewportPlotly),
  { ssr: false },
);

const DEFAULT_ORBIT: Orbit = { yaw: 0.9, pitch: 0.25, radius: 2.8 };

export function CrisprWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.crispr;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [jobExec, setJobExec] = useState("");
  const [models, setModels] = useState<Record<string, unknown>>({});
  const [casOffinder, setCasOffinder] = useState<Record<string, unknown>>({});
  const [catalog, setCatalog] = useState<Record<string, unknown>>({});
  const [complex, setComplex] = useState<Record<string, unknown> | null>(null);
  const [orbit, setOrbit] = useState<Orbit>(DEFAULT_ORBIT);
  const [complexId, setComplexId] = useState<"4UN3" | "4OO8">("4UN3");
  const [selectedGuide, setSelectedGuide] = useState(0);
  const [systems, setSystems] = useState<Array<Record<string, unknown>>>([
    { canonical_key: "SpCas9", display_name: "SpCas9", target_molecule: "DNA" },
  ]);
  const selectedSystem = systems.find((row) => String(row.canonical_key) === draft.cas) ?? systems[0];
  const targetMolecule = String(selectedSystem?.target_molecule || "DNA");
  const envelope = results.crispr;
  const result = asRecord(envelope?.result);
  const guides = asList(result.guides).map((row) => asRecord(row));
  const availability = asRecord(result.model_availability);
  const targetLength = Number(result.target_length ?? 0);
  const azimuth = asRecord(models.on_target_azimuth ?? asRecord(availability.on_target_azimuth));
  const advancedAvailable = azimuth.available === true;

  useEffect(() => {
    void helixApi.crisprAvailability().then(({ data }) => setModels(asRecord(data.result))).catch(setError);
    void helixApi.casOffinderAvailability().then(({ data }) => setCasOffinder(asRecord(data.result))).catch(() => null);
    void helixApi.structureCatalog().then(({ data }) => setCatalog(asRecord(data.result))).catch(() => null);
    void helixApi
      .crisprSystems()
      .then(({ data }) => {
        const rows = asList(asRecord(data.result).systems).map((row) => asRecord(row));
        if (rows.length > 0) setSystems(rows);
      })
      .catch(setError);
  }, []);

  async function run() {
    const gate = next();
    setError(null);
    setResult("crispr", undefined);
    try {
      const { data } = await helixApi.crisprGuides(
        {
          sequence: draft.sequence,
          cas_system: draft.cas,
          include_explanation: true,
          run_off_target: draft.offTarget,
          max_mismatches: 3,
        },
        gate.signal,
      );
      if (gate.isCurrent()) {
        setResult("crispr", data);
        setSelectedGuide(0);
      }
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  const first = guides[selectedGuide] ?? guides[0] ?? {};
  const otSummary = asRecord(first.off_target_summary);
  const otRows = asRecordRows(otSummary, "mismatches");
  const otHits = asList(first.off_targets).map((row) => asRecord(row));
  const pamCounts = countByLabel(otHits, "pam_class");
  const contigCounts = countByLabel(otHits, "chromosome_or_contig");
  const hsuScores = otHits.map((hit) => hit.hsu_hit_score).filter((value) => value != null);
  const poly = asRecord(first.poly_t);
  const hairpin = asRecord(first.self_complementarity);
  const primers = asRecord(result.validation_primers);
  const guideGc = chartNumeric(first.gc_content);
  const complexViewer = asRecord(complex?.viewer);
  const complexes = asList(catalog.complexes).map((row) => asRecord(row));
  const scopeLabel = offTargetScopeLabel(result.off_target_scope, result.not_genome_wide);
  const targetSequence = String(result.target_sequence || draft.sequence);
  const pamLabel = String(asRecord(result.system).pam_display || selectedSystem?.pam_display || "PAM/PFS from registry");
  const cutLabel =
    String(asRecord(result.system).mode || selectedSystem?.mode) === "rna_nuclease"
      ? "No DNA cut site. Cas13 is RNA-targeting."
      : String(asRecord(result.system).editor || selectedSystem?.editor) === "base_editor"
        ? "Edit window from API, not a DSB cut."
        : "Expected cut/edit position from the API when defined.";
  const cutTracks = [
    {
      name: "Guides",
      features: guides
        .filter((guide) => guide.start_0based != null && guide.end_0based != null)
        .map((guide) => ({
          start: guide.start_0based,
          end: guide.end_0based,
          label: String(guide.guide_sequence || "guide"),
        })),
    },
  ];

  return (
    <div className="hs-page" data-testid="hs-page">
      <p className="kicker">Discovery</p>
      <h1>CRISPR Design</h1>
      <ScientificWarning>
        Cas-system identity comes from GET /api/v1/crispr/systems. Informal names such as Cas12a, CBE,
        ABE, or Prime Editor are invalid. Cas13 requires RNA. DNA-target nucleases require DNA including
        PAM. Cas-OFFinder on HELIXSCOPE_TEST_REF is TEST REFERENCE scoped. It is not genome-wide.
      </ScientificWarning>
      <p data-testid="crispr-models">
        Advanced on-target models (Azimuth / Rule Set 2 / DeepHF):{" "}
        {advancedAvailable ? "Available in this installation" : "Unavailable in this installation"}
        . Cas-OFFinder {formatScientificValue(casOffinder.available ?? casOffinder.status)}.
      </p>
      <CompactInput compact={guides.length > 0} summary="Edit CRISPR input">
        <div className="field">
          <label htmlFor="seq-name">Sequence name</label>
          <input id="seq-name" value={draft.sequenceName} onChange={(e) => patchDraft("crispr", { sequenceName: e.target.value })} data-testid="crispr-name" />
        </div>
        <SequenceInput
          id="crispr-dna"
          label={targetMolecule === "RNA" ? "Target RNA" : "Target DNA"}
          hint={
            targetMolecule === "RNA"
              ? "Cas13 requires RNA (A/U/C/G) including PFS context. DNA is not transcribed automatically."
              : "Paste the genomic or plasmid locus to scan, including PAM. Single-record FASTA is accepted."
          }
          value={draft.sequence}
          onChange={(sequence) => patchDraft("crispr", { sequence })}
          disabled={scienceLocked}
        />
        <FastaFileField id="crispr-file" disabled={scienceLocked} onText={(sequence) => patchDraft("crispr", { sequence })} />
        <div className="field">
          <label htmlFor="cas">Cas system</label>
          <select
            id="cas"
            data-testid="crispr-cas"
            value={draft.cas}
            onChange={(e) => patchDraft("crispr", { cas: e.target.value })}
          >
            {systems.map((row) => (
              <option key={String(row.canonical_key)} value={String(row.canonical_key)}>
                {String(row.display_name || row.canonical_key)}
              </option>
            ))}
          </select>
        </div>
        <label>
          <input
            type="checkbox"
            checked={draft.offTarget}
            onChange={(e) => patchDraft("crispr", { offTarget: e.target.checked })}
          />{" "}
          Sequence-local off-target scan (not genome-wide)
        </label>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void run()}>
          Design guides
        </button>
      </CompactInput>
      <CasSystemCard system={selectedSystem} />
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? <EmptyState title="No CRISPR result. Target DNA starts empty." /> : (
        <>
          <AnalysisNav
            items={[
              { id: "crispr-summary", label: "Overview" },
              { id: "crispr-locus", label: "Locus", available: guides.length > 0 },
              { id: "crispr-guides-section", label: "Guides", available: guides.length > 0 },
              { id: "crispr-inspector", label: "Guide inspector", available: guides.length > 0 },
              { id: "crispr-ot", label: "Off-target", available: otRows.length > 0 || Boolean(result.off_target_scope) },
              { id: "crispr-meaning", label: "Interpretation" },
              { id: "crispr-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="crispr-summary">
          <MetricGrid
            metrics={[
              { id: "cas", label: "Cas", value: result.cas_system },
              { id: "n", label: "Guides", value: result.n_guides },
              { id: "tlen", label: "Target length", value: result.target_length, unit: "nt" },
              { id: "scope", label: "Off-target scope", value: result.off_target_scope },
              { id: "gw", label: "Not genome-wide", value: result.not_genome_wide },
            ]}
          />
          <p className="hs-topbar-meta">{String(result.input_note || result.spacer_is_not_target || "")}</p>
          </ResultAnchor>
          {guides.length === 0 ? (
            <EmptyState title="No PAM sites in this target">
              {targetLength > 0 && targetLength < 23
                ? "A 20-nt spacer is not a valid substitute for Target DNA. Paste the locus including PAM (NGG for SpCas9), typically at least 23 nt."
                : "No guide candidates were found on the provided sequence. This is a CRISPR search result with 0 guides, not a molecule-classification error."}
            </EmptyState>
          ) : null}
          <p data-testid="crispr-scope">
            scope {String(result.off_target_scope)} · visible label {scopeLabel}. TEST REFERENCE labeling is
            required whenever the search is not a public assembly. This is not genome-wide.
          </p>
          {guides.length > 0 ? (
            <ResultAnchor id="crispr-locus">
              <CrisprLocusViewer
                sequence={targetSequence}
                guides={guides}
                selectedIndex={Math.min(selectedGuide, guides.length - 1)}
                onSelect={setSelectedGuide}
                pamLabel={pamLabel}
                cutLabel={cutLabel}
              />
            </ResultAnchor>
          ) : null}
          {guides[0] ? (
            <ResultAnchor id="crispr-inspector">
            <section className="hs-panel science-surface" data-testid="crispr-guide-inspector">
              <h2>Guide inspector</h2>
              <MetricGrid
                metrics={[
                  { id: "guide", label: "Guide", value: first.guide_sequence },
                  { id: "pam", label: "PAM/PFS", value: first.pam_sequence ?? first.pfs },
                  { id: "strand", label: "Strand", value: first.strand },
                  { id: "start", label: "Start 0-based", value: first.start_0based ?? first.start },
                  { id: "end", label: "End 0-based", value: first.end_0based ?? first.end },
                  { id: "cut", label: "Cut/edit", value: first.cut_site ?? first.nick_site ?? first.edit_window_start },
                  { id: "cfd", label: "CFD", value: asRecord(availability.specificity_cfd).available === false ? null : first.cfd_score },
                  { id: "mit", label: "MIT", value: asRecord(availability.specificity_mit).available === false ? null : first.mit_score },
                  { id: "eff", label: "Efficiency (heuristic)", value: first.efficiency_score },
                  { id: "risk", label: "Local risk", value: first.risk },
                  { id: "gc-guide", label: "GC", value: first.gc_content, unit: "%" },
                ]}
              />
              <p className="hs-topbar-meta">
                Off-target scope for this result: {scopeLabel}. Mapping to deposited Cas complexes remains
                UNCERTAIN unless an exact polymer match exists.
              </p>
            </section>
            </ResultAnchor>
          ) : null}
          {guides.length > 0 ? (
            <ResultAnchor id="crispr-guides-section">
            <section className="hs-panel science-surface">
              <h2>Guide ranking</h2>
              <InfographicFrame
                title="Heuristic efficiency versus position"
                dataSource="HelixScope Core CRISPR guide table"
                method="Local heuristic inspired by Rule Set 2 positional terms; not the published model"
                units="heuristic score; coordinates"
                status="HEURISTIC"
                shows="Where candidate guides sit and how the local heuristic ranked them"
                doesNotShow="Azimuth, DeepHF, or genome-wide activity"
              >
                <ScienceScatter
                  x={guides.map((guide) => guide.position ?? guide.start_0based)}
                  y={guides.map((guide) => guide.efficiency_score)}
                  text={guides.map((guide) => guide.guide_sequence)}
                  xTitle="position"
                  yTitle="heuristic efficiency"
                  title="Guide heuristic efficiency versus position (API scores)"
                  testId="crispr-guide-scatter"
                />
              </InfographicFrame>
              {targetLength > 0 && cutTracks[0].features.length > 0 ? (
                <InfographicFrame
                  title="Guide spans on the target"
                  dataSource="HelixScope Core CRISPR coordinates"
                  method="PAM search on the pasted DNA locus"
                  units="0-based coordinates"
                  status={result.status}
                  shows="Guide and PAM intervals on the submitted target DNA"
                  doesNotShow="A Cas-DNA complex or experimental cleavage"
                >
                  <ScienceFeatureMap
                    length={targetLength}
                    tracks={cutTracks}
                    title="Guide spans on the target (API coordinates)"
                    testId="crispr-cut-map"
                  />
                </InfographicFrame>
              ) : null}
            </section>
            </ResultAnchor>
          ) : null}
          {otRows.length > 0 ? (
            <ResultAnchor id="crispr-ot">
            <section className="hs-panel science-surface">
              <h2>Sequence-local off-target (top guide)</h2>
              <InfographicFrame
                title="Local mismatch histogram"
                dataSource="HelixScope Core sequence-local off-target summary"
                method="Scan of the pasted sequence only"
                units="site counts by mismatch number"
                status={result.off_target_scope}
                shows="How many local candidate sites have 0-N mismatches"
                doesNotShow="Genome-wide off-targets, CFD, or MIT Specificity Score"
              >
                <ScienceBarChart
                  labels={otRows.map((row) => row.mismatches)}
                  values={otRows.map((row) => row.value ?? row.count)}
                  xTitle="mismatches"
                  yTitle="sites"
                  title="Local off-target mismatch histogram (API summary)"
                  testId="crispr-ot-hist"
                />
              </InfographicFrame>
              <p className="hs-topbar-meta">
                Sequence-local scan of the pasted target. Not CFD, not MIT, not genome-wide.
                Specificity proxy {formatScientificValue(first.specificity_proxy)}.
              </p>
              {contigCounts.labels.filter((label) => label && label !== "unnamed").length > 1 ? (
                <ScienceBarChart
                  labels={contigCounts.labels}
                  values={contigCounts.values}
                  xTitle="Contig in provided reference"
                  yTitle="Verified hits"
                  title="Off-target hits by contig (API hit labels)"
                  testId="crispr-ot-contig"
                />
              ) : null}
              {pamCounts.labels.length > 0 ? (
                <ScienceBarChart
                  labels={pamCounts.labels}
                  values={pamCounts.values}
                  xTitle="PAM class in verified hits"
                  yTitle="Verified hits"
                  title="Off-target PAM classes (API hit labels)"
                  testId="crispr-ot-pam"
                />
              ) : null}
              {hsuScores.length > 0 ? (
                <ScienceHistogram
                  values={hsuScores}
                  xTitle="Hsu hit score (API; not MIT Specificity Score)"
                  yTitle="Verified hits"
                  title="Hsu single-hit scores on verified local sites"
                  testId="crispr-ot-hsu"
                />
              ) : null}
              <DataTable
                rows={otHits}
                caption="Local off-target sites (pasted sequence only)"
                testId="crispr-ot-hits"
                emptyLabel="0 verified hits under selected method/reference."
              />
            </section>
            </ResultAnchor>
          ) : null}
          <DataTable
            rows={guides.map((g, index) => ({
              index,
              guide_sequence: g.guide_sequence,
              pam_sequence: g.pam_sequence,
              position: g.position,
              strand: g.strand,
              start_0based: g.start_0based,
              end_0based: g.end_0based,
              cut_site: g.cut_site,
              gc: g.gc_content ?? g.gc_percent,
              cfd: g.cfd_score,
              mit: g.mit_score,
              ruleset2: g.ruleset2_score ?? g.doench_score,
              deephf: g.deephf_score,
              efficiency: g.efficiency_score,
              risk: g.risk,
              specificity_proxy: g.specificity_proxy,
            }))}
            testId="crispr-guides"
            emptyLabel="0 guides."
            selectedKey={String(selectedGuide)}
            onRowSelect={(row) => setSelectedGuide(Number(row.index) || 0)}
          />
          {guides.length > 0 ? (
            <button type="button" className="btn-secondary" onClick={() => downloadCsv("crispr_guides.csv", guides)}>
              Download guides (CSV)
            </button>
          ) : null}
          <p className="hs-topbar-meta">
            Azimuth, DeepHF, CFD and MIT unavailable values stay N/A, never 0. Efficiency method:{" "}
            {String(first.efficiency_method || "")}. This is a local heuristic inspired by Rule Set 2
            positional terms, not the published model.
          </p>
          {guides[0] ? (
            <section className="hs-panel science-surface">
              <h2>Guide sequence viewer</h2>
              <SequenceViewer
                sequence={`${String(first.guide_sequence || "")}${String(first.pam_sequence || "")}`}
                testId="crispr-guide-seq"
              />
              <p className="hs-topbar-meta">
                Cut site {formatScientificValue(first.cut_site)} · GC {formatScientificValue(first.gc_content)}% ·
                hairpin stem {formatScientificValue(hairpin.max_stem)} · poly-T {formatScientificValue(poly.max_run)}
              </p>
              {guideGc !== null && (guideGc < 40 || guideGc > 80) ? (
                <ScientificWarning>
                  GC of this spacer is outside the 40-80% range usually recommended for reliable activity.
                </ScientificWarning>
              ) : null}
              {poly.has_signal ? (
                <ScientificWarning>
                  Poly-T run may terminate Pol III transcription and truncate the guide.
                </ScientificWarning>
              ) : null}
              {hairpin.risk === "high" ? (
                <ScientificWarning>
                  Self-complementary stem may form a hairpin that competes with sgRNA scaffold folding.
                </ScientificWarning>
              ) : null}
            </section>
          ) : null}
          {primers.forward_primer || primers.status === "UNAVAILABLE" ? (
            <section className="hs-panel science-surface" data-testid="crispr-primers">
              <h2>Primer design for validation</h2>
              {primers.status === "UNAVAILABLE" ? (
                <p>{String(primers.reason || "Primers unavailable for this locus.")}</p>
              ) : (
                <>
                  <p>
                    Forward Tm {formatScientificValue(primers.forward_tm)} C · Reverse Tm{" "}
                    {formatScientificValue(primers.reverse_tm)} C · amplicon{" "}
                    {formatScientificValue(primers.amplicon_size)} bp
                  </p>
                  <p className="mono">F {String(primers.forward_primer || "")}</p>
                  <p className="mono">R {String(primers.reverse_primer || "")}</p>
                  <p className="hs-topbar-meta">
                    These primers are designed for validation PCR around the predicted cut, not as a cloning kit.
                  </p>
                </>
              )}
            </section>
          ) : null}
          <section className="hs-panel science-surface">
            <h2>Cas-OFFinder job</h2>
            <button
              type="button"
              className="btn-secondary"
              onClick={async () => {
                if (!guides[0]) return;
                setError(null);
                try {
                  const { data } = await helixApi.casOffinderSubmit({
                    guide_sequence: String(guides[0].guide_sequence),
                    assembly_id: "HELIXSCOPE_TEST_REF",
                    pam: String(guides[0].pam_sequence || "NGG"),
                    cas_system: String(result.cas_system || draft.cas),
                  });
                  recordTransfer({
                    kind: "crispr-to-casoffinder",
                    note: "Guide submitted to Cas-OFFinder TEST REFERENCE, not GRCh38.",
                    payload: { assembly_id: "HELIXSCOPE_TEST_REF", cas_system: String(result.cas_system || draft.cas) },
                  });
                  setJobExec(data.execution_status);
                  const done = await pollJob(data);
                  setJobExec(String(done.execution_status || "COMPLETED"));
                } catch (err) {
                  setError(err);
                }
              }}
            >
              Run Cas-OFFinder on TEST REFERENCE
            </button>
            {jobExec ? <JobState execution={jobExec} /> : null}
          </section>
          <section className="hs-panel science-surface hs-viz-hero">
            <h2>CRISPR 3D (deposited complexes)</h2>
            <ScientificWarning>
              HelixScope does not fabricate a Cas protein or a DNA-Cas complex from a spacer. 4UN3 and
              4OO8 are experimental RCSB depositions. Guide/PAM highlighting remains UNCERTAIN unless an
              exact polymer subsequence match exists.
            </ScientificWarning>
            <ExternalDisclosure service="RCSB Files API">
              Fetch is explicit. Coordinates are experimental mmCIF, not a predicted Cas complex.
            </ExternalDisclosure>
            <p className="hs-topbar-meta">
              Catalog mapping contract: {formatScientificValue(asRecord(catalog.mapping_contract).structure_mapping ?? catalog.status)}.
            </p>
            <DataTable rows={complexes} caption="Experimental Cas-guide-DNA pointers" testId="crispr-complex-catalog" />
            <div className="field">
              <label htmlFor="complex-id">Catalog PDB</label>
              <select id="complex-id" value={complexId} onChange={(e) => setComplexId(e.target.value as "4UN3" | "4OO8")}>
                <option value="4UN3">4UN3</option>
                <option value="4OO8">4OO8</option>
              </select>
            </div>
            <button
              type="button"
              className="btn-secondary"
              data-testid="crispr-fetch-4un3"
              onClick={async () => {
                setError(null);
                try {
                  const { data } = await helixApi.structureExperimentalComplex({ structure_id: complexId });
                  setComplex(asRecord(data.result));
                } catch (err) {
                  setError(err);
                }
              }}
            >
              Fetch deposited Cas complex from RCSB
            </button>
            {complex ? (
              <>
                <p data-testid="crispr-3d-status">
                  {String(complex.structure_id)} · <ScientificStatus status={complex.status} /> · mapping{" "}
                  <ScientificStatus status={complex.mapping_status} />
                </p>
                <p className="hs-topbar-meta">{String(complex.mapping_note || complex.disclaimer || "")}</p>
                {complexViewer.coordinates || complexViewer.atoms ? (
                  <StructureViewportPlotly
                    scene={{ ...complexViewer, kind: complex.kind, status: complex.status, structure_id: complex.structure_id }}
                    orbit={orbit}
                    onOrbit={setOrbit}
                  />
                ) : null}
              </>
            ) : null}
          </section>
          <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-crispr.json", result)}>
            Download CRISPR JSON
          </button>
          <ResultAnchor id="crispr-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="crispr-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} />
        </>
      )}
    </div>
  );
}
