"use client";

import dynamic from "next/dynamic";
import { useRouter } from "next/navigation";
import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, asRecordRows, chartNumeric } from "@/lib/api/numeric";
import { dinucleotideOeMatrix } from "@/lib/charts/dinucleotide";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import {
  EmptyState,
  InvalidMoleculeResult,
  MetricGrid,
  ScientificError,
  SequenceInput,
} from "@/components/science/SciencePrimitives";
import { InterpretationPanel, MetricExplanations, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput, FastaFileField, SequenceViewer } from "@/components/science/ResultChrome";
import { AnalysisNav, AnonymousSequenceNotice, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { GenomicTrack } from "@/components/tracks/GenomicTrack";
import { BeginnerOnly } from "@/lib/visual/presentation";
import { ScienceChart } from "@/components/charts/ScienceChart";
import { ScienceBarChart } from "@/components/charts/ScienceBarChart";
import { ScienceHeatmap } from "@/components/charts/ScienceHeatmap";
import { ScienceHistogram } from "@/components/charts/ScienceHistogram";
import { ScienceFeatureMap } from "@/components/charts/ScienceFeatureMap";
import { downloadCsv, downloadJson, downloadText, fastaRecord } from "@/lib/export/download";
import type { Orbit } from "@/lib/gesture/math";

const StructureViewportPlotly = dynamic(
  () => import("@/components/structure/StructureViewportPlotly").then((m) => m.StructureViewportPlotly),
  { ssr: false },
);

const DEFAULT_ORBIT: Orbit = { yaw: 0.7, pitch: 0.35, radius: 2.4 };

export function DnaWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.dna;
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
  const [selectedFeature, setSelectedFeature] = useState<{ start?: number; end?: number; kind?: string } | null>(null);
  const envelope = results.dna;
  const result = asRecord(envelope?.result);
  const status = String(result.status || "");
  const invalid = status === "ERROR";
  const computed = status === "COMPUTED";

  async function analyze() {
    const gate = next();
    setLoading(true);
    setError(null);
    setResult("dna", undefined);
    setHelix(null);
    try {
      const { data } = await helixApi.dnaAnalyze(
        {
          sequence: draft.sequence,
          include_orfs: draft.includeOrfs,
          include_profiles: draft.includeProfiles,
          include_restriction: draft.includeRestriction,
          include_cpg: draft.includeCpg,
          include_kmers: draft.includeKmers,
          include_santalucia: draft.includeSantalucia,
          include_explanation: true,
          min_orf_length: draft.minOrfLength,
          drop_overlapping_orfs: draft.dropNested,
          profile_window: draft.profileWindow,
          profile_step: draft.profileStep,
          kmer_k: draft.kmerK,
          kmer_top_n: 30,
          santalucia_na_m: Number(draft.santaluciaNa) || 0.05,
          santalucia_oligo_nm: Number(draft.santaluciaOligo) || 250,
        },
        gate.signal,
      );
      if (!gate.isCurrent()) return;
      setResult("dna", data);
      const nextLen = Number(asRecord(data.result).length) || 0;
      setHelixStart(0);
      setHelixEnd(nextLen || 1);
      setHelixInspect(0);
    } catch (err) {
      if (!gate.isCurrent()) return;
      setError(err);
    } finally {
      if (gate.isCurrent()) setLoading(false);
    }
  }

  const tm = asRecord(result.tm_report);
  const profiles = asList(result.windowed_profiles).map((row) => asRecord(row));
  const compositionRows = asRecordRows(result.composition, "symbol").filter((row) => Number(row.count) > 0);
  const dinucRows = asRecordRows(result.dinucleotides, "dinucleotide");
  const dinucHeat = dinucleotideOeMatrix(result.dinucleotides);
  const santalucia = asRecord(result.santalucia_tm);
  const kmers = asRecord(result.kmer_summary);
  const kmerTop = asList(kmers.top).map((row) => asRecord(row));
  const orfs = asList(result.orfs).map((row) => asRecord(row));
  const cpg = asList(result.cpg_islands).map((row) => asRecord(row));
  const enzymeRows = asList(result.restriction_enzyme_rows).map((row) => asRecord(row));
  const restrictionHits = asList(result.restriction_hits).map((row) => asRecord(row));
  const tracks = asList(result.sequence_feature_tracks);
  const x = profiles.map((row, index) => chartNumeric(row.midpoint ?? row.center ?? row.position) ?? index);
  const helixViewer = asRecord(helix?.viewer);
  const busy = loading || scienceLocked;
  const seqLen = Number(result.length) || 0;
  const cpgOe = asRecord(asRecord(result.dinucleotides).CG).observed_expected;
  const compositionBar = compositionRows.filter((row) => ["A", "T", "C", "G"].includes(String(row.symbol)));

  return (
    <div className="hs-page" data-testid="hs-page">
      <p className="kicker">Analysis</p>
      <h1>DNA</h1>
      <p className="hs-lede">
        GC, Tm, molecular weight, entropy, ORF, restriction, dinucleotide O/E and profile values are
        computed by HelixScope API v1. This page only renders the response. Single-record FASTA is
        parsed the same way as the validated legacy reference; multi-FASTA is refused.
      </p>
      <CompactInput compact={computed} summary="Edit DNA input">
        <h2>Input</h2>
        <SequenceInput
          id="dna-sequence"
          label="DNA sequence"
          hint="Paste a canonical DNA string (A, T, C, G) or a single-record FASTA. Line wraps are removed. Numbers and punctuation in a non-FASTA paste remain invalid."
          value={draft.sequence}
          onChange={(sequence) => patchDraft("dna", { sequence })}
          disabled={busy}
        />
        <FastaFileField id="dna-file" disabled={busy} onText={(sequence) => patchDraft("dna", { sequence })} />
        <div className="hs-option-row">
          <label>
            <input
              type="checkbox"
              checked={draft.includeOrfs}
              onChange={(e) => patchDraft("dna", { includeOrfs: e.target.checked })}
            />
            Include ORFs
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.dropNested}
              onChange={(e) => patchDraft("dna", { dropNested: e.target.checked })}
            />
            Discard nested/overlapping ORFs
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.includeProfiles}
              onChange={(e) => patchDraft("dna", { includeProfiles: e.target.checked })}
            />
            Include windowed profiles
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.includeRestriction}
              onChange={(e) => patchDraft("dna", { includeRestriction: e.target.checked })}
            />
            Include restriction sites
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.includeCpg}
              onChange={(e) => patchDraft("dna", { includeCpg: e.target.checked })}
            />
            Include CpG islands
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.includeKmers}
              onChange={(e) => patchDraft("dna", { includeKmers: e.target.checked })}
            />
            Include k-mers
          </label>
          <label>
            <input
              type="checkbox"
              checked={draft.includeSantalucia}
              onChange={(e) => patchDraft("dna", { includeSantalucia: e.target.checked })}
            />
            Include SantaLucia Tm
          </label>
        </div>
        <div className="hs-numeric-row">
          <div className="field">
            <label htmlFor="dna-orf-min">Minimum ORF length (nt)</label>
            <input
              id="dna-orf-min"
              type="number"
              min={30}
              max={3000}
              value={draft.minOrfLength}
              onChange={(e) => patchDraft("dna", { minOrfLength: Number(e.target.value) || 150 })}
            />
          </div>
          <div className="field">
            <label htmlFor="dna-window">Profile window (bp)</label>
            <input
              id="dna-window"
              type="number"
              min={10}
              max={2000}
              value={draft.profileWindow}
              onChange={(e) => patchDraft("dna", { profileWindow: Number(e.target.value) || 100 })}
            />
          </div>
          <div className="field">
            <label htmlFor="dna-step">Profile step (bp)</label>
            <input
              id="dna-step"
              type="number"
              min={1}
              max={2000}
              value={draft.profileStep}
              onChange={(e) => patchDraft("dna", { profileStep: Number(e.target.value) || 50 })}
            />
          </div>
          <div className="field">
            <label htmlFor="dna-k">k-mer length</label>
            <input
              id="dna-k"
              type="number"
              min={1}
              max={5}
              value={draft.kmerK}
              onChange={(e) => patchDraft("dna", { kmerK: Number(e.target.value) || 3 })}
            />
          </div>
          <div className="field">
            <label htmlFor="dna-na">SantaLucia Na+ (M)</label>
            <input
              id="dna-na"
              value={draft.santaluciaNa}
              onChange={(e) => patchDraft("dna", { santaluciaNa: e.target.value })}
            />
          </div>
          <div className="field">
            <label htmlFor="dna-oligo">SantaLucia oligo (nM)</label>
            <input
              id="dna-oligo"
              value={draft.santaluciaOligo}
              onChange={(e) => patchDraft("dna", { santaluciaOligo: e.target.value })}
            />
          </div>
        </div>
        <div className="hs-actions">
          <button type="button" className="btn-primary" onClick={() => void analyze()} disabled={busy}>
            Analyze
          </button>
          <button type="button" className="btn-ghost" onClick={() => patchDraft("dna", { sequence: "" })}>
            Clear
          </button>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              recordTransfer({
                kind: "dna-to-crispr",
                note: "Explicit DNA to CRISPR sequence copy",
                payload: { sequence: draft.sequence },
              });
              patchDraft("crispr", { sequence: draft.sequence });
              router.push("/discovery/crispr");
            }}
          >
            Send sequence to CRISPR
          </button>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              recordTransfer({
                kind: "dna-to-blast",
                note: "Explicit DNA to BLAST query copy",
                payload: { sequence: draft.sequence },
              });
              patchDraft("blast", { query: draft.sequence });
              router.push("/data/blast");
            }}
          >
            Send sequence to BLAST
          </button>
          {computed ? (
            <>
              <button
                type="button"
                className="btn-secondary"
                onClick={() =>
                  downloadText(
                    "helixscope-dna.fasta",
                    fastaRecord(String(result.sequence_hash || "dna"), String(result.sequence || draft.sequence)),
                  )
                }
              >
                Download FASTA
              </button>
              <button
                type="button"
                className="btn-secondary"
                onClick={() => downloadJson("helixscope-dna.json", result)}
              >
                Download JSON
              </button>
            </>
          ) : null}
        </div>
      </CompactInput>
      {loading ? <p data-testid="dna-loading">Analyzing via API…</p> : null}
      {error ? <ScientificError error={error} /> : null}
      {!envelope && !loading && !error ? (
        <EmptyState title="No DNA result yet">Paste a sequence and analyze.</EmptyState>
      ) : null}
      {invalid ? <InvalidMoleculeResult expected="canonical DNA (A, T, C, G), optionally as single-record FASTA" result={result} /> : null}
      {computed ? (
        <>
          <AnalysisNav
            items={[
              { id: "dna-summary", label: "Overview" },
              { id: "dna-composition", label: "Composition", available: compositionBar.length > 0 },
              { id: "dna-properties", label: "Properties" },
              { id: "dna-thermo", label: "Thermodynamics" },
              { id: "dna-profiles", label: "Profiles", available: profiles.length > 0 },
              { id: "dna-features", label: "Features" },
              { id: "dna-orfs-section", label: "ORFs", available: draft.includeOrfs },
              { id: "dna-cpg", label: "CpG", available: draft.includeCpg },
              { id: "dna-restriction", label: "Restriction", available: draft.includeRestriction },
              { id: "dna-kmers", label: "k-mers", available: kmerTop.length > 0 },
              { id: "dna-meaning", label: "Interpretation" },
              { id: "dna-raw", label: "Raw data" },
              { id: "dna-provenance", label: "Provenance" },
            ]}
          />
          <AnonymousSequenceNotice />
          <BeginnerOnly>
            <p className="hs-lede">
              After analysis, the sequence recedes and the computed DNA result becomes the investigation
              surface: composition, profiles, predicted features, then raw values and provenance.
            </p>
          </BeginnerOnly>
          <ResultAnchor id="dna-summary">
          <section className="hs-panel science-surface">
            <h2>Summary</h2>
            <MetricGrid
              metrics={[
                { id: "length", label: "Length", value: result.length, unit: "bp" },
                { id: "gc", label: "GC", value: result.gc_content, unit: "%" },
                { id: "at", label: "AT", value: result.at_content, unit: "%" },
                { id: "tm", label: "Tm", value: tm.value_c, unit: "C" },
                { id: "mw", label: "MW", value: result.molecular_weight_dna, unit: "g/mol" },
                { id: "entropy", label: "Entropy", value: result.shannon_entropy, unit: "bits" },
                { id: "status", label: "Status", value: result.status },
                { id: "hash", label: "Sequence hash", value: String(result.sequence_hash || "").slice(0, 16) },
                { id: "gc_skew", label: "GC skew", value: result.gc_skew },
                { id: "at_skew", label: "AT skew", value: result.at_skew },
                { id: "cpg_oe", label: "CpG O/E", value: cpgOe },
              ]}
            />
            <p className="hs-topbar-meta">
              Tm method {String(tm.method || "from API")}. Wallace for fewer than 14 canonical bp and
              salt-adjusted (0.05 M Na+) from 14 to 200 bp. That value is not nearest-neighbor.
              Sequences longer than 200 bp or without canonical bases are N/A, not 0. GC skew is
              (G-C)/(G+C); AT skew is (A-T)/(A+T). A zero denominator is Undefined, not 0. Composition
              consistent: {String(result.composition_consistent)}.
            </p>
            {result.reverse_complement ? (
              <p data-testid="dna-reverse-complement" className="mono">
                Reverse complement {String(result.reverse_complement)}
              </p>
            ) : null}
            <p className="hs-topbar-meta">{String(result.identity_note || "")}</p>
          </section>
          </ResultAnchor>
          {compositionBar.length > 0 ? (
            <ResultAnchor id="dna-composition">
            <section className="hs-panel science-surface">
              <h2>Composition</h2>
              <InfographicFrame
                title="Nucleotide composition"
                dataSource="HelixScope Core composition counts"
                method="Direct base counts on the normalized DNA string"
                units="count"
                status={result.status}
                shows="How many A, C, G and T residues are in the analyzed sequence"
                doesNotShow="Organism, gene identity, or sequencing quality"
              >
                <ScienceBarChart
                  labels={compositionBar.map((row) => row.symbol)}
                  values={compositionBar.map((row) => row.count)}
                  xTitle="base"
                  yTitle="count"
                  title="Nucleotide counts (API composition)"
                  testId="dna-composition-bar"
                />
              </InfographicFrame>
              <DataTable rows={compositionRows} caption="Nucleotide composition" testId="dna-composition" />
            </section>
            </ResultAnchor>
          ) : null}
          {dinucRows.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Dinucleotide frequencies and CpG odds ratio</h2>
              <InfographicFrame
                title="Dinucleotide observed/expected"
                dataSource="HelixScope Core dinucleotide_frequencies"
                method="Observed pair counts divided by independent-base expectation"
                units="O/E ratio"
                status={result.status}
                shows="Which adjacent base pairs are over- or under-represented"
                doesNotShow="Methylation, chromatin, or experimental CpG island assays. Zero expected frequency is Undefined, not 0."
              >
                <ScienceHeatmap
                  z={dinucHeat.z}
                  x={dinucHeat.x}
                  y={dinucHeat.y}
                  zmid={1}
                  title="Dinucleotide observed/expected (API O/E; NaN is a gap, not 0)"
                  testId="dna-dinuc-heatmap"
                />
              </InfographicFrame>
              <DataTable rows={dinucRows} caption="Dinucleotides" testId="dna-dinucleotides" />
              <button type="button" className="btn-secondary" onClick={() => downloadCsv("dna_dinucleotides.csv", dinucRows)}>
                Download dinucleotides (CSV)
              </button>
              <p className="hs-topbar-meta">
                Observed/Expected near 1.0 means adjacent bases behave independently. When the expected
                frequency is zero the ratio is Undefined, not 0.
              </p>
            </section>
          ) : null}
          <ResultAnchor id="dna-thermo">
          {Object.keys(santalucia).length > 0 || tm.value_c !== undefined ? (
            <section className="hs-panel science-surface">
              <h2>Thermodynamics</h2>
              <p className="hs-topbar-meta">
                Wallace or salt-adjusted Tm applies only inside the implemented length domain. SantaLucia
                nearest-neighbor is a separate optional method. A missing Tm is N/A, not 0 C.
              </p>
              {Object.keys(santalucia).length > 0 ? (
                <MetricGrid
                  metrics={[
                    { id: "sl_status", label: "SantaLucia status", value: santalucia.status },
                    { id: "sl_tm", label: "SantaLucia Tm", value: santalucia.value_c, unit: "C" },
                    { id: "sl_method", label: "Method", value: santalucia.method },
                  ]}
                />
              ) : null}
              <p className="hs-topbar-meta">{String(santalucia.reason || santalucia.parameter_set || "")}</p>
            </section>
          ) : null}
          </ResultAnchor>
          {draft.includeProfiles && profiles.length === 0 ? (
            <p data-testid="dna-profiles-too-short">
              Windowed GC, AT, skew and entropy profiles require a sequence at least as long as the API
              window ({draft.profileWindow} nt). This result has 0 windowed profile rows.
            </p>
          ) : null}
          {profiles.length > 0 ? (
            <ResultAnchor id="dna-profiles">
            <section className="hs-panel science-surface">
              <h2>Windowed profiles</h2>
              <InfographicFrame
                title="Windowed GC"
                dataSource="HelixScope Core windowed_profiles"
                method={`Sliding window ${draft.profileWindow} bp, step ${draft.profileStep} bp`}
                units="% GC"
                status={result.status}
                shows="Local GC along the analyzed coordinates"
                doesNotShow="Isochores, experimental melting maps, or interpolated missing windows"
              >
                <ScienceChart
                  x={x}
                  y={profiles.map((row) => row.gc_percent)}
                  xTitle="position"
                  yTitle="GC %"
                  title="Windowed GC (API arrays)"
                  testId="dna-gc-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Windowed AT"
                dataSource="HelixScope Core windowed_profiles"
                method={`Sliding window ${draft.profileWindow} bp, step ${draft.profileStep} bp`}
                units="% AT"
                status={result.status}
                shows="Local AT along the analyzed coordinates"
                doesNotShow="Strand-specific transcription bias"
              >
                <ScienceChart
                  x={x}
                  y={profiles.map((row) => row.at_percent)}
                  xTitle="position"
                  yTitle="AT %"
                  title="Windowed AT (API arrays)"
                  testId="dna-at-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Windowed GC skew"
                dataSource="HelixScope Core windowed_profiles"
                method="(G-C)/(G+C) per window; zero denominator is Undefined"
                units="dimensionless"
                status={result.status}
                shows="Local G versus C excess"
                doesNotShow="Confirmed replication origin"
              >
                <ScienceChart
                  x={x}
                  y={profiles.map((row) => row.gc_skew)}
                  xTitle="position"
                  yTitle="GC skew"
                  title="Windowed GC skew (API arrays)"
                  testId="dna-gc-skew-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Windowed AT skew"
                dataSource="HelixScope Core windowed_profiles"
                method="(A-T)/(A+T) per window; zero denominator is Undefined"
                units="dimensionless"
                status={result.status}
                shows="Local A versus T excess"
                doesNotShow="Confirmed transcription direction"
              >
                <ScienceChart
                  x={x}
                  y={profiles.map((row) => row.at_skew)}
                  xTitle="position"
                  yTitle="AT skew"
                  title="Windowed AT skew (API arrays)"
                  testId="dna-at-skew-chart"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Windowed Shannon entropy"
                dataSource="HelixScope Core windowed_profiles"
                method="Shannon entropy of window base frequencies"
                units="bits"
                status={result.status}
                shows="Local sequence complexity of the window alphabet"
                doesNotShow="Biological information content or motif function"
              >
                <ScienceChart
                  x={x}
                  y={profiles.map((row) => row.entropy)}
                  xTitle="position"
                  yTitle="entropy (bits)"
                  title="Windowed Shannon entropy (API arrays)"
                  testId="dna-entropy-chart"
                />
              </InfographicFrame>
              <DataTable rows={profiles} caption="Windowed profile table" testId="dna-window-table" />
              <button type="button" className="btn-secondary" onClick={() => downloadCsv("dna_windowed_profiles.csv", profiles)}>
                Download windowed profiles (CSV)
              </button>
            </section>
            </ResultAnchor>
          ) : null}
          {draft.includeKmers && kmerTop.length > 0 ? (
            <ResultAnchor id="dna-kmers">
            <section className="hs-panel science-surface">
              <h2>Entropy, k-mers and AT skew</h2>
              <InfographicFrame
                title="Most abundant k-mers"
                dataSource="HelixScope Core kmer_summary"
                method={`Exact k-mer counts, k=${String(kmers.k)}`}
                units="count"
                status={kmers.status}
                shows="Which oligomers occur most often in the normalized sequence"
                doesNotShow="Biological motif function or genome-wide enrichment"
              >
                <ScienceBarChart
                  labels={kmerTop.map((row) => row.kmer)}
                  values={kmerTop.map((row) => row.count)}
                  xTitle="k-mer"
                  yTitle="count"
                  title="Most abundant k-mers (API ranking)"
                  testId="dna-kmer-bar"
                />
              </InfographicFrame>
              <p className="hs-topbar-meta">
                k={String(kmers.k)} · valid windows {String(kmers.valid_windows)} · unique{" "}
                {String(kmers.unique_kmers)} / {String(kmers.possible_kmers)} · diversity{" "}
                {String(kmers.diversity)} · status {String(kmers.status)}
              </p>
              <DataTable rows={kmerTop} caption="k-mer summary" testId="dna-kmers" />
            </section>
            </ResultAnchor>
          ) : null}
          {tracks.length > 0 || orfs.length > 0 || cpg.length > 0 || restrictionHits.length > 0 ? (
            <ResultAnchor id="dna-features">
            <section className="hs-panel science-surface">
              <h2>Linked sequence track</h2>
              <p className="hs-topbar-meta">
                ORF, CpG and restriction intervals use API coordinates. Click a table row to highlight the
                matching interval when start/end identity is exact.
              </p>
              <GenomicTrack
                length={Number(result.length) || 0}
                sequence={String(result.sequence || "")}
                tracks={tracks}
                profiles={profiles}
                selected={selectedFeature}
                onSelect={(feature) => setSelectedFeature(feature)}
              />
            </section>
            </ResultAnchor>
          ) : null}
          {tracks.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Sequence map (predicted ORFs)</h2>
              <p className="hs-topbar-meta">
                Predicted ORFs, not confirmed genes. Coordinates are mapped onto the input strand.
              </p>
              <InfographicFrame
                title="Predicted ORF track"
                dataSource="HelixScope Core sequence_feature_tracks"
                method="ATG start, in-frame stop, minimum length, optional non-overlap filter"
                units="bp coordinates"
                status="PREDICTED"
                shows="Where predicted open reading frames sit on the pasted sequence"
                doesNotShow="Confirmed genes, promoters, or spliced transcripts"
              >
                <ScienceFeatureMap
                  length={Number(result.length) || 0}
                  tracks={tracks}
                  testId="dna-feature-map"
                />
              </InfographicFrame>
            </section>
          ) : null}
          <section className="hs-panel science-surface">
            <h2>Detailed analysis</h2>
            {draft.includeOrfs ? (
              <ResultAnchor id="dna-orfs-section">
              <>
                <p className="hs-topbar-meta">
                  Raw ORF candidates {String(result.orfs_raw_count ?? orfs.length)}. Reported {orfs.length}.
                  Nested filter {String(result.orfs_drop_overlapping)}. Status of every row is PREDICTED,
                  not gene.
                </p>
                {orfs.length > 0 ? (
                  <InfographicFrame
                    title="Predicted ORF lengths"
                    dataSource="HelixScope Core orfs"
                    method="Length of each predicted ORF in nucleotides"
                    units="nt"
                    status="PREDICTED"
                    shows="Distribution of predicted ORF lengths"
                    doesNotShow="Protein-coding proof or expression"
                  >
                    <ScienceHistogram
                      values={orfs.map((row) => row.length_bp ?? row.length)}
                      xTitle="ORF length (nt)"
                      yTitle="count"
                      title="Predicted ORF length histogram (API lengths)"
                      testId="dna-orf-hist"
                    />
                  </InfographicFrame>
                ) : null}
                <DataTable
                  rows={orfs}
                  caption="ORFs"
                  testId="dna-orfs"
                  emptyLabel="0 ORFs."
                  selectedKey={selectedFeature?.kind === "orf" ? String(selectedFeature.start) : null}
                  onRowSelect={(row) =>
                    setSelectedFeature({
                      start: Number(row.start),
                      end: Number(row.end),
                      kind: "orf",
                    })
                  }
                />
                {orfs.length > 0 ? (
                  <button type="button" className="btn-secondary" onClick={() => downloadCsv("dna_orfs.csv", orfs)}>
                    Download ORFs (CSV)
                  </button>
                ) : null}
              </>
              </ResultAnchor>
            ) : null}
            {draft.includeRestriction ? (
              <ResultAnchor id="dna-restriction">
              <>
                <DataTable rows={enzymeRows} caption="Restriction enzymes" testId="dna-restriction" emptyLabel="0 restriction sites." />
                <DataTable
                  rows={restrictionHits}
                  caption="Restriction hits"
                  testId="dna-restriction-hits"
                  emptyLabel="0 restriction hits."
                  onRowSelect={(row) =>
                    setSelectedFeature({
                      start: Number(row.start ?? row.position),
                      end: Number(row.end ?? row.position) + 1,
                      kind: "restriction",
                    })
                  }
                />
              </>
              </ResultAnchor>
            ) : null}
            {draft.includeCpg ? (
              <ResultAnchor id="dna-cpg">
              <DataTable
                rows={cpg}
                caption="CpG islands"
                testId="dna-cpg"
                emptyLabel="0 CpG islands."
                onRowSelect={(row) =>
                  setSelectedFeature({
                    start: Number(row.start),
                    end: Number(row.end),
                    kind: "cpg",
                  })
                }
              />
              </ResultAnchor>
            ) : null}
          </section>
          <ResultAnchor id="dna-properties">
          <SequenceViewer sequence={String(result.sequence || "")} hash={String(result.sequence_hash || "")} testId="dna-sequence-viewer" />
          </ResultAnchor>
          <section className="hs-panel science-surface hs-viz-hero">
            <h2>DNA 3D</h2>
            <p className="hs-lede">
              A DNA sequence is not a deposited structure. The illustrative B-DNA helix uses fibre
              parameters and stays ILLUSTRATIVE. Atomic detail is a consecutive window capped by LOD;
              omitted flanks are not interpolated. 1BNA is an experimental duplex from the bundled
              RCSB copy; it is not remapped onto this analysis sequence unless hashes match.
            </p>
            <div className="hs-numeric-row">
              <div className="field">
                <label htmlFor="dna-3d-lod">Visual LOD</label>
                <select id="dna-3d-lod" value={helixLod} onChange={(e) => setHelixLod(e.target.value)}>
                  <option value="high">high</option>
                  <option value="medium">medium</option>
                  <option value="low">low</option>
                </select>
              </div>
              <div className="field">
                <label htmlFor="dna-3d-start">Helix start (0-based)</label>
                <input
                  id="dna-3d-start"
                  type="number"
                  min={0}
                  max={Math.max(0, seqLen - 1)}
                  value={helixStart}
                  onChange={(e) => setHelixStart(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="dna-3d-end">Helix end (exclusive)</label>
                <input
                  id="dna-3d-end"
                  type="number"
                  min={1}
                  max={Math.max(1, seqLen)}
                  value={helixEnd}
                  onChange={(e) => setHelixEnd(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="dna-3d-pos">Inspect DNA position (0-based)</label>
                <input
                  id="dna-3d-pos"
                  type="number"
                  min={0}
                  max={Math.max(0, seqLen - 1)}
                  value={helixInspect}
                  onChange={(e) => setHelixInspect(Number(e.target.value))}
                />
              </div>
              <div className="field">
                <label htmlFor="dna-3d-prefer">DNA 3D path</label>
                <select id="dna-3d-prefer" value={helixPrefer} onChange={(e) => setHelixPrefer(e.target.value)}>
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
                    const { data } = await helixApi.dnaHelix3d({
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
                Load illustrative B-DNA helix
              </button>
              <button
                type="button"
                className="btn-secondary"
                onClick={() => {
                  patchDraft("viewer", { structureId: "1BNA" });
                  router.push("/structure/viewer");
                }}
              >
                Open experimental 1BNA
              </button>
            </div>
            {helix ? (
              <>
                <p data-testid="dna-3d-kind">
                  kind {String(helix.kind || "")} · status {String(helix.status || "")}
                  {helix.partial_view ? " · PARTIAL_VIEW" : ""}
                </p>
                <p data-testid="dna-3d-mapping">
                  Sequence position {String(helix.inspect_position)} base {String(helix.inspect_base ?? "N/A")} mapping{" "}
                  {String(asRecord(helix.position_mapping).status || "")}
                </p>
                <p className="hs-topbar-meta">{String(helix.disclaimer || helix.kind_label || "")}</p>
                {helixViewer.coordinates || helixViewer.atoms ? (
                  <StructureViewportPlotly
                    scene={{ ...helixViewer, kind: helix.kind, status: helix.status, structure_id: "illustrative-bdna" }}
                    orbit={orbit}
                    onOrbit={setOrbit}
                  />
                ) : null}
              </>
            ) : null}
          </section>
          <ResultAnchor id="dna-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="dna-raw">
          <RawDataPanel
            result={result}
            fasta={fastaRecord(String(result.sequence_hash || "dna"), String(result.sequence || draft.sequence))}
            fastaName="helixscope-dna.fasta"
          />
          </ResultAnchor>
          <ResultAnchor id="dna-provenance">
          <ProvenancePanel envelope={envelope} />
          </ResultAnchor>
        </>
      ) : null}
    </div>
  );
}
