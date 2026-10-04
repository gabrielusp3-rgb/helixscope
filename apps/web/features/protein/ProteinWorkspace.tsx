"use client";

import dynamic from "next/dynamic";
import { useRouter } from "next/navigation";
import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, asRecordRows, chartNumeric } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import {
  EmptyState,
  ExternalDisclosure,
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
import { ScienceBarChart } from "@/components/charts/ScienceBarChart";
import { ScienceChart } from "@/components/charts/ScienceChart";
import { ScienceRadar } from "@/components/charts/ScienceRadar";
import { downloadJson, downloadText, fastaRecord } from "@/lib/export/download";
import type { Orbit } from "@/lib/gesture/math";

const StructureViewport = dynamic(
  () => import("@/components/structure/StructureViewport").then((m) => m.StructureViewport),
  { ssr: false },
);

const DEFAULT_ORBIT: Orbit = { yaw: 0.7, pitch: 0.35, radius: 2.4 };

export function ProteinWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.protein;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [loading, setLoading] = useState(false);
  const [caps, setCaps] = useState<Record<string, unknown>>({});
  const [orbit, setOrbit] = useState<Orbit>(DEFAULT_ORBIT);
  const [loadedScene, setLoadedScene] = useState<Record<string, unknown> | null>(null);
  const envelope = results.protein;
  const uniprotEnvelope = results.proteinUniprot;
  const structureEnvelope = results.proteinStructure;
  const result = asRecord(envelope?.result);
  const uniprot = asRecord(uniprotEnvelope?.result);
  const structure = asRecord(structureEnvelope?.result);
  const ss = asRecord(result.ss_prediction);
  const hydropathy = asRecordRows(result.hydropathy_profile, "midpoint");
  const charge = asRecordRows(result.charge_profile, "midpoint");
  const composition = asRecordRows(result.composition, "symbol");
  const categories = asRecordRows(result.amino_acid_categories, "category");
  const status = String(result.status || "");
  const invalid = status === "ERROR";
  const computed = Boolean(envelope) && !error && !invalid;
  const busy = loading || scienceLocked;

  useEffect(() => {
    void helixApi.proteinStructureCapabilities().then(({ data }) => setCaps(asRecord(data.result))).catch(() => null);
  }, []);

  async function run() {
    const gate = next();
    setLoading(true);
    setError(null);
    setResult("protein", undefined);
    try {
      const { data } = await helixApi.proteinAnalyze(
        { sequence: draft.sequence, ph: Number(draft.ph) || 7, include_explanation: true },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("protein", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    } finally {
      if (gate.isCurrent()) setLoading(false);
    }
  }

  async function retrieveUniprot() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.proteinUniprot(
        {
          accession: draft.uniprotAccession,
          include_domains: draft.includeDomains,
          include_explanation: true,
        },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("proteinUniprot", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  return (
    <div className="hs-page" data-testid="hs-page">
      <p className="kicker">Analysis</p>
      <h1>Protein</h1>
      <p className="hs-lede">
        Physicochemical values come from HelixScope API v1. Single-record FASTA headers are stripped
        before the 20 canonical amino-acid alphabet is applied. Sequence secondary-structure predictors
        remain UNAVAILABLE.
      </p>
      <CompactInput compact={computed} summary="Edit protein input">
        <SequenceInput
          id="prot-seq"
          label="Protein sequence"
          hint="20 canonical amino-acid symbols (ACDEFGHIKLMNPQRSTVWY), or a single-record FASTA. O and U are noncanonical and are rejected. FASTA description lines are not residues."
          value={draft.sequence}
          onChange={(sequence) => patchDraft("protein", { sequence })}
          disabled={busy}
        />
        <FastaFileField id="prot-file" disabled={busy} onText={(sequence) => patchDraft("protein", { sequence })} />
        <div className="field">
          <label htmlFor="ph">pH for net charge</label>
          <input id="ph" value={draft.ph} onChange={(e) => patchDraft("protein", { ph: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="uniprot-acc">UniProtKB accession (optional, retrieved annotation)</label>
          <input
            id="uniprot-acc"
            value={draft.uniprotAccession}
            onChange={(e) => patchDraft("protein", { uniprotAccession: e.target.value })}
            data-testid="protein-uniprot-accession"
          />
        </div>
        <label>
          <input
            type="checkbox"
            checked={draft.includeDomains}
            onChange={(e) => patchDraft("protein", { includeDomains: e.target.checked })}
          />{" "}
          Include InterPro domains for that accession
        </label>
        {draft.uniprotAccession.trim() ? (
          <ExternalDisclosure service="UniProt REST (rest.uniprot.org)">
            Only the accession is sent. An anonymous pasted sequence is not mapped automatically.
          </ExternalDisclosure>
        ) : null}
        <div className="hs-actions">
          <button type="button" className="btn-primary" disabled={busy} onClick={() => void run()}>
            Analyze
          </button>
          <button
            type="button"
            className="btn-secondary"
            disabled={busy || !draft.uniprotAccession.trim()}
            onClick={() => void retrieveUniprot()}
            data-testid="protein-uniprot-retrieve"
          >
            Retrieve UniProt entry
          </button>
          {computed ? (
            <>
              <button
                type="button"
                className="btn-secondary"
                onClick={() =>
                  downloadText(
                    "helixscope-protein.fasta",
                    fastaRecord(String(result.sequence_hash || "protein"), String(result.sequence || draft.sequence)),
                  )
                }
              >
                Download FASTA
              </button>
              <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-protein.json", result)}>
                Download JSON
              </button>
            </>
          ) : null}
        </div>
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      <section className="hs-panel science-surface">
        <h2>Structure sources</h2>
        <p className="hs-lede">
          Structure lookup stays visible after sequence analysis. Bundled 1CRN is experimental. HelixScope
          does not invent a fold from the pasted sequence.
        </p>
        <div className="hs-actions">
          <button
            type="button"
            className="btn-secondary"
            onClick={() => {
              recordTransfer({
                kind: "protein-to-structure",
                note: "Open 3D viewer after protein analysis",
                payload: { sequence: draft.sequence },
              });
              patchDraft("viewer", { structureId: "1CRN" });
              router.push("/structure/viewer");
            }}
          >
            Open experimental 1CRN
          </button>
          <button
            type="button"
            className="btn-secondary"
            data-testid="protein-structure-search"
            onClick={async () => {
              setError(null);
              try {
                const { data } = await helixApi.proteinStructureSearch({
                  sequence: draft.sequence,
                  molecule: "PROTEIN",
                  pdb_id: "1CRN",
                  search_pdb_by_sequence: false,
                });
                setResult("proteinStructure", data);
                recordTransfer({
                  kind: "protein-to-structure",
                  note: "Structure source lookup from Protein sequence",
                  payload: { pdb_id: "1CRN" },
                });
              } catch (err) {
                setError(err);
              }
            }}
          >
            Search structure sources
          </button>
          <button
            type="button"
            className="btn-secondary"
            data-testid="protein-structure-analyze"
            onClick={async () => {
              setError(null);
              try {
                const { data } = await helixApi.proteinStructureAnalyze({
                  structure_id: "1CRN",
                  run_sasa: true,
                  run_dssp: false,
                  run_edtsurf: false,
                });
                setResult("proteinStructure", data);
              } catch (err) {
                setError(err);
              }
            }}
          >
            Run bundled 1CRN SASA
          </button>
        </div>
      </section>
      <p data-testid="protein-structure-caps">
        DSSP local {String(asRecord(caps.dssp_local).status || asRecord(caps.dssp_local).available || "pending")} ·
        PDB-REDO {String(asRecord(caps.dssp_pdb_redo).status || "pending")} · SASA {String(asRecord(caps.sasa).status || "pending")} ·
        EDTSurf {String(asRecord(caps.edtsurf).status || asRecord(caps.edtsurf).available || "pending")}
      </p>
      {structureEnvelope ? (
        <section className="hs-panel science-surface" data-testid="protein-structure-result">
          <h2>Structure sources and analysis</h2>
          <p>
            status {String(structure.status)} · kind {String(structure.kind || "")} · source{" "}
            {String(structure.source || "")}
          </p>
          <p className="hs-lede">
            Hits are retrieved records. HelixScope does not pick a best structure. Experimental and
            predicted candidates stay labeled separately.
          </p>
          <DataTable
            rows={asList(structure.hits).map((row) => asRecord(row))}
            caption="Structure candidates"
            testId="protein-structure-hits"
            emptyLabel="0 structure candidates."
            onRowSelect={(row) => {
              const id = String(row.structure_id || row.pdb_id || "");
              if (!id) return;
              void (async () => {
                setError(null);
                try {
                  const { data } = await helixApi.proteinStructureLoad({
                    structure_id: id,
                    sequence: draft.sequence,
                    molecule: "PROTEIN",
                    bundled: ["1CRN", "1BNA", "1RNA"].includes(id.toUpperCase()),
                  });
                  setResult("proteinStructure", data);
                  const scene = asRecord(asRecord(data.result).viewer);
                  setLoadedScene({
                    ...scene,
                    kind: asRecord(data.result).kind,
                    status: asRecord(data.result).status,
                    structure_id: asRecord(data.result).structure_id || id,
                    source: asRecord(data.result).source,
                  });
                  recordTransfer({
                    kind: "protein-to-structure",
                    note: "Loaded structure candidate into 3D viewport",
                    payload: { structure_id: id },
                  });
                } catch (err) {
                  setError(err);
                }
              })();
            }}
          />
          <MetricGrid
            metrics={[
              { id: "sasa-status", label: "SASA", value: asRecord(structure.sasa).status },
              { id: "dssp-status", label: "DSSP", value: asRecord(structure.dssp ?? asRecord(caps.dssp_local)).status },
              { id: "dssp-redo", label: "PDB-REDO DSSP", value: asRecord(structure.dssp_pdb_redo ?? asRecord(caps.dssp_pdb_redo)).status },
              { id: "edtsurf-status", label: "EDTSurf", value: asRecord(structure.edtsurf ?? asRecord(caps.edtsurf)).status },
            ]}
          />
          <p>SASA {String(asRecord(structure.sasa).status || "")}</p>
          {loadedScene && (loadedScene.coordinates || loadedScene.atoms) ? (
            <StructureViewport scene={loadedScene} orbit={orbit} onOrbit={setOrbit} />
          ) : null}
        </section>
      ) : null}
      {invalid ? <InvalidMoleculeResult expected="20 canonical amino acids, optionally as single-record FASTA" result={result} /> : null}
      {!envelope && !uniprotEnvelope && !error ? <EmptyState title="No protein result yet" /> : null}
      {computed ? (
        <>
          <AnalysisNav
            items={[
              { id: "protein-summary", label: "Overview" },
              { id: "protein-composition", label: "Composition", available: composition.length > 0 },
              { id: "protein-profiles", label: "Profiles", available: hydropathy.length > 0 || charge.length > 0 },
              { id: "protein-uniprot", label: "UniProt", available: Boolean(uniprotEnvelope) },
              { id: "protein-structure", label: "Structure", available: Boolean(structureEnvelope) },
              { id: "protein-meaning", label: "Interpretation" },
              { id: "protein-raw", label: "Raw data" },
            ]}
          />
          <AnonymousSequenceNotice />
          <ResultAnchor id="protein-summary">
          <MetricGrid
            metrics={[
              { id: "len", label: "Length", value: result.length, unit: "aa" },
              { id: "mwk", label: "MW", value: result.molecular_weight_kda, unit: "kDa" },
              { id: "pi", label: "pI", value: result.isoelectric_point },
              { id: "gravy", label: "GRAVY", value: result.gravy },
              { id: "charge", label: "Net charge", value: result.net_charge },
              { id: "inst", label: "Instability", value: result.instability_index },
              { id: "ali", label: "Aliphatic index", value: result.aliphatic_index },
              { id: "ss", label: "SS prediction", value: ss.available === false ? "UNAVAILABLE" : ss.status },
            ]}
          />
          <p className="hs-topbar-meta">{String(result.identity_note || "")}</p>
          </ResultAnchor>
          {composition.length > 0 ? (
            <ResultAnchor id="protein-composition">
            <section className="hs-panel science-surface">
              <h2>Amino-acid composition</h2>
              <InfographicFrame
                title="Amino-acid counts"
                dataSource="HelixScope Core amino_acid_composition"
                method="Direct residue counts on the canonical protein string"
                units="count"
                status={result.status}
                shows="How often each of the 20 canonical residues occurs"
                doesNotShow="Domains, organism, or UniProt identity"
              >
                <ScienceBarChart
                  labels={composition.map((row) => row.symbol)}
                  values={composition.map((row) => row.count ?? row.frequency)}
                  title="Amino-acid counts (API composition)"
                  testId="protein-aa-bar"
                />
              </InfographicFrame>
              <InfographicFrame
                title="Amino-acid frequencies"
                dataSource="HelixScope Core amino_acid_composition"
                method="Residue count divided by length"
                units="frequency"
                status={result.status}
                shows="Relative residue frequencies"
                doesNotShow="Structure or function"
              >
                <ScienceRadar
                  labels={composition.map((row) => row.symbol)}
                  values={composition.map((row) => row.frequency)}
                  title="Amino-acid frequency radar (API frequencies)"
                  testId="protein-aa-radar"
                />
              </InfographicFrame>
              <DataTable rows={composition} caption="Amino-acid composition" testId="protein-composition" />
            </section>
            </ResultAnchor>
          ) : null}
          {categories.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Side-chain categories</h2>
              <ScienceBarChart
                labels={categories.map((row) => row.category)}
                values={categories.map((row) => row.frequency)}
                title="Amino-acid categories (API grouping)"
                testId="protein-category-bar"
              />
              <DataTable rows={categories} caption="Categories" testId="protein-categories" />
            </section>
          ) : null}
          {hydropathy.length > 0 ? (
            <ResultAnchor id="protein-profiles">
            <section className="hs-panel science-surface">
              <h2>Kyte-Doolittle hydropathy</h2>
              <InfographicFrame
                title="Kyte-Doolittle profile"
                dataSource="HelixScope Core hydropathy_profile"
                method="Kyte-Doolittle windowed scores from Core"
                units="hydropathy score"
                status={result.status}
                shows="Local hydrophobicity along the residue string"
                doesNotShow="Transmembrane helix proof or experimental topology"
              >
                <ScienceChart
                  x={hydropathy.map((row) => row.midpoint)}
                  y={hydropathy.map((row) => row.score)}
                  xTitle="midpoint"
                  yTitle="hydropathy"
                  title="Kyte-Doolittle profile (API arrays)"
                  testId="protein-hydropathy"
                />
              </InfographicFrame>
            </section>
            </ResultAnchor>
          ) : null}
          {charge.length > 0 ? (
            <section className="hs-panel science-surface">
              <h2>Formal charge profile</h2>
              <ScienceChart
                x={charge.map((row) => row.midpoint)}
                y={charge.map((row) => row.score)}
                xTitle="midpoint"
                yTitle="formal charge / window"
                title="Side-chain formal charge profile (API arrays; His excluded)"
                testId="protein-charge"
              />
              <p className="hs-topbar-meta">
                This is a residue-count profile, not Henderson-Hasselbalch net charge. Window{" "}
                {String(charge[0]?.window)}. Status {String(result.charge_profile_status)}.
              </p>
            </section>
          ) : null}
          {chartNumeric(result.length) ? (
            <SequenceViewer sequence={String(result.sequence || "")} hash={String(result.sequence_hash || "")} testId="protein-sequence-viewer" />
          ) : null}
          <ResultAnchor id="protein-meaning">
          <InterpretationPanel explanation={result.explanation} />
          <MetricExplanations explanations={result.metric_explanations} />
          </ResultAnchor>
          <ResultAnchor id="protein-raw">
          <RawDataPanel
            result={result}
            fasta={fastaRecord(String(result.sequence_hash || "protein"), String(result.sequence || draft.sequence))}
            fastaName="helixscope-protein.fasta"
          />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} />
        </>
      ) : null}
      {uniprotEnvelope ? (
        <ResultAnchor id="protein-uniprot">
        <section className="hs-panel science-surface" data-testid="protein-uniprot-panel">
          <h2>UniProt retrieved annotation</h2>
          <p className="hs-lede">
            These fields are RETRIEVED from UniProtKB for the requested accession. They are not
            computed from the pasted sequence, and they are not assigned to an anonymous paste.
          </p>
          <MetricGrid
            metrics={[
              { id: "up-acc", label: "Accession", value: uniprot.accession },
              { id: "up-name", label: "Protein name", value: uniprot.protein_name },
              { id: "up-gene", label: "Gene names", value: Array.isArray(uniprot.gene_names) ? uniprot.gene_names.join(", ") : uniprot.gene_names },
              { id: "up-org", label: "Organism", value: uniprot.organism },
              { id: "up-status", label: "Status", value: uniprot.status ?? uniprot.evidence_status },
            ]}
          />
          <DataTable rows={asRecordRows(uniprot.features, "type")} caption="UniProt features" testId="protein-uniprot-features" />
          <DataTable rows={asRecordRows(uniprot.cross_references, "database")} caption="Cross-references" />
          <DataTable rows={asRecordRows(uniprot.literature, "title")} caption="Literature" />
          {uniprot.retrieved_sequence ? (
            <button
              type="button"
              className="btn-secondary"
              onClick={() => {
                recordTransfer({
                  kind: "uniprot-to-protein",
                  note: "Explicit UniProt retrieved sequence load",
                  payload: { accession: String(uniprot.accession || "") },
                });
                patchDraft("protein", { sequence: String(uniprot.retrieved_sequence) });
              }}
            >
              Load retrieved UniProt sequence into Protein analysis
            </button>
          ) : null}
          <InterpretationPanel explanation={uniprot.explanation} />
          <RawDataPanel result={uniprot} />
          <ProvenancePanel envelope={uniprotEnvelope} extra={uniprot.entry_audit} />
        </section>
        </ResultAnchor>
      ) : null}
    </div>
  );
}
