"use client";

import { useRouter } from "next/navigation";
import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace, type TransferKind } from "@/lib/workspace/WorkstationProvider";
import {
  ncbiAccession,
  ncbiIsPubmed,
  ncbiRecord,
  ncbiRetrievedSequence,
} from "@/lib/workspace/ncbiTransfer";
import {
  EmptyState,
  ExternalDisclosure,
  MetricGrid,
  ScientificError,
  SequenceInput,
} from "@/components/science/SciencePrimitives";
import { InterpretationPanel, ProvenancePanel } from "@/components/science/Provenance";
import { DataTable } from "@/components/science/DataTable";
import { CompactInput, SequenceViewer } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { downloadJson, downloadText, fastaRecord } from "@/lib/export/download";

export function NcbiWorkspace() {
  const router = useRouter();
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult, recordTransfer } = useWorkspace();
  const draft = drafts.ncbi;
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const envelope = results.ncbi;
  const result = asRecord(envelope?.result);
  const record = ncbiRecord(result);
  const features = asList(record.features ?? result.features).map((row) => asRecord(row));
  const references = asList(record.references ?? result.references).map((row) => asRecord(row));
  const hits = asList(result.hits ?? result.records).map((row) => asRecord(row));
  const sequence = ncbiRetrievedSequence(result);
  const accession = ncbiAccession(result);
  const pubmed = ncbiIsPubmed(result, draft.db);
  const authors = asList(record.authors).map((item) => String(item));
  const computed = Boolean(envelope) && !error;

  async function fetchAcc() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.ncbiFetch(
        { accession: draft.accession, email: draft.email, db: draft.db },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("ncbi", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  async function search() {
    const gate = next();
    setError(null);
    try {
      const { data } = await helixApi.ncbiSearch(
        { term: draft.term, email: draft.email, db: draft.db, retmax: 15 },
        gate.signal,
      );
      if (gate.isCurrent()) setResult("ncbi", data);
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  function sendSequence(kind: "dna" | "rna" | "protein" | "alignment" | "blast" | "msa") {
    if (!sequence) return;
    const fasta = fastaRecord(accession || "ncbi", sequence);
    recordTransfer({
      kind: (`ncbi-to-${kind}`) as TransferKind,
      note: "Explicit NCBI retrieved sequence copy. HSP strings are not used.",
      payload: { accession, db: draft.db },
    });
    if (kind === "dna") {
      patchDraft("dna", { sequence });
      router.push("/analysis/dna");
      return;
    }
    if (kind === "rna") {
      patchDraft("rna", { sequence });
      router.push("/analysis/rna");
      return;
    }
    if (kind === "protein") {
      patchDraft("protein", { sequence });
      router.push("/analysis/protein");
      return;
    }
    if (kind === "alignment") {
      patchDraft("alignment", { seq1: sequence });
      router.push("/analysis/alignment");
      return;
    }
    if (kind === "blast") {
      patchDraft("blast", { query: sequence, email: draft.email });
      router.push("/data/blast");
      return;
    }
    patchDraft("msa", { fasta, email: draft.email });
    router.push("/evolution/msa");
  }

  return (
    <div className="hs-page">
      <p className="kicker">Data</p>
      <h1>NCBI Fetch</h1>
      <ExternalDisclosure service="NCBI Entrez">The browser does not contact NCBI directly.</ExternalDisclosure>
      <CompactInput compact={computed} summary="Edit NCBI query">
        <div className="field">
          <label htmlFor="email">Contact email (NCBI policy)</label>
          <input id="email" value={draft.email} onChange={(e) => patchDraft("ncbi", { email: e.target.value })} />
        </div>
        <div className="field">
          <label htmlFor="db">Database</label>
          <select id="db" value={draft.db} onChange={(e) => patchDraft("ncbi", { db: e.target.value })}>
            <option value="nucleotide">nucleotide</option>
            <option value="protein">protein</option>
            <option value="gene">gene</option>
            <option value="pubmed">pubmed</option>
          </select>
        </div>
        <div className="field">
          <label htmlFor="acc">Accession or identifier</label>
          <input id="acc" value={draft.accession} onChange={(e) => patchDraft("ncbi", { accession: e.target.value })} />
        </div>
        <button type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void fetchAcc()}>
          Fetch record
        </button>
        <SequenceInput id="term" label="Search term" value={draft.term} onChange={(term) => patchDraft("ncbi", { term })} rows={3} />
        <button type="button" className="btn-secondary" disabled={scienceLocked} onClick={() => void search()}>
          Search
        </button>
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {!envelope ? (
        <EmptyState title="No NCBI result in this session" />
      ) : (
        <>
          <AnalysisNav
            items={[
              { id: "ncbi-record", label: "Record" },
              { id: "ncbi-pubmed", label: "PubMed", available: pubmed },
              { id: "ncbi-hits", label: "Hits", available: hits.length > 0 },
              { id: "ncbi-features", label: "Features", available: features.length > 0 },
              { id: "ncbi-refs", label: "References", available: references.length > 0 },
              { id: "ncbi-transfer", label: "Send to analysis", available: Boolean(sequence) },
              { id: "ncbi-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="ncbi-record">
          <section className="hs-panel science-surface">
            <h2>Record</h2>
            <MetricGrid
              metrics={[
                { id: "acc", label: "Accession", value: record.accession ?? record.Accession ?? record.accession_version },
                { id: "org", label: "Organism", value: record.organism ?? record.Organism },
                { id: "tax", label: "Taxonomy", value: record.taxonomy ?? record.Taxonomy },
                { id: "mol", label: "Molecule", value: record.molecule ?? record.mol_type },
                { id: "topo", label: "Topology", value: record.topology },
                { id: "len", label: "Length", value: record.length ?? record.Length },
                { id: "def", label: "Definition", value: record.definition ?? record.Definition ?? record.title },
              ]}
            />
          </section>
          </ResultAnchor>
          {pubmed ? (
            <ResultAnchor id="ncbi-pubmed">
            <section className="hs-panel science-surface" data-testid="ncbi-pubmed">
              <h2>PubMed metadata</h2>
              <p data-testid="pubmed-title">{String(record.title || record.definition || "")}</p>
              <p data-testid="pubmed-authors">{authors.length ? authors.join(", ") : "Authors unavailable"}</p>
              <p>
                {String(record.journal || "")} {String(record.publication_year || record.pubdate || "")}
              </p>
              <p data-testid="pubmed-abstract-status">
                Abstract {record.abstract_available ? "available from NCBI" : String(record.abstract_status || "not returned")}
              </p>
              {record.abstract ? <p>{String(record.abstract)}</p> : null}
              <p className="hs-topbar-meta">
                HelixScope retrieves bibliographic metadata and abstract text when NCBI returns them. It
                does not copy copyrighted full papers.
              </p>
            </section>
            </ResultAnchor>
          ) : null}
          {hits.length > 0 ? (
            <ResultAnchor id="ncbi-hits">
            <DataTable rows={hits} caption="Search hits" testId="ncbi-hits" />
            </ResultAnchor>
          ) : null}
          {features.length > 0 ? (
            <ResultAnchor id="ncbi-features">
            <DataTable rows={features} caption="Features" testId="ncbi-features" />
            </ResultAnchor>
          ) : null}
          {references.length > 0 ? (
            <ResultAnchor id="ncbi-refs">
            <DataTable rows={references} caption="References" testId="ncbi-references" />
            </ResultAnchor>
          ) : null}
          {sequence ? (
            <SequenceViewer sequence={sequence} testId="ncbi-sequence" />
          ) : (
            <pre className="mono">{String(record.sequence ?? result.sequence ?? "").slice(0, 2000)}</pre>
          )}
          {sequence ? (
            <ResultAnchor id="ncbi-transfer">
            <section className="hs-panel science-surface" data-testid="ncbi-workspace-transfer">
              <h2>Send retrieved sequence to a HelixScope workspace</h2>
              <p className="hs-lede">
                These actions copy the retrieved NCBI sequence. They do not mutate drafts until you click.
                BLAST HSP strings are never used here.
              </p>
              <div className="hs-actions">
                <button type="button" className="btn-secondary" data-testid="ncbi-to-dna" onClick={() => sendSequence("dna")}>
                  Analyze DNA
                </button>
                <button type="button" className="btn-secondary" data-testid="ncbi-to-rna" onClick={() => sendSequence("rna")}>
                  Analyze RNA
                </button>
                <button type="button" className="btn-secondary" data-testid="ncbi-to-protein" onClick={() => sendSequence("protein")}>
                  Analyze Protein
                </button>
                <button type="button" className="btn-secondary" onClick={() => sendSequence("alignment")}>
                  Send to Alignment
                </button>
                <button type="button" className="btn-secondary" onClick={() => sendSequence("blast")}>
                  Send to BLAST
                </button>
                <button type="button" className="btn-secondary" onClick={() => sendSequence("msa")}>
                  Add to MSA
                </button>
              </div>
            </section>
            </ResultAnchor>
          ) : (
            <p className="hs-topbar-meta">
              No retrieved sequence is present, so workspace transfer is unavailable. PubMed and some
              Gene summaries do not include a sequence.
            </p>
          )}
          {sequence ? (
            <button
              type="button"
              className="btn-secondary"
              onClick={() => downloadText("helixscope-ncbi.fasta", fastaRecord(accession || "ncbi", sequence))}
            >
              Download FASTA
            </button>
          ) : null}
          <button type="button" className="btn-secondary" onClick={() => downloadJson("helixscope-ncbi.json", result)}>
            Download JSON
          </button>
          <InterpretationPanel explanation={result.explanation} />
          <ResultAnchor id="ncbi-raw">
          <RawDataPanel result={result} extra={record} fasta={sequence ? fastaRecord(accession || "ncbi", sequence) : undefined} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} extra={record.provenance} />
        </>
      )}
    </div>
  );
}
