"use client";

import { useState } from "react";
import { helixApi } from "@/lib/api/client";
import { jobScientificResult, pollJob } from "@/lib/api/jobs";
import { asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useRequestGate } from "@/lib/hooks/useRequestGate";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { EmptyState, JobState, ScientificError, SequenceInput } from "@/components/science/SciencePrimitives";
import { InterpretationPanel, ProvenancePanel } from "@/components/science/Provenance";
import { CompactInput } from "@/components/science/ResultChrome";
import { AnalysisNav, ResultAnchor } from "@/components/science/AnalysisNav";
import { InfographicFrame } from "@/components/science/Infographic";
import { RawDataPanel } from "@/components/science/RawDataPanel";
import { TreeViewer } from "@/components/evolution/TreeViewer";
import { downloadText } from "@/lib/export/download";

const SYNC = new Set(["neighbor_joining", "upgma"]);

export function PhylogenyWorkspace() {
  const { scienceLocked } = useApiContract();
  const { drafts, patchDraft, results, setResult } = useWorkspace();
  const draft = drafts.phylogeny;
  const msa = asRecord(results.msa?.result);
  const { next } = useRequestGate();
  const [error, setError] = useState<unknown>(null);
  const [execution, setExecution] = useState("");
  const envelope = results.phylogeny;
  const result = asRecord(envelope?.result);
  const readyFasta = draft.fasta.trim() || String(msa.raw_alignment || "");

  async function infer() {
    const gate = next();
    setError(null);
    try {
      if (SYNC.has(draft.method)) {
        const { data } = await helixApi.phylogenyInfer(
          {
            fasta: readyFasta,
            method: draft.method,
            distance_model: draft.distance,
            alignment_hash: draft.alignmentHash || String(msa.alignment_hash || "") || null,
          },
          gate.signal,
        );
        if (gate.isCurrent()) setResult("phylogeny", data);
        return;
      }
      const { data: accepted } = await helixApi.phylogenySubmit(
        {
          fasta: readyFasta,
          method: draft.method,
          distance_model: draft.distance,
          alignment_hash: draft.alignmentHash || String(msa.alignment_hash || "") || null,
        },
        gate.signal,
      );
      setExecution(accepted.execution_status);
      const done = await pollJob(accepted, { signal: gate.signal });
      if (gate.isCurrent()) {
        setExecution(String(done.execution_status || ""));
        setResult("phylogeny", { ...done, result: jobScientificResult(done) });
      }
    } catch (err) {
      if (gate.isCurrent()) setError(err);
    }
  }

  const hasMsa = Boolean(readyFasta.trim()) || Boolean(msa.alignment_hash);

  return (
    <div className="hs-page">
      <p className="kicker">Evolution</p>
      <h1>Phylogeny</h1>
      <p data-testid="phy-source-hash">
        transferred {String(draft.alignmentHash || msa.alignment_hash || "")}
      </p>
      {!hasMsa ? (
        <EmptyState title="No completed MSA is loaded.">
          Run an MSA or load a valid prealigned FASTA before inferring a phylogenetic tree.
        </EmptyState>
      ) : null}
      <CompactInput compact={Boolean(envelope)} summary="Edit phylogeny input">
        <SequenceInput id="tree-fasta" label="Prealigned FASTA" value={draft.fasta} onChange={(fasta) => patchDraft("phylogeny", { fasta })} rows={10} />
        <div className="field">
          <label htmlFor="method">Method</label>
          <select id="method" value={draft.method} onChange={(e) => patchDraft("phylogeny", { method: e.target.value })}>
            <option value="neighbor_joining">NJ</option>
            <option value="upgma">UPGMA</option>
            <option value="iqtree_ml">IQ-TREE</option>
            <option value="fasttree_ml">FastTree</option>
          </select>
        </div>
        <button type="button" className="btn-primary" onClick={() => void infer()} disabled={scienceLocked || (!hasMsa && !draft.fasta)}>
          Infer tree
        </button>
        {execution ? <JobState execution={execution} /> : null}
      </CompactInput>
      {error ? <ScientificError error={error} /> : null}
      {envelope ? (
        <>
          <AnalysisNav
            items={[
              { id: "phy-summary", label: "Overview" },
              { id: "phy-tree", label: "Tree" },
              { id: "phy-meaning", label: "Interpretation" },
              { id: "phy-raw", label: "Raw data" },
            ]}
          />
          <ResultAnchor id="phy-summary">
          <p data-testid="phy-method">method {String(result.method)} · not taxonomic {String(result.not_taxonomic_tree)}</p>
          <p data-testid="phy-engine">
            engine {formatScientificValue(result.engine ?? result.software ?? result.tool)} · version{" "}
            {formatScientificValue(result.engine_version ?? result.version ?? asRecord(result.software).version)} ·
            model {formatScientificValue(result.model ?? result.substitution_model)}
          </p>
          <p>
            tree_hash {formatScientificValue(result.tree_hash)} · source MSA hash{" "}
            {formatScientificValue(result.source_msa_hash ?? result.alignment_hash)}
          </p>
          <p data-testid="newick" className="mono">
            {String(result.newick || "")}
          </p>
          <button
            type="button"
            className="btn-secondary"
            onClick={() => downloadText("helixscope-tree.nwk", `${String(result.newick || "").trim()}\n`)}
          >
            Download Newick
          </button>
          <p>source MSA hash {formatScientificValue(result.source_msa_hash ?? result.alignment_hash)}</p>
          </ResultAnchor>
          <ResultAnchor id="phy-tree">
          <InfographicFrame
            title="Inferred tree"
            dataSource="HelixScope Core phylogeny layout"
            method={String(result.method || "see provenance")}
            units="branch length as returned by the engine"
            status={result.status}
            shows="The inferred topology and branch lengths from Core"
            doesNotShow="A taxonomic classification, a rooted species tree, or fabricated support"
          >
            <TreeViewer result={result} />
          </InfographicFrame>
          </ResultAnchor>
          <ResultAnchor id="phy-meaning">
          <InterpretationPanel explanation={result.explanation} />
          </ResultAnchor>
          <ResultAnchor id="phy-raw">
          <RawDataPanel result={result} />
          </ResultAnchor>
          <ProvenancePanel envelope={envelope} extra={result} />
          <section className="hs-panel glass-medium">
            <h2>Optional NCBI taxonomy</h2>
            <p className="hs-lede" data-testid="taxonomy-not-phylogeny">
              TAXONOMIC ANNOTATION. Taxonomy colors or groups leaves from declared organism names. It
              does not change Newick or tree_hash and is not phylogenetic inference.
            </p>
            <div className="field">
              <label htmlFor="tax-email">NCBI email</label>
              <input id="tax-email" value={draft.taxonomyEmail} onChange={(e) => patchDraft("phylogeny", { taxonomyEmail: e.target.value })} />
            </div>
            <SequenceInput
              id="tax-orgs"
              label="Declared organisms (leaf_id=scientific name, one per line)"
              value={draft.declaredOrganisms}
              onChange={(declaredOrganisms) => patchDraft("phylogeny", { declaredOrganisms })}
              rows={4}
            />
            <button
              type="button"
              className="btn-secondary"
              data-testid="phy-taxonomy"
              onClick={async () => {
                setError(null);
                try {
                  const declared: Record<string, string> = {};
                  for (const line of draft.declaredOrganisms.split("\n")) {
                    const [id, ...rest] = line.split("=");
                    const name = rest.join("=").trim();
                    if (id.trim() && name) declared[id.trim()] = name;
                  }
                  const { data } = await helixApi.phylogenyTaxonomy({
                    tree: result,
                    email: draft.taxonomyEmail,
                    declared_organisms: declared,
                  });
                  setResult("phylogeny", data);
                } catch (err) {
                  setError(err);
                }
              }}
            >
              Attach NCBI taxonomy
            </button>
          </section>
        </>
      ) : null}
    </div>
  );
}
