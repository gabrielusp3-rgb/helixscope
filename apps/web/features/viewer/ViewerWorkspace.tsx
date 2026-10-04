"use client";

import dynamic from "next/dynamic";
import { useCallback, useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord } from "@/lib/api/numeric";
import { useApiContract } from "@/lib/api/contractIdentity";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { ScientificError, ScientificStatus, ScientificWarning } from "@/components/science/SciencePrimitives";
import { ProvenancePanel } from "@/components/science/Provenance";
import type { Orbit } from "@/lib/gesture/math";

const StructureViewport = dynamic(
  () => import("@/components/structure/StructureViewport").then((m) => m.StructureViewport),
  { ssr: false },
);
const HandControl = dynamic(
  () => import("@/components/structure/HandControl").then((m) => m.HandControl),
  { ssr: false },
);

const DEFAULT_ORBIT: Orbit = { yaw: 0.7, pitch: 0.35, radius: 2.4 };

export function ViewerWorkspace() {
  const { scienceLocked } = useApiContract();
  const { patchDraft, results, setResult } = useWorkspace();
  const [error, setError] = useState<unknown>(null);
  const [orbit, setOrbit] = useState<Orbit>(DEFAULT_ORBIT);
  const [selected, setSelected] = useState<number | null>(null);
  const [mappingStatus, setMappingStatus] = useState<string | null>(null);
  const [catalog, setCatalog] = useState<Record<string, unknown> | null>(null);
  const envelope = results.viewer;
  const result = asRecord(envelope?.result);
  const viewer = asRecord(result.viewer);

  useEffect(() => {
    void helixApi.structureCatalog().then(({ data }) => setCatalog(asRecord(data.result))).catch(() => null);
  }, []);

  const load = useCallback(async (id: "1CRN" | "1BNA" | "1RNA") => {
    setError(null);
    try {
      const { data } = await helixApi.structureFixture({ structure_id: id });
      setResult("viewer", data);
      patchDraft("viewer", { structureId: id });
    } catch (err) {
      setError(err);
    }
  }, [patchDraft, setResult]);

  const classify4un3 = useCallback(async () => {
    setError(null);
    try {
      const { data } = await helixApi.structureClassifyMapping({
        sequences_identical: true,
        polymer_index_mode: false,
        observed_coordinate_residues: 0,
        polymer_length: 0,
      });
      setMappingStatus(String(asRecord(data.result).mapping_status || ""));
    } catch (err) {
      setError(err);
    }
  }, []);

  return (
    <div className="hs-page">
      <p className="kicker">Structure</p>
      <h1>3D Viewer</h1>
      <p className="hs-lede">
        Coordinates come from the API. Camera, zoom, orbit and selection presentation are owned here.
        Source coordinates are not mutated.
      </p>
      <div className="hs-actions">
        {(["1BNA", "1RNA", "1CRN"] as const).map((id) => (
          <button key={id} type="button" className="btn-primary" disabled={scienceLocked} onClick={() => void load(id)}>
            Load {id}
          </button>
        ))}
        <button type="button" className="btn-secondary" disabled={scienceLocked} onClick={() => void classify4un3()} data-testid="classify-4un3">
          Classify 4UN3-class mapping
        </button>
        <button
          type="button"
          className="btn-secondary"
          data-testid="fetch-4un3"
          disabled={scienceLocked}
          onClick={() => {
            setError(null);
            void helixApi
              .structureExperimentalComplex({ structure_id: "4UN3" })
              .then(({ data }) => setResult("viewer", data))
              .catch(setError);
          }}
        >
          Fetch experimental 4UN3
        </button>
        <button
          type="button"
          className="btn-secondary"
          data-testid="load-mapped-1crn"
          onClick={() => {
            setError(null);
            void helixApi
              .structureFixtureMapped({
                structure_id: "1CRN",
                query_sequence: "TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN",
                molecule: "PROTEIN",
              })
              .then(({ data }) => {
                setResult("viewer", data);
                patchDraft("viewer", { structureId: "1CRN" });
              })
              .catch(setError);
          }}
        >
          Load mapped 1CRN
        </button>
      </div>
      {error ? <ScientificError error={error} /> : null}
      <ScientificWarning>
        4UN3 remains a catalog pointer. Mapping with polymer_length 0 stays UNCERTAIN. Coordinates are not fabricated.
      </ScientificWarning>
      {catalog ? (
        <section className="hs-panel science-surface">
          <h2>Experimental complex pointers</h2>
          {asList(catalog.complexes).map((item, index) => {
            const row = asRecord(item);
            return (
              <p key={index}>
                {String(row.pdb_id || row.structure_id || JSON.stringify(row).slice(0, 80))}
              </p>
            );
          })}
        </section>
      ) : null}
      {mappingStatus ? (
        <p data-testid="mapping-status">
          mapping_status <ScientificStatus status={mappingStatus} /> {mappingStatus}
        </p>
      ) : null}
      {viewer.coordinates || result.kind ? (
        <div className="split">
          <HandControl orbit={orbit} onOrbit={setOrbit} onReset={() => setOrbit(DEFAULT_ORBIT)} />
          <StructureViewport
            scene={{ ...viewer, ...result }}
            selectedIndex={selected}
            onSelect={setSelected}
            orbit={orbit}
            onOrbit={setOrbit}
          />
        </div>
      ) : null}
      {envelope ? <ProvenancePanel envelope={envelope} extra={result} /> : null}
    </div>
  );
}
