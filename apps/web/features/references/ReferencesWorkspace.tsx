"use client";

import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { ScientificError, ScientificStatus, ScientificWarning } from "@/components/science/SciencePrimitives";

export function ReferencesWorkspace() {
  const { results, setResult } = useWorkspace();
  const [error, setError] = useState<unknown>(null);
  const [selected, setSelected] = useState("GRCh38.p14");
  const envelope = results.references;
  const admission = asRecord(results.referencesAdmission?.result);

  useEffect(() => {
    void helixApi.references().then(({ data }) => setResult("references", data)).catch(setError);
  }, [setResult]);

  const rows = asList(asRecord(envelope?.result).references).map((row) => asRecord(row));
  const decision = String(admission.decision || admission.status || "");
  const canInstall = decision === "PROCEED" && admission.not_a_public_assembly !== true;

  return (
    <div className="hs-page">
      <p className="kicker">Data</p>
      <h1>References</h1>
      <p className="hs-lede">
        Assembly readiness is whatever the API reports. This browser does not infer GRCh38 or T2T status.
      </p>
      <ScientificWarning>
        Install is catalogued assemblies only. Resource admission runs first. GRCh38 is not downloaded when
        admission returns RESOURCE_LIMIT. TEST REFERENCE is not a public genome.
      </ScientificWarning>
      {error ? <ScientificError error={error} /> : null}
      <div className="field">
        <label htmlFor="ref-id">Catalog assembly</label>
        <select id="ref-id" data-testid="references-select" value={selected} onChange={(e) => setSelected(e.target.value)}>
          {rows.map((row) => (
            <option key={String(row.assembly_id)} value={String(row.assembly_id)}>
              {String(row.assembly_id)}
            </option>
          ))}
        </select>
      </div>
      <button
        type="button"
        className="btn-secondary"
        data-testid="references-admit"
        onClick={async () => {
          setError(null);
          try {
            const { data } = await helixApi.referencesAdmission({ assembly_id: selected });
            setResult("referencesAdmission", data);
          } catch (err) {
            setError(err);
          }
        }}
      >
        Check resource admission
      </button>
      {results.referencesAdmission ? (
        <p data-testid="references-admission">
          decision {decision} · <ScientificStatus status={decision || admission.status} />
          {decision === "RESOURCE_LIMIT" ? ` · ${String(admission.reason || admission.message || "resource admission refused")}` : ""}
        </p>
      ) : null}
      <button
        type="button"
        className="btn-primary"
        data-testid="references-install"
        disabled={!canInstall}
        onClick={async () => {
          setError(null);
          try {
            const { data } = await helixApi.referencesInstall({ assembly_id: selected, confirm: true });
            setResult("referencesAdmission", data);
          } catch (err) {
            setError(err);
          }
        }}
      >
        Install after PROCEED
      </button>
      <div className="hs-table-wrap">
        <table className="hs-table" data-testid="references-table">
          <thead>
            <tr>
              <th>Assembly</th>
              <th>Status</th>
              <th>Kind</th>
              <th>Ready</th>
              <th>Version</th>
              <th>Checksum</th>
              <th>Size</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((row) => (
              <tr key={String(row.assembly_id)}>
                <td>{String(row.assembly_id)}</td>
                <td>
                  <ScientificStatus status={row.status} />
                </td>
                <td>{String(row.kind)}</td>
                <td data-testid={`ref-ready-${String(row.assembly_id)}`}>
                  {row.ready === true || String(row.status).includes("READY") ? "READY" : "NOT READY"}
                </td>
                <td>{formatScientificValue(row.version ?? row.assembly_version)}</td>
                <td className="mono">{formatScientificValue(row.checksum ?? row.sha256)}</td>
                <td>{formatScientificValue(row.size_bytes ?? row.size)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
