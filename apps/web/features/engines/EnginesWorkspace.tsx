"use client";

import { useEffect, useState } from "react";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { ScientificError, ScientificStatus } from "@/components/science/SciencePrimitives";
import { DataTable } from "@/components/science/DataTable";

export function EnginesWorkspace() {
  const { results, setResult } = useWorkspace();
  const [error, setError] = useState<unknown>(null);
  const envelope = results.engines;

  useEffect(() => {
    void helixApi.capabilities().then(({ data }) => setResult("engines", data)).catch(setError);
  }, [setResult]);

  const payload = asRecord(envelope?.result);
  const tools = asRecord(payload.tools);
  const providers = asList(payload.remote_providers).map((row) => asRecord(row));

  return (
    <div className="hs-page">
      <p className="kicker">Workspace</p>
      <h1>Scientific Engines</h1>
      <p className="hs-lede">{String(payload.note || "")}</p>
      {error ? <ScientificError error={error} /> : null}
      <div className="hs-table-wrap">
        <table className="hs-table">
          <thead>
            <tr>
              <th>Name</th>
              <th>Status</th>
              <th>Backend</th>
              <th>Version</th>
              <th>Purpose</th>
              <th>Validation</th>
              <th>Limitation</th>
            </tr>
          </thead>
          <tbody>
            {Object.values(tools).map((tool) => {
              const row = asRecord(tool);
              const live = asRecord(row.live);
              return (
                <tr key={String(row.name)}>
                  <td>{String(row.name)}</td>
                  <td>
                    <ScientificStatus status={row.status} />
                  </td>
                  <td>{String(row.source ?? row.local_or_remote ?? live.source)}</td>
                  <td>{formatScientificValue(row.version || live.version)}</td>
                  <td>{String(row.purpose || row.scientific_purpose || asList(row.capabilities)[0] || "")}</td>
                  <td>{formatScientificValue(live.status ?? row.last_validated)}</td>
                  <td>{String(row.reason || row.limitation || "")}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      <section className="hs-panel science-surface" data-testid="remote-providers">
        <h2>Remote scientific providers</h2>
        <p className="hs-lede">
          Status CONFIGURED means the client exists. REMOTE_VALIDATED is only used after a controlled
          live interaction. Secrets and API keys are not displayed.
        </p>
        <DataTable
          rows={providers.map((row) => ({
            provider: row.provider,
            service: row.service,
            local_or_remote: row.local_or_remote,
            status: row.status,
            api_key_configured: row.api_key_configured,
            limitation: row.limitation,
            privacy: row.privacy,
          }))}
          caption="Remote providers"
          testId="remote-providers-table"
        />
      </section>
      <details className="disclosure hs-panel science-surface">
        <summary>Advanced details</summary>
        <DataTable
          rows={Object.values(tools).map((tool) => {
            const row = asRecord(tool);
            return {
              name: row.name,
              live: asRecord(row.live).status,
              capabilities: asList(row.capabilities).join(", "),
            };
          })}
        />
      </details>
    </div>
  );
}
