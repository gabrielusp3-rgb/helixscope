"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { helixApi } from "@/lib/api/client";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { useWorkspace } from "@/lib/workspace/WorkstationProvider";
import { ScientificStatus } from "@/components/science/SciencePrimitives";
import { NAV_ITEMS } from "@/lib/nav";

export function OverviewWorkspace() {
  const { results, setResult, transfers } = useWorkspace();
  const [error, setError] = useState("");

  const [health, setHealth] = useState({ live: "", ready: "" });

  useEffect(() => {
    void helixApi
      .capabilities()
      .then(({ data }) => setResult("engines", data))
      .catch((err) => setError(err instanceof Error ? err.message : "API unavailable"));
    void helixApi.healthLive().then(({ data }) => setHealth((prev) => ({ ...prev, live: String(asRecord(data).status || "LIVE") }))).catch(() => setHealth((prev) => ({ ...prev, live: "UNREACHABLE" })));
    void helixApi.healthReady().then(({ data }) => setHealth((prev) => ({ ...prev, ready: String(asRecord(data).status || "READY") }))).catch(() => setHealth((prev) => ({ ...prev, ready: "UNREACHABLE" })));
  }, [setResult]);

  const envelope = results.engines;
  const payload = asRecord(envelope?.result);
  const tools = asRecord(payload.tools);
  const refs = asList(payload.references).map((row) => asRecord(row));

  return (
    <div className="hs-page">
      <p className="kicker">Workspace</p>
      <h1>HelixScope</h1>
      <p className="hs-lede">
        Independent scientific workstation. Biology is calculated by HelixScope Core through FastAPI v1.
        This interface presents results, provenance, and camera/layout interaction only.
      </p>
      <div className="chip-row">
        <span>API {envelope?.api_version || "v1"}</span>
        <span>product {envelope?.product_version || "from backend"}</span>
        <span data-testid="health-live">live {health.live || "pending"}</span>
        <span data-testid="health-ready">ready {health.ready || "pending"}</span>
      </div>
      {error ? <p className="error-box">{error}</p> : null}
      <section className="hs-panel science-surface">
        <h2>Engine status summary</h2>
        <div className="chip-row">
          {Object.values(tools)
            .slice(0, 12)
            .map((tool) => {
              const row = asRecord(tool);
              return (
                <span key={String(row.name)}>
                  {String(row.name)} <ScientificStatus status={row.status} />
                </span>
              );
            })}
        </div>
        <p className="hs-topbar-meta">{String(payload.note || "")}</p>
      </section>
      <section className="hs-panel science-surface">
        <h2>Reference snapshot</h2>
        {refs.map((row) => (
          <p key={String(row.assembly_id)}>
            {String(row.assembly_id)} <ScientificStatus status={row.status} />
          </p>
        ))}
      </section>
      <section className="hs-panel glass-medium">
        <h2>Modules</h2>
        <div className="metric-grid">
          {NAV_ITEMS.filter((item) => item.href !== "/overview").map((item) => (
            <Link key={item.href} href={item.href} className="metric-card">
              <div className="label">{item.group}</div>
              <div className="value">{item.label}</div>
            </Link>
          ))}
        </div>
      </section>
      <section className="hs-panel science-surface">
        <h2>Session transfers</h2>
        {transfers.length === 0 ? (
          <p className="hs-topbar-meta">No explicit object transfers yet. Refresh clears in-memory drafts.</p>
        ) : (
          transfers.map((item) => (
            <p key={item.at}>
              {item.kind} · {item.note} · {formatScientificValue(item.payload.alignment_hash || item.payload.sequence?.slice(0, 12))}
            </p>
          ))
        )}
      </section>
    </div>
  );
}
