"use client";

import dynamic from "next/dynamic";
import { useCallback, useState } from "react";
import { asRecord, formatScientificValue } from "@/lib/api/numeric";
import { ScientificStatus } from "@/components/science/SciencePrimitives";
import { StructureViewportPlotly } from "@/components/structure/StructureViewportPlotly";
import { BUNDLED_MOLSTAR_IDS, evidenceBanner, type StructureScene } from "@/components/structure/structureScene";
import type { Orbit } from "@/lib/gesture/math";

const MolstarViewport = dynamic(
  () => import("@/components/structure/MolstarViewport").then((m) => m.MolstarViewport),
  { ssr: false },
);

export type { StructureScene };

export function StructureViewport({
  scene,
  selectedIndex,
  onSelect,
  orbit,
  onOrbit,
}: {
  scene: StructureScene | Record<string, unknown> | null;
  selectedIndex?: number | null;
  onSelect?: (index: number | null) => void;
  orbit?: Orbit;
  onOrbit?: (orbit: Orbit) => void;
}) {
  const rec = asRecord(scene);
  const id = String(rec.structure_id || "").toUpperCase();
  const kind = String(rec.kind || rec.status || "");
  const illustrative = kind.toLowerCase().includes("illustrative") || String(rec.status).includes("ILLUSTRATIVE");
  const molstarEligible = BUNDLED_MOLSTAR_IDS.has(id) && !illustrative;
  const [engine, setEngine] = useState<"molstar" | "plotly">(molstarEligible ? "molstar" : "plotly");
  const [plugin, setPlugin] = useState<{ dispose: () => void; managers?: { camera?: { reset?: () => void } }; helpers?: { viewportScreenshot?: { getImage?: () => Promise<{ dataUri?: string }> } } } | null>(null);
  const [full, setFull] = useState(false);

  const onReady = useCallback((next: typeof plugin) => setPlugin(next), []);

  async function screenshot() {
    const image = await plugin?.helpers?.viewportScreenshot?.getImage?.();
    if (!image?.dataUri) return;
    const anchor = document.createElement("a");
    anchor.href = image.dataUri;
    anchor.download = `helixscope_${id || "structure"}.png`;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
  }

  return (
    <div className={`viewport-3d ${full ? "viewport-3d-full" : ""}`} data-testid="structure-viewport">
      <div className="floating-control glass-thin structure-banner">
        <div className="chip-row">
          <ScientificStatus status={rec.kind || rec.status || rec.mapping_status} />
          <span data-testid="structure-evidence-banner">{evidenceBanner(rec)}</span>
        </div>
        <p className="hs-topbar-meta">
          source {String(rec.source ?? "API")} · mapping {formatScientificValue(rec.mapping_status)}
        </p>
        <div className="hs-actions">
          {molstarEligible ? (
            <button type="button" className="btn-ghost" onClick={() => setEngine(engine === "molstar" ? "plotly" : "molstar")}>
              {engine === "molstar" ? "Legacy Plotly atoms" : "Mol* cartoon"}
            </button>
          ) : null}
          <button type="button" className="btn-ghost" onClick={() => plugin?.managers?.camera?.reset?.()}>
            Reset camera
          </button>
          <button type="button" className="btn-ghost" onClick={() => setFull((value) => !value)}>
            {full ? "Exit fullscreen" : "Fullscreen"}
          </button>
          <button type="button" className="btn-ghost" onClick={() => void screenshot()}>
            Screenshot
          </button>
        </div>
        {engine === "molstar" ? (
          <p className="hs-topbar-meta">
            Mol* 5.11.0 · Hand Control is PARTIAL on Mol* (mouse/touch camera is authoritative).
          </p>
        ) : (
          <p className="hs-topbar-meta">Legacy Plotly atom scatter remains available for illustrative helices and fallback.</p>
        )}
      </div>
      {engine === "molstar" && molstarEligible ? (
        <MolstarViewport structureId={id} onReady={onReady} />
      ) : (
        <StructureViewportPlotly
          scene={scene}
          selectedIndex={selectedIndex}
          onSelect={onSelect}
          orbit={orbit}
          onOrbit={onOrbit}
        />
      )}
    </div>
  );
}
