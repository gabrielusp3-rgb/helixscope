"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { asList, asRecord, formatScientificValue } from "@/lib/api/numeric";
import { tracesXyzFingerprint, type Orbit, orbitToEye } from "@/lib/gesture/math";
import { loadPlotly } from "@/lib/plotly/load";
import { ScientificStatus } from "@/components/science/SciencePrimitives";
import type { StructureScene } from "@/components/structure/structureScene";

const DEFAULT_ORBIT: Orbit = { yaw: 0.7, pitch: 0.35, radius: 2.4 };

export function StructureViewportPlotly({
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
  const host = useRef<HTMLDivElement>(null);
  const [localOrbit, setLocalOrbit] = useState<Orbit>(DEFAULT_ORBIT);
  const usedOrbit = orbit ?? localOrbit;
  const rec = asRecord(scene);
  const coords = asList(rec.coordinates ?? rec.atoms);
  const traces = useMemo(() => {
    const xs: number[] = [];
    const ys: number[] = [];
    const zs: number[] = [];
    const text: string[] = [];
    coords.forEach((item, index) => {
      const atom = asRecord(item);
      const x = Number(atom.x);
      const y = Number(atom.y);
      const z = Number(atom.z);
      if (![x, y, z].every(Number.isFinite)) return;
      xs.push(x);
      ys.push(y);
      zs.push(z);
      text.push(
        `${String(atom.residue ?? atom.comp_id ?? "")} ${String(atom.auth_seq_id ?? atom.label_seq_id ?? index)} ${String(atom.atom ?? atom.atom_name ?? "")}`,
      );
    });
    return [{ x: xs, y: ys, z: zs, text, type: "scatter3d" as const, mode: "markers" as const }];
  }, [coords]);

  const fingerprint = useMemo(() => tracesXyzFingerprint(traces), [traces]);

  useEffect(() => {
    let disposed = false;
    const el = host.current;
    if (!el) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed || !host.current) return;
      const plotEl = host.current;
      const eye = orbitToEye(usedOrbit.yaw, usedOrbit.pitch, usedOrbit.radius);
      await Plotly.newPlot(
        plotEl,
        [
          {
            ...traces[0],
            marker: {
              size: 3,
              color: traces[0].x.map((_, i) => (i === selectedIndex ? "#D9E4FC" : "#6893D0")),
            },
            hoverinfo: "text",
          },
        ],
        {
          paper_bgcolor: "#01040B",
          scene: {
            bgcolor: "#01040B",
            xaxis: { visible: false },
            yaxis: { visible: false },
            zaxis: { visible: false },
            camera: { eye },
          },
          margin: { t: 8, r: 8, b: 8, l: 8 },
        },
        { displayModeBar: false, responsive: true },
      );
      const graph = plotEl as unknown as HTMLElement & {
        on: (event: string, cb: (payload: Record<string, unknown>) => void) => void;
      };
      graph.on("plotly_click", (event) => {
        const points = event.points as Array<{ pointNumber?: number }> | undefined;
        const idx = points?.[0]?.pointNumber;
        onSelect?.(typeof idx === "number" ? idx : null);
      });
      graph.on("plotly_relayout", (event) => {
        const ex = event["scene.camera.eye.x"];
        const ey = event["scene.camera.eye.y"];
        const ez = event["scene.camera.eye.z"];
        if ([ex, ey, ez].every((v) => typeof v === "number")) {
          const next = {
            yaw: Math.atan2(Number(ey), Number(ex)),
            pitch: Math.asin(
              Math.max(-1, Math.min(1, Number(ez) / Math.hypot(Number(ex), Number(ey), Number(ez) || 1))),
            ),
            radius: Math.hypot(Number(ex), Number(ey), Number(ez)),
          };
          if (onOrbit) onOrbit(next);
          else setLocalOrbit(next);
        }
      });
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => Plotly.purge(el));
    };
  }, [traces, usedOrbit, selectedIndex, onOrbit, onSelect]);

  return (
    <>
      <div className="hs-topbar-meta" data-testid="coord-fingerprint">
        coord fingerprint n={fingerprint.count} · viewer Plotly 3.3.0 (legacy atom scatter)
      </div>
      <div className="hs-topbar-meta">
        atoms {formatScientificValue(rec.n_atoms ?? coords.length)} · hash {String(rec.structure_hash ?? "").slice(0, 12)}
      </div>
      <ScientificStatus status={rec.kind || rec.status || rec.mapping_status} />
      <div ref={host} className="plot-host" />
    </>
  );
}
