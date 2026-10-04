"use client";

import { useEffect, useRef } from "react";
import { asList, asRecord, chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  length: number;
  pairs: unknown;
  title?: string;
  testId?: string;
  height?: number;
};

/** Draw MFE pair arcs from API base_pairs. Does not infer pairing. */
export function ScienceArcDiagram({ length, pairs, title, testId, height = 240 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const traces: unknown[] = [
        {
          type: "scatter",
          mode: "lines",
          x: [0, Math.max(length - 1, 0)],
          y: [0, 0],
          line: { color: "rgba(175, 203, 241, 0.35)", width: 2 },
          showlegend: false,
          hoverinfo: "skip",
        },
      ];
      asList(pairs).forEach((row) => {
        const rec = asRecord(row);
        const left = chartNumeric(rec.position_0based ?? rec.i ?? rec.left);
        const right = chartNumeric(rec.paired_position_0based ?? rec.j ?? rec.right);
        if (left === null || right === null || right <= left) return;
        const mid = (left + right) / 2;
        const radius = (right - left) / 2;
        const xs: number[] = [];
        const ys: number[] = [];
        for (let step = 0; step <= 24; step += 1) {
          const angle = Math.PI * (step / 24);
          xs.push(mid + radius * Math.cos(Math.PI - angle));
          ys.push(radius * Math.sin(angle));
        }
        traces.push({
          type: "scatter",
          mode: "lines",
          x: xs,
          y: ys,
          line: { color: HELIX_TRACE, width: 1.2 },
          showlegend: false,
          hovertemplate: `${left}–${right}<extra></extra>`,
        });
      });
      await Plotly.newPlot(
        node,
        traces,
        helixLayout({
          title: title ?? "Predicted base-pair arcs (MFE pairs from API)",
          xTitle: "Position (0-based)",
          yTitle: "Arc height (presentation)",
          height,
        }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "arcs")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [length, pairs, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-arcs"} style={{ height, width: "100%" }} />;
}
