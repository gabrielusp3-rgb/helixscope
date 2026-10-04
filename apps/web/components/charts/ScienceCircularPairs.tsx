"use client";

import { useEffect, useRef } from "react";
import { asList, asRecord, chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  layout: unknown;
  title?: string;
  testId?: string;
  height?: number;
};

/** Plot Core circular_layout coordinates. Does not infer base pairs. */
export function ScienceCircularPairs({ layout, title, testId, height = 360 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const rec = asRecord(layout);
      const xs = asList(rec.x).map((value) => chartNumeric(value));
      const ys = asList(rec.y).map((value) => chartNumeric(value));
      const bases = asList(rec.bases).map((value) => String(value));
      const traces: unknown[] = [];
      asList(rec.chords).forEach((chord) => {
        const item = asRecord(chord);
        traces.push({
          type: "scatter",
          mode: "lines",
          x: asList(item.x).map((value) => chartNumeric(value)),
          y: asList(item.y).map((value) => chartNumeric(value)),
          line: { color: "rgba(104, 147, 208, 0.75)", width: 1.4 },
          hovertemplate: `${item.base_i ?? ""}${item.i ?? ""}–${item.base_j ?? ""}${item.j ?? ""}<extra></extra>`,
          showlegend: false,
        });
      });
      traces.push({
        type: "scatter",
        mode: "markers+text",
        x: xs,
        y: ys,
        text: bases,
        textposition: "top center",
        marker: { size: 8, color: HELIX_TRACE },
        hovertext: bases.map((base, index) => `pos ${index} (0-based)<br>${base}`),
        hoverinfo: "text",
        showlegend: false,
      });
      await Plotly.newPlot(
        node,
        traces,
        {
          ...helixLayout({
            title: title ?? "Predicted pairing topology (circular). Not RNAPlot and not 3D.",
            height,
          }),
          xaxis: { visible: false, scaleanchor: "y", constrain: "domain" },
          yaxis: { visible: false },
        },
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "circular-pairs")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [layout, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-circular-pairs"} style={{ height, width: "100%" }} />;
}
