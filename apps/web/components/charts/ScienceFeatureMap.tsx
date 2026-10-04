"use client";

import { useEffect, useRef } from "react";
import { asList, asRecord, chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, HELIX_TRACE_SOFT, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  length: number;
  tracks: unknown;
  title?: string;
  testId?: string;
  height?: number;
};

const COLORS = [HELIX_TRACE, HELIX_TRACE_SOFT, "#2355A0", "#91B4E4"];

export function ScienceFeatureMap({ length, tracks, title, testId, height = 180 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const traces: unknown[] = [];
      asList(tracks).forEach((track, index) => {
        const rec = asRecord(track);
        const name = String(rec.name || `track ${index + 1}`);
        const color = COLORS[index % COLORS.length];
        asList(rec.features).forEach((feature) => {
          const item = asRecord(feature);
          const start = chartNumeric(item.start);
          const end = chartNumeric(item.end);
          if (start === null || end === null || end <= start) return;
          traces.push({
            type: "scatter",
            mode: "lines",
            x: [start, end],
            y: [name, name],
            line: { color, width: 14 },
            hovertemplate: `${String(item.label || name)}<br>start ${start}<br>end ${end}<extra></extra>`,
            showlegend: false,
          });
        });
      });
      await Plotly.newPlot(
        node,
        traces,
        {
          ...helixLayout({
            title: title ?? "Sequence map (coordinates from computed features)",
            xTitle: "Position on input strand (0-based)",
            yTitle: "Track",
            height,
          }),
          xaxis: {
            ...(helixLayout({}).xaxis as object),
            range: [0, Math.max(length, 1)],
            title: "Position on input strand (0-based)",
          },
        },
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "feature-map")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [length, tracks, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-feature-map"} style={{ height, width: "100%" }} />;
}
