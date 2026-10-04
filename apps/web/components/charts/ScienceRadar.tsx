"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_PLOT, HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  labels: unknown[];
  values: unknown[];
  title?: string;
  testId?: string;
  height?: number;
};

export function ScienceRadar({ labels, values, title, testId, height = 320 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const theta = labels.map((label) => String(label));
      const r = values.map((value) => chartNumeric(value));
      await Plotly.newPlot(
        node,
        [
          {
            type: "scatterpolar",
            r,
            theta,
            fill: "toself",
            line: { color: HELIX_TRACE },
          },
        ],
        {
          ...helixLayout({ title, height }),
          polar: {
            bgcolor: HELIX_PLOT,
            radialaxis: { gridcolor: "rgba(145, 180, 228, 0.12)", tickfont: { size: 10 } },
            angularaxis: { gridcolor: "rgba(145, 180, 228, 0.12)", tickfont: { size: 10 } },
          },
        },
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "radar")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [labels, values, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-radar"} style={{ height, width: "100%" }} />;
}
