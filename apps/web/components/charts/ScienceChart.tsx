"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  x: unknown[];
  y: unknown[];
  xTitle?: string;
  yTitle?: string;
  title?: string;
  testId?: string;
  height?: number;
};

export function ScienceChart({ x, y, xTitle, yTitle, title, testId, height = 280 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const xs = x.map((v, i) => (typeof v === "number" ? v : i));
      const ys = y.map((v) => chartNumeric(v));
      await Plotly.newPlot(
        node,
        [
          {
            type: "scatter",
            mode: "lines",
            x: xs,
            y: ys,
            connectgaps: false,
            line: { color: HELIX_TRACE, width: 1.6 },
          },
        ],
        helixLayout({ title, xTitle, yTitle, height }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "chart")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [x, y, xTitle, yTitle, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-chart"} style={{ height, width: "100%" }} />;
}
