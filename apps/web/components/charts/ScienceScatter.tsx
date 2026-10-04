"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  x: unknown[];
  y: unknown[];
  text?: unknown[];
  xTitle?: string;
  yTitle?: string;
  title?: string;
  testId?: string;
  height?: number;
};

export function ScienceScatter({ x, y, text, xTitle, yTitle, title, testId, height = 280 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      await Plotly.newPlot(
        node,
        [
          {
            type: "scatter",
            mode: "markers",
            x: x.map((v, i) => chartNumeric(v) ?? i),
            y: y.map((v) => chartNumeric(v)),
            text: text?.map((item) => String(item)),
            marker: { color: HELIX_TRACE, size: 8 },
          },
        ],
        helixLayout({ title, xTitle, yTitle, height }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "scatter")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [x, y, text, xTitle, yTitle, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-scatter"} style={{ height, width: "100%" }} />;
}
