"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  values: unknown[];
  xTitle?: string;
  yTitle?: string;
  title?: string;
  testId?: string;
  height?: number;
};

export function ScienceHistogram({ values, xTitle, yTitle, title, testId, height = 260 }: Props) {
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
            type: "histogram",
            x: values.map((value) => chartNumeric(value)).filter((value): value is number => value !== null),
            marker: { color: HELIX_TRACE },
          },
        ],
        helixLayout({ title, xTitle, yTitle, height }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "histogram")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [values, xTitle, yTitle, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-histogram"} style={{ height, width: "100%" }} />;
}
