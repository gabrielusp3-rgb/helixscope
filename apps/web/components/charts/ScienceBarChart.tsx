"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_TRACE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  labels: unknown[];
  values: unknown[];
  xTitle?: string;
  yTitle?: string;
  title?: string;
  testId?: string;
  height?: number;
  orientation?: "v" | "h";
  logY?: boolean;
  yRange?: [number, number];
};

export function ScienceBarChart({
  labels,
  values,
  xTitle,
  yTitle,
  title,
  testId,
  height = 280,
  orientation = "v",
  logY = false,
  yRange,
}: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const numeric = values.map((value) => chartNumeric(value));
      const names = labels.map((label) => String(label));
      const horizontal = orientation === "h";
      await Plotly.newPlot(
        node,
        [
          {
            type: "bar",
            orientation: horizontal ? "h" : "v",
            x: horizontal ? numeric : names,
            y: horizontal ? names : numeric,
            marker: { color: HELIX_TRACE },
          },
        ],
        helixLayout({
          title,
          xTitle,
          yTitle,
          height,
          yType: !horizontal && logY ? "log" : "linear",
          yRange: horizontal ? undefined : yRange,
        }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "bar")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [labels, values, xTitle, yTitle, title, height, orientation, logY, yRange]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-bar"} style={{ height, width: "100%" }} />;
}
