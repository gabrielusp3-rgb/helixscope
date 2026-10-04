"use client";

import { useEffect, useRef } from "react";
import { loadPlotly } from "@/lib/plotly/load";
import { figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  classes: unknown[];
  title?: string;
  testId?: string;
  height?: number;
};

const COLOR: Record<string, string> = {
  match: "#6893D0",
  mismatch: "#91B4E4",
  gap: "#093986",
};

export function ScienceColumnMap({ classes, title, testId, height = 120 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      const traces = (["match", "mismatch", "gap"] as const).map((category) => {
        const xs: number[] = [];
        const ys: number[] = [];
        classes.forEach((item, index) => {
          if (String(item) === category) {
            xs.push(index);
            ys.push(1);
          }
        });
        return {
          type: "bar",
          x: xs,
          y: ys,
          name: category,
          marker: { color: COLOR[category] },
          width: 1,
        };
      });
      await Plotly.newPlot(
        node,
        traces,
        {
          ...helixLayout({ title, xTitle: "Alignment column", yTitle: "", height, legend: true }),
          barmode: "stack",
          yaxis: { visible: false },
        },
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "column-map")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [classes, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-column-map"} style={{ height, width: "100%" }} />;
}
