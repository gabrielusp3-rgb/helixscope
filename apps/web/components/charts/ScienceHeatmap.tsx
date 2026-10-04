"use client";

import { useEffect, useRef } from "react";
import { chartNumeric } from "@/lib/api/numeric";
import { loadPlotly } from "@/lib/plotly/load";
import { HELIX_HEATMAP_SCALE, figureExportName, helixLayout, helixPlotConfigWithSvg } from "@/lib/plotly/theme";

type Props = {
  matrix?: unknown;
  z?: Array<Array<number | null>>;
  x?: string[];
  y?: string[];
  zmid?: number;
  title?: string;
  testId?: string;
  height?: number;
};

export function ScienceHeatmap({ matrix, z, x, y, zmid, title, testId, height = 320 }: Props) {
  const host = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let disposed = false;
    const node = host.current;
    if (!node) return;
    void (async () => {
      const Plotly = await loadPlotly();
      if (disposed) return;
      let values = z;
      if (!values) {
        const packed = matrix && typeof matrix === "object" ? (matrix as { data?: unknown }) : null;
        const raw = Array.isArray(packed?.data) ? packed.data : Array.isArray(matrix) ? matrix : [];
        values = (raw as unknown[]).map((row) =>
          Array.isArray(row) ? row.map((cell) => chartNumeric(cell)) : [],
        );
      }
      await Plotly.newPlot(
        node,
        [
          {
            type: "heatmap",
            z: values,
            x,
            y,
            zmid,
            colorscale: HELIX_HEATMAP_SCALE,
            showscale: true,
            hoverongaps: false,
          },
        ],
        helixLayout({ title, height }),
        helixPlotConfigWithSvg(Plotly, figureExportName(title || "heatmap")),
      );
    })();
    return () => {
      disposed = true;
      void loadPlotly().then((Plotly) => {
        Plotly.purge(node);
      });
    };
  }, [matrix, z, x, y, zmid, title, height]);

  return <div ref={host} className="hs-chart" data-testid={testId ?? "science-heatmap"} style={{ height, width: "100%" }} />;
}
