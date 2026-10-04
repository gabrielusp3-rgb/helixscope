/** Shared HelixScope Plotly presentation. Does not alter scientific values. */

import { PLOT_AXIS, PLOT_FONT, PLOT_GRID, PLOT_INNER, PLOT_PAPER } from "@/lib/visual/nucleotides";

export const HELIX_TRACE = "#6893D0";
export const HELIX_TRACE_SOFT = "#91B4E4";
export const HELIX_TRACE_ICE = "#AFCBF1";
export const HELIX_PAPER = PLOT_PAPER;
export const HELIX_PLOT = PLOT_INNER;
export const HELIX_GRID = PLOT_GRID;
export const HELIX_FONT = PLOT_FONT;

export function helixLayout(options: {
  title?: string;
  xTitle?: string;
  yTitle?: string;
  height?: number;
  legend?: boolean;
  yType?: "linear" | "log";
  yRange?: [number, number];
} = {}): Record<string, unknown> {
  return {
    paper_bgcolor: HELIX_PAPER,
    plot_bgcolor: HELIX_PLOT,
    font: { color: HELIX_FONT, family: "IBM Plex Sans, sans-serif", size: 12 },
    margin: { t: options.title ? 40 : 18, r: 18, b: 52, l: 56 },
    title: options.title ? { text: options.title, font: { size: 13, color: "#D9E4FC" } } : undefined,
    xaxis: {
      title: options.xTitle,
      gridcolor: HELIX_GRID,
      zeroline: false,
      linecolor: PLOT_AXIS,
      tickfont: { family: "JetBrains Mono, monospace", size: 11 },
    },
    yaxis: {
      title: options.yTitle,
      type: options.yType ?? "linear",
      range: options.yRange,
      gridcolor: HELIX_GRID,
      zeroline: false,
      linecolor: PLOT_AXIS,
      tickfont: { family: "JetBrains Mono, monospace", size: 11 },
    },
    showlegend: Boolean(options.legend),
    height: options.height ?? 280,
    hoverlabel: {
      bgcolor: "rgba(2, 11, 22, 0.92)",
      bordercolor: "rgba(175, 203, 241, 0.35)",
      font: { color: "#D9E4FC" },
    },
  };
}

export function figureExportName(label: string): string {
  const slug = label
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_|_$/g, "")
    .slice(0, 56);
  return `helixscope_${slug || "figure"}`;
}

export function helixPlotConfig(filename = "helixscope-figure"): Record<string, unknown> {
  return {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
    scrollZoom: true,
    toImageButtonOptions: {
      format: "png",
      filename,
      height: 900,
      width: 1400,
      scale: 2,
    },
  };
}

type PlotlyExport = {
  Icons?: { camera?: unknown };
  downloadImage?: (gd: HTMLElement, opts: Record<string, unknown>) => Promise<unknown>;
};

export function helixPlotConfigWithSvg(Plotly: PlotlyExport, filename = "helixscope-figure"): Record<string, unknown> {
  const extra =
    Plotly.downloadImage && Plotly.Icons?.camera
      ? [
          {
            name: "downloadSvg",
            title: "Download SVG",
            icon: Plotly.Icons.camera,
            click: (gd: HTMLElement) => {
              void Plotly.downloadImage?.(gd, { format: "svg", filename, height: 900, width: 1400 });
            },
          },
        ]
      : [];
  return {
    ...helixPlotConfig(filename),
    modeBarButtonsToAdd: extra,
  };
}

export const HELIX_PLOT_CONFIG = helixPlotConfig();

export const HELIX_HEATMAP_SCALE: Array<[number, string]> = [
  [0, "#01040B"],
  [0.35, "#001F59"],
  [0.65, "#2355A0"],
  [1, "#91B4E4"],
];
