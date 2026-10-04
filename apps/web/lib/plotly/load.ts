type PlotlyModule = {
  newPlot: (
    el: HTMLElement,
    data: unknown[],
    layout?: Record<string, unknown>,
    config?: Record<string, unknown>,
  ) => Promise<HTMLElement> | HTMLElement;
  purge: (el: HTMLElement) => void;
  downloadImage?: (el: HTMLElement, opts: Record<string, unknown>) => Promise<unknown>;
  Icons?: { camera?: unknown };
};

export async function loadPlotly(): Promise<PlotlyModule> {
  const mod = (await import("plotly.js/dist/plotly.js")) as { default?: PlotlyModule } & PlotlyModule;
  return mod.default ?? mod;
}
