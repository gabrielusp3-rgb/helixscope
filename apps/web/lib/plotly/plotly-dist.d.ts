declare module "plotly.js/dist/plotly.js" {
  const Plotly: {
    newPlot: (
      el: HTMLElement,
      data: unknown[],
      layout?: Record<string, unknown>,
      config?: Record<string, unknown>,
    ) => Promise<HTMLElement> | HTMLElement;
    purge: (el: HTMLElement) => void;
  };
  export default Plotly;
}
