/// <reference types="vite/client" />

declare module 'plotly.js-dist-min' {
  interface PlotlyApi {
    react: (
      element: HTMLElement,
      data: unknown[],
      layout?: Record<string, unknown>,
      config?: Record<string, unknown>,
    ) => Promise<void>
    purge: (element: HTMLElement) => void
    Plots: { resize: (element: HTMLElement) => void }
  }
  const Plotly: PlotlyApi
  export default Plotly
}
