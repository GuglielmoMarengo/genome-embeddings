import { useEffect, useMemo, useRef } from 'react'
import Plotly from 'plotly.js-dist-min'
import { useThemeContext } from './ThemeProvider'

interface PlotProps {
  data: unknown[]
  layout?: Record<string, unknown>
  className?: string
  height?: number
}

export function Plot({ data, layout = {}, className = '', height = 440 }: PlotProps) {
  const container = useRef<HTMLDivElement>(null)
  const { theme } = useThemeContext()

  const themedLayout = useMemo(() => {
    const dark = theme === 'dark'
    return {
      autosize: true,
      height,
      margin: { l: 64, r: 28, t: 58, b: 58 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: {
        family: 'Inter, ui-sans-serif, system-ui, sans-serif',
        color: dark ? '#b8c8d8' : '#415267',
        size: 12,
      },
      xaxis: {
        gridcolor: dark ? 'rgba(148,190,222,.10)' : 'rgba(31,61,90,.09)',
        zerolinecolor: dark ? 'rgba(148,190,222,.14)' : 'rgba(31,61,90,.12)',
        linecolor: 'transparent',
        ...((layout.xaxis as object | undefined) ?? {}),
      },
      yaxis: {
        gridcolor: dark ? 'rgba(148,190,222,.10)' : 'rgba(31,61,90,.09)',
        zerolinecolor: dark ? 'rgba(148,190,222,.14)' : 'rgba(31,61,90,.12)',
        linecolor: 'transparent',
        ...((layout.yaxis as object | undefined) ?? {}),
      },
      hoverlabel: {
        bgcolor: dark ? '#102238' : '#ffffff',
        bordercolor: dark ? 'rgba(148,190,222,.18)' : 'rgba(31,61,90,.14)',
        font: { color: dark ? '#eff8ff' : '#0d1b2a' },
      },
      ...layout,
    }
  }, [height, layout, theme])

  useEffect(() => {
    const element = container.current
    if (!element) return

    void Plotly.react(element, data, themedLayout, {
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      modeBarButtonsToRemove: ['lasso2d', 'select2d'],
      toImageButtonOptions: {
        format: 'png',
        filename: 'genome_embeddings_chart',
        scale: 2,
      },
    })

    const observer = new ResizeObserver(() => Plotly.Plots.resize(element))
    observer.observe(element)
    return () => {
      observer.disconnect()
      Plotly.purge(element)
    }
  }, [data, themedLayout])

  return <div ref={container} className={`plot-container w-full ${className}`} />
}
