import { Binary, Grid3X3, Layers3, Maximize2, Network } from 'lucide-react'
import { useMemo, useState } from 'react'
import type { AnalysisData, MatrixData } from '../types'
import { ChartCard } from '../components/ChartCard'
import { Plot } from '../components/Plot'
import { Badge, Card, SectionHeading, cn } from '../components/ui'

const matrixMeta: Record<string, { label: string; description: string; family: string }> = {
  legacy_euclidean: {
    label: 'Legacy Euclidean',
    description: 'Six normalized aggregate descriptors',
    family: 'Legacy',
  },
  legacy_cosine: {
    label: 'Legacy cosine',
    description: 'Directional similarity in the baseline descriptor space',
    family: 'Legacy',
  },
  descriptor_v2_euclidean: {
    label: 'Descriptor V2 Euclidean',
    description: 'Finite-sample, dependency, and sparsity-aware descriptor geometry',
    family: 'Descriptor V2',
  },
  descriptor_v2_cosine: {
    label: 'Descriptor V2 cosine',
    description: 'Directional similarity in the expanded descriptor space',
    family: 'Descriptor V2',
  },
  embedding_v2_euclidean: {
    label: 'Multiscale embedding V2',
    description: 'Global features once plus scale-specific diagnostic blocks',
    family: 'Embedding V2',
  },
  jensen_shannon: {
    label: 'Jensen–Shannon distance',
    description: 'Complete normalized k-mer probability distributions',
    family: 'Distribution',
  },
}

function heatmap(matrix: MatrixData) {
  const isSimilarity = matrix.metric.includes('cosine')
  return [
    {
      type: 'heatmap',
      z: matrix.values,
      x: matrix.labels,
      y: matrix.labels,
      colorscale: isSimilarity
        ? [[0, '#312e81'], [0.5, '#0284c7'], [1, '#f0fdfa']]
        : [[0, '#06101d'], [0.36, '#0f766e'], [0.72, '#0ea5e9'], [1, '#fbbf24']],
      reversescale: false,
      colorbar: { thickness: 13, len: 0.82, title: { text: isSimilarity ? 'Similarity' : 'Distance', side: 'right' } },
      hovertemplate: '%{y} vs %{x}<br>Value=%{z:.6f}<extra></extra>',
    },
  ]
}

export function MatricesPage({ analysis }: { analysis: AnalysisData }) {
  const keys = Object.keys(analysis.matrices).filter((key) => matrixMeta[key])
  const [selected, setSelected] = useState(keys.includes('jensen_shannon') ? 'jensen_shannon' : keys[0])
  const matrix = analysis.matrices[selected]
  const meta = matrixMeta[selected]

  const ranking = useMemo(() => {
    const reference = analysis.summary.reference_label
    const index = matrix.labels.indexOf(reference)
    const values = matrix.labels
      .map((label, candidate) => ({ label, value: matrix.values[index][candidate] }))
      .filter((row) => row.label !== reference)
    return values.sort((first, second) =>
      matrix.metric.includes('cosine') ? second.value - first.value : first.value - second.value,
    )
  }, [analysis.summary.reference_label, matrix])

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Pairwise geometry"
        title="Comparison matrices"
        description="Switch representations, inspect exact values, zoom into the heatmap, and follow the selected reference ranking."
      />

      <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-3">
        {keys.map((key) => {
          const item = matrixMeta[key]
          const active = key === selected
          const Icon = key === 'jensen_shannon' ? Binary : key.includes('embedding') ? Layers3 : key.includes('v2') ? Network : Grid3X3
          return (
            <button
              key={key}
              type="button"
              onClick={() => setSelected(key)}
              className={cn(
                'focus-ring glass-panel flex items-start gap-3 rounded-[1.25rem] p-4 text-left transition-all duration-200',
                active ? 'border-[var(--primary)] bg-[var(--surface-hover)] shadow-lg shadow-sky-500/5' : 'hover:border-[var(--border-strong)]',
              )}
            >
              <span className={cn('grid h-10 w-10 shrink-0 place-items-center rounded-xl', active ? 'bg-sky-500/15 text-sky-500' : 'bg-[var(--bg-subtle)] text-[var(--text-muted)]')}>
                <Icon size={18} />
              </span>
              <span>
                <span className="flex flex-wrap items-center gap-2">
                  <span className="text-sm font-semibold text-[var(--text)]">{item.label}</span>
                  {active ? <Badge tone="info">Active</Badge> : null}
                </span>
                <span className="mt-1 block text-xs leading-5 text-[var(--text-muted)]">{item.description}</span>
              </span>
            </button>
          )
        })}
      </div>

      <div className="grid gap-6 2xl:grid-cols-[1fr_22rem]">
        <ChartCard
          title={`${meta.label} · k=${matrix.kmer_length}`}
          description={meta.description}
          toolbar={<Badge>{meta.family}</Badge>}
        >
          <Plot
            data={heatmap(matrix)}
            height={640}
            layout={{
              margin: { l: 125, r: 55, t: 24, b: 105 },
              xaxis: { tickangle: -35, side: 'bottom' },
              yaxis: { autorange: 'reversed' },
            }}
          />
        </ChartCard>

        <Card className="h-fit overflow-hidden">
          <div className="border-b border-[var(--border)] px-5 py-5">
            <div className="flex items-center justify-between gap-3">
              <div>
                <p className="text-xs font-semibold uppercase tracking-[0.15em] text-[var(--primary)]">Reference ranking</p>
                <h3 className="mt-1 font-semibold">{analysis.summary.reference_label}</h3>
              </div>
              <span className="grid h-10 w-10 place-items-center rounded-xl bg-violet-500/10 text-violet-500"><Maximize2 size={18} /></span>
            </div>
          </div>
          <div className="p-3">
            {ranking.map((row, index) => (
              <div key={row.label} className="flex items-center gap-3 rounded-xl px-3 py-3 hover:bg-[var(--surface-hover)]">
                <span className="grid h-8 w-8 place-items-center rounded-lg bg-[var(--bg-subtle)] text-xs font-semibold text-[var(--text-soft)]">{index + 1}</span>
                <div className="min-w-0 flex-1">
                  <p className="truncate text-sm font-semibold">{row.label}</p>
                  <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-[var(--bg-subtle)]">
                    <div
                      className="h-full rounded-full bg-[linear-gradient(90deg,var(--primary),var(--secondary))]"
                      style={{ width: `${Math.max(4, matrix.metric.includes('cosine') ? row.value * 100 : (1 - Math.min(1, row.value / Math.max(...ranking.map((item) => item.value)))) * 100)}%` }}
                    />
                  </div>
                </div>
                <span className="font-mono text-xs font-semibold text-[var(--text-soft)]">{row.value.toFixed(5)}</span>
              </div>
            ))}
          </div>
        </Card>
      </div>

      <Card className="overflow-hidden">
        <div className="border-b border-[var(--border)] px-5 py-4">
          <h3 className="font-semibold">Numerical matrix</h3>
          <p className="mt-1 text-xs text-[var(--text-muted)]">Exact values used by the heatmap and ranking.</p>
        </div>
        <div className="overflow-x-auto table-scroll p-3">
          <table className="w-full min-w-[60rem] border-separate border-spacing-1 text-center text-xs">
            <thead>
              <tr>
                <th className="p-2 text-left text-[var(--text-muted)]">Sequence</th>
                {matrix.labels.map((label) => <th key={label} className="max-w-28 p-2 font-semibold text-[var(--text-soft)]">{label}</th>)}
              </tr>
            </thead>
            <tbody>
              {matrix.values.map((values, rowIndex) => (
                <tr key={matrix.labels[rowIndex]}>
                  <th className="p-2 text-left font-semibold">{matrix.labels[rowIndex]}</th>
                  {values.map((value, columnIndex) => (
                    <td key={matrix.labels[columnIndex]} className="rounded-lg bg-[var(--bg-subtle)] p-2 font-mono text-[var(--text-soft)]">
                      {value.toFixed(4)}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </Card>
    </div>
  )
}
