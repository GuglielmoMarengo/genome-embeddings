import { Activity, ArrowRight, Layers3, TrendingUp } from 'lucide-react'
import type { AnalysisData } from '../types'
import { ChartCard } from '../components/ChartCard'
import { Plot } from '../components/Plot'
import { Badge, Card, SectionHeading } from '../components/ui'

function orderedTrajectory(values: Record<string, number>) {
  return Object.entries(values)
    .map(([k, value]) => [Number(k), value] as const)
    .sort((first, second) => first[0] - second[0])
}

export function MultiscalePage({ analysis }: { analysis: AnalysisData }) {
  const legacy = orderedTrajectory(analysis.pair_trajectories.legacy_euclidean)
  const v2 = orderedTrajectory(analysis.pair_trajectories.descriptor_v2_euclidean)
  const js = orderedTrajectory(analysis.pair_trajectories.jensen_shannon)
  const stepRows = Object.entries(analysis.jensen_shannon_step_distances)
  const stabilityRows = Object.entries(analysis.jensen_shannon_ranking_stability)

  const pairData = [
    {
      type: 'scatter', mode: 'lines+markers', name: 'Legacy Euclidean',
      x: legacy.map(([k]) => k), y: legacy.map(([, value]) => value),
      line: { color: '#38bdf8', width: 3 }, marker: { size: 9 },
      hovertemplate: 'k=%{x}<br>Legacy=%{y:.6f}<extra></extra>',
    },
    {
      type: 'scatter', mode: 'lines+markers', name: 'Descriptor V2 Euclidean',
      x: v2.map(([k]) => k), y: v2.map(([, value]) => value),
      line: { color: '#a78bfa', width: 3 }, marker: { size: 9 },
      hovertemplate: 'k=%{x}<br>Descriptor V2=%{y:.6f}<extra></extra>',
    },
    {
      type: 'scatter', mode: 'lines+markers', name: 'Jensen–Shannon',
      x: js.map(([k]) => k), y: js.map(([, value]) => value),
      line: { color: '#2dd4bf', width: 3 }, marker: { size: 9 },
      hovertemplate: 'k=%{x}<br>Jensen–Shannon=%{y:.6f}<extra></extra>',
    },
  ]

  const stabilityData = [
    {
      type: 'scatter', mode: 'lines+markers', name: 'Kendall tau',
      x: stabilityRows.map(([transition]) => transition),
      y: stabilityRows.map(([, row]) => row.kendall_tau),
      line: { color: '#38bdf8', width: 3 }, marker: { size: 10 },
      fill: 'tozeroy', fillcolor: 'rgba(56,189,248,.10)',
      hovertemplate: '%{x}<br>Kendall τ=%{y:.3f}<extra></extra>',
    },
  ]

  const stepData = [
    {
      type: 'bar',
      x: stepRows.map(([transition]) => transition),
      y: stepRows.map(([, value]) => value),
      marker: {
        color: stepRows.map(([, value]) => value),
        colorscale: [[0, '#0f766e'], [0.5, '#0ea5e9'], [1, '#8b5cf6']],
      },
      hovertemplate: '%{x}<br>Matrix step=%{y:.6f}<extra></extra>',
    },
  ]

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Cross-scale behavior"
        title="Multiscale trajectories"
        description="Follow one selected pair, global Jensen–Shannon deformation, and reference-based rank stability as k increases."
      />

      <div className="grid gap-4 sm:grid-cols-3">
        <Card className="p-5">
          <span className="grid h-10 w-10 place-items-center rounded-xl bg-sky-500/10 text-sky-500"><Layers3 size={19} /></span>
          <p className="mt-4 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Scales analyzed</p>
          <p className="mt-2 text-3xl font-semibold">{analysis.summary.k_values.length}</p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">k={analysis.summary.k_values.join(', ')}</p>
        </Card>
        <Card className="p-5">
          <span className="grid h-10 w-10 place-items-center rounded-xl bg-teal-500/10 text-teal-500"><TrendingUp size={19} /></span>
          <p className="mt-4 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Largest JS step</p>
          <p className="mt-2 text-3xl font-semibold">{Math.max(...stepRows.map(([, value]) => value)).toFixed(4)}</p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">Global matrix deformation</p>
        </Card>
        <Card className="p-5">
          <span className="grid h-10 w-10 place-items-center rounded-xl bg-violet-500/10 text-violet-500"><Activity size={19} /></span>
          <p className="mt-4 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Mean Kendall tau</p>
          <p className="mt-2 text-3xl font-semibold">
            {(stabilityRows.reduce((total, [, row]) => total + row.kendall_tau, 0) / stabilityRows.length).toFixed(3)}
          </p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">Ordinal stability across transitions</p>
        </Card>
      </div>

      <ChartCard
        title={`${analysis.summary.reference_label} vs ${analysis.summary.comparison_label}`}
        description="The metrics have different scales; compare trajectory shape rather than absolute height alone."
        toolbar={<Badge tone="info">Selected pair</Badge>}
      >
        <Plot
          data={pairData}
          height={500}
          layout={{
            xaxis: { title: 'k-mer length', dtick: 1 },
            yaxis: { title: 'Pair value' },
            legend: { orientation: 'h', y: 1.15 },
          }}
        />
      </ChartCard>

      <div className="grid gap-6 xl:grid-cols-2">
        <ChartCard title="Jensen–Shannon ranking stability" description="Kendall correlation between consecutive reference rankings">
          <Plot
            data={stabilityData}
            layout={{
              yaxis: { title: 'Kendall tau', range: [-1.05, 1.05] },
              xaxis: { title: 'Scale transition' },
              shapes: [
                { type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 1, y1: 1, line: { color: '#2dd4bf', dash: 'dot', width: 1 } },
                { type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: '#94a3b8', dash: 'dot', width: 1 } },
              ],
            }}
          />
        </ChartCard>
        <ChartCard title="Jensen–Shannon matrix step distances" description="Euclidean movement of the full upper-triangle vector">
          <Plot data={stepData} layout={{ yaxis: { title: 'Step distance' }, xaxis: { title: 'Scale transition' } }} />
        </ChartCard>
      </div>

      <Card className="overflow-hidden">
        <div className="border-b border-[var(--border)] px-5 py-5 sm:px-6">
          <h3 className="font-semibold">Jensen–Shannon ranking flow</h3>
          <p className="mt-1 text-xs text-[var(--text-muted)]">Reference: {analysis.summary.reference_label}</p>
        </div>
        <div className="space-y-3 p-4 sm:p-5">
          {Object.entries(analysis.jensen_shannon_rankings)
            .sort(([first], [second]) => Number(first) - Number(second))
            .map(([k, ranking]) => (
              <div key={k} className="grid gap-3 rounded-2xl border border-[var(--border)] bg-[var(--bg-solid)]/55 p-4 lg:grid-cols-[4rem_1fr] lg:items-center">
                <Badge tone="info">k={k}</Badge>
                <div className="flex flex-wrap items-center gap-2">
                  {ranking.map(([label, value], index) => (
                    <div key={label} className="flex items-center gap-2">
                      <span className="rounded-xl bg-[var(--bg-subtle)] px-3 py-2 text-xs font-semibold text-[var(--text-soft)]">
                        {index + 1}. {label} <span className="ml-1 font-mono text-[var(--text-muted)]">{value.toFixed(4)}</span>
                      </span>
                      {index < ranking.length - 1 ? <ArrowRight size={14} className="text-[var(--text-muted)]" /> : null}
                    </div>
                  ))}
                </div>
              </div>
            ))}
        </div>
      </Card>

      <Card className="overflow-hidden">
        <div className="border-b border-[var(--border)] px-5 py-4 sm:px-6">
          <h3 className="font-semibold">Transition diagnostics</h3>
        </div>
        <div className="overflow-x-auto table-scroll">
          <table className="w-full min-w-[48rem] border-collapse text-left">
            <thead><tr className="border-b border-[var(--border)] text-xs uppercase tracking-[0.13em] text-[var(--text-muted)]">
              {['Transition', 'Kendall tau', 'Inversions', 'Mean shift', 'Max shift', 'Exact match'].map((heading) => <th key={heading} className="px-5 py-4 font-semibold first:sm:px-6 last:sm:px-6">{heading}</th>)}
            </tr></thead>
            <tbody>
              {stabilityRows.map(([transition, row]) => (
                <tr key={transition} className="border-b border-[var(--border)] last:border-0 hover:bg-[var(--surface-hover)]">
                  <td className="px-5 py-4 font-semibold sm:px-6">{transition}</td>
                  <td className="px-5 py-4 font-mono">{row.kendall_tau.toFixed(3)}</td>
                  <td className="px-5 py-4 font-mono">{row.discordant_pairs}</td>
                  <td className="px-5 py-4 font-mono">{row.mean_absolute_rank_shift.toFixed(3)}</td>
                  <td className="px-5 py-4 font-mono">{row.max_rank_shift}</td>
                  <td className="px-5 py-4 sm:px-6"><Badge tone={row.exact_match ? 'success' : 'warning'}>{row.exact_match ? 'Yes' : 'No'}</Badge></td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </Card>
    </div>
  )
}
