import { Atom, CircleDot, Gauge, Repeat2 } from 'lucide-react'
import type { AnalysisData } from '../types'
import { ChartCard } from '../components/ChartCard'
import { Plot } from '../components/Plot'
import { Badge, Card, SectionHeading } from '../components/ui'

export function DescriptorsPage({ analysis }: { analysis: AnalysisData }) {
  const rows = analysis.descriptor_v2
  const labels = rows.map((row) => row.label)

  const entropyData = [
    {
      type: 'bar',
      x: labels,
      y: rows.map((row) => row.conditional_entropy),
      name: 'Conditional nucleotide entropy',
      marker: { color: '#38bdf8', opacity: 0.88 },
      hovertemplate: '%{x}<br>Conditional H=%{y:.4f} bits<extra></extra>',
    },
    {
      type: 'bar',
      x: labels,
      y: rows.map((row) => row.finite_sample_entropy),
      name: 'Finite-sample normalized k-mer entropy',
      marker: { color: '#2dd4bf', opacity: 0.82 },
      hovertemplate: '%{x}<br>Finite-sample H=%{y:.4f}<extra></extra>',
    },
  ]

  const coverageData = [
    {
      type: 'scatter',
      mode: 'lines+markers',
      x: labels,
      y: rows.map((row) => row.theoretical_coverage),
      name: 'Theoretical space',
      line: { color: '#a78bfa', width: 3 },
      marker: { size: 8 },
      hovertemplate: '%{x}<br>Theoretical coverage=%{y:.4f}<extra></extra>',
    },
    {
      type: 'scatter',
      mode: 'lines+markers',
      x: labels,
      y: rows.map((row) => row.observable_coverage),
      name: 'Observable space',
      line: { color: '#38bdf8', width: 3 },
      marker: { size: 8 },
      hovertemplate: '%{x}<br>Observable coverage=%{y:.4f}<extra></extra>',
    },
  ]

  const effectiveData = [
    {
      type: 'bar',
      x: labels,
      y: rows.map((row) => row.effective_kmer_count),
      marker: {
        color: rows.map((row) => row.effective_kmer_count),
        colorscale: [[0, '#0f766e'], [0.55, '#0ea5e9'], [1, '#8b5cf6']],
      },
      hovertemplate: '%{x}<br>Effective k-mers=%{y:.2f}<extra></extra>',
    },
  ]

  const dinucleotides = Object.keys(rows[0]?.dinucleotide_odds_ratios ?? {})
  const dinucleotideData = [
    {
      type: 'heatmap',
      x: dinucleotides,
      y: labels,
      z: rows.map((row) => dinucleotides.map((key) => row.dinucleotide_odds_ratios[key])),
      colorscale: [[0, '#312e81'], [0.5, '#f8fafc'], [1, '#0f766e']],
      zmid: 1,
      colorbar: { thickness: 13, len: 0.8, title: { text: 'Observed / expected', side: 'right' } },
      hovertemplate: '%{y}<br>%{x}: %{z:.4f}<extra></extra>',
    },
  ]

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Descriptor Foundation V2"
        title="Finite-sample and dependency diagnostics"
        description="Inspect complexity, coverage, repeat structure, and local dependency without confusing sequence length with biological signal."
      />

      <div className="grid gap-4 md:grid-cols-2 xl:grid-cols-4">
        <Card className="p-5">
          <div className="flex items-center justify-between">
            <span className="grid h-10 w-10 place-items-center rounded-xl bg-sky-500/10 text-sky-500"><Atom size={19} /></span>
            <Badge tone="info">bits</Badge>
          </div>
          <p className="mt-5 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Conditional entropy</p>
          <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">Uncertainty of the next nucleotide given the current nucleotide.</p>
        </Card>
        <Card className="p-5">
          <div className="flex items-center justify-between">
            <span className="grid h-10 w-10 place-items-center rounded-xl bg-teal-500/10 text-teal-500"><Gauge size={19} /></span>
            <Badge tone="success">normalized</Badge>
          </div>
          <p className="mt-5 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Finite-sample entropy</p>
          <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">Entropy divided by the maximum that the available windows can actually observe.</p>
        </Card>
        <Card className="p-5">
          <div className="flex items-center justify-between">
            <span className="grid h-10 w-10 place-items-center rounded-xl bg-violet-500/10 text-violet-500"><CircleDot size={19} /></span>
            <Badge>2ᴴ</Badge>
          </div>
          <p className="mt-5 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Effective k-mers</p>
          <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">The diversity-equivalent number of equally probable k-mer categories.</p>
        </Card>
        <Card className="p-5">
          <div className="flex items-center justify-between">
            <span className="grid h-10 w-10 place-items-center rounded-xl bg-amber-500/10 text-amber-500"><Repeat2 size={19} /></span>
            <Badge tone="warning">sparsity</Badge>
          </div>
          <p className="mt-5 text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Repeat diagnostics</p>
          <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">Singleton and repeated-window fractions expose high-k saturation and periodicity.</p>
        </Card>
      </div>

      <div className="grid gap-6 xl:grid-cols-2">
        <ChartCard title="Entropy diagnostics" description={`Primary scale k=${analysis.summary.selected_k}`}>
          <Plot
            data={entropyData}
            layout={{
              barmode: 'group',
              yaxis: { title: 'Entropy / normalized entropy' },
              xaxis: { tickangle: -22 },
              legend: { orientation: 'h', y: 1.16 },
            }}
          />
        </ChartCard>
        <ChartCard title="k-mer space coverage" description="Theoretical versus finite observable support">
          <Plot
            data={coverageData}
            layout={{
              yaxis: { title: 'Coverage', range: [0, 1.05] },
              xaxis: { tickangle: -22 },
              legend: { orientation: 'h', y: 1.16 },
            }}
          />
        </ChartCard>
      </div>

      <div className="grid gap-6 xl:grid-cols-[.8fr_1.2fr]">
        <ChartCard title="Effective k-mer count" description="Entropy converted into an intuitive diversity-equivalent count">
          <Plot data={effectiveData} layout={{ yaxis: { title: 'Effective categories' }, xaxis: { tickangle: -22 } }} />
        </ChartCard>
        <ChartCard title="Dinucleotide dependency signature" description="Observed/expected enrichment ratios; values above 1 indicate enrichment">
          <Plot
            data={dinucleotideData}
            layout={{
              margin: { l: 125, r: 55, t: 20, b: 55 },
              xaxis: { title: 'Dinucleotide' },
              yaxis: { autorange: 'reversed' },
            }}
          />
        </ChartCard>
      </div>

      <Card className="overflow-hidden">
        <div className="border-b border-[var(--border)] px-5 py-5 sm:px-6">
          <h3 className="font-semibold">Descriptor table</h3>
          <p className="mt-1 text-xs text-[var(--text-muted)]">Values are generated at the selected single-scale k.</p>
        </div>
        <div className="overflow-x-auto table-scroll">
          <table className="w-full min-w-[88rem] border-collapse text-left">
            <thead>
              <tr className="border-b border-[var(--border)] text-xs uppercase tracking-[0.12em] text-[var(--text-muted)]">
                {['Sequence', 'Length', 'GC', 'Conditional H', 'Finite-sample H', 'Effective k-mers', 'Theoretical coverage', 'Observable coverage', 'Singleton', 'Repeated'].map((heading) => (
                  <th key={heading} className="px-5 py-4 font-semibold first:sm:px-6 last:sm:px-6">{heading}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr key={row.label} className="border-b border-[var(--border)] last:border-0 hover:bg-[var(--surface-hover)]">
                  <td className="px-5 py-4 font-semibold sm:px-6">{row.label}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.length}</td>
                  <td className="px-5 py-4 font-mono text-sm">{(row.gc_content * 100).toFixed(2)}%</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.conditional_entropy.toFixed(4)}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.finite_sample_entropy.toFixed(4)}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.effective_kmer_count.toFixed(2)}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.theoretical_coverage.toFixed(4)}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.observable_coverage.toFixed(4)}</td>
                  <td className="px-5 py-4 font-mono text-sm">{row.singleton_fraction.toFixed(4)}</td>
                  <td className="px-5 py-4 font-mono text-sm sm:px-6">{row.repeated_window_fraction.toFixed(4)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </Card>
    </div>
  )
}
