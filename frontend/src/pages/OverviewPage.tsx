import { motion } from 'motion/react'
import {
  ArrowRight,
  Binary,
  Dna,
  Gauge,
  GitCompareArrows,
  Layers3,
  Network,
  Sparkles,
} from 'lucide-react'
import type { AnalysisData } from '../types'
import { DnaVisual } from '../components/DnaVisual'
import { Badge, Button, Card, MetricCard, SectionHeading } from '../components/ui'

function format(value: number, digits = 5): string {
  return value.toFixed(digits)
}

export function OverviewPage({
  analysis,
  openSetup,
  goToMatrices,
}: {
  analysis: AnalysisData
  openSetup: () => void
  goToMatrices: () => void
}) {
  const { summary } = analysis
  const comparisonRows = [
    ['Legacy Euclidean', summary.legacy_euclidean_distance, 'Baseline descriptor geometry'],
    ['Legacy cosine', summary.legacy_cosine_similarity, 'Directional similarity'],
    ['Descriptor V2 Euclidean', summary.descriptor_v2_euclidean_distance, 'Finite-sample + dependencies'],
    ['Descriptor V2 cosine', summary.descriptor_v2_cosine_similarity, 'Directional V2 comparison'],
    ['Embedding V2 Euclidean', summary.embedding_v2_euclidean_distance, 'Global + multiscale blocks'],
    ['Jensen–Shannon', summary.jensen_shannon_distance, 'Complete k-mer distributions'],
  ] as const

  const biologicalRows = analysis.dataset.filter((row) => !row.label.toLowerCase().includes('periodic'))
  const averageLength = biologicalRows.reduce((total, row) => total + row.length, 0) / biologicalRows.length
  const stableTransitions = Object.values(analysis.jensen_shannon_ranking_stability).filter(
    (row) => row.kendall_tau >= 0.8,
  ).length

  return (
    <div className="space-y-7">
      <section className="grid gap-6 xl:grid-cols-[1.15fr_.85fr]">
        <Card className="relative overflow-hidden p-6 sm:p-8 lg:p-10">
          <div className="absolute right-0 top-0 h-64 w-64 rounded-full bg-sky-500/10 blur-3xl" />
          <motion.div
            initial={{ opacity: 0, y: 14 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.55 }}
            className="relative z-10"
          >
            <div className="flex flex-wrap items-center gap-2">
              <Badge tone="info">Descriptor Foundation V2</Badge>
              <Badge tone="success">Analysis complete</Badge>
              <Badge>k={summary.k_values.join(', ')}</Badge>
            </div>
            <h1 className="mt-7 max-w-4xl text-4xl font-semibold leading-[1.04] tracking-[-0.045em] text-[var(--text)] sm:text-5xl lg:text-6xl">
              Explore genomic structure as an{' '}
              <span className="text-gradient">interpretable mathematical space.</span>
            </h1>
            <p className="mt-6 max-w-2xl text-base leading-7 text-[var(--text-soft)] sm:text-lg">
              Compare composition, finite-sample complexity, dependency signatures, multiscale embeddings,
              and complete k-mer distributions without hiding the mathematics behind a black box.
            </p>

            <div className="mt-8 flex flex-wrap gap-3">
              <Button icon={Sparkles} onClick={openSetup}>Configure analysis</Button>
              <Button variant="secondary" icon={Network} onClick={goToMatrices}>
                Explore matrices
              </Button>
            </div>

            <div className="mt-9 flex flex-wrap items-center gap-x-7 gap-y-3 border-t border-[var(--border)] pt-6 text-sm text-[var(--text-muted)]">
              <span className="flex items-center gap-2"><Dna size={16} className="text-sky-500" />{summary.dataset_size} sequences</span>
              <span className="flex items-center gap-2"><Layers3 size={16} className="text-teal-500" />{summary.k_values.length} scales</span>
              <span className="flex items-center gap-2"><GitCompareArrows size={16} className="text-violet-500" />{summary.reference_label} → {summary.comparison_label}</span>
            </div>
          </motion.div>
        </Card>
        <DnaVisual />
      </section>

      <section className="grid gap-4 sm:grid-cols-2 xl:grid-cols-4">
        <MetricCard
          label="Sequences"
          value={summary.dataset_size}
          detail={`${Math.round(averageLength).toLocaleString()} nt mean biological length`}
          icon={Dna}
          accent="blue"
        />
        <MetricCard
          label="Descriptor V2 distance"
          value={format(summary.descriptor_v2_euclidean_distance)}
          detail={`Primary scale k=${summary.selected_k}`}
          icon={Gauge}
          accent="teal"
        />
        <MetricCard
          label="Embedding V2 distance"
          value={format(summary.embedding_v2_euclidean_distance)}
          detail="Global and scale-specific blocks"
          icon={Layers3}
          accent="violet"
        />
        <MetricCard
          label="Jensen–Shannon"
          value={format(summary.jensen_shannon_distance)}
          detail={`${stableTransitions} stable scale transitions`}
          icon={Binary}
          accent="amber"
        />
      </section>

      <section className="grid gap-6 xl:grid-cols-[1.2fr_.8fr]">
        <Card className="overflow-hidden">
          <div className="border-b border-[var(--border)] px-5 py-5 sm:px-6">
            <SectionHeading
              eyebrow="Selected pair"
              title="Representation comparison"
              description="The same sequence pair viewed through six complementary mathematical representations."
            />
          </div>
          <div className="overflow-x-auto table-scroll">
            <table className="w-full min-w-[46rem] border-collapse text-left">
              <thead>
                <tr className="border-b border-[var(--border)] text-xs uppercase tracking-[0.14em] text-[var(--text-muted)]">
                  <th className="px-5 py-4 font-semibold sm:px-6">Representation</th>
                  <th className="px-5 py-4 font-semibold">Interpretation</th>
                  <th className="px-5 py-4 text-right font-semibold sm:px-6">Pair value</th>
                </tr>
              </thead>
              <tbody>
                {comparisonRows.map(([label, value, detail], index) => (
                  <motion.tr
                    key={label}
                    initial={{ opacity: 0, x: -8 }}
                    animate={{ opacity: 1, x: 0 }}
                    transition={{ delay: index * 0.045 }}
                    className="border-b border-[var(--border)] last:border-0 hover:bg-[var(--surface-hover)]"
                  >
                    <td className="px-5 py-4 sm:px-6">
                      <div className="flex items-center gap-3">
                        <span className="grid h-9 w-9 place-items-center rounded-xl bg-[var(--bg-subtle)] text-[var(--primary)]">
                          {index + 1}
                        </span>
                        <span className="font-semibold text-[var(--text)]">{label}</span>
                      </div>
                    </td>
                    <td className="px-5 py-4 text-sm text-[var(--text-soft)]">{detail}</td>
                    <td className="px-5 py-4 text-right font-mono text-sm font-semibold text-[var(--text)] sm:px-6">
                      {format(value, 6)}
                    </td>
                  </motion.tr>
                ))}
              </tbody>
            </table>
          </div>
        </Card>

        <div className="space-y-4">
          <Card className="p-6">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-teal-500/10 text-teal-500">
                <Binary size={20} />
              </span>
              <div>
                <p className="text-xs font-semibold uppercase tracking-[0.16em] text-teal-500">Signal expansion</p>
                <h3 className="mt-2 text-lg font-semibold">Jensen–Shannon avoids cosine compression</h3>
                <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                  Complete k-mer probability profiles preserve identity information that entropy and diversity summaries necessarily discard.
                </p>
              </div>
            </div>
          </Card>
          <Card className="p-6">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-violet-500/10 text-violet-500">
                <Layers3 size={20} />
              </span>
              <div>
                <p className="text-xs font-semibold uppercase tracking-[0.16em] text-violet-500">Interpretable embedding</p>
                <h3 className="mt-2 text-lg font-semibold">Global features are represented once</h3>
                <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                  The V2 embedding separates global composition from scale-specific diagnostics, preventing repeated weighting of invariant features.
                </p>
              </div>
            </div>
          </Card>
          <button
            type="button"
            onClick={goToMatrices}
            className="focus-ring group flex w-full items-center justify-between rounded-[1.4rem] border border-[var(--border-strong)] bg-[var(--surface-hover)] p-5 text-left transition hover:border-[var(--primary)]"
          >
            <div>
              <p className="font-semibold text-[var(--text)]">Inspect pairwise geometry</p>
              <p className="mt-1 text-sm text-[var(--text-muted)]">Switch metrics and explore interactive heatmaps.</p>
            </div>
            <ArrowRight className="text-[var(--primary)] transition-transform group-hover:translate-x-1" size={20} />
          </button>
        </div>
      </section>
    </div>
  )
}
