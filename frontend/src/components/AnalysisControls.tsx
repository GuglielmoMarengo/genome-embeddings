import { AnimatePresence, motion } from 'motion/react'
import { ChevronRight, Dna, FlaskConical, RotateCcw, Settings2, X } from 'lucide-react'
import type { AnalysisConfig, SequenceRecord } from '../types'
import { Button, Card, cn } from './ui'
import { FileDropzone } from './FileDropzone'

const allK = [1, 2, 3, 4, 5, 6, 7, 8]

export function AnalysisControls({
  open,
  close,
  records,
  config,
  setConfig,
  onAnalyze,
  onReset,
  onFiles,
  onRemove,
  busy,
}: {
  open: boolean
  close: () => void
  records: SequenceRecord[]
  config: AnalysisConfig
  setConfig: (config: AnalysisConfig) => void
  onAnalyze: () => Promise<void>
  onReset: () => Promise<void>
  onFiles: (files: File[]) => Promise<void>
  onRemove: (index: number) => void
  busy: boolean
}) {
  const labels = records.map((record) => record.label)

  const toggleK = (value: number) => {
    const selected = config.k_values.includes(value)
    if (selected && config.k_values.length <= 2) return
    const next = selected
      ? config.k_values.filter((item) => item !== value)
      : [...config.k_values, value].sort((a, b) => a - b)
    const selectedK = next.includes(config.selected_k) ? config.selected_k : next[0]
    setConfig({ ...config, k_values: next, selected_k: selectedK })
  }

  return (
    <AnimatePresence>
      {open ? (
        <>
          <motion.button
            type="button"
            className="fixed inset-0 z-[60] bg-slate-950/55 backdrop-blur-md"
            aria-label="Close analysis setup"
            onClick={close}
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
          />
          <motion.aside
            initial={{ x: '100%' }}
            animate={{ x: 0 }}
            exit={{ x: '100%' }}
            transition={{ type: 'spring', stiffness: 340, damping: 34 }}
            className="glass-panel fixed inset-y-0 right-0 z-[70] w-full max-w-[46rem] overflow-y-auto rounded-none border-y-0 border-r-0"
          >
            <div className="sticky top-0 z-10 flex items-center justify-between border-b border-[var(--border)] bg-[color-mix(in_srgb,var(--bg-solid)_84%,transparent)] px-5 py-4 backdrop-blur-2xl sm:px-7">
              <div className="flex items-center gap-3">
                <span className="grid h-11 w-11 place-items-center rounded-2xl bg-sky-500/10 text-sky-500">
                  <Settings2 size={20} />
                </span>
                <div>
                  <p className="text-base font-semibold">Analysis setup</p>
                  <p className="text-xs text-[var(--text-muted)]">Dataset, scale, and pair configuration</p>
                </div>
              </div>
              <button
                type="button"
                onClick={close}
                className="focus-ring rounded-xl p-2.5 text-[var(--text-muted)] hover:bg-[var(--surface-hover)] hover:text-[var(--text)]"
                aria-label="Close analysis setup"
              >
                <X size={20} />
              </button>
            </div>

            <div className="space-y-6 p-5 sm:p-7">
              <section>
                <div className="mb-4 flex items-center gap-3">
                  <span className="grid h-8 w-8 place-items-center rounded-xl bg-teal-500/10 text-teal-500">
                    <Dna size={17} />
                  </span>
                  <div>
                    <h3 className="text-sm font-semibold">1. Sequence dataset</h3>
                    <p className="text-xs text-[var(--text-muted)]">At least two validated sequences are required.</p>
                  </div>
                </div>
                <FileDropzone records={records} onFiles={onFiles} onRemove={onRemove} disabled={busy} />
              </section>

              <section className="border-t border-[var(--border)] pt-6">
                <div className="mb-4 flex items-center gap-3">
                  <span className="grid h-8 w-8 place-items-center rounded-xl bg-violet-500/10 text-violet-500">
                    <FlaskConical size={17} />
                  </span>
                  <div>
                    <h3 className="text-sm font-semibold">2. Multiscale resolution</h3>
                    <p className="text-xs text-[var(--text-muted)]">Select standard or exploratory k-mer scales.</p>
                  </div>
                </div>

                <div className="grid grid-cols-4 gap-2 sm:grid-cols-8">
                  {allK.map((value) => {
                    const selected = config.k_values.includes(value)
                    return (
                      <button
                        type="button"
                        key={value}
                        onClick={() => toggleK(value)}
                        className={cn(
                          'focus-ring rounded-xl border px-3 py-2.5 text-sm font-semibold transition',
                          selected
                            ? 'border-[var(--primary)] bg-[var(--primary)] text-white shadow-lg shadow-sky-500/15'
                            : 'border-[var(--border)] bg-[var(--bg-solid)]/55 text-[var(--text-soft)] hover:border-[var(--border-strong)] hover:text-[var(--primary)]',
                        )}
                      >
                        k={value}
                      </button>
                    )
                  })}
                </div>
                {config.k_values.some((value) => value > 5) ? (
                  <p className="mt-3 rounded-xl border border-amber-500/20 bg-amber-500/8 px-3 py-2 text-xs leading-5 text-amber-500">
                    k ≥ 6 is an exploratory high-sparsity regime for the current demonstration sequences.
                  </p>
                ) : null}

                <label className="mt-4 block text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">
                  Primary single-scale view
                </label>
                <select
                  value={config.selected_k}
                  onChange={(event) => setConfig({ ...config, selected_k: Number(event.target.value) })}
                  className="focus-ring mt-2 h-12 w-full rounded-xl border border-[var(--border)] bg-[var(--bg-solid)] px-3 text-sm text-[var(--text)]"
                >
                  {config.k_values.map((value) => (
                    <option key={value} value={value}>k={value}</option>
                  ))}
                </select>
              </section>

              <section className="border-t border-[var(--border)] pt-6">
                <div className="mb-4">
                  <h3 className="text-sm font-semibold">3. Selected pair</h3>
                  <p className="mt-1 text-xs text-[var(--text-muted)]">The pair drives trajectories and representation comparisons.</p>
                </div>
                <div className="grid gap-4 sm:grid-cols-[1fr_auto_1fr] sm:items-end">
                  <label className="block">
                    <span className="text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Reference</span>
                    <select
                      value={config.reference_label}
                      onChange={(event) => {
                        const reference = event.target.value
                        const comparison = reference === config.comparison_label
                          ? labels.find((label) => label !== reference) ?? ''
                          : config.comparison_label
                        setConfig({ ...config, reference_label: reference, comparison_label: comparison })
                      }}
                      className="focus-ring mt-2 h-12 w-full rounded-xl border border-[var(--border)] bg-[var(--bg-solid)] px-3 text-sm text-[var(--text)]"
                    >
                      {labels.map((label) => <option key={label}>{label}</option>)}
                    </select>
                  </label>
                  <ChevronRight className="mb-3 hidden text-[var(--text-muted)] sm:block" size={18} />
                  <label className="block">
                    <span className="text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Comparison</span>
                    <select
                      value={config.comparison_label}
                      onChange={(event) => setConfig({ ...config, comparison_label: event.target.value })}
                      className="focus-ring mt-2 h-12 w-full rounded-xl border border-[var(--border)] bg-[var(--bg-solid)] px-3 text-sm text-[var(--text)]"
                    >
                      {labels.filter((label) => label !== config.reference_label).map((label) => <option key={label}>{label}</option>)}
                    </select>
                  </label>
                </div>
              </section>

              <Card className="p-4">
                <div className="flex flex-wrap items-center justify-between gap-3">
                  <div>
                    <p className="text-sm font-semibold">Ready to recompute</p>
                    <p className="mt-1 text-xs text-[var(--text-muted)]">
                      {records.length} sequences · k={config.k_values.join(', ')} · primary k={config.selected_k}
                    </p>
                  </div>
                  <div className="flex gap-2">
                    <Button variant="secondary" icon={RotateCcw} onClick={() => void onReset()} disabled={busy}>
                      Demo
                    </Button>
                    <Button
                      icon={FlaskConical}
                      loading={busy}
                      disabled={records.length < 2 || !config.reference_label || !config.comparison_label}
                      onClick={() => void onAnalyze()}
                    >
                      Run analysis
                    </Button>
                  </div>
                </div>
              </Card>
            </div>
          </motion.aside>
        </>
      ) : null}
    </AnimatePresence>
  )
}
