import { motion } from 'motion/react'
import { Dna, FileText, FlaskConical, Search, Settings2 } from 'lucide-react'
import { useMemo, useState } from 'react'
import type { AnalysisData, SequenceRecord } from '../types'
import { Badge, Button, Card, SectionHeading } from '../components/ui'

export function DatasetPage({
  analysis,
  records,
  openSetup,
}: {
  analysis: AnalysisData
  records: SequenceRecord[]
  openSetup: () => void
}) {
  const [query, setQuery] = useState('')
  const rows = useMemo(() => {
    const needle = query.trim().toLowerCase()
    if (!needle) return analysis.dataset
    return analysis.dataset.filter((row) =>
      [row.label, row.source, row.header].some((value) => value.toLowerCase().includes(needle)),
    )
  }, [analysis.dataset, query])

  const totalBases = analysis.dataset.reduce((total, row) => total + row.length, 0)
  const gcValues = analysis.dataset.map((row) => row.gc_content)
  const minGc = Math.min(...gcValues)
  const maxGc = Math.max(...gcValues)

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Sequence inventory"
        title="Dataset"
        description="Inspect provenance, sequence length, composition, and FASTA metadata before interpreting any downstream geometry."
        action={<Button icon={Settings2} onClick={openSetup}>Manage dataset</Button>}
      />

      <div className="grid gap-4 sm:grid-cols-3">
        <Card className="p-5">
          <p className="text-xs font-semibold uppercase tracking-[0.16em] text-[var(--text-muted)]">Sequence records</p>
          <p className="mt-3 text-3xl font-semibold">{analysis.dataset.length}</p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">Single-record FASTA inputs</p>
        </Card>
        <Card className="p-5">
          <p className="text-xs font-semibold uppercase tracking-[0.16em] text-[var(--text-muted)]">Total nucleotides</p>
          <p className="mt-3 text-3xl font-semibold">{totalBases.toLocaleString()}</p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">Across the active dataset</p>
        </Card>
        <Card className="p-5">
          <p className="text-xs font-semibold uppercase tracking-[0.16em] text-[var(--text-muted)]">GC range</p>
          <p className="mt-3 text-3xl font-semibold">{(minGc * 100).toFixed(1)}–{(maxGc * 100).toFixed(1)}%</p>
          <p className="mt-2 text-sm text-[var(--text-muted)]">Minimum to maximum composition</p>
        </Card>
      </div>

      <Card className="overflow-hidden">
        <div className="flex flex-col gap-4 border-b border-[var(--border)] p-5 sm:flex-row sm:items-center sm:justify-between sm:px-6">
          <div>
            <h3 className="font-semibold">Active sequence records</h3>
            <p className="mt-1 text-xs text-[var(--text-muted)]">Search labels, sources, or FASTA headers.</p>
          </div>
          <label className="relative block w-full sm:max-w-xs">
            <Search className="pointer-events-none absolute left-3 top-1/2 -translate-y-1/2 text-[var(--text-muted)]" size={17} />
            <input
              value={query}
              onChange={(event) => setQuery(event.target.value)}
              placeholder="Search dataset"
              className="focus-ring h-11 w-full rounded-xl border border-[var(--border)] bg-[var(--bg-solid)] pl-10 pr-3 text-sm text-[var(--text)] placeholder:text-[var(--text-muted)]"
            />
          </label>
        </div>

        <div className="overflow-x-auto table-scroll">
          <table className="w-full min-w-[72rem] border-collapse text-left">
            <thead>
              <tr className="border-b border-[var(--border)] text-xs uppercase tracking-[0.14em] text-[var(--text-muted)]">
                <th className="px-5 py-4 font-semibold sm:px-6">Sequence</th>
                <th className="px-5 py-4 font-semibold">Source</th>
                <th className="px-5 py-4 text-right font-semibold">Length</th>
                <th className="px-5 py-4 text-right font-semibold">GC</th>
                <th className="px-5 py-4 font-semibold sm:px-6">FASTA header</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((row, index) => (
                <motion.tr
                  key={row.label}
                  initial={{ opacity: 0, y: 5 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: Math.min(index * 0.035, 0.25) }}
                  className="border-b border-[var(--border)] last:border-0 hover:bg-[var(--surface-hover)]"
                >
                  <td className="px-5 py-4 sm:px-6">
                    <div className="flex items-center gap-3">
                      <span className="grid h-10 w-10 place-items-center rounded-xl bg-sky-500/10 text-sky-500"><Dna size={18} /></span>
                      <div>
                        <p className="font-semibold text-[var(--text)]">{row.label}</p>
                        <Badge tone={row.label.toLowerCase().includes('periodic') ? 'warning' : 'info'}>
                          {row.label.toLowerCase().includes('periodic') ? 'Synthetic control' : 'Biological sequence'}
                        </Badge>
                      </div>
                    </div>
                  </td>
                  <td className="max-w-[18rem] truncate px-5 py-4 text-sm text-[var(--text-soft)]" title={row.source}>{row.source}</td>
                  <td className="px-5 py-4 text-right font-mono text-sm">{row.length.toLocaleString()}</td>
                  <td className="px-5 py-4 text-right font-mono text-sm">{(row.gc_content * 100).toFixed(2)}%</td>
                  <td className="max-w-[32rem] truncate px-5 py-4 text-xs text-[var(--text-muted)] sm:px-6" title={row.header}>{row.header}</td>
                </motion.tr>
              ))}
            </tbody>
          </table>
        </div>
      </Card>

      <div className="grid gap-4 lg:grid-cols-2">
        {records.map((record, index) => (
          <Card key={`${record.label}-${index}`} className="p-5">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-teal-500/10 text-teal-500"><FileText size={20} /></span>
              <div className="min-w-0 flex-1">
                <div className="flex flex-wrap items-center justify-between gap-2">
                  <h3 className="font-semibold">{record.label}</h3>
                  <span className="font-mono text-xs text-[var(--text-muted)]">{record.sequence.length.toLocaleString()} nt</span>
                </div>
                <p className="mt-2 truncate text-xs text-[var(--text-muted)]" title={record.header}>{record.header || 'No FASTA header'}</p>
                <div className="mt-4 overflow-hidden rounded-xl border border-[var(--border)] bg-[var(--bg-subtle)] p-3 font-mono text-xs leading-5 text-[var(--text-soft)]">
                  {record.sequence.slice(0, 120)}{record.sequence.length > 120 ? '…' : ''}
                </div>
              </div>
            </div>
          </Card>
        ))}
      </div>

      <Card className="flex flex-col items-start justify-between gap-4 p-5 sm:flex-row sm:items-center">
        <div className="flex items-center gap-3">
          <span className="grid h-10 w-10 place-items-center rounded-xl bg-violet-500/10 text-violet-500"><FlaskConical size={18} /></span>
          <div>
            <p className="font-semibold">Dataset changes require recomputation</p>
            <p className="mt-1 text-sm text-[var(--text-muted)]">Open the analysis setup to upload, remove, or reset sequence records.</p>
          </div>
        </div>
        <Button variant="secondary" onClick={openSetup}>Open setup</Button>
      </Card>
    </div>
  )
}
