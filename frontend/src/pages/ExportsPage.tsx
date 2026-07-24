import { Braces, CheckCircle2, Download, FileJson2, FileSpreadsheet, ShieldCheck } from 'lucide-react'
import type { AnalysisConfig, AnalysisData, SequenceRecord } from '../types'
import { Badge, Button, Card, SectionHeading } from '../components/ui'

export function ExportsPage({
  analysis,
  records,
  config,
  onDownload,
  busy,
}: {
  analysis: AnalysisData
  records: SequenceRecord[]
  config: AnalysisConfig
  onDownload: (kind: 'json' | 'csv') => Promise<void>
  busy: boolean
}) {
  const jsonPreview = JSON.stringify({
    summary: analysis.summary,
    matrices: Object.keys(analysis.matrices),
    trajectories: Object.keys(analysis.pair_trajectories),
  }, null, 2)

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Reusable outputs"
        title="Exports"
        description="Download complete structured analysis data or a compact summary suitable for spreadsheets and reporting workflows."
      />

      <div className="grid gap-6 xl:grid-cols-2">
        <Card className="relative overflow-hidden p-6 sm:p-7">
          <div className="absolute right-0 top-0 h-48 w-48 rounded-full bg-sky-500/12 blur-3xl" />
          <div className="relative z-10">
            <span className="grid h-13 w-13 place-items-center rounded-2xl bg-sky-500/10 text-sky-500"><FileJson2 size={23} /></span>
            <div className="mt-5 flex flex-wrap items-center gap-2">
              <h3 className="text-xl font-semibold">Complete JSON analysis</h3>
              <Badge tone="info">Structured</Badge>
            </div>
            <p className="mt-3 text-sm leading-6 text-[var(--text-soft)]">
              Includes summary metrics, sequence inventory, Descriptor V2 diagnostics, all matrices, pair trajectories, step distances, rankings, and stability statistics.
            </p>
            <ul className="mt-5 space-y-2 text-sm text-[var(--text-soft)]">
              {['Machine-readable scientific output', 'Exact floating-point values', 'Matrix labels and analytical context'].map((item) => (
                <li key={item} className="flex items-center gap-2"><CheckCircle2 size={16} className="text-emerald-500" />{item}</li>
              ))}
            </ul>
            <Button className="mt-6" icon={Download} loading={busy} onClick={() => void onDownload('json')}>
              Download JSON
            </Button>
          </div>
        </Card>

        <Card className="relative overflow-hidden p-6 sm:p-7">
          <div className="absolute right-0 top-0 h-48 w-48 rounded-full bg-teal-500/12 blur-3xl" />
          <div className="relative z-10">
            <span className="grid h-13 w-13 place-items-center rounded-2xl bg-teal-500/10 text-teal-500"><FileSpreadsheet size={23} /></span>
            <div className="mt-5 flex flex-wrap items-center gap-2">
              <h3 className="text-xl font-semibold">Summary CSV</h3>
              <Badge tone="success">Spreadsheet-ready</Badge>
            </div>
            <p className="mt-3 text-sm leading-6 text-[var(--text-soft)]">
              A compact two-column summary of the active configuration and selected pair values, encoded with UTF-8 BOM for reliable spreadsheet opening.
            </p>
            <ul className="mt-5 space-y-2 text-sm text-[var(--text-soft)]">
              {['Reference and comparison labels', 'Selected k and multiscale range', 'Legacy, V2, embedding, and JS values'].map((item) => (
                <li key={item} className="flex items-center gap-2"><CheckCircle2 size={16} className="text-emerald-500" />{item}</li>
              ))}
            </ul>
            <Button variant="secondary" className="mt-6" icon={Download} loading={busy} onClick={() => void onDownload('csv')}>
              Download summary CSV
            </Button>
          </div>
        </Card>
      </div>

      <div className="grid gap-6 xl:grid-cols-[1.1fr_.9fr]">
        <Card className="overflow-hidden">
          <div className="flex items-center gap-3 border-b border-[var(--border)] px-5 py-4 sm:px-6">
            <span className="grid h-9 w-9 place-items-center rounded-xl bg-violet-500/10 text-violet-500"><Braces size={18} /></span>
            <div>
              <h3 className="font-semibold">JSON structure preview</h3>
              <p className="mt-0.5 text-xs text-[var(--text-muted)]">Abbreviated client-side preview</p>
            </div>
          </div>
          <pre className="max-h-[30rem] overflow-auto p-5 font-mono text-xs leading-6 text-[var(--text-soft)] sm:p-6">{jsonPreview}</pre>
        </Card>

        <div className="space-y-4">
          <Card className="p-5">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-emerald-500/10 text-emerald-500"><ShieldCheck size={20} /></span>
              <div>
                <h3 className="font-semibold">Stateless export workflow</h3>
                <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                  The browser sends the active dataset and configuration to the API. Uploaded filenames are sanitized and never interpreted as server filesystem paths.
                </p>
              </div>
            </div>
          </Card>
          <Card className="p-5">
            <p className="text-xs font-semibold uppercase tracking-[0.15em] text-[var(--text-muted)]">Active export context</p>
            <dl className="mt-4 space-y-3 text-sm">
              <div className="flex justify-between gap-4"><dt className="text-[var(--text-muted)]">Sequences</dt><dd className="font-semibold">{records.length}</dd></div>
              <div className="flex justify-between gap-4"><dt className="text-[var(--text-muted)]">k values</dt><dd className="font-mono">{config.k_values.join(', ')}</dd></div>
              <div className="flex justify-between gap-4"><dt className="text-[var(--text-muted)]">Primary k</dt><dd className="font-mono">{config.selected_k}</dd></div>
              <div className="flex justify-between gap-4"><dt className="text-[var(--text-muted)]">Reference</dt><dd className="max-w-[12rem] truncate font-semibold" title={config.reference_label}>{config.reference_label}</dd></div>
              <div className="flex justify-between gap-4"><dt className="text-[var(--text-muted)]">Comparison</dt><dd className="max-w-[12rem] truncate font-semibold" title={config.comparison_label}>{config.comparison_label}</dd></div>
            </dl>
          </Card>
        </div>
      </div>
    </div>
  )
}
