import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { AnimatePresence, motion } from 'motion/react'
import { AlertTriangle, RefreshCcw } from 'lucide-react'
import { analyze, downloadExport, loadDemo, parseFasta } from './lib/api'
import type {
  AnalysisConfig,
  AnalysisData,
  AnalysisRequest,
  PageId,
  SequenceRecord,
  ToastMessage,
} from './types'
import { AppShell } from './components/AppShell'
import { AnalysisControls } from './components/AnalysisControls'
import { LoadingOverlay } from './components/LoadingOverlay'
import { ThemeProvider } from './components/ThemeProvider'
import { ToastStack } from './components/ToastStack'
import { Button, Card } from './components/ui'
import { OverviewPage } from './pages/OverviewPage'
import { DatasetPage } from './pages/DatasetPage'
import { DescriptorsPage } from './pages/DescriptorsPage'
import { MatricesPage } from './pages/MatricesPage'
import { MultiscalePage } from './pages/MultiscalePage'
import { ExportsPage } from './pages/ExportsPage'
import { MethodologyPage } from './pages/MethodologyPage'

function uniqueLabel(records: SequenceRecord[], requested: string): string {
  const existing = new Set(records.map((record) => record.label))
  if (!existing.has(requested)) return requested
  let index = 2
  while (existing.has(`${requested} (${index})`)) index += 1
  return `${requested} (${index})`
}

function AppContent() {
  const [page, setPage] = useState<PageId>('overview')
  const [records, setRecords] = useState<SequenceRecord[]>([])
  const [analysisData, setAnalysisData] = useState<AnalysisData | null>(null)
  const [config, setConfig] = useState<AnalysisConfig>({
    k_values: [1, 2, 3, 4, 5],
    selected_k: 3,
    reference_label: '',
    comparison_label: '',
  })
  const [setupOpen, setSetupOpen] = useState(false)
  const [busy, setBusy] = useState(true)
  const [exporting, setExporting] = useState(false)
  const [loadingStage, setLoadingStage] = useState(0)
  const [fatalError, setFatalError] = useState<string | null>(null)
  const [toasts, setToasts] = useState<ToastMessage[]>([])
  const toastId = useRef(0)

  const notify = useCallback((tone: ToastMessage['tone'], title: string, description?: string) => {
    const id = ++toastId.current
    setToasts((current) => [...current, { id, tone, title, description }])
    window.setTimeout(() => {
      setToasts((current) => current.filter((message) => message.id !== id))
    }, tone === 'error' ? 6500 : 4200)
  }, [])

  const applyDemo = useCallback(async () => {
    setBusy(true)
    setLoadingStage(0)
    setFatalError(null)
    try {
      const response = await loadDemo()
      setRecords(response.records)
      setAnalysisData(response.analysis)
      setConfig({
        k_values: response.analysis.summary.k_values,
        selected_k: response.analysis.summary.selected_k,
        reference_label: response.analysis.summary.reference_label,
        comparison_label: response.analysis.summary.comparison_label,
      })
      notify('success', 'Demonstration dataset ready', 'Six validated sequences and all V2 analyses are available.')
    } catch (error) {
      const message = error instanceof Error ? error.message : 'Unable to load the scientific API.'
      setFatalError(message)
    } finally {
      setBusy(false)
    }
  }, [notify])

  useEffect(() => {
    void applyDemo()
  }, [applyDemo])

  useEffect(() => {
    if (!busy) return
    const interval = window.setInterval(() => {
      setLoadingStage((current) => Math.min(3, current + 1))
    }, 520)
    return () => window.clearInterval(interval)
  }, [busy])

  const request = useMemo<AnalysisRequest>(() => ({ records, config }), [records, config])

  const runAnalysis = async () => {
    setBusy(true)
    setLoadingStage(0)
    try {
      const result = await analyze(request)
      setAnalysisData(result)
      setSetupOpen(false)
      notify(
        'success',
        'Analysis complete',
        `${records.length} sequences · k=${config.k_values.join(', ')} · primary k=${config.selected_k}`,
      )
    } catch (error) {
      notify('error', 'Analysis failed', error instanceof Error ? error.message : 'Unexpected API error.')
    } finally {
      setBusy(false)
    }
  }

  const addFiles = async (files: File[]) => {
    setBusy(true)
    setLoadingStage(0)
    const added: SequenceRecord[] = []
    const failures: string[] = []
    try {
      for (const file of files) {
        try {
          const parsed = await parseFasta(file)
          const current = [...records, ...added]
          added.push({ ...parsed, label: uniqueLabel(current, parsed.label) })
        } catch (error) {
          failures.push(`${file.name}: ${error instanceof Error ? error.message : 'Invalid FASTA file.'}`)
        }
      }

      if (added.length) {
        const nextRecords = [...records, ...added]
        setRecords(nextRecords)
        if (!config.reference_label && nextRecords[0]) {
          setConfig((current) => ({
            ...current,
            reference_label: nextRecords[0].label,
            comparison_label: nextRecords[1]?.label ?? '',
          }))
        }
        notify('success', `${added.length} sequence${added.length === 1 ? '' : 's'} added`, 'Run the analysis to include the new records.')
      }
      if (failures.length) notify('error', 'Some files were rejected', failures.join(' '))
    } finally {
      setBusy(false)
    }
  }

  const removeRecord = (index: number) => {
    const next = records.filter((_, candidate) => candidate !== index)
    setRecords(next)
    setConfig((current) => {
      let reference = current.reference_label
      let comparison = current.comparison_label
      if (!next.some((record) => record.label === reference)) reference = next[0]?.label ?? ''
      if (!next.some((record) => record.label === comparison) || comparison === reference) {
        comparison = next.find((record) => record.label !== reference)?.label ?? ''
      }
      return { ...current, reference_label: reference, comparison_label: comparison }
    })
  }

  const download = async (kind: 'json' | 'csv') => {
    setExporting(true)
    try {
      await downloadExport(kind, request)
      notify('success', `${kind.toUpperCase()} export generated`, 'The download contains the current dataset and analysis configuration.')
    } catch (error) {
      notify('error', 'Export failed', error instanceof Error ? error.message : 'Unexpected export error.')
    } finally {
      setExporting(false)
    }
  }

  if (fatalError) {
    return (
      <div className="grid min-h-screen place-items-center p-5">
        <Card className="max-w-xl p-8 text-center">
          <span className="mx-auto grid h-16 w-16 place-items-center rounded-2xl bg-rose-500/10 text-rose-500"><AlertTriangle size={28} /></span>
          <h1 className="mt-6 text-2xl font-semibold">Unable to connect to the scientific backend</h1>
          <p className="mt-3 text-sm leading-6 text-[var(--text-soft)]">{fatalError}</p>
          <p className="mt-3 rounded-xl bg-[var(--bg-subtle)] p-3 font-mono text-xs text-[var(--text-muted)]">
            Start the API with: python app.py
          </p>
          <Button className="mt-6" icon={RefreshCcw} onClick={() => void applyDemo()}>Retry connection</Button>
        </Card>
      </div>
    )
  }

  const pageContent = analysisData ? {
    overview: <OverviewPage analysis={analysisData} openSetup={() => setSetupOpen(true)} goToMatrices={() => setPage('matrices')} />,
    dataset: <DatasetPage analysis={analysisData} records={records} openSetup={() => setSetupOpen(true)} />,
    descriptors: <DescriptorsPage analysis={analysisData} />,
    matrices: <MatricesPage analysis={analysisData} />,
    multiscale: <MultiscalePage analysis={analysisData} />,
    exports: <ExportsPage analysis={analysisData} records={records} config={config} onDownload={download} busy={exporting} />,
    methodology: <MethodologyPage />,
  }[page] : (
    <Card className="shimmer min-h-[36rem] p-8" />
  )

  return (
    <>
      <AppShell
        page={page}
        setPage={setPage}
        analysisReady={Boolean(analysisData) && !busy}
        openSetup={() => setSetupOpen(true)}
      >
        {pageContent}
      </AppShell>
      <AnalysisControls
        open={setupOpen}
        close={() => setSetupOpen(false)}
        records={records}
        config={config}
        setConfig={setConfig}
        onAnalyze={runAnalysis}
        onReset={applyDemo}
        onFiles={addFiles}
        onRemove={removeRecord}
        busy={busy}
      />
      <LoadingOverlay visible={busy} stage={loadingStage} />
      <ToastStack
        messages={toasts}
        dismiss={(id) => setToasts((current) => current.filter((message) => message.id !== id))}
      />
      <AnimatePresence>
        {analysisData && !setupOpen && page === 'overview' ? (
          <motion.button
            type="button"
            onClick={() => setSetupOpen(true)}
            className="focus-ring fixed bottom-5 right-5 z-30 flex h-13 items-center gap-2 rounded-2xl bg-[var(--primary)] px-4 text-sm font-semibold text-white shadow-2xl shadow-sky-500/30 sm:hidden"
            initial={{ y: 30, opacity: 0 }}
            animate={{ y: 0, opacity: 1 }}
            exit={{ y: 30, opacity: 0 }}
          >
            Configure
          </motion.button>
        ) : null}
      </AnimatePresence>
    </>
  )
}

export default function App() {
  return (
    <ThemeProvider>
      <AppContent />
    </ThemeProvider>
  )
}
