import { AnimatePresence, motion } from 'motion/react'
import { FilePlus2, FileText, UploadCloud, X } from 'lucide-react'
import { useRef, useState } from 'react'
import type { SequenceRecord } from '../types'
import { Button, cn } from './ui'

const accepted = '.fasta,.fa,.fna,.txt'

export function FileDropzone({
  records,
  onFiles,
  onRemove,
  disabled = false,
}: {
  records: SequenceRecord[]
  onFiles: (files: File[]) => Promise<void>
  onRemove: (index: number) => void
  disabled?: boolean
}) {
  const input = useRef<HTMLInputElement>(null)
  const [dragging, setDragging] = useState(false)

  const receive = async (files: FileList | null) => {
    if (!files?.length || disabled) return
    await onFiles(Array.from(files))
    if (input.current) input.current.value = ''
  }

  return (
    <div className="space-y-4">
      <button
        type="button"
        disabled={disabled}
        onClick={() => input.current?.click()}
        onDragEnter={(event) => {
          event.preventDefault()
          setDragging(true)
        }}
        onDragOver={(event) => event.preventDefault()}
        onDragLeave={(event) => {
          event.preventDefault()
          if (event.currentTarget === event.target) setDragging(false)
        }}
        onDrop={(event) => {
          event.preventDefault()
          setDragging(false)
          void receive(event.dataTransfer.files)
        }}
        className={cn(
          'focus-ring relative grid min-h-52 w-full place-items-center overflow-hidden rounded-[1.4rem] border border-dashed p-6 text-center transition-all duration-200',
          dragging
            ? 'scale-[1.01] border-[var(--primary)] bg-[color-mix(in_srgb,var(--primary)_10%,transparent)]'
            : 'border-[var(--border-strong)] bg-[var(--bg-subtle)]/55 hover:border-[var(--primary)] hover:bg-[var(--surface-hover)]',
          disabled && 'cursor-not-allowed opacity-55',
        )}
      >
        <input
          ref={input}
          type="file"
          accept={accepted}
          multiple
          className="hidden"
          onChange={(event) => void receive(event.target.files)}
        />
        <div className="pointer-events-none relative z-10 max-w-md">
          <motion.span
            animate={dragging ? { scale: 1.1, y: -4 } : { scale: 1, y: 0 }}
            className="mx-auto grid h-16 w-16 place-items-center rounded-[1.3rem] bg-[linear-gradient(135deg,var(--primary),var(--secondary))] text-white shadow-xl shadow-sky-500/20"
          >
            {dragging ? <FilePlus2 size={27} /> : <UploadCloud size={27} />}
          </motion.span>
          <h3 className="mt-5 text-base font-semibold text-[var(--text)]">
            {dragging ? 'Release to add sequences' : 'Drop FASTA files here'}
          </h3>
          <p className="mt-2 text-sm leading-6 text-[var(--text-muted)]">
            Single-record FASTA files are validated by the Python backend before they enter the workspace.
          </p>
          <span className="mt-4 inline-flex rounded-full border border-[var(--border)] bg-[var(--bg-solid)]/65 px-3 py-1 text-xs font-medium text-[var(--text-soft)]">
            FASTA · FA · FNA · TXT
          </span>
        </div>
        <div className="dna-grid absolute inset-0 opacity-30" />
      </button>

      <AnimatePresence initial={false}>
        {records.length ? (
          <motion.div layout className="grid gap-2 sm:grid-cols-2">
            {records.map((record, index) => (
              <motion.div
                layout
                key={`${record.label}-${record.source}-${index}`}
                initial={{ opacity: 0, y: 8 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, scale: 0.96 }}
                className="flex items-center gap-3 rounded-2xl border border-[var(--border)] bg-[var(--bg-solid)]/65 p-3"
              >
                <span className="grid h-10 w-10 shrink-0 place-items-center rounded-xl bg-sky-500/10 text-sky-500">
                  <FileText size={18} />
                </span>
                <div className="min-w-0 flex-1 text-left">
                  <p className="truncate text-sm font-semibold text-[var(--text)]">{record.label}</p>
                  <p className="truncate text-xs text-[var(--text-muted)]">
                    {record.sequence.length.toLocaleString()} nt · {record.source}
                  </p>
                </div>
                <Button
                  variant="ghost"
                  className="h-9 min-h-9 w-9 px-0"
                  onClick={(event) => {
                    event.stopPropagation()
                    onRemove(index)
                  }}
                  aria-label={`Remove ${record.label}`}
                >
                  <X size={16} />
                </Button>
              </motion.div>
            ))}
          </motion.div>
        ) : null}
      </AnimatePresence>
    </div>
  )
}
