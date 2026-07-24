import { AnimatePresence, motion } from 'motion/react'
import { Maximize2, X } from 'lucide-react'
import { useState, type ReactNode } from 'react'
import { Card } from './ui'

export function ChartCard({
  title,
  description,
  children,
  toolbar,
  fullscreen = true,
}: {
  title: string
  description?: string
  children: ReactNode
  toolbar?: ReactNode
  fullscreen?: boolean
}) {
  const [expanded, setExpanded] = useState(false)
  const header = (
    <div className="flex flex-col gap-3 border-b border-[var(--border)] px-5 py-4 sm:flex-row sm:items-center sm:justify-between">
      <div>
        <h3 className="font-semibold text-[var(--text)]">{title}</h3>
        {description ? <p className="mt-1 text-xs text-[var(--text-muted)]">{description}</p> : null}
      </div>
      <div className="flex items-center gap-2">
        {toolbar}
        {fullscreen ? (
          <button
            type="button"
            onClick={() => setExpanded(true)}
            className="focus-ring grid h-9 w-9 place-items-center rounded-xl border border-[var(--border)] bg-[var(--bg-solid)]/65 text-[var(--text-muted)] transition hover:border-[var(--border-strong)] hover:text-[var(--primary)]"
            aria-label={`Expand ${title}`}
          >
            <Maximize2 size={16} />
          </button>
        ) : null}
      </div>
    </div>
  )

  return (
    <>
      <Card className="overflow-hidden">
        {header}
        <div className="p-2 sm:p-4">{children}</div>
      </Card>
      <AnimatePresence>
        {expanded ? (
          <motion.div
            className="fixed inset-0 z-[80] flex flex-col bg-[var(--bg)]"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
          >
            <div className="flex items-center justify-between border-b border-[var(--border)] px-5 py-4">
              <div>
                <h2 className="text-lg font-semibold">{title}</h2>
                {description ? <p className="mt-1 text-xs text-[var(--text-muted)]">{description}</p> : null}
              </div>
              <button
                type="button"
                onClick={() => setExpanded(false)}
                className="focus-ring grid h-11 w-11 place-items-center rounded-2xl border border-[var(--border)] bg-[var(--bg-solid)] text-[var(--text-muted)] hover:text-[var(--primary)]"
                aria-label="Close fullscreen chart"
              >
                <X size={20} />
              </button>
            </div>
            <div className="min-h-0 flex-1 overflow-auto p-4 sm:p-7">
              <Card className="min-h-full p-4">{children}</Card>
            </div>
          </motion.div>
        ) : null}
      </AnimatePresence>
    </>
  )
}
