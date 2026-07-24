import { AnimatePresence, motion } from 'motion/react'
import { CheckCircle2, CircleAlert, Info, X } from 'lucide-react'
import type { ToastMessage } from '../types'

export function ToastStack({
  messages,
  dismiss,
}: {
  messages: ToastMessage[]
  dismiss: (id: number) => void
}) {
  const icons = {
    success: CheckCircle2,
    error: CircleAlert,
    info: Info,
  }
  const styles = {
    success: 'text-emerald-500 bg-emerald-500/10',
    error: 'text-rose-500 bg-rose-500/10',
    info: 'text-sky-500 bg-sky-500/10',
  }

  return (
    <div className="pointer-events-none fixed right-4 top-4 z-[100] flex w-[min(25rem,calc(100vw-2rem))] flex-col gap-3">
      <AnimatePresence initial={false}>
        {messages.map((message) => {
          const Icon = icons[message.tone]
          return (
            <motion.div
              key={message.id}
              initial={{ opacity: 0, x: 35, scale: 0.96 }}
              animate={{ opacity: 1, x: 0, scale: 1 }}
              exit={{ opacity: 0, x: 35, scale: 0.96 }}
              transition={{ type: 'spring', stiffness: 420, damping: 34 }}
              className="glass-panel pointer-events-auto flex items-start gap-3 rounded-2xl p-4"
            >
              <span className={`mt-0.5 grid h-9 w-9 shrink-0 place-items-center rounded-xl ${styles[message.tone]}`}>
                <Icon size={18} />
              </span>
              <div className="min-w-0 flex-1">
                <p className="font-semibold text-[var(--text)]">{message.title}</p>
                {message.description ? (
                  <p className="mt-1 text-sm leading-5 text-[var(--text-soft)]">
                    {message.description}
                  </p>
                ) : null}
              </div>
              <button
                type="button"
                className="focus-ring rounded-lg p-1.5 text-[var(--text-muted)] hover:bg-[var(--surface-hover)] hover:text-[var(--text)]"
                onClick={() => dismiss(message.id)}
                aria-label="Dismiss notification"
              >
                <X size={16} />
              </button>
            </motion.div>
          )
        })}
      </AnimatePresence>
    </div>
  )
}
