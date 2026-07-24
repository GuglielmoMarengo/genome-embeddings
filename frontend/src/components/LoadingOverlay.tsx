import { AnimatePresence, motion } from 'motion/react'
import { Binary, Braces, ChartSpline, Dna } from 'lucide-react'

const steps = [
  { label: 'Validating sequences', icon: Dna },
  { label: 'Building descriptors', icon: Braces },
  { label: 'Comparing k-mer profiles', icon: Binary },
  { label: 'Rendering geometry', icon: ChartSpline },
]

export function LoadingOverlay({ visible, stage }: { visible: boolean; stage: number }) {
  return (
    <AnimatePresence>
      {visible ? (
        <motion.div
          className="fixed inset-0 z-[90] grid place-items-center bg-[color-mix(in_srgb,var(--bg)_84%,transparent)] p-5 backdrop-blur-xl"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
        >
          <motion.div
            initial={{ y: 18, scale: 0.96, opacity: 0 }}
            animate={{ y: 0, scale: 1, opacity: 1 }}
            exit={{ y: 14, scale: 0.97, opacity: 0 }}
            className="glass-panel w-full max-w-xl rounded-[2rem] p-7 sm:p-9"
          >
            <div className="flex items-center gap-4">
              <div className="relative grid h-14 w-14 place-items-center rounded-2xl bg-[linear-gradient(135deg,var(--primary),var(--secondary))] text-white shadow-lg shadow-sky-500/20">
                <Dna size={25} />
                <span className="absolute -inset-2 animate-ping rounded-[1.35rem] border border-sky-400/20" />
              </div>
              <div>
                <p className="text-xs font-semibold uppercase tracking-[0.18em] text-[var(--primary)]">
                  Scientific workflow
                </p>
                <h2 className="mt-1 text-xl font-semibold">Computing genome embeddings</h2>
              </div>
            </div>

            <div className="mt-8 space-y-3">
              {steps.map(({ label, icon: Icon }, index) => {
                const active = index === stage
                const complete = index < stage
                return (
                  <motion.div
                    key={label}
                    animate={{ opacity: index <= stage ? 1 : 0.42, x: active ? 5 : 0 }}
                    className="flex items-center gap-3 rounded-xl border border-[var(--border)] bg-[var(--bg-solid)]/55 px-4 py-3"
                  >
                    <span
                      className={`grid h-9 w-9 place-items-center rounded-xl ${
                        complete
                          ? 'bg-emerald-500/12 text-emerald-500'
                          : active
                            ? 'bg-sky-500/12 text-sky-500'
                            : 'bg-[var(--bg-subtle)] text-[var(--text-muted)]'
                      }`}
                    >
                      {complete ? '✓' : <Icon size={17} />}
                    </span>
                    <span className={active ? 'font-semibold text-[var(--text)]' : 'text-sm text-[var(--text-soft)]'}>
                      {label}
                    </span>
                    {active ? (
                      <span className="ml-auto flex gap-1">
                        {[0, 1, 2].map((dot) => (
                          <motion.span
                            key={dot}
                            className="h-1.5 w-1.5 rounded-full bg-sky-500"
                            animate={{ y: [0, -5, 0] }}
                            transition={{ delay: dot * 0.13, duration: 0.7, repeat: Infinity }}
                          />
                        ))}
                      </span>
                    ) : null}
                  </motion.div>
                )
              })}
            </div>
          </motion.div>
        </motion.div>
      ) : null}
    </AnimatePresence>
  )
}
