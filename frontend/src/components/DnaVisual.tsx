import { motion } from 'motion/react'

const rungs = Array.from({ length: 13 }, (_, index) => index)

export function DnaVisual({ compact = false }: { compact?: boolean }) {
  return (
    <div
      className={`relative isolate overflow-hidden rounded-[2rem] border border-[var(--border)] bg-[linear-gradient(145deg,color-mix(in_srgb,var(--bg-solid)_92%,transparent),color-mix(in_srgb,var(--primary)_7%,var(--bg-solid)))] ${
        compact ? 'h-52' : 'h-[25rem]'
      }`}
      aria-hidden="true"
    >
      <div className="dna-grid absolute inset-0 opacity-55" />
      <motion.div
        className="absolute left-[13%] top-[8%] h-40 w-40 rounded-full bg-sky-400/20 blur-3xl"
        animate={{ x: [0, 24, -10, 0], y: [0, 12, 28, 0], scale: [1, 1.18, 0.95, 1] }}
        transition={{ duration: 10, repeat: Infinity, ease: 'easeInOut' }}
      />
      <motion.div
        className="absolute bottom-[7%] right-[10%] h-44 w-44 rounded-full bg-violet-500/15 blur-3xl"
        animate={{ x: [0, -18, 12, 0], y: [0, -25, 8, 0], scale: [1, 0.9, 1.2, 1] }}
        transition={{ duration: 12, repeat: Infinity, ease: 'easeInOut' }}
      />

      <div className="absolute inset-0 grid place-items-center">
        <motion.svg
          viewBox="0 0 420 420"
          className={compact ? 'h-48 w-48' : 'h-[22rem] w-[22rem]'}
          initial={{ opacity: 0, scale: 0.9, rotate: -7 }}
          animate={{ opacity: 1, scale: 1, rotate: 0 }}
          transition={{ duration: 0.9, ease: [0.22, 1, 0.36, 1] }}
        >
          <defs>
            <linearGradient id="helix-a" x1="0" y1="0" x2="1" y2="1">
              <stop offset="0%" stopColor="#38bdf8" />
              <stop offset="55%" stopColor="#2dd4bf" />
              <stop offset="100%" stopColor="#a78bfa" />
            </linearGradient>
            <linearGradient id="helix-b" x1="1" y1="0" x2="0" y2="1">
              <stop offset="0%" stopColor="#a78bfa" />
              <stop offset="52%" stopColor="#38bdf8" />
              <stop offset="100%" stopColor="#2dd4bf" />
            </linearGradient>
            <filter id="glow">
              <feGaussianBlur stdDeviation="4" result="coloredBlur" />
              <feMerge>
                <feMergeNode in="coloredBlur" />
                <feMergeNode in="SourceGraphic" />
              </feMerge>
            </filter>
          </defs>

          <motion.path
            d="M112 38 C330 115 88 300 308 382"
            fill="none"
            stroke="url(#helix-a)"
            strokeWidth="12"
            strokeLinecap="round"
            filter="url(#glow)"
            pathLength={1}
            initial={{ pathLength: 0 }}
            animate={{ pathLength: 1 }}
            transition={{ duration: 1.4, ease: 'easeInOut' }}
          />
          <motion.path
            d="M308 38 C90 115 332 300 112 382"
            fill="none"
            stroke="url(#helix-b)"
            strokeWidth="12"
            strokeLinecap="round"
            filter="url(#glow)"
            pathLength={1}
            initial={{ pathLength: 0 }}
            animate={{ pathLength: 1 }}
            transition={{ duration: 1.4, delay: 0.12, ease: 'easeInOut' }}
          />

          {rungs.map((rung) => {
            const y = 62 + rung * 25
            const phase = rung / (rungs.length - 1)
            const wave = Math.sin(phase * Math.PI * 4.1)
            const x1 = 210 - wave * 92
            const x2 = 210 + wave * 92
            return (
              <motion.line
                key={rung}
                x1={x1}
                y1={y}
                x2={x2}
                y2={y}
                stroke="url(#helix-a)"
                strokeWidth="5"
                strokeLinecap="round"
                opacity="0.72"
                initial={{ opacity: 0, pathLength: 0 }}
                animate={{ opacity: 0.72, pathLength: 1 }}
                transition={{ delay: 0.3 + rung * 0.045, duration: 0.5 }}
              />
            )
          })}
        </motion.svg>
      </div>

      <motion.div
        className="absolute right-6 top-6 rounded-full border border-sky-400/20 bg-sky-400/10 px-3 py-1.5 text-[11px] font-semibold uppercase tracking-[0.18em] text-sky-500"
        animate={{ y: [0, -4, 0] }}
        transition={{ duration: 3, repeat: Infinity, ease: 'easeInOut' }}
      >
        interpretable vectors
      </motion.div>
      <motion.div
        className="absolute bottom-6 left-6 rounded-full border border-teal-400/20 bg-teal-400/10 px-3 py-1.5 text-[11px] font-semibold uppercase tracking-[0.18em] text-teal-500"
        animate={{ y: [0, 4, 0] }}
        transition={{ duration: 3.6, repeat: Infinity, ease: 'easeInOut' }}
      >
        multiscale geometry
      </motion.div>
    </div>
  )
}
