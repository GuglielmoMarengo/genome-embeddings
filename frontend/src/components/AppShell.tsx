import { AnimatePresence, motion } from 'motion/react'
import {
  Activity,
  Atom,
  BarChart3,
  ChevronLeft,
  ChevronRight,
  Database,
  Download,
  Dna,
  FileBarChart,
  Menu,
  Moon,
  Network,
  PanelLeftClose,
  PanelLeftOpen,
  Settings2,
  Sun,
  X,
} from 'lucide-react'
import { useState, type ReactNode } from 'react'
import type { PageId } from '../types'
import { useThemeContext } from './ThemeProvider'
import { cn } from './ui'

const navigation: Array<{
  id: PageId
  label: string
  description: string
  icon: typeof Activity
}> = [
  { id: 'overview', label: 'Overview', description: 'Analysis summary', icon: Activity },
  { id: 'dataset', label: 'Dataset', description: 'Sequence inventory', icon: Database },
  { id: 'descriptors', label: 'Descriptors V2', description: 'Finite-sample diagnostics', icon: Atom },
  { id: 'matrices', label: 'Matrices', description: 'Pairwise geometry', icon: Network },
  { id: 'multiscale', label: 'Multiscale', description: 'Trajectories and stability', icon: BarChart3 },
  { id: 'exports', label: 'Exports', description: 'Reusable analysis data', icon: Download },
  { id: 'methodology', label: 'Methodology', description: 'Definitions and limits', icon: FileBarChart },
]

function Brand({ compact = false }: { compact?: boolean }) {
  return (
    <div className="flex min-w-0 items-center gap-3">
      <div className="relative grid h-11 w-11 shrink-0 place-items-center overflow-hidden rounded-2xl bg-[linear-gradient(135deg,var(--primary),var(--secondary)_58%,var(--accent))] text-white shadow-lg shadow-sky-500/20">
        <Dna size={22} strokeWidth={2.2} />
        <span className="absolute inset-0 bg-[linear-gradient(120deg,transparent_20%,rgba(255,255,255,.28)_50%,transparent_80%)] opacity-60 [transform:translateX(-100%)] [animation:shimmer_3.6s_infinite]" />
      </div>
      {!compact ? (
        <div className="min-w-0">
          <p className="truncate text-[15px] font-bold tracking-tight text-[var(--text)]">
            Genome Embeddings
          </p>
          <p className="truncate text-[11px] text-[var(--text-muted)]">
            Turning genomes into mathematics.
          </p>
        </div>
      ) : null}
    </div>
  )
}

export function AppShell({
  page,
  setPage,
  children,
  analysisReady,
  openSetup,
}: {
  page: PageId
  setPage: (page: PageId) => void
  children: ReactNode
  analysisReady: boolean
  openSetup: () => void
}) {
  const [collapsed, setCollapsed] = useState(false)
  const [mobileOpen, setMobileOpen] = useState(false)
  const { theme, toggleTheme } = useThemeContext()

  const selectPage = (next: PageId) => {
    setPage(next)
    setMobileOpen(false)
  }

  const sidebar = (mobile = false) => (
    <div className="flex h-full flex-col">
      <div className={cn('flex h-20 items-center border-b border-[var(--border)] px-4', collapsed && !mobile ? 'justify-center' : 'justify-between')}>
        <Brand compact={collapsed && !mobile} />
        {!mobile && !collapsed ? (
          <button
            type="button"
            className="focus-ring rounded-xl p-2 text-[var(--text-muted)] transition hover:bg-[var(--surface-hover)] hover:text-[var(--text)]"
            onClick={() => setCollapsed(true)}
            aria-label="Collapse sidebar"
          >
            <PanelLeftClose size={18} />
          </button>
        ) : null}
        {mobile ? (
          <button
            type="button"
            onClick={() => setMobileOpen(false)}
            className="focus-ring rounded-xl p-2 text-[var(--text-muted)] hover:bg-[var(--surface-hover)]"
            aria-label="Close navigation"
          >
            <X size={19} />
          </button>
        ) : null}
      </div>

      <nav className="flex-1 space-y-1.5 overflow-y-auto px-3 py-5" aria-label="Main navigation">
        {navigation.map((item) => {
          const Icon = item.icon
          const active = page === item.id
          return (
            <button
              key={item.id}
              type="button"
              onClick={() => selectPage(item.id)}
              title={collapsed && !mobile ? item.label : undefined}
              className={cn(
                'focus-ring group relative flex w-full items-center rounded-2xl px-3 py-3 text-left transition-all duration-200',
                collapsed && !mobile ? 'justify-center' : 'gap-3.5',
                active
                  ? 'bg-[linear-gradient(110deg,color-mix(in_srgb,var(--primary)_16%,transparent),color-mix(in_srgb,var(--secondary)_8%,transparent))] text-[var(--text)] shadow-[inset_0_0_0_1px_var(--border-strong)]'
                  : 'text-[var(--text-soft)] hover:bg-[var(--surface-hover)] hover:text-[var(--text)]',
              )}
            >
              {active ? (
                <motion.span
                  layoutId={mobile ? 'mobile-nav-indicator' : 'desktop-nav-indicator'}
                  className="absolute left-0 h-7 w-1 rounded-r-full bg-[var(--primary)]"
                  transition={{ type: 'spring', stiffness: 380, damping: 32 }}
                />
              ) : null}
              <span
                className={cn(
                  'grid h-9 w-9 shrink-0 place-items-center rounded-xl transition-colors',
                  active
                    ? 'bg-[color-mix(in_srgb,var(--primary)_15%,transparent)] text-[var(--primary)]'
                    : 'bg-[var(--bg-subtle)] text-[var(--text-muted)] group-hover:text-[var(--primary)]',
                )}
              >
                <Icon size={18} />
              </span>
              {!(collapsed && !mobile) ? (
                <span className="min-w-0 flex-1">
                  <span className="block truncate text-sm font-semibold">{item.label}</span>
                  <span className="mt-0.5 block truncate text-[11px] text-[var(--text-muted)]">
                    {item.description}
                  </span>
                </span>
              ) : null}
            </button>
          )
        })}
      </nav>

      <div className="border-t border-[var(--border)] p-3">
        <button
          type="button"
          onClick={openSetup}
          className={cn(
            'focus-ring flex w-full items-center rounded-2xl bg-[var(--bg-subtle)] px-3 py-3 text-[var(--text-soft)] transition hover:bg-[var(--surface-hover)] hover:text-[var(--text)]',
            collapsed && !mobile ? 'justify-center' : 'gap-3',
          )}
        >
          <Settings2 size={18} className="text-[var(--primary)]" />
          {!(collapsed && !mobile) ? (
            <span className="text-sm font-semibold">Analysis setup</span>
          ) : null}
        </button>
        {!mobile && collapsed ? (
          <button
            type="button"
            onClick={() => setCollapsed(false)}
            className="focus-ring mt-2 flex w-full justify-center rounded-xl p-2 text-[var(--text-muted)] hover:bg-[var(--surface-hover)]"
            aria-label="Expand sidebar"
          >
            <PanelLeftOpen size={18} />
          </button>
        ) : null}
      </div>
    </div>
  )

  return (
    <div className="min-h-screen">
      <aside
        className={cn(
          'glass-panel fixed inset-y-0 left-0 z-40 hidden rounded-none border-y-0 border-l-0 transition-[width] duration-300 lg:block',
          collapsed ? 'w-[5.25rem]' : 'w-[17.5rem]',
        )}
      >
        {sidebar()}
      </aside>

      <AnimatePresence>
        {mobileOpen ? (
          <>
            <motion.button
              type="button"
              aria-label="Close navigation overlay"
              className="fixed inset-0 z-40 bg-black/45 backdrop-blur-sm lg:hidden"
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              exit={{ opacity: 0 }}
              onClick={() => setMobileOpen(false)}
            />
            <motion.aside
              className="glass-panel fixed inset-y-0 left-0 z-50 w-[min(21rem,88vw)] rounded-none border-y-0 border-l-0 lg:hidden"
              initial={{ x: '-100%' }}
              animate={{ x: 0 }}
              exit={{ x: '-100%' }}
              transition={{ type: 'spring', stiffness: 360, damping: 36 }}
            >
              {sidebar(true)}
            </motion.aside>
          </>
        ) : null}
      </AnimatePresence>

      <div className={cn('transition-[padding] duration-300', collapsed ? 'lg:pl-[5.25rem]' : 'lg:pl-[17.5rem]')}>
        <header className="sticky top-0 z-30 border-b border-[var(--border)] bg-[color-mix(in_srgb,var(--bg)_78%,transparent)] backdrop-blur-2xl">
          <div className="flex h-20 items-center justify-between gap-4 px-4 sm:px-6 xl:px-8">
            <div className="flex min-w-0 items-center gap-3">
              <button
                type="button"
                onClick={() => setMobileOpen(true)}
                className="focus-ring rounded-xl p-2.5 text-[var(--text-soft)] hover:bg-[var(--surface-hover)] lg:hidden"
                aria-label="Open navigation"
              >
                <Menu size={20} />
              </button>
              <div className="lg:hidden">
                <Brand />
              </div>
              <div className="hidden min-w-0 lg:block">
                <div className="flex items-center gap-2 text-xs font-medium text-[var(--text-muted)]">
                  <span>Workspace</span>
                  <ChevronRight size={13} />
                  <span className="text-[var(--primary)]">
                    {navigation.find((item) => item.id === page)?.label}
                  </span>
                </div>
                <p className="mt-1 truncate text-lg font-semibold text-[var(--text)]">
                  {analysisReady ? 'Interactive sequence analysis' : 'Preparing scientific workspace'}
                </p>
              </div>
            </div>

            <div className="flex items-center gap-2 sm:gap-3">
              <span className="hidden items-center gap-2 rounded-full border border-[var(--border)] bg-[var(--bg-solid)]/65 px-3 py-2 text-xs font-medium text-[var(--text-soft)] sm:flex">
                <span className={cn('h-2 w-2 rounded-full', analysisReady ? 'bg-emerald-400 shadow-[0_0_12px_rgba(52,211,153,.65)]' : 'bg-amber-400')} />
                {analysisReady ? 'Analysis ready' : 'Connecting'}
              </span>
              <button
                type="button"
                onClick={toggleTheme}
                className="focus-ring relative grid h-11 w-11 place-items-center overflow-hidden rounded-2xl border border-[var(--border)] bg-[var(--bg-solid)]/75 text-[var(--text-soft)] transition hover:border-[var(--border-strong)] hover:text-[var(--primary)]"
                aria-label={`Switch to ${theme === 'dark' ? 'light' : 'dark'} theme`}
              >
                <AnimatePresence mode="wait" initial={false}>
                  <motion.span
                    key={theme}
                    initial={{ y: 16, opacity: 0, rotate: -35 }}
                    animate={{ y: 0, opacity: 1, rotate: 0 }}
                    exit={{ y: -16, opacity: 0, rotate: 35 }}
                    transition={{ duration: 0.2 }}
                  >
                    {theme === 'dark' ? <Sun size={19} /> : <Moon size={19} />}
                  </motion.span>
                </AnimatePresence>
              </button>
              <button
                type="button"
                onClick={openSetup}
                className="focus-ring hidden min-h-11 items-center gap-2 rounded-2xl bg-[var(--primary)] px-4 text-sm font-semibold text-white shadow-lg shadow-sky-500/20 transition hover:brightness-110 sm:flex"
              >
                <Settings2 size={17} />
                Configure
              </button>
            </div>
          </div>
        </header>

        <main className="mx-auto w-full max-w-[104rem] px-4 py-6 sm:px-6 sm:py-8 xl:px-8">
          <AnimatePresence mode="wait">
            <motion.div
              key={page}
              initial={{ opacity: 0, y: 12 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -8 }}
              transition={{ duration: 0.28, ease: [0.22, 1, 0.36, 1] }}
            >
              {children}
            </motion.div>
          </AnimatePresence>
        </main>
      </div>
    </div>
  )
}
