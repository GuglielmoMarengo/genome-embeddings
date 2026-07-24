import type { HTMLAttributes, ReactNode } from 'react'
import { motion, type HTMLMotionProps } from 'motion/react'
import type { LucideIcon } from 'lucide-react'

export function cn(...values: Array<string | false | null | undefined>): string {
  return values.filter(Boolean).join(' ')
}

interface ButtonProps extends HTMLMotionProps<'button'> {
  variant?: 'primary' | 'secondary' | 'ghost' | 'danger'
  icon?: LucideIcon
  loading?: boolean
}

export function Button({
  className,
  variant = 'primary',
  icon: Icon,
  loading = false,
  children,
  disabled,
  ...props
}: ButtonProps) {
  const styles = {
    primary:
      'bg-[var(--primary)] text-white shadow-[0_12px_30px_color-mix(in_srgb,var(--primary)_28%,transparent)] hover:brightness-110',
    secondary:
      'border border-[var(--border-strong)] bg-[var(--surface-hover)] text-[var(--primary)] hover:bg-[color-mix(in_srgb,var(--primary)_13%,transparent)]',
    ghost:
      'text-[var(--text-soft)] hover:bg-[var(--surface-hover)] hover:text-[var(--text)]',
    danger:
      'border border-[color-mix(in_srgb,var(--danger)_28%,transparent)] bg-[color-mix(in_srgb,var(--danger)_10%,transparent)] text-[var(--danger)] hover:bg-[color-mix(in_srgb,var(--danger)_17%,transparent)]',
  }

  return (
    <motion.button
      whileHover={disabled || loading ? undefined : { y: -1 }}
      whileTap={disabled || loading ? undefined : { scale: 0.98 }}
      className={cn(
        'focus-ring inline-flex min-h-10 items-center justify-center gap-2 rounded-xl px-4 py-2 text-sm font-semibold transition-all duration-200 disabled:cursor-not-allowed disabled:opacity-45',
        styles[variant],
        className,
      )}
      disabled={disabled || loading}
      {...props}
    >
      {loading ? (
        <span className="h-4 w-4 animate-spin rounded-full border-2 border-current border-t-transparent" />
      ) : Icon ? (
        <Icon size={17} strokeWidth={2} />
      ) : null}
      {children}
    </motion.button>
  )
}

export function Card({ className, children, ...props }: HTMLAttributes<HTMLDivElement>) {
  return (
    <div
      className={cn('glass-panel rounded-[1.4rem]', className)}
      {...props}
    >
      {children}
    </div>
  )
}

interface MetricCardProps {
  label: string
  value: string | number
  detail?: string
  icon: LucideIcon
  accent?: 'blue' | 'teal' | 'violet' | 'amber'
}

export function MetricCard({
  label,
  value,
  detail,
  icon: Icon,
  accent = 'blue',
}: MetricCardProps) {
  const accents = {
    blue: 'text-sky-500 bg-sky-500/10',
    teal: 'text-teal-500 bg-teal-500/10',
    violet: 'text-violet-500 bg-violet-500/10',
    amber: 'text-amber-500 bg-amber-500/10',
  }

  return (
    <motion.div
      layout
      whileHover={{ y: -4 }}
      transition={{ type: 'spring', stiffness: 360, damping: 28 }}
      className="metric-glow glass-panel rounded-[1.35rem] p-5"
    >
      <div className="flex items-start justify-between gap-4">
        <div>
          <p className="text-xs font-semibold uppercase tracking-[0.16em] text-[var(--text-muted)]">
            {label}
          </p>
          <p className="mt-3 text-3xl font-semibold tracking-tight text-[var(--text)]">
            {value}
          </p>
          {detail ? (
            <p className="mt-2 text-sm text-[var(--text-muted)]">{detail}</p>
          ) : null}
        </div>
        <span className={cn('grid h-11 w-11 place-items-center rounded-2xl', accents[accent])}>
          <Icon size={21} />
        </span>
      </div>
    </motion.div>
  )
}

export function SectionHeading({
  eyebrow,
  title,
  description,
  action,
}: {
  eyebrow?: string
  title: string
  description?: string
  action?: ReactNode
}) {
  return (
    <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
      <div className="max-w-3xl">
        {eyebrow ? (
          <p className="mb-2 text-xs font-semibold uppercase tracking-[0.19em] text-[var(--primary)]">
            {eyebrow}
          </p>
        ) : null}
        <h2 className="text-2xl font-semibold tracking-tight text-[var(--text)] sm:text-3xl">
          {title}
        </h2>
        {description ? (
          <p className="mt-2 text-sm leading-6 text-[var(--text-soft)] sm:text-base">
            {description}
          </p>
        ) : null}
      </div>
      {action}
    </div>
  )
}

export function Badge({
  children,
  tone = 'neutral',
}: {
  children: ReactNode
  tone?: 'neutral' | 'success' | 'info' | 'warning'
}) {
  const styles = {
    neutral: 'bg-[var(--bg-subtle)] text-[var(--text-soft)] border-[var(--border)]',
    success:
      'bg-emerald-500/10 text-emerald-500 border-emerald-500/20',
    info: 'bg-sky-500/10 text-sky-500 border-sky-500/20',
    warning: 'bg-amber-500/10 text-amber-500 border-amber-500/20',
  }
  return (
    <span className={cn('inline-flex items-center rounded-full border px-2.5 py-1 text-xs font-medium', styles[tone])}>
      {children}
    </span>
  )
}

export function EmptyState({
  icon: Icon,
  title,
  description,
  action,
}: {
  icon: LucideIcon
  title: string
  description: string
  action?: ReactNode
}) {
  return (
    <Card className="grid min-h-72 place-items-center p-8 text-center">
      <div className="max-w-md">
        <span className="mx-auto grid h-14 w-14 place-items-center rounded-2xl bg-[var(--surface-hover)] text-[var(--primary)]">
          <Icon size={25} />
        </span>
        <h3 className="mt-5 text-lg font-semibold">{title}</h3>
        <p className="mt-2 text-sm leading-6 text-[var(--text-muted)]">{description}</p>
        {action ? <div className="mt-5">{action}</div> : null}
      </div>
    </Card>
  )
}
