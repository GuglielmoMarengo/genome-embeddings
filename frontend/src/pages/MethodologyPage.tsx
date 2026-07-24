import { Atom, Binary, BookOpen, BrainCircuit, CircleDot, Layers3, ShieldAlert } from 'lucide-react'
import { Badge, Card, SectionHeading } from '../components/ui'

const methods = [
  {
    title: 'Legacy descriptor baseline',
    icon: Atom,
    tone: 'sky',
    body: 'A compact six-feature normalized vector covering composition, entropy, skew, purine content, k-mer diversity, and normalized k-mer entropy.',
  },
  {
    title: 'Descriptor Foundation V2',
    icon: CircleDot,
    tone: 'teal',
    body: 'Adds finite-sample entropy normalization, effective k-mer counts, coverage, repeat diagnostics, conditional entropy, and dinucleotide dependency ratios.',
  },
  {
    title: 'Multiscale embedding V2',
    icon: Layers3,
    tone: 'violet',
    body: 'Represents global composition and dependency features once, then appends only scale-specific diagnostic blocks for each selected k.',
  },
  {
    title: 'Jensen–Shannon distance',
    icon: Binary,
    tone: 'amber',
    body: 'Compares complete normalized k-mer probability vectors. The square root of the Jensen–Shannon divergence is used as a bounded symmetric distance.',
  },
]

export function MethodologyPage() {
  const toneClasses: Record<string, string> = {
    sky: 'bg-sky-500/10 text-sky-500',
    teal: 'bg-teal-500/10 text-teal-500',
    violet: 'bg-violet-500/10 text-violet-500',
    amber: 'bg-amber-500/10 text-amber-500',
  }

  return (
    <div className="space-y-7">
      <SectionHeading
        eyebrow="Transparent mathematics"
        title="Methodology"
        description="Every result in the application maps to an explicit descriptor, distribution, matrix coordinate, or rank statistic."
      />

      <div className="grid gap-5 md:grid-cols-2">
        {methods.map(({ title, icon: Icon, tone, body }) => (
          <Card key={title} className="p-6">
            <span className={`grid h-12 w-12 place-items-center rounded-2xl ${toneClasses[tone]}`}><Icon size={22} /></span>
            <h3 className="mt-5 text-lg font-semibold">{title}</h3>
            <p className="mt-3 text-sm leading-7 text-[var(--text-soft)]">{body}</p>
          </Card>
        ))}
      </div>

      <div className="grid gap-6 xl:grid-cols-[1.1fr_.9fr]">
        <Card className="p-6 sm:p-7">
          <div className="flex items-start gap-4">
            <span className="grid h-12 w-12 shrink-0 place-items-center rounded-2xl bg-sky-500/10 text-sky-500"><BookOpen size={22} /></span>
            <div>
              <h3 className="text-lg font-semibold">Interpretation hierarchy</h3>
              <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                The representations answer different questions and should not be treated as interchangeable scores.
              </p>
            </div>
          </div>
          <div className="mt-6 space-y-3">
            {[
              ['Composition', 'Which broad nucleotide balances characterize the sequence?'],
              ['Dependency', 'How predictable is the next base and which dinucleotides are enriched?'],
              ['Complexity', 'How many effective k-mer categories are supported at the selected scale?'],
              ['Distribution', 'Which exact k-mers are used and with what probabilities?'],
              ['Geometry', 'How do all pairwise relationships change across scales?'],
            ].map(([label, description], index) => (
              <div key={label} className="grid gap-3 rounded-2xl border border-[var(--border)] bg-[var(--bg-solid)]/55 p-4 sm:grid-cols-[2.5rem_8rem_1fr] sm:items-center">
                <span className="grid h-9 w-9 place-items-center rounded-xl bg-[var(--bg-subtle)] text-sm font-semibold text-[var(--primary)]">{index + 1}</span>
                <span className="font-semibold">{label}</span>
                <span className="text-sm leading-6 text-[var(--text-soft)]">{description}</span>
              </div>
            ))}
          </div>
        </Card>

        <div className="space-y-5">
          <Card className="p-6">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-violet-500/10 text-violet-500"><BrainCircuit size={20} /></span>
              <div>
                <h3 className="font-semibold">Interpretability before complexity</h3>
                <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                  The current project deliberately avoids learned latent vectors. Future machine-learning integration should preserve descriptor provenance and ablation testing.
                </p>
              </div>
            </div>
          </Card>
          <Card className="border-amber-500/20 bg-amber-500/[0.055] p-6">
            <div className="flex items-start gap-4">
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-2xl bg-amber-500/12 text-amber-500"><ShieldAlert size={20} /></span>
              <div>
                <div className="flex flex-wrap items-center gap-2"><h3 className="font-semibold">Research limitation</h3><Badge tone="warning">Not validated</Badge></div>
                <p className="mt-2 text-sm leading-6 text-[var(--text-soft)]">
                  Distances and rankings are properties of the implemented descriptor spaces. They are not direct measurements of sequence identity, evolutionary distance, protein function, phenotype, or clinical relevance.
                </p>
              </div>
            </div>
          </Card>
        </div>
      </div>
    </div>
  )
}
