export type Theme = 'light' | 'dark'

export type PageId =
  | 'overview'
  | 'dataset'
  | 'descriptors'
  | 'matrices'
  | 'multiscale'
  | 'exports'
  | 'methodology'

export interface SequenceRecord {
  label: string
  sequence: string
  header: string
  source: string
}

export interface AnalysisConfig {
  k_values: number[]
  selected_k: number
  reference_label: string
  comparison_label: string
}

export interface MatrixData {
  labels: string[]
  values: number[][]
  metric: string
  kmer_length: number
}

export interface DatasetRow {
  label: string
  source: string
  length: number
  gc_content: number
  header: string
}

export interface DescriptorRow {
  label: string
  length: number
  gc_content: number
  conditional_entropy: number
  finite_sample_entropy: number
  kmer_window_count: number
  possible_kmer_count: number
  observable_kmer_count: number
  distinct_kmer_count: number
  effective_kmer_count: number
  dinucleotide_odds_ratios: Record<string, number>
  theoretical_coverage: number
  observable_coverage: number
  singleton_fraction: number
  repeated_window_fraction: number
}

export interface SummaryData {
  dataset_size: number
  selected_k: number
  k_values: number[]
  reference_label: string
  comparison_label: string
  legacy_euclidean_distance: number
  legacy_cosine_similarity: number
  descriptor_v2_euclidean_distance: number
  descriptor_v2_cosine_similarity: number
  embedding_v2_euclidean_distance: number
  jensen_shannon_distance: number
}

export interface RankingStabilityRow {
  kendall_tau: number
  concordant_pairs: number
  discordant_pairs: number
  total_pairs: number
  mean_absolute_rank_shift: number
  max_rank_shift: number
  exact_match: boolean
}

export interface AnalysisData {
  summary: SummaryData
  dataset: DatasetRow[]
  descriptor_v2: DescriptorRow[]
  matrices: Record<string, MatrixData>
  pair_trajectories: Record<string, Record<string, number>>
  jensen_shannon_step_distances: Record<string, number>
  jensen_shannon_rankings: Record<string, [string, number][]>
  jensen_shannon_ranking_stability: Record<string, RankingStabilityRow>
}

export interface DemoResponse {
  records: SequenceRecord[]
  analysis: AnalysisData
}

export interface AnalysisRequest {
  records: SequenceRecord[]
  config: AnalysisConfig
}

export interface ToastMessage {
  id: number
  title: string
  description?: string
  tone: 'success' | 'error' | 'info'
}
