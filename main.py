"""Command-line demonstration for Genome Embeddings.

The report preserves the legacy descriptor analysis and adds Descriptor
Foundation V2, a first interpretable multiscale embedding, and full k-mer
Jensen-Shannon distribution comparisons.
"""

from __future__ import annotations

from pathlib import Path

from src.dashboard import DashboardConfig, analyze_records, load_demo_records
from src.descriptor_v2 import (
    DINUCLEOTIDE_ORDER,
    DescriptorV2Collection,
    GenomeDescriptorV2,
)
from src.genome import Genome, GenomeCollection, GenomeMatrix
from src.kmer_distribution import KmerDistributionCollection
from src.visualization import (
    plot_matrix_distribution,
    plot_matrix_heatmap,
    plot_pair_trajectory,
    plot_ranking_stability,
    plot_trajectory_distributions,
    save_figure,
)


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"

FULL_GENOME_LABELS = [
    "Aequorea GFP",
    "Acropora GFP",
    "Discosoma FP583",
    "S. aureus catA",
    "S. cerevisiae TPI1",
    "Periodic control",
]
PERIODIC_CONTROL_LABEL = "Periodic control"
BIOLOGICAL_GENOME_LABELS = [
    label for label in FULL_GENOME_LABELS if label != PERIODIC_CONTROL_LABEL
]

DEFAULT_KMER_LENGTH = 3
KMER_SENSITIVITY_LENGTHS = [1, 2, 3, 4, 5]
DEFAULT_KMER_LIMIT = 10
REFERENCE_LABEL = "Aequorea GFP"
COMPARISON_LABEL = "Acropora GFP"
REPORT_WIDTH = 96


def report_title(title: str, subtitle: str | None = None) -> None:
    print("=" * REPORT_WIDTH)
    print(title.upper())
    if subtitle:
        print(subtitle)
    print("=" * REPORT_WIDTH)


def section(number: int, title: str) -> None:
    print()
    print(f"{number}. {title}".upper())
    print("-" * REPORT_WIDTH)


def subsection(title: str) -> None:
    print()
    print(title)
    print("." * min(len(title), REPORT_WIDTH))


def key_value(label: str, value: object, indent: int = 2) -> None:
    print(f"{' ' * indent}{label}: {value}")


def load_genomes() -> list[Genome]:
    return [record.genome for record in load_demo_records(PROJECT_ROOT)]


def print_genome_summary(genome: Genome) -> None:
    key_value("Header", genome.header)
    key_value("Length", f"{genome.length()} bp")
    key_value("Sequence preview", f"{genome.sequence[:64]}...")
    key_value(
        "Reverse-complement preview",
        f"{genome.reverse_complement()[:64]}...",
    )


def print_legacy_descriptor(genome: Genome, k: int) -> None:
    descriptor = genome.descriptor(k)
    rows = [
        ("GC content", f"{descriptor.gc_content * 100:.2f}%"),
        ("AT content", f"{descriptor.at_content * 100:.2f}%"),
        ("Shannon entropy", f"{descriptor.shannon_entropy:.4f} bits"),
        ("GC skew", f"{descriptor.gc_skew:.4f}"),
        ("Purine content", f"{descriptor.purine_content * 100:.2f}%"),
        ("k-mer diversity", f"{descriptor.kmer_diversity:.4f}"),
        ("k-mer entropy", f"{descriptor.kmer_entropy:.4f} bits"),
    ]
    for label, value in rows:
        key_value(label, value)


def print_v2_descriptor(descriptor: GenomeDescriptorV2) -> None:
    kmer = descriptor.kmer
    rows = [
        (
            "Conditional nucleotide entropy",
            f"{descriptor.conditional_nucleotide_entropy:.4f} bits",
        ),
        ("k-mer windows", kmer.window_count),
        ("Theoretical k-mer space", kmer.possible_kmer_count),
        ("Observable maximum", kmer.observable_kmer_count),
        ("Distinct observed k-mers", kmer.distinct_kmer_count),
        (
            "Finite-sample normalized entropy",
            f"{kmer.finite_sample_normalized_kmer_entropy:.4f}",
        ),
        ("Effective k-mer count", f"{kmer.effective_kmer_count:.2f}"),
        (
            "Theoretical-space coverage",
            f"{kmer.theoretical_space_coverage:.4f}",
        ),
        (
            "Observable-space coverage",
            f"{kmer.observable_space_coverage:.4f}",
        ),
        ("Singleton fraction", f"{kmer.singleton_fraction:.4f}"),
        (
            "Repeated-window fraction",
            f"{kmer.repeated_window_fraction:.4f}",
        ),
    ]
    for label, value in rows:
        key_value(label, value)

    print("  Strongest dinucleotide enrichments:")
    enriched = sorted(
        descriptor.dinucleotide_odds_ratios.items(),
        key=lambda item: (-item[1], item[0]),
    )[:5]
    for dinucleotide, ratio in enriched:
        print(f"    - {dinucleotide}: {ratio:.4f}")


def print_top_kmers(genome: Genome, k: int, limit: int) -> None:
    frequencies = sorted(
        genome.kmer_frequencies(k).items(),
        key=lambda item: (-item[1], item[0]),
    )
    for index, (kmer, count) in enumerate(frequencies[:limit], start=1):
        print(f"  {index:>2}. {kmer}: {count} occurrences")


def print_matrix(matrix: GenomeMatrix) -> None:
    label_width = max(len(label) for label in matrix.labels)
    value_width = max(16, max(len(label) + 2 for label in matrix.labels))
    header = " " * (label_width + 2) + "".join(
        f"{label:>{value_width}}" for label in matrix.labels
    )
    print(header)
    for label, row in zip(matrix.labels, matrix.values, strict=True):
        values = "".join(f"{value:>{value_width}.4f}" for value in row)
        print(f"{label:<{label_width}}  {values}")


def print_ranking(matrix: GenomeMatrix, reference_label: str) -> None:
    for index, (label, value) in enumerate(
        matrix.rank_by_label(reference_label),
        start=1,
    ):
        print(f"  {index}. {label}: {value:.6f}")


def print_pair_comparison(
    legacy_collection: GenomeCollection,
    v2_collection: DescriptorV2Collection,
    distribution_collection: KmerDistributionCollection,
    labels: list[str],
) -> None:
    legacy_euclidean = legacy_collection.euclidean_distance_matrix(
        labels, DEFAULT_KMER_LENGTH
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)
    legacy_cosine = legacy_collection.cosine_similarity_matrix(
        labels, DEFAULT_KMER_LENGTH
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)
    v2_euclidean = v2_collection.euclidean_distance_matrix(
        labels, DEFAULT_KMER_LENGTH
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)
    v2_cosine = v2_collection.cosine_similarity_matrix(
        labels, DEFAULT_KMER_LENGTH
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)
    embedding_distance = v2_collection.multiscale_embedding_distance_matrix(
        labels,
        KMER_SENSITIVITY_LENGTHS,
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)
    jensen_shannon = distribution_collection.distance_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    ).get_value(REFERENCE_LABEL, COMPARISON_LABEL)

    rows = [
        ("Legacy Euclidean", legacy_euclidean),
        ("Legacy cosine", legacy_cosine),
        ("Descriptor V2 Euclidean", v2_euclidean),
        ("Descriptor V2 cosine", v2_cosine),
        ("Embedding V2 Euclidean", embedding_distance),
        ("Jensen-Shannon distance", jensen_shannon),
    ]
    for label, value in rows:
        key_value(label, f"{value:.6f}")


def print_step_comparison(
    full: dict[tuple[int, int], float],
    biological: dict[tuple[int, int], float],
) -> None:
    print("  Transition     Full       Biological    Reduction")
    for transition, full_value in full.items():
        biological_value = biological[transition]
        reduction = (
            (full_value - biological_value) / full_value
            if full_value != 0.0
            else 0.0
        )
        first_k, second_k = transition
        print(
            f"  {first_k}->{second_k:<8} "
            f"{full_value:>9.6f}   "
            f"{biological_value:>9.6f}   "
            f"{reduction * 100:>8.2f}%"
        )


def print_partition(
    partition: dict[tuple[int, int], dict[str, float | int]],
) -> None:
    print("  Transition   Control share   Biological share")
    for (first_k, second_k), row in partition.items():
        print(
            f"  {first_k}->{second_k:<8} "
            f"{float(row['selected_share']) * 100:>11.2f}%   "
            f"{float(row['remaining_share']) * 100:>14.2f}%"
        )


def print_trajectory(trajectory: dict[int, float]) -> None:
    for k, value in trajectory.items():
        print(f"  k={k}: {value:.6f}")


def print_step_distances(distances: dict[tuple[int, int], float]) -> None:
    for (first_k, second_k), value in distances.items():
        print(f"  k={first_k} -> k={second_k}: {value:.6f}")


def print_ranking_trajectory(
    rankings: dict[int, list[tuple[str, float]]],
) -> None:
    for k, ranking in rankings.items():
        print(f"  k={k}: {' > '.join(label for label, _ in ranking)}")


def print_ranking_stability(
    stability: dict[tuple[int, int], dict[str, float | int | bool]],
) -> None:
    print("  Transition   Kendall tau   Inversions   Mean shift   Max shift")
    for (first_k, second_k), row in stability.items():
        print(
            f"  {first_k}->{second_k:<8} "
            f"{float(row['kendall_tau']):>10.3f}   "
            f"{int(row['discordant_pairs']):>10}   "
            f"{float(row['mean_absolute_rank_shift']):>10.3f}   "
            f"{int(row['max_rank_shift']):>9}"
        )


def matrix_trajectory(
    matrices: dict[int, GenomeMatrix],
) -> dict[int, list[float]]:
    return {
        k: matrix.to_upper_triangle_vector()
        for k, matrix in matrices.items()
    }


def save_visualizations(
    *,
    legacy_euclidean: GenomeMatrix,
    legacy_cosine: GenomeMatrix,
    descriptor_v2_euclidean: GenomeMatrix,
    embedding_v2_euclidean: GenomeMatrix,
    jensen_shannon: GenomeMatrix,
    legacy_euclidean_trajectory: dict[int, list[float]],
    legacy_cosine_trajectory: dict[int, list[float]],
    jensen_shannon_trajectory: dict[int, list[float]],
    legacy_pair_euclidean: dict[int, float],
    legacy_pair_cosine: dict[int, float],
    jensen_shannon_pair: dict[int, float],
    legacy_euclidean_stability: dict,
    legacy_cosine_stability: dict,
    jensen_shannon_stability: dict,
) -> list[Path]:
    figures = [
        ("euclidean_heatmap.png", plot_matrix_heatmap(legacy_euclidean)),
        ("cosine_heatmap.png", plot_matrix_heatmap(legacy_cosine)),
        (
            "descriptor_v2_euclidean_heatmap.png",
            plot_matrix_heatmap(descriptor_v2_euclidean),
        ),
        (
            "embedding_v2_euclidean_heatmap.png",
            plot_matrix_heatmap(embedding_v2_euclidean),
        ),
        (
            "jensen_shannon_heatmap.png",
            plot_matrix_heatmap(jensen_shannon),
        ),
        (
            "euclidean_distribution.png",
            plot_matrix_distribution(legacy_euclidean),
        ),
        (
            "cosine_distribution.png",
            plot_matrix_distribution(legacy_cosine),
        ),
        (
            "jensen_shannon_distribution.png",
            plot_matrix_distribution(jensen_shannon),
        ),
        (
            "aequorea_acropora_euclidean_trajectory.png",
            plot_pair_trajectory(
                legacy_pair_euclidean,
                REFERENCE_LABEL,
                COMPARISON_LABEL,
                "euclidean",
            ),
        ),
        (
            "aequorea_acropora_cosine_trajectory.png",
            plot_pair_trajectory(
                legacy_pair_cosine,
                REFERENCE_LABEL,
                COMPARISON_LABEL,
                "cosine",
            ),
        ),
        (
            "aequorea_acropora_jensen_shannon_trajectory.png",
            plot_pair_trajectory(
                jensen_shannon_pair,
                REFERENCE_LABEL,
                COMPARISON_LABEL,
                "jensen_shannon",
            ),
        ),
        (
            "euclidean_multi_k_distribution.png",
            plot_trajectory_distributions(
                legacy_euclidean_trajectory,
                "euclidean",
            ),
        ),
        (
            "cosine_multi_k_distribution.png",
            plot_trajectory_distributions(
                legacy_cosine_trajectory,
                "cosine",
            ),
        ),
        (
            "jensen_shannon_multi_k_distribution.png",
            plot_trajectory_distributions(
                jensen_shannon_trajectory,
                "jensen_shannon",
            ),
        ),
        (
            "euclidean_ranking_stability.png",
            plot_ranking_stability(legacy_euclidean_stability, "euclidean"),
        ),
        (
            "cosine_ranking_stability.png",
            plot_ranking_stability(legacy_cosine_stability, "cosine"),
        ),
        (
            "jensen_shannon_ranking_stability.png",
            plot_ranking_stability(
                jensen_shannon_stability,
                "jensen_shannon",
            ),
        ),
    ]

    paths = []
    for filename, figure in figures:
        paths.append(
            save_figure(
                figure,
                OUTPUT_DIR / filename,
                close=True,
            )
        )
    return paths


def main() -> None:
    records = load_demo_records(PROJECT_ROOT)
    genomes = [record.genome for record in records]
    labels = [record.label for record in records]
    biological_genomes = [
        record.genome
        for record in records
        if record.label != PERIODIC_CONTROL_LABEL
    ]

    legacy = GenomeCollection(genomes)
    biological_legacy = GenomeCollection(biological_genomes)
    descriptor_v2 = DescriptorV2Collection(genomes)
    distribution = KmerDistributionCollection(genomes)

    reference = genomes[0]
    reference_v2 = GenomeDescriptorV2.from_genome(
        reference,
        DEFAULT_KMER_LENGTH,
    )

    legacy_euclidean = legacy.euclidean_distance_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    )
    legacy_cosine = legacy.cosine_similarity_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    )
    descriptor_v2_euclidean = descriptor_v2.euclidean_distance_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    )
    descriptor_v2_cosine = descriptor_v2.cosine_similarity_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    )
    embedding_v2_euclidean = (
        descriptor_v2.multiscale_embedding_distance_matrix(
            labels,
            KMER_SENSITIVITY_LENGTHS,
        )
    )
    jensen_shannon = distribution.distance_matrix(
        labels,
        DEFAULT_KMER_LENGTH,
    )

    legacy_euclidean_matrices = legacy.euclidean_distance_matrices(
        labels,
        KMER_SENSITIVITY_LENGTHS,
    )
    legacy_cosine_matrices = legacy.cosine_similarity_matrices(
        labels,
        KMER_SENSITIVITY_LENGTHS,
    )
    js_matrices = distribution.distance_matrices(
        labels,
        KMER_SENSITIVITY_LENGTHS,
    )
    legacy_euclidean_trajectory = matrix_trajectory(legacy_euclidean_matrices)
    legacy_cosine_trajectory = matrix_trajectory(legacy_cosine_matrices)
    js_trajectory = matrix_trajectory(js_matrices)

    biological_euclidean_trajectory = biological_legacy.euclidean_matrix_trajectory(
        BIOLOGICAL_GENOME_LABELS,
        KMER_SENSITIVITY_LENGTHS,
    )
    biological_cosine_trajectory = biological_legacy.cosine_matrix_trajectory(
        BIOLOGICAL_GENOME_LABELS,
        KMER_SENSITIVITY_LENGTHS,
    )

    full_euclidean_steps = GenomeCollection.matrix_trajectory_step_distances(
        legacy_euclidean_trajectory
    )
    full_cosine_steps = GenomeCollection.matrix_trajectory_step_distances(
        legacy_cosine_trajectory
    )
    biological_euclidean_steps = (
        GenomeCollection.matrix_trajectory_step_distances(
            biological_euclidean_trajectory
        )
    )
    biological_cosine_steps = GenomeCollection.matrix_trajectory_step_distances(
        biological_cosine_trajectory
    )

    euclidean_contributions = GenomeCollection.matrix_trajectory_pair_contributions(
        labels,
        legacy_euclidean_trajectory,
    )
    cosine_contributions = GenomeCollection.matrix_trajectory_pair_contributions(
        labels,
        legacy_cosine_trajectory,
    )
    euclidean_partition = GenomeCollection.matrix_trajectory_deformation_partition(
        euclidean_contributions,
        PERIODIC_CONTROL_LABEL,
    )
    cosine_partition = GenomeCollection.matrix_trajectory_deformation_partition(
        cosine_contributions,
        PERIODIC_CONTROL_LABEL,
    )

    legacy_euclidean_rankings = GenomeCollection.matrix_ranking_trajectory(
        legacy_euclidean_matrices,
        REFERENCE_LABEL,
    )
    legacy_cosine_rankings = GenomeCollection.matrix_ranking_trajectory(
        legacy_cosine_matrices,
        REFERENCE_LABEL,
    )
    js_rankings = distribution.ranking_trajectory(
        labels,
        REFERENCE_LABEL,
        KMER_SENSITIVITY_LENGTHS,
    )
    legacy_euclidean_stability = GenomeCollection.ranking_trajectory_stability(
        legacy_euclidean_rankings
    )
    legacy_cosine_stability = GenomeCollection.ranking_trajectory_stability(
        legacy_cosine_rankings
    )
    js_stability = GenomeCollection.ranking_trajectory_stability(js_rankings)

    legacy_pair_euclidean = legacy.euclidean_pair_trajectory(
        labels,
        REFERENCE_LABEL,
        COMPARISON_LABEL,
        KMER_SENSITIVITY_LENGTHS,
    )
    legacy_pair_cosine = legacy.cosine_pair_trajectory(
        labels,
        REFERENCE_LABEL,
        COMPARISON_LABEL,
        KMER_SENSITIVITY_LENGTHS,
    )
    js_pair = distribution.pair_trajectory(
        labels,
        REFERENCE_LABEL,
        COMPARISON_LABEL,
        KMER_SENSITIVITY_LENGTHS,
    )

    output_paths = save_visualizations(
        legacy_euclidean=legacy_euclidean,
        legacy_cosine=legacy_cosine,
        descriptor_v2_euclidean=descriptor_v2_euclidean,
        embedding_v2_euclidean=embedding_v2_euclidean,
        jensen_shannon=jensen_shannon,
        legacy_euclidean_trajectory=legacy_euclidean_trajectory,
        legacy_cosine_trajectory=legacy_cosine_trajectory,
        jensen_shannon_trajectory=js_trajectory,
        legacy_pair_euclidean=legacy_pair_euclidean,
        legacy_pair_cosine=legacy_pair_cosine,
        jensen_shannon_pair=js_pair,
        legacy_euclidean_stability=legacy_euclidean_stability,
        legacy_cosine_stability=legacy_cosine_stability,
        jensen_shannon_stability=js_stability,
    )

    report_title(
        "Genome Embeddings Analysis Report",
        (
            f"Reference: {REFERENCE_LABEL} | comparison: {COMPARISON_LABEL} | "
            f"k-mer scales: {KMER_SENSITIVITY_LENGTHS}"
        ),
    )

    section(1, "Dataset")
    key_value("Full dataset", f"{len(genomes)} sequences")
    key_value("Biological-only dataset", f"{len(biological_genomes)} sequences")
    key_value("Synthetic control", PERIODIC_CONTROL_LABEL)
    subsection("Reference sequence")
    print_genome_summary(reference)

    section(2, f"Legacy descriptor baseline (k={DEFAULT_KMER_LENGTH})")
    subsection("Reference descriptor")
    print_legacy_descriptor(reference, DEFAULT_KMER_LENGTH)
    subsection(f"Top {DEFAULT_KMER_LIMIT} reference k-mers by frequency")
    print_top_kmers(reference, DEFAULT_KMER_LENGTH, DEFAULT_KMER_LIMIT)
    subsection("Legacy Euclidean ranking")
    print_ranking(legacy_euclidean, REFERENCE_LABEL)
    subsection("Legacy cosine ranking")
    print_ranking(legacy_cosine, REFERENCE_LABEL)

    section(3, "Descriptor Foundation V2")
    subsection("Reference finite-sample and dependency descriptors")
    print_v2_descriptor(reference_v2)
    subsection(f"{REFERENCE_LABEL} vs {COMPARISON_LABEL}")
    print_pair_comparison(legacy, descriptor_v2, distribution, labels)
    subsection("Descriptor V2 Euclidean ranking")
    print_ranking(descriptor_v2_euclidean, REFERENCE_LABEL)
    subsection("Multiscale embedding V2 Euclidean ranking")
    print_ranking(embedding_v2_euclidean, REFERENCE_LABEL)

    section(4, "K-mer distribution comparison")
    subsection("Jensen-Shannon distance ranking")
    print_ranking(jensen_shannon, REFERENCE_LABEL)
    subsection("Jensen-Shannon pair trajectory")
    print_trajectory(js_pair)
    subsection("Jensen-Shannon matrix step distances")
    print_step_distances(
        GenomeCollection.matrix_trajectory_step_distances(js_trajectory)
    )
    subsection("Jensen-Shannon ranking trajectory")
    print_ranking_trajectory(js_rankings)
    subsection("Jensen-Shannon ranking stability")
    print_ranking_stability(js_stability)

    section(5, "Legacy multiscale geometry and control influence")
    subsection("Euclidean full vs biological-only step distances")
    print_step_comparison(full_euclidean_steps, biological_euclidean_steps)
    subsection("Cosine full vs biological-only step distances")
    print_step_comparison(full_cosine_steps, biological_cosine_steps)
    subsection("Euclidean squared-deformation partition")
    print_partition(euclidean_partition)
    subsection("Cosine squared-deformation partition")
    print_partition(cosine_partition)

    section(6, "Selected pair trajectories")
    subsection("Legacy Euclidean")
    print_trajectory(legacy_pair_euclidean)
    subsection("Legacy cosine")
    print_trajectory(legacy_pair_cosine)
    subsection("Jensen-Shannon")
    print_trajectory(js_pair)

    section(7, "Comparison matrices")
    subsection("Legacy Euclidean distance matrix")
    print_matrix(legacy_euclidean)
    subsection("Descriptor V2 Euclidean distance matrix")
    print_matrix(descriptor_v2_euclidean)
    subsection("Multiscale embedding V2 distance matrix")
    print_matrix(embedding_v2_euclidean)
    subsection("Jensen-Shannon distance matrix")
    print_matrix(jensen_shannon)

    section(8, "Generated artifacts and dashboard")
    for path in output_paths:
        print(f"  - {path.relative_to(PROJECT_ROOT)}")
    key_value("Scientific dashboard", "python app.py")

    dashboard_analysis = analyze_records(
        records,
        DashboardConfig(
            k_values=tuple(KMER_SENSITIVITY_LENGTHS),
            selected_k=DEFAULT_KMER_LENGTH,
            reference_label=REFERENCE_LABEL,
            comparison_label=COMPARISON_LABEL,
        ),
    )
    subsection("Serialization summary")
    key_value("Dashboard JSON characters", len(dashboard_analysis.to_json()))
    key_value("Summary CSV rows", len(dashboard_analysis.summary_csv().splitlines()))

    print()
    print("=" * REPORT_WIDTH)
    print("ANALYSIS COMPLETE")
    print("=" * REPORT_WIDTH)


if __name__ == "__main__":
    main()
