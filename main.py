from dataclasses import dataclass
from pathlib import Path

from src.genome import (
    Genome,
    GenomeCollection,
    GenomeComparison,
    GenomeDescriptor,
    GenomeMatrix,
)
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

AEQUOREA_GFP_PATH = (
    PROJECT_ROOT
    / "data"
    / "fluorescent_proteins"
    / "aequorea_victoria_gfp_cds.fasta"
)

ACROPORA_GFP_PATH = (
    PROJECT_ROOT
    / "data"
    / "fluorescent_proteins"
    / "acropora_millepora_gfp_cds.fasta"
)

DISCOSOMA_FP583_PATH = (
    PROJECT_ROOT
    / "data"
    / "fluorescent_proteins"
    / "discosoma_fp583_cds.fasta"
)

STAPHYLOCOCCUS_AUREUS_CATA_PATH = (
    PROJECT_ROOT
    / "data"
    / "controls"
    / "biological"
    / "staphylococcus_aureus_cata_cds.fasta"
)

SACCHAROMYCES_CEREVISIAE_TPI1_PATH = (
    PROJECT_ROOT
    / "data"
    / "controls"
    / "biological"
    / "saccharomyces_cerevisiae_tpi1_cds.fasta"
)

PERIODIC_CONTROL_PATH = (
    PROJECT_ROOT
    / "data"
    / "controls"
    / "periodic_sequence.fasta"
)

PERIODIC_CONTROL_LABEL = "Periodic control"

FULL_GENOME_LABELS = [
    "Aequorea GFP",
    "Acropora GFP",
    "Discosoma FP583",
    "S. aureus catA",
    "S. cerevisiae TPI1",
    PERIODIC_CONTROL_LABEL,
]

BIOLOGICAL_GENOME_LABELS = [
    label
    for label in FULL_GENOME_LABELS
    if label != PERIODIC_CONTROL_LABEL
]

DEFAULT_KMER_LENGTH = 3
DEFAULT_KMER_LIMIT = 10

KMER_SENSITIVITY_LENGTHS = [
    1,
    2,
    3,
    4,
    5,
]

REFERENCE_LABEL = "Aequorea GFP"
COMPARISON_LABEL = "Acropora GFP"
REPORT_WIDTH = 88


RankingTrajectory = dict[
    int,
    list[tuple[str, float]],
]

RankingStability = dict[
    tuple[int, int],
    dict[str, float | int | bool],
]


@dataclass(slots=True)
class DatasetMultiscaleAnalysis:
    labels: list[str]

    euclidean_trajectory: dict[
        int,
        list[float],
    ]

    cosine_trajectory: dict[
        int,
        list[float],
    ]

    euclidean_step_distances: dict[
        tuple[int, int],
        float,
    ]

    cosine_step_distances: dict[
        tuple[int, int],
        float,
    ]

    euclidean_pair_contributions: dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ]

    cosine_pair_contributions: dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ]

    euclidean_deformation_partition: dict[
        tuple[int, int],
        dict[str, float | int],
    ]

    cosine_deformation_partition: dict[
        tuple[int, int],
        dict[str, float | int],
    ]

    euclidean_ranking_trajectory: RankingTrajectory
    cosine_ranking_trajectory: RankingTrajectory
    euclidean_ranking_stability: RankingStability
    cosine_ranking_stability: RankingStability


def print_report_title(
    title: str,
    subtitle: str | None = None,
) -> None:
    print("=" * REPORT_WIDTH)
    print(title.upper())

    if subtitle is not None:
        print(subtitle)

    print("=" * REPORT_WIDTH)


def print_section(
    title: str,
) -> None:
    print()
    print(title.upper())
    print("-" * REPORT_WIDTH)


def print_subsection(
    title: str,
) -> None:
    print()
    print(title)
    print("." * min(len(title), REPORT_WIDTH))


def print_key_value(
    label: str,
    value: object,
    indent: int = 0,
) -> None:
    prefix = " " * indent
    print(f"{prefix}{label}: {value}")


def print_genome_summary(
    genome: Genome,
) -> None:
    print_key_value("Header", genome.header)
    print_key_value("Length", f"{genome.length()} bp")
    print_key_value(
        "Sequence preview",
        f"{genome.sequence[:60]}...",
    )
    print_key_value(
        "Reverse-complement preview",
        f"{genome.reverse_complement()[:60]}...",
    )


def print_descriptor(
    descriptor: GenomeDescriptor,
) -> None:
    rows = [
        ("GC content", f"{descriptor.gc_content * 100:.2f}%"),
        ("AT content", f"{descriptor.at_content * 100:.2f}%"),
        ("Shannon entropy", f"{descriptor.shannon_entropy:.4f} bits"),
        ("GC skew", f"{descriptor.gc_skew:.4f}"),
        ("Purine content", f"{descriptor.purine_content * 100:.2f}%"),
        (
            "Pyrimidine content",
            f"{descriptor.pyrimidine_content * 100:.2f}%",
        ),
        ("k-mer length", descriptor.kmer_length),
        ("k-mer diversity", f"{descriptor.kmer_diversity:.4f}"),
        ("k-mer entropy", f"{descriptor.kmer_entropy:.4f} bits"),
    ]

    for label, value in rows:
        print_key_value(label, value, indent=2)


def print_kmer_frequencies(
    genome: Genome,
    k: int = DEFAULT_KMER_LENGTH,
    limit: int = DEFAULT_KMER_LIMIT,
) -> None:
    frequencies = genome.kmer_frequencies(k)

    for index, (kmer, frequency) in enumerate(
        frequencies.items(),
        start=1,
    ):
        print(
            f"  {index:>2}. {kmer}: "
            f"{frequency} occurrences"
        )

        if index == limit:
            break


def print_genome_comparison(
    comparison: GenomeComparison,
) -> None:
    print_key_value(
        "Euclidean distance",
        f"{comparison.euclidean_distance:.6f}",
        indent=2,
    )

    print_key_value(
        "Cosine similarity",
        f"{comparison.cosine_similarity:.6f}",
        indent=2,
    )

    print("  Feature differences:")

    for feature_name, difference in (
        comparison.sorted_feature_differences()
    ):
        print(
            f"    - {feature_name}: "
            f"{difference:.6f}"
        )


def print_genome_matrix(
    matrix: GenomeMatrix,
) -> None:
    label_width = max(
        len(label)
        for label in matrix.labels
    )

    value_width = max(
        12,
        max(
            len(label) + 2
            for label in matrix.labels
        ),
    )

    header = " " * (label_width + 2)

    for label in matrix.labels:
        header += f"{label:>{value_width}}"

    print(header)

    for label, row in zip(
        matrix.labels,
        matrix.values,
        strict=True,
    ):
        formatted_values = "".join(
            f"{value:>{value_width}.4f}"
            for value in row
        )

        print(
            f"{label:<{label_width}}  "
            f"{formatted_values}"
        )


def print_ranking(
    ranking: list[tuple[str, float]],
) -> None:
    for position, (label, value) in enumerate(
        ranking,
        start=1,
    ):
        print(
            f"  {position}. {label}: "
            f"{value:.6f}"
        )


def print_dataset_dimensions(
    full_analysis: DatasetMultiscaleAnalysis,
    biological_analysis: DatasetMultiscaleAnalysis,
) -> None:
    for k, full_vector in (
        full_analysis.euclidean_trajectory.items()
    ):
        biological_vector = (
            biological_analysis
            .euclidean_trajectory[k]
        )

        print(
            f"  k={k}: "
            f"full={len(full_vector)} coordinates, "
            f"biological-only="
            f"{len(biological_vector)} coordinates"
        )


def print_dataset_step_comparison(
    full_distances: dict[
        tuple[int, int],
        float,
    ],
    biological_distances: dict[
        tuple[int, int],
        float,
    ],
) -> None:
    if (
        full_distances.keys()
        != biological_distances.keys()
    ):
        raise ValueError(
            "Dataset step-distance transitions "
            "must match."
        )

    print(
        "  Transition   Full       Biological   "
        "Reduction"
    )

    for transition, full_distance in (
        full_distances.items()
    ):
        biological_distance = (
            biological_distances[transition]
        )

        relative_reduction = (
            (
                full_distance
                - biological_distance
            )
            / full_distance
            if full_distance != 0
            else 0.0
        )

        first_k, second_k = transition

        print(
            f"  {first_k}->{second_k:<7}"
            f"{full_distance:>10.6f}"
            f"{biological_distance:>13.6f}"
            f"{relative_reduction * 100:>11.2f}%"
        )


def print_deformation_partition(
    partitions: dict[
        tuple[int, int],
        dict[str, float | int],
    ],
) -> None:
    print(
        "  Transition   Control share   "
        "Biological share"
    )

    for (first_k, second_k), partition in (
        partitions.items()
    ):
        selected_share = float(
            partition["selected_share"]
        )

        remaining_share = float(
            partition["remaining_share"]
        )

        print(
            f"  {first_k}->{second_k:<7}"
            f"{selected_share * 100:>12.2f}%"
            f"{remaining_share * 100:>17.2f}%"
        )


def print_pair_trajectory(
    trajectory: dict[int, float],
) -> None:
    for k, value in trajectory.items():
        print(
            f"  k={k}: {value:.6f}"
        )


def print_pair_step_differences(
    differences: dict[
        tuple[int, int],
        float,
    ],
) -> None:
    for (first_k, second_k), difference in (
        differences.items()
    ):
        print(
            f"  k={first_k} -> k={second_k}: "
            f"{difference:+.6f}"
        )


def print_top_pair_contributions(
    contributions: dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ],
    limit: int = 3,
) -> None:
    for (first_k, second_k), rows in (
        contributions.items()
    ):
        print(
            f"  k={first_k} -> k={second_k}"
        )

        for position, row in enumerate(
            rows[:limit],
            start=1,
        ):
            print(
                f"    {position}. "
                f"{row['row_label']} -> "
                f"{row['column_label']}: "
                f"{float(row['difference']):+.6f}"
            )


def print_ranking_trajectory(
    rankings: RankingTrajectory,
) -> None:
    for k, ranking in rankings.items():
        order = " > ".join(
            label
            for label, _ in ranking
        )

        print(
            f"  k={k}: {order}"
        )


def print_ranking_stability(
    stability: RankingStability,
) -> None:
    print(
        "  Transition   Kendall tau   "
        "Inversions   Mean shift   Max shift"
    )

    for (first_k, second_k), values in (
        stability.items()
    ):
        print(
            f"  {first_k}->{second_k:<7}"
            f"{float(values['kendall_tau']):>12.3f}"
            f"{int(values['discordant_pairs']):>13}"
            f"{float(values['mean_absolute_rank_shift']):>13.3f}"
            f"{int(values['max_rank_shift']):>12}"
        )


def load_genomes() -> list[Genome]:
    return [
        Genome.from_fasta(AEQUOREA_GFP_PATH),
        Genome.from_fasta(ACROPORA_GFP_PATH),
        Genome.from_fasta(DISCOSOMA_FP583_PATH),
        Genome.from_fasta(
            STAPHYLOCOCCUS_AUREUS_CATA_PATH
        ),
        Genome.from_fasta(
            SACCHAROMYCES_CEREVISIAE_TPI1_PATH
        ),
        Genome.from_fasta(PERIODIC_CONTROL_PATH),
    ]


def analyze_multiscale_dataset(
    collection: GenomeCollection,
    labels: list[str],
    reference_label: str,
) -> DatasetMultiscaleAnalysis:
    euclidean_trajectory = (
        collection.euclidean_matrix_trajectory(
            labels=labels,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    cosine_trajectory = (
        collection.cosine_matrix_trajectory(
            labels=labels,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    euclidean_step_distances = (
        collection.matrix_trajectory_step_distances(
            euclidean_trajectory
        )
    )

    cosine_step_distances = (
        collection.matrix_trajectory_step_distances(
            cosine_trajectory
        )
    )

    euclidean_pair_contributions = (
        collection
        .matrix_trajectory_pair_contributions(
            labels=labels,
            trajectory=euclidean_trajectory,
        )
    )

    cosine_pair_contributions = (
        collection
        .matrix_trajectory_pair_contributions(
            labels=labels,
            trajectory=cosine_trajectory,
        )
    )

    euclidean_deformation_partition = (
        collection
        .matrix_trajectory_deformation_partition(
            contributions=(
                euclidean_pair_contributions
            ),
            selected_label=PERIODIC_CONTROL_LABEL,
        )
        if PERIODIC_CONTROL_LABEL in labels
        else {}
    )

    cosine_deformation_partition = (
        collection
        .matrix_trajectory_deformation_partition(
            contributions=(
                cosine_pair_contributions
            ),
            selected_label=PERIODIC_CONTROL_LABEL,
        )
        if PERIODIC_CONTROL_LABEL in labels
        else {}
    )

    euclidean_ranking_trajectory = (
        collection.euclidean_ranking_trajectory(
            labels=labels,
            reference_label=reference_label,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    cosine_ranking_trajectory = (
        collection.cosine_ranking_trajectory(
            labels=labels,
            reference_label=reference_label,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    return DatasetMultiscaleAnalysis(
        labels=labels.copy(),
        euclidean_trajectory=(
            euclidean_trajectory
        ),
        cosine_trajectory=(
            cosine_trajectory
        ),
        euclidean_step_distances=(
            euclidean_step_distances
        ),
        cosine_step_distances=(
            cosine_step_distances
        ),
        euclidean_pair_contributions=(
            euclidean_pair_contributions
        ),
        cosine_pair_contributions=(
            cosine_pair_contributions
        ),
        euclidean_deformation_partition=(
            euclidean_deformation_partition
        ),
        cosine_deformation_partition=(
            cosine_deformation_partition
        ),
        euclidean_ranking_trajectory=(
            euclidean_ranking_trajectory
        ),
        cosine_ranking_trajectory=(
            cosine_ranking_trajectory
        ),
        euclidean_ranking_stability=(
            collection.ranking_trajectory_stability(
                euclidean_ranking_trajectory
            )
        ),
        cosine_ranking_stability=(
            collection.ranking_trajectory_stability(
                cosine_ranking_trajectory
            )
        ),
    )


def save_visualizations(
    euclidean_matrix: GenomeMatrix,
    cosine_matrix: GenomeMatrix,
    full_analysis: DatasetMultiscaleAnalysis,
    euclidean_pair_trajectory: dict[int, float],
    cosine_pair_trajectory: dict[int, float],
) -> list[Path]:
    figures = [
        (
            "euclidean_heatmap.png",
            plot_matrix_heatmap(
                euclidean_matrix
            ),
        ),
        (
            "cosine_heatmap.png",
            plot_matrix_heatmap(
                cosine_matrix
            ),
        ),
        (
            "aequorea_acropora_"
            "euclidean_trajectory.png",
            plot_pair_trajectory(
                trajectory=(
                    euclidean_pair_trajectory
                ),
                row_label=REFERENCE_LABEL,
                column_label=COMPARISON_LABEL,
                metric="euclidean",
            ),
        ),
        (
            "aequorea_acropora_"
            "cosine_trajectory.png",
            plot_pair_trajectory(
                trajectory=(
                    cosine_pair_trajectory
                ),
                row_label=REFERENCE_LABEL,
                column_label=COMPARISON_LABEL,
                metric="cosine",
            ),
        ),
        (
            "euclidean_distribution.png",
            plot_matrix_distribution(
                euclidean_matrix
            ),
        ),
        (
            "cosine_distribution.png",
            plot_matrix_distribution(
                cosine_matrix
            ),
        ),
        (
            "euclidean_multi_k_distribution.png",
            plot_trajectory_distributions(
                trajectory=(
                    full_analysis
                    .euclidean_trajectory
                ),
                metric="euclidean",
            ),
        ),
        (
            "cosine_multi_k_distribution.png",
            plot_trajectory_distributions(
                trajectory=(
                    full_analysis
                    .cosine_trajectory
                ),
                metric="cosine",
            ),
        ),
        (
            "euclidean_ranking_stability.png",
            plot_ranking_stability(
                stability=(
                    full_analysis
                    .euclidean_ranking_stability
                ),
                metric="euclidean",
            ),
        ),
        (
            "cosine_ranking_stability.png",
            plot_ranking_stability(
                stability=(
                    full_analysis
                    .cosine_ranking_stability
                ),
                metric="cosine",
            ),
        ),
    ]

    output_paths: list[Path] = []

    for filename, figure in figures:
        output_paths.append(
            save_figure(
                figure=figure,
                output_path=(
                    OUTPUT_DIR
                    / filename
                ),
                close=True,
            )
        )

    return output_paths


def main() -> None:
    genomes = load_genomes()

    biological_genomes = [
        genome
        for genome, label in zip(
            genomes,
            FULL_GENOME_LABELS,
            strict=True,
        )
        if label != PERIODIC_CONTROL_LABEL
    ]

    full_collection = GenomeCollection(genomes)

    biological_collection = GenomeCollection(
        biological_genomes
    )

    reference_genome = genomes[0]
    comparison_genome = genomes[1]

    reference_descriptor = (
        reference_genome.descriptor(
            k=DEFAULT_KMER_LENGTH
        )
    )

    comparison_descriptor = (
        comparison_genome.descriptor(
            k=DEFAULT_KMER_LENGTH
        )
    )

    comparison = reference_descriptor.compare(
        comparison_descriptor
    )

    euclidean_matrix = (
        full_collection.euclidean_distance_matrix(
            labels=FULL_GENOME_LABELS,
            k=DEFAULT_KMER_LENGTH,
        )
    )

    cosine_matrix = (
        full_collection.cosine_similarity_matrix(
            labels=FULL_GENOME_LABELS,
            k=DEFAULT_KMER_LENGTH,
        )
    )

    full_analysis = analyze_multiscale_dataset(
        collection=full_collection,
        labels=FULL_GENOME_LABELS,
        reference_label=REFERENCE_LABEL,
    )

    biological_analysis = (
        analyze_multiscale_dataset(
            collection=biological_collection,
            labels=BIOLOGICAL_GENOME_LABELS,
            reference_label=REFERENCE_LABEL,
        )
    )

    euclidean_pair_trajectory = (
        full_collection.euclidean_pair_trajectory(
            labels=FULL_GENOME_LABELS,
            row_label=REFERENCE_LABEL,
            column_label=COMPARISON_LABEL,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    cosine_pair_trajectory = (
        full_collection.cosine_pair_trajectory(
            labels=FULL_GENOME_LABELS,
            row_label=REFERENCE_LABEL,
            column_label=COMPARISON_LABEL,
            k_values=KMER_SENSITIVITY_LENGTHS,
        )
    )

    euclidean_pair_steps = (
        GenomeCollection
        .pair_trajectory_step_differences(
            euclidean_pair_trajectory
        )
    )

    cosine_pair_steps = (
        GenomeCollection
        .pair_trajectory_step_differences(
            cosine_pair_trajectory
        )
    )

    visualization_paths = save_visualizations(
        euclidean_matrix=euclidean_matrix,
        cosine_matrix=cosine_matrix,
        full_analysis=full_analysis,
        euclidean_pair_trajectory=(
            euclidean_pair_trajectory
        ),
        cosine_pair_trajectory=(
            cosine_pair_trajectory
        ),
    )

    print_report_title(
        "Genome Embeddings Analysis Report",
        (
            f"Reference: {REFERENCE_LABEL} | "
            f"k-mer scales: "
            f"{KMER_SENSITIVITY_LENGTHS}"
        ),
    )

    print_section("1. Dataset")
    print_key_value(
        "Full dataset",
        f"{len(full_collection.genomes)} genomes",
    )
    print_key_value(
        "Biological-only dataset",
        (
            f"{len(biological_collection.genomes)} "
            "genomes"
        ),
    )
    print_key_value(
        "Synthetic control",
        PERIODIC_CONTROL_LABEL,
    )

    print_subsection("Reference sequence")
    print_genome_summary(reference_genome)

    print_section(
        f"2. Single-scale analysis (k={DEFAULT_KMER_LENGTH})"
    )

    print_subsection("Reference descriptor")
    print_descriptor(reference_descriptor)

    print_subsection(
        f"Top {DEFAULT_KMER_LIMIT} reference k-mers"
    )
    print_kmer_frequencies(
        reference_genome,
        k=DEFAULT_KMER_LENGTH,
        limit=DEFAULT_KMER_LIMIT,
    )

    print_subsection(
        f"{REFERENCE_LABEL} vs {COMPARISON_LABEL}"
    )
    print_genome_comparison(comparison)

    print_subsection(
        f"Euclidean ranking from {REFERENCE_LABEL}"
    )
    print_ranking(
        euclidean_matrix.rank_by_label(
            label=REFERENCE_LABEL,
        )
    )

    print_subsection(
        f"Cosine ranking from {REFERENCE_LABEL}"
    )
    print_ranking(
        cosine_matrix.rank_by_label(
            label=REFERENCE_LABEL,
        )
    )

    print_subsection("Euclidean distance matrix")
    print_genome_matrix(euclidean_matrix)

    print_subsection("Cosine similarity matrix")
    print_genome_matrix(cosine_matrix)

    print_section("3. Multiscale geometry")
    print_subsection("Trajectory dimensions")
    print_dataset_dimensions(
        full_analysis=full_analysis,
        biological_analysis=(
            biological_analysis
        ),
    )

    print_subsection(
        "Euclidean full vs biological-only step distances"
    )
    print_dataset_step_comparison(
        full_distances=(
            full_analysis
            .euclidean_step_distances
        ),
        biological_distances=(
            biological_analysis
            .euclidean_step_distances
        ),
    )

    print_subsection(
        "Cosine full vs biological-only step distances"
    )
    print_dataset_step_comparison(
        full_distances=(
            full_analysis
            .cosine_step_distances
        ),
        biological_distances=(
            biological_analysis
            .cosine_step_distances
        ),
    )

    print_subsection(
        "Euclidean squared-deformation partition"
    )
    print_deformation_partition(
        full_analysis
        .euclidean_deformation_partition
    )

    print_subsection(
        "Cosine squared-deformation partition"
    )
    print_deformation_partition(
        full_analysis
        .cosine_deformation_partition
    )

    print_section("4. Ranking stability across k")
    print_subsection("Euclidean ranking trajectory")
    print_ranking_trajectory(
        full_analysis
        .euclidean_ranking_trajectory
    )

    print_subsection("Euclidean stability")
    print_ranking_stability(
        full_analysis
        .euclidean_ranking_stability
    )

    print_subsection("Cosine ranking trajectory")
    print_ranking_trajectory(
        full_analysis
        .cosine_ranking_trajectory
    )

    print_subsection("Cosine stability")
    print_ranking_stability(
        full_analysis
        .cosine_ranking_stability
    )

    print_section("5. Selected pair trajectory")
    print_subsection("Euclidean trajectory")
    print_pair_trajectory(
        euclidean_pair_trajectory
    )

    print_subsection("Euclidean step differences")
    print_pair_step_differences(
        euclidean_pair_steps
    )

    print_subsection("Cosine trajectory")
    print_pair_trajectory(
        cosine_pair_trajectory
    )

    print_subsection("Cosine step differences")
    print_pair_step_differences(
        cosine_pair_steps
    )

    print_section("6. Main deformation drivers")
    print_subsection(
        "Top Euclidean pair contributions"
    )
    print_top_pair_contributions(
        full_analysis
        .euclidean_pair_contributions
    )

    print_subsection(
        "Top cosine pair contributions"
    )
    print_top_pair_contributions(
        full_analysis
        .cosine_pair_contributions
    )

    print_section("7. Generated artifacts")

    for path in visualization_paths:
        print(
            f"  - {path.relative_to(PROJECT_ROOT)}"
        )

    print_subsection("Serialization summary")
    print_key_value(
        "Matrix metric",
        euclidean_matrix.metric,
        indent=2,
    )
    print_key_value(
        "Matrix k-mer length",
        euclidean_matrix.kmer_length,
        indent=2,
    )
    print_key_value(
        "JSON characters",
        len(euclidean_matrix.to_json(indent=2)),
        indent=2,
    )
    print_key_value(
        "CSV rows",
        len(
            euclidean_matrix
            .to_csv()
            .splitlines()
        ),
        indent=2,
    )

    print()
    print("=" * REPORT_WIDTH)
    print("ANALYSIS COMPLETE")
    print("=" * REPORT_WIDTH)


if __name__ == "__main__":
    main()
