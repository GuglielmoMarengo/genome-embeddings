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
    plot_trajectory_distributions,
    save_figure,
)


PROJECT_ROOT = (
    Path(__file__).resolve().parent
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs"
)

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

GENOME_LABELS = [
    "Aequorea GFP",
    "Acropora GFP",
    "Discosoma FP583",
    "S. aureus catA",
    "S. cerevisiae TPI1",
    "Periodic control",
]

DEFAULT_KMER_LENGTH = 3
DEFAULT_KMER_LIMIT = 10

KMER_SENSITIVITY_LENGTHS = [
    1,
    2,
    3,
    4,
]

REFERENCE_LABEL = "Aequorea GFP"
COMPARISON_LABEL = "Acropora GFP"


def print_genome_summary(
    genome: Genome,
) -> None:
    print(f"Header: {genome.header}")

    print(
        "Sequence (first 100 bp): "
        f"{genome.sequence[:100]}..."
    )

    print(
        "Reverse complement "
        "(first 100 bp): "
        f"{genome.reverse_complement()[:100]}..."
    )


def print_descriptor(
    descriptor: GenomeDescriptor,
) -> None:
    print("\nGenome Descriptor:")
    print(
        f"Length: {descriptor.length} bp"
    )

    print(
        "GC content: "
        f"{descriptor.gc_content * 100:.2f}%"
    )

    print(
        "AT content: "
        f"{descriptor.at_content * 100:.2f}%"
    )

    print(
        "Shannon entropy: "
        f"{descriptor.shannon_entropy:.4f} bits"
    )

    print(
        f"GC skew: "
        f"{descriptor.gc_skew:.4f}"
    )

    print(
        "Purine content: "
        f"{descriptor.purine_content * 100:.2f}%"
    )

    print(
        "Pyrimidine content: "
        f"{descriptor.pyrimidine_content * 100:.2f}%"
    )

    print(
        f"k-mer length: "
        f"{descriptor.kmer_length}"
    )

    print(
        "k-mer diversity: "
        f"{descriptor.kmer_diversity:.4f}"
    )

    print(
        "k-mer entropy: "
        f"{descriptor.kmer_entropy:.4f} bits"
    )


def print_descriptor_vectors(
    descriptor: GenomeDescriptor,
) -> None:
    print("\nRaw Descriptor Vector:")
    print(descriptor.to_vector())

    print(
        "\nNormalized Descriptor Vector:"
    )

    print(
        descriptor.to_normalized_vector()
    )


def print_kmer_frequencies(
    genome: Genome,
    k: int = DEFAULT_KMER_LENGTH,
    limit: int = DEFAULT_KMER_LIMIT,
) -> None:
    print(
        f"\nFirst {limit} "
        f"k-mer frequencies (k={k}):"
    )

    frequencies = (
        genome.kmer_frequencies(k)
    )

    for index, (
        kmer,
        frequency,
    ) in enumerate(
        frequencies.items(),
        start=1,
    ):
        print(
            f"{kmer}: "
            f"{frequency} times"
        )

        if index == limit:
            break


def print_genome_comparison(
    comparison: GenomeComparison,
) -> None:
    print("\nGenome Comparison:")

    print(
        "Euclidean distance: "
        f"{comparison.euclidean_distance:.4f}"
    )

    print(
        "Cosine similarity: "
        f"{comparison.cosine_similarity:.4f}"
    )

    print("\nFeature Differences:")

    for (
        feature_name,
        difference,
    ) in (
        comparison
        .sorted_feature_differences()
    ):
        print(
            f"{feature_name}: "
            f"{difference:.4f}"
        )


def print_genome_matrix(
    matrix: GenomeMatrix,
) -> None:
    label_width = max(
        len(label)
        for label in matrix.labels
    )

    value_width = max(
        20,
        max(
            len(label) + 2
            for label in matrix.labels
        ),
    )

    metric_title = {
        "euclidean": (
            "Euclidean Distance Matrix"
        ),
        "cosine": (
            "Cosine Similarity Matrix"
        ),
    }[matrix.metric]

    print(
        f"\n{metric_title} "
        f"(k={matrix.kmer_length}):"
    )

    header = (
        " " * (label_width + 2)
    )

    for label in matrix.labels:
        header += (
            f"{label:>{value_width}}"
        )

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


def print_matrix_trajectory(
    trajectory: dict[int, list[float]],
    metric_name: str,
    preview_limit: int = 3,
) -> None:
    print(
        f"\n{metric_name} "
        "Matrix-Geometry Trajectory:"
    )

    for k, vector in trajectory.items():
        preview = (
            vector[:preview_limit]
        )

        print(
            f"k={k}: "
            f"dimensions={len(vector)}, "
            f"first values={preview}"
        )


def print_pair_trajectory(
    trajectory: dict[int, float],
    row_label: str,
    column_label: str,
    metric_name: str,
) -> None:
    print(
        f"\n{row_label} -> "
        f"{column_label} "
        f"{metric_name} Pair Trajectory:"
    )

    for k, value in trajectory.items():
        print(
            f"k={k}: {value:.4f}"
        )


def load_genomes() -> list[Genome]:
    return [
        Genome.from_fasta(
            AEQUOREA_GFP_PATH
        ),
        Genome.from_fasta(
            ACROPORA_GFP_PATH
        ),
        Genome.from_fasta(
            DISCOSOMA_FP583_PATH
        ),
        Genome.from_fasta(
            STAPHYLOCOCCUS_AUREUS_CATA_PATH
        ),
        Genome.from_fasta(
            SACCHAROMYCES_CEREVISIAE_TPI1_PATH
        ),
        Genome.from_fasta(
            PERIODIC_CONTROL_PATH
        ),
    ]


def save_visualizations(
    euclidean_matrix: GenomeMatrix,
    cosine_matrix: GenomeMatrix,
    euclidean_trajectory: dict[
        int,
        list[float],
    ],
    cosine_trajectory: dict[
        int,
        list[float],
    ],
    euclidean_pair_trajectory: dict[
        int,
        float,
    ],
    cosine_pair_trajectory: dict[
        int,
        float,
    ],
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
            (
                "aequorea_acropora_"
                "euclidean_trajectory.png"
            ),
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
            (
                "aequorea_acropora_"
                "cosine_trajectory.png"
            ),
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
            (
                "euclidean_multi_k_"
                "distribution.png"
            ),
            plot_trajectory_distributions(
                trajectory=(
                    euclidean_trajectory
                ),
                metric="euclidean",
            ),
        ),
        (
            (
                "cosine_multi_k_"
                "distribution.png"
            ),
            plot_trajectory_distributions(
                trajectory=(
                    cosine_trajectory
                ),
                metric="cosine",
            ),
        ),
    ]

    output_paths = []

    for filename, figure in figures:
        output_path = save_figure(
            figure=figure,
            output_path=(
                OUTPUT_DIR
                / filename
            ),
            close=True,
        )

        output_paths.append(
            output_path
        )

    return output_paths


def main() -> None:
    genomes = load_genomes()

    collection = GenomeCollection(
        genomes
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

    comparison = (
        reference_descriptor.compare(
            comparison_descriptor
        )
    )

    print_genome_summary(
        reference_genome
    )

    print_descriptor(
        reference_descriptor
    )

    print_descriptor_vectors(
        reference_descriptor
    )

    print_kmer_frequencies(
        reference_genome,
        k=DEFAULT_KMER_LENGTH,
        limit=DEFAULT_KMER_LIMIT,
    )

    print_genome_comparison(
        comparison
    )

    euclidean_matrix = (
        collection
        .euclidean_distance_matrix(
            labels=GENOME_LABELS,
            k=DEFAULT_KMER_LENGTH,
        )
    )

    cosine_matrix = (
        collection
        .cosine_similarity_matrix(
            labels=GENOME_LABELS,
            k=DEFAULT_KMER_LENGTH,
        )
    )

    print_genome_matrix(
        euclidean_matrix
    )

    print_genome_matrix(
        cosine_matrix
    )

    euclidean_ranking = (
        euclidean_matrix.rank_by_label(
            label=REFERENCE_LABEL,
        )
    )

    cosine_ranking = (
        cosine_matrix.rank_by_label(
            label=REFERENCE_LABEL,
        )
    )

    print(
        "\nEuclidean Ranking "
        f"from {REFERENCE_LABEL}:"
    )

    for label, value in (
        euclidean_ranking
    ):
        print(
            f"{label}: {value:.4f}"
        )

    print(
        "\nCosine Ranking "
        f"from {REFERENCE_LABEL}:"
    )

    for label, value in (
        cosine_ranking
    ):
        print(
            f"{label}: {value:.4f}"
        )

    euclidean_trajectory = (
        collection
        .euclidean_matrix_trajectory(
            labels=GENOME_LABELS,
            k_values=(
                KMER_SENSITIVITY_LENGTHS
            ),
        )
    )

    cosine_trajectory = (
        collection
        .cosine_matrix_trajectory(
            labels=GENOME_LABELS,
            k_values=(
                KMER_SENSITIVITY_LENGTHS
            ),
        )
    )

    euclidean_pair_trajectory = (
        collection
        .euclidean_pair_trajectory(
            labels=GENOME_LABELS,
            row_label=REFERENCE_LABEL,
            column_label=COMPARISON_LABEL,
            k_values=(
                KMER_SENSITIVITY_LENGTHS
            ),
        )
    )

    cosine_pair_trajectory = (
        collection
        .cosine_pair_trajectory(
            labels=GENOME_LABELS,
            row_label=REFERENCE_LABEL,
            column_label=COMPARISON_LABEL,
            k_values=(
                KMER_SENSITIVITY_LENGTHS
            ),
        )
    )

    print_matrix_trajectory(
        trajectory=(
            euclidean_trajectory
        ),
        metric_name="Euclidean",
    )

    print_matrix_trajectory(
        trajectory=cosine_trajectory,
        metric_name="Cosine",
    )

    print_pair_trajectory(
        trajectory=(
            euclidean_pair_trajectory
        ),
        row_label=REFERENCE_LABEL,
        column_label=COMPARISON_LABEL,
        metric_name="Euclidean",
    )

    print_pair_trajectory(
        trajectory=(
            cosine_pair_trajectory
        ),
        row_label=REFERENCE_LABEL,
        column_label=COMPARISON_LABEL,
        metric_name="Cosine",
    )

    visualization_paths = (
        save_visualizations(
            euclidean_matrix=(
                euclidean_matrix
            ),
            cosine_matrix=(
                cosine_matrix
            ),
            euclidean_trajectory=(
                euclidean_trajectory
            ),
            cosine_trajectory=(
                cosine_trajectory
            ),
            euclidean_pair_trajectory=(
                euclidean_pair_trajectory
            ),
            cosine_pair_trajectory=(
                cosine_pair_trajectory
            ),
        )
    )

    print(
        "\nGenerated Visualizations:"
    )

    for path in visualization_paths:
        print(
            path.relative_to(
                PROJECT_ROOT
            )
        )

    euclidean_json = (
        euclidean_matrix.to_json(
            indent=2,
        )
    )

    euclidean_csv = (
        euclidean_matrix.to_csv()
    )

    print(
        "\nEuclidean Matrix "
        "JSON Preview:"
    )

    print(
        euclidean_json[:300]
        + "..."
    )

    print(
        "\nEuclidean Matrix "
        "CSV Preview:"
    )

    print(
        "\n".join(
            euclidean_csv
            .splitlines()[:3]
        )
    )

    selected_distance = (
        euclidean_matrix.get_value(
            row_label=REFERENCE_LABEL,
            column_label=COMPARISON_LABEL,
        )
    )

    print(
        "\nSelected Euclidean Distance:"
    )

    print(
        f"{REFERENCE_LABEL} -> "
        f"{COMPARISON_LABEL}: "
        f"{selected_distance:.4f}"
    )

    matrix_rows = (
        euclidean_matrix.to_rows()
    )

    matrix_dict = (
        euclidean_matrix.to_dict()
    )

    print(
        "\nFirst Euclidean Matrix Row:"
    )

    print(
        matrix_rows[0]
    )

    print(
        "\nEuclidean Matrix Metadata:"
    )

    print(
        f"Metric: "
        f"{matrix_dict['metric']}"
    )

    print(
        "k-mer length: "
        f"{matrix_dict['kmer_length']}"
    )

    print(
        f"Labels: "
        f"{matrix_dict['labels']}"
    )


if __name__ == "__main__":
    main()