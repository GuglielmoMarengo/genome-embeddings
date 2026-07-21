from pathlib import Path

from src.genome import Genome, GenomeCollection


PROJECT_ROOT = Path(__file__).resolve().parent

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
    "Periodic control",
]

DEFAULT_KMER_LENGTH = 3
DEFAULT_KMER_LIMIT = 10


def print_genome_summary(genome):
    print(f"Header: {genome.header}")
    print(f"Sequence (first 100 bp): {genome.sequence[:100]}...")
    print(
        "Reverse complement (first 100 bp): "
        f"{genome.reverse_complement()[:100]}..."
    )


def print_descriptor(descriptor):
    print("\nGenome Descriptor:")
    print(f"Length: {descriptor.length} bp")
    print(f"GC content: {descriptor.gc_content * 100:.2f}%")
    print(f"AT content: {descriptor.at_content * 100:.2f}%")
    print(f"Shannon entropy: {descriptor.shannon_entropy:.4f} bits")
    print(f"GC skew: {descriptor.gc_skew:.4f}")
    print(f"Purine content: {descriptor.purine_content * 100:.2f}%")
    print(f"Pyrimidine content: {descriptor.pyrimidine_content * 100:.2f}%")
    print(f"k-mer length: {descriptor.kmer_length}")
    print(f"k-mer diversity: {descriptor.kmer_diversity:.4f}")
    print(f"k-mer entropy: {descriptor.kmer_entropy:.4f} bits")


def print_descriptor_vectors(descriptor):
    print("\nRaw Descriptor Vector:")
    print(descriptor.to_vector())

    print("\nNormalized Descriptor Vector:")
    print(descriptor.to_normalized_vector())


def print_kmer_frequencies(
    genome,
    k=DEFAULT_KMER_LENGTH,
    limit=DEFAULT_KMER_LIMIT,
):
    print(f"\nFirst {limit} k-mer frequencies (k={k}):")

    frequencies = genome.kmer_frequencies(k)

    for index, (kmer, frequency) in enumerate(
        frequencies.items(),
        start=1,
    ):
        print(f"{kmer}: {frequency} times")

        if index == limit:
            break


def print_genome_comparison(comparison):
    print("\nGenome Comparison:")
    print(f"Euclidean distance: {comparison.euclidean_distance:.4f}")
    print(f"Cosine similarity: {comparison.cosine_similarity:.4f}")

    print("\nFeature Differences:")

    for feature_name, difference in comparison.sorted_feature_differences():
        print(f"{feature_name}: {difference:.4f}")


def print_matrix(
    title: str,
    labels: list[str],
    matrix: list[list[float]],
) -> None:
    label_width = max(len(label) for label in labels)
    value_width = 20

    print(f"\n{title}:")

    header = " " * (label_width + 2)

    for label in labels:
        header += f"{label:>{value_width}}"

    print(header)

    for label, row in zip(labels, matrix, strict=True):
        formatted_values = "".join(
            f"{value:>{value_width}.4f}"
            for value in row
        )

        print(
            f"{label:<{label_width}}  "
            f"{formatted_values}"
        )


def main():
    genomes = [
        Genome.from_fasta(AEQUOREA_GFP_PATH),
        Genome.from_fasta(ACROPORA_GFP_PATH),
        Genome.from_fasta(DISCOSOMA_FP583_PATH),
        Genome.from_fasta(PERIODIC_CONTROL_PATH),
    ]

    collection = GenomeCollection(genomes)

    reference_genome = genomes[0]
    reference_descriptor = reference_genome.descriptor(
        k=DEFAULT_KMER_LENGTH
    )

    comparison_descriptor = genomes[1].descriptor(
        k=DEFAULT_KMER_LENGTH
    )

    comparison = reference_descriptor.compare(
        comparison_descriptor
    )

    print_genome_summary(reference_genome)
    print_descriptor(reference_descriptor)
    print_descriptor_vectors(reference_descriptor)

    print_kmer_frequencies(
        reference_genome,
        k=DEFAULT_KMER_LENGTH,
        limit=DEFAULT_KMER_LIMIT,
    )

    print_genome_comparison(comparison)

    euclidean_matrix = collection.euclidean_distance_matrix(
        k=DEFAULT_KMER_LENGTH
    )

    cosine_matrix = collection.cosine_similarity_matrix(
        k=DEFAULT_KMER_LENGTH
    )

    print_matrix(
        title="Euclidean Distance Matrix",
        labels=GENOME_LABELS,
        matrix=euclidean_matrix,
    )

    print_matrix(
        title="Cosine Similarity Matrix",
        labels=GENOME_LABELS,
        matrix=cosine_matrix,
    )


if __name__ == "__main__":
    main()