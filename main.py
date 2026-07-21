from pathlib import Path

from src.genome import Genome


PROJECT_ROOT = Path(__file__).resolve().parent
GFP_FASTA_PATH = PROJECT_ROOT / "data" / "gfp.fasta"

COMPARISON_SEQUENCE = "ACGT" * 230 + "AC"
COMPARISON_HEADER = ">Synthetic balanced comparison sequence"

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

    for index, (kmer, frequency) in enumerate(frequencies.items(), start=1):
        print(f"{kmer}: {frequency} times")

        if index == limit:
            break


def print_genome_comparison(first_descriptor, second_descriptor):
    euclidean_distance = first_descriptor.euclidean_distance(second_descriptor)
    cosine_similarity = first_descriptor.cosine_similarity(second_descriptor)
    feature_differences = first_descriptor.feature_differences(second_descriptor)

    print("\nGenome Comparison:")
    print(f"Euclidean distance: {euclidean_distance:.4f}")
    print(f"Cosine similarity: {cosine_similarity:.4f}")

    print("\nFeature Differences:")

    for feature_name, difference in sorted(
        feature_differences.items(),
        key=lambda item: item[1],
        reverse=True,
    ):
        print(f"{feature_name}: {difference:.4f}")


def main():
    genome = Genome.from_fasta(GFP_FASTA_PATH)
    descriptor = genome.descriptor(k=DEFAULT_KMER_LENGTH)

    comparison_genome = Genome(
        sequence=COMPARISON_SEQUENCE,
        header=COMPARISON_HEADER,
    )
    comparison_descriptor = comparison_genome.descriptor(
        k=DEFAULT_KMER_LENGTH
    )

    print_genome_summary(genome)
    print_descriptor(descriptor)
    print_descriptor_vectors(descriptor)
    print_kmer_frequencies(
        genome,
        k=DEFAULT_KMER_LENGTH,
        limit=DEFAULT_KMER_LIMIT,
    )
    print_genome_comparison(
        descriptor,
        comparison_descriptor,
    )


if __name__ == "__main__":
    main()