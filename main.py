from pathlib import Path

from src.genome import Genome


PROJECT_ROOT = Path(__file__).resolve().parent
GFP_FASTA_PATH = PROJECT_ROOT / "data" / "gfp.fasta"

DEFAULT_KMER_LENGTH = 3
DEFAULT_KMER_LIMIT = 10


def print_genome_summary(genome):
    print(f"Header: {genome.header}")
    print(f"Sequence (first 100 bp): {genome.sequence[:100]}...")
    print(f"Reverse complement (first 100 bp): {genome.reverse_complement()[:100]}...")


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


def print_descriptor_vector(descriptor):
    print("\nDescriptor Vector:")
    print(descriptor.to_vector())


def print_kmer_frequencies(genome, k=DEFAULT_KMER_LENGTH, limit=DEFAULT_KMER_LIMIT):
    print(f"\nFirst {limit} k-mer frequencies (k={k}):")

    frequencies = genome.kmer_frequencies(k)

    for index, (kmer, frequency) in enumerate(frequencies.items(), start=1):
        print(f"{kmer}: {frequency} times")

        if index == limit:
            break


def main():
    genome = Genome.from_fasta(GFP_FASTA_PATH)
    descriptor = genome.descriptor(k=DEFAULT_KMER_LENGTH)

    print_genome_summary(genome)
    print_descriptor(descriptor)
    print_descriptor_vector(descriptor)
    print_kmer_frequencies(
        genome,
        k=DEFAULT_KMER_LENGTH,
        limit=DEFAULT_KMER_LIMIT,
    )


if __name__ == "__main__":
    main()