from pathlib import Path

from src.genome import Genome


PROJECT_ROOT = Path(__file__).resolve().parent
GFP_FASTA_PATH = PROJECT_ROOT / "data" / "gfp.fasta"


def print_kmer_frequencies(genome, k=3, limit=10):
    print(f"First {limit} k-mer frequencies (k={k}):")

    for i, (kmer, frequency) in enumerate(genome.kmer_frequencies(k).items(), start=1):
        print(f"{kmer}: {frequency} times")

        if i == limit:
            break


def main():
    genome = Genome.from_fasta(GFP_FASTA_PATH)
    descriptor = genome.descriptor()

    print(f"Header: {genome.header}")
    print(f"Sequence (first 100 bp): {genome.sequence[:100]}...")
    print(f"Sequence Length: {genome.length()} bp")
    print(f"GC Content: {genome.gc_content() * 100:.2f}%")
    print(f"Shannon entropy: {genome.shannon_entropy():.4f} bits")
    print(f"Reverse Complement (first 100 bp): {genome.reverse_complement()[:100]}...")
    print(f"Genome Descriptor: {descriptor.to_dict()}")
    print(f"Genome Descriptor Vector: {descriptor.to_vector()}")
    print_kmer_frequencies(genome, k=3, limit=10)


if __name__ == "__main__":
    main()