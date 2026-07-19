from src.genome import Genome

genome = Genome.from_fasta("data/gfp.fasta")

print(f"Header: {genome.header}")
print(f"Sequence (first 100 bp): {genome.sequence[:100]}...")
print(f"Sequence Length: {genome.length()} bp")
print(f"GC Content: {genome.gc_content()*100:.2f}%")
print(f"Reverse Complement (first 100 bp): {genome.reverse_complement()[:100]}")
print(f"First 10 K-mers (k=3): {genome.kmers(3)[:10]}...")

frequencies = genome.kmer_frequencies(3)

print("First 10 k-mer frequencies (k=3):")

for i, (kmer, frequency) in enumerate(frequencies.items(), start=1):
    print(f"{kmer}: {frequency} times")

    if i == 10:
        break