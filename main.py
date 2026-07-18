from src.genome import Genome

genome = Genome.from_fasta("data/gfp.fasta")

print(f"Header: {genome.header}")
print(f"Sequence (first 100 bp): {genome.sequence[:100]}...")
print(f"Sequence Length: {genome.length()} bp")
print(f"GC Content: {genome.gc_content()*100:.2f}%")
print(f"Reverse Complement: {genome.reverse_complement()}")