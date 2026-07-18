from src.genome import Genome

genome = Genome.from_fasta("data/example.fasta")

print(f"Header: {genome.header}")
print(f"Sequence: {genome.sequence}")
print(f"Length: {genome.length()}")
print(f"GC Content: {genome.gc_content()*100:.2f}%")
print(f"Reverse Complement: {genome.reverse_complement()}")