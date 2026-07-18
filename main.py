from src.genome import Genome

human = Genome(
    organism="Human",
    chromosome="1",
    sequence="ACGTACGT"
)

print(f"Organism: {human.organism}")
print(f"Chromosome: {human.chromosome}")
print(f"Sequence: {human.sequence}")
print(f"Length: {human.length()}")
print(f"GC Content: {human.gc_content()*100:.2f}%")