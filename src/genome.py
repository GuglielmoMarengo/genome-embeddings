class Genome:
    def __init__(self, organism, chromosome ,sequence):

        if type(sequence) != str:
            raise TypeError("Sequence must be a string.")
        if len(sequence) == 0:
            raise ValueError("Sequence cannot be empty.")

        sequence = sequence.upper()
        valid_nucleotides = {"A", "C", "G", "T"}

        for i, character in enumerate(sequence, start=1):
            if character not in valid_nucleotides:
                raise ValueError(f"Invalid character {character} in sequence at position {i}.")

        if type(organism) != str:
            raise TypeError("Organism must be a string.")
        if len(organism) == 0:
            raise ValueError("Organism cannot be empty.")

        if type(chromosome) != str:
            raise TypeError("Chromosome must be a string.")
        if len(chromosome) == 0:
            raise ValueError("Chromosome cannot be empty.")

        self.sequence = sequence
        self.organism = organism
        self.chromosome = chromosome

    def length(self):
        return len(self.sequence)