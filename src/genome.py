class Genome:
    def __init__(self, organism, chromosome , sequence):
        self.validate_string("Organism", organism)
        self.validate_string("Chromosome", chromosome)
        self.validate_string("Sequence", sequence)

        sequence = sequence.upper()
        valid_nucleotides = {"A", "C", "G", "T"}

        for i, character in enumerate(sequence, start=1):
            if character not in valid_nucleotides:
                raise ValueError(f"Invalid character {character} in sequence at position {i}.")

        self.sequence = sequence
        self.organism = organism
        self.chromosome = chromosome

    @staticmethod
    def validate_string(name, value):
        if type(value) != str:
            raise TypeError(f"{name} must be a string.")
        if len(value) == 0:
            raise ValueError(f"{name} cannot be empty.")

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        gc_bases = {"C", "G"}
        gc_count = 0

        for character in self.sequence:
            if character in gc_bases:
                gc_count += 1
        
        return gc_count / self.length()

    def reverse_complement(self):
        reverse_complement = []

        reversed_sequence = self.sequence[::-1]

        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
        } 
        
        for character in reversed_sequence:
            reverse_complement.append(complement[character])
        
        return "".join(reverse_complement)