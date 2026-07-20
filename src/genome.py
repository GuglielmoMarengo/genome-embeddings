import math

class Genome:
    VALID_NUCLEOTIDES = {"A", "C", "G", "T"}
    GC_BASES = {"C", "G"}
    COMPLEMENT = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    def __init__(self, sequence, header=None):
        self.validate_string("Sequence", sequence)

        sequence = sequence.upper()
        self._validate_sequence(sequence)

        self.sequence = sequence
        self.header = header

    @classmethod
    def from_fasta(cls, filepath):
        with open(filepath, encoding="utf-8") as file:
            lines = [line.strip() for line in file if line.strip()]

        if not lines:
            raise ValueError("FASTA file is empty.")

        header = lines[0]
        if not header.startswith(">"):
            raise ValueError("FASTA header must start with '>'.")

        sequence = "".join(lines[1:])

        return cls(sequence=sequence, header=header)

    @staticmethod
    def validate_string(name, value):
        if type(value) != str:
            raise TypeError(f"{name} must be a string.")
        if len(value) == 0:
            raise ValueError(f"{name} cannot be empty.")

    @classmethod
    def _validate_sequence(cls, sequence):
        for i, character in enumerate(sequence, start=1):
            if character not in cls.VALID_NUCLEOTIDES:
                raise ValueError(
                    f"Invalid character {character} in sequence at position {i}."
                )

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        gc_count = sum(1 for character in self.sequence if character in self.GC_BASES)
        return gc_count / self.length()

    def reverse_complement(self):
        return "".join(self.COMPLEMENT[character] for character in reversed(self.sequence))

    def kmers(self, k):
        if type(k) != int:
            raise TypeError(f"k must be an integer, got {type(k).__name__}.")
        if k <= 0:
            raise ValueError(f"k must be positive. Got {k}.")
        if k > self.length():
            raise ValueError(f"k cannot exceed the sequence length. Got {k}.")

        sequence_length = self.length()
        return [self.sequence[i:i + k] for i in range(sequence_length - k + 1)]

    def kmer_frequencies(self, k):
        frequencies = {}

        for kmer in self.kmers(k):
            frequencies[kmer] = frequencies.get(kmer, 0) + 1

        return frequencies
    
    def nucleotide_frequencies(self):
        frequencies = {}

        for character in self.sequence:
            frequencies[character] = frequencies.get(character, 0) + 1
        
        return frequencies

    def shannon_entropy(self):
        frequencies = self.nucleotide_frequencies()
        total = self.length()
        entropy = 0

        for frequency in frequencies.values():
            probability = frequency / total
            entropy -= probability * math.log2(probability)

        return entropy