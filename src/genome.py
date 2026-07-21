import math
from dataclasses import dataclass


@dataclass(slots=True)
class GenomeDescriptor:
    length: int
    gc_content: float
    at_content: float
    shannon_entropy: float
    gc_skew: float
    purine_content: float
    pyrimidine_content: float
    kmer_length: int
    kmer_diversity: float
    kmer_entropy: float

    def to_dict(self) -> dict[str, int | float]:
        return {
            "length": self.length,
            "gc_content": self.gc_content,
            "at_content": self.at_content,
            "shannon_entropy": self.shannon_entropy,
            "gc_skew": self.gc_skew,
            "purine_content": self.purine_content,
            "pyrimidine_content": self.pyrimidine_content,
            "kmer_length": self.kmer_length,
            "kmer_diversity": self.kmer_diversity,
            "kmer_entropy": self.kmer_entropy,
        }

    def to_vector(self) -> list[float]:
        return [
            float(self.length),
            float(self.gc_content),
            float(self.at_content),
            float(self.shannon_entropy),
            float(self.gc_skew),
            float(self.purine_content),
            float(self.pyrimidine_content),
            float(self.kmer_length),
            float(self.kmer_diversity),
            float(self.kmer_entropy),
        ]


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

    @staticmethod
    def _entropy_from_frequencies(frequencies, total):
        entropy = 0.0

        for frequency in frequencies.values():
            probability = frequency / total
            entropy -= probability * math.log2(probability)

        return entropy

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

        return self._entropy_from_frequencies(frequencies, total)

    def descriptor(self, k=3):
        return GenomeDescriptor(
            length=self.length(),
            gc_content=self.gc_content(),
            at_content=self.at_content(),
            shannon_entropy=self.shannon_entropy(),
            gc_skew=self.gc_skew(),
            purine_content=self.purine_content(),
            pyrimidine_content=self.pyrimidine_content(),
            kmer_length=k,
            kmer_diversity=self.kmer_diversity(k),
            kmer_entropy=self.kmer_entropy(k),
        )
    
    def at_content(self):
        return 1.0 - self.gc_content()
    
    def gc_skew(self):
        frequencies = self.nucleotide_frequencies()

        g_count = frequencies.get("G", 0)
        c_count = frequencies.get("C", 0)
        total_gc = g_count + c_count

        if total_gc == 0:
            return 0.0

        return (g_count - c_count) / total_gc
    
    def purine_content(self):
        frequencies = self.nucleotide_frequencies()

        a_count = frequencies.get("A", 0)
        g_count = frequencies.get("G", 0)

        return (a_count + g_count) / self.length()
    
    def pyrimidine_content(self):
        frequencies = self.nucleotide_frequencies()

        c_count = frequencies.get("C", 0)
        t_count = frequencies.get("T", 0)

        return (c_count + t_count) / self.length()
    
    def kmer_diversity(self, k):
        frequencies = self.kmer_frequencies(k)

        distinct_kmers = len(frequencies)
        total_kmers = len(self.kmers(k))

        return distinct_kmers / total_kmers
    
    def kmer_entropy(self, k):
        frequencies = self.kmer_frequencies(k)
        total_kmers = len(self.kmers(k))

        return self._entropy_from_frequencies(frequencies, total_kmers)
    
    @staticmethod
    def _entropy_from_frequencies(frequencies, total):
        entropy = 0.0

        for frequency in frequencies.values():
            probability = frequency / total
            entropy -= probability * math.log2(probability)

        return entropy