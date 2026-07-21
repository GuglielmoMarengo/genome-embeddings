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
    
    def to_normalized_vector(self) -> list[float]:
        normalized_shannon_entropy = self.shannon_entropy / 2
        normalized_gc_skew = (self.gc_skew + 1) / 2
        normalized_kmer_entropy = (
            self.kmer_entropy / (2 * self.kmer_length)
        )

        return [
            self.gc_content,
            normalized_shannon_entropy,
            normalized_gc_skew,
            self.purine_content,
            self.kmer_diversity,
            normalized_kmer_entropy,
        ]
    
    def euclidean_distance(self, other) -> float:
        if not isinstance(other, GenomeDescriptor):
            raise TypeError("other must be a GenomeDescriptor.")

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        squared_differences = [
            (first_value - second_value) ** 2
            for first_value, second_value in zip(first_vector, second_vector)
        ]

        return math.sqrt(sum(squared_differences))
    
    def cosine_similarity(self, other) -> float:
        if not isinstance(other, GenomeDescriptor):
            raise TypeError("other must be a GenomeDescriptor.")

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        dot_product = sum(
            first_value * second_value
            for first_value, second_value in zip(first_vector, second_vector)
        )

        first_magnitude = math.sqrt(
            sum(value ** 2 for value in first_vector)
        )

        second_magnitude = math.sqrt(
            sum(value ** 2 for value in second_vector)
        )

        if first_magnitude == 0 or second_magnitude == 0:
            return 0.0

        return dot_product / (first_magnitude * second_magnitude)
    
    def feature_differences(self, other) -> dict[str, float]:
        if not isinstance(other, GenomeDescriptor):
            raise TypeError("other must be a GenomeDescriptor.")

        feature_names = [
            "gc_content",
            "normalized_shannon_entropy",
            "normalized_gc_skew",
            "purine_content",
            "kmer_diversity",
            "normalized_kmer_entropy",
        ]

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        return {
            name: abs(first_value - second_value)
            for name, first_value, second_value in zip(
                feature_names,
                first_vector,
                second_vector,
            )
        }
    
    def compare(self, other) -> "GenomeComparison":
        if not isinstance(other, GenomeDescriptor):
            raise TypeError("other must be a GenomeDescriptor.")

        return GenomeComparison(
            euclidean_distance=self.euclidean_distance(other),
            cosine_similarity=self.cosine_similarity(other),
            feature_differences=self.feature_differences(other),
        )

@dataclass(slots=True)
class GenomeComparison:
    euclidean_distance: float
    cosine_similarity: float
    feature_differences: dict[str, float]

    def sorted_feature_differences(self) -> list[tuple[str, float]]:
        return sorted(
            self.feature_differences.items(),
            key=lambda item: item[1],
            reverse=True,
        )

class GenomeCollection:
    def __init__(self, genomes: list["Genome"]) -> None:
        if not isinstance(genomes, list):
            raise TypeError(
                "genomes must be a list of Genome objects."
            )

        if not genomes:
            raise ValueError(
                "Genome collection cannot be empty."
            )

        if not all(isinstance(genome, Genome) for genome in genomes):
            raise TypeError(
                "All items must be Genome objects."
            )

        self.genomes = genomes
    
    def descriptors(self, k: int) -> list[GenomeDescriptor]:
        return [
            genome.descriptor(k=k)
            for genome in self.genomes
        ]
    
    def euclidean_distance_matrix(
        self,
        k: int,
    ) -> list[list[float]]:
        descriptors = self.descriptors(k=k)
        matrix: list[list[float]] = []

        for first_descriptor in descriptors:
            row = [
                first_descriptor.euclidean_distance(second_descriptor)
                for second_descriptor in descriptors
            ]
            matrix.append(row)

        return matrix
    
    def cosine_similarity_matrix(
        self,
        k: int,
    ) -> list[list[float]]:
        descriptors = self.descriptors(k=k)
        matrix: list[list[float]] = []

        for first_descriptor in descriptors:
            row = [
                first_descriptor.cosine_similarity(second_descriptor)
                for second_descriptor in descriptors
            ]
            matrix.append(row)

        return matrix

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
        maximum_distinct_kmers = min(total_kmers, len(self.VALID_NUCLEOTIDES) ** k)

        return distinct_kmers / maximum_distinct_kmers
    
    def kmer_entropy(self, k):
        frequencies = self.kmer_frequencies(k)
        total_kmers = len(self.kmers(k))

        return self._entropy_from_frequencies(frequencies, total_kmers)
