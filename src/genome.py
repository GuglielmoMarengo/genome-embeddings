import math
from dataclasses import dataclass
from pathlib import Path


NORMALIZED_FEATURE_NAMES = [
    "gc_content",
    "normalized_shannon_entropy",
    "normalized_gc_skew",
    "purine_content",
    "kmer_diversity",
    "normalized_kmer_entropy",
]

SUPPORTED_MATRIX_METRICS = {
    "euclidean",
    "cosine",
}


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
        normalized_shannon_entropy = (
            self.shannon_entropy / 2
        )

        normalized_gc_skew = (
            self.gc_skew + 1
        ) / 2

        normalized_kmer_entropy = (
            self.kmer_entropy
            / (2 * self.kmer_length)
        )

        return [
            self.gc_content,
            normalized_shannon_entropy,
            normalized_gc_skew,
            self.purine_content,
            self.kmer_diversity,
            normalized_kmer_entropy,
        ]

    def euclidean_distance(
        self,
        other: "GenomeDescriptor",
    ) -> float:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        squared_differences = [
            (first_value - second_value) ** 2
            for first_value, second_value in zip(
                first_vector,
                second_vector,
                strict=True,
            )
        ]

        return math.sqrt(sum(squared_differences))

    def cosine_similarity(
        self,
        other: "GenomeDescriptor",
    ) -> float:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        dot_product = sum(
            first_value * second_value
            for first_value, second_value in zip(
                first_vector,
                second_vector,
                strict=True,
            )
        )

        first_magnitude = math.sqrt(
            sum(
                value ** 2
                for value in first_vector
            )
        )

        second_magnitude = math.sqrt(
            sum(
                value ** 2
                for value in second_vector
            )
        )

        if (
            first_magnitude == 0
            or second_magnitude == 0
        ):
            return 0.0

        return (
            dot_product
            / (first_magnitude * second_magnitude)
        )

    def feature_differences(
        self,
        other: "GenomeDescriptor",
    ) -> dict[str, float]:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        return {
            feature_name: abs(
                first_value - second_value
            )
            for (
                feature_name,
                first_value,
                second_value,
            ) in zip(
                NORMALIZED_FEATURE_NAMES,
                first_vector,
                second_vector,
                strict=True,
            )
        }

    def compare(
        self,
        other: "GenomeDescriptor",
    ) -> "GenomeComparison":
        self._validate_other_descriptor(other)

        return GenomeComparison(
            euclidean_distance=(
                self.euclidean_distance(other)
            ),
            cosine_similarity=(
                self.cosine_similarity(other)
            ),
            feature_differences=(
                self.feature_differences(other)
            ),
        )

    @staticmethod
    def _validate_other_descriptor(
        other: object,
    ) -> None:
        if not isinstance(
            other,
            GenomeDescriptor,
        ):
            raise TypeError(
                "other must be a GenomeDescriptor."
            )


@dataclass(slots=True)
class GenomeComparison:
    euclidean_distance: float
    cosine_similarity: float
    feature_differences: dict[str, float]

    def sorted_feature_differences(
        self,
    ) -> list[tuple[str, float]]:
        return sorted(
            self.feature_differences.items(),
            key=lambda item: item[1],
            reverse=True,
        )


@dataclass(slots=True)
class GenomeMatrix:
    labels: list[str]
    values: list[list[float]]
    metric: str
    kmer_length: int

    def __post_init__(self) -> None:
        if not self.labels:
            raise ValueError(
                "Genome matrix labels cannot be empty."
            )

        matrix_size = len(self.values)

        if any(
            len(row) != matrix_size
            for row in self.values
        ):
            raise ValueError(
                "Genome matrix values must be square."
            )

        if len(self.labels) != matrix_size:
            raise ValueError(
                "Genome matrix labels must match matrix size."
            )

        if self.metric not in SUPPORTED_MATRIX_METRICS:
            raise ValueError(
                "Unsupported genome matrix metric."
            )
    
    def get_value(
        self,
        row_label: str,
        column_label: str,
    ) -> float:
        row_index = self._label_index(row_label)
        column_index = self._label_index(column_label)

        return self.values[row_index][column_index]
    
    def to_rows(
        self,
    ) -> list[dict[str, str | list[float]]]:
        return [
            {
                "label": label,
                "values": row.copy(),
            }
            for label, row in zip(
                self.labels,
                self.values,
                strict=True,
            )
        ]
    
    def to_dict(
        self,
    ) -> dict[
        str,
        list[str]
        | list[list[float]]
        | str
        | int,
    ]:
        return {
            "labels": self.labels.copy(),
            "values": [
                row.copy()
                for row in self.values
            ],
            "metric": self.metric,
            "kmer_length": self.kmer_length,
        }

    def _label_index(
        self,
        label: str,
    ) -> int:
        try:
            return self.labels.index(label)
        except ValueError as error:
            raise ValueError(
                f"Unknown genome matrix label: {label}."
            ) from error


class GenomeCollection:
    def __init__(
        self,
        genomes: list["Genome"],
    ) -> None:
        if not isinstance(genomes, list):
            raise TypeError(
                "genomes must be a list of Genome objects."
            )

        if not genomes:
            raise ValueError(
                "Genome collection cannot be empty."
            )

        if not all(
            isinstance(genome, Genome)
            for genome in genomes
        ):
            raise TypeError(
                "All items must be Genome objects."
            )

        self.genomes = genomes

    def descriptors(
        self,
        k: int,
    ) -> list[GenomeDescriptor]:
        return [
            genome.descriptor(k=k)
            for genome in self.genomes
        ]

    def euclidean_distance_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k=k)

        values = self._build_matrix(
            descriptors=descriptors,
            metric="euclidean",
        )

        return GenomeMatrix(
            labels=labels,
            values=values,
            metric="euclidean",
            kmer_length=k,
        )

    def cosine_similarity_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k=k)

        values = self._build_matrix(
            descriptors=descriptors,
            metric="cosine",
        )

        return GenomeMatrix(
            labels=labels,
            values=values,
            metric="cosine",
            kmer_length=k,
        )

    @staticmethod
    def _build_matrix(
        descriptors: list[GenomeDescriptor],
        metric: str,
    ) -> list[list[float]]:
        matrix: list[list[float]] = []

        for first_descriptor in descriptors:
            row: list[float] = []

            for second_descriptor in descriptors:
                if metric == "euclidean":
                    value = (
                        first_descriptor.euclidean_distance(
                            second_descriptor
                        )
                    )
                elif metric == "cosine":
                    value = (
                        first_descriptor.cosine_similarity(
                            second_descriptor
                        )
                    )
                else:
                    raise ValueError(
                        "Unsupported genome matrix metric."
                    )

                row.append(value)

            matrix.append(row)

        return matrix


class Genome:
    VALID_NUCLEOTIDES = {
        "A",
        "C",
        "G",
        "T",
    }

    GC_BASES = {
        "C",
        "G",
    }

    COMPLEMENT = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    def __init__(
        self,
        sequence: str,
        header: str | None = None,
    ) -> None:
        self.validate_string(
            "Sequence",
            sequence,
        )

        normalized_sequence = sequence.upper()

        self._validate_sequence(
            normalized_sequence
        )

        self.sequence = normalized_sequence
        self.header = header

    @classmethod
    def from_fasta(
        cls,
        filepath: str | Path,
    ) -> "Genome":
        with open(
            filepath,
            encoding="utf-8",
        ) as file:
            lines = [
                line.strip()
                for line in file
                if line.strip()
            ]

        if not lines:
            raise ValueError(
                "FASTA file is empty."
            )

        header = lines[0]

        if not header.startswith(">"):
            raise ValueError(
                "FASTA header must start with '>'."
            )

        sequence = "".join(lines[1:])

        return cls(
            sequence=sequence,
            header=header,
        )

    @staticmethod
    def validate_string(
        name: str,
        value: object,
    ) -> None:
        if type(value) is not str:
            raise TypeError(
                f"{name} must be a string."
            )

        if len(value) == 0:
            raise ValueError(
                f"{name} cannot be empty."
            )

    @classmethod
    def _validate_sequence(
        cls,
        sequence: str,
    ) -> None:
        for position, character in enumerate(
            sequence,
            start=1,
        ):
            if character not in cls.VALID_NUCLEOTIDES:
                raise ValueError(
                    f"Invalid character {character} "
                    f"in sequence at position {position}."
                )

    @staticmethod
    def _entropy_from_frequencies(
        frequencies: dict[str, int],
        total: int,
    ) -> float:
        entropy = 0.0

        for frequency in frequencies.values():
            probability = frequency / total
            entropy -= (
                probability
                * math.log2(probability)
            )

        return entropy

    @staticmethod
    def _validate_k(
        k: int,
        sequence_length: int,
    ) -> None:
        if type(k) is not int:
            raise TypeError(
                "k must be an integer, "
                f"got {type(k).__name__}."
            )

        if k <= 0:
            raise ValueError(
                f"k must be positive. Got {k}."
            )

        if k > sequence_length:
            raise ValueError(
                "k cannot exceed the sequence length. "
                f"Got {k}."
            )

    def length(self) -> int:
        return len(self.sequence)

    def gc_content(self) -> float:
        gc_count = sum(
            1
            for character in self.sequence
            if character in self.GC_BASES
        )

        return gc_count / self.length()

    def at_content(self) -> float:
        return 1.0 - self.gc_content()

    def reverse_complement(self) -> str:
        return "".join(
            self.COMPLEMENT[character]
            for character in reversed(self.sequence)
        )

    def kmers(
        self,
        k: int,
    ) -> list[str]:
        self._validate_k(
            k=k,
            sequence_length=self.length(),
        )

        return [
            self.sequence[index:index + k]
            for index in range(
                self.length() - k + 1
            )
        ]

    def kmer_frequencies(
        self,
        k: int,
    ) -> dict[str, int]:
        frequencies: dict[str, int] = {}

        for kmer in self.kmers(k):
            frequencies[kmer] = (
                frequencies.get(kmer, 0) + 1
            )

        return frequencies

    def nucleotide_frequencies(
        self,
    ) -> dict[str, int]:
        frequencies: dict[str, int] = {}

        for character in self.sequence:
            frequencies[character] = (
                frequencies.get(character, 0) + 1
            )

        return frequencies

    def shannon_entropy(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        return self._entropy_from_frequencies(
            frequencies=frequencies,
            total=self.length(),
        )

    def gc_skew(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        g_count = frequencies.get("G", 0)
        c_count = frequencies.get("C", 0)
        total_gc = g_count + c_count

        if total_gc == 0:
            return 0.0

        return (
            g_count - c_count
        ) / total_gc

    def purine_content(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        a_count = frequencies.get("A", 0)
        g_count = frequencies.get("G", 0)

        return (
            a_count + g_count
        ) / self.length()

    def pyrimidine_content(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        c_count = frequencies.get("C", 0)
        t_count = frequencies.get("T", 0)

        return (
            c_count + t_count
        ) / self.length()

    def kmer_diversity(
        self,
        k: int,
    ) -> float:
        frequencies = self.kmer_frequencies(k)
        total_kmers = self.length() - k + 1

        distinct_kmers = len(frequencies)

        maximum_distinct_kmers = min(
            total_kmers,
            len(self.VALID_NUCLEOTIDES) ** k,
        )

        return (
            distinct_kmers
            / maximum_distinct_kmers
        )

    def kmer_entropy(
        self,
        k: int,
    ) -> float:
        frequencies = self.kmer_frequencies(k)
        total_kmers = self.length() - k + 1

        return self._entropy_from_frequencies(
            frequencies=frequencies,
            total=total_kmers,
        )

    def descriptor(
        self,
        k: int = 3,
    ) -> GenomeDescriptor:
        return GenomeDescriptor(
            length=self.length(),
            gc_content=self.gc_content(),
            at_content=self.at_content(),
            shannon_entropy=self.shannon_entropy(),
            gc_skew=self.gc_skew(),
            purine_content=self.purine_content(),
            pyrimidine_content=(
                self.pyrimidine_content()
            ),
            kmer_length=k,
            kmer_diversity=(
                self.kmer_diversity(k)
            ),
            kmer_entropy=(
                self.kmer_entropy(k)
            ),
        )