"""Versioned, interpretable descriptors for finite-sample sequence analysis.

The legacy :class:`src.genome.GenomeDescriptor` remains unchanged.  This
module introduces a second descriptor family with explicit finite-sample,
dependency, and sparsity diagnostics and a first interpretable multiscale
embedding.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from itertools import product

from src.genome import Genome, GenomeCollection, GenomeMatrix


NUCLEOTIDE_ORDER = ("A", "C", "G", "T")
DINUCLEOTIDE_ORDER = tuple(
    "".join(symbols)
    for symbols in product(NUCLEOTIDE_ORDER, repeat=2)
)

GLOBAL_V2_FEATURE_NAMES = (
    "gc_content",
    "normalized_shannon_entropy",
    "normalized_gc_skew",
    "purine_content",
    "normalized_conditional_nucleotide_entropy",
)

DINUCLEOTIDE_V2_FEATURE_NAMES = tuple(
    f"dinucleotide_{dinucleotide}_bounded_odds_ratio"
    for dinucleotide in DINUCLEOTIDE_ORDER
)

KMER_V2_FEATURE_NAMES = (
    "finite_sample_normalized_kmer_entropy",
    "effective_kmer_fraction",
    "theoretical_space_coverage",
    "observable_space_coverage",
    "singleton_fraction",
)

DESCRIPTOR_V2_FEATURE_NAMES = (
    *GLOBAL_V2_FEATURE_NAMES,
    *DINUCLEOTIDE_V2_FEATURE_NAMES,
    *KMER_V2_FEATURE_NAMES,
)


@dataclass(frozen=True, slots=True)
class KmerScaleDescriptorV2:
    """Finite-sample and sparsity diagnostics for one k-mer scale."""

    kmer_length: int
    window_count: int
    possible_kmer_count: int
    observable_kmer_count: int
    distinct_kmer_count: int
    kmer_entropy: float
    theoretical_normalized_kmer_entropy: float
    finite_sample_normalized_kmer_entropy: float
    effective_kmer_count: float
    theoretical_space_coverage: float
    observable_space_coverage: float
    singleton_fraction: float
    repeated_window_fraction: float

    def to_dict(self) -> dict[str, int | float]:
        return {
            "kmer_length": self.kmer_length,
            "window_count": self.window_count,
            "possible_kmer_count": self.possible_kmer_count,
            "observable_kmer_count": self.observable_kmer_count,
            "distinct_kmer_count": self.distinct_kmer_count,
            "kmer_entropy": self.kmer_entropy,
            "theoretical_normalized_kmer_entropy": (
                self.theoretical_normalized_kmer_entropy
            ),
            "finite_sample_normalized_kmer_entropy": (
                self.finite_sample_normalized_kmer_entropy
            ),
            "effective_kmer_count": self.effective_kmer_count,
            "theoretical_space_coverage": self.theoretical_space_coverage,
            "observable_space_coverage": self.observable_space_coverage,
            "singleton_fraction": self.singleton_fraction,
            "repeated_window_fraction": self.repeated_window_fraction,
        }

    @property
    def effective_kmer_fraction(self) -> float:
        if self.observable_kmer_count == 0:
            return 0.0
        return self.effective_kmer_count / self.observable_kmer_count

    def to_normalized_vector(self) -> list[float]:
        return [
            self.finite_sample_normalized_kmer_entropy,
            self.effective_kmer_fraction,
            self.theoretical_space_coverage,
            self.observable_space_coverage,
            self.singleton_fraction,
        ]


@dataclass(frozen=True, slots=True)
class GenomeDescriptorV2:
    """Interpretable descriptor with global, dependency, and k-mer blocks."""

    length: int
    gc_content: float
    shannon_entropy: float
    gc_skew: float
    purine_content: float
    conditional_nucleotide_entropy: float
    dinucleotide_odds_ratios: dict[str, float]
    kmer: KmerScaleDescriptorV2

    @classmethod
    def from_genome(
        cls,
        genome: Genome,
        k: int,
    ) -> "GenomeDescriptorV2":
        if not isinstance(genome, Genome):
            raise TypeError("genome must be a Genome object.")

        return cls(
            length=genome.length(),
            gc_content=genome.gc_content(),
            shannon_entropy=genome.shannon_entropy(),
            gc_skew=genome.gc_skew(),
            purine_content=genome.purine_content(),
            conditional_nucleotide_entropy=(
                conditional_nucleotide_entropy(genome)
            ),
            dinucleotide_odds_ratios=(
                dinucleotide_odds_ratios(genome)
            ),
            kmer=kmer_scale_descriptor(genome, k),
        )

    @property
    def normalized_shannon_entropy(self) -> float:
        return self.shannon_entropy / 2.0

    @property
    def normalized_gc_skew(self) -> float:
        return (self.gc_skew + 1.0) / 2.0

    @property
    def normalized_conditional_nucleotide_entropy(self) -> float:
        return self.conditional_nucleotide_entropy / 2.0

    def bounded_dinucleotide_odds_ratios(self) -> dict[str, float]:
        return {
            dinucleotide: _bounded_nonnegative_value(
                self.dinucleotide_odds_ratios[dinucleotide]
            )
            for dinucleotide in DINUCLEOTIDE_ORDER
        }

    def to_dict(self) -> dict[str, object]:
        return {
            "version": 2,
            "length": self.length,
            "gc_content": self.gc_content,
            "shannon_entropy": self.shannon_entropy,
            "normalized_shannon_entropy": self.normalized_shannon_entropy,
            "gc_skew": self.gc_skew,
            "normalized_gc_skew": self.normalized_gc_skew,
            "purine_content": self.purine_content,
            "conditional_nucleotide_entropy": (
                self.conditional_nucleotide_entropy
            ),
            "normalized_conditional_nucleotide_entropy": (
                self.normalized_conditional_nucleotide_entropy
            ),
            "dinucleotide_odds_ratios": {
                key: self.dinucleotide_odds_ratios[key]
                for key in DINUCLEOTIDE_ORDER
            },
            "kmer": self.kmer.to_dict(),
        }

    def to_feature_dict(self) -> dict[str, float]:
        bounded_ratios = self.bounded_dinucleotide_odds_ratios()
        feature_values = [
            self.gc_content,
            self.normalized_shannon_entropy,
            self.normalized_gc_skew,
            self.purine_content,
            self.normalized_conditional_nucleotide_entropy,
            *(
                bounded_ratios[dinucleotide]
                for dinucleotide in DINUCLEOTIDE_ORDER
            ),
            *self.kmer.to_normalized_vector(),
        ]

        return dict(
            zip(
                DESCRIPTOR_V2_FEATURE_NAMES,
                feature_values,
                strict=True,
            )
        )

    def to_normalized_vector(self) -> list[float]:
        return list(self.to_feature_dict().values())

    def euclidean_distance(self, other: "GenomeDescriptorV2") -> float:
        self._validate_other(other)
        return math.dist(
            self.to_normalized_vector(),
            other.to_normalized_vector(),
        )

    def cosine_similarity(self, other: "GenomeDescriptorV2") -> float:
        self._validate_other(other)
        first = self.to_normalized_vector()
        second = other.to_normalized_vector()
        dot_product = sum(
            first_value * second_value
            for first_value, second_value in zip(
                first,
                second,
                strict=True,
            )
        )
        first_magnitude = math.sqrt(sum(value ** 2 for value in first))
        second_magnitude = math.sqrt(sum(value ** 2 for value in second))
        if first_magnitude == 0.0 or second_magnitude == 0.0:
            return 0.0
        return dot_product / (first_magnitude * second_magnitude)

    def feature_differences(
        self,
        other: "GenomeDescriptorV2",
    ) -> dict[str, float]:
        self._validate_other(other)
        first = self.to_feature_dict()
        second = other.to_feature_dict()
        return {
            feature_name: abs(first[feature_name] - second[feature_name])
            for feature_name in DESCRIPTOR_V2_FEATURE_NAMES
        }

    @staticmethod
    def _validate_other(other: object) -> None:
        if not isinstance(other, GenomeDescriptorV2):
            raise TypeError("other must be a GenomeDescriptorV2.")


@dataclass(frozen=True, slots=True)
class GenomeEmbeddingV2:
    """Global features once plus scale-specific diagnostics for each k."""

    k_values: tuple[int, ...]
    feature_names: tuple[str, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.k_values:
            raise ValueError("Embedding k-mer lengths cannot be empty.")
        if len(self.k_values) != len(set(self.k_values)):
            raise ValueError("Embedding k-mer lengths must be unique.")
        if len(self.feature_names) != len(self.values):
            raise ValueError("Embedding feature names must match values.")

    def to_dict(self) -> dict[str, object]:
        return {
            "version": 2,
            "k_values": list(self.k_values),
            "features": dict(
                zip(self.feature_names, self.values, strict=True)
            ),
        }

    def to_json(self, indent: int | None = None) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def euclidean_distance(self, other: "GenomeEmbeddingV2") -> float:
        self._validate_compatible(other)
        return math.dist(self.values, other.values)

    def cosine_similarity(self, other: "GenomeEmbeddingV2") -> float:
        self._validate_compatible(other)
        dot_product = sum(
            first * second
            for first, second in zip(
                self.values,
                other.values,
                strict=True,
            )
        )
        first_magnitude = math.sqrt(sum(value ** 2 for value in self.values))
        second_magnitude = math.sqrt(sum(value ** 2 for value in other.values))
        if first_magnitude == 0.0 or second_magnitude == 0.0:
            return 0.0
        return dot_product / (first_magnitude * second_magnitude)

    def _validate_compatible(self, other: object) -> None:
        if not isinstance(other, GenomeEmbeddingV2):
            raise TypeError("other must be a GenomeEmbeddingV2.")
        if self.feature_names != other.feature_names:
            raise ValueError("Embedding features must match.")


class DescriptorV2Collection:
    """Collection-level matrices for V2 descriptors and embeddings."""

    def __init__(self, genomes: list[Genome]) -> None:
        self.collection = GenomeCollection(genomes)

    @property
    def genomes(self) -> list[Genome]:
        return self.collection.genomes

    def descriptors(self, k: int) -> list[GenomeDescriptorV2]:
        return [
            GenomeDescriptorV2.from_genome(genome, k)
            for genome in self.genomes
        ]

    def embeddings(self, k_values: list[int]) -> list[GenomeEmbeddingV2]:
        _validate_k_values(k_values)
        return [
            multiscale_embedding(genome, k_values)
            for genome in self.genomes
        ]

    def euclidean_distance_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k)
        return _comparison_matrix(
            labels=labels,
            items=descriptors,
            metric="euclidean_v2",
            kmer_length=k,
            comparison=lambda first, second: first.euclidean_distance(second),
        )

    def cosine_similarity_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k)
        return _comparison_matrix(
            labels=labels,
            items=descriptors,
            metric="cosine_v2",
            kmer_length=k,
            comparison=lambda first, second: first.cosine_similarity(second),
        )

    def multiscale_embedding_distance_matrix(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> GenomeMatrix:
        embeddings = self.embeddings(k_values)
        return _comparison_matrix(
            labels=labels,
            items=embeddings,
            metric="embedding_v2_euclidean",
            kmer_length=max(k_values),
            comparison=lambda first, second: first.euclidean_distance(second),
        )

    def multiscale_embedding_cosine_matrix(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> GenomeMatrix:
        embeddings = self.embeddings(k_values)
        return _comparison_matrix(
            labels=labels,
            items=embeddings,
            metric="embedding_v2_cosine",
            kmer_length=max(k_values),
            comparison=lambda first, second: first.cosine_similarity(second),
        )

    def euclidean_distance_matrices(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, GenomeMatrix]:
        _validate_k_values(k_values)
        return {
            k: self.euclidean_distance_matrix(labels, k)
            for k in k_values
        }

    def cosine_similarity_matrices(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, GenomeMatrix]:
        _validate_k_values(k_values)
        return {
            k: self.cosine_similarity_matrix(labels, k)
            for k in k_values
        }


def conditional_nucleotide_entropy(genome: Genome) -> float:
    """Return H(next nucleotide | current nucleotide) in bits."""

    if not isinstance(genome, Genome):
        raise TypeError("genome must be a Genome object.")
    if genome.length() < 2:
        return 0.0

    transition_counts = {
        current: {next_base: 0 for next_base in NUCLEOTIDE_ORDER}
        for current in NUCLEOTIDE_ORDER
    }
    current_counts = {base: 0 for base in NUCLEOTIDE_ORDER}

    for current, next_base in zip(
        genome.sequence[:-1],
        genome.sequence[1:],
        strict=True,
    ):
        transition_counts[current][next_base] += 1
        current_counts[current] += 1

    total_transitions = genome.length() - 1
    entropy = 0.0

    for current in NUCLEOTIDE_ORDER:
        current_count = current_counts[current]
        if current_count == 0:
            continue

        conditional_entropy = 0.0
        for next_base in NUCLEOTIDE_ORDER:
            count = transition_counts[current][next_base]
            if count == 0:
                continue
            probability = count / current_count
            conditional_entropy -= probability * math.log2(probability)

        entropy += (current_count / total_transitions) * conditional_entropy

    return entropy


def dinucleotide_odds_ratios(genome: Genome) -> dict[str, float]:
    """Return observed/expected dinucleotide abundance ratios."""

    if not isinstance(genome, Genome):
        raise TypeError("genome must be a Genome object.")

    nucleotide_counts = genome.nucleotide_frequencies()
    nucleotide_probabilities = {
        base: nucleotide_counts.get(base, 0) / genome.length()
        for base in NUCLEOTIDE_ORDER
    }

    if genome.length() < 2:
        return {dinucleotide: 0.0 for dinucleotide in DINUCLEOTIDE_ORDER}

    pair_counts = {dinucleotide: 0 for dinucleotide in DINUCLEOTIDE_ORDER}
    for index in range(genome.length() - 1):
        pair_counts[genome.sequence[index:index + 2]] += 1

    total_pairs = genome.length() - 1
    ratios: dict[str, float] = {}

    for dinucleotide in DINUCLEOTIDE_ORDER:
        first, second = dinucleotide
        observed = pair_counts[dinucleotide] / total_pairs
        expected = (
            nucleotide_probabilities[first]
            * nucleotide_probabilities[second]
        )
        ratios[dinucleotide] = observed / expected if expected > 0.0 else 0.0

    return ratios


def kmer_scale_descriptor(
    genome: Genome,
    k: int,
) -> KmerScaleDescriptorV2:
    """Build finite-sample diagnostics for a genome and k-mer scale."""

    if not isinstance(genome, Genome):
        raise TypeError("genome must be a Genome object.")

    frequencies = genome.kmer_frequencies(k)
    window_count = genome.length() - k + 1
    possible_kmer_count = len(NUCLEOTIDE_ORDER) ** k
    observable_kmer_count = min(possible_kmer_count, window_count)
    distinct_kmer_count = len(frequencies)
    entropy = genome.kmer_entropy(k)

    theoretical_maximum_entropy = math.log2(possible_kmer_count)
    observable_maximum_entropy = (
        math.log2(observable_kmer_count)
        if observable_kmer_count > 1
        else 0.0
    )

    theoretical_normalized_entropy = (
        entropy / theoretical_maximum_entropy
        if theoretical_maximum_entropy > 0.0
        else 0.0
    )
    finite_sample_normalized_entropy = (
        entropy / observable_maximum_entropy
        if observable_maximum_entropy > 0.0
        else 0.0
    )

    singleton_windows = sum(
        1
        for count in frequencies.values()
        if count == 1
    )
    singleton_fraction = singleton_windows / window_count

    return KmerScaleDescriptorV2(
        kmer_length=k,
        window_count=window_count,
        possible_kmer_count=possible_kmer_count,
        observable_kmer_count=observable_kmer_count,
        distinct_kmer_count=distinct_kmer_count,
        kmer_entropy=entropy,
        theoretical_normalized_kmer_entropy=(
            theoretical_normalized_entropy
        ),
        finite_sample_normalized_kmer_entropy=(
            finite_sample_normalized_entropy
        ),
        effective_kmer_count=2.0 ** entropy,
        theoretical_space_coverage=(
            distinct_kmer_count / possible_kmer_count
        ),
        observable_space_coverage=(
            distinct_kmer_count / observable_kmer_count
        ),
        singleton_fraction=singleton_fraction,
        repeated_window_fraction=1.0 - singleton_fraction,
    )


def multiscale_embedding(
    genome: Genome,
    k_values: list[int],
) -> GenomeEmbeddingV2:
    """Create an interpretable embedding without duplicating global features."""

    _validate_k_values(k_values)
    descriptors = [
        GenomeDescriptorV2.from_genome(genome, k)
        for k in k_values
    ]
    base_descriptor = descriptors[0]
    bounded_ratios = base_descriptor.bounded_dinucleotide_odds_ratios()

    feature_names = [
        *GLOBAL_V2_FEATURE_NAMES,
        *DINUCLEOTIDE_V2_FEATURE_NAMES,
    ]
    values = [
        base_descriptor.gc_content,
        base_descriptor.normalized_shannon_entropy,
        base_descriptor.normalized_gc_skew,
        base_descriptor.purine_content,
        base_descriptor.normalized_conditional_nucleotide_entropy,
        *(
            bounded_ratios[dinucleotide]
            for dinucleotide in DINUCLEOTIDE_ORDER
        ),
    ]

    for descriptor in descriptors:
        for feature_name, value in zip(
            KMER_V2_FEATURE_NAMES,
            descriptor.kmer.to_normalized_vector(),
            strict=True,
        ):
            feature_names.append(f"k{descriptor.kmer.kmer_length}_{feature_name}")
            values.append(value)

    return GenomeEmbeddingV2(
        k_values=tuple(k_values),
        feature_names=tuple(feature_names),
        values=tuple(values),
    )


def _comparison_matrix(
    *,
    labels: list[str],
    items: list[object],
    metric: str,
    kmer_length: int,
    comparison,
) -> GenomeMatrix:
    if len(labels) != len(items):
        raise ValueError("Genome labels must match collection size.")

    values = [
        [comparison(first, second) for second in items]
        for first in items
    ]
    return GenomeMatrix(
        labels=labels.copy(),
        values=values,
        metric=metric,
        kmer_length=kmer_length,
    )


def _bounded_nonnegative_value(value: float) -> float:
    if value < 0.0:
        raise ValueError("Expected a non-negative descriptor value.")
    return value / (1.0 + value)


def _validate_k_values(k_values: list[int]) -> None:
    if not isinstance(k_values, list):
        raise TypeError("k_values must be a list of integers.")
    if not k_values:
        raise ValueError("k-mer lengths cannot be empty.")
    if len(k_values) != len(set(k_values)):
        raise ValueError("k-mer lengths must be unique.")
    for k in k_values:
        if type(k) is not int:
            raise TypeError("k-mer lengths must be integers.")
        if k <= 0:
            raise ValueError("k-mer lengths must be positive.")
