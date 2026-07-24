"""Full k-mer probability profiles and Jensen-Shannon comparisons."""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product

from src.genome import Genome, GenomeCollection, GenomeMatrix


NUCLEOTIDE_ORDER = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class KmerDistributionProfile:
    kmer_length: int
    vocabulary: tuple[str, ...]
    probabilities: tuple[float, ...]

    def __post_init__(self) -> None:
        if self.kmer_length <= 0:
            raise ValueError("k-mer length must be positive.")
        if not self.vocabulary:
            raise ValueError("k-mer vocabulary cannot be empty.")
        if len(self.vocabulary) != len(self.probabilities):
            raise ValueError("k-mer vocabulary must match probabilities.")
        if any(probability < 0.0 for probability in self.probabilities):
            raise ValueError("k-mer probabilities cannot be negative.")
        if not math.isclose(
            sum(self.probabilities),
            1.0,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise ValueError("k-mer probabilities must sum to one.")

    @classmethod
    def from_genome(
        cls,
        genome: Genome,
        k: int,
    ) -> "KmerDistributionProfile":
        if not isinstance(genome, Genome):
            raise TypeError("genome must be a Genome object.")

        vocabulary = kmer_vocabulary(k)
        frequencies = genome.kmer_frequencies(k)
        total = genome.length() - k + 1
        probabilities = tuple(
            frequencies.get(kmer, 0) / total
            for kmer in vocabulary
        )
        return cls(
            kmer_length=k,
            vocabulary=vocabulary,
            probabilities=probabilities,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "kmer_length": self.kmer_length,
            "probabilities": dict(
                zip(
                    self.vocabulary,
                    self.probabilities,
                    strict=True,
                )
            ),
        }

    def jensen_shannon_divergence(
        self,
        other: "KmerDistributionProfile",
    ) -> float:
        self._validate_compatible(other)
        return jensen_shannon_divergence(
            self.probabilities,
            other.probabilities,
        )

    def jensen_shannon_distance(
        self,
        other: "KmerDistributionProfile",
    ) -> float:
        return math.sqrt(self.jensen_shannon_divergence(other))

    def _validate_compatible(self, other: object) -> None:
        if not isinstance(other, KmerDistributionProfile):
            raise TypeError("other must be a KmerDistributionProfile.")
        if self.vocabulary != other.vocabulary:
            raise ValueError("k-mer profile vocabularies must match.")


class KmerDistributionCollection:
    """Collection-level Jensen-Shannon matrices and trajectories."""

    def __init__(self, genomes: list[Genome]) -> None:
        self.collection = GenomeCollection(genomes)

    @property
    def genomes(self) -> list[Genome]:
        return self.collection.genomes

    def profiles(self, k: int) -> list[KmerDistributionProfile]:
        return [
            KmerDistributionProfile.from_genome(genome, k)
            for genome in self.genomes
        ]

    def distance_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        if len(labels) != len(self.genomes):
            raise ValueError("Genome labels must match collection size.")
        profiles = self.profiles(k)
        values = [
            [
                first.jensen_shannon_distance(second)
                for second in profiles
            ]
            for first in profiles
        ]
        return GenomeMatrix(
            labels=labels.copy(),
            values=values,
            metric="jensen_shannon",
            kmer_length=k,
        )

    def distance_matrices(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, GenomeMatrix]:
        _validate_k_values(k_values)
        return {
            k: self.distance_matrix(labels=labels, k=k)
            for k in k_values
        }

    def matrix_trajectory(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, list[float]]:
        matrices = self.distance_matrices(labels, k_values)
        return {
            k: matrix.to_upper_triangle_vector()
            for k, matrix in matrices.items()
        }

    def pair_trajectory(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
    ) -> dict[int, float]:
        matrices = self.distance_matrices(labels, k_values)
        return {
            k: matrix.get_value(row_label, column_label)
            for k, matrix in matrices.items()
        }

    def matrix_trajectory_step_distances(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[tuple[int, int], float]:
        return GenomeCollection.matrix_trajectory_step_distances(
            self.matrix_trajectory(labels, k_values)
        )

    def ranking_trajectory(
        self,
        labels: list[str],
        reference_label: str,
        k_values: list[int],
    ) -> dict[int, list[tuple[str, float]]]:
        return GenomeCollection.matrix_ranking_trajectory(
            matrices=self.distance_matrices(labels, k_values),
            reference_label=reference_label,
        )

    def ranking_stability(
        self,
        labels: list[str],
        reference_label: str,
        k_values: list[int],
    ) -> dict[tuple[int, int], dict[str, float | int | bool]]:
        return GenomeCollection.ranking_trajectory_stability(
            self.ranking_trajectory(
                labels=labels,
                reference_label=reference_label,
                k_values=k_values,
            )
        )


def kmer_vocabulary(k: int) -> tuple[str, ...]:
    if type(k) is not int:
        raise TypeError("k must be an integer.")
    if k <= 0:
        raise ValueError("k must be positive.")
    return tuple(
        "".join(symbols)
        for symbols in product(NUCLEOTIDE_ORDER, repeat=k)
    )


def jensen_shannon_divergence(
    first: tuple[float, ...] | list[float],
    second: tuple[float, ...] | list[float],
) -> float:
    """Return base-2 Jensen-Shannon divergence in the interval [0, 1]."""

    if not first or not second:
        raise ValueError("Probability vectors cannot be empty.")
    if len(first) != len(second):
        raise ValueError("Probability vectors must have the same length.")
    if (
        any(value < 0.0 for value in first)
        or any(value < 0.0 for value in second)
    ):
        raise ValueError("Probability vectors cannot contain negatives.")
    if not math.isclose(sum(first), 1.0, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("First probability vector must sum to one.")
    if not math.isclose(sum(second), 1.0, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("Second probability vector must sum to one.")

    midpoint = tuple(
        (first_value + second_value) / 2.0
        for first_value, second_value in zip(first, second, strict=True)
    )

    return 0.5 * _kl_divergence(first, midpoint) + 0.5 * _kl_divergence(
        second,
        midpoint,
    )


def _kl_divergence(
    distribution: tuple[float, ...] | list[float],
    reference: tuple[float, ...],
) -> float:
    return sum(
        probability * math.log2(probability / reference_probability)
        for probability, reference_probability in zip(
            distribution,
            reference,
            strict=True,
        )
        if probability > 0.0
    )


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
