import math

import pytest

from src.genome import Genome
from src.kmer_distribution import (
    KmerDistributionCollection,
    KmerDistributionProfile,
    jensen_shannon_divergence,
    kmer_vocabulary,
)


def test_kmer_vocabulary_is_lexicographic_and_complete():
    vocabulary = kmer_vocabulary(2)

    assert vocabulary[:4] == ("AA", "AC", "AG", "AT")
    assert vocabulary[-1] == "TT"
    assert len(vocabulary) == 16


def test_profile_contains_complete_probability_vector():
    profile = KmerDistributionProfile.from_genome(
        Genome("AAAA"),
        k=2,
    )

    assert len(profile.probabilities) == 16
    assert sum(profile.probabilities) == pytest.approx(1.0)
    assert profile.probabilities[profile.vocabulary.index("AA")] == pytest.approx(
        1.0
    )


def test_jensen_shannon_divergence_is_zero_for_identical_distributions():
    distribution = [0.25, 0.25, 0.25, 0.25]

    assert jensen_shannon_divergence(distribution, distribution) == pytest.approx(
        0.0
    )


def test_jensen_shannon_divergence_is_one_for_disjoint_distributions():
    assert jensen_shannon_divergence([1.0, 0.0], [0.0, 1.0]) == pytest.approx(
        1.0
    )


def test_jensen_shannon_divergence_is_symmetric():
    first = [0.5, 0.5, 0.0]
    second = [0.0, 0.5, 0.5]

    assert jensen_shannon_divergence(first, second) == pytest.approx(
        jensen_shannon_divergence(second, first)
    )


def test_jensen_shannon_distance_is_square_root_of_divergence():
    first = KmerDistributionProfile.from_genome(Genome("AAAAAA"), k=1)
    second = KmerDistributionProfile.from_genome(Genome("CCCCCC"), k=1)

    divergence = first.jensen_shannon_divergence(second)

    assert first.jensen_shannon_distance(second) == pytest.approx(
        math.sqrt(divergence)
    )
    assert first.jensen_shannon_distance(second) == pytest.approx(1.0)


def test_probability_vectors_must_sum_to_one():
    with pytest.raises(
        ValueError,
        match=r"First probability vector must sum to one\.",
    ):
        jensen_shannon_divergence([0.2, 0.2], [0.5, 0.5])


def test_distribution_collection_builds_symmetric_distance_matrix():
    collection = KmerDistributionCollection(
        [Genome("ACGTACGT"), Genome("AAAACCCC")]
    )
    matrix = collection.distance_matrix(
        labels=["First", "Second"],
        k=2,
    )

    assert matrix.metric == "jensen_shannon"
    assert matrix.values[0][0] == pytest.approx(0.0)
    assert matrix.values[0][1] == pytest.approx(matrix.values[1][0])
    assert matrix.values[0][1] > 0.0


def test_distribution_matrix_ranking_uses_smaller_distance_first():
    collection = KmerDistributionCollection(
        [
            Genome("ACGTACGT"),
            Genome("ACGTACGA"),
            Genome("AAAAAAAC"),
        ]
    )
    matrix = collection.distance_matrix(
        labels=["Reference", "Near", "Far"],
        k=2,
    )

    assert matrix.rank_by_label("Reference")[0][0] == "Near"


def test_distribution_collection_preserves_scale_order():
    collection = KmerDistributionCollection(
        [Genome("ACGTACGT"), Genome("AAAACCCC")]
    )

    matrices = collection.distance_matrices(
        labels=["First", "Second"],
        k_values=[3, 1, 2],
    )

    assert list(matrices) == [3, 1, 2]


def test_distribution_trajectory_step_distances_are_available():
    collection = KmerDistributionCollection(
        [Genome("ACGTACGT"), Genome("AAAACCCC"), Genome("GGGGTTTT")]
    )

    distances = collection.matrix_trajectory_step_distances(
        labels=["A", "B", "C"],
        k_values=[1, 2, 3],
    )

    assert list(distances) == [(1, 2), (2, 3)]
    assert all(value >= 0.0 for value in distances.values())


def test_distribution_ranking_stability_uses_existing_core_analysis():
    collection = KmerDistributionCollection(
        [Genome("ACGTACGT"), Genome("ACGTACGA"), Genome("AAAACCCC")]
    )

    stability = collection.ranking_stability(
        labels=["Reference", "Near", "Far"],
        reference_label="Reference",
        k_values=[1, 2, 3],
    )

    assert list(stability) == [(1, 2), (2, 3)]
    assert all(-1.0 <= row["kendall_tau"] <= 1.0 for row in stability.values())
