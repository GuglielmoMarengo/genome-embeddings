import math

import pytest

from src.descriptor_v2 import (
    DESCRIPTOR_V2_FEATURE_NAMES,
    DINUCLEOTIDE_ORDER,
    DescriptorV2Collection,
    GenomeDescriptorV2,
    conditional_nucleotide_entropy,
    dinucleotide_odds_ratios,
    kmer_scale_descriptor,
    multiscale_embedding,
)
from src.genome import Genome


def test_periodic_sequence_has_zero_first_order_conditional_entropy():
    genome = Genome("ACGTACGTACGT")

    assert conditional_nucleotide_entropy(genome) == pytest.approx(0.0)


def test_conditional_entropy_detects_branching_transitions():
    genome = Genome("AACAGATAC")

    assert conditional_nucleotide_entropy(genome) > 0.0
    assert conditional_nucleotide_entropy(genome) <= 2.0


def test_dinucleotide_odds_ratios_have_deterministic_coordinates():
    ratios = dinucleotide_odds_ratios(Genome("ACGTACGT"))

    assert tuple(ratios) == DINUCLEOTIDE_ORDER
    assert ratios["AC"] > 1.0
    assert ratios["AA"] == pytest.approx(0.0)


def test_finite_sample_entropy_uses_observable_maximum():
    genome = Genome("ACGTACGT")
    descriptor = kmer_scale_descriptor(genome, k=5)

    assert descriptor.window_count == 4
    assert descriptor.possible_kmer_count == 1024
    assert descriptor.observable_kmer_count == 4
    assert descriptor.distinct_kmer_count == 4
    assert descriptor.finite_sample_normalized_kmer_entropy == pytest.approx(
        1.0
    )
    assert descriptor.theoretical_normalized_kmer_entropy == pytest.approx(
        0.2
    )


def test_singleton_and_repeated_window_fractions_are_complementary():
    descriptor = kmer_scale_descriptor(Genome("AAAAAC"), k=2)

    assert descriptor.singleton_fraction + descriptor.repeated_window_fraction == (
        pytest.approx(1.0)
    )
    assert descriptor.repeated_window_fraction > 0.0


def test_effective_kmer_count_matches_entropy_definition():
    descriptor = kmer_scale_descriptor(Genome("ACGTACGT"), k=2)

    assert descriptor.effective_kmer_count == pytest.approx(
        2.0 ** descriptor.kmer_entropy
    )
    assert descriptor.effective_kmer_fraction <= 1.0


def test_descriptor_v2_vector_has_stable_feature_order():
    descriptor = GenomeDescriptorV2.from_genome(Genome("ACGTACGT"), k=2)

    assert tuple(descriptor.to_feature_dict()) == DESCRIPTOR_V2_FEATURE_NAMES
    assert len(descriptor.to_normalized_vector()) == len(
        DESCRIPTOR_V2_FEATURE_NAMES
    )
    assert all(0.0 <= value <= 1.0 for value in descriptor.to_normalized_vector())


def test_descriptor_v2_distances_are_zero_for_identical_sequences():
    first = GenomeDescriptorV2.from_genome(Genome("ACGTACGT"), k=3)
    second = GenomeDescriptorV2.from_genome(Genome("ACGTACGT"), k=3)

    assert first.euclidean_distance(second) == pytest.approx(0.0)
    assert first.cosine_similarity(second) == pytest.approx(1.0)


def test_descriptor_v2_detects_order_with_equal_nucleotide_composition():
    periodic = GenomeDescriptorV2.from_genome(Genome("ACGT" * 8), k=2)
    blocked = GenomeDescriptorV2.from_genome(
        Genome("A" * 8 + "C" * 8 + "G" * 8 + "T" * 8),
        k=2,
    )

    assert periodic.gc_content == pytest.approx(blocked.gc_content)
    assert periodic.euclidean_distance(blocked) > 0.0
    assert periodic.feature_differences(blocked)[
        "normalized_conditional_nucleotide_entropy"
    ] > 0.0


def test_multiscale_embedding_includes_global_features_once():
    embedding = multiscale_embedding(Genome("ACGT" * 20), [1, 2, 3])

    assert embedding.k_values == (1, 2, 3)
    assert embedding.feature_names.count("gc_content") == 1
    assert len(embedding.values) == 21 + 5 * 3


def test_multiscale_embedding_rejects_duplicate_scales():
    with pytest.raises(
        ValueError,
        match=r"k-mer lengths must be unique\.",
    ):
        multiscale_embedding(Genome("ACGTACGT"), [1, 2, 2])


def test_descriptor_v2_collection_builds_supported_matrices():
    collection = DescriptorV2Collection(
        [Genome("ACGTACGT"), Genome("AAAACCCC")]
    )
    labels = ["First", "Second"]

    matrix = collection.euclidean_distance_matrix(labels=labels, k=2)
    embedding_matrix = collection.multiscale_embedding_distance_matrix(
        labels=labels,
        k_values=[1, 2, 3],
    )

    assert matrix.metric == "euclidean_v2"
    assert embedding_matrix.metric == "embedding_v2_euclidean"
    assert matrix.values[0][1] > 0.0
    assert embedding_matrix.values[0][1] > 0.0


def test_descriptor_v2_cosine_matrix_ranks_larger_similarity_first():
    collection = DescriptorV2Collection(
        [
            Genome("ACGTACGT"),
            Genome("ACGTACGA"),
            Genome("AAAACCCC"),
        ]
    )
    matrix = collection.cosine_similarity_matrix(
        labels=["Reference", "Near", "Far"],
        k=2,
    )

    assert matrix.rank_by_label("Reference")[0][0] == "Near"


def test_embedding_json_is_serializable():
    embedding = multiscale_embedding(Genome("ACGT" * 10), [1, 2])

    output = embedding.to_json(indent=2)

    assert '"version": 2' in output
    assert '"k_values"' in output
