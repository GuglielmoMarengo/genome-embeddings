from pathlib import Path

import pytest

from src.genome import Genome, GenomeComparison, GenomeDescriptor


TEST_DATA_DIR = Path(__file__).parent / "data"


def make_genome(sequence, header=None):
    return Genome(sequence=sequence, header=header)


def test_length():
    assert make_genome("ACGTACGT").length() == 8


@pytest.mark.parametrize(
    "sequence, expected_gc",
    [
        ("ACGTACGT", 0.5),
        ("AAAAAAAA", 0.0),
        ("GCCGGGCC", 1.0),
    ],
)
def test_gc_content(sequence, expected_gc):
    assert make_genome(sequence).gc_content() == expected_gc


def test_reverse_complement():
    assert make_genome("AAGC").reverse_complement() == "GCTT"


@pytest.mark.parametrize(
    "sequence, message",
    [
        ("ACGX", "Invalid character X in sequence at position 4."),
        ("", "Sequence cannot be empty."),
    ],
)
def test_invalid_sequence(sequence, message):
    with pytest.raises(ValueError, match=message):
        Genome(sequence=sequence)


def test_from_fasta_reads_sequence():
    genome = Genome.from_fasta(TEST_DATA_DIR / "example.fasta")
    assert genome.sequence == "ACGTTGCA"
    assert genome.header == ">Example sequence"


def test_kmers():
    assert make_genome("ACGTAC").kmers(3) == ["ACG", "CGT", "GTA", "TAC"]


def test_kmer_frequencies():
    assert make_genome("ACGTAC").kmer_frequencies(3) == {
        "ACG": 1,
        "CGT": 1,
        "GTA": 1,
        "TAC": 1,
    }


def test_nucleotide_frequencies():
    assert make_genome("ACGTAC").nucleotide_frequencies() == {
        "A": 2,
        "C": 2,
        "G": 1,
        "T": 1,
    }


def test_shannon_entropy_single_nucleotide():
    genome = make_genome("AAAAAAAAAA")
    assert genome.shannon_entropy() == pytest.approx(0.0)


def test_shannon_entropy_uniform_distribution():
    genome = make_genome("ACGT")
    assert genome.shannon_entropy() == pytest.approx(2.0)


def test_descriptor():
    genome = make_genome("ACGTACGT")
    descriptor = genome.descriptor()

    assert isinstance(descriptor, GenomeDescriptor)
    assert descriptor.length == 8
    assert descriptor.gc_content == 0.5
    assert descriptor.at_content == 0.5
    assert descriptor.shannon_entropy == pytest.approx(2.0)
    assert descriptor.gc_skew == 0.0
    assert descriptor.purine_content == 0.5
    assert descriptor.pyrimidine_content == 0.5
    assert descriptor.kmer_length == 3
    assert descriptor.kmer_diversity == pytest.approx(4 / 6)
    assert descriptor.kmer_entropy == pytest.approx(1.9182958340544896)

    assert descriptor.to_vector() == pytest.approx(
        [8.0, 0.5, 0.5, 2.0, 0.0, 0.5, 0.5, 3.0, 4 / 6, 1.9182958340544896]
    )

def test_at_content():
    assert make_genome("ACGTACGT").at_content() == 0.5

def test_gc_skew():
    assert make_genome("GGGC").gc_skew() == pytest.approx(0.5)

def test_gc_skew_without_gc_bases():
    assert make_genome("AAAA").gc_skew() == pytest.approx(0.0)

def test_purine_content():
    assert make_genome("AAGC").purine_content() == pytest.approx(0.75)

def test_pyrimidine_content():
    assert make_genome("ACTT").pyrimidine_content() == pytest.approx(0.75)

def test_purine_and_pyrimidine_content_sum_to_one():
    genome = make_genome("ACGTACGT")

    assert (
        genome.purine_content() + genome.pyrimidine_content()
        == pytest.approx(1.0)
    )

def test_kmer_diversity():
    genome = make_genome("ACGTAC")

    assert genome.kmer_diversity(3) == pytest.approx(1.0)

def test_kmer_diversity_with_repeated_kmers():
    genome = make_genome("AAAAA")

    assert genome.kmer_diversity(2) == pytest.approx(0.25)

def test_kmer_diversity_uses_theoretical_maximum():
    genome = make_genome("ACGT" * 20)

    assert genome.kmer_diversity(1) == pytest.approx(1.0)

def test_kmer_entropy_uniform_distribution():
    genome = make_genome("ACGTAC")

    assert genome.kmer_entropy(3) == pytest.approx(2.0)

def test_kmer_entropy_single_kmer():
    genome = make_genome("AAAAA")

    assert genome.kmer_entropy(2) == pytest.approx(0.0)

def test_descriptor_normalized_vector():
    genome = make_genome("ACGTACGT")
    descriptor = genome.descriptor(k=3)

    normalized_vector = descriptor.to_normalized_vector()

    assert len(normalized_vector) == 6
    assert all(0.0 <= value <= 1.0 for value in normalized_vector)

def test_euclidean_distance_between_identical_descriptors():
    genome = make_genome("ACGTACGT")
    descriptor = genome.descriptor(k=3)

    assert descriptor.euclidean_distance(descriptor) == pytest.approx(0.0)

def test_euclidean_distance_between_different_descriptors():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    assert first.euclidean_distance(second) > 0.0

def test_euclidean_distance_is_symmetric():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    assert first.euclidean_distance(second) == pytest.approx(
        second.euclidean_distance(first)
    )

def test_euclidean_distance_rejects_invalid_type():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    with pytest.raises(TypeError, match="other must be a GenomeDescriptor."):
        descriptor.euclidean_distance("invalid")
    
def test_cosine_similarity_identical_descriptors():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    assert descriptor.cosine_similarity(descriptor) == pytest.approx(1.0)

def test_cosine_similarity_is_symmetric():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    assert first.cosine_similarity(second) == pytest.approx(
        second.cosine_similarity(first)
    )

def test_cosine_similarity_rejects_invalid_type():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    with pytest.raises(TypeError, match="other must be a GenomeDescriptor."):
        descriptor.cosine_similarity("invalid")

def test_cosine_similarity_is_between_zero_and_one():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    similarity = first.cosine_similarity(second)

    assert 0.0 <= similarity <= 1.0

def test_feature_differences_identical_descriptors_are_zero():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    differences = descriptor.feature_differences(descriptor)

    assert all(value == pytest.approx(0.0) for value in differences.values())

def test_feature_differences_contains_normalized_features():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    differences = first.feature_differences(second)

    assert set(differences) == {
        "gc_content",
        "normalized_shannon_entropy",
        "normalized_gc_skew",
        "purine_content",
        "kmer_diversity",
        "normalized_kmer_entropy",
    }

def test_compare_returns_genome_comparison():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    comparison = first.compare(second)

    assert isinstance(comparison, GenomeComparison)
    assert comparison.euclidean_distance > 0.0
    assert 0.0 <= comparison.cosine_similarity <= 1.0
    assert comparison.feature_differences == first.feature_differences(second)

def test_compare_rejects_invalid_type():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    with pytest.raises(TypeError, match="other must be a GenomeDescriptor."):
        descriptor.compare("invalid")
    
def test_sorted_feature_differences_are_descending():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    comparison = first.compare(second)
    sorted_differences = comparison.sorted_feature_differences()

    values = [value for _, value in sorted_differences]

    assert values == sorted(values, reverse=True)