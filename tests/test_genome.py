from pathlib import Path

import pytest

from src.genome import Genome, GenomeDescriptor


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
    assert descriptor.to_dict() == {
        "length": 8,
        "gc_content": 0.5,
        "shannon_entropy": 2.0,
    }
    assert descriptor.to_vector() == [8.0, 0.5, 2.0]