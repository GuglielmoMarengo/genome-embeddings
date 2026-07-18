import pytest
from src.genome import Genome

@pytest.fixture
def genome():
    return Genome (sequence="ACGTACGT")

@pytest.fixture
def aagc_genome():
    return Genome (sequence="AAGC")

@pytest.fixture
def all_a_genome():
    return Genome (sequence="AAAAAAAA")

@pytest.fixture
def all_gc_genome():
    return Genome (sequence="GCCGGGCC")

def test_length(genome):
    assert genome.length() == 8

def test_gc_content(genome):
    assert genome.gc_content() == 0.5

def test_reverse_complement(aagc_genome):
    assert aagc_genome.reverse_complement() == "GCTT"

def test_gc_content_all_a(all_a_genome):
    assert all_a_genome.gc_content() == 0

def test_gc_content_all_gc(all_gc_genome):
    assert all_gc_genome.gc_content() == 1

def test_value_error():
    with pytest.raises(ValueError, match="Invalid character X in sequence at position 4."):
        Genome(sequence="ACGX")

def test_empty_sequence_raises_value_error():
    with pytest.raises(ValueError, match="Sequence cannot be empty."):
        Genome(sequence="")

def test_from_fasta_reads_sequence():
    genome = Genome.from_fasta("tests/data/example.fasta")
    assert genome.sequence == "ACGTTGCA"
    assert genome.header == ">Example sequence"