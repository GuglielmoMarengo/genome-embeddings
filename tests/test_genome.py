from pathlib import Path

import pytest, json, csv, io

from src.genome import (
    Genome,
    GenomeCollection,
    GenomeComparison,
    GenomeDescriptor,
    GenomeMatrix,
)


TEST_DATA_DIR = Path(__file__).parent / "data"

THREE_GENOME_LABELS = [
    "Balanced",
    "Adenine",
    "Cytosine",
]

TWO_GENOME_LABELS = [
    "First",
    "Second",
]


def make_genome(
    sequence: str,
    header: str | None = None,
) -> Genome:
    return Genome(
        sequence=sequence,
        header=header,
    )


@pytest.fixture
def two_genomes() -> list[Genome]:
    return [
        make_genome("ACGTACGT"),
        make_genome("AAAAAAAA"),
    ]


@pytest.fixture
def three_genomes() -> list[Genome]:
    return [
        make_genome("ACGTACGT"),
        make_genome("AAAAAAAA"),
        make_genome("CCCCCCCC"),
    ]


@pytest.fixture
def two_genome_collection(
    two_genomes: list[Genome],
) -> GenomeCollection:
    return GenomeCollection(two_genomes)


@pytest.fixture
def three_genome_collection(
    three_genomes: list[Genome],
) -> GenomeCollection:
    return GenomeCollection(three_genomes)


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
def test_gc_content(
    sequence: str,
    expected_gc: float,
):
    assert make_genome(sequence).gc_content() == expected_gc


def test_at_content():
    assert make_genome("ACGTACGT").at_content() == pytest.approx(0.5)


def test_reverse_complement():
    assert make_genome("AAGC").reverse_complement() == "GCTT"


@pytest.mark.parametrize(
    "sequence, message",
    [
        (
            "ACGX",
            r"Invalid character X in sequence at position 4\.",
        ),
        (
            "",
            r"Sequence cannot be empty\.",
        ),
    ],
)
def test_invalid_sequence(
    sequence: str,
    message: str,
):
    with pytest.raises(ValueError, match=message):
        Genome(sequence=sequence)


def test_from_fasta_reads_sequence():
    genome = Genome.from_fasta(
        TEST_DATA_DIR / "example.fasta"
    )

    assert genome.sequence == "ACGTTGCA"
    assert genome.header == ">Example sequence"


def test_kmers():
    assert make_genome("ACGTAC").kmers(3) == [
        "ACG",
        "CGT",
        "GTA",
        "TAC",
    ]


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


def test_gc_skew():
    assert make_genome("GGGC").gc_skew() == pytest.approx(0.5)


def test_gc_skew_without_gc_bases():
    assert make_genome("AAAA").gc_skew() == pytest.approx(0.0)


def test_purine_content():
    assert make_genome("AAGC").purine_content() == pytest.approx(
        0.75
    )


def test_pyrimidine_content():
    assert make_genome("ACTT").pyrimidine_content() == pytest.approx(
        0.75
    )


def test_purine_and_pyrimidine_content_sum_to_one():
    genome = make_genome("ACGTACGT")

    total = (
        genome.purine_content()
        + genome.pyrimidine_content()
    )

    assert total == pytest.approx(1.0)


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


def test_descriptor():
    genome = make_genome("ACGTACGT")
    descriptor = genome.descriptor()

    assert isinstance(descriptor, GenomeDescriptor)
    assert descriptor.length == 8
    assert descriptor.gc_content == pytest.approx(0.5)
    assert descriptor.at_content == pytest.approx(0.5)
    assert descriptor.shannon_entropy == pytest.approx(2.0)
    assert descriptor.gc_skew == pytest.approx(0.0)
    assert descriptor.purine_content == pytest.approx(0.5)
    assert descriptor.pyrimidine_content == pytest.approx(0.5)
    assert descriptor.kmer_length == 3
    assert descriptor.kmer_diversity == pytest.approx(4 / 6)
    assert descriptor.kmer_entropy == pytest.approx(
        1.9182958340544896
    )

    assert descriptor.to_vector() == pytest.approx(
        [
            8.0,
            0.5,
            0.5,
            2.0,
            0.0,
            0.5,
            0.5,
            3.0,
            4 / 6,
            1.9182958340544896,
        ]
    )


def test_descriptor_normalized_vector():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    normalized_vector = descriptor.to_normalized_vector()

    assert len(normalized_vector) == 6
    assert all(
        0.0 <= value <= 1.0
        for value in normalized_vector
    )


def test_euclidean_distance_between_identical_descriptors():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    assert descriptor.euclidean_distance(
        descriptor
    ) == pytest.approx(0.0)


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

    with pytest.raises(
        TypeError,
        match=r"other must be a GenomeDescriptor\.",
    ):
        descriptor.euclidean_distance("invalid")


def test_cosine_similarity_identical_descriptors():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    assert descriptor.cosine_similarity(
        descriptor
    ) == pytest.approx(1.0)


def test_cosine_similarity_is_symmetric():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    assert first.cosine_similarity(second) == pytest.approx(
        second.cosine_similarity(first)
    )


def test_cosine_similarity_rejects_invalid_type():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    with pytest.raises(
        TypeError,
        match=r"other must be a GenomeDescriptor\.",
    ):
        descriptor.cosine_similarity("invalid")


def test_cosine_similarity_is_between_zero_and_one():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    similarity = first.cosine_similarity(second)

    assert 0.0 <= similarity <= 1.0


def test_feature_differences_identical_descriptors_are_zero():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    differences = descriptor.feature_differences(descriptor)

    assert all(
        value == pytest.approx(0.0)
        for value in differences.values()
    )


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
    assert comparison.feature_differences == (
        first.feature_differences(second)
    )


def test_compare_rejects_invalid_type():
    descriptor = make_genome("ACGTACGT").descriptor(k=3)

    with pytest.raises(
        TypeError,
        match=r"other must be a GenomeDescriptor\.",
    ):
        descriptor.compare("invalid")


def test_sorted_feature_differences_are_descending():
    first = make_genome("ACGTACGT").descriptor(k=3)
    second = make_genome("AAAAAAAA").descriptor(k=3)

    comparison = first.compare(second)
    sorted_differences = (
        comparison.sorted_feature_differences()
    )

    values = [
        value
        for _, value in sorted_differences
    ]

    assert values == sorted(values, reverse=True)


def test_genome_collection_stores_genomes(
    two_genomes: list[Genome],
):
    collection = GenomeCollection(two_genomes)

    assert collection.genomes == two_genomes


def test_genome_collection_rejects_empty_collection():
    with pytest.raises(
        ValueError,
        match=r"Genome collection cannot be empty\.",
    ):
        GenomeCollection([])


def test_genome_collection_rejects_non_list_input():
    with pytest.raises(
        TypeError,
        match=r"genomes must be a list of Genome objects\.",
    ):
        GenomeCollection("invalid")


def test_genome_collection_rejects_invalid_items():
    genomes = [
        make_genome("ACGTACGT"),
        "invalid",
    ]

    with pytest.raises(
        TypeError,
        match=r"All items must be Genome objects\.",
    ):
        GenomeCollection(genomes)


def test_genome_collection_generates_descriptors(
    two_genome_collection: GenomeCollection,
):
    descriptors = two_genome_collection.descriptors(k=3)

    assert len(descriptors) == 2
    assert all(
        isinstance(descriptor, GenomeDescriptor)
        for descriptor in descriptors
    )


def test_genome_collection_preserves_descriptor_order(
    two_genomes: list[Genome],
):
    collection = GenomeCollection(two_genomes)

    descriptors = collection.descriptors(k=3)

    assert descriptors[0] == two_genomes[0].descriptor(k=3)
    assert descriptors[1] == two_genomes[1].descriptor(k=3)


def test_genome_collection_descriptors_reject_invalid_k(
    two_genome_collection: GenomeCollection,
):
    with pytest.raises(
        ValueError,
        match=r"k must be positive\. Got 0\.",
    ):
        two_genome_collection.descriptors(k=0)


def test_genome_matrix_stores_values_and_metadata():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.labels == ["Genome A", "Genome B"]
    assert matrix.values == [
        [0.0, 0.5],
        [0.5, 0.0],
    ]
    assert matrix.metric == "euclidean"
    assert matrix.kmer_length == 3


def test_genome_matrix_rejects_empty_labels():
    with pytest.raises(
        ValueError,
        match=r"Genome matrix labels cannot be empty\.",
    ):
        GenomeMatrix(
            labels=[],
            values=[],
            metric="euclidean",
            kmer_length=3,
        )


def test_genome_matrix_rejects_non_square_values():
    with pytest.raises(
        ValueError,
        match=r"Genome matrix values must be square\.",
    ):
        GenomeMatrix(
            labels=["Genome A", "Genome B"],
            values=[
                [0.0, 0.5],
            ],
            metric="euclidean",
            kmer_length=3,
        )


def test_genome_matrix_rejects_label_value_size_mismatch():
    with pytest.raises(
        ValueError,
        match=(
            r"Genome matrix labels must match matrix size\."
        ),
    ):
        GenomeMatrix(
            labels=["Genome A"],
            values=[
                [0.0, 0.5],
                [0.5, 0.0],
            ],
            metric="euclidean",
            kmer_length=3,
        )


def test_genome_matrix_rejects_invalid_metric():
    with pytest.raises(
        ValueError,
        match=r"Unsupported genome matrix metric\.",
    ):
        GenomeMatrix(
            labels=["Genome A"],
            values=[[0.0]],
            metric="invalid",
            kmer_length=3,
        )


def test_euclidean_distance_matrix_returns_genome_matrix(
    two_genome_collection: GenomeCollection,
):
    matrix = (
        two_genome_collection.euclidean_distance_matrix(
            labels=TWO_GENOME_LABELS,
            k=3,
        )
    )

    assert isinstance(matrix, GenomeMatrix)
    assert matrix.labels == TWO_GENOME_LABELS
    assert matrix.metric == "euclidean"
    assert matrix.kmer_length == 3


def test_euclidean_distance_matrix_has_correct_shape(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.euclidean_distance_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    assert len(matrix.values) == 3
    assert all(
        len(row) == 3
        for row in matrix.values
    )


def test_euclidean_distance_matrix_has_zero_diagonal(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.euclidean_distance_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    for index in range(len(matrix.values)):
        assert matrix.values[index][index] == pytest.approx(
            0.0
        )


def test_euclidean_distance_matrix_is_symmetric(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.euclidean_distance_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    for row_index in range(len(matrix.values)):
        for column_index in range(len(matrix.values)):
            assert (
                matrix.values[row_index][column_index]
                == pytest.approx(
                    matrix.values[column_index][row_index]
                )
            )


def test_euclidean_distance_matrix_matches_descriptor_distance(
    two_genomes: list[Genome],
):
    collection = GenomeCollection(two_genomes)

    matrix = collection.euclidean_distance_matrix(
        labels=TWO_GENOME_LABELS,
        k=3,
    )

    expected_distance = (
        two_genomes[0]
        .descriptor(k=3)
        .euclidean_distance(
            two_genomes[1].descriptor(k=3)
        )
    )

    assert matrix.values[0][1] == pytest.approx(
        expected_distance
    )
    assert matrix.values[1][0] == pytest.approx(
        expected_distance
    )


def test_cosine_similarity_matrix_returns_genome_matrix(
    two_genome_collection: GenomeCollection,
):
    matrix = (
        two_genome_collection.cosine_similarity_matrix(
            labels=TWO_GENOME_LABELS,
            k=3,
        )
    )

    assert isinstance(matrix, GenomeMatrix)
    assert matrix.labels == TWO_GENOME_LABELS
    assert matrix.metric == "cosine"
    assert matrix.kmer_length == 3


def test_cosine_similarity_matrix_has_correct_shape(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.cosine_similarity_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    assert len(matrix.values) == 3
    assert all(
        len(row) == 3
        for row in matrix.values
    )


def test_cosine_similarity_matrix_has_one_diagonal(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.cosine_similarity_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    for index in range(len(matrix.values)):
        assert matrix.values[index][index] == pytest.approx(
            1.0
        )


def test_cosine_similarity_matrix_is_symmetric(
    three_genome_collection: GenomeCollection,
):
    matrix = (
        three_genome_collection.cosine_similarity_matrix(
            labels=THREE_GENOME_LABELS,
            k=3,
        )
    )

    for row_index in range(len(matrix.values)):
        for column_index in range(len(matrix.values)):
            assert (
                matrix.values[row_index][column_index]
                == pytest.approx(
                    matrix.values[column_index][row_index]
                )
            )


def test_cosine_similarity_matrix_matches_descriptor_similarity(
    two_genomes: list[Genome],
):
    collection = GenomeCollection(two_genomes)

    matrix = collection.cosine_similarity_matrix(
        labels=TWO_GENOME_LABELS,
        k=3,
    )

    expected_similarity = (
        two_genomes[0]
        .descriptor(k=3)
        .cosine_similarity(
            two_genomes[1].descriptor(k=3)
        )
    )

    assert matrix.values[0][1] == pytest.approx(
        expected_similarity
    )
    assert matrix.values[1][0] == pytest.approx(
        expected_similarity
    )


@pytest.mark.parametrize(
    "matrix_method_name",
    [
        "euclidean_distance_matrix",
        "cosine_similarity_matrix",
    ],
)
def test_collection_matrix_rejects_label_count_mismatch(
    two_genome_collection: GenomeCollection,
    matrix_method_name: str,
):
    matrix_method = getattr(
        two_genome_collection,
        matrix_method_name,
    )

    with pytest.raises(
        ValueError,
        match=(
            r"Genome matrix labels must match matrix size\."
        ),
    ):
        matrix_method(
            labels=["Only one label"],
            k=3,
        )

def test_genome_matrix_get_value_by_labels():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    value = matrix.get_value(
        row_label="Genome A",
        column_label="Genome B",
    )

    assert value == pytest.approx(0.5)

def test_genome_matrix_get_value_preserves_direction():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [1.0, 0.8],
            [0.7, 1.0],
        ],
        metric="cosine",
        kmer_length=3,
    )

    assert matrix.get_value(
        row_label="Genome A",
        column_label="Genome B",
    ) == pytest.approx(0.8)

    assert matrix.get_value(
        row_label="Genome B",
        column_label="Genome A",
    ) == pytest.approx(0.7)

def test_genome_matrix_get_value_rejects_unknown_row_label():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    with pytest.raises(
        ValueError,
        match=r"Unknown genome matrix label: Unknown\.",
    ):
        matrix.get_value(
            row_label="Unknown",
            column_label="Genome B",
        )

def test_genome_matrix_get_value_rejects_unknown_column_label():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    with pytest.raises(
        ValueError,
        match=r"Unknown genome matrix label: Unknown\.",
    ):
        matrix.get_value(
            row_label="Genome A",
            column_label="Unknown",
        )

def test_genome_matrix_to_rows():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.to_rows() == [
        {
            "label": "Genome A",
            "values": [0.0, 0.5],
        },
        {
            "label": "Genome B",
            "values": [0.5, 0.0],
        },
    ]

def test_genome_matrix_to_rows_returns_copies():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    rows = matrix.to_rows()
    rows[0]["values"][0] = 999.0

    assert matrix.values[0][0] == pytest.approx(0.0)

def test_genome_matrix_to_dict():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.to_dict() == {
        "labels": ["Genome A", "Genome B"],
        "values": [
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        "metric": "euclidean",
        "kmer_length": 3,
    }

def test_genome_matrix_to_dict_returns_copies():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    matrix_dict = matrix.to_dict()

    matrix_dict["labels"][0] = "Modified"
    matrix_dict["values"][0][0] = 999.0

    assert matrix.labels[0] == "Genome A"
    assert matrix.values[0][0] == pytest.approx(0.0)

def test_genome_matrix_rank_by_label_for_euclidean():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.1],
            [0.2, 0.0, 0.3],
            [0.1, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    ranking = matrix.rank_by_label(
        label="Genome A",
    )

    assert ranking == [
        ("Genome C", pytest.approx(0.1)),
        ("Genome B", pytest.approx(0.2)),
    ]

def test_genome_matrix_rank_by_label_for_cosine():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [1.0, 0.8, 0.95],
            [0.8, 1.0, 0.7],
            [0.95, 0.7, 1.0],
        ],
        metric="cosine",
        kmer_length=3,
    )

    ranking = matrix.rank_by_label(
        label="Genome A",
    )

    assert ranking == [
        ("Genome C", pytest.approx(0.95)),
        ("Genome B", pytest.approx(0.8)),
    ]

def test_genome_matrix_rank_by_label_excludes_reference():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
        ],
        values=[
            [0.0, 0.2],
            [0.2, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    ranking = matrix.rank_by_label(
        label="Genome A",
    )

    assert ranking == [
        ("Genome B", pytest.approx(0.2)),
    ]

def test_genome_matrix_rank_by_label_rejects_unknown_label():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
        ],
        values=[
            [0.0, 0.2],
            [0.2, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    with pytest.raises(
        ValueError,
        match=r"Unknown genome matrix label: Unknown\.",
    ):
        matrix.rank_by_label(
            label="Unknown",
        )

def test_genome_matrix_to_json():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    json_output = matrix.to_json()

    assert json.loads(json_output) == {
        "labels": ["Genome A", "Genome B"],
        "values": [
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        "metric": "euclidean",
        "kmer_length": 3,
    }

def test_genome_matrix_to_json_supports_indentation():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [1.0, 0.8],
            [0.8, 1.0],
        ],
        metric="cosine",
        kmer_length=3,
    )

    json_output = matrix.to_json(
        indent=2,
    )

    assert "\n" in json_output
    assert '  "labels"' in json_output

def test_genome_matrix_to_csv():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    csv_output = matrix.to_csv()

    rows = list(
        csv.reader(
            io.StringIO(csv_output)
        )
    )

    assert rows == [
        ["label", "Genome A", "Genome B"],
        ["Genome A", "0.0", "0.5"],
        ["Genome B", "0.5", "0.0"],
    ]

def test_genome_matrix_to_csv_supports_custom_delimiter():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[
            [1.0, 0.8],
            [0.8, 1.0],
        ],
        metric="cosine",
        kmer_length=3,
    )

    csv_output = matrix.to_csv(
        delimiter=";",
    )

    rows = list(
        csv.reader(
            io.StringIO(csv_output),
            delimiter=";",
        )
    )

    assert rows == [
        ["label", "Genome A", "Genome B"],
        ["Genome A", "1.0", "0.8"],
        ["Genome B", "0.8", "1.0"],
    ]

def test_genome_collection_generates_multi_k_euclidean_matrices():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.euclidean_distance_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[1, 2, 3],
    )

    assert list(matrices) == [1, 2, 3]

    assert all(
        isinstance(matrix, GenomeMatrix)
        for matrix in matrices.values()
    )

    assert matrices[1].kmer_length == 1
    assert matrices[2].kmer_length == 2
    assert matrices[3].kmer_length == 3

    assert all(
        matrix.metric == "euclidean"
        for matrix in matrices.values()
    )

def test_multi_k_euclidean_matrices_match_single_k_results():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.euclidean_distance_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[1, 2, 3],
    )

    for k in [1, 2, 3]:
        expected_matrix = (
            collection.euclidean_distance_matrix(
                labels=["Genome A", "Genome B"],
                k=k,
            )
        )

        assert matrices[k] == expected_matrix

def test_multi_k_euclidean_matrices_preserve_requested_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.euclidean_distance_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[3, 1, 2],
    )

    assert list(matrices) == [3, 1, 2]

def test_multi_k_euclidean_matrices_reject_empty_k_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths cannot be empty",
    ):
        collection.euclidean_distance_matrices(
            labels=["Genome A", "Genome B"],
            k_values=[],
        )

def test_genome_collection_generates_multi_k_cosine_matrices():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.cosine_similarity_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[1, 2, 3],
    )

    assert list(matrices) == [1, 2, 3]

    assert all(
        isinstance(matrix, GenomeMatrix)
        for matrix in matrices.values()
    )

    assert matrices[1].kmer_length == 1
    assert matrices[2].kmer_length == 2
    assert matrices[3].kmer_length == 3

    assert all(
        matrix.metric == "cosine"
        for matrix in matrices.values()
    )

def test_multi_k_cosine_matrices_match_single_k_results():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.cosine_similarity_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[1, 2, 3],
    )

    for k in [1, 2, 3]:
        expected_matrix = (
            collection.cosine_similarity_matrix(
                labels=["Genome A", "Genome B"],
                k=k,
            )
        )

        assert matrices[k] == expected_matrix

def test_multi_k_cosine_matrices_preserve_requested_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    matrices = collection.cosine_similarity_matrices(
        labels=["Genome A", "Genome B"],
        k_values=[3, 1, 2],
    )

    assert list(matrices) == [3, 1, 2]

def test_multi_k_cosine_matrices_reject_empty_k_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(

        ValueError,
        match="k-mer lengths cannot be empty",
    ):
        collection.cosine_similarity_matrices(
            labels=["Genome A", "Genome B"],
            k_values=[],
        )

def test_multi_k_euclidean_matrices_reject_duplicate_k_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.euclidean_distance_matrices(
            labels=["Genome A", "Genome B"],
            k_values=[1, 2, 2, 3],
        )

def test_multi_k_cosine_matrices_reject_duplicate_k_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.cosine_similarity_matrices(
            labels=["Genome A", "Genome B"],
            k_values=[1, 2, 2, 3],
        )

def test_genome_matrix_to_upper_triangle_vector():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.to_upper_triangle_vector() == [
        0.2,
        0.4,
        0.3,
    ]

def test_upper_triangle_vector_excludes_diagonal():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
        ],
        values=[
            [1.0, 0.8],
            [0.8, 1.0],
        ],
        metric="cosine",
        kmer_length=3,
    )

    assert matrix.to_upper_triangle_vector() == [
        0.8,
    ]

def test_upper_triangle_vector_returns_independent_list():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    vector = matrix.to_upper_triangle_vector()
    vector[0] = 999.0

    assert matrix.values[0][1] == 0.2

def test_genome_matrix_upper_triangle_pairs():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.upper_triangle_pairs() == [
        ("Genome A", "Genome B"),
        ("Genome A", "Genome C"),
        ("Genome B", "Genome C"),
    ]

def test_upper_triangle_pairs_match_vector_order():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    pairs = matrix.upper_triangle_pairs()
    vector = matrix.to_upper_triangle_vector()

    assert list(zip(pairs, vector, strict=True)) == [
        (("Genome A", "Genome B"), 0.2),
        (("Genome A", "Genome C"), 0.4),
        (("Genome B", "Genome C"), 0.3),
    ]

def test_upper_triangle_pairs_returns_independent_list():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
        ],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    pairs = matrix.upper_triangle_pairs()
    pairs.append(("Modified", "Pair"))

    assert matrix.upper_triangle_pairs() == [
        ("Genome A", "Genome B"),
    ]

def test_genome_matrix_to_upper_triangle_rows():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    assert matrix.to_upper_triangle_rows() == [
        {
            "row_label": "Genome A",
            "column_label": "Genome B",
            "value": 0.2,
        },
        {
            "row_label": "Genome A",
            "column_label": "Genome C",
            "value": 0.4,
        },
        {
            "row_label": "Genome B",
            "column_label": "Genome C",
            "value": 0.3,
        },
    ]

def test_upper_triangle_rows_match_pairs_and_vector_order():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        values=[
            [0.0, 0.2, 0.4],
            [0.2, 0.0, 0.3],
            [0.4, 0.3, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    rows = matrix.to_upper_triangle_rows()
    pairs = matrix.upper_triangle_pairs()
    vector = matrix.to_upper_triangle_vector()

    assert [
        (
            row["row_label"],
            row["column_label"],
            row["value"],
        )
        for row in rows
    ] == [
        (
            row_label,
            column_label,
            value,
        )
        for (
            row_label,
            column_label,
        ), value in zip(
            pairs,
            vector,
            strict=True,
        )
    ]

def test_upper_triangle_rows_returns_independent_dictionaries():
    matrix = GenomeMatrix(
        labels=[
            "Genome A",
            "Genome B",
        ],
        values=[
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        metric="euclidean",
        kmer_length=3,
    )

    rows = matrix.to_upper_triangle_rows()
    rows[0]["value"] = 999.0

    assert matrix.values[0][1] == 0.5

def test_genome_collection_generates_euclidean_matrix_trajectory():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.euclidean_matrix_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        k_values=[1, 2, 3],
    )

    assert list(trajectory) == [1, 2, 3]

    assert all(
        len(vector) == 3
        for vector in trajectory.values()
    )

def test_euclidean_matrix_trajectory_matches_matrix_vectors():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    labels = [
        "Genome A",
        "Genome B",
        "Genome C",
    ]

    trajectory = collection.euclidean_matrix_trajectory(
        labels=labels,
        k_values=[1, 2, 3],
    )

    matrices = collection.euclidean_distance_matrices(
        labels=labels,
        k_values=[1, 2, 3],
    )

    assert trajectory == {
        k: matrix.to_upper_triangle_vector()
        for k, matrix in matrices.items()
    }

def test_euclidean_matrix_trajectory_preserves_k_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.euclidean_matrix_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        k_values=[3, 1, 2],
    )

    assert list(trajectory) == [3, 1, 2]

def test_euclidean_matrix_trajectory_reuses_multi_k_validation():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.euclidean_matrix_trajectory(
            labels=["Genome A", "Genome B"],
            k_values=[1, 2, 2],
        )
    
def test_genome_collection_generates_cosine_matrix_trajectory():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.cosine_matrix_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        k_values=[1, 2, 3],
    )

    assert list(trajectory) == [1, 2, 3]

    assert all(
        len(vector) == 3
        for vector in trajectory.values()
    )

def test_cosine_matrix_trajectory_matches_matrix_vectors():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    labels = [
        "Genome A",
        "Genome B",
        "Genome C",
    ]

    trajectory = collection.cosine_matrix_trajectory(
        labels=labels,
        k_values=[1, 2, 3],
    )

    matrices = collection.cosine_similarity_matrices(
        labels=labels,
        k_values=[1, 2, 3],
    )

    assert trajectory == {
        k: matrix.to_upper_triangle_vector()
        for k, matrix in matrices.items()
    }

def test_cosine_matrix_trajectory_preserves_k_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.cosine_matrix_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        k_values=[3, 1, 2],
    )

    assert list(trajectory) == [3, 1, 2]

def test_cosine_matrix_trajectory_reuses_multi_k_validation():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.cosine_matrix_trajectory(
            labels=["Genome A", "Genome B"],
            k_values=[1, 2, 2],
        )

def test_genome_collection_generates_euclidean_pair_trajectory():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.euclidean_pair_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        row_label="Genome A",
        column_label="Genome B",
        k_values=[1, 2, 3],
    )

    assert list(trajectory) == [1, 2, 3]

    assert all(
        isinstance(value, float)
        for value in trajectory.values()
    )

def test_euclidean_pair_trajectory_matches_matrix_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    labels = [
        "Genome A",
        "Genome B",
        "Genome C",
    ]

    trajectory = collection.euclidean_pair_trajectory(
        labels=labels,
        row_label="Genome A",
        column_label="Genome C",
        k_values=[1, 2, 3],
    )

    matrices = collection.euclidean_distance_matrices(
        labels=labels,
        k_values=[1, 2, 3],
    )

    assert trajectory == {
        k: matrix.get_value(
            row_label="Genome A",
            column_label="Genome C",
        )
        for k, matrix in matrices.items()
    }

def test_euclidean_pair_trajectory_preserves_k_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    trajectory = collection.euclidean_pair_trajectory(
        labels=["Genome A", "Genome B"],
        row_label="Genome A",
        column_label="Genome B",
        k_values=[3, 1, 2],
    )

    assert list(trajectory) == [3, 1, 2]

def test_euclidean_pair_trajectory_rejects_unknown_label():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="Unknown genome matrix label: Unknown",
    ):
        collection.euclidean_pair_trajectory(
            labels=["Genome A", "Genome B"],
            row_label="Genome A",
            column_label="Unknown",
            k_values=[1, 2, 3],
        )

def test_euclidean_pair_trajectory_reuses_multi_k_validation():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.euclidean_pair_trajectory(
            labels=["Genome A", "Genome B"],
            row_label="Genome A",
            column_label="Genome B",
            k_values=[1, 2, 2],
        )

def test_genome_collection_generates_cosine_pair_trajectory():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    trajectory = collection.cosine_pair_trajectory(
        labels=[
            "Genome A",
            "Genome B",
            "Genome C",
        ],
        row_label="Genome A",
        column_label="Genome B",
        k_values=[1, 2, 3],
    )

    assert list(trajectory) == [1, 2, 3]

    assert all(
        isinstance(value, float)
        for value in trajectory.values()
    )

def test_cosine_pair_trajectory_matches_matrix_values():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    labels = [
        "Genome A",
        "Genome B",
        "Genome C",
    ]

    trajectory = collection.cosine_pair_trajectory(
        labels=labels,
        row_label="Genome A",
        column_label="Genome C",
        k_values=[1, 2, 3],
    )

    matrices = collection.cosine_similarity_matrices(
        labels=labels,
        k_values=[1, 2, 3],
    )

    assert trajectory == {
        k: matrix.get_value(
            row_label="Genome A",
            column_label="Genome C",
        )
        for k, matrix in matrices.items()
    }

def test_cosine_pair_trajectory_preserves_k_order():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    trajectory = collection.cosine_pair_trajectory(
        labels=["Genome A", "Genome B"],
        row_label="Genome A",
        column_label="Genome B",
        k_values=[3, 1, 2],
    )

    assert list(trajectory) == [3, 1, 2]


def test_cosine_pair_trajectory_rejects_unknown_label():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="Unknown genome matrix label: Unknown",
    ):
        collection.cosine_pair_trajectory(
            labels=["Genome A", "Genome B"],
            row_label="Genome A",
            column_label="Unknown",
            k_values=[1, 2, 3],
        )

def test_cosine_pair_trajectory_reuses_multi_k_validation():
    collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    with pytest.raises(
        ValueError,
        match="k-mer lengths must be unique",
    ):
        collection.cosine_pair_trajectory(
            labels=["Genome A", "Genome B"],
            row_label="Genome A",
            column_label="Genome B",
            k_values=[1, 2, 2],
        )