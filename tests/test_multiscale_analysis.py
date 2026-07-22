import math

import pytest

from src.genome import (
    Genome,
    GenomeCollection,
)


LABELS = [
    "Genome A",
    "Genome B",
    "Genome C",
]


@pytest.fixture
def collection() -> GenomeCollection:
    return GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )


def test_pair_trajectory_step_differences_are_signed():
    differences = (
        GenomeCollection
        .pair_trajectory_step_differences(
            {
                1: 0.10,
                2: 0.15,
                3: 0.12,
            }
        )
    )

    assert differences == {
        (1, 2): pytest.approx(0.05),
        (2, 3): pytest.approx(-0.03),
    }


def test_pair_trajectory_step_differences_preserve_order():
    differences = (
        GenomeCollection
        .pair_trajectory_step_differences(
            {
                3: 0.30,
                1: 0.10,
                2: 0.20,
            }
        )
    )

    assert list(differences) == [
        (3, 1),
        (1, 2),
    ]


def test_pair_trajectory_step_differences_reject_empty():
    with pytest.raises(
        ValueError,
        match=r"Pair trajectory cannot be empty\.",
    ):
        GenomeCollection.pair_trajectory_step_differences({})


def test_pair_trajectory_step_differences_require_two_scales():
    with pytest.raises(
        ValueError,
        match=(
            r"Pair trajectory must contain "
            r"at least two k-mer scales\."
        ),
    ):
        GenomeCollection.pair_trajectory_step_differences(
            {1: 0.10}
        )


def test_matrix_trajectory_step_distances():
    distances = (
        GenomeCollection.matrix_trajectory_step_distances(
            {
                1: [0.0, 0.0],
                2: [3.0, 4.0],
                3: [3.0, 4.0],
            }
        )
    )

    assert distances == {
        (1, 2): pytest.approx(5.0),
        (2, 3): pytest.approx(0.0),
    }


def test_matrix_trajectory_step_distances_match_euclidean_formula():
    trajectory = {
        1: [0.10, 0.20, 0.30],
        2: [0.15, 0.18, 0.35],
    }

    expected = math.sqrt(
        (0.15 - 0.10) ** 2
        + (0.18 - 0.20) ** 2
        + (0.35 - 0.30) ** 2
    )

    distances = (
        GenomeCollection
        .matrix_trajectory_step_distances(
            trajectory
        )
    )

    assert distances[(1, 2)] == pytest.approx(
        expected
    )


def test_matrix_trajectory_step_distances_reject_unequal_vectors():
    with pytest.raises(
        ValueError,
        match=(
            r"Matrix trajectory vectors must have "
            r"the same length\."
        ),
    ):
        GenomeCollection.matrix_trajectory_step_distances(
            {
                1: [0.1, 0.2],
                2: [0.1],
            }
        )


def test_matrix_pair_contributions_preserve_labels_and_sign():
    contributions = (
        GenomeCollection.matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory={
                1: [0.10, 0.20, 0.30],
                2: [0.15, 0.18, 0.40],
            },
        )
    )

    rows = contributions[(1, 2)]

    assert rows[0] == {
        "row_label": "Genome B",
        "column_label": "Genome C",
        "difference": pytest.approx(0.10),
        "absolute_difference": pytest.approx(0.10),
    }

    negative_row = next(
        row
        for row in rows
        if row["row_label"] == "Genome A"
        and row["column_label"] == "Genome C"
    )

    assert negative_row["difference"] == (
        pytest.approx(-0.02)
    )


def test_matrix_pair_contributions_are_sorted_by_absolute_change():
    contributions = (
        GenomeCollection.matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory={
                1: [0.10, 0.20, 0.30],
                2: [0.11, 0.25, 0.28],
            },
        )
    )

    absolute_differences = [
        row["absolute_difference"]
        for row in contributions[(1, 2)]
    ]

    assert absolute_differences == sorted(
        absolute_differences,
        reverse=True,
    )


def test_matrix_pair_contributions_reject_vector_label_mismatch():
    with pytest.raises(
        ValueError,
        match=(
            r"Matrix trajectory vector length must "
            r"match the number of unique label pairs\."
        ),
    ):
        GenomeCollection.matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory={
                1: [0.10, 0.20],
                2: [0.15, 0.25],
            },
        )


def test_euclidean_pair_step_wrapper_matches_generic_method(
    collection: GenomeCollection,
):
    trajectory = collection.euclidean_pair_trajectory(
        labels=LABELS,
        row_label="Genome A",
        column_label="Genome B",
        k_values=[1, 2, 3],
    )

    assert (
        collection
        .euclidean_pair_trajectory_step_differences(
            labels=LABELS,
            row_label="Genome A",
            column_label="Genome B",
            k_values=[1, 2, 3],
        )
        == collection
        .pair_trajectory_step_differences(
            trajectory
        )
    )


def test_cosine_pair_step_wrapper_matches_generic_method(
    collection: GenomeCollection,
):
    trajectory = collection.cosine_pair_trajectory(
        labels=LABELS,
        row_label="Genome A",
        column_label="Genome B",
        k_values=[1, 2, 3],
    )

    assert (
        collection
        .cosine_pair_trajectory_step_differences(
            labels=LABELS,
            row_label="Genome A",
            column_label="Genome B",
            k_values=[1, 2, 3],
        )
        == collection
        .pair_trajectory_step_differences(
            trajectory
        )
    )


def test_euclidean_matrix_step_wrapper_matches_generic_method(
    collection: GenomeCollection,
):
    trajectory = collection.euclidean_matrix_trajectory(
        labels=LABELS,
        k_values=[1, 2, 3],
    )

    assert (
        collection
        .euclidean_matrix_trajectory_step_distances(
            labels=LABELS,
            k_values=[1, 2, 3],
        )
        == collection
        .matrix_trajectory_step_distances(
            trajectory
        )
    )


def test_cosine_matrix_contribution_wrapper_matches_generic_method(
    collection: GenomeCollection,
):
    trajectory = collection.cosine_matrix_trajectory(
        labels=LABELS,
        k_values=[1, 2, 3],
    )

    assert (
        collection.cosine_matrix_pair_contributions(
            labels=LABELS,
            k_values=[1, 2, 3],
        )
        == collection
        .matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory=trajectory,
        )
    )
