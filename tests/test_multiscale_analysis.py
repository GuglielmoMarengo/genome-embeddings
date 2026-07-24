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


def test_deformation_partition_is_additive():
    contributions = (
        GenomeCollection.matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory={
                1: [0.10, 0.20, 0.30],
                2: [0.15, 0.18, 0.40],
            },
        )
    )

    partition = (
        GenomeCollection
        .matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label="Genome C",
        )[(1, 2)]
    )

    assert partition[
        "total_squared_deformation"
    ] == pytest.approx(
        partition["selected_squared_deformation"]
        + partition["remaining_squared_deformation"]
    )

    assert partition["selected_share"] + partition[
        "remaining_share"
    ] == pytest.approx(1.0)


def test_deformation_partition_reports_expected_values():
    contributions = {
        (1, 2): [
            {
                "row_label": "Genome A",
                "column_label": "Genome B",
                "difference": 3.0,
                "absolute_difference": 3.0,
            },
            {
                "row_label": "Genome A",
                "column_label": "Genome C",
                "difference": 4.0,
                "absolute_difference": 4.0,
            },
            {
                "row_label": "Genome B",
                "column_label": "Genome C",
                "difference": 12.0,
                "absolute_difference": 12.0,
            },
        ]
    }

    partition = (
        GenomeCollection
        .matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label="Genome C",
        )[(1, 2)]
    )

    assert partition[
        "total_squared_deformation"
    ] == pytest.approx(169.0)
    assert partition[
        "selected_squared_deformation"
    ] == pytest.approx(160.0)
    assert partition[
        "remaining_squared_deformation"
    ] == pytest.approx(9.0)
    assert partition["total_distance"] == pytest.approx(13.0)
    assert partition["selected_distance"] == pytest.approx(
        math.sqrt(160.0)
    )
    assert partition["remaining_distance"] == pytest.approx(3.0)
    assert partition["selected_pair_count"] == 2
    assert partition["remaining_pair_count"] == 1


def test_deformation_partition_preserves_transition_order():
    contributions = {
        (3, 1): [
            {
                "row_label": "Genome A",
                "column_label": "Genome B",
                "difference": 1.0,
                "absolute_difference": 1.0,
            }
        ],
        (1, 2): [
            {
                "row_label": "Genome A",
                "column_label": "Genome B",
                "difference": 2.0,
                "absolute_difference": 2.0,
            }
        ],
    }

    partitions = (
        GenomeCollection
        .matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label="Genome A",
        )
    )

    assert list(partitions) == [
        (3, 1),
        (1, 2),
    ]


def test_deformation_partition_rejects_empty_contributions():
    with pytest.raises(
        ValueError,
        match=(
            r"Matrix trajectory contributions "
            r"cannot be empty\."
        ),
    ):
        GenomeCollection.matrix_trajectory_deformation_partition(
            contributions={},
            selected_label="Genome A",
        )


def test_deformation_partition_rejects_empty_rows():
    with pytest.raises(
        ValueError,
        match=(
            r"Matrix trajectory contribution "
            r"rows cannot be empty\."
        ),
    ):
        GenomeCollection.matrix_trajectory_deformation_partition(
            contributions={(1, 2): []},
            selected_label="Genome A",
        )


def test_deformation_partition_rejects_unknown_label():
    contributions = (
        GenomeCollection.matrix_trajectory_pair_contributions(
            labels=LABELS,
            trajectory={
                1: [0.10, 0.20, 0.30],
                2: [0.15, 0.18, 0.40],
            },
        )
    )

    with pytest.raises(
        ValueError,
        match=(
            r"Selected label is not present in "
            r"the contribution rows\."
        ),
    ):
        GenomeCollection.matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label="Unknown",
        )


def test_euclidean_deformation_partition_wrapper_matches_generic(
    collection: GenomeCollection,
):
    contributions = collection.euclidean_matrix_pair_contributions(
        labels=LABELS,
        k_values=[1, 2, 3],
    )

    assert collection.euclidean_matrix_deformation_partition(
        labels=LABELS,
        k_values=[1, 2, 3],
        selected_label="Genome C",
    ) == collection.matrix_trajectory_deformation_partition(
        contributions=contributions,
        selected_label="Genome C",
    )


def test_remaining_deformation_matches_subset_step_distance():
    full_collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
            Genome("GGGGTTTT"),
        ]
    )

    subset_collection = GenomeCollection(
        [
            Genome("ACGTACGT"),
            Genome("AAAACCCC"),
        ]
    )

    full_partition = (
        full_collection
        .euclidean_matrix_deformation_partition(
            labels=LABELS,
            k_values=[1, 2, 3],
            selected_label="Genome C",
        )
    )

    subset_distances = (
        subset_collection
        .euclidean_matrix_trajectory_step_distances(
            labels=["Genome A", "Genome B"],
            k_values=[1, 2, 3],
        )
    )

    for transition, subset_distance in subset_distances.items():
        assert full_partition[transition][
            "remaining_squared_deformation"
        ] == pytest.approx(subset_distance ** 2)


def test_cosine_deformation_partition_wrapper_matches_generic(
    collection: GenomeCollection,
):
    contributions = collection.cosine_matrix_pair_contributions(
        labels=LABELS,
        k_values=[1, 2, 3],
    )

    assert collection.cosine_matrix_deformation_partition(
        labels=LABELS,
        k_values=[1, 2, 3],
        selected_label="Genome C",
    ) == collection.matrix_trajectory_deformation_partition(
        contributions=contributions,
        selected_label="Genome C",
    )

def test_matrix_ranking_trajectory_preserves_k_order(
    collection: GenomeCollection,
):
    matrices = collection.euclidean_distance_matrices(
        labels=LABELS,
        k_values=[3, 1, 2],
    )

    trajectory = (
        GenomeCollection.matrix_ranking_trajectory(
            matrices=matrices,
            reference_label="Genome A",
        )
    )

    assert list(trajectory) == [3, 1, 2]

    assert trajectory[3] == matrices[3].rank_by_label(
        label="Genome A"
    )


def test_ranking_stability_reports_perfect_match():
    stability = (
        GenomeCollection.ranking_trajectory_stability(
            {
                1: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                    ("Genome D", 0.3),
                ],
                2: [
                    ("Genome B", 0.4),
                    ("Genome C", 0.5),
                    ("Genome D", 0.6),
                ],
            }
        )[(1, 2)]
    )

    assert stability["kendall_tau"] == pytest.approx(1.0)
    assert stability["concordant_pairs"] == 3
    assert stability["discordant_pairs"] == 0
    assert stability["mean_absolute_rank_shift"] == pytest.approx(0.0)
    assert stability["max_rank_shift"] == 0
    assert stability["exact_match"] is True


def test_ranking_stability_reports_complete_reversal():
    stability = (
        GenomeCollection.ranking_trajectory_stability(
            {
                1: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                    ("Genome D", 0.3),
                ],
                2: [
                    ("Genome D", 0.1),
                    ("Genome C", 0.2),
                    ("Genome B", 0.3),
                ],
            }
        )[(1, 2)]
    )

    assert stability["kendall_tau"] == pytest.approx(-1.0)
    assert stability["concordant_pairs"] == 0
    assert stability["discordant_pairs"] == 3
    assert stability["max_rank_shift"] == 2
    assert stability["exact_match"] is False


def test_ranking_stability_reports_single_inversion():
    stability = (
        GenomeCollection.ranking_trajectory_stability(
            {
                1: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                    ("Genome D", 0.3),
                ],
                2: [
                    ("Genome C", 0.1),
                    ("Genome B", 0.2),
                    ("Genome D", 0.3),
                ],
            }
        )[(1, 2)]
    )

    assert stability["kendall_tau"] == pytest.approx(1 / 3)
    assert stability["concordant_pairs"] == 2
    assert stability["discordant_pairs"] == 1
    assert stability["mean_absolute_rank_shift"] == pytest.approx(
        2 / 3
    )


def test_ranking_stability_preserves_transition_order():
    stability = (
        GenomeCollection.ranking_trajectory_stability(
            {
                3: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                ],
                1: [
                    ("Genome C", 0.1),
                    ("Genome B", 0.2),
                ],
                2: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                ],
            }
        )
    )

    assert list(stability) == [
        (3, 1),
        (1, 2),
    ]


def test_ranking_stability_rejects_empty_rankings():
    with pytest.raises(
        ValueError,
        match=r"Ranking trajectory cannot be empty\.",
    ):
        GenomeCollection.ranking_trajectory_stability({})


def test_ranking_stability_requires_two_scales():
    with pytest.raises(
        ValueError,
        match=(
            r"Ranking trajectory must contain "
            r"at least two k-mer scales\."
        ),
    ):
        GenomeCollection.ranking_trajectory_stability(
            {
                1: [
                    ("Genome B", 0.1),
                ]
            }
        )


def test_ranking_stability_rejects_mismatched_labels():
    with pytest.raises(
        ValueError,
        match=(
            r"Ranking trajectory entries must "
            r"contain the same labels\."
        ),
    ):
        GenomeCollection.ranking_trajectory_stability(
            {
                1: [
                    ("Genome B", 0.1),
                    ("Genome C", 0.2),
                ],
                2: [
                    ("Genome B", 0.1),
                    ("Genome D", 0.2),
                ],
            }
        )


def test_euclidean_ranking_stability_wrapper_matches_generic(
    collection: GenomeCollection,
):
    rankings = collection.euclidean_ranking_trajectory(
        labels=LABELS,
        reference_label="Genome A",
        k_values=[1, 2, 3],
    )

    assert collection.euclidean_ranking_stability(
        labels=LABELS,
        reference_label="Genome A",
        k_values=[1, 2, 3],
    ) == collection.ranking_trajectory_stability(
        rankings
    )


def test_cosine_ranking_stability_wrapper_matches_generic(
    collection: GenomeCollection,
):
    rankings = collection.cosine_ranking_trajectory(
        labels=LABELS,
        reference_label="Genome A",
        k_values=[1, 2, 3],
    )

    assert collection.cosine_ranking_stability(
        labels=LABELS,
        reference_label="Genome A",
        k_values=[1, 2, 3],
    ) == collection.ranking_trajectory_stability(
        rankings
    )
