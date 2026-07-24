from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from matplotlib.figure import Figure

from src.genome import GenomeMatrix
from src.visualization import (
    plot_matrix_distribution,
    plot_matrix_heatmap,
    plot_pair_trajectory,
    plot_trajectory_distributions,
    save_figure,
)


@pytest.fixture
def euclidean_matrix() -> GenomeMatrix:
    return GenomeMatrix(
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


@pytest.fixture
def cosine_matrix() -> GenomeMatrix:
    return GenomeMatrix(
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


def test_plot_matrix_heatmap_returns_figure(
    euclidean_matrix: GenomeMatrix,
):
    figure = plot_matrix_heatmap(
        euclidean_matrix
    )

    assert isinstance(
        figure,
        Figure,
    )


def test_plot_matrix_heatmap_uses_metric_title(
    cosine_matrix: GenomeMatrix,
):
    figure = plot_matrix_heatmap(
        cosine_matrix,
        annotate=False,
    )

    axis = figure.axes[0]

    assert axis.get_title() == (
        "Cosine Similarity Matrix (k=3)"
    )


def test_plot_matrix_heatmap_adds_annotations(
    euclidean_matrix: GenomeMatrix,
):
    figure = plot_matrix_heatmap(
        euclidean_matrix,
        annotate=True,
    )

    axis = figure.axes[0]

    assert len(axis.texts) == 4


def test_plot_pair_trajectory_uses_trajectory_values():
    trajectory = {
        1: 0.10,
        2: 0.12,
        3: 0.15,
    }

    figure = plot_pair_trajectory(
        trajectory=trajectory,
        row_label="Genome A",
        column_label="Genome B",
        metric="euclidean",
    )

    axis = figure.axes[0]
    line = axis.lines[0]

    assert list(
        line.get_xdata()
    ) == [
        1,
        2,
        3,
    ]

    assert list(
        line.get_ydata()
    ) == pytest.approx(
        [
            0.10,
            0.12,
            0.15,
        ]
    )


def test_plot_pair_trajectory_rejects_empty_trajectory():
    with pytest.raises(
        ValueError,
        match=(
            r"Trajectory cannot be empty\."
        ),
    ):
        plot_pair_trajectory(
            trajectory={},
            row_label="Genome A",
            column_label="Genome B",
            metric="euclidean",
        )


def test_plot_matrix_distribution_uses_unique_values(
    euclidean_matrix: GenomeMatrix,
):
    figure = plot_matrix_distribution(
        euclidean_matrix
    )

    axis = figure.axes[0]

    assert axis.get_title() == (
        "Euclidean Distance Distribution "
        "(k=3)"
    )

    assert axis.get_xlabel() == (
        "Euclidean Distance"
    )


def test_plot_trajectory_distributions_uses_k_labels():
    trajectory = {
        1: [
            0.10,
            0.20,
            0.30,
        ],
        2: [
            0.12,
            0.22,
            0.32,
        ],
        3: [
            0.15,
            0.25,
            0.35,
        ],
    }

    figure = (
        plot_trajectory_distributions(
            trajectory=trajectory,
            metric="euclidean",
        )
    )

    axis = figure.axes[0]

    tick_labels = [
        label.get_text()
        for label in axis.get_xticklabels()
    ]

    assert tick_labels == [
        "1",
        "2",
        "3",
    ]


def test_plot_trajectory_distributions_rejects_empty_trajectory():
    with pytest.raises(
        ValueError,
        match=(
            r"Trajectory cannot be empty\."
        ),
    ):
        plot_trajectory_distributions(
            trajectory={},
            metric="cosine",
        )


def test_save_figure_creates_output_file(
    tmp_path: Path,
    euclidean_matrix: GenomeMatrix,
):
    figure = plot_matrix_heatmap(
        euclidean_matrix
    )

    output_path = (
        tmp_path
        / "visualizations"
        / "heatmap.png"
    )

    saved_path = save_figure(
        figure=figure,
        output_path=output_path,
    )

    assert saved_path.exists()
    assert saved_path.is_file()


def test_save_figure_returns_requested_path(
    tmp_path: Path,
    euclidean_matrix: GenomeMatrix,
):
    figure = plot_matrix_distribution(
        euclidean_matrix
    )

    output_path = (
        tmp_path
        / "distribution.png"
    )

    saved_path = save_figure(
        figure=figure,
        output_path=output_path,
    )

    assert saved_path == output_path

def test_plot_jensen_shannon_heatmap_uses_metric_title():
    matrix = GenomeMatrix(
        labels=["Genome A", "Genome B"],
        values=[[0.0, 0.4], [0.4, 0.0]],
        metric="jensen_shannon",
        kmer_length=3,
    )

    figure = plot_matrix_heatmap(matrix, annotate=False)

    assert figure.axes[0].get_title() == (
        "Jensen-Shannon Distance Matrix (k=3)"
    )


def test_plot_jensen_shannon_pair_trajectory_uses_distance_label():
    figure = plot_pair_trajectory(
        trajectory={1: 0.1, 2: 0.2},
        row_label="Genome A",
        column_label="Genome B",
        metric="jensen_shannon",
    )

    assert figure.axes[0].get_ylabel() == "Jensen-Shannon Distance"
