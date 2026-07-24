from pathlib import Path
from typing import Literal

from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from src.genome import GenomeMatrix


MetricName = Literal["euclidean", "cosine"]


def _metric_title(
    metric: str,
) -> str:
    titles = {
        "euclidean": "Euclidean Distance",
        "cosine": "Cosine Similarity",
    }

    try:
        return titles[metric]
    except KeyError as error:
        raise ValueError(
            f"Unsupported visualization metric: {metric}."
        ) from error


def _validate_trajectory(
    trajectory: dict[int, float] | dict[int, list[float]],
) -> None:
    if not trajectory:
        raise ValueError(
            "Trajectory cannot be empty."
        )


def plot_matrix_heatmap(
    matrix: GenomeMatrix,
    annotate: bool = True,
    decimal_places: int = 3,
) -> Figure:
    if decimal_places < 0:
        raise ValueError(
            "decimal_places cannot be negative."
        )

    metric_title = _metric_title(matrix.metric)

    matrix_size = len(matrix.labels)

    figure_width = max(
        7.0,
        matrix_size * 1.35,
    )

    figure_height = max(
        6.0,
        matrix_size * 1.15,
    )

    figure, axis = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        )
    )

    image = axis.imshow(
        matrix.values,
        aspect="auto",
    )

    colorbar = figure.colorbar(
        image,
        ax=axis,
    )

    colorbar.set_label(metric_title)

    tick_positions = list(
        range(matrix_size)
    )

    axis.set_xticks(
        tick_positions,
        labels=matrix.labels,
        rotation=45,
        ha="right",
    )

    axis.set_yticks(
        tick_positions,
        labels=matrix.labels,
    )

    axis.set_xlabel("Genome")
    axis.set_ylabel("Genome")

    axis.set_title(
        f"{metric_title} Matrix "
        f"(k={matrix.kmer_length})"
    )

    if annotate:
        number_format = (
            f".{decimal_places}f"
        )

        for row_index, row in enumerate(
            matrix.values
        ):
            for column_index, value in enumerate(
                row
            ):
                axis.text(
                    column_index,
                    row_index,
                    format(
                        value,
                        number_format,
                    ),
                    ha="center",
                    va="center",
                )

    figure.tight_layout()

    return figure


def plot_pair_trajectory(
    trajectory: dict[int, float],
    row_label: str,
    column_label: str,
    metric: MetricName,
) -> Figure:
    _validate_trajectory(trajectory)

    metric_title = _metric_title(metric)

    k_values = list(
        trajectory.keys()
    )

    values = list(
        trajectory.values()
    )

    figure, axis = plt.subplots(
        figsize=(8.0, 5.0)
    )

    axis.plot(
        k_values,
        values,
        marker="o",
    )

    axis.set_xticks(k_values)
    axis.set_xlabel("k-mer length")
    axis.set_ylabel(metric_title)

    axis.set_title(
        f"{row_label} vs {column_label}\n"
        f"{metric_title} Across k-mer Scales"
    )

    axis.grid(
        visible=True,
        alpha=0.3,
    )

    figure.tight_layout()

    return figure


def plot_matrix_distribution(
    matrix: GenomeMatrix,
    bins: int | None = None,
) -> Figure:
    values = (
        matrix.to_upper_triangle_vector()
    )

    if not values:
        raise ValueError(
            "Matrix distribution requires "
            "at least two genome labels."
        )

    if bins is None:
        bins = min(
            10,
            max(
                1,
                len(values),
            ),
        )

    if bins <= 0:
        raise ValueError(
            "bins must be positive."
        )

    metric_title = _metric_title(
        matrix.metric
    )

    figure, axis = plt.subplots(
        figsize=(8.0, 5.0)
    )

    axis.hist(
        values,
        bins=bins,
        edgecolor="black",
    )

    axis.set_xlabel(metric_title)
    axis.set_ylabel(
        "Number of genome pairs"
    )

    axis.set_title(
        f"{metric_title} Distribution "
        f"(k={matrix.kmer_length})"
    )

    axis.grid(
        visible=True,
        axis="y",
        alpha=0.3,
    )

    figure.tight_layout()

    return figure


def plot_trajectory_distributions(
    trajectory: dict[int, list[float]],
    metric: MetricName,
) -> Figure:
    _validate_trajectory(trajectory)

    for k, values in trajectory.items():
        if not values:
            raise ValueError(
                "Trajectory vectors cannot "
                f"be empty. Got empty vector for k={k}."
            )

    metric_title = _metric_title(metric)

    k_values = list(
        trajectory.keys()
    )

    distributions = list(
        trajectory.values()
    )

    figure, axis = plt.subplots(
        figsize=(9.0, 5.5)
    )

    axis.boxplot(
        distributions,
        tick_labels=[
            str(k)
            for k in k_values
        ],
    )

    axis.set_xlabel("k-mer length")
    axis.set_ylabel(metric_title)

    axis.set_title(
        f"{metric_title} Distribution "
        "Across k-mer Scales"
    )

    axis.grid(
        visible=True,
        axis="y",
        alpha=0.3,
    )

    figure.tight_layout()

    return figure



def plot_ranking_stability(
    stability: dict[
        tuple[int, int],
        dict[str, float | int | bool],
    ],
    metric: MetricName,
) -> Figure:
    if not stability:
        raise ValueError(
            "Ranking stability cannot be empty."
        )

    metric_title = _metric_title(metric)

    transitions = [
        f"{first_k}->{second_k}"
        for first_k, second_k in stability
    ]

    tau_values = [
        float(values["kendall_tau"])
        for values in stability.values()
    ]

    figure, axis = plt.subplots(
        figsize=(8.5, 5.0)
    )

    axis.plot(
        transitions,
        tau_values,
        marker="o",
    )

    axis.axhline(
        1.0,
        linewidth=1.0,
        linestyle="--",
        alpha=0.4,
    )

    axis.axhline(
        0.0,
        linewidth=1.0,
        alpha=0.4,
    )

    axis.set_ylim(-1.05, 1.05)
    axis.set_xlabel("k-mer scale transition")
    axis.set_ylabel("Kendall rank correlation")

    axis.set_title(
        f"{metric_title} Ranking Stability "
        "Across k-mer Scales"
    )

    axis.grid(
        visible=True,
        axis="y",
        alpha=0.3,
    )

    figure.tight_layout()

    return figure


def save_figure(
    figure: Figure,
    output_path: str | Path,
    dpi: int = 150,
    close: bool = False,
) -> Path:
    if dpi <= 0:
        raise ValueError(
            "dpi must be positive."
        )

    path = Path(output_path)

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if close:
        plt.close(figure)

    return path