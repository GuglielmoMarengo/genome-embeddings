"""Application-facing dataset and analysis services.

The module contains no NiceGUI imports, which keeps the scientific workflow
unit-testable and reusable from the command-line report and future clients.
"""

from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path

from src.descriptor_v2 import DescriptorV2Collection, GenomeDescriptorV2
from src.genome import Genome, GenomeCollection, GenomeMatrix
from src.kmer_distribution import KmerDistributionCollection


DEMO_LABELS = [
    "Aequorea GFP",
    "Acropora GFP",
    "Discosoma FP583",
    "S. aureus catA",
    "S. cerevisiae TPI1",
    "Periodic control",
]

DEMO_RELATIVE_PATHS = [
    "data/fluorescent_proteins/aequorea_victoria_gfp_cds.fasta",
    "data/fluorescent_proteins/acropora_millepora_gfp_cds.fasta",
    "data/fluorescent_proteins/discosoma_fp583_cds.fasta",
    "data/controls/biological/staphylococcus_aureus_cata_cds.fasta",
    "data/controls/biological/saccharomyces_cerevisiae_tpi1_cds.fasta",
    "data/controls/periodic_sequence.fasta",
]


@dataclass(frozen=True, slots=True)
class DatasetRecord:
    label: str
    genome: Genome
    source: str

    def to_row(self) -> dict[str, object]:
        return {
            "label": self.label,
            "source": self.source,
            "length": self.genome.length(),
            "gc_content": self.genome.gc_content(),
            "header": self.genome.header or "",
        }


@dataclass(frozen=True, slots=True)
class DashboardConfig:
    k_values: tuple[int, ...]
    selected_k: int
    reference_label: str
    comparison_label: str

    def __post_init__(self) -> None:
        if not self.k_values:
            raise ValueError("k-mer lengths cannot be empty.")
        if len(self.k_values) != len(set(self.k_values)):
            raise ValueError("k-mer lengths must be unique.")
        if self.selected_k not in self.k_values:
            raise ValueError("selected_k must be included in k_values.")
        if self.reference_label == self.comparison_label:
            raise ValueError("Reference and comparison labels must differ.")


@dataclass(slots=True)
class DashboardAnalysis:
    config: DashboardConfig
    labels: list[str]
    dataset_rows: list[dict[str, object]]
    descriptor_v2_rows: list[dict[str, object]]
    legacy_euclidean_matrix: GenomeMatrix
    legacy_cosine_matrix: GenomeMatrix
    descriptor_v2_euclidean_matrix: GenomeMatrix
    descriptor_v2_cosine_matrix: GenomeMatrix
    embedding_v2_euclidean_matrix: GenomeMatrix
    jensen_shannon_matrix: GenomeMatrix
    legacy_euclidean_pair_trajectory: dict[int, float]
    descriptor_v2_euclidean_pair_trajectory: dict[int, float]
    jensen_shannon_pair_trajectory: dict[int, float]
    jensen_shannon_step_distances: dict[tuple[int, int], float]
    jensen_shannon_ranking_trajectory: dict[int, list[tuple[str, float]]]
    jensen_shannon_ranking_stability: dict[
        tuple[int, int],
        dict[str, float | int | bool],
    ]

    def summary(self) -> dict[str, object]:
        reference = self.config.reference_label
        comparison = self.config.comparison_label
        return {
            "dataset_size": len(self.labels),
            "selected_k": self.config.selected_k,
            "k_values": list(self.config.k_values),
            "reference_label": reference,
            "comparison_label": comparison,
            "legacy_euclidean_distance": (
                self.legacy_euclidean_matrix.get_value(reference, comparison)
            ),
            "legacy_cosine_similarity": (
                self.legacy_cosine_matrix.get_value(reference, comparison)
            ),
            "descriptor_v2_euclidean_distance": (
                self.descriptor_v2_euclidean_matrix.get_value(
                    reference,
                    comparison,
                )
            ),
            "descriptor_v2_cosine_similarity": (
                self.descriptor_v2_cosine_matrix.get_value(
                    reference,
                    comparison,
                )
            ),
            "embedding_v2_euclidean_distance": (
                self.embedding_v2_euclidean_matrix.get_value(
                    reference,
                    comparison,
                )
            ),
            "jensen_shannon_distance": (
                self.jensen_shannon_matrix.get_value(reference, comparison)
            ),
        }

    def to_dict(self) -> dict[str, object]:
        return {
            "summary": self.summary(),
            "dataset": self.dataset_rows,
            "descriptor_v2": self.descriptor_v2_rows,
            "matrices": {
                "legacy_euclidean": self.legacy_euclidean_matrix.to_dict(),
                "legacy_cosine": self.legacy_cosine_matrix.to_dict(),
                "descriptor_v2_euclidean": (
                    self.descriptor_v2_euclidean_matrix.to_dict()
                ),
                "descriptor_v2_cosine": (
                    self.descriptor_v2_cosine_matrix.to_dict()
                ),
                "embedding_v2_euclidean": (
                    self.embedding_v2_euclidean_matrix.to_dict()
                ),
                "jensen_shannon": self.jensen_shannon_matrix.to_dict(),
            },
            "pair_trajectories": {
                "legacy_euclidean": self.legacy_euclidean_pair_trajectory,
                "descriptor_v2_euclidean": (
                    self.descriptor_v2_euclidean_pair_trajectory
                ),
                "jensen_shannon": self.jensen_shannon_pair_trajectory,
            },
            "jensen_shannon_step_distances": {
                f"{first_k}->{second_k}": value
                for (first_k, second_k), value in (
                    self.jensen_shannon_step_distances.items()
                )
            },
            "jensen_shannon_rankings": {
                str(k): ranking
                for k, ranking in self.jensen_shannon_ranking_trajectory.items()
            },
            "jensen_shannon_ranking_stability": {
                f"{first_k}->{second_k}": row
                for (first_k, second_k), row in (
                    self.jensen_shannon_ranking_stability.items()
                )
            },
        }

    def to_json(self, indent: int | None = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def summary_csv(self) -> str:
        output = io.StringIO()
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for name, value in self.summary().items():
            if isinstance(value, list):
                value = ",".join(str(item) for item in value)
            writer.writerow([name, value])
        return output.getvalue()


def load_demo_records(project_root: str | Path) -> list[DatasetRecord]:
    root = Path(project_root)
    return [
        DatasetRecord(
            label=label,
            genome=Genome.from_fasta(root / relative_path),
            source=relative_path,
        )
        for label, relative_path in zip(
            DEMO_LABELS,
            DEMO_RELATIVE_PATHS,
            strict=True,
        )
    ]


def parse_fasta_bytes(
    filename: str,
    data: bytes,
    label: str | None = None,
) -> DatasetRecord:
    safe_name = Path(filename).name
    if not safe_name:
        raise ValueError("Uploaded filename cannot be empty.")

    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("FASTA files must be UTF-8 encoded.") from error

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("FASTA file is empty.")
    if not lines[0].startswith(">"):
        raise ValueError("FASTA header must start with '>'.")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Dashboard uploads currently support one FASTA record.")

    genome = Genome(sequence="".join(lines[1:]), header=lines[0])
    resolved_label = label or Path(safe_name).stem.replace("_", " ")
    Genome.validate_string("Dataset label", resolved_label)

    return DatasetRecord(
        label=resolved_label,
        genome=genome,
        source=safe_name,
    )


def analyze_records(
    records: list[DatasetRecord],
    config: DashboardConfig,
) -> DashboardAnalysis:
    if not records:
        raise ValueError("Dataset cannot be empty.")

    labels = [record.label for record in records]
    if len(labels) != len(set(labels)):
        raise ValueError("Dataset labels must be unique.")
    if config.reference_label not in labels:
        raise ValueError("Reference label is not present in the dataset.")
    if config.comparison_label not in labels:
        raise ValueError("Comparison label is not present in the dataset.")

    genomes = [record.genome for record in records]
    legacy = GenomeCollection(genomes)
    descriptor_v2 = DescriptorV2Collection(genomes)
    distribution = KmerDistributionCollection(genomes)
    k_values = list(config.k_values)

    descriptor_rows = []
    for record in records:
        descriptor = GenomeDescriptorV2.from_genome(
            record.genome,
            config.selected_k,
        )
        descriptor_rows.append(
            {
                "label": record.label,
                "length": descriptor.length,
                "gc_content": descriptor.gc_content,
                "conditional_entropy": (
                    descriptor.conditional_nucleotide_entropy
                ),
                "finite_sample_entropy": (
                    descriptor.kmer.finite_sample_normalized_kmer_entropy
                ),
                "effective_kmer_count": (
                    descriptor.kmer.effective_kmer_count
                ),
                "theoretical_coverage": (
                    descriptor.kmer.theoretical_space_coverage
                ),
                "observable_coverage": (
                    descriptor.kmer.observable_space_coverage
                ),
                "singleton_fraction": descriptor.kmer.singleton_fraction,
                "repeated_window_fraction": (
                    descriptor.kmer.repeated_window_fraction
                ),
            }
        )

    v2_matrices = descriptor_v2.euclidean_distance_matrices(labels, k_values)
    legacy_trajectory = legacy.euclidean_pair_trajectory(
        labels=labels,
        row_label=config.reference_label,
        column_label=config.comparison_label,
        k_values=k_values,
    )
    v2_pair_trajectory = {
        k: matrix.get_value(config.reference_label, config.comparison_label)
        for k, matrix in v2_matrices.items()
    }

    return DashboardAnalysis(
        config=config,
        labels=labels,
        dataset_rows=[record.to_row() for record in records],
        descriptor_v2_rows=descriptor_rows,
        legacy_euclidean_matrix=legacy.euclidean_distance_matrix(
            labels=labels,
            k=config.selected_k,
        ),
        legacy_cosine_matrix=legacy.cosine_similarity_matrix(
            labels=labels,
            k=config.selected_k,
        ),
        descriptor_v2_euclidean_matrix=v2_matrices[config.selected_k],
        descriptor_v2_cosine_matrix=(
            descriptor_v2.cosine_similarity_matrix(
                labels=labels,
                k=config.selected_k,
            )
        ),
        embedding_v2_euclidean_matrix=(
            descriptor_v2.multiscale_embedding_distance_matrix(
                labels=labels,
                k_values=k_values,
            )
        ),
        jensen_shannon_matrix=distribution.distance_matrix(
            labels=labels,
            k=config.selected_k,
        ),
        legacy_euclidean_pair_trajectory=legacy_trajectory,
        descriptor_v2_euclidean_pair_trajectory=v2_pair_trajectory,
        jensen_shannon_pair_trajectory=distribution.pair_trajectory(
            labels=labels,
            row_label=config.reference_label,
            column_label=config.comparison_label,
            k_values=k_values,
        ),
        jensen_shannon_step_distances=(
            distribution.matrix_trajectory_step_distances(labels, k_values)
        ),
        jensen_shannon_ranking_trajectory=(
            distribution.ranking_trajectory(
                labels=labels,
                reference_label=config.reference_label,
                k_values=k_values,
            )
        ),
        jensen_shannon_ranking_stability=(
            distribution.ranking_stability(
                labels=labels,
                reference_label=config.reference_label,
                k_values=k_values,
            )
        ),
    )
