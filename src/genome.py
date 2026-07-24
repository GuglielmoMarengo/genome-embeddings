import math, json, csv, io
from dataclasses import dataclass
from pathlib import Path


NORMALIZED_FEATURE_NAMES = [
    "gc_content",
    "normalized_shannon_entropy",
    "normalized_gc_skew",
    "purine_content",
    "kmer_diversity",
    "normalized_kmer_entropy",
]

SUPPORTED_MATRIX_METRICS = {
    "euclidean",
    "cosine",
}


@dataclass(slots=True)
class GenomeDescriptor:
    length: int
    gc_content: float
    at_content: float
    shannon_entropy: float
    gc_skew: float
    purine_content: float
    pyrimidine_content: float
    kmer_length: int
    kmer_diversity: float
    kmer_entropy: float

    def to_dict(self) -> dict[str, int | float]:
        return {
            "length": self.length,
            "gc_content": self.gc_content,
            "at_content": self.at_content,
            "shannon_entropy": self.shannon_entropy,
            "gc_skew": self.gc_skew,
            "purine_content": self.purine_content,
            "pyrimidine_content": self.pyrimidine_content,
            "kmer_length": self.kmer_length,
            "kmer_diversity": self.kmer_diversity,
            "kmer_entropy": self.kmer_entropy,
        }

    def to_vector(self) -> list[float]:
        return [
            float(self.length),
            float(self.gc_content),
            float(self.at_content),
            float(self.shannon_entropy),
            float(self.gc_skew),
            float(self.purine_content),
            float(self.pyrimidine_content),
            float(self.kmer_length),
            float(self.kmer_diversity),
            float(self.kmer_entropy),
        ]

    def to_normalized_vector(self) -> list[float]:
        normalized_shannon_entropy = (
            self.shannon_entropy / 2
        )

        normalized_gc_skew = (
            self.gc_skew + 1
        ) / 2

        normalized_kmer_entropy = (
            self.kmer_entropy
            / (2 * self.kmer_length)
        )

        return [
            self.gc_content,
            normalized_shannon_entropy,
            normalized_gc_skew,
            self.purine_content,
            self.kmer_diversity,
            normalized_kmer_entropy,
        ]

    def euclidean_distance(
        self,
        other: "GenomeDescriptor",
    ) -> float:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        squared_differences = [
            (first_value - second_value) ** 2
            for first_value, second_value in zip(
                first_vector,
                second_vector,
                strict=True,
            )
        ]

        return math.sqrt(sum(squared_differences))

    def cosine_similarity(
        self,
        other: "GenomeDescriptor",
    ) -> float:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        dot_product = sum(
            first_value * second_value
            for first_value, second_value in zip(
                first_vector,
                second_vector,
                strict=True,
            )
        )

        first_magnitude = math.sqrt(
            sum(
                value ** 2
                for value in first_vector
            )
        )

        second_magnitude = math.sqrt(
            sum(
                value ** 2
                for value in second_vector
            )
        )

        if (
            first_magnitude == 0
            or second_magnitude == 0
        ):
            return 0.0

        return (
            dot_product
            / (first_magnitude * second_magnitude)
        )

    def feature_differences(
        self,
        other: "GenomeDescriptor",
    ) -> dict[str, float]:
        self._validate_other_descriptor(other)

        first_vector = self.to_normalized_vector()
        second_vector = other.to_normalized_vector()

        return {
            feature_name: abs(
                first_value - second_value
            )
            for (
                feature_name,
                first_value,
                second_value,
            ) in zip(
                NORMALIZED_FEATURE_NAMES,
                first_vector,
                second_vector,
                strict=True,
            )
        }

    def compare(
        self,
        other: "GenomeDescriptor",
    ) -> "GenomeComparison":
        self._validate_other_descriptor(other)

        return GenomeComparison(
            euclidean_distance=(
                self.euclidean_distance(other)
            ),
            cosine_similarity=(
                self.cosine_similarity(other)
            ),
            feature_differences=(
                self.feature_differences(other)
            ),
        )

    @staticmethod
    def _validate_other_descriptor(
        other: object,
    ) -> None:
        if not isinstance(
            other,
            GenomeDescriptor,
        ):
            raise TypeError(
                "other must be a GenomeDescriptor."
            )


@dataclass(slots=True)
class GenomeComparison:
    euclidean_distance: float
    cosine_similarity: float
    feature_differences: dict[str, float]

    def sorted_feature_differences(
        self,
    ) -> list[tuple[str, float]]:
        return sorted(
            self.feature_differences.items(),
            key=lambda item: item[1],
            reverse=True,
        )


@dataclass(slots=True)
class GenomeMatrix:
    labels: list[str]
    values: list[list[float]]
    metric: str
    kmer_length: int

    def __post_init__(self) -> None:
        if not self.labels:
            raise ValueError(
                "Genome matrix labels cannot be empty."
            )

        matrix_size = len(self.values)

        if any(
            len(row) != matrix_size
            for row in self.values
        ):
            raise ValueError(
                "Genome matrix values must be square."
            )

        if len(self.labels) != matrix_size:
            raise ValueError(
                "Genome matrix labels must match matrix size."
            )

        if self.metric not in SUPPORTED_MATRIX_METRICS:
            raise ValueError(
                "Unsupported genome matrix metric."
            )
    
    def get_value(
        self,
        row_label: str,
        column_label: str,
    ) -> float:
        row_index = self._label_index(row_label)
        column_index = self._label_index(column_label)

        return self.values[row_index][column_index]
    
    def to_rows(
        self,
    ) -> list[dict[str, str | list[float]]]:
        return [
            {
                "label": label,
                "values": row.copy(),
            }
            for label, row in zip(
                self.labels,
                self.values,
                strict=True,
            )
        ]
    
    def to_dict(
        self,
    ) -> dict[
        str,
        list[str]
        | list[list[float]]
        | str
        | int,
    ]:
        return {
            "labels": self.labels.copy(),
            "values": [
                row.copy()
                for row in self.values
            ],
            "metric": self.metric,
            "kmer_length": self.kmer_length,
        }
    
    def to_json(
        self,
        indent: int | None = None,
    ) -> str:
        return json.dumps(
            self.to_dict(),
            indent=indent,
        )

    def to_csv(
        self,
        delimiter: str = ",",
    ) -> str:
        output = io.StringIO()

        writer = csv.writer(
            output,
            delimiter=delimiter,
            lineterminator="\n",
        )

        writer.writerow(
            [
                "label",
                *self.labels,
            ]
        )

        for label, row in zip(
            self.labels,
            self.values,
            strict=True,
        ):
            writer.writerow(
                [
                    label,
                    *row,
                ]
            )

        return output.getvalue()
    
    def to_upper_triangle_vector(
        self,
    ) -> list[float]:
        return [
            self.values[row_index][column_index]
            for row_index in range(len(self.values))
            for column_index in range(
                row_index + 1,
                len(self.values),
            )
        ]
    
    def upper_triangle_pairs(
        self,
    ) -> list[tuple[str, str]]:
        return [
            (
                self.labels[row_index],
                self.labels[column_index],
            )
            for row_index in range(len(self.labels))
            for column_index in range(
                row_index + 1,
                len(self.labels),
            )
        ]
    
    def to_upper_triangle_rows(
        self,
    ) -> list[dict[str, str | float]]:
        return [
            {
                "row_label": self.labels[row_index],
                "column_label": self.labels[column_index],
                "value": self.values[row_index][column_index],
            }
            for row_index in range(len(self.labels))
            for column_index in range(
                row_index + 1,
                len(self.labels),
            )
        ]
    
    def rank_by_label(
        self,
        label: str,
    ) -> list[tuple[str, float]]:
        reference_index = self._label_index(label)

        ranking = [
            (
                candidate_label,
                self.values[reference_index][candidate_index],
            )
            for candidate_index, candidate_label in enumerate(
                self.labels
            )
            if candidate_index != reference_index
        ]

        reverse = self.metric == "cosine"

        return sorted(
            ranking,
            key=lambda item: item[1],
            reverse=reverse,
        )

    def _label_index(
        self,
        label: str,
    ) -> int:
        try:
            return self.labels.index(label)
        except ValueError as error:
            raise ValueError(
                f"Unknown genome matrix label: {label}."
            ) from error


class GenomeCollection:
    def __init__(
        self,
        genomes: list["Genome"],
    ) -> None:
        if not isinstance(genomes, list):
            raise TypeError(
                "genomes must be a list of Genome objects."
            )

        if not genomes:
            raise ValueError(
                "Genome collection cannot be empty."
            )

        if not all(
            isinstance(genome, Genome)
            for genome in genomes
        ):
            raise TypeError(
                "All items must be Genome objects."
            )

        self.genomes = genomes

    def descriptors(
        self,
        k: int,
    ) -> list[GenomeDescriptor]:
        return [
            genome.descriptor(k=k)
            for genome in self.genomes
        ]
    
    def _build_matrices(
        self,
        labels: list[str],
        k_values: list[int],
        matrix_builder,
    ) -> dict[int, GenomeMatrix]:
        if not k_values:
            raise ValueError(
                "k-mer lengths cannot be empty."
            )

        if len(k_values) != len(set(k_values)):
            raise ValueError(
                "k-mer lengths must be unique."
            )

        return {
            k: matrix_builder(
                labels=labels,
                k=k,
            )
            for k in k_values
        }
    
    def _build_matrix_trajectory(
        self,
        labels: list[str],
        k_values: list[int],
        matrices_builder,
    ) -> dict[int, list[float]]:
        matrices = matrices_builder(
            labels=labels,
            k_values=k_values,
        )

        return {
            k: matrix.to_upper_triangle_vector()
            for k, matrix in matrices.items()
        }
    
    def _build_pair_trajectory(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
        matrices_builder,
    ) -> dict[int, float]:
        matrices = matrices_builder(
            labels=labels,
            k_values=k_values,
        )

        return {
            k: matrix.get_value(
                row_label=row_label,
                column_label=column_label,
            )
            for k, matrix in matrices.items()
        }

    @staticmethod
    def _validate_trajectory_steps(
        trajectory: dict[int, object],
        trajectory_name: str,
    ) -> None:
        if not trajectory:
            raise ValueError(
                f"{trajectory_name} cannot be empty."
            )

        if len(trajectory) < 2:
            raise ValueError(
                f"{trajectory_name} must contain "
                "at least two k-mer scales."
            )

    @staticmethod
    def pair_trajectory_step_differences(
        trajectory: dict[int, float],
    ) -> dict[tuple[int, int], float]:
        GenomeCollection._validate_trajectory_steps(
            trajectory=trajectory,
            trajectory_name="Pair trajectory",
        )

        trajectory_items = list(trajectory.items())

        return {
            (first_k, second_k): (
                second_value - first_value
            )
            for (
                first_k,
                first_value,
            ), (
                second_k,
                second_value,
            ) in zip(
                trajectory_items,
                trajectory_items[1:],
            )
        }

    @staticmethod
    def matrix_trajectory_step_distances(
        trajectory: dict[int, list[float]],
    ) -> dict[tuple[int, int], float]:
        GenomeCollection._validate_trajectory_steps(
            trajectory=trajectory,
            trajectory_name="Matrix trajectory",
        )

        vector_lengths = {
            len(vector)
            for vector in trajectory.values()
        }

        if 0 in vector_lengths:
            raise ValueError(
                "Matrix trajectory vectors cannot be empty."
            )

        if len(vector_lengths) != 1:
            raise ValueError(
                "Matrix trajectory vectors must have "
                "the same length."
            )

        trajectory_items = list(trajectory.items())

        return {
            (first_k, second_k): math.sqrt(
                sum(
                    (second_value - first_value) ** 2
                    for (
                        first_value,
                        second_value,
                    ) in zip(
                        first_vector,
                        second_vector,
                        strict=True,
                    )
                )
            )
            for (
                first_k,
                first_vector,
            ), (
                second_k,
                second_vector,
            ) in zip(
                trajectory_items,
                trajectory_items[1:],
            )
        }

    @staticmethod
    def matrix_trajectory_pair_contributions(
        labels: list[str],
        trajectory: dict[int, list[float]],
    ) -> dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ]:
        GenomeCollection._validate_trajectory_steps(
            trajectory=trajectory,
            trajectory_name="Matrix trajectory",
        )

        if not labels:
            raise ValueError(
                "Genome labels cannot be empty."
            )

        if len(labels) != len(set(labels)):
            raise ValueError(
                "Genome labels must be unique."
            )

        pairs = [
            (labels[row_index], labels[column_index])
            for row_index in range(len(labels))
            for column_index in range(
                row_index + 1,
                len(labels),
            )
        ]

        expected_vector_length = len(pairs)

        if any(
            len(vector) != expected_vector_length
            for vector in trajectory.values()
        ):
            raise ValueError(
                "Matrix trajectory vector length must "
                "match the number of unique label pairs."
            )

        trajectory_items = list(trajectory.items())
        contributions: dict[
            tuple[int, int],
            list[dict[str, str | float]],
        ] = {}

        for (
            first_k,
            first_vector,
        ), (
            second_k,
            second_vector,
        ) in zip(
            trajectory_items,
            trajectory_items[1:],
        ):
            step_rows = [
                {
                    "row_label": row_label,
                    "column_label": column_label,
                    "difference": (
                        second_value - first_value
                    ),
                    "absolute_difference": abs(
                        second_value - first_value
                    ),
                }
                for (
                    row_label,
                    column_label,
                ), (
                    first_value,
                    second_value,
                ) in zip(
                    pairs,
                    zip(
                        first_vector,
                        second_vector,
                        strict=True,
                    ),
                    strict=True,
                )
            ]

            contributions[(first_k, second_k)] = sorted(
                step_rows,
                key=lambda row: row[
                    "absolute_difference"
                ],
                reverse=True,
            )

        return contributions

    @staticmethod
    def matrix_trajectory_deformation_partition(
        contributions: dict[
            tuple[int, int],
            list[dict[str, str | float]],
        ],
        selected_label: str,
    ) -> dict[
        tuple[int, int],
        dict[str, float | int],
    ]:
        Genome.validate_string(
            "Selected label",
            selected_label,
        )

        if not contributions:
            raise ValueError(
                "Matrix trajectory contributions "
                "cannot be empty."
            )

        partitions: dict[
            tuple[int, int],
            dict[str, float | int],
        ] = {}

        for transition, rows in contributions.items():
            if not rows:
                raise ValueError(
                    "Matrix trajectory contribution "
                    "rows cannot be empty."
                )

            total_squared_deformation = 0.0
            selected_squared_deformation = 0.0
            selected_pair_count = 0

            for row in rows:
                difference = float(row["difference"])
                squared_difference = difference ** 2

                total_squared_deformation += (
                    squared_difference
                )

                if (
                    row["row_label"] == selected_label
                    or row["column_label"]
                    == selected_label
                ):
                    selected_squared_deformation += (
                        squared_difference
                    )
                    selected_pair_count += 1

            if selected_pair_count == 0:
                raise ValueError(
                    "Selected label is not present in "
                    "the contribution rows."
                )

            remaining_squared_deformation = (
                total_squared_deformation
                - selected_squared_deformation
            )

            selected_share = (
                selected_squared_deformation
                / total_squared_deformation
                if total_squared_deformation != 0
                else 0.0
            )

            remaining_share = (
                remaining_squared_deformation
                / total_squared_deformation
                if total_squared_deformation != 0
                else 0.0
            )

            partitions[transition] = {
                "total_squared_deformation": (
                    total_squared_deformation
                ),
                "selected_squared_deformation": (
                    selected_squared_deformation
                ),
                "remaining_squared_deformation": (
                    remaining_squared_deformation
                ),
                "selected_share": selected_share,
                "remaining_share": remaining_share,
                "total_distance": math.sqrt(
                    total_squared_deformation
                ),
                "selected_distance": math.sqrt(
                    selected_squared_deformation
                ),
                "remaining_distance": math.sqrt(
                    remaining_squared_deformation
                ),
                "selected_pair_count": (
                    selected_pair_count
                ),
                "remaining_pair_count": (
                    len(rows) - selected_pair_count
                ),
            }

        return partitions

    def euclidean_distance_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k=k)

        values = self._build_matrix(
            descriptors=descriptors,
            metric="euclidean",
        )

        return GenomeMatrix(
            labels=labels,
            values=values,
            metric="euclidean",
            kmer_length=k,
        )
    
    def euclidean_distance_matrices(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, GenomeMatrix]:
        return self._build_matrices(
            labels=labels,
            k_values=k_values,
            matrix_builder=self.euclidean_distance_matrix,
        )

    def euclidean_matrix_trajectory(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, list[float]]:
        return self._build_matrix_trajectory(
            labels=labels,
            k_values=k_values,
            matrices_builder=self.euclidean_distance_matrices,
        )
    
    def euclidean_pair_trajectory(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
    ) -> dict[int, float]:
        return self._build_pair_trajectory(
            labels=labels,
            row_label=row_label,
            column_label=column_label,
            k_values=k_values,
            matrices_builder=self.euclidean_distance_matrices,
        )

    def euclidean_pair_trajectory_step_differences(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
    ) -> dict[tuple[int, int], float]:
        trajectory = self.euclidean_pair_trajectory(
            labels=labels,
            row_label=row_label,
            column_label=column_label,
            k_values=k_values,
        )

        return self.pair_trajectory_step_differences(
            trajectory
        )

    def euclidean_matrix_trajectory_step_distances(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[tuple[int, int], float]:
        trajectory = self.euclidean_matrix_trajectory(
            labels=labels,
            k_values=k_values,
        )

        return self.matrix_trajectory_step_distances(
            trajectory
        )

    def euclidean_matrix_pair_contributions(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ]:
        trajectory = self.euclidean_matrix_trajectory(
            labels=labels,
            k_values=k_values,
        )

        return self.matrix_trajectory_pair_contributions(
            labels=labels,
            trajectory=trajectory,
        )

    def euclidean_matrix_deformation_partition(
        self,
        labels: list[str],
        k_values: list[int],
        selected_label: str,
    ) -> dict[
        tuple[int, int],
        dict[str, float | int],
    ]:
        contributions = (
            self.euclidean_matrix_pair_contributions(
                labels=labels,
                k_values=k_values,
            )
        )

        return self.matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label=selected_label,
        )

    def cosine_similarity_matrix(
        self,
        labels: list[str],
        k: int,
    ) -> GenomeMatrix:
        descriptors = self.descriptors(k=k)

        values = self._build_matrix(
            descriptors=descriptors,
            metric="cosine",
        )

        return GenomeMatrix(
            labels=labels,
            values=values,
            metric="cosine",
            kmer_length=k,
        )
    
    def cosine_similarity_matrices(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, GenomeMatrix]:
        return self._build_matrices(
            labels=labels,
            k_values=k_values,
            matrix_builder=self.cosine_similarity_matrix,
        )

    def cosine_matrix_trajectory(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[int, list[float]]:
        return self._build_matrix_trajectory(
            labels=labels,
            k_values=k_values,
            matrices_builder=self.cosine_similarity_matrices,
        )
    
    def cosine_pair_trajectory(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
    ) -> dict[int, float]:
        return self._build_pair_trajectory(
            labels=labels,
            row_label=row_label,
            column_label=column_label,
            k_values=k_values,
            matrices_builder=self.cosine_similarity_matrices,
        )

    def cosine_pair_trajectory_step_differences(
        self,
        labels: list[str],
        row_label: str,
        column_label: str,
        k_values: list[int],
    ) -> dict[tuple[int, int], float]:
        trajectory = self.cosine_pair_trajectory(
            labels=labels,
            row_label=row_label,
            column_label=column_label,
            k_values=k_values,
        )

        return self.pair_trajectory_step_differences(
            trajectory
        )

    def cosine_matrix_trajectory_step_distances(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[tuple[int, int], float]:
        trajectory = self.cosine_matrix_trajectory(
            labels=labels,
            k_values=k_values,
        )

        return self.matrix_trajectory_step_distances(
            trajectory
        )

    def cosine_matrix_pair_contributions(
        self,
        labels: list[str],
        k_values: list[int],
    ) -> dict[
        tuple[int, int],
        list[dict[str, str | float]],
    ]:
        trajectory = self.cosine_matrix_trajectory(
            labels=labels,
            k_values=k_values,
        )

        return self.matrix_trajectory_pair_contributions(
            labels=labels,
            trajectory=trajectory,
        )

    def cosine_matrix_deformation_partition(
        self,
        labels: list[str],
        k_values: list[int],
        selected_label: str,
    ) -> dict[
        tuple[int, int],
        dict[str, float | int],
    ]:
        contributions = (
            self.cosine_matrix_pair_contributions(
                labels=labels,
                k_values=k_values,
            )
        )

        return self.matrix_trajectory_deformation_partition(
            contributions=contributions,
            selected_label=selected_label,
        )

    @staticmethod
    def _build_matrix(
        descriptors: list[GenomeDescriptor],
        metric: str,
    ) -> list[list[float]]:
        matrix: list[list[float]] = []

        for first_descriptor in descriptors:
            row: list[float] = []

            for second_descriptor in descriptors:
                if metric == "euclidean":
                    value = (
                        first_descriptor.euclidean_distance(
                            second_descriptor
                        )
                    )
                elif metric == "cosine":
                    value = (
                        first_descriptor.cosine_similarity(
                            second_descriptor
                        )
                    )
                else:
                    raise ValueError(
                        "Unsupported genome matrix metric."
                    )

                row.append(value)

            matrix.append(row)

        return matrix


class Genome:
    VALID_NUCLEOTIDES = {
        "A",
        "C",
        "G",
        "T",
    }

    GC_BASES = {
        "C",
        "G",
    }

    COMPLEMENT = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    def __init__(
        self,
        sequence: str,
        header: str | None = None,
    ) -> None:
        self.validate_string(
            "Sequence",
            sequence,
        )

        normalized_sequence = sequence.upper()

        self._validate_sequence(
            normalized_sequence
        )

        self.sequence = normalized_sequence
        self.header = header

    @classmethod
    def from_fasta(
        cls,
        filepath: str | Path,
    ) -> "Genome":
        with open(
            filepath,
            encoding="utf-8",
        ) as file:
            lines = [
                line.strip()
                for line in file
                if line.strip()
            ]

        if not lines:
            raise ValueError(
                "FASTA file is empty."
            )

        header = lines[0]

        if not header.startswith(">"):
            raise ValueError(
                "FASTA header must start with '>'."
            )

        sequence = "".join(lines[1:])

        return cls(
            sequence=sequence,
            header=header,
        )

    @staticmethod
    def validate_string(
        name: str,
        value: object,
    ) -> None:
        if type(value) is not str:
            raise TypeError(
                f"{name} must be a string."
            )

        if len(value) == 0:
            raise ValueError(
                f"{name} cannot be empty."
            )

    @classmethod
    def _validate_sequence(
        cls,
        sequence: str,
    ) -> None:
        for position, character in enumerate(
            sequence,
            start=1,
        ):
            if character not in cls.VALID_NUCLEOTIDES:
                raise ValueError(
                    f"Invalid character {character} "
                    f"in sequence at position {position}."
                )

    @staticmethod
    def _entropy_from_frequencies(
        frequencies: dict[str, int],
        total: int,
    ) -> float:
        entropy = 0.0

        for frequency in frequencies.values():
            probability = frequency / total
            entropy -= (
                probability
                * math.log2(probability)
            )

        return entropy

    @staticmethod
    def _validate_k(
        k: int,
        sequence_length: int,
    ) -> None:
        if type(k) is not int:
            raise TypeError(
                "k must be an integer, "
                f"got {type(k).__name__}."
            )

        if k <= 0:
            raise ValueError(
                f"k must be positive. Got {k}."
            )

        if k > sequence_length:
            raise ValueError(
                "k cannot exceed the sequence length. "
                f"Got {k}."
            )

    def length(self) -> int:
        return len(self.sequence)

    def gc_content(self) -> float:
        gc_count = sum(
            1
            for character in self.sequence
            if character in self.GC_BASES
        )

        return gc_count / self.length()

    def at_content(self) -> float:
        return 1.0 - self.gc_content()

    def reverse_complement(self) -> str:
        return "".join(
            self.COMPLEMENT[character]
            for character in reversed(self.sequence)
        )

    def kmers(
        self,
        k: int,
    ) -> list[str]:
        self._validate_k(
            k=k,
            sequence_length=self.length(),
        )

        return [
            self.sequence[index:index + k]
            for index in range(
                self.length() - k + 1
            )
        ]

    def kmer_frequencies(
        self,
        k: int,
    ) -> dict[str, int]:
        frequencies: dict[str, int] = {}

        for kmer in self.kmers(k):
            frequencies[kmer] = (
                frequencies.get(kmer, 0) + 1
            )

        return frequencies

    def nucleotide_frequencies(
        self,
    ) -> dict[str, int]:
        frequencies: dict[str, int] = {}

        for character in self.sequence:
            frequencies[character] = (
                frequencies.get(character, 0) + 1
            )

        return frequencies

    def shannon_entropy(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        return self._entropy_from_frequencies(
            frequencies=frequencies,
            total=self.length(),
        )

    def gc_skew(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        g_count = frequencies.get("G", 0)
        c_count = frequencies.get("C", 0)
        total_gc = g_count + c_count

        if total_gc == 0:
            return 0.0

        return (
            g_count - c_count
        ) / total_gc

    def purine_content(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        a_count = frequencies.get("A", 0)
        g_count = frequencies.get("G", 0)

        return (
            a_count + g_count
        ) / self.length()

    def pyrimidine_content(self) -> float:
        frequencies = (
            self.nucleotide_frequencies()
        )

        c_count = frequencies.get("C", 0)
        t_count = frequencies.get("T", 0)

        return (
            c_count + t_count
        ) / self.length()

    def kmer_diversity(
        self,
        k: int,
    ) -> float:
        frequencies = self.kmer_frequencies(k)
        total_kmers = self.length() - k + 1

        distinct_kmers = len(frequencies)

        maximum_distinct_kmers = min(
            total_kmers,
            len(self.VALID_NUCLEOTIDES) ** k,
        )

        return (
            distinct_kmers
            / maximum_distinct_kmers
        )

    def kmer_entropy(
        self,
        k: int,
    ) -> float:
        frequencies = self.kmer_frequencies(k)
        total_kmers = self.length() - k + 1

        return self._entropy_from_frequencies(
            frequencies=frequencies,
            total=total_kmers,
        )

    def descriptor(
        self,
        k: int = 3,
    ) -> GenomeDescriptor:
        return GenomeDescriptor(
            length=self.length(),
            gc_content=self.gc_content(),
            at_content=self.at_content(),
            shannon_entropy=self.shannon_entropy(),
            gc_skew=self.gc_skew(),
            purine_content=self.purine_content(),
            pyrimidine_content=(
                self.pyrimidine_content()
            ),
            kmer_length=k,
            kmer_diversity=(
                self.kmer_diversity(k)
            ),
            kmer_entropy=(
                self.kmer_entropy(k)
            ),
        )