from pathlib import Path

import pytest

from src.dashboard import (
    DashboardConfig,
    DatasetRecord,
    analyze_records,
    load_demo_records,
    parse_fasta_bytes,
)
from src.genome import Genome


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_load_demo_records_returns_named_sequences():
    records = load_demo_records(PROJECT_ROOT)

    assert len(records) == 6
    assert records[0].label == "Aequorea GFP"
    assert records[-1].label == "Periodic control"


def test_parse_fasta_bytes_sanitizes_filename_and_builds_label():
    record = parse_fasta_bytes(
        "../example_sequence.fasta",
        b">Example\nACGTACGT\n",
    )

    assert record.source == "example_sequence.fasta"
    assert record.label == "example sequence"
    assert record.genome.sequence == "ACGTACGT"


def test_parse_fasta_bytes_rejects_multiple_records():
    with pytest.raises(
        ValueError,
        match=r"currently support one FASTA record\.",
    ):
        parse_fasta_bytes(
            "multiple.fasta",
            b">First\nACGT\n>Second\nACGT\n",
        )


def test_dashboard_config_requires_selected_k_in_range():
    with pytest.raises(
        ValueError,
        match=r"selected_k must be included in k_values\.",
    ):
        DashboardConfig(
            k_values=(1, 2, 3),
            selected_k=4,
            reference_label="A",
            comparison_label="B",
        )


def test_analyze_records_returns_all_comparison_families():
    records = [
        DatasetRecord("Reference", Genome("ACGTACGT"), "first.fasta"),
        DatasetRecord("Near", Genome("ACGTACGA"), "second.fasta"),
        DatasetRecord("Far", Genome("AAAACCCC"), "third.fasta"),
    ]
    config = DashboardConfig(
        k_values=(1, 2, 3),
        selected_k=2,
        reference_label="Reference",
        comparison_label="Near",
    )

    analysis = analyze_records(records, config)

    assert analysis.legacy_euclidean_matrix.metric == "euclidean"
    assert analysis.descriptor_v2_euclidean_matrix.metric == "euclidean_v2"
    assert analysis.embedding_v2_euclidean_matrix.metric == (
        "embedding_v2_euclidean"
    )
    assert analysis.jensen_shannon_matrix.metric == "jensen_shannon"
    assert list(analysis.jensen_shannon_pair_trajectory) == [1, 2, 3]


def test_dashboard_analysis_exports_json_and_csv():
    records = [
        DatasetRecord("Reference", Genome("ACGTACGT"), "first.fasta"),
        DatasetRecord("Comparison", Genome("AAAACCCC"), "second.fasta"),
    ]
    config = DashboardConfig(
        k_values=(1, 2),
        selected_k=2,
        reference_label="Reference",
        comparison_label="Comparison",
    )

    analysis = analyze_records(records, config)

    assert '"jensen_shannon"' in analysis.to_json()
    assert analysis.summary_csv().splitlines()[0] == "metric,value"


def test_dashboard_config_requires_two_scales():
    with pytest.raises(
        ValueError,
        match=r"At least two k-mer lengths are required\.",
    ):
        DashboardConfig(
            k_values=(3,),
            selected_k=3,
            reference_label="A",
            comparison_label="B",
        )
