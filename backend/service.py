"""Conversions between API payloads and the scientific analysis layer."""

from __future__ import annotations

from pathlib import Path

from backend.schemas import AnalysisConfigPayload, AnalysisRequest, SequencePayload
from src.dashboard import (
    DashboardAnalysis,
    DashboardConfig,
    DatasetRecord,
    analyze_records,
    load_demo_records,
    parse_fasta_bytes,
)
from src.genome import Genome


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def record_to_payload(record: DatasetRecord) -> SequencePayload:
    return SequencePayload(
        label=record.label,
        sequence=record.genome.sequence,
        header=record.genome.header or "",
        source=record.source,
    )


def payload_to_record(payload: SequencePayload) -> DatasetRecord:
    return DatasetRecord(
        label=payload.label,
        genome=Genome(
            sequence=payload.sequence,
            header=payload.header or None,
        ),
        source=Path(payload.source).name or "uploaded",
    )


def config_to_dashboard(payload: AnalysisConfigPayload) -> DashboardConfig:
    return DashboardConfig(
        k_values=tuple(payload.k_values),
        selected_k=payload.selected_k,
        reference_label=payload.reference_label,
        comparison_label=payload.comparison_label,
    )


def run_request(request: AnalysisRequest) -> DashboardAnalysis:
    records = [payload_to_record(record) for record in request.records]
    return analyze_records(records, config_to_dashboard(request.config))


def demo_request() -> AnalysisRequest:
    records = load_demo_records(PROJECT_ROOT)
    payloads = [record_to_payload(record) for record in records]
    return AnalysisRequest(
        records=payloads,
        config=AnalysisConfigPayload(
            k_values=[1, 2, 3, 4, 5],
            selected_k=3,
            reference_label=payloads[0].label,
            comparison_label=payloads[1].label,
        ),
    )


def parse_upload(filename: str, data: bytes, label: str | None = None) -> SequencePayload:
    return record_to_payload(parse_fasta_bytes(filename, data, label=label))
