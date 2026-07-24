"""Validated API schemas for the Genome Embeddings dashboard."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class SequencePayload(BaseModel):
    """A single validated DNA sequence used by the web client."""

    model_config = ConfigDict(str_strip_whitespace=True)

    label: str = Field(min_length=1, max_length=160)
    sequence: str = Field(min_length=1)
    header: str = ""
    source: str = "uploaded"

    @field_validator("sequence")
    @classmethod
    def normalize_sequence(cls, value: str) -> str:
        return "".join(value.split()).upper()


class AnalysisConfigPayload(BaseModel):
    """Client-selected multiscale analysis settings."""

    k_values: list[int] = Field(min_length=2)
    selected_k: int = Field(ge=1)
    reference_label: str = Field(min_length=1)
    comparison_label: str = Field(min_length=1)

    @field_validator("k_values")
    @classmethod
    def validate_k_values(cls, values: list[int]) -> list[int]:
        if any(type(value) is not int or value <= 0 for value in values):
            raise ValueError("k_values must contain positive integers.")
        if len(values) != len(set(values)):
            raise ValueError("k_values must be unique.")
        return sorted(values)

    @model_validator(mode="after")
    def validate_relationships(self) -> "AnalysisConfigPayload":
        if self.selected_k not in self.k_values:
            raise ValueError("selected_k must be included in k_values.")
        if self.reference_label == self.comparison_label:
            raise ValueError("Reference and comparison labels must differ.")
        return self


class AnalysisRequest(BaseModel):
    """Stateless request containing a complete dataset and configuration."""

    records: list[SequencePayload] = Field(min_length=2)
    config: AnalysisConfigPayload


class DemoResponse(BaseModel):
    """Demonstration records together with a complete initial analysis."""

    records: list[SequencePayload]
    analysis: dict[str, object]


class ApiMetadata(BaseModel):
    """Discoverable information used by the frontend shell."""

    name: str
    tagline: str
    version: str
    supported_k_min: int
    supported_k_max: int
    upload_extensions: list[str]
