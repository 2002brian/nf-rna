"""Strict Pydantic models for the Milestone 1 project contract."""

from __future__ import annotations

import re
from enum import Enum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, StrictStr, field_validator, model_validator

SUPPORTED_SCHEMA_VERSION = "1.2"
LEGACY_SCHEMA_VERSION = "1.0"
PIPELINE_VERSION = "0.5.0"
NFCORE_RNASEQ_VERSION = "3.26.0"
PROJECT_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class StrictModel(BaseModel):
    """Base model that rejects unrecognized configuration fields."""

    model_config = ConfigDict(extra="forbid")


class Preset(str, Enum):
    QC = "QC"
    L1 = "L1"
    L2 = "L2"


class Species(str, Enum):
    HUMAN = "Homo sapiens"
    MOUSE = "Mus musculus"
    RAT = "Rattus norvegicus"
    OTHER = "Other"


class DesignType(str, Enum):
    TWO_GROUP = "two_group"
    MULTI_GROUP = "multi_group"
    PAIRED = "paired"


class InputType(str, Enum):
    FASTQ = "fastq"
    RAW_COUNTS = "raw_counts"


class SequencingLayout(str, Enum):
    PAIRED_END = "paired_end"
    SINGLE_END = "single_end"


class FastqPreprocessing(str, Enum):
    RAW = "raw"
    PRETRIMMED = "pretrimmed"


class ProjectInfo(StrictModel):
    id: StrictStr
    pipeline: Literal["bulk_rnaseq"]
    preset: Preset

    @field_validator("id")
    @classmethod
    def validate_project_id(cls, value: str) -> str:
        if not PROJECT_ID_PATTERN.fullmatch(value):
            raise ValueError(
                "Project ID must start with a letter or digit and contain only "
                "letters, digits, '.', '_' or '-'."
            )
        return value


class OrganismConfig(StrictModel):
    species: Species


class InputConfig(StrictModel):
    type: InputType
    path: StrictStr
    layout: SequencingLayout | None = None
    preprocessing: FastqPreprocessing | None = None

    @model_validator(mode="after")
    def validate_layout(self) -> "InputConfig":
        if self.type is InputType.FASTQ and self.layout is None:
            raise ValueError("input.layout is required when input.type is fastq.")
        if self.type is InputType.FASTQ and self.preprocessing is None:
            # Legacy FASTQ projects omitted this declaration. Preserve their
            # trimming-enabled behavior while exposing it after parsing.
            self.preprocessing = FastqPreprocessing.RAW
        if self.type is InputType.RAW_COUNTS and self.layout is not None:
            raise ValueError("input.layout is only supported when input.type is fastq.")
        if self.type is InputType.RAW_COUNTS and self.preprocessing is not None:
            raise ValueError("input.preprocessing is only supported when input.type is fastq.")
        return self


class DesignConfig(StrictModel):
    type: DesignType
    formula: StrictStr


class ThresholdsConfig(StrictModel):
    padj: float = Field(ge=0.0, le=1.0)
    abs_log2fc: float = Field(ge=0.0)


class GoConfig(StrictModel):
    pvalue_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    qvalue_cutoff: float = Field(default=0.2, ge=0.0, le=1.0)
    p_adjust_method: Literal["BH"] = "BH"


class GseaConfig(StrictModel):
    """Deliberately narrow GO-only preranked GSEA contract for M4B-2."""

    minimum_ranked_genes: int = Field(default=50, ge=1)
    min_gs_size: int = Field(default=10, ge=1)
    max_gs_size: int = Field(default=500, ge=1)
    pvalue_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    padj_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    p_adjust_method: Literal["BH"] = "BH"
    seed: int = Field(default=1, ge=0)

    @model_validator(mode="after")
    def validate_size_bounds(self) -> "GseaConfig":
        if self.max_gs_size < self.min_gs_size:
            raise ValueError("annotation.enrichment.gsea.max_gs_size must be >= min_gs_size.")
        return self


class KeggOraConfig(StrictModel):
    pvalue_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    qvalue_cutoff: float = Field(default=0.2, ge=0.0, le=1.0)
    p_adjust_method: Literal["BH"] = "BH"
    min_gs_size: int = Field(default=10, ge=1)
    max_gs_size: int = Field(default=500, ge=1)

    @model_validator(mode="after")
    def validate_size_bounds(self) -> "KeggOraConfig":
        if self.max_gs_size < self.min_gs_size:
            raise ValueError("annotation.enrichment.kegg.ora.max_gs_size must be >= min_gs_size.")
        return self


class KeggGseaConfig(StrictModel):
    minimum_ranked_genes: int = Field(default=50, ge=1)
    pvalue_cutoff: float = Field(default=1.0, ge=0.0, le=1.0)
    padj_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    p_adjust_method: Literal["BH"] = "BH"
    min_gs_size: int = Field(default=10, ge=1)
    max_gs_size: int = Field(default=500, ge=1)
    seed: int = Field(default=1, ge=0)

    @model_validator(mode="after")
    def validate_size_bounds(self) -> "KeggGseaConfig":
        if self.max_gs_size < self.min_gs_size:
            raise ValueError("annotation.enrichment.kegg.gsea.max_gs_size must be >= min_gs_size.")
        return self


class KeggConfig(StrictModel):
    """Online-only KEGG contract; no provider fallback is implemented."""

    resource_provider: Literal["online_kegg_rest_via_clusterprofiler"] = "online_kegg_rest_via_clusterprofiler"
    ora: KeggOraConfig = Field(default_factory=KeggOraConfig)
    gsea: KeggGseaConfig = Field(default_factory=KeggGseaConfig)


class EnrichmentConfig(StrictModel):
    go: GoConfig = Field(default_factory=GoConfig)
    gsea: GseaConfig = Field(default_factory=GseaConfig)
    kegg: KeggConfig = Field(default_factory=KeggConfig)


class AnnotationConfig(StrictModel):
    """Deliberately narrow, offline-only identifier contract for M4B-1."""

    organism: Literal["Homo sapiens", "Mus musculus"]
    input_id_type: Literal["ENSEMBL", "ENTREZID", "SYMBOL"]
    target_id_type: Literal["ENTREZID"] = "ENTREZID"
    gene_symbol_output: bool = True
    minimum_mapping_rate: float = Field(default=0.70, ge=0.0, le=1.0)
    minimum_mapped_foreground: int = Field(default=5, ge=1)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)


PUBLIC_ENRICHMENT_METHOD = "gsea"
LEGACY_ENRICHMENT_MODULES = frozenset({"go", "gsea-go", "kegg", "gsea-kegg"})
INTERNAL_GSEA_BACKENDS = ("gsea-go", "gsea-kegg")


def normalize_enrichment_selection(value: object) -> tuple[str, ...]:
    """Normalize only the documented public method and complete legacy scope."""

    if value == PUBLIC_ENRICHMENT_METHOD:
        return (PUBLIC_ENRICHMENT_METHOD,)
    if isinstance(value, (list, tuple)):
        selected = tuple(value)
        if not selected:
            return ()
        if selected == (PUBLIC_ENRICHMENT_METHOD,):
            return selected
        if len(selected) == len(LEGACY_ENRICHMENT_MODULES) and set(selected) == LEGACY_ENRICHMENT_MODULES:
            return (PUBLIC_ENRICHMENT_METHOD,)
        if any(item in LEGACY_ENRICHMENT_MODULES for item in selected):
            raise ValueError(
                "Legacy analysis.enrichment is accepted only when it contains exactly "
                "go, gsea-go, kegg, gsea-kegg; use analysis.enrichment: gsea."
            )
    return value  # Pydantic reports invalid types or unsupported public values.


def production_enrichment_backends(value: object) -> tuple[str, ...]:
    """Expand the one public production method into its two internal backends."""

    selected = normalize_enrichment_selection(value)
    if selected == ():
        return ()
    if selected == (PUBLIC_ENRICHMENT_METHOD,):
        return INTERNAL_GSEA_BACKENDS
    raise ValueError("analysis.enrichment must be empty or the public method: gsea.")


class AnalysisConfig(StrictModel):
    """Execution choices, separated from statistical thresholds.

    An explicit list makes the delivery contract reproducible: a later run never
    gains an enrichment step merely because the CLI happened to use a flag.
    """

    enrichment: tuple[Literal["gsea"], ...] = ()

    @field_validator("enrichment", mode="before")
    @classmethod
    def normalize_legacy_enrichment(cls, value: object) -> object:
        return normalize_enrichment_selection(value)

    @field_validator("enrichment")
    @classmethod
    def validate_unique_enrichment(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        if len(set(value)) != len(value):
            raise ValueError("analysis.enrichment must not contain duplicate modules.")
        return value


class UpstreamConfig(StrictModel):
    engine: Literal["nfcore_rnaseq", "external"]
    pipeline_version: StrictStr | None = None
    aligner: StrictStr | None = None
    strandedness: Literal["auto", "unstranded", "forward", "reverse"] | None = None
    provider: StrictStr | None = None
    quantification_method: StrictStr | None = None
    quantification: "QuantificationConfig | None" = None


class QuantificationConfig(StrictModel):
    """Explicit FASTQ quantification route.

    ``salmon`` remains the compatibility default for projects written before
    this field existed.  Alignment/counting is deliberately a separate route:
    it must never inherit tximport semantics.
    """

    method: Literal["salmon", "hisat2_featurecounts"]


class ReferenceConfig(StrictModel):
    source: Literal["igenomes", "custom", "local"]
    genome: StrictStr | None = None
    fasta: StrictStr | None = None
    gtf: StrictStr | None = None
    transcript_fasta: StrictStr | None = None
    salmon_index: StrictStr | None = None
    hisat2_index: StrictStr | None = None
    hisat2_splice_sites: StrictStr | None = None
    root: StrictStr | None = None
    manifest: StrictStr | None = None

    @model_validator(mode="after")
    def validate_source_contract(self) -> "ReferenceConfig":
        if self.source == "local":
            if self.root is None or not self.root.strip():
                raise ValueError("reference.root is required when reference.source is local.")
            if self.manifest is None or not self.manifest.strip():
                raise ValueError("reference.manifest is required when reference.source is local.")
            if any(value is not None for value in (self.genome, self.fasta, self.gtf, self.transcript_fasta, self.salmon_index, self.hisat2_index, self.hisat2_splice_sites)):
                raise ValueError("reference.source local resolves assets only from reference.manifest; do not set genome/fasta/gtf/transcript_fasta/salmon_index/hisat2_index/hisat2_splice_sites in project.yaml.")
        return self


class ProjectConfig(StrictModel):
    schema_version: StrictStr
    project: ProjectInfo
    organism: OrganismConfig
    input: InputConfig
    design: DesignConfig
    metadata_file: StrictStr
    contrasts_file: StrictStr
    upstream: UpstreamConfig
    reference: ReferenceConfig = Field(
        default_factory=lambda: ReferenceConfig(source="igenomes", genome=None)
    )
    thresholds: ThresholdsConfig
    annotation: AnnotationConfig | None = None
    analysis: AnalysisConfig | None = None

    @field_validator("schema_version")
    @classmethod
    def validate_schema_version(cls, value: str) -> str:
        if value not in {LEGACY_SCHEMA_VERSION, "1.1", SUPPORTED_SCHEMA_VERSION}:
            raise ValueError(
                f"Unsupported project schema version: {value}. "
                f"Supported schema versions: {LEGACY_SCHEMA_VERSION}, 1.1, {SUPPORTED_SCHEMA_VERSION}"
            )
        return value

    @model_validator(mode="after")
    def validate_input_upstream_contract(self) -> "ProjectConfig":
        if self.project.preset is Preset.QC and self.input.type is not InputType.FASTQ:
            raise ValueError("QC preset is available only for FASTQ projects.")
        if self.input.type is InputType.FASTQ:
            if self.upstream.engine != "nfcore_rnaseq":
                raise ValueError("FASTQ projects require upstream.engine: nfcore_rnaseq.")
            if self.upstream.pipeline_version != NFCORE_RNASEQ_VERSION:
                raise ValueError(
                    "FASTQ projects require upstream.pipeline_version: "
                    f"{NFCORE_RNASEQ_VERSION}."
                )
            if self.upstream.strandedness is None:
                raise ValueError("FASTQ projects require upstream.strandedness.")
            method = self.upstream.quantification.method if self.upstream.quantification else "salmon"
            if method == "hisat2_featurecounts" and self.upstream.strandedness == "auto":
                raise ValueError(
                    "HISAT2 + featureCounts requires explicit upstream.strandedness: "
                    "unstranded, forward, or reverse; auto inference is not implemented."
                )
        elif self.upstream.engine != "external":
            raise ValueError("raw_counts projects require upstream.engine: external.")
        if self.annotation is not None and self.annotation.organism != self.organism.species.value:
            raise ValueError("annotation.organism must match organism.species.")
        if self.schema_version == SUPPORTED_SCHEMA_VERSION and self.analysis is None:
            raise ValueError(f"schema_version {SUPPORTED_SCHEMA_VERSION} requires an explicit analysis.enrichment list (it may be empty).")
        selected = self.analysis.enrichment if self.analysis is not None else ()
        if selected and self.annotation is None:
            raise ValueError("analysis.enrichment requires an explicit annotation contract.")
        if selected and self.project.preset is not Preset.L2:
            raise ValueError("analysis.enrichment is only supported by L2 projects.")
        return self
