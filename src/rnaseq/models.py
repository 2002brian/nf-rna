"""Strict Pydantic models for the Milestone 1 project contract."""

from __future__ import annotations

import re
from enum import Enum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, StrictStr, field_validator, model_validator

from rnaseq import __version__

SUPPORTED_SCHEMA_VERSION = "1.3"
LEGACY_SCHEMA_VERSION = "1.0"
PIPELINE_VERSION = __version__
OFFICIAL_EXECUTION_IMAGE_REPOSITORY = "ghcr.io/2002brian/nf-rna"


def execution_image_for_version(version: str) -> str:
    """Return the official execution image paired with one CLI version.

    The package version is the sole release identity for the control plane.
    Keeping this mapping here makes stable and prerelease builds use the same
    explicit tag (for example, ``1.0.1rc1``) rather than falling back to a
    mutable convenience tag.
    """

    return f"{OFFICIAL_EXECUTION_IMAGE_REPOSITORY}:{version}"


DEFAULT_EXECUTION_IMAGE = execution_image_for_version(PIPELINE_VERSION)
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
    PAIRED_TWO_GROUP = "paired_two_group"


class MetadataVariableType(str, Enum):
    """Statistical representation of a design variable in DESeq2."""

    CATEGORICAL = "categorical"
    CONTINUOUS = "continuous"


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
    variables: dict[StrictStr, MetadataVariableType] | None = None
    pair_id: StrictStr | None = None

    @model_validator(mode="before")
    @classmethod
    def normalize_interim_paired_contract(cls, value: object) -> object:
        """Read projects written by the short-lived pre-M5A paired spelling."""

        if not isinstance(value, dict):
            return value
        normalized = dict(value)
        if normalized.get("type") == "paired":
            normalized["type"] = DesignType.PAIRED_TWO_GROUP.value
        if "pairing_column" in normalized:
            if "pair_id" in normalized and normalized["pair_id"] != normalized["pairing_column"]:
                raise ValueError("design.pair_id and legacy design.pairing_column disagree.")
            normalized["pair_id"] = normalized.pop("pairing_column")
        return normalized

    @model_validator(mode="after")
    def validate_pairing_contract(self) -> "DesignConfig":
        if self.variables is not None:
            for name in self.variables:
                if not name.strip():
                    raise ValueError("design.variables cannot contain a blank variable name.")
        if self.type is DesignType.PAIRED_TWO_GROUP:
            if self.pair_id is None or not self.pair_id.strip():
                raise ValueError(
                    "design.pair_id is required for a paired_two_group biological design; "
                    "it is independent of paired-end sequencing layout."
                )
            if self.variables is not None and self.variables.get(self.pair_id) not in {None, MetadataVariableType.CATEGORICAL}:
                raise ValueError("design.pair_id must be declared categorical when listed in design.variables.")
        elif self.pair_id is not None:
            raise ValueError("design.pair_id is supported only when design.type is paired_two_group.")
        return self


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
    """Configured nf-rna significance thresholds for KEGG preranked GSEA.

    clusterProfiler calculation intentionally uses ``pvalueCutoff=1`` so
    nf-rna can retain every evaluated term before applying these thresholds.
    """

    minimum_ranked_genes: int = Field(default=50, ge=1)
    pvalue_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
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
    """Identifier and annotation-QC contract for optional preranked GSEA."""

    organism: Literal["Homo sapiens", "Mus musculus"]
    input_id_type: Literal["ENSEMBL", "ENTREZID", "SYMBOL"]
    target_id_type: Literal["ENTREZID"] = "ENTREZID"
    gene_symbol_output: bool = True
    # ``minimum_mapping_rate`` was the only (blocking) guardrail in schema
    # 1.0/1.1 projects.  It now deliberately names the fail-closed threshold;
    # a distinct warning threshold permits an auditable, non-blocking warning.
    mapping_warning_rate: float = Field(default=0.70, ge=0.0, le=1.0)
    minimum_mapping_rate: float = Field(default=0.50, ge=0.0, le=1.0)
    minimum_mapped_foreground: int = Field(default=5, ge=1)
    enrichment: EnrichmentConfig = Field(default_factory=EnrichmentConfig)

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_mapping_threshold(cls, value: object) -> object:
        """Treat a lone legacy threshold as the old warning-level contract.

        Existing projects expressed only one threshold. Preserve it as their
        warning threshold while applying the documented 0.50 blocking default
        where that does not make a historically permissive project stricter.
        New configurations can set both fields explicitly.
        """

        if not isinstance(value, dict):
            return value
        migrated = dict(value)
        if "mapping_warning_rate" not in migrated and "minimum_mapping_rate" in migrated:
            legacy_threshold = migrated["minimum_mapping_rate"]
            migrated["mapping_warning_rate"] = legacy_threshold
            try:
                # Do not make a historically more-permissive project stricter
                # merely by loading it under the dual-threshold contract.
                migrated["minimum_mapping_rate"] = min(0.50, float(legacy_threshold))
            except (TypeError, ValueError):
                migrated["minimum_mapping_rate"] = 0.50
        return migrated

    @model_validator(mode="after")
    def validate_mapping_thresholds(self) -> "AnnotationConfig":
        if self.minimum_mapping_rate > self.mapping_warning_rate:
            raise ValueError(
                "annotation.minimum_mapping_rate must be <= annotation.mapping_warning_rate."
            )
        return self


PUBLIC_ENRICHMENT_METHODS = ("go", "kegg", "gsea")
# Retained for callers that used the former single-method public constant.
PUBLIC_ENRICHMENT_METHOD = "gsea"
LEGACY_ENRICHMENT_MODULES = frozenset({"go", "gsea-go", "kegg", "gsea-kegg"})
INTERNAL_ENRICHMENT_BACKENDS = {
    "go": ("go",),
    "kegg": ("kegg",),
    "gsea": ("gsea-go", "gsea-kegg"),
}


def normalize_enrichment_selection(value: object) -> tuple[str, ...]:
    """Normalize the public enrichment selection and the former complete scope.

    GO ORA, KEGG ORA, and preranked GSEA are independently selectable public
    methods. Public selections are deduplicated into their canonical order so
    a repeated selection cannot schedule a backend twice. The historical
    four-backend spelling remains a lossless alias.
    """

    if isinstance(value, str) and value in PUBLIC_ENRICHMENT_METHODS:
        return (value,)
    if isinstance(value, (list, tuple)):
        selected = tuple(value)
        if not selected:
            return ()
        if all(isinstance(item, str) and item in PUBLIC_ENRICHMENT_METHODS for item in selected):
            return tuple(method for method in PUBLIC_ENRICHMENT_METHODS if method in selected)
        if len(selected) == len(LEGACY_ENRICHMENT_MODULES) and set(selected) == LEGACY_ENRICHMENT_MODULES:
            return PUBLIC_ENRICHMENT_METHODS
        if any(item in LEGACY_ENRICHMENT_MODULES for item in selected):
            raise ValueError(
                "Legacy analysis.enrichment is accepted only when it contains exactly "
                "go, gsea-go, kegg, gsea-kegg; use public methods go, kegg, and/or gsea."
            )
    return value  # Pydantic reports invalid types or unsupported public values.


def production_enrichment_backends(value: object) -> tuple[str, ...]:
    """Expand selected public methods into their stable runtime backends."""

    selected = normalize_enrichment_selection(value)
    if selected == ():
        return ()
    if all(item in INTERNAL_ENRICHMENT_BACKENDS for item in selected):
        return tuple(backend for item in selected for backend in INTERNAL_ENRICHMENT_BACKENDS[item])
    raise ValueError("analysis.enrichment must contain only go, kegg, and/or gsea.")


class AnalysisConfig(StrictModel):
    """Execution choices, separated from statistical thresholds.

    An explicit list makes the delivery contract reproducible: a later run never
    gains an enrichment step merely because the CLI happened to use a flag.
    """

    enrichment: tuple[Literal["go", "kegg", "gsea"], ...] = ()

    @field_validator("enrichment", mode="before")
    @classmethod
    def normalize_legacy_enrichment(cls, value: object) -> object:
        return normalize_enrichment_selection(value)



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
    acceptance: Literal["standard", "production"] = "standard"

    @model_validator(mode="after")
    def validate_source_contract(self) -> "ReferenceConfig":
        if self.acceptance == "production" and self.source != "local":
            raise ValueError(
                "reference.acceptance: production requires reference.source: local "
                "with a checksum-bound managed manifest."
            )
        if self.source == "local":
            if self.root is None or not self.root.strip():
                raise ValueError("reference.root is required when reference.source is local.")
            if self.manifest is None or not self.manifest.strip():
                raise ValueError("reference.manifest is required when reference.source is local.")
            if any(value is not None for value in (self.genome, self.fasta, self.gtf, self.transcript_fasta, self.salmon_index, self.hisat2_index, self.hisat2_splice_sites)):
                raise ValueError("reference.source local resolves assets only from reference.manifest; do not set genome/fasta/gtf/transcript_fasta/salmon_index/hisat2_index/hisat2_splice_sites in project.yaml.")
        return self


class RuntimeConfig(StrictModel):
    """Requested first-party execution image for downstream Nextflow tasks."""

    execution_image: StrictStr = DEFAULT_EXECUTION_IMAGE

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_control_plane_image(cls, value: object) -> object:
        """Read the pre-1.0 runtime spelling without serializing it again."""

        if not isinstance(value, dict) or "control_plane_image" not in value:
            return value
        if "execution_image" in value:
            raise ValueError(
                "runtime may specify only execution_image; control_plane_image is a legacy read-only alias."
            )
        migrated = dict(value)
        migrated["execution_image"] = migrated.pop("control_plane_image")
        return migrated

    @field_validator("execution_image")
    @classmethod
    def validate_image_reference(cls, value: str) -> str:
        if not value.strip() or any(char.isspace() for char in value):
            raise ValueError("runtime.execution_image must be a non-blank container reference without whitespace.")
        return value


class ExecutionConfig(StrictModel):
    """Portable aggregate local execution ceiling, not a scientific setting."""

    profile: Literal["local"] = "local"
    # ``auto`` sizes the ceiling from this machine at run time (capacity minus
    # an OS reserve); an explicit integer is a deliberate user limit.
    max_cpus: Annotated[int, Field(ge=1)] | Literal["auto"] = "auto"
    max_memory_gb: Annotated[int, Field(ge=1)] | Literal["auto"] = "auto"


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
    runtime: RuntimeConfig = Field(default_factory=RuntimeConfig)
    execution: ExecutionConfig = Field(default_factory=ExecutionConfig)
    thresholds: ThresholdsConfig
    annotation: AnnotationConfig | None = None
    analysis: AnalysisConfig | None = None

    @field_validator("schema_version")
    @classmethod
    def validate_schema_version(cls, value: str) -> str:
        if value not in {LEGACY_SCHEMA_VERSION, "1.1", "1.2", SUPPORTED_SCHEMA_VERSION}:
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
        if self.design.type is DesignType.PAIRED_TWO_GROUP:
            assert self.design.pair_id is not None
            variables = tuple(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", self.design.formula.removeprefix("~")))
            if self.design.pair_id not in variables:
                raise ValueError("design.pair_id must be present in design.formula.")
        if self.schema_version == SUPPORTED_SCHEMA_VERSION:
            formula_variables = tuple(
                re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", self.design.formula.removeprefix("~"))
            )
            declared = self.design.variables
            if declared is None:
                raise ValueError(
                    "schema_version 1.3 requires design.variables to declare the type of every design formula variable."
                )
            missing = sorted(set(formula_variables) - set(declared))
            if missing:
                raise ValueError(
                    "schema_version 1.3 design.variables is missing formula variable(s): "
                    + ", ".join(missing)
                )
        if self.reference.acceptance == "production":
            if self.input.type is not InputType.FASTQ:
                raise ValueError("reference.acceptance: production is supported only for FASTQ projects.")
            image = self.runtime.execution_image
            if image.endswith(":latest") or (":" not in image.rsplit("/", 1)[-1] and "@sha256:" not in image):
                raise ValueError(
                    "Production acceptance requires an immutable runtime.execution_image: "
                    "use an image digest or a versioned tag, never latest or an untagged reference."
                )
        if self.schema_version in {"1.2", SUPPORTED_SCHEMA_VERSION} and self.analysis is None:
            raise ValueError(f"schema_version {self.schema_version} requires an explicit analysis.enrichment list (it may be empty).")
        selected = self.analysis.enrichment if self.analysis is not None else ()
        if selected and self.annotation is None:
            raise ValueError("analysis.enrichment requires an explicit annotation contract.")
        if selected and self.project.preset is not Preset.L2:
            raise ValueError("analysis.enrichment is only supported by L2 projects.")
        return self
