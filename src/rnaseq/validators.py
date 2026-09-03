"""Read-only validation for FASTQ and raw-count RNA-seq projects."""

from __future__ import annotations

import csv
import re
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from enum import Enum
from pathlib import Path
from typing import Iterable

from rnaseq.errors import ProjectConfigError
from rnaseq.models import DesignType, InputType, ProjectConfig, SequencingLayout
from rnaseq.project import LoadedProject, load_project
from rnaseq.references import LocalReference, LocalReferenceError, load_local_reference

CONTRAST_HEADER = ["contrast_id", "factor", "numerator", "denominator"]
FORMULA_PATTERN = re.compile(
    r"^~\s*([A-Za-z_][A-Za-z0-9_.]*)(?:\s*\+\s*([A-Za-z_][A-Za-z0-9_.]*))*\s*$"
)
VARIABLE_PATTERN = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*")
FASTQ_NAME_PATTERN = re.compile(
    r"^(?P<sample>.+?)(?:_L(?P<lane>\d{3}))?_R(?P<read>[12])(?:_\d{3})?\.(?:fastq|fq)\.gz$"
)


class Severity(str, Enum):
    ERROR = "ERROR"
    WARNING = "WARNING"


@dataclass(frozen=True)
class ValidationIssue:
    severity: Severity
    code: str
    message: str


@dataclass(frozen=True)
class CountsSummary:
    path: Path
    sample_ids: tuple[str, ...]
    gene_count: int
    all_zero_genes: int


@dataclass(frozen=True)
class FastqRecord:
    sample_id: str
    lane: str
    fastq_1: Path
    fastq_2: Path | None


@dataclass(frozen=True)
class FastqSummary:
    path: Path
    layout: SequencingLayout
    records: tuple[FastqRecord, ...]

    @property
    def sample_ids(self) -> tuple[str, ...]:
        return tuple(sorted({record.sample_id for record in self.records}))


@dataclass(frozen=True)
class MetadataSummary:
    path: Path
    columns: tuple[str, ...]
    rows: tuple[dict[str, str], ...]
    sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class ContrastDefinition:
    contrast_id: str
    factor: str
    numerator: str
    denominator: str


@dataclass(frozen=True)
class ContrastsSummary:
    path: Path
    contrasts: tuple[ContrastDefinition, ...]


@dataclass
class ValidationReport:
    project_dir: Path
    loaded: LoadedProject | None = None
    formula_variables: tuple[str, ...] = ()
    counts: CountsSummary | None = None
    fastq: FastqSummary | None = None
    local_reference: LocalReference | None = None
    metadata: MetadataSummary | None = None
    contrasts: ContrastsSummary | None = None
    groups: OrderedDict[str, OrderedDict[str, int]] = field(default_factory=OrderedDict)
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def config(self) -> ProjectConfig | None:
        return self.loaded.config if self.loaded else None

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity is Severity.ERROR]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity is Severity.WARNING]

    @property
    def is_valid(self) -> bool:
        return not self.errors

    @property
    def state(self) -> str:
        if self.errors:
            return "FAIL"
        if self.warnings:
            return "PASS WITH WARNINGS"
        return "PASS"

    @property
    def execution_ready(self) -> bool:
        if not self.is_valid or self.config is None:
            return False
        if self.config.input.type is InputType.RAW_COUNTS:
            return True
        reference = self.config.reference
        if reference is None:
            return False
        if reference.source == "igenomes":
            return reference.genome is not None
        if reference.source == "local":
            return self.local_reference is not None and self.local_reference.salmon_index is not None
        return reference.fasta is not None and reference.gtf is not None

    @property
    def execution_blockers(self) -> tuple[str, ...]:
        if self.config is None or self.config.input.type is InputType.RAW_COUNTS:
            return ()
        if self.execution_ready:
            return ()
        if self.config.reference.source == "local":
            if self.local_reference is not None and self.local_reference.salmon_index is None:
                return (f"local Salmon index is not built. Run: rnaseq reference prepare {self.local_reference.root}",)
            return ("local reference configuration or integrity validation failed",)
        return ("reference genome not configured",)

    def error(self, code: str, message: str) -> None:
        self.issues.append(ValidationIssue(Severity.ERROR, code, message))

    def warning(self, code: str, message: str) -> None:
        self.issues.append(ValidationIssue(Severity.WARNING, code, message))


def _open_csv(path: Path):
    return path.open("r", encoding="utf-8", newline="")


def _is_blank_row(row: list[str]) -> bool:
    return not row or all(value == "" for value in row)


def parse_formula(formula: str, report: ValidationReport) -> tuple[str, ...]:
    """Parse the deliberately restricted additive formula grammar."""

    if any(operator in formula for operator in (":", "*", "/", "^")):
        report.error(
            "unsupported_formula",
            "Interaction/factorial formulas are not supported in the current pipeline version.",
        )
        return ()
    if not FORMULA_PATTERN.fullmatch(formula):
        report.error(
            "unsupported_formula",
            "Formula must use only additive variables, for example '~ condition' or "
            "'~ subject_id + condition'.",
        )
        return ()
    variables = tuple(VARIABLE_PATTERN.findall(formula.removeprefix("~")))
    if len(set(variables)) != len(variables):
        report.error("duplicate_formula_variable", "Design formula contains duplicate variables.")
        return ()
    return variables


def _validate_header(
    header: list[str],
    *,
    report: ValidationReport,
    code_prefix: str,
    label: str,
) -> bool:
    valid = True
    blank_positions = [str(index + 1) for index, name in enumerate(header) if name.strip() == ""]
    if blank_positions:
        report.error(
            f"blank_{code_prefix}_column",
            f"{label} contains blank column names at positions: "
            + ", ".join(blank_positions),
        )
        valid = False
    duplicates = sorted(name for name, count in Counter(header).items() if count > 1)
    if duplicates:
        report.error(
            f"duplicate_{code_prefix}_column",
            f"{label} contains duplicate column names: " + ", ".join(duplicates),
        )
        valid = False
    return valid


def validate_counts(path: Path, report: ValidationReport) -> CountsSummary | None:
    if not path.exists():
        report.error("missing_counts", f"Count matrix not found: {path}")
        return None
    if not path.is_file():
        report.error("invalid_counts_path", f"Count matrix is not a file: {path}")
        return None

    with _open_csv(path) as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            report.error("empty_counts", "Count matrix is empty.")
            return None

        header_valid = _validate_header(
            header, report=report, code_prefix="count_matrix", label="Count matrix"
        )
        if len(header) < 2:
            report.error(
                "missing_count_samples",
                "Count matrix must contain gene_id and at least one sample column.",
            )
            return None
        if header[0] != "gene_id":
            report.error("invalid_gene_column", "First count-matrix column must be named gene_id.")
            header_valid = False

        sample_ids = header[1:]
        seen_genes: set[str] = set()
        seen_records: set[tuple[str, ...]] = set()
        gene_count = 0
        all_zero_genes = 0

        for line_number, row in enumerate(reader, start=2):
            if _is_blank_row(row):
                continue
            gene_count += 1
            if len(row) != len(header):
                report.error(
                    "malformed_count_row",
                    f"Count matrix row {line_number} has {len(row)} fields; expected {len(header)}.",
                )
                continue

            record = tuple(row)
            if record in seen_records:
                report.error(
                    "duplicate_count_row", f"Count matrix contains a duplicate row at line {line_number}."
                )
            seen_records.add(record)

            gene_id = row[0]
            if gene_id.strip() == "":
                report.error("blank_gene_id", f"Count matrix has a blank gene_id at line {line_number}.")
            elif gene_id in seen_genes:
                report.error("duplicate_gene_id", f"Duplicate gene_id at line {line_number}: {gene_id}")
            seen_genes.add(gene_id)

            parsed: list[Decimal] = []
            row_valid = True
            for column_number, raw_value in enumerate(row[1:], start=2):
                try:
                    value = Decimal(raw_value)
                except InvalidOperation:
                    report.error(
                        "non_numeric_count",
                        f"Non-numeric count at row {line_number}, column {column_number}: {raw_value!r}",
                    )
                    row_valid = False
                    continue
                if not value.is_finite():
                    report.error(
                        "non_finite_count",
                        f"Non-finite count at row {line_number}, column {column_number}: {raw_value!r}",
                    )
                    row_valid = False
                elif value < 0:
                    report.error(
                        "negative_count",
                        f"Negative count at row {line_number}, column {column_number}: {raw_value!r}",
                    )
                    row_valid = False
                elif value != value.to_integral_value():
                    report.error(
                        "fractional_count",
                        f"Non-integer count at row {line_number}, column {column_number}: {raw_value!r}",
                    )
                    row_valid = False
                parsed.append(value)
            if row_valid and parsed and all(value == 0 for value in parsed):
                all_zero_genes += 1

    if gene_count == 0:
        report.error("missing_genes", "Count matrix must contain at least one gene row.")
    if all_zero_genes:
        report.warning(
            "all_zero_genes",
            f"Count matrix contains {all_zero_genes} all-zero gene(s); no genes were filtered.",
        )
    if not header_valid:
        return None
    return CountsSummary(path, tuple(sample_ids), gene_count, all_zero_genes)


def validate_fastq(
    path: Path, layout: SequencingLayout, report: ValidationReport
) -> FastqSummary | None:
    """Discover conventional gzipped FASTQs without reading their biological content."""

    if not path.exists():
        report.error("missing_fastq_directory", f"FASTQ input directory not found: {path}")
        return None
    if not path.is_dir():
        report.error("invalid_fastq_path", f"FASTQ input path is not a directory: {path}")
        return None

    # AppleDouble sidecars are filesystem metadata, not FASTQ inputs.
    files = sorted(item for item in path.iterdir() if item.is_file() and not item.name.startswith("._"))
    if not files:
        report.error("missing_fastq_files", "FASTQ input directory contains no files.")
        return None
    parsed: dict[tuple[str, str], dict[str, Path]] = {}
    for file_path in files:
        match = FASTQ_NAME_PATTERN.fullmatch(file_path.name)
        if match is None:
            report.error(
                "invalid_fastq_filename",
                "FASTQ files must use .fastq.gz or .fq.gz and an unambiguous "
                f"_R1/_R2 name marker: {file_path.name}",
            )
            continue
        sample_id = match.group("sample")
        lane = match.group("lane") or "L000"
        read = match.group("read")
        key = (sample_id, lane)
        reads = parsed.setdefault(key, {})
        if read in reads:
            report.error(
                "duplicate_fastq_assignment",
                f"Multiple R{read} FASTQ files assign to sample {sample_id!r}, lane {lane}.",
            )
            continue
        reads[read] = file_path

    records: list[FastqRecord] = []
    for (sample_id, lane), reads in sorted(parsed.items()):
        r1, r2 = reads.get("1"), reads.get("2")
        if layout is SequencingLayout.PAIRED_END:
            if r1 is None:
                report.error("orphan_r2", f"Sample {sample_id!r}, lane {lane} has R2 but no R1.")
            if r2 is None:
                report.error("orphan_r1", f"Sample {sample_id!r}, lane {lane} has R1 but no R2.")
            if r1 is not None and r2 is not None:
                records.append(FastqRecord(sample_id, lane, r1, r2))
        else:
            if r1 is None:
                report.error("orphan_r2", f"Sample {sample_id!r}, lane {lane} has R2 but no R1.")
            elif r2 is not None:
                report.error(
                    "unexpected_r2", f"Single-end sample {sample_id!r}, lane {lane} includes R2."
                )
            else:
                records.append(FastqRecord(sample_id, lane, r1, None))
    if not records and not report.errors:
        report.error("missing_fastq_files", "No valid FASTQ assignments were discovered.")
    return FastqSummary(path, layout, tuple(records))


def validate_metadata(
    path: Path,
    required_variables: Iterable[str],
    report: ValidationReport,
) -> MetadataSummary | None:
    if not path.exists():
        report.error("missing_metadata", f"Metadata file not found: {path}")
        return None
    if not path.is_file():
        report.error("invalid_metadata_path", f"Metadata path is not a file: {path}")
        return None

    with _open_csv(path) as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            report.error("empty_metadata", "Metadata file is empty.")
            return None

        header_valid = _validate_header(
            header, report=report, code_prefix="metadata", label="Metadata"
        )
        if "sample_id" not in header:
            report.error("missing_sample_id_column", "Metadata must contain a sample_id column.")
            header_valid = False
        missing_variables = sorted(set(required_variables) - set(header))
        if missing_variables:
            report.error(
                "missing_design_variables",
                "Metadata is missing design variable(s): " + ", ".join(missing_variables),
            )
            header_valid = False
        if not header_valid:
            return None

        rows: list[dict[str, str]] = []
        sample_ids: list[str] = []
        seen_samples: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            if _is_blank_row(row):
                continue
            if len(row) != len(header):
                report.error(
                    "malformed_metadata_row",
                    f"Metadata row {line_number} has {len(row)} fields; expected {len(header)}.",
                )
                continue
            record = dict(zip(header, row, strict=True))
            sample_id = record["sample_id"]
            sample_ids.append(sample_id)
            if sample_id.strip() == "":
                report.error("blank_sample_id", f"Metadata has a blank sample_id at line {line_number}.")
            elif sample_id in seen_samples:
                report.error(
                    "duplicate_sample_id", f"Duplicate metadata sample_id at line {line_number}: {sample_id}"
                )
            seen_samples.add(sample_id)
            for variable in required_variables:
                if record[variable].strip() == "":
                    report.error(
                        "missing_design_value",
                        f"Metadata value for required design variable {variable!r} is missing "
                        f"at line {line_number}.",
                    )
            rows.append(record)

    if not rows:
        report.error("missing_metadata_rows", "Metadata must contain at least one sample row.")
    return MetadataSummary(path, tuple(header), tuple(rows), tuple(sample_ids))


def validate_contrasts(
    path: Path,
    metadata: MetadataSummary | None,
    formula_variables: tuple[str, ...],
    report: ValidationReport,
) -> ContrastsSummary | None:
    if not path.exists():
        report.error("missing_contrasts", f"Contrasts file not found: {path}")
        return None
    if not path.is_file():
        report.error("invalid_contrasts_path", f"Contrasts path is not a file: {path}")
        return None

    with _open_csv(path) as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            report.error("empty_contrasts", "Contrasts file is empty.")
            return None
        if header != CONTRAST_HEADER:
            report.error(
                "invalid_contrast_schema",
                "Contrasts columns must be exactly, in order: " + ",".join(CONTRAST_HEADER),
            )
            return None

        seen_ids: set[str] = set()
        definitions: list[ContrastDefinition] = []
        for line_number, row in enumerate(reader, start=2):
            if _is_blank_row(row):
                continue
            if len(row) != len(CONTRAST_HEADER):
                report.error(
                    "malformed_contrast_row",
                    f"Contrasts row {line_number} has {len(row)} fields; expected 4.",
                )
                continue
            contrast_id, factor, numerator, denominator = row
            if any(value.strip() == "" for value in row):
                report.error(
                    "blank_contrast_value", f"Contrasts row {line_number} contains a blank required value."
                )
                continue
            if contrast_id in seen_ids:
                report.error(
                    "duplicate_contrast_id", f"Duplicate contrast_id at line {line_number}: {contrast_id}"
                )
            seen_ids.add(contrast_id)
            if numerator == denominator:
                report.error(
                    "identical_contrast_levels",
                    f"Contrast {contrast_id} uses the same numerator and denominator: {numerator}",
                )
            if factor not in formula_variables:
                report.error(
                    "contrast_factor_not_in_design",
                    f"Contrast factor {factor!r} is not present in the design formula.",
                )
            if metadata is not None:
                if factor not in metadata.columns:
                    report.error(
                        "missing_contrast_factor",
                        f"Contrast factor {factor!r} does not exist in metadata.",
                    )
                else:
                    levels = {record[factor] for record in metadata.rows}
                    if numerator not in levels:
                        report.error(
                            "missing_contrast_numerator",
                            f"Contrast {contrast_id} numerator level not found in {factor}: {numerator}",
                        )
                    if denominator not in levels:
                        report.error(
                            "missing_contrast_denominator",
                            f"Contrast {contrast_id} denominator level not found in {factor}: {denominator}",
                        )
            definitions.append(
                ContrastDefinition(contrast_id, factor, numerator, denominator)
            )

    if not definitions:
        report.error("missing_contrast_rows", "Contrasts file must contain at least one contrast row.")
    return ContrastsSummary(path, tuple(definitions))


def _validate_sample_agreement(
    counts: CountsSummary,
    metadata: MetadataSummary,
    report: ValidationReport,
) -> None:
    count_set = set(counts.sample_ids)
    metadata_set = set(metadata.sample_ids)
    count_only = sorted(count_set - metadata_set)
    metadata_only = sorted(metadata_set - count_set)
    if count_only:
        report.error(
            "count_only_samples",
            "Samples present in count matrix but absent from metadata: " + ", ".join(count_only),
        )
    if metadata_only:
        report.error(
            "metadata_only_samples",
            "Samples present in metadata but absent from count matrix: " + ", ".join(metadata_only),
        )


def _validate_fastq_sample_agreement(
    fastq: FastqSummary, metadata: MetadataSummary, report: ValidationReport
) -> None:
    fastq_set = set(fastq.sample_ids)
    metadata_set = set(metadata.sample_ids)
    fastq_only = sorted(fastq_set - metadata_set)
    metadata_only = sorted(metadata_set - fastq_set)
    if fastq_only:
        report.error(
            "fastq_only_samples",
            "Samples present in FASTQ input but absent from metadata: " + ", ".join(fastq_only),
        )
    if metadata_only:
        report.error(
            "metadata_only_fastq_samples",
            "Samples present in metadata but absent from FASTQ input: " + ", ".join(metadata_only),
        )


def _build_groups_and_validate_design(report: ValidationReport) -> None:
    if report.metadata is None or report.contrasts is None or report.config is None:
        return
    metadata = report.metadata
    design_type = report.config.design.type
    factors = list(dict.fromkeys(contrast.factor for contrast in report.contrasts.contrasts))

    for factor in factors:
        if factor not in metadata.columns:
            continue
        counts = Counter(record[factor] for record in metadata.rows)
        ordered = OrderedDict((level, counts[level]) for level in sorted(counts))
        report.groups[factor] = ordered
        for level, count in ordered.items():
            if count in (1, 2):
                report.warning(
                    "limited_replication",
                    f"Factor {factor}, level {level} has n={count}; biological replication is limited.",
                )

        level_count = len(ordered)
        if design_type is DesignType.TWO_GROUP and level_count != 2:
            report.error(
                "invalid_two_group_design",
                f"two_group design requires exactly 2 levels for {factor}; found {level_count}.",
            )
        elif design_type is DesignType.MULTI_GROUP and level_count < 3:
            report.error(
                "invalid_multi_group_design",
                f"multi_group design requires at least 3 levels for {factor}; found {level_count}.",
            )
        elif design_type is DesignType.PAIRED and level_count != 2:
            report.error(
                "invalid_paired_levels",
                f"paired design requires exactly 2 levels for {factor}; found {level_count}.",
            )

    if design_type is not DesignType.PAIRED:
        return
    if "subject_id" not in report.formula_variables:
        report.error(
            "missing_pairing_variable", "Paired design formula must include subject_id."
        )
        return
    if "subject_id" not in metadata.columns:
        return

    for factor in factors:
        if factor not in metadata.columns:
            continue
        expected_levels = set(record[factor] for record in metadata.rows)
        by_subject: dict[str, list[str]] = {}
        for record in metadata.rows:
            by_subject.setdefault(record["subject_id"], []).append(record[factor])
        for subject in sorted(by_subject):
            observed = by_subject[subject]
            if len(observed) != 2 or set(observed) != expected_levels or len(set(observed)) != 2:
                report.error(
                    "incomplete_pair",
                    f"Subject {subject!r} must have exactly one sample from each {factor} level; "
                    f"observed: {', '.join(observed)}",
                )


def validate_project(project_dir: Path | str) -> ValidationReport:
    """Validate a project without modifying any project input."""

    report = ValidationReport(project_dir=Path(project_dir).resolve())
    try:
        report.loaded = load_project(project_dir)
    except ProjectConfigError as exc:
        report.error("invalid_project_config", str(exc))
        return report

    report.formula_variables = parse_formula(report.config.design.formula, report)
    if report.config.input.type is InputType.RAW_COUNTS:
        report.counts = validate_counts(report.loaded.input_path, report)
    else:
        assert report.config.input.layout is not None
        report.fastq = validate_fastq(
            report.loaded.input_path, report.config.input.layout, report
        )
        if report.config.reference.source == "local":
            try:
                report.local_reference = load_local_reference(
                    report.config.reference, report.config.organism.species.value
                )
            except LocalReferenceError as exc:
                report.error("invalid_local_reference", str(exc))
    report.metadata = validate_metadata(
        report.loaded.metadata_path, report.formula_variables, report
    )
    report.contrasts = validate_contrasts(
        report.loaded.contrasts_path,
        report.metadata,
        report.formula_variables,
        report,
    )
    if report.counts is not None and report.metadata is not None:
        _validate_sample_agreement(report.counts, report.metadata, report)
    if report.fastq is not None and report.metadata is not None:
        _validate_fastq_sample_agreement(report.fastq, report.metadata, report)
    _build_groups_and_validate_design(report)
    return report
