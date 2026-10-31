"""Server-side dataset browsing for the Datasets page (no Streamlit here).

A *dataset* is any folder inside the allowed roots. Scanning only lists file
names, sizes and permissions and reads the first line of small tables; FASTQ
contents are never read. Sample recognition reuses the project validator
(FASTQ naming rule, count-matrix checks) and the samplesheet reader of
``rnaseq new --fastq-samplesheet``. Nothing here writes anywhere.

Paths are resolved (symlinks included) and must stay inside an available
allowed root: a symlink that points elsewhere is reported and never followed,
and symlinked folders are not descended into.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from rnaseq.gui.config import GuiConfig, PathNotAllowed
from rnaseq.models import SequencingLayout
from rnaseq.project import FASTQ_SAMPLESHEET_HEADER
from rnaseq.validators import CountsSummary, FastqSummary, ValidationIssue, ValidationReport, validate_counts, validate_fastq

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
TABLE_SUFFIXES = (".csv", ".tsv", ".txt")
# Folders that hold run outputs or engine state, not input data.
SKIPPED_DIRS = {"runs", "work", "downstream", "delivery", "upstream", "frozen", "planning", "provenance", ".nextflow"}
MAX_DEPTH = 3
MAX_FILES = 5000
HEADER_BYTES = 65536

KIND_LABELS = {
    "fastq": "FASTQ",
    "count_matrix": "count matrix",
    "metadata": "sample metadata",
    "samplesheet": "FASTQ samplesheet",
    "contrasts": "contrasts",
    "table": "other table",
    "other": "other",
}


@dataclass(frozen=True)
class DatasetFile:
    path: Path
    relative: str
    kind: str
    size: int | None
    readable: bool
    symlink: bool = False


@dataclass
class DatasetScan:
    root: Path
    files: list[DatasetFile] = field(default_factory=list)
    subfolders: list[str] = field(default_factory=list)
    outside_links: list[str] = field(default_factory=list)
    unlisted: list[str] = field(default_factory=list)
    truncated: bool = False

    def of_kind(self, kind: str) -> list[DatasetFile]:
        return [item for item in self.files if item.kind == kind]

    @property
    def total_bytes(self) -> int:
        return sum(item.size or 0 for item in self.files)

    @property
    def unreadable(self) -> list[DatasetFile]:
        return [item for item in self.files if not item.readable]

    def kind_summary(self) -> list[dict[str, object]]:
        rows = []
        for kind, label in KIND_LABELS.items():
            members = self.of_kind(kind)
            if members:
                rows.append({"kind": label, "files": len(members), "size": human_size(sum(item.size or 0 for item in members)),
                             "unreadable": sum(1 for item in members if not item.readable)})
        return rows

    @property
    def fastq_folders(self) -> list[Path]:
        """Folders that directly contain FASTQ files (each is one recognisable FASTQ set)."""

        return sorted({item.path.parent for item in self.of_kind("fastq")})


def human_size(size: int | None) -> str:
    if size is None:
        return "?"
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def _first_line(path: Path) -> str | None:
    try:
        with path.open("rb") as handle:
            data = handle.read(HEADER_BYTES)
    except OSError:
        return None
    text = data.decode("utf-8-sig", errors="replace")
    return text.splitlines()[0] if text else ""


def _second_line(path: Path) -> str | None:
    try:
        with path.open("rb") as handle:
            data = handle.read(HEADER_BYTES)
    except OSError:
        return None
    lines = data.decode("utf-8-sig", errors="replace").splitlines()
    return lines[1] if len(lines) > 1 else None


def _cells(line: str) -> list[str]:
    delimiter = "\t" if "\t" in line else ","
    return [cell.strip().strip('"') for cell in line.split(delimiter)]


def classify(path: Path) -> str:
    """The kind of a file, from its name and (for tables) its header line only."""

    name = path.name.lower()
    if name.endswith(FASTQ_SUFFIXES):
        return "fastq"
    if not name.endswith(TABLE_SUFFIXES):
        return "other"
    header_line = _first_line(path)
    if not header_line:
        return "table"
    header = _cells(header_line)
    if tuple(header) == FASTQ_SAMPLESHEET_HEADER or "fastq_1" in header:
        return "samplesheet"
    if "contrast_id" in header:
        return "contrasts"
    if header and header[0] == "sample_id":
        return "metadata"
    second = _second_line(path)
    if len(header) >= 2 and second:
        values = _cells(second)[1:]
        if values and all(value.isdigit() for value in values):
            return "count_matrix"
    return "table"


def scan_dataset(config: GuiConfig, folder: Path | str) -> DatasetScan:
    """List a folder (up to MAX_DEPTH levels, MAX_FILES files) inside the allowed roots."""

    root = config.resolve_allowed(folder)
    if not root.is_dir():
        raise PathNotAllowed(f"Not a folder: {root}")
    scan = DatasetScan(root)
    pending: list[tuple[Path, int]] = [(root, 0)]
    while pending:
        directory, depth = pending.pop(0)
        try:
            entries = sorted(directory.iterdir(), key=lambda item: item.name.lower())
        except OSError as exc:
            scan.unlisted.append(f"{directory}: {exc.strerror or exc}")
            continue
        for entry in entries:
            if entry.name.startswith("."):
                continue
            relative = entry.relative_to(root).as_posix()
            is_link = entry.is_symlink()
            if is_link:
                try:
                    config.resolve_allowed(entry)
                except PathNotAllowed:
                    scan.outside_links.append(relative)
                    continue
            if entry.is_dir():
                if depth == 0:
                    scan.subfolders.append(entry.name)
                if not is_link and entry.name not in SKIPPED_DIRS and depth + 1 < MAX_DEPTH:
                    pending.append((entry, depth + 1))
                continue
            if not entry.is_file():
                continue
            if len(scan.files) >= MAX_FILES:
                scan.truncated = True
                return scan
            try:
                size: int | None = entry.stat().st_size
            except OSError:
                size = None
            readable = os.access(entry, os.R_OK)
            scan.files.append(DatasetFile(entry, relative, classify(entry) if readable else "other", size, readable, is_link))
    return scan


# --------------------------------------------------------------------- recognition (existing rules)

def guess_layout(files: list[DatasetFile]) -> SequencingLayout:
    names = [item.path.name for item in files]
    return SequencingLayout.PAIRED_END if any("_R2" in name for name in names) else SequencingLayout.SINGLE_END


def recognise_fastq(folder: Path, layout: SequencingLayout) -> tuple[FastqSummary | None, list[ValidationIssue]]:
    """Samples and R1/R2 pairs by the validator's FASTQ naming rule (file names only)."""

    report = ValidationReport(project_dir=folder)
    summary = validate_fastq(folder, layout, report)
    # A dataset folder may also hold tables (metadata, samplesheets); only
    # FASTQ-named files are copied by the importer, so naming findings about
    # other files are not FASTQ problems. Misnamed FASTQ files stay errors.
    errors = [
        issue for issue in report.errors
        if not (issue.code == "invalid_fastq_filename" and not issue.message.rsplit(": ", 1)[-1].lower().endswith(FASTQ_SUFFIXES))
    ]
    return summary, errors


def recognise_counts(path: Path) -> tuple[CountsSummary | None, list[ValidationIssue]]:
    """The project validator's count-matrix checks (comma-separated, gene_id first, integer counts)."""

    report = ValidationReport(project_dir=path.parent)
    summary = validate_counts(path, report)
    return summary, report.issues


def table_preview(path: Path, rows: int = 20) -> tuple[list[str], list[list[str]]]:
    """Header and the first rows of a small table, as text."""

    try:
        with path.open(encoding="utf-8-sig", errors="replace") as handle:
            lines = [line.rstrip("\n\r") for _, line in zip(range(rows + 1), handle)]
    except OSError:
        return [], []
    if not lines:
        return [], []
    header = _cells(lines[0])
    body = [_cells(line) for line in lines[1:] if line.strip()]
    width = len(header)
    return header, [(row + [""] * width)[:width] for row in body]


@dataclass(frozen=True)
class CopyPlan:
    """What importing a FASTQ set copies, and whether the target filesystem has room."""

    files: int
    bytes: int
    target: Path
    free_bytes: int | None

    @property
    def fits(self) -> bool:
        return self.free_bytes is None or self.bytes < self.free_bytes


def fastq_copy_plan(records, target: Path) -> CopyPlan:
    """``records``: FastqRecord or ImportedFastq values (anything with fastq_1/fastq_2)."""

    paths = [path for record in records for path in (record.fastq_1, record.fastq_2) if path is not None]
    total = 0
    for path in paths:
        try:
            total += path.stat().st_size
        except OSError:
            pass
    probe = target
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    try:
        free = shutil.disk_usage(probe).free
    except OSError:
        free = None
    return CopyPlan(len(paths), total, target, free)


@dataclass(frozen=True)
class SamplesheetView:
    records: tuple
    layout: SequencingLayout | None
    strandedness: str | None
    error: str | None
    outside: tuple[str, ...] = ()


def recognise_samplesheet(config: GuiConfig, path: Path) -> SamplesheetView:
    """The rows ``rnaseq new --fastq-samplesheet`` would accept, and any FASTQ outside the allowed roots."""

    from rnaseq.errors import ProjectCreationError
    from rnaseq.project import _read_fastq_samplesheet

    try:
        records, layout, strandedness = _read_fastq_samplesheet(path)
    except ProjectCreationError as exc:
        return SamplesheetView((), None, None, str(exc))
    paths = [item for record in records for item in (record.fastq_1, record.fastq_2) if item is not None]
    return SamplesheetView(records, layout, strandedness, None, tuple(outside_allowed(config, paths)))


def outside_allowed(config: GuiConfig, paths: list[Path]) -> list[str]:
    """Paths that resolve (symlinks included) outside the available allowed roots."""

    found = []
    for path in paths:
        try:
            config.resolve_allowed(path)
        except PathNotAllowed:
            found.append(str(path))
    return found


def fastq_paths(records) -> list[Path]:
    return [path for record in records for path in (record.fastq_1, record.fastq_2) if path is not None]
