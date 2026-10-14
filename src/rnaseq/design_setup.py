"""Template-first L1 -> L2 transition: derive and confirm contrasts from completed metadata.

``rnaseq new`` writes only templates.  After L1 QC, ``rnaseq design`` reads the
completed ``metadata.csv`` and proposes the contrasts for the project's design.
Nothing here guesses which comparisons are scientifically meaningful beyond the
single documented rule below, and nothing is written until the caller confirms.

Contrast rule (deterministic, order-preserving, never alphabetical):

* levels are the ``condition`` values in order of first appearance in
  ``metadata.csv`` (or the order frozen in ``design.condition_order``);
* every pair of levels yields one contrast in which the *later* level is the
  numerator and the *earlier* level the denominator, so ``log2FC > 0`` means
  higher in the later level;
* for ``Control, DrugA, DrugB`` the rows are, in this order,
  ``DrugA_vs_Control``, ``DrugB_vs_Control``, ``DrugB_vs_DrugA``.
"""

from __future__ import annotations

import csv
import io
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

import yaml

from rnaseq.models import DesignType, Preset, ProjectConfig
from rnaseq.project import derive_design
from rnaseq.validators import CONTRAST_HEADER, ContrastDefinition, MetadataSummary, ValidationReport, condition_factor

CONTRAST_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class DesignSetupError(ValueError):
    """A design proposal cannot be produced or written; the message is user-facing."""


@dataclass(frozen=True)
class DesignProposal:
    factor: str
    levels: tuple[str, ...]
    counts: dict[str, int]
    contrasts: tuple[ContrastDefinition, ...]
    formula: str
    standard_formula: bool
    has_batch: bool
    swapped: bool = False


def discover_levels(metadata: MetadataSummary, factor: str) -> tuple[str, ...]:
    """Return ``factor`` levels in order of first appearance in the metadata rows."""

    return tuple(dict.fromkeys(record[factor] for record in metadata.rows))


def generate_contrasts(factor: str, levels: tuple[str, ...]) -> tuple[ContrastDefinition, ...]:
    """Return every pairwise contrast, later level over earlier level.

    Ordered by denominator (the earlier level) and then numerator, giving
    ``k * (k - 1) / 2`` rows for ``k`` levels.
    """

    contrasts: list[ContrastDefinition] = []
    for earlier_index, denominator in enumerate(levels[:-1]):
        for numerator in levels[earlier_index + 1:]:
            contrasts.append(ContrastDefinition(f"{numerator}_vs_{denominator}", factor, numerator, denominator))
    unsafe = sorted({level for level in levels if not CONTRAST_ID_PATTERN.fullmatch(level)})
    if unsafe:
        raise DesignSetupError(
            "Condition level(s) cannot be used in a contrast ID (use letters, digits, '.', '_' or '-', "
            "starting with a letter or digit): " + ", ".join(repr(level) for level in unsafe)
            + ". Rename the level(s) in metadata.csv, or write contrasts.csv by hand."
        )
    ids = [contrast.contrast_id for contrast in contrasts]
    if len(set(ids)) != len(ids):
        raise DesignSetupError("Condition levels produce colliding contrast IDs: " + ", ".join(sorted({i for i in ids if ids.count(i) > 1})))
    return tuple(contrasts)


def build_proposal(
    report: ValidationReport, *, condition_order: tuple[str, ...] | None = None, swap: bool = False
) -> DesignProposal:
    """Derive the proposal for a validated project; raises :class:`DesignSetupError` if none applies."""

    config, metadata = report.config, report.metadata
    if config is None or metadata is None:
        raise DesignSetupError("The project configuration or metadata could not be read.")
    if config.project.preset is Preset.QC:
        raise DesignSetupError("This project's scope is QC-only; it has no statistical design to configure.")
    factor = condition_factor(report.formula_variables, report.design_variable_types)
    if factor is None or factor not in metadata.columns:
        raise DesignSetupError("The design formula has no categorical condition factor present in metadata.csv.")
    observed = discover_levels(metadata, factor)
    if len(observed) < 2:
        raise DesignSetupError(f"Metadata column {factor!r} has fewer than two levels; there is nothing to compare.")
    levels = observed
    frozen = config.design.condition_order
    requested = condition_order if condition_order is not None else (tuple(frozen) if frozen else None)
    if requested is not None:
        if set(requested) != set(observed) or len(requested) != len(observed):
            raise DesignSetupError(
                f"The requested condition order ({', '.join(requested)}) must list exactly the levels in "
                f"metadata ({', '.join(observed)})."
            )
        levels = tuple(requested)
    if swap and len(levels) != 2:
        raise DesignSetupError("Only a two-level design can be swapped; use --condition-order for more levels.")
    counts = {level: sum(1 for record in metadata.rows if record[factor] == level) for level in levels}
    variables = config.design.variables
    has_batch = "batch" in report.formula_variables
    standard = (
        variables is not None
        and config.design.formula == derive_design(config.design.type, config.design.pair_id, has_batch)[0]
        and all(kind.value == "categorical" for kind in variables.values())
    )
    contrasts = generate_contrasts(factor, levels)
    if swap:
        # Swapping reverses only the contrast direction.  `levels` stays the
        # condition order (first appearance, or explicit) that is frozen as
        # design.condition_order; contrasts.csv is the authority on direction.
        contrasts = tuple(
            ContrastDefinition(f"{c.denominator}_vs_{c.numerator}", c.factor, c.denominator, c.numerator) for c in contrasts
        )
    return DesignProposal(factor, levels, counts, contrasts, config.design.formula, standard, has_batch, swap)


def contrasts_csv_text(contrasts: tuple[ContrastDefinition, ...]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(CONTRAST_HEADER)
    for contrast in contrasts:
        writer.writerow([contrast.contrast_id, contrast.factor, contrast.numerator, contrast.denominator])
    return buffer.getvalue()


def _indent_block(text: str, spaces: int = 2) -> str:
    return "".join(" " * spaces + line + "\n" for line in text.splitlines())


def _set_top_level_block(text: str, key: str, value: object, *, before: str | None = None) -> str:
    """Replace (or insert) one top-level YAML mapping entry, leaving other text intact."""

    block = yaml.safe_dump({key: value}, sort_keys=False, allow_unicode=True, default_flow_style=False)
    lines = text.splitlines(keepends=True)
    if lines and not lines[-1].endswith("\n"):
        lines[-1] += "\n"
    start = next((i for i, line in enumerate(lines) if re.match(rf"{re.escape(key)}:(\s|$)", line)), None)
    if start is None:
        anchor = next((i for i, line in enumerate(lines) if before and re.match(rf"{re.escape(before)}:(\s|$)", line)), len(lines))
        lines[anchor:anchor] = [block]
    else:
        end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() and not lines[i][0].isspace()), len(lines))
        lines[start:end] = [block]
    return "".join(lines)


def updated_project_yaml(
    text: str, levels: tuple[str, ...], *,
    enrichment: tuple[str, ...] | None = None, annotation: dict[str, str] | None = None,
) -> str:
    """Return L2 ``project.yaml`` text, optionally with an enrichment selection and annotation.

    ``enrichment`` ``None`` leaves ``analysis`` untouched; ``annotation`` is
    written only when given (the caller never replaces an existing contract).
    """

    edited = _updated_design_yaml(text, levels)
    if annotation is not None:
        edited = _set_top_level_block(edited, "annotation", dict(annotation), before="analysis")
    if enrichment is not None:
        edited = _set_top_level_block(edited, "analysis", {"enrichment": list(enrichment)})
    return edited


def _updated_design_yaml(text: str, levels: tuple[str, ...]) -> str:
    """Return ``project.yaml`` text with ``preset: L2`` and the frozen condition order.

    Edits only those two entries so the user's other content, ordering and
    comments are preserved; if the expected layout is not found the file is
    re-serialised instead.  The result is always re-parsed and checked.
    """

    order_block = "  condition_order:\n" + _indent_block(yaml.safe_dump(list(levels), default_flow_style=False), 4)
    lines = text.splitlines(keepends=True)
    try:
        project_start = next(i for i, line in enumerate(lines) if re.fullmatch(r"project:\s*", line.rstrip("\n")))
        preset_index = next(
            i for i in range(project_start + 1, len(lines))
            if re.match(r"\s+preset:\s", lines[i]) and not re.match(r"\S", lines[i])
        )
        lines[preset_index] = re.sub(r"(preset:\s*)\S+", r"\g<1>L2", lines[preset_index], count=1)
        design_start = next(i for i, line in enumerate(lines) if re.fullmatch(r"design:\s*", line.rstrip("\n")))
        design_end = next(
            (i for i in range(design_start + 1, len(lines)) if lines[i].strip() and not lines[i][0].isspace()), len(lines)
        )
        block = lines[design_start:design_end]
        keep: list[str] = []
        skipping = False
        for line in block:
            if re.match(r"  condition_order:", line):
                skipping = True
                continue
            if skipping and (line.startswith("    ") or line.lstrip().startswith("- ")):
                continue
            skipping = False
            keep.append(line)
        while keep and not keep[-1].strip():
            keep.pop()
        if keep and not keep[-1].endswith("\n"):
            keep[-1] += "\n"
        lines[design_start:design_end] = [*keep, order_block]
        edited = "".join(lines)
        parsed = yaml.safe_load(edited)
        if parsed["project"]["preset"] == "L2" and parsed["design"]["condition_order"] == list(levels):
            return edited
    except (StopIteration, KeyError, TypeError, yaml.YAMLError):
        pass
    parsed = yaml.safe_load(text)
    parsed["project"]["preset"] = "L2"
    parsed["design"]["condition_order"] = list(levels)
    return yaml.safe_dump(parsed, sort_keys=False, allow_unicode=True)


def _atomic_write(path: Path, content: str) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def write_design(
    project_dir: Path, config_path: Path, contrasts_path: Path, proposal: DesignProposal, *,
    enrichment: tuple[str, ...] | None = None, annotation: dict[str, str] | None = None,
) -> tuple[str, bytes]:
    """Write contrasts and the frozen design; return the prior ``(project.yaml, contrasts.csv)`` for rollback."""

    previous_yaml = config_path.read_text(encoding="utf-8")
    previous_contrasts = contrasts_path.read_bytes()
    new_yaml = updated_project_yaml(previous_yaml, proposal.levels, enrichment=enrichment, annotation=annotation)
    # Validate the model before touching anything on disk.
    ProjectConfig.model_validate(yaml.safe_load(new_yaml))
    _atomic_write(contrasts_path, contrasts_csv_text(proposal.contrasts))
    try:
        _atomic_write(config_path, new_yaml)
    except BaseException:
        contrasts_path.write_bytes(previous_contrasts)
        raise
    return previous_yaml, previous_contrasts


def rollback_design(config_path: Path, contrasts_path: Path, previous: tuple[str, bytes]) -> None:
    _atomic_write(config_path, previous[0])
    contrasts_path.write_bytes(previous[1])


def design_is_paired(config: ProjectConfig) -> bool:
    return config.design.type is DesignType.PAIRED_TWO_GROUP
