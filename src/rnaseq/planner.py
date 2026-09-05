"""Deterministic Milestone 1 planning-artifact generation."""

from __future__ import annotations

import csv
import hashlib
import os
import tempfile
from io import StringIO
from pathlib import Path

import yaml

from rnaseq.models import FastqPreprocessing, InputType, PIPELINE_VERSION, Preset
from rnaseq.hisat2_featurecounts import HISAT2_VERSION, SAMTOOLS_VERSION, SUBREAD_VERSION
from rnaseq.validators import ValidationReport

PLANNED_STATUS = "PLANNED — awaiting immutable case/run execution"
DISABLED_STATUS = "DISABLED — not selected in the current project configuration"
L1_MODULES = (
    ("count_import", "Count import"),
    ("normalization", "Normalization"),
    ("vst", "Variance-stabilizing transformation (VST)"),
    ("library_size_assessment", "Library-size assessment"),
    ("pca", "Principal component analysis (PCA)"),
    ("sample_correlation", "Sample correlation"),
    ("expression_qc", "Expression-level QC"),
)
L2_ADDITIONAL_MODULES = (
    ("differential_expression", "DESeq2 differential expression"),
    ("deg_tables", "Differentially expressed gene tables"),
    ("volcano", "Volcano plot"),
    ("deg_heatmap", "DEG heatmap"),
    ("functional_enrichment", "Functional enrichment"),
    ("standardized_report", "Standardized report"),
)


def planned_modules(preset: Preset):
    if preset is Preset.QC:
        return ()
    return L1_MODULES if preset is Preset.L1 else L1_MODULES + L2_ADDITIONAL_MODULES


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_valid_report(report: ValidationReport) -> None:
    if not report.is_valid:
        raise ValueError("Cannot generate a plan for an invalid project.")
    if not report.loaded or not report.config:
        raise ValueError("Validated project is missing required planning data.")
    if report.config.project.preset is not Preset.QC and not (report.metadata and report.contrasts):
        raise ValueError("Validated analytical project is missing metadata or contrasts.")
    if report.config.input.type is InputType.RAW_COUNTS and report.counts is None:
        raise ValueError("Validated raw-count project is missing count data.")
    if report.config.input.type is InputType.FASTQ and report.fastq is None:
        raise ValueError("Validated FASTQ project is missing FASTQ data.")


def _relative_to_project(report: ValidationReport, path: Path) -> str:
    return path.resolve().relative_to(report.project_dir.resolve()).as_posix()


def render_samplesheet(report: ValidationReport) -> str:
    """Render the four-column nf-core/rnaseq 3.26.0 samplesheet."""

    _require_valid_report(report)
    if report.fastq is None or report.config is None:
        raise ValueError("Samplesheets are only generated for FASTQ projects.")
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
    for record in report.fastq.records:
        writer.writerow(
            [
                record.sample_id,
                _relative_to_project(report, record.fastq_1),
                _relative_to_project(report, record.fastq_2) if record.fastq_2 else "",
                report.config.upstream.strandedness,
            ]
        )
    return output.getvalue()


def _contrast_lines(report: ValidationReport) -> list[str]:
    assert report.contrasts is not None
    lines: list[str] = []
    for contrast in report.contrasts.contrasts:
        lines.extend(
            [
                f"### {contrast.contrast_id}", "", f"- Factor: `{contrast.factor}`",
                f"- Comparison: {contrast.numerator} vs {contrast.denominator}",
                f"- log2FC > 0: higher in {contrast.numerator}",
                f"- log2FC < 0: higher in {contrast.denominator}", "",
            ]
        )
    return lines


def render_analysis_plan(report: ValidationReport) -> str:
    _require_valid_report(report)
    config = report.config
    assert config is not None
    lines = [
        "# Bulk RNA-seq Analysis Plan", "", "## Project", "",
        f"- Project ID: `{config.project.id}`",
        f"- Project schema version: `{config.schema_version}`",
        f"- Pipeline software version: `{PIPELINE_VERSION}`",
        f"- Preset: `{config.project.preset.value}`", "", "## Input", "",
        f"- Input type: `{config.input.type.value}`",
    ]
    if report.counts is not None:
        lines.extend([f"- Count matrix: `{config.input.path}`", f"- Genes: {report.counts.gene_count}", f"- Biological samples: {len(report.counts.sample_ids)}", f"- All-zero genes: {report.counts.all_zero_genes}"])
    else:
        assert report.fastq is not None
        files_per_record = 2 if report.fastq.layout.value == "paired_end" else 1
        preprocessing_arguments = ["--skip_trimming"] if config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
        lines.extend([f"- FASTQ directory: `{config.input.path}`", f"- FASTQ files: {len(report.fastq.records) * files_per_record}", f"- Biological samples: {len(report.fastq.sample_ids)}", f"- Sequencing layout: `{report.fastq.layout.value}`", f"- FASTQ preprocessing: `{config.input.preprocessing.value}`", f"- nf-core skip_trimming: `{str(config.input.preprocessing is FastqPreprocessing.PRETRIMMED).lower()}`", f"- nf-core preprocessing arguments: `{', '.join(preprocessing_arguments) if preprocessing_arguments else 'none'}`", "- nf-core samplesheet: `planning/samplesheet.csv`"])
    if config.project.preset is Preset.QC:
        lines.extend([
            "", "## Analysis scope", "", "- Quantification and technical FASTQ QC only.",
            "- No downstream normalization, differential expression, metadata design, or contrasts are run.", "",
        ])
    else:
        lines.extend([
        "", "## Species", "", config.organism.species.value, "", "## Experimental design", "",
        f"- Design type: `{config.design.type.value}`", f"- Formula: `{config.design.formula}`",
        f"- Design variables: {', '.join(report.formula_variables)}", "",
        "Only variables explicitly present in the formula are part of the planned design. Additional metadata columns are preserved but are not added automatically.",
        "", "## Contrasts", "",
        ])
        lines.extend(_contrast_lines(report))
    if config.input.type is InputType.FASTQ:
        reference = config.reference
        assert reference is not None
        method = config.upstream.quantification.method if config.upstream.quantification else "salmon"
        lines.extend([
            "## Upstream processing", "", f"- Backend: `{method}`",
            f"- Strandedness: `{config.upstream.strandedness}`", f"- Reference source: `{reference.source}`",
        ])
        if method == "salmon":
            lines.append(f"- Pipeline: `nf-core/rnaseq {config.upstream.pipeline_version}`")
        else:
            lines.extend([
                f"- Pinned tools: HISAT2 {HISAT2_VERSION}; SAMtools {SAMTOOLS_VERSION}; Subread/featureCounts {SUBREAD_VERSION}.",
                "- Counting: exon/gene_id; multimappers, multi-gene overlaps, fractional counts and duplicate removal are disabled; MAPQ is 0.",
                "- Paired-end uses fragment counting with both mates mapped and chimeras excluded; single-end uses read counting.",
            ])
        if report.local_reference is not None:
            local = report.local_reference
            lines.extend([
                f"- Local reference identity: `{local.species}` / `{local.provider}` release `{local.release}` / `{local.assembly}` patch `{local.assembly_patch}`",
                f"- Manifest: `{local.manifest_path}` (SHA256 `{local.manifest_sha256}`)",
                f"- Genome FASTA: `{local.genome_fasta.path}` (SHA256 `{local.genome_fasta.sha256}`)",
                f"- GTF: `{local.annotation_gtf.path}` (SHA256 `{local.annotation_gtf.sha256}`)",
                f"- Transcript FASTA asset: `{local.transcript_fasta.path}` (SHA256 `{local.transcript_fasta.sha256}`; runtime used: `{str(local.external_transcript_fasta_used).lower()}`)",
                f"- Transcriptome strategy: `{local.transcriptome_strategy}`",
                f"- {'Salmon' if method == 'salmon' else 'HISAT2'} index strategy: {local.salmon_strategy if method == 'salmon' else local.hisat2_status}",
            ])
            if local.salmon_transcriptome is not None:
                lines.append(
                    f"- Adopted GTF-derived transcriptome: `{local.salmon_transcriptome.path}` "
                    f"(SHA256 `{local.salmon_transcriptome.sha256}`; provenance only, not an nf-core `--transcript_fasta` argument)"
                )
            if method == "salmon" and local.salmon_index is not None:
                lines.append("- nf-core reference arguments: " + " ".join(
                    f"`{option} {path}`" for option, path in local.nfcore_arguments()
                ))
            elif method == "salmon":
                lines.append(f"- Reference preparation required: `rnaseq reference prepare {local.root}`")
            elif local.hisat2_index is None:
                lines.append(f"- Reference preparation required: `rnaseq reference prepare-hisat2 {local.root}`")
        else:
            lines.append(f"- Reference genome: `{reference.genome}`")
        lines.extend([
            f"- FASTQ QC, alignment / quantification, and MultiQC: **PLANNED — executed by the immutable case/run {'nf-core/rnaseq' if method == 'salmon' else 'first-party HISAT2 + featureCounts'} route, not by `rnaseq plan`**", "",
            "## Execution readiness", "", f"- Status: `{'READY' if report.execution_ready else 'NOT READY'}`",
        ])
        lines.extend(f"- Blocking requirement: {item}" for item in report.execution_blockers)
        lines.append("")
    lines.extend(["## Downstream module status", ""])
    selected_enrichment = set(config.analysis.enrichment) if config.analysis is not None else set()
    for module, label in planned_modules(config.project.preset):
        status = DISABLED_STATUS if module == "functional_enrichment" and not selected_enrichment else PLANNED_STATUS
        lines.append(f"- {label}: **{status}**")
    if selected_enrichment:
        lines.extend([
            "", "## Functional enrichment", "", "- Enrichment method: `GSEA`",
            "- Gene-set resources: `GO BP`, `GO MF`, `GO CC`, `KEGG`",
            "- Internal backends: `gsea-go`, `gsea-kegg`", "",
        ])
    lines.extend(["", "## Warnings", ""])
    if report.warnings:
        lines.extend(f"- {issue.message}" for issue in report.warnings)
    else:
        lines.append("None")
    lines.extend([
        "", "## Limitations", "",
        "- This planning artifact is non-executing; execution is initiated explicitly through the immutable case/run service.",
        "- Operational Nextflow cache and work data are local execution scratch, while frozen inputs, outputs, logs and provenance persist in the case/run directory.",
        "- A module marked PLANNED has not been run for this project; DISABLED means it is not selected by the current configuration.", "",
    ])
    return "\n".join(lines)


def _fastq_manifest_files(report: ValidationReport) -> list[dict[str, str]]:
    assert report.fastq is not None
    files = [record.fastq_1 for record in report.fastq.records]
    files.extend(record.fastq_2 for record in report.fastq.records if record.fastq_2)
    return [{"relative_path": _relative_to_project(report, path), "sha256": sha256_file(path)} for path in sorted(files)]


def render_manifest(report: ValidationReport) -> str:
    _require_valid_report(report)
    loaded, config, metadata, contrasts = report.loaded, report.config, report.metadata, report.contrasts
    assert loaded and config
    manifest: dict[str, object] = {
        "schema": {"project_schema_version": config.schema_version},
        "pipeline": {"name": config.project.pipeline, "version": PIPELINE_VERSION, "preset": config.project.preset.value},
        "project": {"id": config.project.id},
        "configuration": {"file": "project.yaml", "sha256": sha256_file(loaded.config_path)},
    }
    if report.counts is not None:
        manifest["input"] = {"type": config.input.type.value, "file": config.input.path, "sha256": sha256_file(report.counts.path), "genes": report.counts.gene_count, "samples": len(report.counts.sample_ids), "all_zero_genes": report.counts.all_zero_genes}
    else:
        assert report.fastq is not None
        manifest["input"] = {"type": "fastq", "files": _fastq_manifest_files(report)}
        manifest["upstream"] = {
            "engine": config.upstream.engine,
            "pipeline_version": config.upstream.pipeline_version,
            "sequencing_layout": report.fastq.layout.value,
            "strandedness": config.upstream.strandedness,
            "preprocessing": config.input.preprocessing.value,
            "skip_trimming": config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
            "nfcore_preprocessing_arguments": (
                ["--skip_trimming"]
                if config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
            ),
            "nfcore_runtime_params": (
                {"skip_alignment": True, "skip_trimming": True}
                if config.input.preprocessing is FastqPreprocessing.PRETRIMMED
                else {"skip_alignment": True}
            ),
            "quantification": (
                config.upstream.quantification.model_dump()
                if config.upstream.quantification is not None
                else None
            ),
        }
        assert config.reference is not None
        manifest["reference"] = (
            report.local_reference.provenance()
            if report.local_reference is not None
            else config.reference.model_dump()
        )
    manifest.update({
        "metadata": ({"file": config.metadata_file, "sha256": sha256_file(metadata.path), "design_variables": list(report.formula_variables)} if metadata is not None else None),
        "contrasts": ({"file": config.contrasts_file, "sha256": sha256_file(contrasts.path), "definitions": [item.__dict__ for item in contrasts.contrasts]} if contrasts is not None else None),
        "organism": {"species": config.organism.species.value},
        "design": {"type": config.design.type.value, "formula": config.design.formula},
        "planned_downstream": {
            "preset": config.project.preset.value,
            "enrichment": list(config.analysis.enrichment) if config.analysis is not None else [],
        },
        "thresholds": {"padj": config.thresholds.padj, "abs_log2fc": config.thresholds.abs_log2fc},
    })
    return yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True)


def render_upstream_preview(report: ValidationReport) -> str:
    _require_valid_report(report)
    if report.config is None or report.config.input.type is not InputType.FASTQ or report.fastq is None:
        raise ValueError("Upstream previews are only generated for FASTQ projects.")
    config = report.config
    assert config.reference is not None
    method = config.upstream.quantification.method if config.upstream.quantification else "salmon"
    preview = {
        "engine": {"name": "nextflow"},
        "pipeline": ({"name": "nf-core/rnaseq", "version": config.upstream.pipeline_version} if method == "salmon" else {"name": "nf-rna/hisat2_featurecounts", "versions": {"hisat2": HISAT2_VERSION, "samtools": SAMTOOLS_VERSION, "subread": SUBREAD_VERSION}}),
        "input": {"samplesheet": "planning/samplesheet.csv", "project_working_directory": "."},
        "sequencing": {"layout": report.fastq.layout.value, "strandedness": config.upstream.strandedness},
        "preprocessing": config.input.preprocessing.value,
        "skip_trimming": config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
        "nfcore_preprocessing_arguments": (
            ["--skip_trimming"]
            if config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
        ),
        "nfcore_runtime_params": (
            {"skip_alignment": True, "skip_trimming": True}
            if config.input.preprocessing is FastqPreprocessing.PRETRIMMED
            else {"skip_alignment": True}
        ),
        "quantification": (
            config.upstream.quantification.model_dump()
            if config.upstream.quantification is not None
            else None
        ),
        "reference": (
            report.local_reference.provenance()
            if report.local_reference is not None
            else config.reference.model_dump()
        ),
        "nfcore_reference_arguments": (
            [item for option, path in report.local_reference.nfcore_arguments() for item in (option, str(path))]
            if report.local_reference is not None and report.local_reference.salmon_index is not None
            else (["--genome", config.reference.genome] if config.reference.source == "igenomes" and config.reference.genome else [])
        ),
        "salmon_index_strategy": (
            report.local_reference.salmon_strategy
            if report.local_reference is not None
            else "configured by the iGenomes/custom reference route"
        ),
        "nfcore_dynamic_salmon_index_build": (
            False
            if report.local_reference is not None and report.local_reference.salmon_index is not None
            else None
        ),
        "transcriptome_strategy": (
            report.local_reference.transcriptome_strategy
            if report.local_reference is not None
            else "configured by the iGenomes/custom reference route"
        ),
        "execution": {"profile": "local"},
        "status": {"execution_ready": report.execution_ready},
        "blocking_requirements": list(report.execution_blockers),
    }
    if method == "hisat2_featurecounts":
        preview["counting"] = {"feature_type": "exon", "grouping_attribute": "gene_id", "paired_end": "-p --countReadPairs -B -C", "single_end": "read", "strand_mapping": {"unstranded": 0, "forward": 1, "reverse": 2}}
        preview["reference_readiness"] = "HISAT2 index required; iGenomes is not an execution route for this backend"
    return yaml.safe_dump(preview, sort_keys=False, allow_unicode=True)


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def generate_plan(report: ValidationReport) -> tuple[Path, ...]:
    """Write deterministic planning files for a valid project."""

    _require_valid_report(report)
    planning_dir = report.project_dir / "planning"
    files = {planning_dir / "analysis_plan.md": render_analysis_plan(report), planning_dir / "manifest.preview.yaml": render_manifest(report)}
    if report.config is not None and report.config.input.type is InputType.FASTQ:
        files[planning_dir / "samplesheet.csv"] = render_samplesheet(report)
        files[planning_dir / "upstream_run.preview.yaml"] = render_upstream_preview(report)
    for path, content in files.items():
        _atomic_write(path, content)
    return tuple(files)
