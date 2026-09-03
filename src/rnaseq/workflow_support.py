"""Small non-statistical helpers used inside the downstream Nextflow processes."""

from __future__ import annotations

import argparse
import base64
import csv
import json
from html import escape
from pathlib import Path
from typing import Any

import yaml

from rnaseq.models import normalize_enrichment_selection


FILTER = {"rule": "remove genes with total count < 10 after all-zero removal", "minimum_total_count": 10}


REPORT_STYLES = """<style>
figure { max-width: 100%; margin: 1rem 0; }
img.report-figure {
  display: block;
  max-width: 100%;
  width: auto;
  height: auto;
  object-fit: contain;
}
</style>"""


_ENRICHMENT_REQUIRED_PATHS = {
    "gsea-go": (
        "annotation.organism", "annotation.input_id_type", "annotation.minimum_mapping_rate",
        "annotation.enrichment.gsea.minimum_ranked_genes", "annotation.enrichment.gsea.min_gs_size",
        "annotation.enrichment.gsea.max_gs_size", "annotation.enrichment.gsea.pvalue_cutoff",
        "annotation.enrichment.gsea.padj_cutoff", "annotation.enrichment.gsea.p_adjust_method",
        "annotation.enrichment.gsea.seed",
    ),
    "gsea-kegg": (
        "annotation.organism", "annotation.input_id_type", "annotation.minimum_mapping_rate",
        "annotation.enrichment.kegg.resource_provider", "annotation.enrichment.kegg.gsea.minimum_ranked_genes",
        "annotation.enrichment.kegg.gsea.pvalue_cutoff", "annotation.enrichment.kegg.gsea.padj_cutoff",
        "annotation.enrichment.kegg.gsea.p_adjust_method", "annotation.enrichment.kegg.gsea.min_gs_size",
        "annotation.enrichment.kegg.gsea.max_gs_size", "annotation.enrichment.kegg.gsea.seed",
    ),
}


def _read_contract(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("downstream contract must be a JSON object")
    return value


def _read_execution_inputs(inputs: Path) -> tuple[Path, dict[str, Any]]:
    """Read a staged-only execution manifest without resolving provenance paths."""

    manifest_path = inputs / "execution_inputs.json"
    value = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("staged execution inputs must contain a JSON object")
    return inputs, value


def _staged_file(root: Path, configured: object, label: str) -> Path:
    if not isinstance(configured, str) or not configured or Path(configured).is_absolute():
        raise ValueError(f"staged execution inputs has invalid {label}")
    root = root.resolve()
    candidate = root / configured
    try:
        candidate = candidate.resolve(strict=True)
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"staged execution inputs {label} escapes its task directory") from exc
    except OSError as exc:
        raise ValueError(f"staged execution inputs cannot resolve {label}: {exc}") from exc
    if candidate.name.startswith("._"):
        raise ValueError(f"staged execution inputs {label} selected an AppleDouble artifact")
    if not candidate.is_file():
        raise ValueError(f"staged execution inputs is missing {label}: {candidate}")
    try:
        with candidate.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise ValueError(f"staged execution inputs cannot read {label}: {candidate}: {exc}") from exc
    return candidate


def _project(inputs: Path) -> dict[str, Any]:
    inputs_root, manifest = _read_execution_inputs(inputs)
    value = yaml.safe_load(_staged_file(inputs_root, manifest.get("project_config"), "project_config").read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("staged project configuration must be a YAML mapping")
    return value


def _samples(inputs: Path) -> list[str]:
    inputs_root, manifest = _read_execution_inputs(inputs)
    samples = manifest.get("samples")
    if not isinstance(samples, list) or not samples or not all(isinstance(item, str) and item for item in samples):
        raise ValueError("staged execution inputs has no valid sample list")
    if samples != sorted(samples) or len(set(samples)) != len(samples):
        raise ValueError("staged execution inputs sample list must be sorted and unique")
    metadata = _staged_file(inputs_root, manifest.get("metadata"), "metadata")
    with metadata.open(encoding="utf-8", newline="") as handle:
        observed = [row.get("sample_id") for row in csv.DictReader(handle)]
    if sorted(observed) != samples or len(set(observed)) != len(observed):
        raise ValueError("staged metadata sample IDs disagree with staged execution inputs sample list")
    return list(samples)


def _source_config(contract: dict[str, Any], inputs: Path) -> dict[str, Any]:
    root, manifest = _read_execution_inputs(inputs)
    source = manifest.get("source")
    if not isinstance(source, dict) or source.get("type") != contract.get("source", {}).get("type"):
        raise ValueError("staged execution source disagrees with frozen provenance contract")
    if source["type"] == "raw_counts":
        return {"source_type": "raw_counts", "counts": str(_staged_file(root, source.get("counts"), "raw-count matrix"))}
    samples = _samples(inputs)
    quant = source.get("quant_sf")
    if not isinstance(quant, dict) or set(quant) != set(samples):
        raise ValueError("staged Salmon quant.sf sample mapping disagrees with metadata samples")
    return {
        "source_type": "salmon_tximport",
        "quant_sf": {sample: str(_staged_file(root, quant[sample], f"salmon.quant_sf.{sample}")) for sample in samples},
        "tx2gene": str(_staged_file(root, source.get("tx2gene"), "salmon.tx2gene")),
    }


def _contrasts(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _required_frozen_annotation(contract: dict[str, Any], kind: str) -> dict[str, Any]:
    """Read only the immutable annotation snapshot required by one backend."""

    if kind not in _ENRICHMENT_REQUIRED_PATHS:
        raise ValueError(f"unsupported enrichment module: {kind}")
    annotation = contract.get("annotation")
    if not isinstance(annotation, dict):
        raise ValueError("frozen enrichment configuration is missing annotation")
    for dotted_path in _ENRICHMENT_REQUIRED_PATHS[kind]:
        current: Any = contract
        for part in dotted_path.split("."):
            if not isinstance(current, dict) or part not in current or current[part] is None:
                raise ValueError(f"frozen enrichment configuration is missing {dotted_path}")
            current = current[part]
    return annotation


def l1_config(contract_path: Path, inputs: Path, output: Path) -> dict[str, Any]:
    contract = _read_contract(contract_path)
    root, manifest = _read_execution_inputs(inputs)
    project = _project(inputs)
    return {
        **_source_config(contract, inputs),
        "metadata": str(_staged_file(root, manifest.get("metadata"), "metadata")),
        "formula": project["design"]["formula"],
        "samples": _samples(inputs),
        "output_dir": str(output),
        "filter": FILTER,
    }


def l2_config(contract_path: Path, inputs: Path, l1: Path, output: Path) -> dict[str, Any]:
    contract = _read_contract(contract_path)
    root, manifest = _read_execution_inputs(inputs)
    project = _project(inputs)
    return {
        **_source_config(contract, inputs),
        "metadata": str(_staged_file(root, manifest.get("metadata"), "metadata")),
        "formula": project["design"]["formula"],
        "samples": _samples(inputs),
        "output_dir": str(output),
        "filter": FILTER,
        "thresholds": project["thresholds"],
        "contrasts": _contrasts(_staged_file(root, manifest.get("contrasts"), "contrasts")),
        "l1_vst": str(l1 / "vst.tsv"),
        "heatmap_top_n": 50,
    }


def _read_json_artifact(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"final report requires {label}: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"final report requires {label} to be a JSON object: {path}")
    return value


def _image_html(path: Path, caption: str) -> str:
    if not path.is_file():
        return f"<p>{escape(caption)}: not produced.</p>"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        f"<figure><figcaption>{escape(caption)}</figcaption>"
        f"<img class='report-figure' alt='{escape(path.name)}' "
        f"src='data:image/png;base64,{encoded}'></figure>"
    )


def _top_table_html(path: Path, columns: tuple[str, ...], *, top_n: int = 10) -> str:
    if not path.is_file():
        return "<p>Complete TSV was not produced.</p>"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in columns if column not in fieldnames]
        if missing_columns:
            raise ValueError(
                f"final report TSV {path} is missing required report columns: "
                + ", ".join(missing_columns)
            )
        visible: list[dict[str, str]] = []
        row_count = 0
        for row in reader:
            row_count += 1
            if not any((row.get(column) or "").strip() for column in columns):
                raise ValueError(
                    f"final report TSV {path} contains a row with no values for requested report columns: "
                    + ", ".join(columns)
                )
            if len(visible) < top_n:
                visible.append(row)
    headers = "".join(f"<th>{escape(column)}</th>" for column in columns)
    body = "".join(
        "<tr>" + "".join(f'<td>{escape(row.get(column, "") or "")}</td>' for column in columns) + "</tr>"
        for row in visible
    )
    table = f"<table><thead><tr>{headers}</tr></thead><tbody>{body}</tbody></table>" if visible else "<p>No returned rows.</p>"
    return f"<p>Showing {len(visible)} of {row_count} rows. Complete TSV: <code>{escape(str(path))}</code>.</p>{table}"


def _tsv_row_count(path: Path) -> int:
    """Count data rows in an analysis TSV using its explicit tab contract."""

    with path.open(encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def _pipeline_provenance_line(execution: object, project_pipeline: object) -> str:
    """Render one clear pipeline provenance value without duplicate fallbacks."""

    pipeline = execution.get("pipeline", {}) if isinstance(execution, dict) else {}
    if not isinstance(pipeline, dict):
        pipeline = {}
    name = pipeline.get("name")
    version = pipeline.get("version")
    if name and version:
        return f"Pipeline: {escape(str(name))} {escape(str(version))}."
    if name:
        return f"Pipeline: {escape(str(name))}; version: not available."
    if version:
        return f"Pipeline version: {escape(str(version))}."
    if project_pipeline:
        return f"Pipeline: {escape(str(project_pipeline))}; version: not available."
    return "Pipeline version: not available."


def _value(mapping: dict[str, Any], name: str) -> str:
    value = mapping.get(name)
    return "not available" if value is None else escape(str(value))


def _require_usable_gsea_summary(summary: dict[str, Any], label: str) -> None:
    status = summary.get("status")
    if status not in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}:
        raise ValueError(f"final report cannot claim enabled {label} completed: backend status is {status!r}")


def _report_l1_only(
    contract: dict[str, Any], inputs_root: Path, manifest: dict[str, Any], project: dict[str, Any],
    l1: Path, output: Path,
) -> None:
    """Render an L1-only technical report without creating an L2 dependency."""

    if normalize_enrichment_selection(contract.get("analysis", {}).get("enrichment", [])):
        raise ValueError("L1 final report cannot include enrichment.")
    case = contract["case"]
    samples = _samples(inputs_root)
    l1_summary = _read_json_artifact(l1 / "backend_summary.json", "L1 backend summary")
    contrasts = _contrasts(_staged_file(inputs_root, manifest.get("contrasts"), "contrasts"))
    execution_path = inputs_root / "execution_manifest.yaml"
    execution = yaml.safe_load(execution_path.read_text(encoding="utf-8")) if execution_path.is_file() else {}
    package_versions = l1_summary.get("package_versions", {})
    sections = [
        "<!doctype html><html><head><meta charset='utf-8'><title>RNA-seq final analysis report</title>",
        REPORT_STYLES,
        "</head><body>",
        "<h1>RNA-seq final analysis report</h1>",
        "<h2>Run and analysis overview</h2>",
        "<ul>"
        f"<li>Case ID: <code>{escape(str(case['id']))}</code></li>"
        f"<li>Run ID: <code>{escape(str(case['run_id']))}</code></li>"
        f"<li>Input type: {escape(str(contract['source']['type']))}</li>"
        f"<li>Organism: {escape(str(project['organism']['species']))}</li>"
        f"<li>Samples: {len(samples)}</li>"
        f"<li>Design: <code>{escape(str(project['design']['formula']))}</code></li>"
        f"<li>Configured contrasts: {len(contrasts)}</li>"
        "<li>Requested analysis level: L1</li>"
        "</ul>",
        "<h2>L1 — quality control</h2>",
        "<ul>"
        f"<li>Genes input: {_value(l1_summary, 'genes_input')}</li>"
        f"<li>All-zero genes removed: {_value(l1_summary, 'genes_removed_all_zero')}</li>"
        f"<li>Low-total genes removed: {_value(l1_summary, 'genes_removed_low_total')}</li>"
        f"<li>Genes retained: {_value(l1_summary, 'genes_retained')}</li>"
        f"<li>Filtering: {escape(str(l1_summary.get('filter', {}).get('rule', 'not available')))}</li>"
        f"<li>Normalization: {escape(str(l1_summary.get('normalization', {}).get('method', 'not available')))}</li>"
        f"<li>VST: {escape(str(l1_summary.get('vst', {}).get('method', 'not available')))}</li>"
        "<li>Samples are never automatically excluded by this pipeline.</li>"
        "</ul>",
        _image_html(l1 / "pca.png", "PCA"),
        _top_table_html(l1 / "pca_variance.tsv", ("component", "proportion_variance"), top_n=2),
        _top_table_html(l1 / "library_size_qc.tsv", ("sample_id", "input_total", "retained_total", "size_factor", "normalized_total"), top_n=20),
        _image_html(l1 / "sample_correlation.png", "Sample correlation (blind VST)"),
        "<h2>L2 — differential expression</h2><p>Not requested for this L1 project.</p>",
        "<h2>L2 — GSEA</h2><p>Not requested for this L1 project.</p>",
        "<h2>Methods and reproducibility</h2>",
        "<ul>"
        f"<li>Import: {escape('DESeqDataSetFromMatrix raw-count import' if contract['source']['type'] == 'raw_counts' else 'Salmon/tximport import')}</li>"
        f"<li>Filtering: {escape(str(l1_summary.get('filter', {}).get('rule', 'not available')))}</li>"
        f"<li>Normalization: {escape(str(l1_summary.get('normalization', {}).get('method', 'not available')))}</li>"
        f"<li>VST: {escape(str(l1_summary.get('vst', {}).get('method', 'not available')))}</li>"
        "<li>Differential expression and GSEA were not requested for this L1 project.</li>"
        f"<li>R: {escape(str(package_versions.get('R', 'not available')))}; DESeq2: {escape(str(package_versions.get('DESeq2', 'not available')))}.</li>"
        f"<li>{_pipeline_provenance_line(execution, project.get('project', {}).get('pipeline') if isinstance(project.get('project'), dict) else None)}</li>"
        "</ul>",
        "<p>This automated report contains technical/statistical results only and no biological interpretation.</p>",
        "</body></html>",
    ]
    (output / "report.html").write_text("\n".join(sections), encoding="utf-8")


def report(contract_path: Path, inputs: Path, l1: Path, l2: Path | None, output: Path, enrichment_dirs: list[Path] | tuple[Path, ...] = ()) -> None:
    """Render the final client report from completed backend artifacts only."""

    contract = _read_contract(contract_path)
    inputs_root, manifest = _read_execution_inputs(inputs)
    project = _project(inputs)
    output.mkdir(parents=True, exist_ok=True)
    analysis_level = contract.get("analysis_level")
    if analysis_level not in {"L1", "L2"}:
        raise ValueError("final report requires frozen analysis_level L1 or L2")
    if analysis_level == "L1":
        if l2 is not None:
            raise ValueError("L1 final report must not receive L2 artifacts")
        _report_l1_only(contract, inputs_root, manifest, project, l1, output)
        return
    if l2 is None:
        raise ValueError("L2 final report requires L2 artifacts")
    case = contract["case"]
    selected = normalize_enrichment_selection(contract.get("analysis", {}).get("enrichment", []))
    samples = _samples(inputs)
    l1_summary = _read_json_artifact(l1 / "backend_summary.json", "L1 backend summary")
    contrasts = _contrasts(_staged_file(inputs_root, manifest.get("contrasts"), "contrasts"))
    thresholds = project["thresholds"]

    sections = [
        "<!doctype html><html><head><meta charset='utf-8'><title>RNA-seq final analysis report</title>",
        REPORT_STYLES,
        "</head><body>",
        "<h1>RNA-seq final analysis report</h1>",
        "<h2>Run and analysis overview</h2>",
        "<ul>"
        f"<li>Case ID: <code>{escape(str(case['id']))}</code></li>"
        f"<li>Run ID: <code>{escape(str(case['run_id']))}</code></li>"
        f"<li>Input type: {escape(str(contract['source']['type']))}</li>"
        f"<li>Organism: {escape(str(project['organism']['species']))}</li>"
        f"<li>Samples: {len(samples)}</li>"
        f"<li>Design: <code>{escape(str(project['design']['formula']))}</code></li>"
        f"<li>Configured contrasts: {len(contrasts)}</li>"
        "</ul>",
        "<h2>L1 — quality control</h2>",
        "<ul>"
        f"<li>Genes input: {_value(l1_summary, 'genes_input')}</li>"
        f"<li>All-zero genes removed: {_value(l1_summary, 'genes_removed_all_zero')}</li>"
        f"<li>Low-total genes removed: {_value(l1_summary, 'genes_removed_low_total')}</li>"
        f"<li>Genes retained: {_value(l1_summary, 'genes_retained')}</li>"
        f"<li>Filtering: {escape(str(l1_summary.get('filter', {}).get('rule', 'not available')))}</li>"
        f"<li>Normalization: {escape(str(l1_summary.get('normalization', {}).get('method', 'not available')))}</li>"
        f"<li>VST: {escape(str(l1_summary.get('vst', {}).get('method', 'not available')))}</li>"
        "<li>Samples are never automatically excluded by this pipeline.</li>"
        "</ul>",
        _image_html(l1 / "pca.png", "PCA"),
        _top_table_html(l1 / "pca_variance.tsv", ("component", "proportion_variance"), top_n=2),
        _top_table_html(l1 / "library_size_qc.tsv", ("sample_id", "input_total", "retained_total", "size_factor", "normalized_total"), top_n=20),
        _image_html(l1 / "sample_correlation.png", "Sample correlation (blind VST)"),
        "<h2>L2 — differential expression</h2>",
    ]

    for contrast in contrasts:
        contrast_id = contrast["contrast_id"]
        root = l2 / "contrasts" / contrast_id
        summary = _read_json_artifact(root / "backend_summary.json", f"L2 summary for contrast {contrast_id}")
        for required in ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv"):
            if not (root / required).is_file():
                raise ValueError(f"final report requires L2 {required} for contrast {contrast_id}")
        sections.extend([
            f"<h3>{escape(contrast_id)}</h3>",
            "<ul>"
            f"<li>Direction: {escape(contrast['numerator'])} / {escape(contrast['denominator'])}</li>"
            f"<li>Factor: {escape(contrast['factor'])}</li>"
            f"<li>Design: <code>{escape(str(project['design']['formula']))}</code></li>"
            f"<li>Tested genes: {_value(summary, 'tested_genes')}</li>"
            f"<li>Significant genes: {_tsv_row_count(root / 'significant.tsv')}</li>"
            f"<li>Upregulated genes: {_tsv_row_count(root / 'upregulated.tsv')}</li>"
            f"<li>Downregulated genes: {_tsv_row_count(root / 'downregulated.tsv')}</li>"
            f"<li>Thresholds: padj &lt; {escape(str(thresholds['padj']))}; |log2FC| ≥ {escape(str(thresholds['abs_log2fc']))}</li>"
            f"<li>Client tables: <code>tables/l2/contrasts/{escape(contrast_id)}/</code></li>"
            "</ul>",
            _image_html(root / "volcano.png", f"Volcano plot: {contrast_id}"),
            _image_html(root / "heatmap.png", f"DEG heatmap: {contrast_id}"),
        ])

    if selected == ("gsea",):
        indexed = {path.name: path for path in enrichment_dirs}
        missing = sorted({"gsea_go", "gsea_kegg"} - set(indexed))
        if missing:
            raise ValueError("final report requires enabled GSEA backend artifacts: " + ", ".join(missing))
        go_root, kegg_root = indexed["gsea_go"], indexed["gsea_kegg"]
        go_summary = _read_json_artifact(go_root / "gsea_backend_summary.json", "GO GSEA backend summary")
        kegg_summary = _read_json_artifact(kegg_root / "gsea_kegg_backend_summary.json", "KEGG GSEA backend summary")
        _require_usable_gsea_summary(go_summary, "GO GSEA")
        _require_usable_gsea_summary(kegg_summary, "KEGG GSEA")
        go_by_contrast = {str(item.get("contrast_id")): item for item in go_summary.get("contrasts", []) if isinstance(item, dict)}
        kegg_by_contrast = {str(item.get("contrast_id")): item for item in kegg_summary.get("contrasts", []) if isinstance(item, dict)}
        sections.extend([
            "<h2>L2 — preranked GSEA</h2>",
            "<p>Enrichment method: GSEA. Gene-set resources: GO BP, GO MF, GO CC, KEGG.</p>",
            "<p>GSEA uses the unmodified DESeq2 statistic-ranked gene list; no DEG, adjusted-p-value, p-value, or fold-change prefilter is applied.</p>",
        ])
        for contrast in contrasts:
            contrast_id = contrast["contrast_id"]
            go_item, kegg_item = go_by_contrast.get(contrast_id), kegg_by_contrast.get(contrast_id)
            if not isinstance(go_item, dict) or not isinstance(kegg_item, dict):
                raise ValueError(f"final report requires GO and KEGG GSEA summaries for contrast {contrast_id}")
            _require_usable_gsea_summary(go_item, f"GO GSEA contrast {contrast_id}")
            _require_usable_gsea_summary(kegg_item, f"KEGG GSEA contrast {contrast_id}")
            ranking = go_item.get("ranking", {})
            if not isinstance(ranking, dict):
                raise ValueError(f"final report requires GO GSEA ranking summary for contrast {contrast_id}")
            sections.extend([
                f"<h3>GSEA ranking: {escape(contrast_id)}</h3>",
                "<ul>"
                f"<li>Source genes: {_value(ranking, 'finite_stat_source_genes')}</li>"
                f"<li>Mapped genes: {_value(ranking, 'mapped_source_genes')}</li>"
                f"<li>Unmapped genes: {_value(ranking, 'unmapped_source_genes')}</li>"
                f"<li>Mapping rate: {_value(ranking, 'mapping_rate')}</li>"
                f"<li>Final ranked genes: {_value(ranking, 'final_ranked_genes')}</li>"
                f"<li>Positive / negative statistics: {_value(ranking, 'positive_stats')} / {_value(ranking, 'negative_stats')}</li>"
                f"<li>Tie handling: {escape(str(ranking.get('tie_handling', 'not available')))}</li>"
                "</ul>",
                f"<h3>GO GSEA: {escape(contrast_id)}</h3>",
            ])
            ontologies = go_item.get("ontologies", {})
            if not isinstance(ontologies, dict):
                raise ValueError(f"final report requires GO ontology summaries for contrast {contrast_id}")
            for ontology in ("BP", "MF", "CC"):
                outcome = ontologies.get(ontology)
                if not isinstance(outcome, dict):
                    raise ValueError(f"final report requires GO {ontology} summary for contrast {contrast_id}")
                _require_usable_gsea_summary(outcome, f"GO {ontology} GSEA contrast {contrast_id}")
                root = go_root / contrast_id / ontology
                terms = root / "all_terms.tsv"
                if not terms.is_file():
                    raise ValueError(f"final report requires GO {ontology} terms for contrast {contrast_id}")
                sections.extend([
                    f"<h4>{ontology}</h4>",
                    "<ul>"
                    f"<li>Returned terms: {_value(outcome, 'all_terms')}</li>"
                    f"<li>Significant terms: {_value(outcome, 'significant_terms')}</li>"
                    f"<li>Positive / negative terms: {_value(outcome, 'positive_terms')} / {_value(outcome, 'negative_terms')}</li>"
                    f"<li>Client tables: <code>tables/l2/enrichment/gsea_go/{escape(contrast_id)}/{ontology}/</code></li>"
                    "</ul>",
                    _image_html(root / "dotplot.png", f"GO {ontology} dotplot: {contrast_id}"),
                    _top_table_html(terms, ("ID", "Description", "NES", "pvalue", "p.adjust")),
                ])
            kegg_terms = kegg_root / contrast_id / "all_terms.tsv"
            if not kegg_terms.is_file():
                raise ValueError(f"final report requires KEGG GSEA terms for contrast {contrast_id}")
            sections.extend([
                f"<h3>KEGG GSEA: {escape(contrast_id)}</h3>",
                "<ul>"
                f"<li>Returned pathways: {_value(kegg_item, 'all_terms')}</li>"
                f"<li>Significant pathways: {_value(kegg_item, 'significant_terms')}</li>"
                f"<li>Positive / negative pathways: {_value(kegg_item, 'positive_terms')} / {_value(kegg_item, 'negative_terms')}</li>"
                f"<li>Resource provider: {escape(str(kegg_summary.get('resource', {}).get('provider', 'not available')))}</li>"
                f"<li>Client tables: <code>tables/l2/enrichment/gsea_kegg/{escape(contrast_id)}/</code></li>"
                "</ul>",
                _image_html(kegg_root / contrast_id / "dotplot.png", f"KEGG GSEA dotplot: {contrast_id}"),
                _top_table_html(kegg_terms, ("ID", "Description", "NES", "pvalue", "p.adjust")),
            ])
    else:
        sections.append("<h2>L2 — GSEA</h2><p>No enrichment module was enabled for this project.</p>")

    execution_path = inputs_root / "execution_manifest.yaml"
    execution = yaml.safe_load(execution_path.read_text(encoding="utf-8")) if execution_path.is_file() else {}
    package_versions = l1_summary.get("package_versions", {})
    sections.extend([
        "<h2>Methods and reproducibility</h2>",
        "<ul>"
        f"<li>Import: {escape('DESeqDataSetFromMatrix raw-count import' if contract['source']['type'] == 'raw_counts' else 'Salmon/tximport import')}</li>"
        f"<li>Filtering: {escape(str(l1_summary.get('filter', {}).get('rule', 'not available')))}</li>"
        f"<li>Normalization: {escape(str(l1_summary.get('normalization', {}).get('method', 'not available')))}</li>"
        f"<li>VST: {escape(str(l1_summary.get('vst', {}).get('method', 'not available')))}</li>"
        f"<li>DESeq2 model: <code>{escape(str(project['design']['formula']))}</code>; multiple testing: Benjamini-Hochberg.</li>"
        f"<li>DEG thresholds: padj &lt; {escape(str(thresholds['padj']))}; |log2FC| ≥ {escape(str(thresholds['abs_log2fc']))}.</li>"
        f"<li>R: {escape(str(package_versions.get('R', 'not available')))}; DESeq2: {escape(str(package_versions.get('DESeq2', 'not available')))}.</li>"
        f"<li>{_pipeline_provenance_line(execution, project.get('project', {}).get('pipeline') if isinstance(project.get('project'), dict) else None)}</li>"
        "<li>GO GSEA evaluates BP, MF, and CC; KEGG GSEA uses the configured KEGG provider.</li>"
        "</ul>",
        "<p>This automated report contains technical/statistical results only and no biological interpretation.</p>",
        "</body></html>",
    ])
    (output / "report.html").write_text("\n".join(sections), encoding="utf-8")


def enrichment_config(contract_path: Path, inputs: Path, l2: Path, kind: str, output: Path) -> dict[str, Any]:
    """Materialize the existing R-backend config without reimplementing its maths."""

    contract = _read_contract(contract_path)
    inputs_root, manifest = _read_execution_inputs(inputs)
    annotation = _required_frozen_annotation(contract, kind)
    organism = annotation["organism"]
    orgdb = {"Homo sapiens": "org.Hs.eg.db", "Mus musculus": "org.Mm.eg.db"}.get(organism)
    if orgdb is None:
        raise ValueError(f"unsupported enrichment organism: {organism}")
    contrasts = []
    for item in _contrasts(_staged_file(inputs_root, manifest.get("contrasts"), "contrasts")):
        root = l2 / "contrasts" / item["contrast_id"]
        contrasts.append({
            "contrast_id": item["contrast_id"],
            "all_genes": str(root / "all_genes.tsv"),
        })
    names = {"gsea-go": "gsea_go", "gsea-kegg": "gsea_kegg"}
    cfg: dict[str, Any] = {"output_dir": str(output / "enrichment" / names[kind]), "annotation": annotation, "orgdb_package": orgdb, "contrasts": contrasts}
    if kind == "gsea-kegg":
        code = {"Homo sapiens": "hsa", "Mus musculus": "mmu"}[organism]
        cfg.update({"mode": "gsea", "kegg": {"organism_code": code, "provider": annotation["enrichment"]["kegg"]["resource_provider"], "probe_endpoint": f"https://rest.kegg.jp/list/pathway/{code}"}})
    return cfg


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("l1-config", "l2-config", "enrichment-config", "report"))
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--l1", type=Path)
    parser.add_argument("--l2", type=Path)
    parser.add_argument("--enrichment", type=Path, nargs="*")
    parser.add_argument("--kind")
    args = parser.parse_args()
    if args.action == "l1-config":
        value = l1_config(args.contract, args.inputs, Path("l1"))
        args.out.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif args.action == "l2-config":
        if args.l1 is None:
            parser.error("--l1 is required for l2-config")
        value = l2_config(args.contract, args.inputs, args.l1, Path("l2"))
        args.out.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif args.action == "enrichment-config":
        if args.l2 is None or args.kind is None:
            parser.error("--l2 and --kind are required for enrichment-config")
        value = enrichment_config(args.contract, args.inputs, args.l2, args.kind, Path("."))
        args.out.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        if args.l1 is None:
            parser.error("--l1 is required for report")
        report(args.contract, args.inputs, args.l1, args.l2, args.out, args.enrichment or [])


if __name__ == "__main__":
    main()
