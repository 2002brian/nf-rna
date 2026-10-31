from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from conftest import BASE_COUNTS, BASE_METADATA, base_config, pdf_text_origins, require_downstream_r_runtime
from rnaseq.downstream import L1Result
from rnaseq.errors import DownstreamExecutionError
from rnaseq.l2 import execute_l2, prepare_l2
from rnaseq.validators import validate_project


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _example_copy(tmp_path: Path) -> Path:
    source = Path(__file__).parents[1] / "examples" / "simple_two_group"
    destination = tmp_path / "replicated"
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("downstream", "planning"))
    return destination


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _typed_design_project(
    tmp_path: Path, *, formula: str, variables: dict[str, str], paired: bool = False,
) -> Path:
    """Create a deterministic real-DESeq2 fixture from the 100-gene example."""

    root = _example_copy(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    design: dict[str, object] = {"type": "paired_two_group" if paired else "two_group", "formula": formula, "variables": variables}
    if paired:
        design["pair_id"] = "pair_id"
        (root / "metadata.csv").write_text(
            "sample_id,pair_id,condition\n"
            "C1,P1,Control\nC2,P2,Control\nC3,P3,Control\n"
            "T1,P1,Treatment\nT2,P2,Treatment\nT3,P3,Treatment\n",
            encoding="utf-8",
        )
    config["design"] = design
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return root


def test_l2_scientific_provenance_records_actual_typed_model(tmp_path):
    require_downstream_r_runtime("jsonlite", "DESeq2", "ggplot2", "pheatmap")
    root = _typed_design_project(
        tmp_path,
        formula="~ batch + age + condition",
        variables={"batch": "categorical", "age": "continuous", "condition": "categorical"},
    )
    result = execute_l2(prepare_l2(validate_project(root), run_id=None))
    provenance = json.loads((result.output_dir / "scientific_provenance.json").read_text())
    design = provenance["design"]
    assert design["formula"] == "~ batch + age + condition"
    assert design["declared_variable_types"] == {
        "batch": "categorical", "age": "continuous", "condition": "categorical",
    }
    assert design["variables"]["batch"] == {"type": "categorical", "levels": ["B1", "B2"]}
    assert design["variables"]["condition"] == {"type": "categorical", "levels": ["Control", "Treatment"]}
    assert design["variables"]["age"] == {
        "type": "continuous", "n": 6, "min": 8, "max": 10, "mean": pytest.approx(8.666666666666666, abs=1e-4),
    }
    assert design["model_matrix"] == {
        "rank": 4,
        "column_count": 4,
        "full_rank": True,
        "column_names": ["(Intercept)", "batchB2", "age", "conditionTreatment"],
    }
    backend = json.loads((result.output_dir / "contrasts" / "Treatment_vs_Control" / "backend_summary.json").read_text())
    assert provenance["result"]["fit_method"] == backend["fit_method"]


def test_real_l2_three_batch_levels_emits_only_the_biological_contrast(tmp_path):
    require_downstream_r_runtime("jsonlite", "DESeq2", "ggplot2", "pheatmap")
    root = _example_copy(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {
        "type": "two_group",
        "formula": "~ batch + condition",
        "variables": {"batch": "categorical", "condition": "categorical"},
    }
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    (root / "metadata.csv").write_text(
        "sample_id,batch,condition\n"
        "C1,B1,Control\nC2,B2,Control\nC3,B3,Control\n"
        "T1,B1,Treatment\nT2,B2,Treatment\nT3,B3,Treatment\n",
        encoding="utf-8",
    )

    result = execute_l2(prepare_l2(validate_project(root), run_id=None))
    contrast_dirs = sorted(path.name for path in (result.output_dir / "contrasts").iterdir() if path.is_dir())
    assert contrast_dirs == ["Treatment_vs_Control"]
    assert (result.output_dir / "contrasts" / "Treatment_vs_Control" / "volcano.png").is_file()


@pytest.mark.parametrize(
    ("formula", "variables", "paired", "expected_columns"),
    (
        ("~ condition", {"condition": "categorical"}, False, ["(Intercept)", "conditionTreatment"]),
        ("~ pair_id + condition", {"pair_id": "categorical", "condition": "categorical"}, True, ["(Intercept)", "pair_idP2", "pair_idP3", "conditionTreatment"]),
    ),
)
def test_l2_scientific_provenance_covers_simple_and_paired_models(
    tmp_path, formula, variables, paired, expected_columns,
):
    require_downstream_r_runtime("jsonlite", "DESeq2", "ggplot2", "pheatmap")
    root = _typed_design_project(tmp_path, formula=formula, variables=variables, paired=paired)
    result = execute_l2(prepare_l2(validate_project(root), run_id=None))
    design = json.loads((result.output_dir / "scientific_provenance.json").read_text())["design"]
    assert design["formula"] == formula
    assert design["model_matrix"] == {
        "rank": len(expected_columns),
        "column_count": len(expected_columns),
        "full_rank": True,
        "column_names": expected_columns,
    }
    assert design["variables"]["condition"] == {"type": "categorical", "levels": ["Control", "Treatment"]}
    if paired:
        assert design["pair_id"] == "pair_id"
        assert design["variables"]["pair_id"] == {"type": "categorical", "levels": ["P1", "P2", "P3"]}


def test_real_l2_replicated_raw_counts_outputs_and_determinism(tmp_path):
    require_downstream_r_runtime("jsonlite", "DESeq2", "ggplot2", "pheatmap")
    root = _example_copy(tmp_path)
    prepared = prepare_l2(validate_project(root), run_id=None)
    first = execute_l2(prepared)
    assert first.status == "SUCCESS"
    contrast = first.output_dir / "contrasts" / "Treatment_vs_Control"
    all_rows = _rows(contrast / "all_genes.tsv")
    significant, up, down = (_rows(contrast / name) for name in ("significant.tsv", "upregulated.tsv", "downregulated.tsv"))
    assert len(all_rows) == 100 and len(significant) == 20 and len(up) == len(down) == 10
    assert {row["gene_id"] for row in up} | {row["gene_id"] for row in down} == {row["gene_id"] for row in significant}
    assert not ({row["gene_id"] for row in up} & {row["gene_id"] for row in down})
    by_gene = {row["gene_id"]: row for row in all_rows}
    assert float(by_gene["Gene001"]["log2FoldChange"]) > 0
    assert float(by_gene["Gene011"]["log2FoldChange"]) < 0
    assert (contrast / "volcano.png").stat().st_size > 0
    assert (contrast / "volcano.tiff").stat().st_size > 0
    assert (contrast / "heatmap.png").stat().st_size > 0
    assert (contrast / "heatmap.tiff").stat().st_size > 0
    assert not list(contrast.glob("*.svg"))
    assert {path.name for path in contrast.glob("*.pdf")} == {"volcano.pdf", "ma_plot.pdf", "heatmap.pdf"}
    assert {"volcano_300dpi.png", "ma_plot_300dpi.png", "heatmap_300dpi.png"} <= {path.name for path in contrast.glob("*.png")}
    shrinkage = json.loads((contrast / "contrast_metadata.json").read_text())["shrinkage"]
    assert shrinkage["status"] == "success" and shrinkage["method"] == "apeglm"
    assert all(row["log2FoldChange_shrunk"] not in {"NA", row["log2FoldChange"]} for row in all_rows)
    heatmap_genes = [row["gene_id"] for row in _rows(contrast / "heatmap_genes.tsv")]
    assert set(heatmap_genes) == {row["gene_id"] for row in significant}
    heatmap = next(item for item in json.loads((first.output_dir / "figure_manifest.json").read_text())["figures"] if item["id"].endswith("/heatmap"))
    assert heatmap["title"] == "Treatment vs Control\nTop 20 DEGs \u00b7 row z-score of VST"
    width_pt = next(item for item in heatmap["files"] if item["format"] == "pdf")["width_pt"]
    assert all(0 <= x <= width_pt for x, _ in pdf_text_origins(contrast / "heatmap.pdf"))
    summary = yaml.safe_load((contrast / "summary.yaml").read_text())
    assert summary["independent_filtering"] is True
    assert summary["independent_filtering_alpha"] == 0.05  # results(alpha) is the declared DEG padj threshold
    filtering = json.loads((contrast / "contrast_metadata.json").read_text())["deg_definition"]["independent_filtering"]
    assert filtering["alpha"] == 0.05 and isinstance(filtering["filter_threshold"], (int, float))
    assert "fallback" in summary["fit_method"]
    # R warnings are recorded (verbatim, by step) without changing the analysis outcome.
    warnings = json.loads((contrast / "contrast_metadata.json").read_text())["r_warnings"]
    recorded = [*warnings["fit"], *warnings["contrast"]]
    assert warnings["status"] == ("WARNINGS" if recorded else "NONE")
    assert all({"step", "message", "call", "count"} <= set(item) and item["count"] >= 1 for item in recorded)
    backend = json.loads((contrast / "backend_summary.json").read_text())
    assert backend["r_warning_status"] == warnings["status"] and backend["r_warning_count"] == sum(item["count"] for item in recorded)
    assert "biological interpretation" in (contrast / "contrast_report.md").read_text()
    assert json.loads((first.output_dir / "l2_state.json").read_text())["status"] == "SUCCESS"
    before = {_path.name: _digest(_path) for _path in (contrast / "all_genes.tsv", contrast / "significant.tsv", contrast / "upregulated.tsv", contrast / "downregulated.tsv", contrast / "summary.yaml", contrast / "heatmap_genes.tsv", contrast / "heatmap_zscores.tsv")}
    execute_l2(prepared)
    assert before == {_path.name: _digest(_path) for _path in (contrast / "all_genes.tsv", contrast / "significant.tsv", contrast / "upregulated.tsv", contrast / "downregulated.tsv", contrast / "summary.yaml", contrast / "heatmap_genes.tsv", contrast / "heatmap_zscores.tsv")}


def test_l2_blocks_n1_inference_without_affecting_l1_state(project_factory):
    config = base_config()
    counts = "gene_id,C1,T1\nGeneA,10,20\nGeneB,30,15\n"
    metadata = "sample_id,condition\nC1,Control\nT1,Treatment\n"
    root = project_factory(config=config, counts=counts, metadata=metadata)
    report = validate_project(root)
    assert report.is_valid
    with pytest.raises(DownstreamExecutionError, match="requires biological replication"):
        prepare_l2(report, run_id=None)
    assert not (root / "downstream" / "l2" / "l2_state.json").exists()


def test_l2_rechecks_contrast_values_immediately_before_inference(project_factory):
    root = project_factory()
    report = validate_project(root)
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nBad,condition,Absent,Control\n", encoding="utf-8")
    with pytest.raises(DownstreamExecutionError, match="absent from metadata"):
        prepare_l2(report, run_id=None)


def test_l2_supports_multiple_named_contrasts_and_paired_metadata(project_factory):
    config = deepcopy(base_config())
    config["design"] = {"type": "multi_group", "formula": "~ condition"}
    metadata = "sample_id,condition\nC1,A\nC2,A\nC3,B\nT1,B\nT2,C\nT3,C\n"
    contrasts = "contrast_id,factor,numerator,denominator\nB_vs_A,condition,B,A\nC_vs_A,condition,C,A\n"
    prepared = prepare_l2(validate_project(project_factory(config=config, metadata=metadata, contrasts=contrasts)), run_id=None)
    assert [item.contrast_id for item in prepared.contrasts] == ["B_vs_A", "C_vs_A"]

    paired = deepcopy(base_config())
    paired["design"] = {"type": "paired_two_group", "formula": "~ subject_id + condition", "pair_id": "subject_id"}
    metadata = "sample_id,subject_id,condition\nC1,S1,Control\nC2,S2,Control\nC3,S3,Control\nT1,S1,Treatment\nT2,S2,Treatment\nT3,S3,Treatment\n"
    paired_prepared = prepare_l2(validate_project(project_factory(config=paired, metadata=metadata)), run_id=None)
    assert paired_prepared.config.design.formula == "~ subject_id + condition"


def test_l2_zero_deg_heatmap_is_not_applicable(tmp_path):
    require_downstream_r_runtime("jsonlite", "DESeq2", "ggplot2", "pheatmap")
    root = _example_copy(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text())
    config["thresholds"]["abs_log2fc"] = 100.0
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    result = execute_l2(prepare_l2(validate_project(root), run_id=None))
    contrast = result.output_dir / "contrasts" / "Treatment_vs_Control"
    assert _rows(contrast / "significant.tsv") == []
    assert not (contrast / "heatmap.png").exists()
    assert yaml.safe_load((contrast / "summary.yaml").read_text())["heatmap_status"] == "NOT_APPLICABLE"


def test_l2_backend_failure_preserves_separate_failed_state(monkeypatch, project_factory):
    prepared = prepare_l2(validate_project(project_factory()), run_id=None)
    monkeypatch.setattr("rnaseq.l2._require_r", lambda: None)
    monkeypatch.setattr("rnaseq.l2.execute_l1", lambda _prepared: L1Result(prepared.l1.output_dir, prepared.l1.output_dir / "l1_state.json", {"package_versions": {}}))
    monkeypatch.setattr("rnaseq.l2.subprocess.run", lambda *_args, **_kwargs: SimpleNamespace(returncode=91))
    with pytest.raises(DownstreamExecutionError, match="return code 91"):
        execute_l2(prepared)
    state = json.loads((prepared.output_dir / "l2_state.json").read_text())
    assert state["status"] == "FAILED"


def test_production_figures_are_exported_only_through_the_central_figure_policy():
    """Every backend draws through figures.R: PNG report preview, 300-dpi LZW TIFF and Cairo PDF; never SVG."""

    root = Path(__file__).parents[1] / "src" / "rnaseq" / "r"
    for name in ("l1_analysis.R", "l2_analysis.R", "go_analysis.R", "gsea_analysis.R", "kegg_analysis.R"):
        text = (root / name).read_text(encoding="utf-8")
        assert '"figures.R"' in text and "provenance.R" in text
        assert any(call in text for call in ("nf_rna_save_figure(", "nf_rna_ora_dotplot(", "nf_rna_gsea_dotplot(")) and "nf_rna_write_figure_manifest(" in text
        for direct in ("ggsave(", "png(", "tiff(", "pdf(", ".svg", "theme_minimal("):
            assert direct not in text, (name, direct)
    for script in sorted(root.glob("*.R")):
        if script.name == "figures.R":
            continue
        code = "\n".join(line.split("#", 1)[0] for line in script.read_text(encoding="utf-8").splitlines())
        assert not re.search(r"\b(?:ggsave|png|tiff|jpeg|bmp|svg|pdf|cairo_pdf|postscript|dev\.copy|dev\.print|savePlot)\s*\(", code), script.name
    policy = (root / "figures.R").read_text(encoding="utf-8")
    assert 'formats=c("png", "tiff", "pdf")' in policy
    assert 'compression="lzw", bg="white"' in policy and "cairo_pdf(" in policy
    assert "general=list(" in policy and "raster_dpi=300" in policy.split("nature=list(")[0]


def _imbalanced_batch_project(tmp_path: Path, name: str, formula: str, variables: dict[str, str]) -> Path:
    """Deterministic 8-sample dataset: genes 1-60 are higher in batch B2 (no condition effect),
    genes 61-80 are truly higher in Treatment; batch is imbalanced against condition (1/3 vs 3/1)."""

    import math
    import random

    rng = random.Random(7)

    def poisson(lam: float) -> int:
        if lam >= 30:
            return max(0, round(rng.gauss(lam, math.sqrt(lam))))
        limit, k, p = math.exp(-lam), 0, 1.0
        while True:
            p *= rng.random()
            if p <= limit:
                return k
            k += 1

    samples = ["C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"]
    condition = ["Control"] * 4 + ["Treatment"] * 4
    batch = ["B1", "B1", "B1", "B2", "B1", "B2", "B2", "B2"]
    lines = ["gene_id," + ",".join(samples)]
    for gene in range(1, 201):
        base = math.exp(rng.gauss(5, 1))
        row = []
        for index in range(8):
            mean = base * (4.0 if gene <= 60 and batch[index] == "B2" else 1.0) * (2.5 if 61 <= gene <= 80 and condition[index] == "Treatment" else 1.0)
            row.append(poisson(mean * rng.gammavariate(40, 1 / 40)))
        lines.append(f"Gene{gene}," + ",".join(map(str, row)))
    root = _example_copy(tmp_path / name)
    (root / "counts.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (root / "metadata.csv").write_text(
        "sample_id,condition,batch\n" + "".join(f"{s},{c},{b}\n" for s, c, b in zip(samples, condition, batch)), encoding="utf-8")
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": formula, "variables": variables}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return root


def test_batch_term_changes_the_fitted_model_and_the_condition_contrast(tmp_path):
    require_downstream_r_runtime("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    results = {}
    for name, formula, variables in (
        ("plain", "~ condition", {"condition": "categorical"}),
        ("adjusted", "~ batch + condition", {"batch": "categorical", "condition": "categorical"}),
    ):
        root = _imbalanced_batch_project(tmp_path, name, formula, variables)
        output = execute_l2(prepare_l2(validate_project(root), run_id=None)).output_dir
        design = json.loads((output / "scientific_provenance.json").read_text(encoding="utf-8"))["design"]
        table = {row["gene_id"]: float(row["log2FoldChange"]) for row in _rows(output / "contrasts" / "Treatment_vs_Control" / "all_genes.tsv")}
        results[name] = (design["model_matrix"]["column_names"], table, json.loads((output / "backend_config.json").read_text())["formula"])
        # Plotted/QC values are never batch-corrected: blind VST does not depend on the design.
        results[name] += (_digest(output.parent / "l1" / "vst.tsv"),)
    assert results["plain"][0] == ["(Intercept)", "conditionTreatment"] and results["plain"][2] == "~ condition"
    assert results["adjusted"][0] == ["(Intercept)", "batchB2", "conditionTreatment"] and results["adjusted"][2] == "~ batch + condition"
    assert results["plain"][3] == results["adjusted"][3]

    plain, adjusted = results["plain"][1], results["adjusted"][1]
    batch_genes = [f"Gene{i}" for i in range(1, 61)]
    median = lambda values: sorted(values)[len(values) // 2]
    # Batch-only genes look condition-responsive without adjustment and lose most of that apparent effect with it.
    assert median([abs(plain[g]) for g in batch_genes]) > 2 * median([abs(adjusted[g]) for g in batch_genes])
    assert sum(abs(plain[g] - adjusted[g]) > 0.01 for g in plain) > 0.9 * len(plain)


def test_l1_boxplot_samples_follow_sample_order_not_metadata_row_order(tmp_path):
    require_downstream_r_runtime("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    root = _imbalanced_batch_project(tmp_path, "shuffled", "~ batch + condition", {"batch": "categorical", "condition": "categorical"})
    lines = (root / "metadata.csv").read_text(encoding="utf-8").splitlines()
    (root / "metadata.csv").write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n", encoding="utf-8")
    prepared = prepare_l2(validate_project(root), run_id=None)
    output = execute_l2(prepared).output_dir.parent / "l1"
    rows = _rows(output / "sample_distribution_qc.tsv")
    assert [row["sample_id"] for row in rows if row["matrix"] == "blind_vst"] == list(prepared.l1.sample_ids)
    manifest = json.loads((output / "figure_manifest.json").read_text(encoding="utf-8"))
    box = next(item for item in manifest["figures"] if item["id"] == "boxplot_vst")
    assert box["encoding"]["fill"] == "condition" and box["encoding"]["batch_marker"] == "black marker shape below each box"
