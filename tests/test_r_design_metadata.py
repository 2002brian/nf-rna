"""Real-R checks that frozen design types reach DESeq2 with textual identity intact."""

from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime
from pathlib import Path

from conftest import base_config, require_r_packages, require_rscript
from rnaseq.planner import generate_plan
from rnaseq.service import create_case_run, freeze_case_inputs, resolve_downstream_inputs
from rnaseq.validators import validate_project
from rnaseq.workflow_support import l1_config, l2_config


ROOT = Path(__file__).resolve().parents[1]
R_SCRIPTS = ROOT / "src" / "rnaseq" / "r"
COUNTS = (ROOT / "examples" / "simple_two_group" / "counts.csv").read_text(encoding="utf-8")
CONTRASTS = "contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n"
R_PACKAGES = ("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap")


def _metadata(column: str, values: list[str]) -> str:
    samples = ("C1", "C2", "C3", "T1", "T2", "T3")
    conditions = ("Control",) * 3 + ("Treatment",) * 3
    rows = "".join(f"{sample},{value},{condition}\n" for sample, value, condition in zip(samples, values, conditions))
    return f"sample_id,{column},condition\n{rows}"


def _run_r(script: str, config: dict, work: Path) -> subprocess.CompletedProcess[str]:
    path = work / f"{script}.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    return subprocess.run(
        ["Rscript", str(R_SCRIPTS / f"{script}.R"), "--config", str(path)], cwd=work, capture_output=True, text=True, check=False,
    )


def _frozen_configs(project_factory, tmp_path: Path, *, formula: str, variables: dict[str, str], metadata: str):
    config = base_config()
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": formula, "variables": variables}
    validation = validate_project(project_factory(config=config, counts=COUNTS, metadata=metadata, contrasts=CONTRASTS))
    assert validation.is_valid, validation.errors
    generate_plan(validation)
    run = create_case_run(validation, "CASE-R-DESIGN", moment=datetime(2026, 9, 23, 12, 0, 0))
    contract = freeze_case_inputs(validation, run, profile="local", command=["rnaseq", "run"]).contract
    inputs = resolve_downstream_inputs(run).root
    work = tmp_path / "task"
    work.mkdir()
    return work, l1_config(contract, inputs, work / "l1"), l2_config(contract, inputs, work / "l1", work / "l2")


def _run_l1_l2(work: Path, l1: dict, l2: dict) -> dict:
    for script, config in (("l1_analysis", l1), ("l2_analysis", l2)):
        result = _run_r(script, config, work)
        assert result.returncode == 0, result.stderr
    return json.loads((work / "l2" / "scientific_provenance.json").read_text(encoding="utf-8"))["design"]


def test_real_r_keeps_leading_zero_categorical_levels_as_distinct_dummies(project_factory, tmp_path):
    require_r_packages(*R_PACKAGES)
    work, l1, l2 = _frozen_configs(
        project_factory, tmp_path, formula="~ batch + condition",
        variables={"batch": "categorical", "condition": "categorical"},
        metadata=_metadata("batch", ["1", "01", "2", "1", "01", "2"]),
    )
    design = _run_l1_l2(work, l1, l2)
    assert design["variables"]["batch"] == {"type": "categorical", "levels": ["01", "1", "2"]}
    assert design["model_matrix"] == {
        "rank": 4, "column_count": 4, "full_rank": True,
        "column_names": ["(Intercept)", "batch1", "batch2", "conditionTreatment"],
    }

    # The fitted model must equal an independent fit of the intended three-level
    # factor, and must not equal the collapsed two-level model of v1.2.0.
    reference = subprocess.run(["Rscript", "-e", f"""
suppressPackageStartupMessages(library(DESeq2))
m <- read.csv({json.dumps(l2["metadata"])}, colClasses="character"); rownames(m) <- m$sample_id
x <- read.csv({json.dumps(l2["counts"])}, check.names=FALSE); rownames(x) <- x$gene_id; x <- as.matrix(x[, m$sample_id]); x <- x[rowSums(x) >= 10, ]
fit <- function(batch) {{ md <- data.frame(batch=batch, condition=factor(m$condition), row.names=m$sample_id)
  dds <- DESeqDataSetFromMatrix(x, md, ~ batch + condition)
  # Same fit procedure as l2_analysis.R, including its documented fallback.
  fitted <- tryCatch(suppressMessages(DESeq(dds, quiet=TRUE)), error=function(e) e)
  if (inherits(fitted, "error")) {{
    if (!grepl("all gene-wise dispersion estimates", conditionMessage(fitted), fixed=TRUE)) stop(fitted)
    dds <- estimateDispersionsGeneEst(estimateSizeFactors(dds)); dispersions(dds) <- mcols(dds)$dispGeneEst
    fitted <- suppressMessages(nbinomWaldTest(dds, quiet=TRUE))
  }}
  results(fitted, contrast=c("condition","Treatment","Control"))$log2FoldChange }}
cat(sprintf("%.15g", fit(factor(m$batch, levels=c("01","1","2")))), sep="\\n"); cat("---\\n")
cat(sprintf("%.15g", fit(factor(as.numeric(m$batch)))), sep="\\n")
"""], capture_output=True, text=True, check=False)
    assert reference.returncode == 0, reference.stderr
    intended, collapsed = (
        [float(value) for value in block.split()] for block in reference.stdout.split("---")
    )
    with (work / "l2" / "contrasts" / "Treatment_vs_Control" / "all_genes.tsv").open(encoding="utf-8") as handle:
        observed = [float(row["log2FoldChange"]) for row in csv.DictReader(handle, delimiter="\t")]
    assert max(abs(a - b) for a, b in zip(observed, intended)) < 1e-10
    assert max(abs(a - b) for a, b in zip(observed, collapsed)) > 1e-4


def test_real_r_continuous_covariate_from_numeric_text(project_factory, tmp_path):
    require_r_packages(*R_PACKAGES)
    work, l1, l2 = _frozen_configs(
        project_factory, tmp_path, formula="~ age + condition",
        variables={"age": "continuous", "condition": "categorical"},
        metadata=_metadata("age", ["8", "9.0", "8e0", " 9 ", "8.5", "10"]),
    )
    design = _run_l1_l2(work, l1, l2)
    assert design["variables"]["age"] == {"type": "continuous", "n": 6, "min": 8, "max": 10, "mean": 8.75}
    assert design["model_matrix"]["column_names"] == ["(Intercept)", "age", "conditionTreatment"]


def test_real_r_rejects_continuous_value_that_r_cannot_parse(project_factory, tmp_path):
    require_r_packages(*R_PACKAGES)
    # Python's float() accepts "1_000"; R must refuse it rather than use NA.
    work, l1, _l2 = _frozen_configs(
        project_factory, tmp_path, formula="~ age + condition",
        variables={"age": "continuous", "condition": "categorical"},
        metadata=_metadata("age", ["8", "9", "1_000", "9", "8", "10"]),
    )
    result = _run_r("l1_analysis", l1, work)
    assert result.returncode != 0
    assert "continuous design variable is not finite: age (sample C3: '1_000')" in result.stderr


def _design_metadata_levels(tmp_path: Path, *, metadata: str, types: str, formula: str, variable: str) -> str:
    require_rscript()
    path = tmp_path / "metadata.csv"
    path.write_text(metadata, encoding="utf-8")
    code = f"""
source({json.dumps(str(R_SCRIPTS / "design_metadata.R"))})
cfg <- list(metadata={json.dumps(str(path))}, formula={json.dumps(formula)}, design_variable_types={types})
m <- nf_rna_design_metadata(cfg, c("C1","C2","C3","T1","T2","T3"))
x <- m[[{json.dumps(variable)}]]
cat(class(x), if (is.factor(x)) levels(x) else x, sep="|")
"""
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_design_metadata_keeps_1_01_001_distinct(tmp_path):
    output = _design_metadata_levels(
        tmp_path, metadata=_metadata("batch", ["1", "01", "001", "1", "01", "001"]),
        types='list(batch="categorical", condition="categorical")', formula="~ batch + condition", variable="batch",
    )
    assert output == "factor|001|01|1"


def test_design_metadata_keeps_previous_categorical_level_order(tmp_path):
    numeric_coded = _design_metadata_levels(
        tmp_path, metadata=_metadata("batch", ["2", "10", "3", "2", "10", "3"]),
        types='list(batch="categorical", condition="categorical")', formula="~ batch + condition", variable="batch",
    )
    labelled = _design_metadata_levels(
        tmp_path, metadata=_metadata("batch", ["B2", "B10", "B1", "B2", "B10", "B1"]),
        types='list(batch="categorical", condition="categorical")', formula="~ batch + condition", variable="batch",
    )
    assert numeric_coded == "factor|2|3|10"
    assert labelled == "factor|B1|B10|B2"


def test_design_metadata_without_frozen_types_keeps_previous_inference(tmp_path):
    output = _design_metadata_levels(
        tmp_path, metadata=_metadata("age", ["8", "9", "8", "9", "8", "10"]),
        types="list()", formula="~ age + condition", variable="age",
    )
    assert output == "integer|8|9|8|9|8|10"
