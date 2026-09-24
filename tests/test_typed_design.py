from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app
from rnaseq.l2 import prepare_l2
from rnaseq.validators import validate_project


runner = CliRunner()


def _config(formula: str, variables: dict[str, str]) -> dict:
    config = deepcopy(base_config())
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": formula, "variables": variables}
    return config


def _codes(report) -> set[str]:
    return {issue.code for issue in report.errors}


def test_explicit_categorical_regression_and_legacy_inference(project_factory):
    explicit = validate_project(project_factory(config=_config("~ condition", {"condition": "categorical"})))
    legacy = validate_project(project_factory())
    assert explicit.is_valid and legacy.is_valid
    assert dict(explicit.design_variable_types) == {"condition": "categorical"}
    assert dict(legacy.design_variable_types) == {"condition": "categorical"}


@pytest.mark.parametrize(
    "variables",
    (None, {"batch": "categorical"}),
)
def test_schema_13_requires_complete_explicit_formula_variable_types(project_factory, variables):
    config = deepcopy(base_config())
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": "~ batch + condition"}
    if variables is not None:
        config["design"]["variables"] = variables
    report = validate_project(project_factory(config=config))
    assert "invalid_project_config" in _codes(report)
    assert "design.variables" in report.errors[0].message


def test_schema_13_numeric_looking_batch_declared_categorical_remains_categorical(project_factory):
    metadata = """sample_id,batch,condition
C1,1,Control
C2,01,Control
C3,2,Control
T1,1,Treatment
T2,01,Treatment
T3,2,Treatment
"""
    config = _config("~ batch + condition", {"batch": "categorical", "condition": "categorical"})
    report = validate_project(project_factory(config=config, metadata=metadata))
    assert report.is_valid, report.errors
    assert dict(report.design_variable_types) == {"batch": "categorical", "condition": "categorical"}


def test_schema_12_without_design_variables_keeps_legacy_inference(project_factory):
    metadata = """sample_id,batch,condition
C1,1,Control
C2,2,Control
C3,3,Control
T1,1,Treatment
T2,2,Treatment
T3,3,Treatment
"""
    config = deepcopy(base_config())
    config["schema_version"] = "1.2"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": "~ batch + condition"}
    report = validate_project(project_factory(config=config, metadata=metadata))
    assert report.is_valid, report.errors
    assert dict(report.design_variable_types) == {"batch": "continuous", "condition": "categorical"}


def _inferred_warnings(report) -> list[str]:
    return [issue.message for issue in report.warnings if issue.code == "inferred_continuous_design_variable"]


def test_legacy_inferred_continuous_covariate_warns_once_without_changing_type(project_factory):
    metadata = """sample_id,batch,condition
C1,1,Control
C2,2,Control
C3,3,Control
T1,1,Treatment
T2,2,Treatment
T3,3,Treatment
"""
    config = deepcopy(base_config())
    config["design"] = {"type": "two_group", "formula": "~ batch + condition"}
    report = validate_project(project_factory(config=config, metadata=metadata))
    assert report.is_valid, report.errors
    assert dict(report.design_variable_types) == {"batch": "continuous", "condition": "categorical"}
    messages = _inferred_warnings(report)
    assert len(messages) == 1
    assert "'batch' was inferred as continuous" in messages[0]


def test_legacy_categorical_inference_and_numeric_contrast_factor_do_not_warn(project_factory):
    labelled = """sample_id,batch,condition
C1,B1,0
C2,B2,0
C3,B1,0
T1,B2,1
T2,B1,1
T3,B2,1
"""
    contrasts = "contrast_id,factor,numerator,denominator\ntreated_vs_control,condition,1,0\n"
    config = deepcopy(base_config())
    config["design"] = {"type": "two_group", "formula": "~ batch + condition"}
    report = validate_project(project_factory(config=config, metadata=labelled, contrasts=contrasts))
    assert report.is_valid, report.errors
    assert dict(report.design_variable_types) == {"batch": "categorical", "condition": "categorical"}
    assert _inferred_warnings(report) == []


def test_explicit_schema_13_continuous_type_does_not_emit_inferred_warning(project_factory):
    metadata = """sample_id,age,condition
C1,30,Control
C2,34,Control
C3,32,Control
T1,31,Treatment
T2,35,Treatment
T3,39,Treatment
"""
    config = _config("~ age + condition", {"age": "continuous", "condition": "categorical"})
    report = validate_project(project_factory(config=config, metadata=metadata))
    assert report.is_valid, report.errors
    assert dict(report.design_variable_types) == {"age": "continuous", "condition": "categorical"}
    assert _inferred_warnings(report) == []


def test_legacy_numeric_looking_contrast_levels_remain_categorical(project_factory):
    metadata = """sample_id,condition
C1,0
C2,0
C3,0
T1,1
T2,1
T3,1
"""
    contrasts = "contrast_id,factor,numerator,denominator\ntreated_vs_control,condition,1,0\n"
    report = validate_project(project_factory(metadata=metadata, contrasts=contrasts))
    assert report.is_valid
    assert report.design_variable_types["condition"] == "categorical"


def test_continuous_covariate_is_typed_and_preserves_decimal_metadata(project_factory):
    metadata = """sample_id,age,condition
C1,30.25,Control
C2,31.5,Control
C3,33.75,Control
T1,29.5,Treatment
T2,35.125,Treatment
T3,38.5,Treatment
"""
    report = validate_project(project_factory(config=_config("~ age + condition", {"age": "continuous", "condition": "categorical"}), metadata=metadata))
    assert report.is_valid
    assert dict(report.design_variable_types) == {"age": "continuous", "condition": "categorical"}
    prepared = prepare_l2(report, run_id=None)
    assert prepared.l1.config["design_variable_types"]["age"] == "continuous"
    assert prepared.metadata_rows[0]["age"] == "30.25"


def test_batch_and_multiple_covariates_are_full_rank(project_factory):
    metadata = """sample_id,batch,age,condition
C1,B1,30,Control
C2,B1,34,Control
C3,B2,32,Control
T1,B1,31,Treatment
T2,B2,35,Treatment
T3,B2,39,Treatment
"""
    report = validate_project(project_factory(config=_config("~ batch + age + condition", {"batch": "categorical", "age": "continuous", "condition": "categorical"}), metadata=metadata))
    assert report.is_valid
    assert dict(report.design_variable_types) == {"batch": "categorical", "age": "continuous", "condition": "categorical"}


def test_paired_regression_remains_categorical(project_factory):
    config = _config("~ pair_id + condition", {"pair_id": "categorical", "condition": "categorical"})
    config["design"].update({"type": "paired_two_group", "pair_id": "pair_id"})
    metadata = """sample_id,pair_id,condition
C1,P1,Control
T1,P1,Treatment
C2,P2,Control
T2,P2,Treatment
C3,P3,Control
T3,P3,Treatment
"""
    counts = "gene_id,C1,T1,C2,T2,C3,T3\nGeneA,10,15,11,16,9,14\nGeneB,30,25,29,24,32,27\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata))
    assert report.is_valid
    assert report.design_variable_types["pair_id"] == "categorical"


def test_full_rank_paired_batch_design_is_valid(project_factory):
    config = _config(
        "~ pair_id + batch + condition",
        {"pair_id": "categorical", "batch": "categorical", "condition": "categorical"},
    )
    config["design"].update({"type": "paired_two_group", "pair_id": "pair_id"})
    metadata = """sample_id,pair_id,batch,condition
C1,P1,B1,Control
T1,P1,B2,Treatment
C2,P2,B2,Control
T2,P2,B1,Treatment
C3,P3,B1,Control
T3,P3,B2,Treatment
"""
    counts = "gene_id,C1,T1,C2,T2,C3,T3\nGeneA,10,15,11,16,9,14\nGeneB,30,25,29,24,32,27\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata))
    assert report.is_valid, report.errors
    prepared = prepare_l2(report, run_id=None)
    assert prepared.config.design.formula == "~ pair_id + batch + condition"


@pytest.mark.parametrize("value", ["bad", "NA", "NaN", "Inf", "-Inf", ""])
def test_invalid_continuous_values_are_rejected(project_factory, value):
    metadata = f"""sample_id,age,condition
C1,30,Control
C2,{value},Control
C3,32,Control
T1,33,Treatment
T2,34,Treatment
T3,35,Treatment
"""
    report = validate_project(project_factory(config=_config("~ age + condition", {"age": "continuous", "condition": "categorical"}), metadata=metadata))
    assert "invalid_continuous_value" in _codes(report) or "missing_design_value" in _codes(report)


def test_zero_variance_continuous_variable_is_rejected(project_factory):
    metadata = """sample_id,age,condition
C1,42,Control
C2,42,Control
C3,42,Control
T1,42,Treatment
T2,42,Treatment
T3,42,Treatment
"""
    report = validate_project(project_factory(config=_config("~ age + condition", {"age": "continuous", "condition": "categorical"}), metadata=metadata))
    assert "zero_variance_continuous_variable" in _codes(report)


def test_continuous_variable_cannot_be_a_contrast_factor(project_factory):
    metadata = """sample_id,age,condition
C1,30,Control
C2,31,Control
C3,32,Control
T1,33,Treatment
T2,34,Treatment
T3,35,Treatment
"""
    contrasts = "contrast_id,factor,numerator,denominator\nage_effect,age,35,30\n"
    report = validate_project(project_factory(config=_config("~ age + condition", {"age": "continuous", "condition": "categorical"}), metadata=metadata, contrasts=contrasts))
    assert "continuous_contrast_factor" in _codes(report)


def test_perfectly_confounded_batch_and_condition_are_rejected(project_factory):
    metadata = """sample_id,batch,condition
C1,B1,Control
C2,B1,Control
C3,B1,Control
T1,B2,Treatment
T2,B2,Treatment
T3,B2,Treatment
"""
    report = validate_project(project_factory(config=_config("~ batch + condition", {"batch": "categorical", "condition": "categorical"}), metadata=metadata))
    assert "rank_deficient_design" in _codes(report)
    assert any("'batch' and 'condition' are confounded" in issue.message for issue in report.errors)


def test_declared_variables_must_match_the_formula(project_factory):
    report = validate_project(project_factory(config=_config("~ condition", {"batch": "categorical", "condition": "categorical"})))
    assert "declared_design_variable_not_in_formula" in _codes(report)


def test_new_cli_writes_declared_continuous_covariates(tmp_path: Path):
    counts = tmp_path / "counts.csv"
    metadata = tmp_path / "metadata.csv"
    contrasts = tmp_path / "contrasts.csv"
    counts.write_text("gene_id,C1,C2,C3,T1,T2,T3\nGeneA,10,12,9,40,45,43\nGeneB,100,110,98,95,102,99\n", encoding="utf-8")
    metadata.write_text("sample_id,batch,age,condition\nC1,B1,30,Control\nC2,B1,34,Control\nC3,B2,32,Control\nT1,B1,31,Treatment\nT2,B2,35,Treatment\nT3,B2,39,Treatment\n", encoding="utf-8")
    contrasts.write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    result = runner.invoke(app, [
        "new", "--name", "typed", "--destination", str(tmp_path), "--species", "mouse",
        "--input-type", "raw_counts", "--preset", "L2", "--design-type", "two_group",
        "--counts", str(counts), "--metadata", str(metadata), "--contrasts", str(contrasts),
        "--condition-column", "condition", "--covariate", "batch", "--continuous-covariate", "age", "--yes",
    ])
    assert result.exit_code == 0, result.output
    report = validate_project(tmp_path / "typed")
    assert report.is_valid
    assert report.config.design.formula == "~ batch + age + condition"
    assert report.config.design.variables is not None
    assert {name: kind.value for name, kind in report.config.design.variables.items()} == {
        "batch": "categorical", "age": "continuous", "condition": "categorical"
    }
