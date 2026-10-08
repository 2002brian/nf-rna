"""Template-first lifecycle: `rnaseq new` -> fill metadata -> L1 -> `rnaseq design` -> L2."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config, require_r_packages
from rnaseq.cli import app
from rnaseq.design_setup import DesignSetupError, generate_contrasts, updated_project_yaml
from rnaseq.models import ProjectConfig
from rnaseq.validators import ContrastDefinition, validate_project

runner = CliRunner()
CONTRAST_HEADER = "contrast_id,factor,numerator,denominator\n"


def _new(tmp_path: Path, name: str, *extra: str, preset: str = "L1", design: str = "two_group", input_type: str = "raw_counts"):
    return runner.invoke(app, [
        "new", "--name", name, "--destination", str(tmp_path), "--species", "mouse", "--input-type", input_type,
        "--scaffold", "--preset", preset, "--design-type", design, "--yes", *extra,
    ])


def _config(root: Path) -> dict:
    return yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))


def _fill(root: Path, metadata: str, genes: int = 2) -> None:
    (root / "metadata.csv").write_text(metadata, encoding="utf-8")
    samples = [line.split(",")[0] for line in metadata.strip().splitlines()[1:]]
    rows = "".join(f"G{g}," + ",".join(str(10 + 3 * g + 7 * index % 11) for index, _ in enumerate(samples)) + "\n" for g in range(genes))
    (root / "input" / "counts.csv").write_text("gene_id," + ",".join(samples) + "\n" + rows, encoding="utf-8")


TWO_GROUP_BATCH = "sample_id,condition,batch\ns1,Control,A\ns2,Control,A\ns3,Control,B\ns4,Treatment,A\ns5,Treatment,B\ns6,Treatment,B\n"
THREE_GROUP = (
    "sample_id,condition\nc1,Control\nc2,Control\na1,DrugA\na2,DrugA\nb1,DrugB\nb2,DrugB\n"
)


# ---------------------------------------------------------------- templates


@pytest.mark.parametrize(("design", "extra", "batch", "header", "formula", "variables"), (
    ("two_group", [], False, "sample_id,condition", "~ condition", {"condition": "categorical"}),
    ("two_group", ["--batch"], True, "sample_id,condition,batch", "~ batch + condition", {"batch": "categorical", "condition": "categorical"}),
    ("multi_group", [], False, "sample_id,condition", "~ condition", {"condition": "categorical"}),
    ("multi_group", ["--batch"], True, "sample_id,condition,batch", "~ batch + condition", {"batch": "categorical", "condition": "categorical"}),
    ("paired_two_group", ["--pair-id", "patient"], False, "sample_id,patient,condition", "~ patient + condition", {"patient": "categorical", "condition": "categorical"}),
    ("paired_two_group", ["--pair-id", "subj", "--batch"], True, "sample_id,subj,condition,batch", "~ subj + batch + condition", {"subj": "categorical", "batch": "categorical", "condition": "categorical"}),
))
def test_new_generates_templates_and_provisional_model(tmp_path, design, extra, batch, header, formula, variables):
    result = _new(tmp_path, "proj", *extra, design=design)
    assert result.exit_code == 0, result.output
    root = tmp_path / "proj"
    assert (root / "metadata.csv").read_text(encoding="utf-8") == header + "\n"
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER
    design_config = _config(root)["design"]
    assert design_config["formula"] == formula
    assert design_config["variables"] == variables
    assert "condition_order" not in design_config
    assert ("Provisional model: " + formula) in result.output
    assert "Contrasts: set later by `rnaseq design`" in result.output
    assert ("Technical batch: yes" in result.output) is batch


def test_no_batch_flag_is_the_default(tmp_path):
    assert _new(tmp_path, "proj", "--no-batch").exit_code == 0
    assert (tmp_path / "proj" / "metadata.csv").read_text(encoding="utf-8") == "sample_id,condition\n"


def test_interactive_wizard_asks_about_batch_and_writes_the_template(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input=(
        "wizard\n" f"{tmp_path}\n" "human\n" "raw_counts\n" "\n"  # scope: default L1
        "paired_two_group\n" "patient\n" "y\n"  # design, pair column, known batch
        "\n" "\n" "y\n"
    ))
    assert result.exit_code == 0, result.output
    assert "Is there a known technical batch" in result.output
    assert "Contrasts: set later by `rnaseq design`" in result.output
    root = tmp_path / "wizard"
    assert (root / "metadata.csv").read_text(encoding="utf-8") == "sample_id,patient,condition,batch\n"
    assert _config(root)["project"]["preset"] == "L1"
    assert _config(root)["design"]["formula"] == "~ patient + batch + condition"
    # The wizard never asks for levels, a reference level, covariates, or contrasts.
    for forbidden in ("level", "Covariate", "Condition/contrast field", "Reference level"):
        assert forbidden not in result.output.split("Is there a known technical batch")[1].split("Local execution")[0]


def test_qc_scope_asks_no_design_questions(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input="qcproj\n" f"{tmp_path}\n" "human\n" "fastq\n" "\n\n\n\n\n" "qc\n" "\n\ny\n")
    assert result.exit_code == 0, result.output
    assert "technical batch" not in result.output


# ---------------------------------------------------------------- no silent flags


@pytest.mark.parametrize("flags", (
    ["--covariate", "batch"], ["--continuous-covariate", "age"], ["--condition-column", "group"],
))
def test_design_flags_without_metadata_are_errors_not_ignored(tmp_path, flags):
    result = _new(tmp_path, "proj", *flags)
    assert result.exit_code == 1
    assert "--metadata" in result.output
    assert not (tmp_path / "proj").exists()


def test_covariate_error_points_to_batch_flag(tmp_path):
    assert "--batch" in _new(tmp_path, "proj", "--covariate", "batch").output


def test_pair_id_requires_a_paired_design(tmp_path):
    result = _new(tmp_path, "proj", "--pair-id", "patient")
    assert result.exit_code == 1 and "paired_two_group" in result.output


@pytest.mark.parametrize("flag", (["--batch"], ["--covariate", "batch"], ["--condition-column", "condition"], ["--pair-id", "p"]))
def test_qc_scope_rejects_inapplicable_design_flags(tmp_path, flag):
    result = _new(tmp_path, "proj", *flag, preset="qc", input_type="fastq")
    assert result.exit_code == 1 and "qc preset" in result.output


def _import(tmp_path: Path, name: str, metadata: str, *extra: str, preset: str = "L1", design: str = "two_group"):
    counts, meta = tmp_path / f"{name}_counts.csv", tmp_path / f"{name}_meta.csv"
    meta.write_text(metadata, encoding="utf-8")
    samples = [line.split(",")[0] for line in metadata.strip().splitlines()[1:]]
    counts.write_text("gene_id," + ",".join(samples) + "\n" + "G1," + ",".join("5" for _ in samples) + "\nG2," + ",".join(str(3 + i) for i in range(len(samples))) + "\n", encoding="utf-8")
    contrasts = tmp_path / f"{name}_ct.csv"
    contrasts.write_text(CONTRAST_HEADER, encoding="utf-8")
    return runner.invoke(app, [
        "new", "--name", name, "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts",
        "--counts", str(counts), "--metadata", str(meta), "--contrasts", str(contrasts), "--preset", preset,
        "--design-type", design, "--yes", *extra,
    ])


def test_imported_metadata_without_condition_column_lists_available_fields(tmp_path):
    result = _import(tmp_path, "imp", "sample_id,group,batch\ns1,a,x\ns2,b,y\ns3,a,y\ns4,b,x\n")
    assert result.exit_code == 1
    assert "no 'condition' column" in result.output and "group, batch" in result.output


def test_imported_metadata_uses_condition_by_default_and_batch_flag(tmp_path):
    result = _import(tmp_path, "imp", TWO_GROUP_BATCH, "--batch")
    assert result.exit_code == 0, result.output
    config = _config(tmp_path / "imp")
    assert config["design"]["formula"] == "~ batch + condition"
    assert (tmp_path / "imp" / "metadata.csv").read_text(encoding="utf-8") == TWO_GROUP_BATCH


def test_imported_covariate_that_is_also_the_condition_is_an_error(tmp_path):
    result = _import(tmp_path, "imp", TWO_GROUP_BATCH, "--covariate", "condition")
    assert result.exit_code == 1 and "cannot also be a covariate" in result.output


def test_no_batch_conflicts_with_covariate_batch(tmp_path):
    result = _import(tmp_path, "imp", TWO_GROUP_BATCH, "--no-batch", "--covariate", "batch")
    assert result.exit_code == 1 and "--no-batch" in result.output


def test_batch_cannot_be_combined_with_an_explicit_formula(tmp_path):
    from rnaseq.errors import ProjectCreationError
    from rnaseq.models import DesignType, Preset, Species
    from rnaseq.project import create_project

    with pytest.raises(ProjectCreationError):
        create_project(project_name="p", destination=tmp_path, species=Species.MOUSE, preset=Preset.L1,
                       design_type=DesignType.TWO_GROUP, batch=True, formula="~ condition")


# ---------------------------------------------------------------- validation without contrasts


def test_l1_validates_with_empty_contrasts_but_l2_does_not(tmp_path):
    root_l1 = tmp_path / "l1"
    _new(tmp_path, "l1", "--batch")
    _fill(root_l1, TWO_GROUP_BATCH)
    report = validate_project(root_l1)
    assert report.is_valid, [i.message for i in report.errors]
    assert report.contrasts is not None and not report.contrasts.contrasts

    _new(tmp_path, "l2", "--batch", preset="L2")
    _fill(tmp_path / "l2", TWO_GROUP_BATCH)
    report = validate_project(tmp_path / "l2")
    assert [i.code for i in report.errors] == ["missing_contrast_rows"]
    assert "rnaseq design" in report.errors[0].message


def test_header_only_template_validates_only_as_missing_rows(tmp_path):
    _new(tmp_path, "proj")
    assert {i.code for i in validate_project(tmp_path / "proj").errors} >= {"missing_metadata_rows"}


@pytest.mark.parametrize(("design", "metadata", "code"), (
    ("two_group", THREE_GROUP, "invalid_two_group_design"),
    ("multi_group", "sample_id,condition\na,X\nb,X\nc,Y\nd,Y\n", "invalid_multi_group_design"),
))
def test_design_type_level_count_is_checked_without_contrasts(tmp_path, design, metadata, code):
    _new(tmp_path, "proj", design=design)
    _fill(tmp_path / "proj", metadata)
    assert code in {i.code for i in validate_project(tmp_path / "proj").errors}


def test_paired_design_blank_pair_and_incomplete_pairs_are_checked_without_contrasts(tmp_path):
    _new(tmp_path, "proj", "--pair-id", "patient", design="paired_two_group")
    _fill(tmp_path / "proj", "sample_id,patient,condition\na1,P1,Control\na2,P1,Treatment\nb1,P2,Control\nb2,,Treatment\n")
    codes = {i.code for i in validate_project(tmp_path / "proj").errors}
    assert "blank_pair_id" in codes and "pair_missing_numerator" in codes


def test_valid_paired_design_without_contrasts_has_no_recorded_pairing_contrast(tmp_path):
    _new(tmp_path, "proj", "--pair-id", "patient", design="paired_two_group")
    _fill(tmp_path / "proj", "sample_id,patient,condition\na1,P1,Control\na2,P1,Treatment\nb1,P2,Control\nb2,P2,Treatment\nc1,P3,Control\nc2,P3,Treatment\n")
    report = validate_project(tmp_path / "proj")
    assert report.is_valid, [i.message for i in report.errors]
    assert report.pairing is not None and report.pairing.contrasts == ()


# ---------------------------------------------------------------- estimability


def _validate_metadata(tmp_path: Path, metadata: str, *flags: str, preset: str = "L1", design: str = "two_group"):
    name = f"p{len(list(tmp_path.iterdir()))}"
    assert _new(tmp_path, name, *flags, preset=preset, design=design).exit_code == 0
    _fill(tmp_path / name, metadata)
    return validate_project(tmp_path / name)


def test_perfectly_confounded_batch_and_condition_is_an_error_before_l1(tmp_path):
    report = _validate_metadata(tmp_path, "sample_id,condition,batch\na,C,A\nb,C,A\nc,T,B\nd,T,B\n", "--batch")
    assert "rank_deficient_design" in {i.code for i in report.errors}
    assert any("'batch' and 'condition' are confounded" in i.message for i in report.errors)


def test_batch_that_is_a_function_of_condition_in_a_multigroup_is_explained(tmp_path):
    report = _validate_metadata(tmp_path, "sample_id,condition,batch\na,X,B1\nb,X,B1\nc,Y,B1\nd,Y,B1\ne,Z,B2\nf,Z,B2\n", "--batch", design="multi_group")
    errors = [i.message for i in report.errors if i.code == "rank_deficient_design"]
    assert errors and "cannot be estimated separately" in errors[0]


def test_paired_batch_constant_within_each_pair_is_absorbed_by_the_pairing(tmp_path):
    report = _validate_metadata(
        tmp_path, "sample_id,patient,condition,batch\na1,P1,C,A\na2,P1,T,A\nb1,P2,C,B\nb2,P2,T,B\nc1,P3,C,A\nc2,P3,T,A\n",
        "--pair-id", "patient", "--batch", design="paired_two_group",
    )
    messages = [i.message for i in report.errors if i.code == "rank_deficient_design"]
    assert messages and "pairing already absorbs it" in messages[0]


def test_single_level_batch_is_an_error(tmp_path):
    report = _validate_metadata(tmp_path, "sample_id,condition,batch\na,C,A\nb,C,A\nc,T,A\nd,T,A\n", "--batch")
    assert "single_level_design_variable" in {i.code for i in report.errors}


def test_blank_batch_value_is_an_error(tmp_path):
    report = _validate_metadata(tmp_path, "sample_id,condition,batch\na,C,A\nb,C,\nc,T,B\nd,T,A\n", "--batch")
    assert "missing_design_value" in {i.code for i in report.errors}


def test_residual_degrees_of_freedom_use_the_expanded_model_matrix(tmp_path):
    # 4 samples, ~ batch + condition = intercept + batch[B] + condition[T] = 3 columns -> 1 residual df.
    ok = _validate_metadata(tmp_path, "sample_id,condition,batch\na,C,A\nb,C,B\nc,T,A\nd,T,B\n", "--batch", preset="L2")
    assert "no_residual_degrees_of_freedom" not in {i.code for i in (*ok.errors, *ok.warnings)}
    # A three-level batch adds a column: 4 samples, 4 columns -> no residual df, though the formula has only two terms.
    metadata = "sample_id,condition,batch\na,C,A\nb,C,B\nc,T,C\nd,T,A\n"
    for preset in ("L1", "L2"):
        report = _validate_metadata(tmp_path, metadata, "--batch", preset=preset)
        # Like limited replication this is a warning, so L1 QC stays available; L2 is refused at design/launch time.
        assert {i.code for i in report.errors} <= {"missing_contrast_rows"}  # the only L2 gap is the not-yet-designed contrasts
        assert "no_residual_degrees_of_freedom" in {i.code for i in report.warnings}


def test_l2_launch_guards_refuse_a_model_without_residual_degrees_of_freedom(tmp_path):
    from rnaseq.errors import DownstreamExecutionError
    from rnaseq.l2 import prepare_l2

    _new(tmp_path, "proj", "--batch", preset="L2")
    root = tmp_path / "proj"
    _fill(root, "sample_id,condition,batch\na,C,A\nb,C,B\nc,T,C\nd,T,A\n")
    (root / "contrasts.csv").write_text(CONTRAST_HEADER + "T_vs_C,condition,T,C\n", encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    with pytest.raises(DownstreamExecutionError, match="no residual degrees of freedom: 4 sample"):
        prepare_l2(report, run_id=None)


def test_frozen_l2_config_builder_refuses_a_model_without_residual_degrees_of_freedom(project_factory, tmp_path):
    from test_workflow_support import _typed_batch_l2_contract

    from rnaseq.workflow_support import l2_config

    contract, inputs, _validation = _typed_batch_l2_contract(
        project_factory, design_type="two_group", formula="~ batch + condition",
        variables={"batch": "categorical", "condition": "categorical"},
        metadata="sample_id,condition,batch\na,C,A\nb,C,B\nc,T,C\nd,T,A\n",
        contrasts="contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n",
    )
    with pytest.raises(ValueError, match="no residual degrees of freedom: 4 sample.*L1 QC remains available"):
        l2_config(contract, inputs, tmp_path / "l1", tmp_path / "l2")


def test_partial_confounding_and_single_sample_batch_are_warned(tmp_path):
    report = _validate_metadata(
        tmp_path, "sample_id,condition,batch\na,C,A\nb,C,B\nc,C,B\nd,T,B\ne,T,B\nf,T,B\n", "--batch",
    )
    codes = {i.code for i in report.warnings}
    assert "partially_confounded_batch" in codes and "batch_level_single_sample" in codes
    assert report.is_valid


def test_unused_and_misnamed_batch_columns_are_surfaced(tmp_path):
    report = _validate_metadata(tmp_path, "sample_id,condition,Batch\na,C,A\nb,C,B\nc,T,A\nd,T,B\n")
    assert "nonstandard_column_name" in {i.code for i in report.warnings}
    report = _validate_metadata(tmp_path, "sample_id,condition,batch\na,C,A\nb,C,B\nc,T,A\nd,T,B\n")
    assert "unused_batch_column" in {i.code for i in report.warnings}


def test_validation_report_shows_the_batch_condition_contingency(tmp_path):
    _new(tmp_path, "proj", "--batch")
    _fill(tmp_path / "proj", TWO_GROUP_BATCH)
    output = runner.invoke(app, ["validate", str(tmp_path / "proj")]).output
    assert "Batch x condition (samples)" in output
    assert "Contrasts: not yet defined" in output


# ---------------------------------------------------------------- contrast generation


def test_two_level_contrast_is_later_over_earlier():
    assert generate_contrasts("condition", ("Control", "Treatment")) == (
        ContrastDefinition("Treatment_vs_Control", "condition", "Treatment", "Control"),
    )


def test_three_levels_generate_the_agreed_pairwise_contrasts_in_order():
    contrasts = generate_contrasts("condition", ("Control", "DrugA", "DrugB"))
    assert [c.contrast_id for c in contrasts] == ["DrugA_vs_Control", "DrugB_vs_Control", "DrugB_vs_DrugA"]
    assert [(c.numerator, c.denominator) for c in contrasts] == [("DrugA", "Control"), ("DrugB", "Control"), ("DrugB", "DrugA")]


def test_four_levels_generate_six_contrasts():
    contrasts = generate_contrasts("condition", ("A", "B", "C", "D"))
    assert [c.contrast_id for c in contrasts] == ["B_vs_A", "C_vs_A", "D_vs_A", "C_vs_B", "D_vs_B", "D_vs_C"]
    assert len(contrasts) == 4 * 3 // 2


def test_unsafe_level_names_are_reported():
    with pytest.raises(DesignSetupError, match="'Drug A'"):
        generate_contrasts("condition", ("Control", "Drug A"))


# ---------------------------------------------------------------- rnaseq design


def _ready(tmp_path: Path, name: str = "proj", metadata: str = TWO_GROUP_BATCH, *flags: str, preset: str = "L1", design: str = "two_group") -> Path:
    assert _new(tmp_path, name, *flags, preset=preset, design=design).exit_code == 0
    root = tmp_path / name
    _fill(root, metadata)
    return root


def _design(root: Path, *args: str):
    return runner.invoke(app, ["design", str(root), *args])


def test_design_two_group_writes_confirmed_contrast_and_freezes_order(tmp_path):
    root = _ready(tmp_path, "proj", TWO_GROUP_BATCH, "--batch")
    result = _design(root, "--yes")
    assert result.exit_code == 0, result.output
    assert "Treatment_vs_Control" in result.output and "log2FC > 0 = higher in Treatment" in result.output
    assert "Batch x condition" in result.output and "Residual degrees of freedom: 6 samples - 3 model coefficients = 3" in result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER + "Treatment_vs_Control,condition,Treatment,Control\n"
    config = _config(root)
    assert config["project"]["preset"] == "L2"
    assert config["design"]["condition_order"] == ["Control", "Treatment"]
    assert config["design"]["formula"] == "~ batch + condition"
    assert validate_project(root).is_valid


def test_design_first_appearance_decides_direction_and_swap_reverses_it(tmp_path):
    treatment_first = TWO_GROUP_BATCH.replace("Control", "X").replace("Treatment", "Control").replace("X", "Treatment")
    root = _ready(tmp_path, "proj", treatment_first)
    assert _design(root, "--yes").exit_code == 0
    assert (root / "contrasts.csv").read_text(encoding="utf-8").splitlines()[1] == "Control_vs_Treatment,condition,Control,Treatment"
    root2 = _ready(tmp_path, "proj2", TWO_GROUP_BATCH)
    assert _design(root2, "--yes", "--swap").exit_code == 0
    assert (root2 / "contrasts.csv").read_text(encoding="utf-8").splitlines()[1] == "Control_vs_Treatment,condition,Control,Treatment"
    # Swap reverses only the contrast; the condition order stays the metadata's first-appearance order.
    assert _config(root2)["design"]["condition_order"] == ["Control", "Treatment"]
    assert _config(root)["design"]["condition_order"] == ["Treatment", "Control"]  # that metadata listed Treatment first


def test_design_interactive_swap_then_accept(monkeypatch, tmp_path):
    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    root = _ready(tmp_path, "proj", TWO_GROUP_BATCH)
    result = runner.invoke(app, ["design", str(root)], input="swap\ny\n")
    assert result.exit_code == 0, result.output
    assert result.output.count("Proposed contrast:") == 2
    assert "Control_vs_Treatment" in (root / "contrasts.csv").read_text(encoding="utf-8")
    assert _config(root)["design"]["condition_order"] == ["Control", "Treatment"]
    # Swapping twice returns to the default direction, still without touching the order.
    root2 = _ready(tmp_path, "proj2", TWO_GROUP_BATCH)
    assert runner.invoke(app, ["design", str(root2)], input="swap\nswap\ny\n").exit_code == 0
    assert (root2 / "contrasts.csv").read_text(encoding="utf-8").splitlines()[1] == "Treatment_vs_Control,condition,Treatment,Control"
    assert _config(root2)["design"]["condition_order"] == ["Control", "Treatment"]


def test_swap_is_only_for_two_level_designs(tmp_path):
    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    result = _design(root, "--yes", "--swap")
    assert result.exit_code == 1 and "two-level" in result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER


def test_design_declined_writes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    root = _ready(tmp_path, "proj")
    before = {name: (root / name).read_bytes() for name in ("project.yaml", "contrasts.csv", "metadata.csv")}
    result = runner.invoke(app, ["design", str(root)], input="n\n")
    assert result.exit_code == 0 and "No changes were written" in result.output
    assert before == {name: (root / name).read_bytes() for name in before}


def test_design_needs_a_terminal_or_yes(tmp_path):
    root = _ready(tmp_path, "proj")
    result = _design(root)
    assert result.exit_code == 1 and "--yes" in result.output


def test_design_multi_group_writes_all_pairwise_contrasts_after_confirmation(monkeypatch, tmp_path):
    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    declined = runner.invoke(app, ["design", str(root)], input="n\n")
    assert declined.exit_code == 0 and (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER
    result = runner.invoke(app, ["design", str(root)], input="y\n")
    assert result.exit_code == 0, result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER + (
        "DrugA_vs_Control,condition,DrugA,Control\nDrugB_vs_Control,condition,DrugB,Control\nDrugB_vs_DrugA,condition,DrugB,DrugA\n"
    )
    assert _config(root)["design"]["condition_order"] == ["Control", "DrugA", "DrugB"]
    assert validate_project(root).is_valid


def test_design_condition_order_override_and_frozen_order_survive_row_resorting(tmp_path):
    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    assert _design(root, "--yes", "--condition-order", "DrugB,Control,DrugA").exit_code == 0
    assert _config(root)["design"]["condition_order"] == ["DrugB", "Control", "DrugA"]
    ids = [line.split(",")[0] for line in (root / "contrasts.csv").read_text(encoding="utf-8").splitlines()[1:]]
    assert ids == ["Control_vs_DrugB", "DrugA_vs_DrugB", "DrugA_vs_Control"]
    # Re-sorting the metadata rows later cannot change the confirmed direction.
    lines = (root / "metadata.csv").read_text(encoding="utf-8").splitlines()
    (root / "metadata.csv").write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n", encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and list(report.groups["condition"]) == ["DrugB", "Control", "DrugA"]
    # A forced re-run reuses the frozen order unless told otherwise.
    assert _design(root, "--yes", "--force").exit_code == 0
    assert _config(root)["design"]["condition_order"] == ["DrugB", "Control", "DrugA"]


def test_design_rejects_a_condition_order_that_does_not_match_metadata(tmp_path):
    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    result = _design(root, "--yes", "--condition-order", "Control,DrugA")
    assert result.exit_code == 1 and "must list exactly the levels" in result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER


def test_stale_frozen_condition_order_is_a_validation_error(tmp_path):
    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    assert _design(root, "--yes").exit_code == 0
    (root / "metadata.csv").write_text(THREE_GROUP.replace("DrugB", "DrugC"), encoding="utf-8")
    assert "condition_order_mismatch" in {i.code for i in validate_project(root).errors}


def test_design_never_overwrites_existing_contrasts_without_force(tmp_path):
    root = _ready(tmp_path, "proj")
    custom = CONTRAST_HEADER + "mine,condition,Treatment,Control\n"
    (root / "contrasts.csv").write_text(custom, encoding="utf-8")
    before = (root / "project.yaml").read_bytes()
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "already defines 1 contrast" in result.output and "--force" in result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == custom
    assert (root / "project.yaml").read_bytes() == before
    forced = _design(root, "--yes", "--force")
    assert forced.exit_code == 0 and "These existing contrasts will be replaced" in forced.output
    assert "Treatment_vs_Control" in (root / "contrasts.csv").read_text(encoding="utf-8")


def test_design_accepts_an_l2_project_with_empty_contrasts(tmp_path):
    root = _ready(tmp_path, "proj", TWO_GROUP_BATCH, "--batch", preset="L2")
    assert not validate_project(root).is_valid
    assert _design(root, "--yes").exit_code == 0
    assert validate_project(root).is_valid


def test_design_refuses_qc_only_projects(tmp_path):
    assert _new(tmp_path, "qc", preset="qc", input_type="fastq").exit_code == 0
    result = _design(tmp_path / "qc", "--yes")
    assert result.exit_code == 1 and "QC-only" in result.output


def test_design_refuses_invalid_metadata_and_writes_nothing(tmp_path):
    root = _ready(tmp_path, "proj", "sample_id,condition,batch\na,C,A\nb,C,A\nc,T,B\nd,T,B\n", "--batch")
    before = {name: (root / name).read_bytes() for name in ("project.yaml", "contrasts.csv")}
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "not valid yet" in result.output and "confounded" in result.output
    assert before == {name: (root / name).read_bytes() for name in before}


def test_design_refuses_a_model_without_residual_degrees_of_freedom_even_from_l1(tmp_path):
    root = _ready(tmp_path, "proj", "sample_id,condition,batch\na,C,A\nb,C,B\nc,T,C\nd,T,A\n", "--batch")
    assert validate_project(root).is_valid  # L1 only warns
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "no residual degrees of freedom" in result.output
    assert _config(root)["project"]["preset"] == "L1"


def test_design_reports_unsafe_level_names_without_writing(tmp_path):
    root = _ready(tmp_path, "proj", "sample_id,condition\na,Drug A\nb,Drug A\nc,Control\nd,Control\n")
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "'Drug A'" in result.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER


def test_design_keeps_a_custom_formula_and_comments(tmp_path):
    root = _ready(tmp_path, "proj", "sample_id,condition,sex\na,C,M\nb,C,F\nc,T,M\nd,T,F\ne,C,M\nf,T,F\n")
    config = _config(root)
    config["design"]["formula"] = "~ sex + condition"
    config["design"]["variables"] = {"sex": "categorical", "condition": "categorical"}
    text = "# my notes\n" + yaml.safe_dump(config, sort_keys=False)
    (root / "project.yaml").write_text(text, encoding="utf-8")
    result = _design(root, "--yes")
    assert result.exit_code == 0, result.output
    assert "custom formula kept as written" in result.output
    updated = (root / "project.yaml").read_text(encoding="utf-8")
    assert updated.startswith("# my notes\n")
    assert _config(root)["design"]["formula"] == "~ sex + condition"


def test_design_works_for_a_paired_two_group_project(tmp_path):
    root = _ready(
        tmp_path, "proj", "sample_id,patient,condition,batch\na1,P1,Control,A\na2,P1,Treatment,A\nb1,P2,Control,A\nb2,P2,Treatment,A\n"
        "c1,P3,Control,B\nc2,P3,Treatment,B\nd1,P4,Control,B\nd2,P4,Treatment,B\n", "--pair-id", "patient", "--batch", design="paired_two_group",
    )
    result = _design(root, "--yes")
    # Batch is constant within each pair, so it is absorbed by the pairing and the design is refused.
    assert result.exit_code == 1 and "pairing already absorbs it" in result.output
    root = _ready(
        tmp_path, "proj2", "sample_id,patient,condition\na1,P1,Control\na2,P1,Treatment\nb1,P2,Control\nb2,P2,Treatment\nc1,P3,Control\nc2,P3,Treatment\n",
        "--pair-id", "patient", design="paired_two_group",
    )
    assert _design(root, "--yes").exit_code == 0
    assert _config(root)["design"]["formula"] == "~ patient + condition"
    assert validate_project(root).pairing.contrasts[0].complete_pairs == 3


def test_design_rolls_back_if_the_final_l2_validation_fails(monkeypatch, tmp_path):
    root = _ready(tmp_path, "proj")
    before = {name: (root / name).read_bytes() for name in ("project.yaml", "contrasts.csv")}
    real = validate_project
    calls = {"count": 0}

    def flaky(path):
        calls["count"] += 1
        report = real(path)
        if calls["count"] == 2:
            report.error("forced", "forced failure")
        return report

    monkeypatch.setattr("rnaseq.cli.validate_project", flaky)
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "nothing was changed" in result.output
    assert before == {name: (root / name).read_bytes() for name in before}


def test_design_refreshes_the_plan_so_run_preflight_stays_fresh(tmp_path):
    """Regression: `run` refuses stale planning artifacts, and `design` changes the planned configuration."""

    from rnaseq.errors import ExecutionPreflightError
    from rnaseq.execution import require_fresh_plan

    root = _ready(tmp_path, "proj", THREE_GROUP, design="multi_group")
    assert runner.invoke(app, ["plan", str(root)]).exit_code == 0
    require_fresh_plan(validate_project(root))  # fresh before design
    result = _design(root, "--yes")
    assert result.exit_code == 0 and "Planning artifacts refreshed" in result.output
    require_fresh_plan(validate_project(root))  # still fresh: design refreshed the plan
    manifest = yaml.safe_load((root / "planning" / "manifest.preview.yaml").read_text(encoding="utf-8"))
    assert [item["contrast_id"] for item in manifest["contrasts"]["definitions"]] == ["DrugA_vs_Control", "DrugB_vs_Control", "DrugB_vs_DrugA"]
    assert manifest["planned_downstream"]["preset"] == "L2"
    # The guard itself is unchanged: a later manual edit still makes the plan stale.
    (root / "contrasts.csv").write_text(CONTRAST_HEADER + "x,condition,DrugA,Control\n", encoding="utf-8")
    with pytest.raises(ExecutionPreflightError, match="changed since planning"):
        require_fresh_plan(validate_project(root))


def test_design_rolls_back_if_the_plan_cannot_be_refreshed(monkeypatch, tmp_path):
    root = _ready(tmp_path, "proj")
    before = {name: (root / name).read_bytes() for name in ("project.yaml", "contrasts.csv")}

    def broken(_report):
        raise OSError("disk full")

    monkeypatch.setattr("rnaseq.cli.generate_plan", broken)
    result = _design(root, "--yes")
    assert result.exit_code == 1 and "nothing was changed" in result.output
    assert before == {name: (root / name).read_bytes() for name in before}


def test_design_never_mutates_existing_run_records(tmp_path):
    root = _ready(tmp_path, "proj")
    frozen = root / "runs" / "L1-001" / "20260101-000000+0000" / "frozen"
    frozen.mkdir(parents=True)
    for name in ("project.yaml", "metadata.csv", "contrasts.csv"):
        (frozen / name).write_bytes((root / name).read_bytes())
    (frozen.parent / "run_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")

    def snapshot():
        return {str(path.relative_to(root / "runs")): path.read_bytes() for path in sorted((root / "runs").rglob("*")) if path.is_file()}

    before = snapshot()
    assert _design(root, "--yes").exit_code == 0
    assert snapshot() == before  # the L1 run keeps its frozen preset: L1 and its empty contrasts
    assert yaml.safe_load((frozen / "project.yaml").read_text(encoding="utf-8"))["project"]["preset"] == "L1"
    assert _config(root)["project"]["preset"] == "L2"


def test_design_points_to_the_latest_successful_run_for_reuse(tmp_path):
    root = _ready(tmp_path, "proj")
    for case, run in (("L1-001", "20260101-000000+0000"), ("L1-002", "20260102-000000+0000")):
        directory = root / "runs" / case / run
        directory.mkdir(parents=True)
        (directory / "run_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    failed = root / "runs" / "L1-003" / "20260103-000000+0000"
    failed.mkdir(parents=True)
    (failed / "run_state.json").write_text(json.dumps({"status": "FAILED"}), encoding="utf-8")
    assert "--reuse-upstream L1-002/20260102-000000+0000" in _design(root, "--yes").output


# ---------------------------------------------------------------- project.yaml editing


def test_updated_project_yaml_preserves_other_content_and_replaces_an_existing_order():
    config = base_config()
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "multi_group", "formula": "~ condition", "variables": {"condition": "categorical"}, "condition_order": ["B", "A", "C"]}
    text = "# keep\n" + yaml.safe_dump(config, sort_keys=False)
    updated = updated_project_yaml(text, ("A", "B", "C"))
    parsed = yaml.safe_load(updated)
    assert parsed["design"]["condition_order"] == ["A", "B", "C"]
    assert parsed["project"]["preset"] == "L2"
    assert updated.startswith("# keep\n")
    parsed["design"]["condition_order"] = config["design"]["condition_order"]
    parsed["project"]["preset"] = config["project"]["preset"]
    assert parsed == config  # nothing else changed
    ProjectConfig.model_validate(yaml.safe_load(updated))


# ---------------------------------------------------------------- end to end


def test_two_group_known_batch_lifecycle_ends_in_the_batch_adjusted_deseq2_model(tmp_path):
    """new --batch -> fill metadata -> L1 validation (empty contrasts) -> design --yes -> L2 validation -> final model."""

    from rnaseq.downstream import execute_l1, prepare_l1
    from rnaseq.l2 import execute_l2, prepare_l2

    require_r_packages("jsonlite", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    created = _new(tmp_path, "e2e", "--batch")
    assert created.exit_code == 0, created.output
    root = tmp_path / "e2e"
    # 1. The generated template.
    assert (root / "metadata.csv").read_text(encoding="utf-8") == "sample_id,condition,batch\n"
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER
    # 2. Completed metadata, with enough genes for real DESeq2.
    metadata = TWO_GROUP_BATCH
    samples = [line.split(",")[0] for line in metadata.strip().splitlines()[1:]]
    (root / "metadata.csv").write_text(metadata, encoding="utf-8")
    import random

    rng = random.Random(7)  # deterministic over-dispersed counts: gene-wise dispersions must span a range for DESeq2's trend fit
    lines = []
    for gene in range(1, 161):
        mean, spread = 30 + gene * 4, 0.05 + (gene % 40) / 40
        lines.append(f"Gene{gene}," + ",".join(
            str(max(int(mean * (1.8 if sample >= "s4" and gene % 5 == 0 else 1.0) * max(rng.gauss(1.0, spread), 0.05)), 1)) for sample in samples
        ))
    counts = "gene_id," + ",".join(samples) + "\n" + "\n".join(lines) + "\n"
    (root / "input" / "counts.csv").write_text(counts, encoding="utf-8")
    # 3. L1 validates and runs with empty contrasts; the provisional design is batch-adjusted but L1 is uncorrected QC.
    report = validate_project(root)
    assert report.is_valid, [issue.message for issue in report.errors]
    assert report.config.project.preset.value == "L1" and not report.contrasts.contrasts
    l1 = execute_l1(prepare_l1(report, run_id=None))
    assert json.loads((l1.output_dir / "backend_config.json").read_text(encoding="utf-8"))["formula"] == "~ batch + condition"
    assert not list(l1.output_dir.glob("*removeBatch*")) and not list(l1.output_dir.glob("*batch_corrected*"))
    # 4. design --yes proposes and writes the contrast and promotes the project to L2.
    designed = _design(root, "--yes")
    assert designed.exit_code == 0, designed.output
    assert (root / "contrasts.csv").read_text(encoding="utf-8") == CONTRAST_HEADER + "Treatment_vs_Control,condition,Treatment,Control\n"
    # 5. L2 validates and the final fitted model is ~ batch + condition.
    l2_report = validate_project(root)
    assert l2_report.is_valid and l2_report.config.project.preset.value == "L2"
    result = execute_l2(prepare_l2(l2_report, run_id=None))
    design = json.loads((result.output_dir / "scientific_provenance.json").read_text(encoding="utf-8"))["design"]
    assert design["formula"] == "~ batch + condition"
    assert design["model_matrix"]["column_names"] == ["(Intercept)", "batchB", "conditionTreatment"]
    assert design["model_matrix"]["full_rank"] is True
    assert (result.output_dir / "contrasts" / "Treatment_vs_Control" / "all_genes.tsv").is_file()
