"""Phase 1 figure-style configuration: validation, resolution and propagation.

R still renders the legacy look in this phase; these tests only pin what is
configured and what reaches each backend / the frozen contract.
"""

from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from pydantic import ValidationError
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app
from rnaseq.downstream import execute_l1, prepare_l1
from rnaseq.errors import DownstreamExecutionError, ExecutionPreflightError
from rnaseq.figure_settings import (
    figure_settings_from_contract,
    figure_settings_from_project,
    resolve_figure_settings,
    resolve_run_figure_settings,
)
from rnaseq.go import execute_go, prepare_go
from rnaseq.gsea import execute_gsea, prepare_gsea
from rnaseq.kegg import execute_kegg, prepare_kegg
from rnaseq.l2 import execute_l2, prepare_l2
from rnaseq.models import FiguresConfig, ProjectConfig
from rnaseq.planner import generate_plan
from rnaseq.service import create_case_run, execute_retry_service_run, freeze_case_inputs, resolve_downstream_inputs
from rnaseq.validators import validate_project
from rnaseq.workflow_support import enrichment_config, l1_config, l2_config

from test_retry import _failed_frozen_run, _mock_runtime, _successful_downstream
from test_service import _REUSE_SOURCE, _UPSTREAM_MANIFEST, _frozen_run, _reuse_pair

runner = CliRunner()
LEGACY = {"preset": "legacy", "profile": "general", "source": "legacy_default"}


def _explicit(preset: str, profile: str = "general") -> dict[str, str]:
    return {"preset": preset, "profile": profile, "source": "explicit"}


# ----------------------------------------------------------------- model validation


@pytest.mark.parametrize("preset", ["standard", "legacy"])
@pytest.mark.parametrize("profile", ["general", "nature", "plos", "elsevier"])
def test_valid_preset_and_profile_values(preset, profile):
    assert FiguresConfig(preset=preset, profile=profile).model_dump() == {"preset": preset, "profile": profile}


def test_profile_defaults_to_general_but_preset_is_required():
    assert FiguresConfig(preset="standard").profile == "general"
    with pytest.raises(ValidationError):
        FiguresConfig.model_validate({"profile": "nature"})
    with pytest.raises(ValidationError):
        FiguresConfig.model_validate({})


@pytest.mark.parametrize(
    "value",
    [
        {"preset": "Bioinfo"},
        {"preset": "modern"},
        {"preset": ""},
        {"preset": None},
        {"preset": 1},
        {"preset": ["standard"]},
        {"preset": "standard", "profile": "science"},
        {"preset": "standard", "profile": "NATURE"},
        {"preset": "standard", "profile": None},
        {"preset": "standard", "font": "Arial"},
    ],
)
def test_invalid_figures_values_are_rejected(value):
    with pytest.raises(ValidationError):
        FiguresConfig.model_validate(value)


@pytest.mark.parametrize("value", ["standard", ["standard"], 3, {"preset": "nope"}])
def test_project_rejects_invalid_figures_section(value):
    config = base_config()
    config["figures"] = value
    with pytest.raises(ValidationError):
        ProjectConfig.model_validate(config)


def test_validate_project_reports_an_invalid_figures_preset(project_factory):
    config = base_config()
    config["figures"] = {"preset": "modern"}
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert any("figures" in issue.message for issue in report.errors)


# ------------------------------------------------------- old inputs resolve to legacy


def test_existing_project_yaml_without_figures_is_valid_and_resolves_to_legacy(project_factory):
    report = validate_project(project_factory())
    assert report.is_valid, report.errors
    assert report.config.figures is None
    assert resolve_figure_settings(report.config.figures) == LEGACY


@pytest.mark.parametrize("schema_version", ["1.0", "1.1", "1.2"])
def test_every_supported_schema_version_without_figures_resolves_to_legacy(schema_version):
    config = base_config()
    config["schema_version"] = schema_version
    if schema_version in {"1.2"}:
        config["analysis"] = {"enrichment": []}
    parsed = ProjectConfig.model_validate(config)
    assert resolve_figure_settings(parsed.figures) == LEGACY


def test_null_figures_key_resolves_to_legacy():
    config = base_config()
    config["figures"] = None
    assert resolve_figure_settings(ProjectConfig.model_validate(config).figures) == LEGACY


@pytest.mark.parametrize(("preset", "profile"), [("legacy", "general"), ("legacy", "nature"), ("standard", "plos")])
def test_explicit_selection_is_always_respected(project_factory, preset, profile):
    config = base_config()
    config["figures"] = {"preset": preset, "profile": profile}
    report = validate_project(project_factory(config=config))
    assert report.is_valid, report.errors
    assert resolve_figure_settings(report.config.figures) == _explicit(preset, profile)
    assert figure_settings_from_project(config) == _explicit(preset, profile)


def test_explicit_legacy_is_distinguishable_from_the_implicit_default():
    assert resolve_figure_settings(FiguresConfig(preset="legacy")) == _explicit("legacy")
    assert resolve_figure_settings(None) == LEGACY


def test_resolver_accepts_a_raw_mapping_and_rejects_a_bad_one():
    assert resolve_figure_settings({"preset": "standard", "profile": "elsevier"}) == _explicit("standard", "elsevier")
    with pytest.raises(ValueError, match="invalid figures configuration"):
        resolve_figure_settings({"preset": "standard", "profile": "cell"})
    with pytest.raises(ValueError, match="invalid figures configuration"):
        resolve_figure_settings("standard")  # type: ignore[arg-type]


def test_contract_without_figures_resolves_to_legacy_and_explicit_is_kept():
    assert figure_settings_from_contract({"analysis_level": "L2"}) == LEGACY
    assert figure_settings_from_contract({"figures": None}) == LEGACY
    assert figure_settings_from_contract({"figures": _explicit("standard", "nature")}) == _explicit("standard", "nature")


@pytest.mark.parametrize(
    "bad",
    [{"preset": "standard"}, {"preset": "x", "profile": "general", "source": "explicit"},
     {"preset": "standard", "profile": "general", "source": "guess"}, "standard", []],
)
def test_malformed_contract_figures_fail_loudly(bad):
    with pytest.raises(ValueError, match="frozen downstream contract"):
        figure_settings_from_contract({"figures": bad})


# ------------------------------------------------------------ new project template


@pytest.mark.parametrize(("input_type", "preset"), [("raw_counts", "L2"), ("fastq", "L1")])
def test_rnaseq_new_writes_an_explicit_standard_selection(tmp_path, input_type, preset):
    arguments = ["new", "--name", "created_project", "--destination", str(tmp_path), "--species", "mouse", "--input-type", input_type, "--preset", preset, "--design-type", "two_group", "--scaffold", "--yes"]
    if preset == "L2" and input_type == "raw_counts":
        arguments += ["--gene-id-type", "ENSEMBL"]  # L2 enables enrichment by default, which needs the count matrix gene ID type
    result = runner.invoke(app, arguments)
    assert result.exit_code == 0, result.output
    raw = yaml.safe_load((tmp_path / "created_project" / "project.yaml").read_text(encoding="utf-8"))
    assert raw["figures"] == {"preset": "standard", "profile": "general"}
    assert ProjectConfig.model_validate(raw).figures == FiguresConfig(preset="standard", profile="general")


def test_create_project_helper_writes_standard_for_every_scope():
    from rnaseq.models import DesignType, InputType, Preset, Species
    from rnaseq.project import _project_yaml

    for preset in Preset:
        raw = yaml.safe_load(_project_yaml("p", Species.MOUSE, preset, DesignType.TWO_GROUP, InputType.RAW_COUNTS, None))
        assert raw["figures"] == {"preset": "standard", "profile": "general"}


def test_shipped_project_yaml_template_selects_standard_and_validates():
    from importlib import resources

    raw = yaml.safe_load(resources.files("rnaseq.templates").joinpath("project.yaml").read_text(encoding="utf-8"))
    assert raw["figures"] == {"preset": "standard", "profile": "general"}
    assert ProjectConfig.model_validate(raw).figures.preset == "standard"


def test_design_command_preserves_an_explicit_figures_selection():
    from rnaseq.design_setup import updated_project_yaml

    text = yaml.safe_dump({**base_config(), "figures": {"preset": "legacy", "profile": "general"}}, sort_keys=False)
    updated = yaml.safe_load(updated_project_yaml(text, ("Control", "Treatment")))
    assert updated["figures"] == {"preset": "legacy", "profile": "general"}
    plain = yaml.safe_load(updated_project_yaml(yaml.safe_dump(base_config(), sort_keys=False), ("Control", "Treatment")))
    assert "figures" not in plain  # nothing is injected into a project that never chose


# ------------------------------------------------------------- direct propagation


def _capture_config(monkeypatch, module: str, execute, prepared) -> dict:
    """Run ``execute`` with a failing fake Rscript; return the backend_config.json it wrote.

    figures.R still accepts the older ``cfg$figures$profile`` hook, which must not be fed
    by the Python side.  The setting therefore travels under
    ``figure_style`` (R reads it through nf_rna_read_figure_style), and ``figures`` must never reach R.
    """

    monkeypatch.setattr(f"{module}.subprocess.run", lambda *_a, **_k: SimpleNamespace(returncode=1))
    with pytest.raises(DownstreamExecutionError):
        execute(prepared)
    config = json.loads((_output_dir(prepared) / "backend_config.json").read_text(encoding="utf-8"))
    assert "figures" not in config
    return config


def _output_dir(prepared) -> Path:
    return prepared.output_dir


def _annotation() -> dict:
    return {
        "organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID",
        "gene_symbol_output": True, "minimum_mapping_rate": 0.70, "minimum_mapped_foreground": 5,
        "enrichment": {
            "go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"},
            "gsea": {"minimum_ranked_genes": 61, "min_gs_size": 11, "max_gs_size": 321, "pvalue_cutoff": 0.041, "padj_cutoff": 0.049, "p_adjust_method": "BH", "seed": 9},
            "kegg": {
                "resource_provider": "online_kegg_rest_via_clusterprofiler",
                "ora": {"pvalue_cutoff": 0.021, "qvalue_cutoff": 0.19, "p_adjust_method": "BH", "min_gs_size": 12, "max_gs_size": 322},
                "gsea": {"minimum_ranked_genes": 62, "pvalue_cutoff": 0.91, "padj_cutoff": 0.039, "p_adjust_method": "BH", "min_gs_size": 13, "max_gs_size": 323, "seed": 10},
            },
        },
    }


def _project(project_factory, figures: dict | None, *, annotation: bool = False) -> Path:
    config = base_config()
    if figures is not None:
        config["figures"] = figures
    if annotation:
        config["annotation"] = _annotation()
    return project_factory(config=config)


CASES = [
    (None, LEGACY),
    ({"preset": "legacy"}, _explicit("legacy")),
    ({"preset": "standard"}, _explicit("standard")),
    ({"preset": "standard", "profile": "nature"}, _explicit("standard", "nature")),
]


@pytest.mark.parametrize(("figures", "expected"), CASES)
def test_direct_l1_carries_the_resolved_figures(monkeypatch, project_factory, figures, expected):
    monkeypatch.setattr("rnaseq.downstream._require_r", lambda: None)
    prepared = prepare_l1(validate_project(_project(project_factory, figures)), run_id=None)
    assert prepared.figures == expected
    # Provenance written from ``config`` must not change shape in Phase 1.
    assert "figures" not in prepared.config
    assert _capture_config(monkeypatch, "rnaseq.downstream", execute_l1, prepared)["figure_style"] == expected


@pytest.mark.parametrize(("figures", "expected"), CASES)
def test_direct_l2_carries_the_resolved_figures_and_matches_its_l1(monkeypatch, project_factory, figures, expected):
    monkeypatch.setattr("rnaseq.l2._require_r", lambda: None)
    prepared = prepare_l2(validate_project(_project(project_factory, figures)), run_id=None)
    monkeypatch.setattr("rnaseq.l2.execute_l1", lambda l1: SimpleNamespace(summary={}))
    config = _capture_config(monkeypatch, "rnaseq.l2", execute_l2, prepared)
    assert config["figure_style"] == expected == prepared.l1.figures


def _ready_l2(root: Path, *, ranked: bool = False) -> None:
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    (contrast / "all_genes.tsv").write_text("gene_id\tstat\n11303\t2.0\n", encoding="utf-8")
    for name in ("significant.tsv", "upregulated.tsv", "downregulated.tsv"):
        (contrast / name).write_text("gene_id\n11303\n", encoding="utf-8")


@pytest.mark.parametrize(("figures", "expected"), CASES)
def test_direct_go_gsea_and_kegg_carry_the_resolved_figures(monkeypatch, project_factory, figures, expected):
    root = _project(project_factory, figures, annotation=True)
    _ready_l2(root)
    monkeypatch.setattr("rnaseq.go._check_go_runtime", lambda _organism: None)
    monkeypatch.setattr("rnaseq.gsea._check_go_runtime", lambda _organism: None)
    monkeypatch.setattr("rnaseq.kegg._check_kegg_runtime", lambda: None)
    seen = {
        "go": _capture_config(monkeypatch, "rnaseq.go", execute_go, prepare_go(validate_project(root), run_id=None)),
        "gsea": _capture_config(monkeypatch, "rnaseq.gsea", execute_gsea, prepare_gsea(validate_project(root), run_id=None)),
        "kegg_ora": _capture_config(monkeypatch, "rnaseq.kegg", execute_kegg, prepare_kegg(validate_project(root), run_id=None, mode="ora")),
        "kegg_gsea": _capture_config(monkeypatch, "rnaseq.kegg", execute_kegg, prepare_kegg(validate_project(root), run_id=None, mode="gsea")),
    }
    assert {name: config["figure_style"] for name, config in seen.items()} == {name: expected for name in seen}


def test_direct_backend_configs_in_one_project_never_disagree(monkeypatch, project_factory):
    root = _project(project_factory, {"preset": "standard", "profile": "plos"}, annotation=True)
    _ready_l2(root)
    monkeypatch.setattr("rnaseq.downstream._require_r", lambda: None)
    monkeypatch.setattr("rnaseq.l2._require_r", lambda: None)
    monkeypatch.setattr("rnaseq.go._check_go_runtime", lambda _organism: None)
    report = validate_project(root)
    l1 = prepare_l1(report, run_id=None)
    l2 = prepare_l2(report, run_id=None)
    go = prepare_go(report, run_id=None)
    assert l1.figures == l2.l1.figures == go.l2.l1.figures == _explicit("standard", "plos")


# ------------------------------------------------------- Nextflow / frozen contract


def _freeze(project_factory, tmp_path, figures: dict | None, *, case: str = "CASE-FIG"):
    config = base_config()
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": "~ condition", "variables": {"condition": "categorical"}}
    config["annotation"] = _annotation()
    if figures is not None:
        config["figures"] = figures
    root = project_factory(config=config)
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    run = create_case_run(report, case, moment=datetime(2026, 10, 14, 12, 0, 0))
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    return root, run, frozen, resolve_downstream_inputs(run).root


@pytest.mark.parametrize(("figures", "expected"), CASES)
def test_frozen_contract_carries_resolved_figures_and_every_nextflow_config_agrees(
    project_factory, tmp_path, figures, expected
):
    _root, run, frozen, inputs = _freeze(project_factory, tmp_path, figures)
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
    assert contract["figures"] == expected
    l2 = run.run_dir / "downstream" / "l2"
    configs = [
        l1_config(frozen.contract, inputs, tmp_path / "l1"),
        l2_config(frozen.contract, inputs, tmp_path / "l1", tmp_path / "l2"),
        *(enrichment_config(frozen.contract, inputs, l2, kind, tmp_path) for kind in ("go", "kegg", "gsea-go", "gsea-kegg")),
    ]
    assert [config["figure_style"] for config in configs] == [expected] * 6
    assert all("figures" not in config for config in configs)  # figures.R would act on cfg$figures$profile


def test_an_old_frozen_contract_without_figures_resolves_to_legacy_in_every_config(project_factory, tmp_path):
    _root, run, frozen, inputs = _freeze(project_factory, tmp_path, None)
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
    del contract["figures"]  # exactly what a pre-figures run froze
    frozen.contract.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    l2 = run.run_dir / "downstream" / "l2"
    configs = [
        l1_config(frozen.contract, inputs, tmp_path / "l1"),
        l2_config(frozen.contract, inputs, tmp_path / "l1", tmp_path / "l2"),
        *(enrichment_config(frozen.contract, inputs, l2, kind, tmp_path) for kind in ("go", "gsea-go")),
    ]
    assert [config["figure_style"] for config in configs] == [LEGACY] * 4


def test_new_standard_project_in_an_old_contract_is_a_mixed_run_not_a_silent_legacy(project_factory, tmp_path):
    _root, run, frozen, inputs = _freeze(project_factory, tmp_path, {"preset": "standard"})
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
    del contract["figures"]
    frozen.contract.write_text(json.dumps(contract), encoding="utf-8")
    with pytest.raises(ValueError, match="mixed figure settings"):
        l1_config(frozen.contract, inputs, tmp_path / "l1")


@pytest.mark.parametrize("builder", ["l1", "l2", "go"])
def test_contract_and_staged_project_that_disagree_are_rejected(project_factory, tmp_path, builder):
    _root, run, frozen, inputs = _freeze(project_factory, tmp_path, {"preset": "standard"})
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
    contract["figures"] = _explicit("legacy")
    frozen.contract.write_text(json.dumps(contract), encoding="utf-8")
    call = {
        "l1": lambda: l1_config(frozen.contract, inputs, tmp_path / "l1"),
        "l2": lambda: l2_config(frozen.contract, inputs, tmp_path / "l1", tmp_path / "l2"),
        "go": lambda: enrichment_config(frozen.contract, inputs, run.run_dir / "downstream" / "l2", "go", tmp_path),
    }[builder]
    with pytest.raises(ValueError, match="mixed figure settings"):
        call()


def test_run_resolver_requires_agreement_on_profile_too():
    contract = {"figures": _explicit("standard", "nature")}
    assert resolve_run_figure_settings(contract, {"figures": {"preset": "standard", "profile": "nature"}}) == _explicit("standard", "nature")
    with pytest.raises(ValueError, match="mixed figure settings"):
        resolve_run_figure_settings(contract, {"figures": {"preset": "standard", "profile": "plos"}})
    with pytest.raises(ValueError, match="mixed figure settings"):
        resolve_run_figure_settings(contract, {})
    assert resolve_run_figure_settings({}, {}) == LEGACY


# --------------------------------------------------------------- retry and clone


def _frozen_files(run) -> dict[str, bytes]:
    frozen = run.run_dir / "frozen"
    return {name: (frozen / name).read_bytes() for name in ("project.yaml", "downstream_contract.json")}


@pytest.mark.usefixtures("production_capable_execution_capacity")
def test_retry_clone_keeps_the_source_figures_even_if_the_live_project_changes(monkeypatch, project_factory):
    config = base_config()
    config["figures"] = {"preset": "standard", "profile": "nature"}
    root = project_factory(config=config)
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-RETRY")
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run", str(root), "--case-id", "CASE-RETRY"])
    (run.run_dir / "provenance" / "run_provenance.yaml").write_text("design: {}\n", encoding="utf-8")
    from rnaseq.service import _write_state

    _write_state(run, "FAILED", phase="downstream", error="synthetic")
    assert json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text())["figures"] == _explicit("standard", "nature")
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)
    live = deepcopy(config)
    live["figures"] = {"preset": "legacy"}
    (root / "project.yaml").write_text(yaml.safe_dump(live, sort_keys=False), encoding="utf-8")

    retry = execute_retry_service_run(root, retry_of=f"{run.case_id}/{run.run_id}")

    cloned = json.loads((retry.run_dir / "frozen" / "downstream_contract.json").read_text())
    assert cloned["figures"] == _explicit("standard", "nature")
    assert (retry.run_dir / "frozen" / "project.yaml").read_bytes() == (run.run_dir / "frozen" / "project.yaml").read_bytes()


@pytest.mark.usefixtures("production_capable_execution_capacity")
def test_retry_of_an_old_run_without_figures_stays_legacy_and_does_not_inject_figures(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)  # base_config: no figures section
    path = source.run_dir / "frozen" / "downstream_contract.json"
    contract = json.loads(path.read_text(encoding="utf-8"))
    del contract["figures"]
    path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)

    retry = execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")

    cloned = json.loads((retry.run_dir / "frozen" / "downstream_contract.json").read_text())
    assert "figures" not in cloned
    assert figure_settings_from_contract(cloned) == LEGACY


@pytest.mark.usefixtures("production_capable_execution_capacity")
def test_retry_refuses_a_source_whose_contract_and_project_disagree(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)
    path = source.run_dir / "frozen" / "downstream_contract.json"
    contract = json.loads(path.read_text(encoding="utf-8"))
    contract["figures"] = _explicit("standard")  # frozen project.yaml has no figures => legacy
    path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)
    with pytest.raises(ExecutionPreflightError, match="mixed figure settings"):
        execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")


# ------------------------------------------------------- upstream reuse invariance


def test_changing_only_the_figures_selection_does_not_change_the_upstream_identity(project_factory, tmp_path):
    """The identity ``--reuse-upstream`` compares excludes project.yaml bytes/hash and ``figures``."""

    from rnaseq.service import _upstream_identity

    manifests = []
    for index, figures in enumerate((None, {"preset": "legacy"}, {"preset": "standard"}, {"preset": "standard", "profile": "elsevier"})):
        config = base_config()
        if figures is not None:
            config["figures"] = figures
        root = project_factory(config=config)
        report = validate_project(root)
        generate_plan(report)
        run = create_case_run(report, f"CASE-INV-{index}", moment=datetime(2026, 10, 14, 12, 0, 0))
        freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
        manifests.append(run.run_dir / "frozen" / "input_manifest.yaml")
    raw = [yaml.safe_load(path.read_text(encoding="utf-8")) for path in manifests]
    # The whole-file configuration hash does move with the selection ...
    assert len({item["configuration"]["sha256"] for item in raw}) > 1
    # ... but the compared identity does not.
    identities = [_upstream_identity(path) for path in manifests]
    assert all(identity == identities[0] for identity in identities)
    assert "figures" not in json.dumps(identities[0])


def test_upstream_reuse_succeeds_across_a_figures_change(project_factory):
    run, frozen, source = _reuse_pair(project_factory, source_manifest=_UPSTREAM_MANIFEST)
    # The earlier run was frozen with legacy figures; the new run asks for standard.
    (source / "frozen" / "project.yaml").write_text(yaml.safe_dump({"figures": {"preset": "legacy"}}), encoding="utf-8")
    (source / "frozen" / "downstream_contract.json").write_text(json.dumps({"figures": _explicit("legacy")}), encoding="utf-8")
    (run.run_dir / "frozen" / "project.yaml").write_text(yaml.safe_dump({"figures": {"preset": "standard"}}), encoding="utf-8")
    (run.run_dir / "frozen" / "downstream_contract.json").write_text(json.dumps({"figures": _explicit("standard")}), encoding="utf-8")
    from rnaseq.service import reuse_upstream_if_compatible

    assert reuse_upstream_if_compatible(run, frozen, _REUSE_SOURCE) == _REUSE_SOURCE


def test_nfcore_inputs_do_not_mention_figures(project_factory):
    """nfcore.params.json / config / samplesheet are byte-compared for reuse; figures must stay out of them."""

    config = base_config()
    config["figures"] = {"preset": "standard"}
    root = project_factory(config=config)
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-NOFIG", moment=datetime(2026, 10, 14, 12, 0, 0))
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    for name in ("nfcore.params.json", "nfcore.local.config", "samplesheet.csv", "nfcore.conda.config"):
        path = run.run_dir / "frozen" / name
        if path.is_file():
            assert "figures" not in path.read_text(encoding="utf-8").lower(), name


# ------------------------------------------------- preset rename: "bioinfo" is an alias of "standard"


@pytest.mark.parametrize("profile", ["general", "nature", "plos", "elsevier"])
def test_bioinfo_is_accepted_and_normalised_to_standard(profile):
    assert FiguresConfig(preset="bioinfo", profile=profile).model_dump() == {"preset": "standard", "profile": profile}
    assert resolve_figure_settings({"preset": "bioinfo", "profile": profile}) == _explicit("standard", profile)
    assert ProjectConfig.model_validate({**base_config(), "figures": {"preset": "bioinfo", "profile": profile}}).figures.preset == "standard"


def test_the_alias_is_the_only_extra_spelling():
    for spelling in ("Standard", "bio-info", "publication", "bioinfo "):
        with pytest.raises(ValidationError):
            FiguresConfig(preset=spelling)


def test_a_frozen_contract_written_with_bioinfo_still_resolves_and_is_never_rewritten():
    old = {"figures": {"preset": "bioinfo", "profile": "elsevier", "source": "explicit"}}
    assert figure_settings_from_contract(old) == _explicit("standard", "elsevier")
    assert old["figures"]["preset"] == "bioinfo"  # resolving reads; it does not rewrite the frozen file
    # contract and project may spell the same preset differently without being a mixed run
    assert resolve_run_figure_settings(old, {"figures": {"preset": "standard", "profile": "elsevier"}}) == _explicit("standard", "elsevier")
    assert resolve_run_figure_settings({"figures": _explicit("standard", "elsevier")}, {"figures": {"preset": "bioinfo", "profile": "elsevier"}}) == _explicit("standard", "elsevier")
    with pytest.raises(ValueError, match="mixed figure settings"):
        resolve_run_figure_settings(old, {"figures": {"preset": "legacy", "profile": "elsevier"}})


def test_r_receives_the_canonical_name_for_a_bioinfo_project(monkeypatch, project_factory):
    config = base_config()
    config["figures"] = {"preset": "bioinfo", "profile": "nature"}
    report = validate_project(project_factory(config=config))
    assert report.config.figures.preset == "standard"
    assert resolve_figure_settings(report.config.figures) == _explicit("standard", "nature")


@pytest.mark.usefixtures("production_capable_execution_capacity")
def test_retry_of_a_run_frozen_with_bioinfo_keeps_its_bytes_and_stays_consistent(monkeypatch, project_factory):
    config = base_config()
    config["figures"] = {"preset": "bioinfo", "profile": "nature"}
    root = project_factory(config=config)
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-ALIAS")
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run", str(root), "--case-id", "CASE-ALIAS"])
    (run.run_dir / "provenance" / "run_provenance.yaml").write_text("design: {}\n", encoding="utf-8")
    from rnaseq.service import _write_state

    _write_state(run, "FAILED", phase="downstream", error="synthetic")
    contract_path = run.run_dir / "frozen" / "downstream_contract.json"
    # Simulate a run frozen before the rename: the contract and the project both say "bioinfo".
    contract = json.loads(contract_path.read_text())
    contract["figures"]["preset"] = "bioinfo"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    frozen_project = run.run_dir / "frozen" / "project.yaml"
    frozen_project.write_text(frozen_project.read_text().replace("preset: standard", "preset: bioinfo"), encoding="utf-8")
    assert "preset: bioinfo" in frozen_project.read_text()
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)

    retry = execute_retry_service_run(root, retry_of=f"{run.case_id}/{run.run_id}")

    cloned = json.loads((retry.run_dir / "frozen" / "downstream_contract.json").read_text())
    assert cloned["figures"]["preset"] == "bioinfo"  # retained verbatim
    assert "preset: bioinfo" in (retry.run_dir / "frozen" / "project.yaml").read_text()
    assert figure_settings_from_contract(cloned) == _explicit("standard", "nature")


def test_new_l2_project_combines_the_standard_figure_default_with_the_enrichment_defaults(tmp_path):
    base = ["new", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts", "--preset", "L2", "--design-type", "two_group", "--scaffold", "--yes"]
    default = runner.invoke(app, [*base, "--name", "with_enrichment", "--gene-id-type", "ENSEMBL"])
    assert default.exit_code == 0, default.output
    raw = yaml.safe_load((tmp_path / "with_enrichment" / "project.yaml").read_text(encoding="utf-8"))
    assert raw["figures"] == {"preset": "standard", "profile": "general"}
    assert raw["analysis"]["enrichment"] == ["go", "kegg", "gsea"] and raw["annotation"]
    opted_out = runner.invoke(app, [*base, "--name", "without_enrichment", "--no-enrichment"])
    assert opted_out.exit_code == 0, opted_out.output
    raw = yaml.safe_load((tmp_path / "without_enrichment" / "project.yaml").read_text(encoding="utf-8"))
    assert raw["figures"] == {"preset": "standard", "profile": "general"} and raw["analysis"]["enrichment"] == []
    # An existing project without the section keeps the legacy style on the general profile.
    del raw["figures"]
    assert ProjectConfig.model_validate(raw).figures is None
