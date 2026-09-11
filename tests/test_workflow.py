from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from rnaseq.workflow_support import _ENRICHMENT_REQUIRED_PATHS
from rnaseq.execution import ResourceContract, render_local_resource_config


ROOT = Path(__file__).parents[1]
WORKFLOW = ROOT / "workflow"


def test_downstream_workflow_allows_independent_gsea_backends_to_share_the_budget():
    text = (WORKFLOW / "main.nf").read_text(encoding="utf-8")
    match = re.search(r"process ENRICHMENT_ANALYSIS \{(?P<body>.*?)^\}", text, flags=re.DOTALL | re.MULTILINE)
    assert match is not None
    body = match.group("body")
    assert "maxForks" not in body
    assert 'publishDir "${params.outdir}/l2", mode: \'copy\', overwrite: false' in body
    assert "publishDir params.outdir, mode: 'copy', overwrite: false" not in body
    assert "tuple val(module), path('enrichment/*')" in body
    assert "path 'enrichment'" not in body
    assert "Channel.fromList(modules).combine(l2)" in text
    assert ".map { module, directory -> directory }" in text
    assert ".collect()" in text
    assert "TECHNICAL_REPORT(l1, l2, enrichment_dirs, contract, inputs)" in text
    assert "TECHNICAL_REPORT_NO_ENRICHMENT(l1, l2, contract, inputs)" in text
    assert "TECHNICAL_REPORT_L1(l1, contract, inputs)" in text
    assert "--enrichment $enrichment_dirs" in text
    assert set(_ENRICHMENT_REQUIRED_PATHS) == {"gsea-go", "gsea-kegg"}
    assert "allowedEnrichment = ['gsea-go', 'gsea-kegg']" in text
    assert "go_analysis.R" not in text


def test_downstream_workflow_keeps_l1_l2_and_report_persistence_contracts():
    text = (WORKFLOW / "main.nf").read_text(encoding="utf-8")
    for name, output in (("L1_ANALYSIS", "l1"), ("L2_ANALYSIS", "l2"), ("TECHNICAL_REPORT", "report"), ("TECHNICAL_REPORT_L1", "report")):
        match = re.search(rf"process {name} \{{(?P<body>.*?)^\}}", text, flags=re.DOTALL | re.MULTILINE)
        assert match is not None
        body = match.group("body")
        assert "publishDir params.outdir, mode: 'copy', overwrite: false" in body
        assert "stageInMode 'copy'" in body
        assert "path inputs" in body
        assert f"path '{output}'" in body


def test_downstream_workflow_stages_a_narrow_execution_input_bundle_without_volume_mounts():
    text = (WORKFLOW / "main.nf").read_text(encoding="utf-8")
    assert "params.inputs = null" in text
    assert "params.analysis_level = null" in text
    assert "Specify --contract, --inputs, --outdir and --analysis_level" in text
    assert "analysis_level must be L1 or L2" in text
    assert "inputs = Channel.value(file(params.inputs))" in text
    assert "--inputs $inputs" in text
    assert "/Volumes/KOXIA" not in text
    for name in ("L1_ANALYSIS", "L2_ANALYSIS", "TECHNICAL_REPORT", "TECHNICAL_REPORT_NO_ENRICHMENT", "TECHNICAL_REPORT_L1", "ENRICHMENT_ANALYSIS"):
        match = re.search(rf"process {name} \{{(?P<body>.*?)^\}}", text, flags=re.DOTALL | re.MULTILINE)
        assert match is not None
        body = match.group("body")
        assert "path inputs" in body
        assert "stageInMode 'copy'" in body


@pytest.mark.parametrize("profile", ("local", "docker", "server"))
def test_downstream_nextflow_profiles_parse(profile: str):
    if shutil.which("nextflow") is None:
        pytest.skip("Nextflow unavailable")
    result = subprocess.run(
        ["nextflow", "config", "-profile", profile, str(WORKFLOW)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_all_downstream_processes_use_the_doctor_checked_runtime_image():
    config = (WORKFLOW / "nextflow.config").read_text(encoding="utf-8")
    assert "process.container = 'rnaseq-control-plane:latest'" in config
    assert "rnaseq-downstream:latest" not in config


def test_downstream_resource_contracts_are_explicit_without_artificial_serialization():
    config = (WORKFLOW / "nextflow.config").read_text(encoding="utf-8")
    assert "executor {\n  cpus = 8\n  memory = '12 GB'\n}" in config
    assert "resourceLimits = [cpus: 8, memory: '12 GB', time: '12 h']" in config
    for process in ("L1_ANALYSIS", "L2_ANALYSIS", "ENRICHMENT_ANALYSIS", "TECHNICAL_REPORT", "TECHNICAL_REPORT_L1"):
        assert f"withName: {process}" in config
    assert "withName: L1_ANALYSIS { cpus = 1; memory = '2 GB'; time = '2 h' }" in config
    assert "withName: ENRICHMENT_ANALYSIS { cpus = 4; memory = '8 GB'; time = '8 h' }" in config
    assert "maxForks" not in config
    main = (WORKFLOW / "main.nf").read_text(encoding="utf-8")
    assert not re.search(r"process ENRICHMENT_ANALYSIS \{(?P<body>.*?)maxForks", main, flags=re.DOTALL)


def test_sample_level_hisat2_workflow_has_no_fixed_parallelism_cap():
    text = (WORKFLOW / "hisat2_featurecounts.nf").read_text(encoding="utf-8")
    assert "maxForks" not in text


def test_frozen_effective_budget_overrides_only_aggregate_nextflow_limits(tmp_path):
    if shutil.which("nextflow") is None:
        pytest.skip("Nextflow unavailable")
    frozen = tmp_path / "effective.config"
    frozen.write_text(render_local_resource_config(ResourceContract("EFFECTIVE_LOCAL", 16, 32, 12)), encoding="utf-8")
    result = subprocess.run(
        ["nextflow", "-c", str(frozen), "config", "-profile", "local", str(WORKFLOW)],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "cpus = 16" in result.stdout and "memory = '32.GB'" in result.stdout
    assert "withName:ENRICHMENT_ANALYSIS" in result.stdout and "cpus = 4" in result.stdout
    assert "maxForks" not in result.stdout


def test_downstream_workflow_uses_the_frozen_analysis_level_to_gate_l2_and_gsea():
    text = (WORKFLOW / "main.nf").read_text(encoding="utf-8")
    assert "if( params.analysis_level == 'L1' )" in text
    l1_branch = text.split("if( params.analysis_level == 'L1' )", 1)[1].split("} else {", 1)[0]
    assert "TECHNICAL_REPORT_L1(l1, contract, inputs)" in l1_branch
    assert "L2_ANALYSIS" not in l1_branch
    assert "ENRICHMENT_ANALYSIS" not in l1_branch
    l2_branch = text.split("} else {", 1)[1]
    assert "l2 = L2_ANALYSIS(l1, contract, inputs)" in l2_branch
    assert "ENRICHMENT_ANALYSIS" in l2_branch
    assert "TECHNICAL_REPORT_NO_ENRICHMENT" in l2_branch
