"""Validate a separately installed distribution without importing this checkout."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def _run(command: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)


def main() -> None:
    import rnaseq
    from importlib.metadata import version
    from rnaseq.execution import resolve_execution_workspace
    from rnaseq.models import RuntimeConfig, execution_image_for_version
    from rnaseq.project import load_project
    from rnaseq.service import _provenance, create_case_run, freeze_case_inputs
    from rnaseq.validators import validate_project
    from rnaseq.workflow_assets import required_workflow_assets

    source_root = os.environ.get("NF_RNA_SOURCE_ROOT")
    package_root = Path(rnaseq.__file__).resolve()
    if source_root and package_root.is_relative_to(Path(source_root).resolve()):
        raise SystemExit(f"Imported nf-rna from source checkout instead of installed artifact: {package_root}")

    temporary = Path(tempfile.mkdtemp(prefix="nf-rna-release-smoke-"))
    try:
        executable = Path(sys.executable).resolve().parent / "rnaseq"
        help_result = _run([str(executable), "--help"], cwd=temporary)
        if help_result.returncode != 0 or "Create, validate, and plan" not in help_result.stdout:
            raise SystemExit(help_result.stderr or help_result.stdout)

        expected_version = rnaseq.__version__
        version_result = _run([str(executable), "--version"], cwd=temporary)
        if version_result.returncode != 0 or version_result.stdout.strip() != expected_version:
            raise SystemExit(version_result.stderr or version_result.stdout)
        if rnaseq.__version__ != expected_version or version("nf-rna") != expected_version:
            raise SystemExit("Installed package and distribution metadata disagree on the release version.")
        expected_image = execution_image_for_version(expected_version)
        if RuntimeConfig().execution_image != expected_image:
            raise SystemExit("Installed package did not retain its version-matched default execution image.")

        assets = required_workflow_assets()
        if not all(path.is_file() for path in assets.values()):
            raise SystemExit("A required first-party workflow asset is missing from the installed distribution.")

        project_result = _run([
            str(executable), "new", "--name", "installed-project", "--destination", str(temporary),
            "--species", "mouse", "--input-type", "raw_counts", "--preset", "L2",
            "--design-type", "two_group", "--scaffold", "--yes",
        ], cwd=temporary)
        if project_result.returncode != 0:
            raise SystemExit(project_result.stderr or project_result.stdout)
        project = temporary / "installed-project"
        (project / "input" / "counts.csv").write_text(
            "gene_id,C1,C2,T1,T2\nGeneA,10,12,20,22\nGeneB,30,31,29,28\n",
            encoding="utf-8",
        )
        (project / "metadata.csv").write_text(
            "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n",
            encoding="utf-8",
        )
        (project / "contrasts.csv").write_text(
            "contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n",
            encoding="utf-8",
        )
        loaded = load_project(project)
        if loaded.config.runtime.execution_image != expected_image:
            raise SystemExit("Installed CLI created a project with the wrong default execution image.")
        report = validate_project(project)
        run = create_case_run(report, "SMOKE-INSTALLED")
        freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
        provenance = _provenance(
            report, run, profile="local", command=["rnaseq", "run"],
            workspace=resolve_execution_workspace(run.case_id, run.run_id),
        )
        if provenance["source_checkout"] is not None or provenance["git_commit"] is not None:
            raise SystemExit("Installed package incorrectly claimed a source checkout in provenance.")

        config_result = _run(["nextflow", "config", "-flat", str(assets["nextflow.config"])], cwd=temporary)
        if config_result.returncode != 0:
            raise SystemExit(config_result.stderr or config_result.stdout)

        for name in ("main.nf", "hisat2_featurecounts.nf"):
            result = _run(
                ["nextflow", "run", str(assets[name]), "-c", str(assets["nextflow.config"])],
                cwd=temporary,
            )
            combined = result.stdout + result.stderr
            if result.returncode == 0 or "Specify" not in combined:
                raise SystemExit(f"{name} did not reach its expected parameter validation:\n{combined}")
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
