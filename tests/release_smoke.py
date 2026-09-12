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

        expected_version = "1.0.0rc1"
        version_result = _run([str(executable), "--version"], cwd=temporary)
        if version_result.returncode != 0 or version_result.stdout.strip() != expected_version:
            raise SystemExit(version_result.stderr or version_result.stdout)
        if rnaseq.__version__ != expected_version or version("nf-rna") != expected_version:
            raise SystemExit("Installed package and distribution metadata disagree on the RC version.")

        assets = required_workflow_assets()
        if not all(path.is_file() for path in assets.values()):
            raise SystemExit("A required first-party workflow asset is missing from the installed distribution.")

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
