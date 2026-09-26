"""Regression coverage for the first-party native downstream runtime source."""

from __future__ import annotations

import hashlib
import re
import subprocess
import tomllib
from pathlib import Path

import yaml

from rnaseq.downstream_runtime import _source_revision


ROOT = Path(__file__).parents[1]
RUNTIME_SOURCE = ROOT / "workflow" / "envs" / "nf-rna-downstream.yml"
LOCKS = (
    ROOT / "workflow" / "envs" / "locks" / "nf-rna-downstream-linux-64.lock.yml",
    ROOT / "workflow" / "envs" / "locks" / "nf-rna-downstream-osx-arm64.lock.yml",
)
CHECKSUMS = ROOT / "workflow" / "envs" / "locks" / "SHA256SUMS"


def test_nextflow_root_temp_is_ignored_without_hiding_tracked_source_changes(tmp_path):
    source = tmp_path / "checkout"
    source.mkdir()
    (source / ".gitignore").write_bytes((ROOT / ".gitignore").read_bytes())
    tracked = source / "source.py"
    tracked.write_text("version = 1\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=source, check=True)
    subprocess.run(["git", "add", ".gitignore", "source.py"], cwd=source, check=True)
    subprocess.run(
        ["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "baseline"],
        cwd=source, check=True,
    )
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()

    (source / "nxf-tmp.pcJ4CY").touch()
    assert _source_revision(source) == revision

    tracked.write_text("version = 2\n", encoding="utf-8")
    assert _source_revision(source) == revision + "+dirty"


def _distribution_name(specification: str) -> str:
    return re.split(r"[<>=!~;\[]", specification, maxsplit=1)[0].strip().lower().replace("_", "-")


def test_downstream_runtime_source_contains_only_first_party_runtime_dependencies():
    runtime = yaml.safe_load(RUNTIME_SOURCE.read_text(encoding="utf-8"))
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    dependencies = runtime["dependencies"]
    supplied = {_distribution_name(item) for item in dependencies}
    declared_python = {_distribution_name(item) for item in pyproject["project"]["dependencies"]}

    assert runtime["name"] == "nf-rna-downstream"
    assert runtime["channels"] == ["conda-forge", "bioconda"]
    assert declared_python <= supplied
    assert {
        "r-base",
        "bioconductor-deseq2",
        "bioconductor-tximport",
        "bioconductor-clusterprofiler",
        "bioconductor-annotationdbi",
        "bioconductor-org.hs.eg.db",
        "bioconductor-org.mm.eg.db",
        "r-ggplot2",
        "r-pheatmap",
        "r-yaml",
        "r-jsonlite",
    } <= supplied

    # These belong to Docker-only support, upstream nf-core execution, or the
    # separate HISAT2/featureCounts workflow, not the downstream R runtime.
    assert not {
        "procps-ng",
        "salmon",
        "hisat2",
        "samtools",
        "subread",
        "fastp",
        "fastqc",
        "multiqc",
        "pytest",
    } & supplied


def test_downstream_runtime_locks_are_exact_platform_resolutions_of_the_source():
    source = yaml.safe_load(RUNTIME_SOURCE.read_text(encoding="utf-8"))
    source_names = {_distribution_name(item) for item in source["dependencies"]}

    for lock_path in LOCKS:
        lock = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
        locked = lock["dependencies"]
        locked_names = {_distribution_name(item.rsplit("::", maxsplit=1)[-1]) for item in locked}

        assert lock["name"] == source["name"]
        assert lock["channels"] == source["channels"]
        assert source_names <= locked_names
        assert all("::" in item and "==" in item for item in locked)


def test_downstream_runtime_lock_checksums_are_current():
    expected = {
        filename: digest
        for digest, filename in (
            line.split(maxsplit=1)
            for line in CHECKSUMS.read_text(encoding="utf-8").splitlines()
            if line
        )
    }

    assert set(expected) == {path.name for path in LOCKS}
    for lock_path in LOCKS:
        actual = hashlib.sha256(lock_path.read_bytes()).hexdigest()
        assert actual == expected[lock_path.name]
