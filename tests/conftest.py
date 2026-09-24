from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

import pytest
import yaml


BASE_COUNTS = """gene_id,C1,C2,C3,T1,T2,T3
GeneA,10,12,9,40,45,43
GeneB,100,110,98,95,102,99
"""
BASE_METADATA = """sample_id,condition
C1,Control
C2,Control
C3,Control
T1,Treatment
T2,Treatment
T3,Treatment
"""
BASE_CONTRASTS = """contrast_id,factor,numerator,denominator
Treatment_vs_Control,condition,Treatment,Control
"""


IDENTIFIED_TEST_REVISION = "0123456789abcdef0123456789abcdef01234567"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "real_source_revision: use the real nf-rna source-revision lookup instead of a clean test revision",
    )


@pytest.fixture(autouse=True)
def identified_nf_rna_source(request, monkeypatch):
    """Keep tests independent of this checkout's git state.

    Runs refuse an unidentified or dirty nf-rna source.  Tests that exercise
    that guard opt out with ``@pytest.mark.real_source_revision``.
    """

    if request.node.get_closest_marker("real_source_revision"):
        return
    monkeypatch.setattr("rnaseq.downstream_runtime.runtime_source_revision", lambda: IDENTIFIED_TEST_REVISION)
    monkeypatch.setattr("rnaseq.service.runtime_source_revision", lambda: IDENTIFIED_TEST_REVISION)


def base_config() -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "project": {"id": "test_project", "pipeline": "bulk_rnaseq", "preset": "L2"},
        "organism": {"species": "Mus musculus"},
        "input": {"type": "raw_counts", "path": "input/counts.csv"},
        "design": {"type": "two_group", "formula": "~ condition"},
        "metadata_file": "metadata.csv",
        "contrasts_file": "contrasts.csv",
        "upstream": {
            "engine": "external",
            "provider": "external_provider",
            "quantification_method": "unknown",
        },
        "reference": {"source": "igenomes", "genome": None},
        "thresholds": {"padj": 0.05, "abs_log2fc": 1.0},
    }


@pytest.fixture
def project_factory(tmp_path: Path) -> Callable[..., Path]:
    counter = 0

    def create(
        *,
        config: dict[str, Any] | None = None,
        counts: str = BASE_COUNTS,
        metadata: str = BASE_METADATA,
        contrasts: str = BASE_CONTRASTS,
    ) -> Path:
        nonlocal counter
        counter += 1
        root = tmp_path / f"project_{counter}"
        (root / "input").mkdir(parents=True)
        (root / "input" / "counts.csv").write_text(counts, encoding="utf-8")
        (root / "metadata.csv").write_text(metadata, encoding="utf-8")
        (root / "contrasts.csv").write_text(contrasts, encoding="utf-8")
        (root / "project.yaml").write_text(
            yaml.safe_dump(config or base_config(), sort_keys=False), encoding="utf-8"
        )
        return root

    return create


@pytest.fixture
def production_capable_execution_capacity(monkeypatch):
    """Make mocked execution tests independent of the host's actual capacity.

    Production code continues to probe the real host. Tests that mock runtime
    tools but are not testing resource calculation opt in to this fixture so a
    four-core CI runner can reach its intended assertion.
    """

    from rnaseq.execution import LocalResourceCapacity, RuntimeSnapshot

    capacity = LocalResourceCapacity(16, 64, 60)
    snapshot = RuntimeSnapshot(
        host_os="Linux",
        host_architecture="amd64",
        logical_cpus=16,
        host_memory_bytes=64 * 1024**3,
        docker_architecture="amd64",
        docker_memory_bytes=64 * 1024**3,
        docker_version="test",
        first_party_image_architecture="amd64",
        docker_cpus=16,
    )
    monkeypatch.setattr("rnaseq.execution.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.service.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.execution.runtime_snapshot", lambda *_args: snapshot)
    monkeypatch.setattr("rnaseq.service.runtime_snapshot", lambda *_args: snapshot)
    return snapshot


def require_rscript() -> str:
    executable = shutil.which("Rscript")
    if executable is None:
        pytest.skip("Rscript unavailable")
    return executable


def require_r_packages(*packages: str) -> str:
    """Skip full statistical-runtime tests while retaining base-R helper tests."""

    executable = require_rscript()
    expression = " && ".join(f'requireNamespace("{package}", quietly=TRUE)' for package in packages)
    result = subprocess.run([executable, "-e", f"quit(status=if ({expression}) 0 else 1)"], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        pytest.skip("R packages unavailable: " + ", ".join(packages))
    return executable
