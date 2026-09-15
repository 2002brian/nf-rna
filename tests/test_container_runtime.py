from __future__ import annotations

import shlex
import tomllib
from pathlib import Path

import yaml


ROOT = Path(__file__).parents[1]


def _copy_sources(dockerfile: str) -> set[str]:
    """Return root-relative sources explicitly supplied to Docker build stages."""

    sources: set[str] = set()
    for line in dockerfile.splitlines():
        if not line.startswith("COPY "):
            continue
        fields = [field for field in shlex.split(line) if not field.startswith("--")]
        sources.update(fields[1:-1])
    return sources


def test_container_runtime_path_contract_uses_non_login_shell():
    environment = (ROOT / "environment.yml").read_text(encoding="utf-8")
    docker_environment = yaml.safe_load((ROOT / "environment.docker.yml").read_text(encoding="utf-8"))
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
    dockerignore = (ROOT / ".dockerignore").read_text(encoding="utf-8")
    # Keep the host environment portable; procps is a Linux image contract.
    assert "procps-ng" not in environment
    assert "-e .[dev]" not in environment
    assert docker_environment["dependencies"] == ["procps-ng"]
    assert "!environment.docker.yml" in dockerignore
    assert "micromamba install --yes --name rnaseq --file environment.docker.yml" in dockerfile
    assert 'ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"' in dockerfile
    assert dockerfile.index('ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"') < dockerfile.index("RUN micromamba create")
    assert "micromamba run --name rnaseq" not in dockerfile
    assert "python -m pip install --no-deps ." in dockerfile
    assert "cp -a src/rnaseq/r/. /opt/nf-rna/r/" in dockerfile
    assert "org.opencontainers.image.revision" in dockerfile
    assert "sh -lc" not in dockerfile
    for command in ("command -v ps", "ps --version", "command -v python", "command -v Rscript"):
        assert command in dockerfile
    for package in ("DESeq2", "tximport", "clusterProfiler", "org.Mm.eg.db"):
        assert package in dockerfile


def test_docker_packaging_context_supplies_every_hatch_force_include_input():
    """Keep image installation equivalent to the repository packaging contract."""

    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
    dockerignore = (ROOT / ".dockerignore").read_text(encoding="utf-8")
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    force_include = pyproject["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    required = {"src", "pyproject.toml", "README.md", *force_include}
    copied = _copy_sources(dockerfile)

    assert required <= copied
    for directory in force_include:
        assert f"!{directory}/" in dockerignore
        assert f"!{directory}/**" in dockerignore
