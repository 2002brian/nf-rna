from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).parents[1]


def test_container_runtime_path_contract_uses_non_login_shell():
    environment = (ROOT / "environment.yml").read_text(encoding="utf-8")
    docker_environment = yaml.safe_load((ROOT / "environment.docker.yml").read_text(encoding="utf-8"))
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
    dockerignore = (ROOT / ".dockerignore").read_text(encoding="utf-8")
    # Keep the host environment portable; procps is a Linux image contract.
    assert "procps-ng" not in environment
    assert docker_environment["dependencies"] == ["procps-ng"]
    assert "!environment.docker.yml" in dockerignore
    assert "micromamba install --yes --name rnaseq --file environment.docker.yml" in dockerfile
    assert 'ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"' in dockerfile
    assert dockerfile.index('ENV PATH="/opt/conda/envs/rnaseq/bin:${PATH}"') < dockerfile.index("RUN micromamba create")
    assert "micromamba run --name rnaseq" not in dockerfile
    assert "sh -c 'command -v ps" in dockerfile
    assert "sh -lc" not in dockerfile
    for command in ("command -v ps", "ps --version", "command -v python", "command -v Rscript"):
        assert command in dockerfile
    for package in ("DESeq2", "tximport", "clusterProfiler", "org.Mm.eg.db"):
        assert package in dockerfile
