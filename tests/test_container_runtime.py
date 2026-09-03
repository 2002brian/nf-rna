from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_runtime_image_declares_procps_and_build_time_prerequisite_probe():
    environment = (ROOT / "environment.yml").read_text(encoding="utf-8")
    dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
    assert "- procps-ng" in environment
    for command in ("command -v ps", "command -v python", "command -v Rscript"):
        assert command in dockerfile
    for package in ("DESeq2", "tximport", "clusterProfiler", "org.Mm.eg.db"):
        assert package in dockerfile
