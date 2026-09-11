from __future__ import annotations

import yaml

from conftest import base_config
from rnaseq.planner import generate_plan
from rnaseq.validators import validate_project


def test_gsea_plan_reports_only_the_public_method_and_its_gene_set_resources(project_factory):
    config = base_config()
    config["schema_version"] = "1.1"
    config["annotation"] = {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}
    config["analysis"] = {"enrichment": "gsea"}
    root = project_factory(config=config)
    report = validate_project(root)
    generate_plan(report)
    plan = (root / "planning" / "analysis_plan.md").read_text(encoding="utf-8")
    manifest = (root / "planning" / "manifest.preview.yaml").read_text(encoding="utf-8")
    assert "Enrichment method: `GSEA`" in plan
    assert "`GO BP`, `GO MF`, `GO CC`, `KEGG`" in plan
    assert "go, gsea-go, kegg, gsea-kegg" not in plan
    assert "enrichment:\n  - gsea\n" in manifest


def test_plan_exposes_requested_and_effective_runtime_budget(monkeypatch, project_factory):
    from rnaseq.execution import RuntimeSnapshot

    config = base_config()
    config["execution"] = {"profile": "local", "max_cpus": 16, "max_memory_gb": 32}
    root = project_factory(config=config)
    monkeypatch.setattr("rnaseq.execution.runtime_snapshot", lambda *_args: RuntimeSnapshot(
        "Darwin", "arm64", 12, 24 * 1024**3, "arm64", 16 * 1024**3, "test", "arm64", 10
    ))
    generate_plan(validate_project(root))
    resources = yaml.safe_load((root / "planning" / "resource_plan.yaml").read_text(encoding="utf-8"))
    assert resources["requested"] == {"cpus": 16, "memory_gib": 32}
    assert resources["effective"] == {"cpus": 10, "memory_gib": 16}
    assert resources["clamped"] is True
    assert "aggregate" in resources["scheduling"]["policy"]
