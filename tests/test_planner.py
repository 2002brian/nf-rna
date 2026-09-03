from __future__ import annotations

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
