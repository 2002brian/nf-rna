"""All four enrichment backends share one dual-threshold annotation-mapping QC.

``mapping_warning_rate`` only warns; ``minimum_mapping_rate`` alone blocks.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

from conftest import require_r_packages


R_ROOT = Path(__file__).resolve().parents[1] / "src" / "rnaseq" / "r"
# Mouse Entrez IDs present in org.Mm.eg.db, and IDs that are not.
MAPPED = ("11287", "11298", "11303", "11304", "11305", "11306", "11307", "11308")
UNMAPPED = tuple(f"99999999{index}" for index in range(10))
# Ten tested genes per contrast; the mapped fraction sets the mapping rate.
RATES = {"pass": (0.80, "PASS"), "warning": (0.60, "WARNING"), "blocked": (0.40, "BLOCKED")}


def test_no_backend_compares_a_mapping_rate_with_a_threshold_directly():
    for script in ("go_analysis.R", "kegg_analysis.R", "gsea_analysis.R"):
        text = (R_ROOT / script).read_text(encoding="utf-8")
        assert 'source(file.path(dirname(normalizePath(script_file)), "annotation_mapping_qc.R"))' in text
        assert not re.search(r"mapping_rate\s*<\s*cfg\$annotation", text), script
    kegg = (R_ROOT / "kegg_analysis.R").read_text(encoding="utf-8")
    # Both KEGG modes (ORA and GSEA) call the shared policy.
    assert kegg.count("annotation_mapping_qc(") == 2


def _contrast_files(root: Path, name: str, mapped: int) -> dict[str, str]:
    genes = MAPPED[:mapped] + UNMAPPED[: 10 - mapped]
    directory = root / name
    directory.mkdir()
    rows = "".join(f"{gene}\t{0.001 * (index + 1)}\t{10 - index}.5\t0.01\n" for index, gene in enumerate(genes))
    (directory / "all_genes.tsv").write_text("gene_id\tpvalue\tstat\tpadj\n" + rows, encoding="utf-8")
    (directory / "significant.tsv").write_text("gene_id\n" + "".join(f"{gene}\n" for gene in genes), encoding="utf-8")
    (directory / "empty.tsv").write_text("gene_id\n", encoding="utf-8")
    return {
        "contrast_id": name, "all_genes": str(directory / "all_genes.tsv"), "significant": str(directory / "significant.tsv"),
        "up": str(directory / "empty.tsv"), "down": str(directory / "empty.tsv"),
    }


def _run_backend(tmp_path: Path, script: str, mode: str | None) -> dict:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    contrasts = [_contrast_files(inputs, name, round(rate * 10)) for name, (rate, _status) in RATES.items()]
    probe = tmp_path / "probe.txt"
    probe.write_text("path:mmu00010\tGlycolysis / Gluconeogenesis\n", encoding="utf-8")
    # Gates placed after mapping QC make "QC passed, backend continued" observable
    # without running a full enrichment or contacting KEGG.
    unreachable = 1000
    config = {
        "output_dir": str(tmp_path / "out"), "orgdb_package": "org.Mm.eg.db", "contrasts": contrasts,
        "annotation": {
            "organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID",
            "mapping_warning_rate": 0.70, "minimum_mapping_rate": 0.50, "minimum_mapped_foreground": unreachable,
            "enrichment": {
                "go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"},
                "gsea": {"minimum_ranked_genes": unreachable, "min_gs_size": 10, "max_gs_size": 500, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "seed": 1},
                "kegg": {
                    "resource_provider": "online_kegg_rest_via_clusterprofiler",
                    "ora": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500},
                    "gsea": {"minimum_ranked_genes": unreachable, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500, "seed": 1},
                },
            },
        },
    }
    if mode is not None:
        config.update({"mode": mode, "kegg": {"organism_code": "mmu", "provider": "online_kegg_rest_via_clusterprofiler", "probe_endpoint": probe.as_uri()}})
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    result = subprocess.run(["Rscript", str(R_ROOT / script), "--config", str(config_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    summary_name = {
        ("go_analysis.R", None): "go_backend_summary.json", ("kegg_analysis.R", "ora"): "kegg_backend_summary.json",
        ("gsea_analysis.R", None): "gsea_backend_summary.json", ("kegg_analysis.R", "gsea"): "gsea_kegg_backend_summary.json",
    }[(script, mode)]
    return json.loads((tmp_path / "out" / summary_name).read_text(encoding="utf-8"))


@pytest.mark.parametrize(("script", "mode"), (("go_analysis.R", None), ("kegg_analysis.R", "ora")), ids=["go-ora", "kegg-ora"])
def test_ora_backends_block_only_below_minimum_mapping_rate(tmp_path, script, mode):
    require_r_packages("jsonlite", "AnnotationDbi", "clusterProfiler", "ggplot2", "org.Mm.eg.db")
    summary = _run_backend(tmp_path, script, mode)
    by_contrast = {item["contrast_id"]: item for item in summary["contrasts"]}
    for name, (rate, status) in RATES.items():
        universe = by_contrast[name]["universe"]
        assert universe["mapping_rate"] == pytest.approx(rate)
        for outcome in by_contrast[name]["foregrounds"].values():
            if status == "BLOCKED":
                assert outcome["status"] == "BLOCKED"
                assert "below blocking threshold 50.0%" in outcome["reason"]
            else:
                # QC passed (PASS or WARNING): ORA continued to its foreground-size gate.
                assert outcome["status"] == "NOT_APPLICABLE", (name, outcome)
    for name, (_rate, status) in RATES.items():
        assert by_contrast[name]["universe"]["annotation_qc"]["status"] == status
    assert summary["annotation"]["warning_threshold"] == 0.70
    assert summary["annotation"]["blocking_threshold"] == 0.50
    assert summary["annotation_qc_status"] == "BLOCKED"


@pytest.mark.parametrize(("script", "mode"), (("gsea_analysis.R", None), ("kegg_analysis.R", "gsea")), ids=["go-gsea", "kegg-gsea"])
def test_gsea_backends_block_only_below_minimum_mapping_rate(tmp_path, script, mode):
    require_r_packages("jsonlite", "AnnotationDbi", "clusterProfiler", "ggplot2", "org.Mm.eg.db")
    summary = _run_backend(tmp_path, script, mode)
    assert summary["annotation"]["warning_threshold"] == 0.70
    assert summary["annotation"]["blocking_threshold"] == 0.50
    assert summary["annotation_qc_status"] == "BLOCKED"
    by_contrast = {item["contrast_id"]: item for item in summary["contrasts"]}
    for name, (rate, status) in RATES.items():
        item = by_contrast[name]
        assert item["ranking"]["mapping_rate"] == pytest.approx(rate)
        assert item["ranking"]["annotation_qc"]["status"] == status
        if status == "BLOCKED":
            assert "below blocking threshold 50.0%" in item["reason"]
        else:
            # QC passed (PASS or WARNING): GSEA continued to its ranked-gene gate.
            assert "unique ranked Entrez genes" in item["reason"], item["reason"]
