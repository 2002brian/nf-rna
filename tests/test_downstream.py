from __future__ import annotations

import hashlib
import json
import shutil
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from conftest import BASE_CONTRASTS, base_config, require_r_packages
from rnaseq.downstream import execute_l1, prepare_l1
from rnaseq.errors import DownstreamExecutionError
from rnaseq.validators import validate_project


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_raw_counts_l1_real_backend_is_structured_and_text_deterministic(project_factory):
    pytest.importorskip("yaml")
    require_r_packages("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap")
    counts = "gene_id,C1,C2,C3,T1,T2,T3\n" + "\n".join(
        f"Gene{index},{10 + index},{12 + index},{9 + index},{40 + index},{45 + index},{43 + index}"
        for index in range(1, 101)
    ) + "\n"
    report = validate_project(project_factory(counts=counts))
    prepared = prepare_l1(report, run_id=None)
    first = execute_l1(prepared)
    assert json.loads(first.state_path.read_text())["status"] == "SUCCESS"
    assert first.summary["source_type"] == "raw_counts"
    assert first.summary["genes_retained"] == 100
    normalized = first.output_dir / "normalized_counts.tsv"
    assert normalized.read_text().splitlines()[0] == "gene_id\tC1\tC2\tC3\tT1\tT2\tT3"
    assert all((first.output_dir / name).is_file() for name in ("vst.tsv", "pca.png", "pca.tiff", "sample_correlation.png", "sample_correlation.tiff", "filtering_summary.yaml", "l1_report.md"))
    assert not list(first.output_dir.glob("*.svg"))
    assert not list(first.output_dir.glob("*.pdf"))
    before = {_path.name: _digest(_path) for _path in (normalized, first.output_dir / "vst.tsv", first.output_dir / "pca_scores.tsv", first.output_dir / "sample_correlation.tsv")}
    execute_l1(prepared)
    assert before == {_path.name: _digest(_path) for _path in (normalized, first.output_dir / "vst.tsv", first.output_dir / "pca_scores.tsv", first.output_dir / "sample_correlation.tsv")}


def test_fastq_l1_requires_explicit_successful_handoff_and_upgrades_legacy(project_factory):
    root = project_factory()
    config = deepcopy(base_config())
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None, "strandedness": "auto", "quantification": {"method": "salmon"}}
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    fastq = root / "input" / "fastq"
    fastq.mkdir()
    for sample in ("C1", "C2", "C3", "T1", "T2", "T3"):
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    report = validate_project(root)
    with pytest.raises(DownstreamExecutionError, match="explicit --run-id"):
        prepare_l1(report, run_id=None)
    run = root / "runs" / "run-fixture"
    for directory in ("frozen", "handoff", "upstream/nfcore_rnaseq/salmon"):
        (run / directory).mkdir(parents=True, exist_ok=True)
    (run / "run_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    shutil.copy2(root / "metadata.csv", run / "frozen" / "metadata.csv")
    shutil.copy2(root / "project.yaml", run / "frozen" / "project.yaml")
    (run / "handoff" / "upstream_manifest.yaml").write_text(yaml.safe_dump({"samples": ["C1", "C2", "C3", "T1", "T2", "T3"]}, sort_keys=False), encoding="utf-8")
    (run / "upstream/nfcore_rnaseq/salmon/salmon.merged.tx2gene.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
    for sample in ("C1", "C2", "C3", "T1", "T2", "T3"):
        item = run / "upstream/nfcore_rnaseq/salmon" / sample
        item.mkdir()
        (item / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t10\n", encoding="utf-8")
    prepared = prepare_l1(report, run_id="run-fixture")
    assert prepared.source_type == "salmon_tximport"
    assert prepared.metadata_path == run / "frozen" / "metadata.csv"
    handoff = yaml.safe_load((run / "handoff" / "upstream_manifest.yaml").read_text())
    assert sorted(handoff["salmon"]["quant_sf"]) == ["C1", "C2", "C3", "T1", "T2", "T3"]
    assert handoff["salmon"]["tx2gene"]["mapping_type"] == "historical_ordinary"
    assert prepared.config["tx2gene_mapping"]["mapping_type"] == "historical_ordinary"


def test_raw_counts_rejects_upstream_run_id(project_factory):
    with pytest.raises(DownstreamExecutionError, match="only valid for FASTQ"):
        prepare_l1(validate_project(project_factory()), run_id="run-1")
