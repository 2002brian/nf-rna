from __future__ import annotations

import hashlib
import json
import shutil
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from conftest import BASE_CONTRASTS, base_config, pdf_text_origins, require_r_packages
from rnaseq.downstream import execute_l1, prepare_l1
from rnaseq.errors import DownstreamExecutionError
from rnaseq.validators import validate_project


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_raw_counts_l1_real_backend_is_structured_and_text_deterministic(project_factory):
    pytest.importorskip("yaml")
    require_r_packages("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
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
    assert all((first.output_dir / name).is_file() for name in ("vst.tsv", "pca.png", "pca.tiff", "sample_correlation.png", "sample_correlation.tiff", "boxplot_log2_raw_counts.png", "boxplot_vst.tiff", "sample_distribution_qc.tsv", "filtering_summary.yaml", "l1_report.md"))
    assert not list(first.output_dir.glob("*.svg"))
    assert {path.name for path in first.output_dir.glob("*.pdf")} == {"pca.pdf", "sample_correlation.pdf", "boxplot_log2_raw_counts.pdf", "boxplot_vst.pdf"}
    figures = json.loads((first.output_dir / "figure_manifest.json").read_text())
    assert {item["id"] for item in figures["figures"]} == {"pca", "sample_correlation", "boxplot_log2_raw_counts", "boxplot_vst"}
    box_stats = [line.split("\t") for line in (first.output_dir / "sample_distribution_qc.tsv").read_text().splitlines()]
    assert box_stats[0][:3] == ["sample_id", "matrix", "genes"]
    assert [row[0] for row in box_stats[1:7]] == ["C1", "C2", "C3", "T1", "T2", "T3"] and [row[0] for row in box_stats[7:]] == ["C1", "C2", "C3", "T1", "T2", "T3"]
    assert {row[1] for row in box_stats[1:]} == {"log2_raw_count_plus_1", "blind_vst"}
    vst_box = next(item for item in figures["figures"] if item["id"] == "boxplot_vst")
    assert vst_box["encoding"]["batch_marker"].startswith("none") and "neither show nor exclude a batch effect" in vst_box["interpretation"]
    pca = next(item for item in figures["figures"] if item["id"] == "pca")
    assert pca["labels"]["text"] == "every sample_id, unabbreviated" and "geom_text_repel, seed" in pca["labels"]["placement"]
    assert pca["axes"] == "PC1/PC2 with proportion of variance explained; equal axis scaling"  # PC2 varies: no constant-axis handling
    assert figures["package_versions"]["ggrepel"]
    width_pt = next(item for item in pca["files"] if item["format"] == "pdf")["width_pt"]
    assert all(0 <= x <= width_pt for x, _ in pdf_text_origins(first.output_dir / "pca.pdf"))
    before = {_path.name: _digest(_path) for _path in (normalized, first.output_dir / "vst.tsv", first.output_dir / "pca_scores.tsv", first.output_dir / "sample_correlation.tsv", first.output_dir / "pca.tiff", first.output_dir / "pca.png")}
    execute_l1(prepared)
    # Repelled PCA labels use a fixed seed and no time limit: rasters are byte-identical.
    assert before == {_path.name: _digest(_path) for _path in (normalized, first.output_dir / "vst.tsv", first.output_dir / "pca_scores.tsv", first.output_dir / "sample_correlation.tsv", first.output_dir / "pca.tiff", first.output_dir / "pca.png")}


def test_l1_boxplots_annotate_batch_when_declared_categorical(project_factory):
    require_r_packages("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    counts = "gene_id,C1,C2,C3,T1,T2,T3\n" + "\n".join(
        f"Gene{index},{10 + index},{12 + index},{9 + index},{40 + index},{45 + index},{43 + index}" for index in range(1, 101)
    ) + "\n"
    metadata = "sample_id,batch,condition\nC1,B1,Control\nC2,B1,Control\nC3,B2,Control\nT1,B1,Treatment\nT2,B2,Treatment\nT3,B2,Treatment\n"
    config = deepcopy(base_config())
    config["schema_version"] = "1.3"
    config["analysis"] = {"enrichment": []}
    config["design"] = {"type": "two_group", "formula": "~ batch + condition", "variables": {"batch": "categorical", "condition": "categorical"}}
    result = execute_l1(prepare_l1(validate_project(project_factory(config=config, counts=counts, metadata=metadata)), run_id=None))
    figures = json.loads((result.output_dir / "figure_manifest.json").read_text())
    for figure_id in ("boxplot_log2_raw_counts", "boxplot_vst"):
        entry = next(item for item in figures["figures"] if item["id"] == figure_id)
        assert entry["encoding"]["fill"] == "condition"
        assert entry["encoding"]["batch_marker"] == "black marker shape below each box"
        assert all((result.output_dir / f"{figure_id}.{extension}").is_file() for extension in ("png", "tiff", "pdf"))
    assert (result.output_dir / "pca_by_batch.png").is_file()


@pytest.mark.parametrize("batch_values, expected_marker, in_design", (
    (("B1", "B1", "B2", "B1", "B2", "B2"), "black marker shape below each box", False),
    (("B1", "", "B2", "B1", "B2", "B2"), "none (batch column has missing values for some samples)", False),
))
def test_l1_boxplot_batch_annotation_is_independent_of_the_de_design(project_factory, batch_values, expected_marker, in_design):
    require_r_packages("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    counts = "gene_id,C1,C2,C3,T1,T2,T3\n" + "\n".join(
        f"Gene{index},{10 + index},{12 + index},{9 + index},{40 + index},{45 + index},{43 + index}" for index in range(1, 101)
    ) + "\n"
    samples, conditions = ("C1", "C2", "C3", "T1", "T2", "T3"), ("Control",) * 3 + ("Treatment",) * 3
    metadata = "sample_id,condition,batch\n" + "".join(f"{s},{c},{b}\n" for s, c, b in zip(samples, conditions, batch_values))
    # base_config() uses "~ condition": batch is in the metadata but not in the DE design.
    result = execute_l1(prepare_l1(validate_project(project_factory(counts=counts, metadata=metadata)), run_id=None))
    assert json.loads((result.output_dir / "backend_config.json").read_text())["formula"] == "~ condition"  # QC annotation never touches the model
    figures = json.loads((result.output_dir / "figure_manifest.json").read_text())
    for figure_id in ("boxplot_log2_raw_counts", "boxplot_vst"):
        entry = next(item for item in figures["figures"] if item["id"] == figure_id)
        assert entry["encoding"]["batch_marker"] == expected_marker
        assert entry["encoding"]["batch_in_de_design"] is in_design
        assert entry["encoding"]["batch_values"] == "plotted values are not batch-corrected"
    # PCA-by-batch remains tied to the design; the boxplot annotation does not create it.
    assert not (result.output_dir / "pca_by_batch.png").exists()
    # The plotted matrices are the same with or without a batch column.
    plain = execute_l1(prepare_l1(validate_project(project_factory(counts=counts)), run_id=None))
    assert _digest(result.output_dir / "vst.tsv") == _digest(plain.output_dir / "vst.tsv")
    assert _digest(result.output_dir / "sample_distribution_qc.tsv") == _digest(plain.output_dir / "sample_distribution_qc.tsv")
    plain_figures = json.loads((plain.output_dir / "figure_manifest.json").read_text())
    assert next(i for i in plain_figures["figures"] if i["id"] == "boxplot_vst")["encoding"]["batch_marker"] == "none (no batch column in metadata)"


@pytest.mark.parametrize("figure_settings", [None, {"preset": "standard", "profile": "general"}], ids=["legacy", "standard"])
def test_l1_runs_for_two_samples_where_pc2_does_not_exist(project_factory, figure_settings):
    require_r_packages("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "ggrepel")
    counts = "gene_id,C1,T1\n" + "\n".join(f"Gene{index},{10 + index},{40 + index}" for index in range(1, 101)) + "\n"
    metadata = "sample_id,condition\nC1,Control\nT1,Treatment\n"
    config = base_config()
    if figure_settings is not None:
        config["figures"] = figure_settings
    result = execute_l1(prepare_l1(validate_project(project_factory(config=config, counts=counts, metadata=metadata)), run_id=None))
    assert all((result.output_dir / name).is_file() for name in ("pca.png", "boxplot_log2_raw_counts.png", "boxplot_vst.png"))
    # The scores and variance are written unchanged: PC1 carries everything, PC2 is constant at 0.
    scores = [line.split("\t") for line in (result.output_dir / "pca_scores.tsv").read_text().splitlines()]
    assert [row[0] for row in scores[1:]] == ["C1", "T1"]
    assert float(scores[1][1]) == pytest.approx(-float(scores[2][1])) and float(scores[1][1]) != 0
    assert all(abs(float(row[2])) < 1e-8 for row in scores[1:])
    variance = dict(line.split("\t") for line in (result.output_dir / "pca_variance.tsv").read_text().splitlines()[1:])
    assert float(variance["PC1"]) == pytest.approx(1) and float(variance["PC2"]) == pytest.approx(0, abs=1e-12)
    figures = json.loads((result.output_dir / "figure_manifest.json").read_text())
    pca = next(item for item in figures["figures"] if item["id"] == "pca")
    assert "PC2 is constant (no variance)" in pca["axes"] and pca["axes"].endswith("single break at 0")
    assert pca["labels"]["text"] == "every sample_id, unabbreviated"
    assert {item["format"] for item in pca["files"]} == {"png", "tiff", "pdf"}
    assert all((result.output_dir / f"pca.{extension}").is_file() for extension in ("png", "tiff", "pdf"))
    # The panel is sized for a readable PC2 extent, not the former 2% sliver at the height floor
    # (50 mm = 142 pt legacy, 40 mm = 113 pt standard).
    pdf = next(item for item in pca["files"] if item["format"] == "pdf")
    assert pdf["height_pt"] > (200 if figure_settings is None else 150)
    assert all(0 <= x <= pdf["width_pt"] and 0 <= y <= pdf["height_pt"] for x, y in pdf_text_origins(result.output_dir / "pca.pdf"))


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
    # The handoff happens to be lexical; frozen metadata order deliberately is
    # not, and remains the order passed to tximport/DESeq2.
    (run / "frozen" / "metadata.csv").write_text(
        "sample_id,condition\nT1,Treatment\nC1,Control\nT2,Treatment\nC2,Control\nT3,Treatment\nC3,Control\n",
        encoding="utf-8",
    )
    prepared = prepare_l1(report, run_id="run-fixture")
    assert prepared.source_type == "salmon_tximport"
    assert prepared.metadata_path == run / "frozen" / "metadata.csv"
    assert prepared.sample_ids == ("T1", "C1", "T2", "C2", "T3", "C3")
    assert list(prepared.config["quant_sf"]) == list(prepared.sample_ids)
    handoff = yaml.safe_load((run / "handoff" / "upstream_manifest.yaml").read_text())
    assert sorted(handoff["salmon"]["quant_sf"]) == ["C1", "C2", "C3", "T1", "T2", "T3"]
    assert handoff["salmon"]["tx2gene"]["mapping_type"] == "historical_ordinary"
    assert prepared.config["tx2gene_mapping"]["mapping_type"] == "historical_ordinary"


def test_raw_counts_rejects_upstream_run_id(project_factory):
    with pytest.raises(DownstreamExecutionError, match="only valid for FASTQ"):
        prepare_l1(validate_project(project_factory()), run_id="run-1")
