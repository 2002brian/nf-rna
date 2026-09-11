from __future__ import annotations

import json
import hashlib
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import pytest
import yaml

from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.execution import RuntimeCheck, load_run_states
from rnaseq.planner import generate_plan
from rnaseq.service import (
    FrozenInputs,
    _sanitize_delivery_appledouble,
    _assert_delivery_appledouble_free,
    assemble_delivery,
    build_downstream_nextflow_command,
    create_case_run,
    delivery_filename,
    execute_service_run,
    finalize_fastq_handoff,
    freeze_case_inputs,
    prepare_service_run,
    resolve_downstream_inputs,
    reuse_upstream_if_compatible,
    sanitize_completed_delivery,
    taipei_run_timestamp,
    validate_case_id,
    write_downstream_docker_user_config,
    write_downstream_observer_config,
)
from rnaseq.validators import validate_project


def _frozen_run(project_factory) -> tuple[Path, object, object]:
    root = project_factory()
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-20260828-001", moment=datetime(2026, 8, 28, 11, 51, 11))
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    return root, report, run


def _fastq_handoff_run(tmp_path: Path):
    """Create only the documented upstream handoff artifacts; no workflow runs."""

    root = tmp_path / "fastq-project"
    (root / "input" / "fastq").mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (root / "input" / "fastq" / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text(
        "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8"
    )
    (root / "contrasts.csv").write_text(
        "contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n",
        encoding="utf-8",
    )
    from conftest import base_config

    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None,
        "strandedness": "auto", "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    run = create_case_run(report, "CASE-STAGING", moment=datetime(2026, 8, 28, 12, 2, 0))
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    salmon = run.run_dir / "upstream" / "nfcore_rnaseq" / "salmon"
    salmon.mkdir(parents=True)
    (salmon / "salmon.merged.gene_counts.tsv").write_text(
        "gene_id\tC1\tC2\tT1\tT2\nGeneA\t1\t2\t3\t4\n", encoding="utf-8"
    )
    (salmon / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
    for sample in ("C1", "C2", "T1", "T2"):
        quant = salmon / sample / "quant.sf"
        quant.parent.mkdir()
        quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
    multiqc = run.run_dir / "upstream" / "nfcore_rnaseq" / "multiqc"
    (multiqc / "multiqc_data").mkdir(parents=True)
    (multiqc / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
    finalize_fastq_handoff(report, run, frozen.contract)
    return run, salmon


def test_case_id_validation_rejects_unsafe_values():
    assert validate_case_id("CASE-20260828-001") == "CASE-20260828-001"
    for value in ("", "..", "case/name", "case\\name", "case\nname", "a..b"):
        with pytest.raises(ExecutionPreflightError):
            validate_case_id(value)


def test_taipei_timestamp_formatting():
    value = taipei_run_timestamp(datetime(2026, 8, 28, 11, 51, 11))
    assert value == "20260828-115111+0800"
    assert delivery_filename("report.html", value) == "report_20260828.html"
    assert delivery_filename("report.html", "20260828-235959+0800-01") == "report_20260828.html"
    with pytest.raises(ExecutionPreflightError, match="unsafe run ID"):
        delivery_filename("report.html", "20260828")


def test_case_runs_are_immutable_and_freeze_raw_count_inputs(project_factory):
    root, _report, first = _frozen_run(project_factory)
    second = create_case_run(validate_project(root), "CASE-20260828-001", moment=datetime(2026, 8, 28, 11, 51, 11))
    assert first.run_dir != second.run_dir
    assert first.run_id == "20260828-115111+0800"
    assert second.run_id == "20260828-115111+0800-01"
    frozen = first.run_dir / "frozen"
    before = (frozen / "metadata.csv").read_text(encoding="utf-8")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,Changed\n", encoding="utf-8")
    assert (frozen / "metadata.csv").read_text(encoding="utf-8") == before
    contract = json.loads((frozen / "downstream_contract.json").read_text(encoding="utf-8"))
    assert Path(contract["source"]["counts"]).is_file()
    assert Path(contract["source"]["counts"]).parent == frozen / "input"


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("missing_tx2gene", "salmon.tx2gene.path is unavailable"),
        ("missing_quant", "salmon.quant_sf.C1 is unavailable"),
        ("appledouble_tx2gene", "selected an AppleDouble artifact"),
        ("outside_upstream_tx2gene", "outside the immutable upstream output"),
    ),
)
def test_fastq_execution_staging_fails_before_downstream_for_missing_or_appledouble_handoff_artifacts(tmp_path, mutation, message):
    run, salmon = _fastq_handoff_run(tmp_path)
    if mutation == "missing_tx2gene":
        (salmon / "salmon.merged.tx2gene_augmented.tsv").unlink()
    elif mutation == "missing_quant":
        (salmon / "C1" / "quant.sf").unlink()
    else:
        handoff_path = run.run_dir / "frozen" / "upstream_handoff_manifest.yaml"
        handoff = yaml.safe_load(handoff_path.read_text(encoding="utf-8"))
        if mutation == "outside_upstream_tx2gene":
            handoff["salmon"]["tx2gene"]["path"] = "frozen/metadata.csv"
            handoff_path.write_text(yaml.safe_dump(handoff, sort_keys=False), encoding="utf-8")
            with pytest.raises(UpstreamExecutionError, match=message):
                resolve_downstream_inputs(run)
            assert not (run.run_dir / "downstream_inputs").exists()
            return
        sidecar = salmon / "._salmon.merged.tx2gene_augmented.tsv"
        sidecar.write_text("not a biological input\n", encoding="utf-8")
        handoff["salmon"]["tx2gene"]["path"] = str(sidecar.relative_to(run.run_dir))
        handoff_path.write_text(yaml.safe_dump(handoff, sort_keys=False), encoding="utf-8")

    with pytest.raises(UpstreamExecutionError, match=message):
        resolve_downstream_inputs(run)
    assert not (run.run_dir / "downstream_inputs").exists()


def test_downstream_command_is_argument_array_and_delivery_is_allowlisted(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    execution_inputs = resolve_downstream_inputs(run)
    observer_config = write_downstream_observer_config(run)
    command = build_downstream_nextflow_command(run, observer_config=observer_config, execution_inputs=execution_inputs)
    assert command[:2] == ["nextflow", "run"]
    assert "--contract" in command and "--inputs" in command and "--outdir" in command
    assert command[command.index("--analysis_level") + 1] == "L2"
    assert command[command.index("--inputs") + 1] == str(execution_inputs.root.resolve())
    assert "/Volumes/KOXIA" not in command
    assert not any(argument.startswith("-with-") for argument in command)
    assert "--enrichment" not in command
    observer_text = observer_config.read_text(encoding="utf-8")
    for observer in ("trace", "report", "timeline", "dag"):
        assert f"{observer} {{" in observer_text
    assert observer_text.count("overwrite = true") == 4
    downstream = run.run_dir / "downstream"
    (downstream / "report").mkdir()
    (downstream / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
    (downstream / "l2" / "contrasts" / "a").mkdir(parents=True)
    (downstream / "l2" / "contrasts" / "a" / "all_genes.tsv").write_text("gene_id\nGeneA\n", encoding="utf-8")
    (downstream / "l2" / "contrasts" / "a" / "volcano.png").write_bytes(b"png")
    (downstream / "l2" / "contrasts" / "a" / "volcano.tiff").write_bytes(b"tiff")
    (downstream / "l2" / "contrasts" / "a" / "volcano.pdf").write_bytes(b"pdf")
    (downstream / "l2" / "contrasts" / "a" / "volcano.svg").write_bytes(b"svg")
    (downstream / "l2" / "logs").mkdir()
    (downstream / "l2" / "logs" / "r.stderr.log").write_text("internal", encoding="utf-8")
    (downstream / "l2" / "enrichment" / "gsea_go" / "BP").mkdir(parents=True)
    (downstream / "l2" / "enrichment" / "gsea_go" / "BP" / "all_terms.tsv").write_text("ID\tDescription\nGO:1\tterm\n", encoding="utf-8")
    (downstream / "l2" / "enrichment" / "gsea_go" / "BP" / "dotplot.png").write_bytes(b"png")
    (downstream / "l2" / "enrichment" / "gsea_go" / "BP" / "dotplot.tiff").write_bytes(b"tiff")
    (downstream / "l2" / "enrichment" / "gsea_kegg" / "a").mkdir(parents=True)
    (downstream / "l2" / "enrichment" / "gsea_kegg" / "a" / "all_terms.tsv").write_text("ID\tDescription\npath:mmu1\tpathway\n", encoding="utf-8")
    (downstream / "l2" / "enrichment" / "gsea_kegg" / "a" / "._ignored.tsv").write_bytes(b"ignored")
    (downstream / "._ignored.png").write_bytes(b"ignored")
    (downstream / "l2" / "contrasts" / "a" / "._ignored.tsv").write_bytes(b"ignored")
    multiqc = run.run_dir / "upstream" / "nfcore_rnaseq" / "multiqc" / "multiqc_report.html"
    multiqc.parent.mkdir(parents=True)
    multiqc.write_text("<html></html>", encoding="utf-8")
    (run.run_dir / "frozen" / "upstream_handoff_manifest.yaml").write_text(
        yaml.safe_dump({"multiqc": {"html": str(multiqc.relative_to(run.run_dir))}}, sort_keys=False),
        encoding="utf-8",
    )
    (run.run_dir / "provenance" / "run_provenance.yaml").write_text("pipeline_version: 0.4.3\n", encoding="utf-8")
    delivery = assemble_delivery(run)
    assert (delivery / "report_20260828.html").is_file()
    assert (delivery / "README_20260828.md").is_file()
    assert not (delivery / "report.html").exists()
    assert not (delivery / "README.md").exists()
    assert (delivery / "methods_and_versions" / "execution_manifest_20260828.yaml").is_file()
    assert (delivery / "methods_and_versions" / "input_manifest_20260828.yaml").is_file()
    assert (delivery / "methods_and_versions" / "upstream_handoff_manifest_20260828.yaml").is_file()
    assert (delivery / "methods_and_versions" / "run_state_20260828.json").is_file()
    assert (delivery / "methods_and_versions" / "run_provenance_20260828.yaml").is_file()
    assert (delivery / "multiqc" / "multiqc_report_20260828.html").is_file()
    assert not any("115111" in path.name or "+0800" in path.name for path in delivery.rglob("*"))
    assert load_run_states(run.run_dir.parents[2])[0]["delivery_available"] is True
    assert (delivery / "tables" / "l2" / "contrasts" / "a" / "all_genes.tsv").is_file()
    assert (delivery / "figures" / "png" / "l2" / "contrasts" / "a" / "volcano.png").is_file()
    assert (delivery / "figures" / "tiff_300dpi" / "l2" / "contrasts" / "a" / "volcano.tiff").is_file()
    assert (delivery / "tables" / "l2" / "enrichment" / "gsea_go" / "BP" / "all_terms.tsv").is_file()
    assert (delivery / "tables" / "l2" / "enrichment" / "gsea_kegg" / "a" / "all_terms.tsv").is_file()
    assert (delivery / "figures" / "png" / "l2" / "enrichment" / "gsea_go" / "BP" / "dotplot.png").is_file()
    assert (delivery / "figures" / "tiff_300dpi" / "l2" / "enrichment" / "gsea_go" / "BP" / "dotplot.tiff").is_file()
    assert not (downstream / "enrichment").exists()
    assert not (downstream / "l2" / "enrichment" / "enrichment").exists()
    assert not (delivery / "figures" / "vector_pdf_svg").exists()
    assert not list(delivery.rglob("*.pdf"))
    assert not list(delivery.rglob("*.svg"))
    assert not any(path.name.endswith(".log") or path.name.startswith("._") for path in delivery.rglob("*"))


def test_imported_raw_counts_delivery_is_exact_and_has_integer_semantics(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    source = run.run_dir / "frozen" / "input" / "counts.csv"
    delivery = assemble_delivery(run)
    target = delivery / "counts" / "raw_counts.csv"
    assert target.read_bytes() == source.read_bytes()
    manifest = json.loads((delivery / "counts" / "artifact_manifest.json").read_text(encoding="utf-8"))
    artifact = manifest["artifacts"][0]
    assert artifact["source_type"] == "raw_counts"
    assert artifact["value_semantics"] == "integer_raw_counts"
    assert artifact["integer_required"] is True
    assert artifact["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest()


def test_featurecounts_delivery_uses_canonical_integer_matrix(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    matrix = run.run_dir / "upstream" / "hisat2_featurecounts" / "counts" / "canonical_counts.csv"
    matrix.parent.mkdir(parents=True)
    matrix.write_text("gene_id,C1,C2,C3,T1,T2,T3\nGeneA,1,2,3,4,5,6\n", encoding="utf-8")
    handoff = {"featurecounts": {"canonical_matrix": str(matrix.relative_to(run.run_dir))}}
    (run.run_dir / "frozen" / "upstream_handoff_manifest.yaml").write_text(yaml.safe_dump(handoff), encoding="utf-8")
    contract_path = run.run_dir / "frozen" / "downstream_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["source"] = {"type": "featurecounts_raw_counts", "construction_method": "DESeqDataSetFromMatrix", "upstream_handoff": str((run.run_dir / "frozen" / "upstream_handoff_manifest.yaml").resolve())}
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    delivery = assemble_delivery(run)
    assert (delivery / "counts" / "raw_counts.csv").read_bytes() == matrix.read_bytes()
    artifact = json.loads((delivery / "counts" / "artifact_manifest.json").read_text())['artifacts'][0]
    assert artifact["source_type"] == "featurecounts_raw_counts"
    assert artifact["deseq2_construction_method"] == "DESeqDataSetFromMatrix"


def test_salmon_delivery_is_estimated_never_raw_and_qc_uses_upstream_matrix(tmp_path):
    run, _salmon = _fastq_handoff_run(tmp_path)
    delivery = assemble_delivery(run)
    estimated = delivery / "counts" / "estimated_counts.csv"
    assert estimated.read_text(encoding="utf-8") == "gene_id,C1,C2,T1,T2\nGeneA,1,2,3,4\n"
    assert not (delivery / "counts" / "raw_counts.csv").exists()
    artifact = json.loads((delivery / "counts" / "artifact_manifest.json").read_text())['artifacts'][0]
    assert artifact["source_type"] == "salmon_tximport"
    assert artifact["value_semantics"] == "salmon_estimated_counts"
    assert artifact["integer_required"] is False


def test_salmon_l1_delivery_preserves_noninteger_tximport_source_counts(tmp_path):
    run, _salmon = _fastq_handoff_run(tmp_path)
    l1 = run.run_dir / "downstream" / "l1"
    l1.mkdir(parents=True)
    source_counts = l1 / "source_counts.csv"
    source_counts.write_text("gene_id,C1,C2,T1,T2\nGeneA,1.25,2.5,3.75,4.125\n", encoding="utf-8")
    delivery = assemble_delivery(run)
    target = delivery / "counts" / "estimated_counts.csv"
    assert target.read_bytes() == source_counts.read_bytes()
    assert not (delivery / "counts" / "raw_counts.csv").exists()
    artifact = json.loads((delivery / "counts" / "artifact_manifest.json").read_text())['artifacts'][0]
    assert artifact["deseq2_construction_method"] == "DESeqDataSetFromTximport"
    assert artifact["integer_required"] is False


def test_vst_delivery_is_explicitly_transformed_and_accepts_negative_values(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    vst = run.run_dir / "downstream" / "l1" / "vst.csv"
    vst.parent.mkdir(parents=True)
    vst.write_text("gene_id,C1,C2,C3,T1,T2,T3\nGeneA,-1.2,0,1.1,2.2,3.3,4.4\n", encoding="utf-8")
    delivery = assemble_delivery(run)
    assert (delivery / "counts" / "vst.csv").read_bytes() == vst.read_bytes()
    artifacts = json.loads((delivery / "counts" / "artifact_manifest.json").read_text())['artifacts']
    vst_artifact = next(item for item in artifacts if item["filename"] == "counts/vst.csv")
    assert vst_artifact["value_semantics"] == "variance_stabilized_expression"
    assert vst_artifact["normalized"] is True


def test_delivery_finalizer_removes_only_appledouble_without_following_symlinks(tmp_path):
    delivery = tmp_path / "run" / "delivery"
    delivery.mkdir(parents=True)
    (delivery / "report.html").write_text("report", encoding="utf-8")
    (delivery / ".keep").write_text("intended hidden file", encoding="utf-8")
    (delivery / "._after_copy.tsv").write_text("sidecar", encoding="utf-8")
    nested_sidecar = delivery / "figures" / "png" / "l1" / "._pca.png"
    nested_sidecar.parent.mkdir(parents=True)
    nested_sidecar.write_text("sidecar", encoding="utf-8")
    directory_sidecar = delivery / "figures" / "._pca"
    directory_sidecar.mkdir(parents=True)
    (directory_sidecar / "metadata").write_text("sidecar", encoding="utf-8")
    external = tmp_path / "outside-delivery.txt"
    external.write_text("must remain", encoding="utf-8")
    (delivery / "._outside").symlink_to(external)
    previous_delivery = tmp_path / "previous-run" / "delivery"
    previous_delivery.mkdir(parents=True)
    (previous_delivery / "._prior.tsv").write_text("prior sidecar", encoding="utf-8")

    _sanitize_delivery_appledouble(delivery)
    _sanitize_delivery_appledouble(delivery)

    assert (delivery / "report.html").read_text(encoding="utf-8") == "report"
    assert (delivery / ".keep").read_text(encoding="utf-8") == "intended hidden file"
    assert not nested_sidecar.exists()
    assert not list(delivery.rglob("._*"))
    assert external.read_text(encoding="utf-8") == "must remain"
    assert (previous_delivery / "._prior.tsv").is_file()

    escaped_root = tmp_path / "external-delivery"
    escaped_root.mkdir()
    (escaped_root / "._must-not-touch").write_text("external", encoding="utf-8")
    redirected_delivery = tmp_path / "redirected-delivery"
    redirected_delivery.symlink_to(escaped_root, target_is_directory=True)
    with pytest.raises(OSError, match="must not be a symlink"):
        _sanitize_delivery_appledouble(redirected_delivery)
    assert (escaped_root / "._must-not-touch").is_file()


def test_delivery_hard_acceptance_check_reports_nested_remaining_appledouble(tmp_path):
    delivery = tmp_path / "run" / "delivery"
    sidecar = delivery / "figures" / "png" / "l2" / "._volcano.png"
    sidecar.parent.mkdir(parents=True)
    sidecar.write_bytes(b"sidecar")

    with pytest.raises(OSError, match="Delivery finalization failed") as error:
        _assert_delivery_appledouble_free(delivery)
    assert "figures/png/l2/._volcano.png" in str(error.value)


def test_sanitize_completed_delivery_touches_only_a_successful_run_delivery(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    delivery = run.run_dir / "delivery"
    nested_sidecar = delivery / "figures" / "png" / "l1" / "._pca.png"
    nested_sidecar.parent.mkdir(parents=True)
    nested_sidecar.write_bytes(b"sidecar")
    (delivery / "report.html").write_text("report", encoding="utf-8")

    run.state_path.write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    assert sanitize_completed_delivery(run.run_dir) == delivery.resolve()
    assert not nested_sidecar.exists()
    assert (delivery / "report.html").read_text(encoding="utf-8") == "report"

    run.state_path.write_text(json.dumps({"status": "FAILED"}), encoding="utf-8")
    with pytest.raises(ExecutionPreflightError, match="only for a successful"):
        sanitize_completed_delivery(run.run_dir)


def test_delivery_sidecar_after_sanitization_prevents_success_state(monkeypatch, project_factory, tmp_path):
    root = project_factory()
    report = validate_project(root)
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "local-nextflow-cache"))
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))

    def fake_nextflow(command, *, cwd, stdout_path, stderr_path):
        outdir = Path(command[command.index("--outdir") + 1])
        (outdir / "report").mkdir(parents=True)
        (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        stdout_path.write_text("mock stdout\n", encoding="utf-8")
        stderr_path.write_text("mock stderr\n", encoding="utf-8")
        return 0

    original_sanitizer = _sanitize_delivery_appledouble

    def delayed_sidecar(delivery):
        original_sanitizer(delivery)
        sidecar = delivery / "figures" / "png" / "l1" / "._late.png"
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        sidecar.write_bytes(b"late sidecar")

    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow)
    monkeypatch.setattr("rnaseq.service._sanitize_delivery_appledouble", delayed_sidecar)
    with pytest.raises(OSError, match="AppleDouble entries remain"):
        execute_service_run(report, case_id="CASE-DELIVERY-SIDECAR")

    state_path = next((root / "runs" / "CASE-DELIVERY-SIDECAR").glob("*/run_state.json"))
    assert json.loads(state_path.read_text(encoding="utf-8"))["status"] == "FAILED"


def test_public_gsea_selection_expands_to_only_the_two_internal_backends(project_factory):
    from conftest import base_config

    config = base_config()
    config["schema_version"] = "1.1"
    config["annotation"] = {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}
    config["analysis"] = {"enrichment": "gsea"}
    root = project_factory(config=config)
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    run = create_case_run(report, "CASE-GSEA", moment=datetime(2026, 8, 28, 12, 1, 0))
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    contract = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    assert contract["analysis_level"] == "L2"
    assert contract["analysis"] == {"enrichment": ["gsea"]}
    command = build_downstream_nextflow_command(run)
    assert command[command.index("--enrichment") + 1] == "gsea-go,gsea-kegg"

    legacy = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    legacy["analysis"] = {"enrichment": ["go", "gsea-go", "kegg", "gsea-kegg"]}
    (run.run_dir / "frozen" / "downstream_contract.json").write_text(
        json.dumps(legacy, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    migrated_command = build_downstream_nextflow_command(run)
    assert migrated_command[migrated_command.index("--enrichment") + 1] == "gsea-go,gsea-kegg"


def test_downstream_command_rejects_an_invalid_frozen_enrichment_module(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    contract_path = run.run_dir / "frozen" / "downstream_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["analysis"] = {"enrichment": ["not-a-production-backend"]}
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="analysis.enrichment"):
        build_downstream_nextflow_command(run)


def test_downstream_command_freezes_l1_and_rejects_enrichment_or_invalid_levels(project_factory):
    from conftest import base_config

    config = base_config()
    config["project"]["preset"] = "L1"
    root = project_factory(config=config)
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    run = create_case_run(report, "CASE-L1", moment=datetime(2026, 8, 28, 12, 2, 0))
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    contract_path = run.run_dir / "frozen" / "downstream_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    assert contract["analysis_level"] == "L1"
    command = build_downstream_nextflow_command(run)
    assert command[command.index("--analysis_level") + 1] == "L1"
    assert "--enrichment" not in command

    contract["analysis"] = {"enrichment": ["gsea"]}
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="analysis_level L1"):
        build_downstream_nextflow_command(run)

    contract["analysis"] = {"enrichment": []}
    contract["analysis_level"] = "L3"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="analysis_level must be L1 or L2"):
        build_downstream_nextflow_command(run)


def test_service_runs_nextflow_from_local_execution_root_and_preserves_case_outputs(monkeypatch, tmp_path):
    """A fake Nextflow process proves cache state never lands under the project."""

    root = tmp_path / "mounted-client-project"
    (root / "input" / "fastq").mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (root / "input" / "fastq" / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    from conftest import base_config

    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None, "strandedness": "auto", "quantification": {"method": "salmon"}}
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    local_root = tmp_path / "local-nextflow-cache"
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(local_root))
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "Docker daemon is available."))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))
    observed: list[tuple[list[str], Path]] = []

    def fake_nextflow(command, *, cwd, stdout_path, stderr_path):
        observed.append((command, cwd))
        (cwd / ".nextflow" / "cache").mkdir(parents=True, exist_ok=True)
        (cwd / ".nextflow" / "cache" / "000003.log").write_text("cache", encoding="utf-8")
        (cwd / ".nextflow" / "cache" / "._000003.log").write_text("appledouble", encoding="utf-8")
        stdout_path.write_text("mock stdout\n", encoding="utf-8")
        stderr_path.write_text("mock stderr\n", encoding="utf-8")
        if "nf-core/rnaseq" in command:
            outdir = Path(command[command.index("--outdir") + 1])
            (outdir / "salmon").mkdir(parents=True)
            (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\tC2\tT1\tT2\nGeneA\t1.0\t2.0\t3.0\t4.0\n", encoding="utf-8")
            (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
            for sample in ("C1", "C2", "T1", "T2"):
                sample_dir = outdir / "salmon" / sample
                sample_dir.mkdir()
                (sample_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
            (outdir / "multiqc" / "multiqc_data").mkdir(parents=True)
            (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
        else:
            outdir = Path(command[command.index("--outdir") + 1])
            (outdir / "report").mkdir(parents=True)
            (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return 0

    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow)
    run = execute_service_run(report, case_id="CASE-20260828-001")
    assert len(observed) == 2
    assert all(cwd == local_root / run.case_id / run.run_id / "launch" for _command, cwd in observed)
    assert all("-work-dir" in command for command, _cwd in observed)
    local_resource_config = run.run_dir / "frozen" / "nfcore.local.config"
    assert all(str(local_resource_config.resolve()) in command for command, _cwd in observed)
    assert (local_root / run.case_id / run.run_id / "launch" / ".nextflow" / "cache" / "000003.log").is_file()
    assert not (run.run_dir / ".nextflow").exists()
    assert not (run.run_dir / "work").exists()
    assert (run.run_dir / "upstream" / "nfcore_rnaseq" / "salmon" / "salmon.merged.gene_counts.tsv").is_file()
    assert (run.run_dir / "downstream" / "report" / "report.html").is_file()
    execution_inputs = run.run_dir / "downstream_inputs"
    staged = json.loads((execution_inputs / "execution_inputs.json").read_text(encoding="utf-8"))
    assert staged["source"]["type"] == "salmon_tximport"
    assert staged["samples"] == ["C1", "C2", "T1", "T2"]
    assert set(staged["source"]["quant_sf"]) == set(staged["samples"])
    assert staged["source"]["tx2gene"] == "source/salmon.merged.tx2gene_augmented.tsv"
    assert all(path.endswith("/quant.sf") for path in staged["source"]["quant_sf"].values())
    assert all((execution_inputs / path).is_file() for path in staged["source"]["quant_sf"].values())
    assert (execution_inputs / staged["source"]["tx2gene"]).is_file()
    assert "/Volumes/KOXIA" not in json.dumps(staged, sort_keys=True)
    frozen_contract = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    assert frozen_contract["source"]["upstream_handoff"] == str((run.run_dir / "frozen" / "upstream_handoff_manifest.yaml").resolve())
    downstream_command = observed[1][0]
    assert downstream_command[downstream_command.index("--inputs") + 1] == str(execution_inputs.resolve())
    runtime_config = run.run_dir / "frozen" / "downstream.runtime.config"
    assert runtime_config.read_text(encoding="utf-8") == 'process.container = "rnaseq-control-plane:latest"\n'
    assert str(runtime_config.resolve()) in downstream_command
    assert "/Volumes/KOXIA" not in downstream_command
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["execution_root"] == str(local_root / run.case_id / run.run_id)
    assert provenance["fastq_preprocessing"] == "raw"
    assert provenance["skip_trimming"] is False
    assert provenance["nfcore_preprocessing_arguments"] == []
    assert provenance["nfcore_runtime_params"] == {"skip_alignment": True}
    assert provenance["salmon_tx2gene"]["mapping_type"] == "nfcore_tx2gene_augmented"
    assert provenance["salmon_tx2gene"]["path"].endswith("salmon.merged.tx2gene_augmented.tsv")
    assert len(provenance["salmon_tx2gene"]["sha256"]) == 64
    assert provenance["container_image"]["reference"] == "rnaseq-control-plane:latest"
    assert provenance["production_intended"] is False
    assert set(provenance["workflow_sha256"]) == {"workflow/main.nf", "workflow/hisat2_featurecounts.nf"}
    assert all(len(value) == 64 for value in provenance["workflow_sha256"].values())
    assert provenance["runtime_resources"]["resource_profile"] == "M5_LOCAL_SMALL_MEDIUM_LARGE"
    assert "host_architecture" in provenance["runtime_resources"]
    assert provenance["runtime_resources"]["requested"] == {"cpus": 8, "memory_gib": 12}
    assert provenance["runtime_resources"]["effective"] == {"cpus": 8, "memory_gib": 12}
    assert provenance["frozen_local_nextflow_config"]["path"] == "frozen/nfcore.local.config"
    assert len(provenance["frozen_local_nextflow_config"]["sha256"]) == 64
    execution_manifest = yaml.safe_load((run.run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    assert execution_manifest["fastq_preprocessing"] == "raw"
    assert execution_manifest["skip_trimming"] is False
    assert execution_manifest["nfcore_preprocessing_arguments"] == []
    assert execution_manifest["nfcore_runtime_params"] == {"skip_alignment": True}
    runtime_params = json.loads((run.run_dir / "frozen" / "nfcore.params.json").read_text(encoding="utf-8"))
    assert runtime_params == {"skip_alignment": True}
    assert runtime_params["skip_alignment"] is True
    assert "skip_trimming" not in runtime_params
    assert execution_manifest["local_nextflow_config"] == provenance["frozen_local_nextflow_config"]
    observer_config = run.run_dir / "frozen" / "downstream.observers.config"
    assert observer_config.is_file()
    assert "overwrite = true" in observer_config.read_text(encoding="utf-8")
    assert json.loads(run.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"


def test_downstream_linux_docker_user_config_is_frozen_and_only_adds_a_nextflow_override(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    config = write_downstream_docker_user_config(run, host_os="Linux", uid=24701, gid=24703)
    assert config == run.run_dir / "frozen" / "downstream.docker-user.config"
    assert config.read_text(encoding="utf-8") == "docker {\n  runOptions = '--user 24701:24703'\n}\n"
    command = build_downstream_nextflow_command(run, docker_user_config=config)
    assert command[:2] == ["nextflow", "run"]
    assert command[command.index("-profile") + 1] == "docker"
    assert command.count("-c") == 2
    docker_config_index = [index for index, value in enumerate(command) if value == "-c"][-1]
    assert command[docker_config_index + 1] == str(config.resolve())
    assert "rnaseq-control-plane:latest" not in config.read_text(encoding="utf-8")
    assert "1000:1000" not in config.read_text(encoding="utf-8")
    if shutil.which("nextflow") is not None:
        probe = config.parent / "docker_user_config_probe.nf"
        probe.write_text("nextflow.enable.dsl=2\nworkflow { }\n", encoding="utf-8")
        result = subprocess.run(
            [
                "nextflow", "run", str(probe),
                "-c", str(Path(__file__).parents[1] / "workflow" / "nextflow.config"),
                "-c", str(config), "-profile", "docker",
            ],
            cwd=config.parent,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr


def test_downstream_macos_keeps_the_existing_container_user_contract(project_factory):
    _root, _report, run = _frozen_run(project_factory)
    assert write_downstream_docker_user_config(run, host_os="Darwin", uid=24701, gid=24703) is None
    command = build_downstream_nextflow_command(run)
    assert command[:2] == ["nextflow", "run"]
    assert command[command.index("-profile") + 1] == "docker"
    assert command.count("-c") == 1


def test_failed_service_run_keeps_frozen_logs_and_provenance(monkeypatch, project_factory, tmp_path):
    root = project_factory()
    report = validate_project(root)
    generate_plan(report)
    local_root = tmp_path / "local-nextflow-cache"
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(local_root))
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "Docker daemon is available."))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))

    def fail_nextflow(_command, *, cwd, stdout_path, stderr_path):
        (cwd / ".nextflow" / "cache").mkdir(parents=True, exist_ok=True)
        stdout_path.write_text("failed stdout\n", encoding="utf-8")
        stderr_path.write_text("failed stderr\n", encoding="utf-8")
        return 23

    monkeypatch.setattr("rnaseq.service._run_command", fail_nextflow)
    with pytest.raises(UpstreamExecutionError, match="return code 23"):
        execute_service_run(report, case_id="CASE-20260828-002")
    state_path = next((root / "runs" / "CASE-20260828-002").glob("*/run_state.json"))
    run_dir = state_path.parent
    assert json.loads(state_path.read_text(encoding="utf-8"))["status"] == "FAILED"
    assert (run_dir / "frozen" / "project.yaml").is_file()
    assert (run_dir / "provenance" / "run_provenance.yaml").is_file()
    assert (run_dir / "logs" / "downstream.stderr.log").is_file()
    assert (local_root / "CASE-20260828-002" / run_dir.name / "launch" / ".nextflow").is_dir()


def test_service_preflight_accepts_a_successful_container_probe(monkeypatch, project_factory):
    root = project_factory()
    report = validate_project(root)
    generate_plan(report)
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))
    prepare_service_run(report, profile="local")


def test_upstream_reuse_requires_matching_frozen_contract(project_factory):
    root, report, run = _frozen_run(project_factory)
    current = run.run_dir / "frozen"
    for name in ("nfcore.params.json", "nfcore.local.config"):
        (current / name).write_text("same\n", encoding="utf-8")
    source = root / "runs" / "CASE-OLD" / "20260828-090000+0800"
    (source / "frozen").mkdir(parents=True)
    (source / "upstream" / "nfcore_rnaseq").mkdir(parents=True)
    (source / "upstream" / "nfcore_rnaseq" / "marker.txt").write_text("upstream", encoding="utf-8")
    (source / "handoff").mkdir()
    (source / "handoff" / "upstream_manifest.yaml").write_text("pipeline: {}\n", encoding="utf-8")
    (source / "run_state.json").write_text('{"status":"SUCCESS"}\n', encoding="utf-8")
    for name in ("input_manifest.yaml", "nfcore.params.json", "nfcore.local.config"):
        (source / "frozen" / name).write_bytes((current / name).read_bytes())
    frozen = FrozenInputs(None, current / "downstream_contract.json", current / "input_manifest.yaml", current / "nfcore.params.json", current / "nfcore.local.config", {})
    assert reuse_upstream_if_compatible(run, frozen, "CASE-OLD/20260828-090000+0800") == "CASE-OLD/20260828-090000+0800"
    assert (run.run_dir / "upstream" / "nfcore_rnaseq" / "marker.txt").is_file()
    (source / "frozen" / "nfcore.params.json").write_text("different\n", encoding="utf-8")
    other = create_case_run(report, "CASE-20260828-001", moment=datetime(2026, 8, 28, 11, 51, 12))
    other_frozen = other.run_dir / "frozen"
    for name in ("input_manifest.yaml", "nfcore.params.json", "nfcore.local.config"):
        (other_frozen / name).write_bytes((current / name).read_bytes())
    with pytest.raises(ExecutionPreflightError, match="incompatible"):
        reuse_upstream_if_compatible(other, frozen, "CASE-OLD/20260828-090000+0800")


def test_raw_counts_delivery_preserves_metadata_sample_order(project_factory):
    _root, _report, run = _frozen_run(project_factory)

    metadata = run.run_dir / "frozen" / "metadata.csv"
    counts = run.run_dir / "frozen" / "input" / "counts.csv"

    # Deliberately use a valid non-alphabetical sample order.
    expected_samples = ["C1", "T1", "C2", "T2", "C3", "T3"]

    metadata.write_text(
        "sample_id,condition\n"
        "C1,Control\n"
        "T1,Treatment\n"
        "C2,Control\n"
        "T2,Treatment\n"
        "C3,Control\n"
        "T3,Treatment\n",
        encoding="utf-8",
    )

    counts.write_text(
        "gene_id,C1,T1,C2,T2,C3,T3\n"
        "GeneA,1,4,2,5,3,6\n",
        encoding="utf-8",
    )

    delivery = assemble_delivery(run)

    target = delivery / "counts" / "raw_counts.csv"
    assert target.read_bytes() == counts.read_bytes()

    manifest = json.loads(
        (delivery / "counts" / "artifact_manifest.json").read_text(encoding="utf-8")
    )
    artifact = manifest["artifacts"][0]

    assert artifact["ordered_sample_ids"] == expected_samples
    assert artifact["columns"] == 6
