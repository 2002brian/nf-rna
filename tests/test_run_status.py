"""``rnaseq status`` and per-run execution logging."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app
from rnaseq.downstream_runtime import DownstreamRuntime
from rnaseq.execution import LocalResourceCapacity, RuntimeCheck
from rnaseq.planner import generate_plan
from rnaseq.run_status import _Cache, build_status, current_process_identity, render_status, select_run
from rnaseq.validators import validate_project

runner = CliRunner()
RUN_ID = "20260925-203032+0800"
TRACE_HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"


# ------------------------------------------------------------------ shared helpers (also used by test_resource_policy)

def salmon_case_project(tmp_path: Path, monkeypatch, *, execution: dict | None = None):
    """A validated, planned Salmon FASTQ project whose runtime checks are mocked."""

    root = tmp_path / "project"
    (root / "input" / "fastq").mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (root / "input" / "fastq" / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None, "strandedness": "auto", "quantification": {"method": "salmon"}}
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    if execution is not None:
        config["execution"] = execution
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution-root"))
    capacity = LocalResourceCapacity(16, 64, 60)
    monkeypatch.setattr("rnaseq.execution.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.service.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    monkeypatch.setattr("rnaseq.service.check_upstream_conda", lambda: RuntimeCheck("Conda", "FOUND", "conda 25.3.1"))
    runtime = DownstreamRuntime(
        prefix=tmp_path / "native-runtime", platform="linux-64",
        lock_filename="nf-rna-downstream-linux-64.lock.yml", lock_sha256="a" * 64,
        wheel_filename="nf_rna-1.3.0-py3-none-any.whl", wheel_sha256="b" * 64,
        nf_rna_version="1.3.0", source_revision="test-revision",
        r_scripts=({"name": "l2_analysis.R", "sha256": "c" * 64},), r_scripts_sha256="d" * 64,
    )
    monkeypatch.setattr("rnaseq.service.downstream_runtime_preflight", lambda: None)
    monkeypatch.setattr("rnaseq.service.ensure_downstream_runtime", lambda: runtime)
    return validate_project(root)


def fake_nextflow_factory(observed: list[list[str]], *, on_upstream=None):
    """A successful fake Nextflow that writes the documented outputs and plain-log lines."""

    def fake(command, *, cwd, stdout_path, stderr_path):
        observed.append(command)
        outdir = Path(command[command.index("--outdir") + 1])
        if "nf-core/rnaseq" in command:
            if on_upstream is not None:
                on_upstream(command)
            stdout_path.write_text("[ab/cdef01] Submitted process > NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (C1)\n", encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            (outdir / "salmon").mkdir(parents=True)
            (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\tC2\tT1\tT2\nGeneA\t1.0\t2.0\t3.0\t4.0\n", encoding="utf-8")
            (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
            for sample in ("C1", "C2", "T1", "T2"):
                (outdir / "salmon" / sample).mkdir()
                (outdir / "salmon" / sample / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
            (outdir / "multiqc" / "multiqc_data").mkdir(parents=True)
            (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
            (outdir / "pipeline_info").mkdir()
            (outdir / "pipeline_info" / "execution_trace_2026-09-25_20-34-41.txt").write_text(
                TRACE_HEADER + "1\tab/cdef01\t1\tNFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (C1)\tCOMPLETED\t0\n", encoding="utf-8",
            )
        else:
            stdout_path.write_text("downstream\n", encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            (outdir / "report").mkdir(parents=True)
            (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return 0

    return fake


def _run_dir(project: Path, case: str = "grcm39-v130-salmon-r6", run_id: str = RUN_ID) -> Path:
    run_dir = project / "runs" / case / run_id
    for child in ("frozen", "logs", "provenance", "upstream", "downstream", "delivery"):
        (run_dir / child).mkdir(parents=True, exist_ok=True)
    return run_dir


def _state(run_dir: Path, **values) -> None:
    state = {"case_id": run_dir.parent.name, "run_id": run_dir.name, "timezone": "Asia/Taipei", "started_at": "2026-09-25T20:30:32+08:00", "completed_at": None}
    state.update(values)
    (run_dir / "run_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")


def _trace(path: Path, rows: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TRACE_HEADER + "".join(f"{index}\t{task_hash}\t{index}\t{name}\t{status}\t0\n" for index, (task_hash, name, status) in enumerate(rows, 1)), encoding="utf-8")


def _upstream_trace(run_dir: Path) -> Path:
    return run_dir / "upstream" / "nfcore_rnaseq" / "pipeline_info" / "execution_trace_2026-09-25_20-34-41.txt"


def _tree_fingerprint(root: Path) -> dict[str, tuple[int, str]]:
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "dir")
        for path in sorted(root.rglob("*"))
    }


def _dead_pid() -> int:
    process = subprocess.Popen([sys.executable, "-c", "pass"])
    process.wait()
    return process.pid


def _norm(text: str) -> str:
    """Collapse column padding so assertions do not depend on the widest stage label."""

    return "\n".join(" ".join(line.split()) for line in text.splitlines()) + ("\n" if text.endswith("\n") else "")


def _status(project: Path, *arguments: str) -> str:
    result = runner.invoke(app, ["status", str(project), *arguments])
    assert result.exit_code == 0, result.output
    return _norm(result.output)


# ------------------------------------------------------------------ route-aware fixtures

SAMPLES = ("S1", "S2", "S3")
NOW = datetime(2026, 9, 25, 21, 0, 0).astimezone()


def _route_run(tmp_path: Path, method: str, *, case: str = "CASE-A", run_id: str = RUN_ID, enrichment=(), preset: str = "L2") -> Path:
    run_dir = _run_dir(tmp_path, case, run_id)
    (run_dir / "frozen" / "input_manifest.yaml").write_text(yaml.safe_dump({
        "project": {"id": "PRJ"}, "input": {"type": "fastq"},
        "upstream": {"engine": "nfcore_rnaseq", "quantification": {"method": method}},
        "planned_downstream": {"preset": preset, "enrichment": list(enrichment)},
    }), encoding="utf-8")
    (run_dir / "frozen" / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2,strandedness\n" + "".join(f"{s},{s}_R1.fq.gz,{s}_R2.fq.gz,unstranded\n" for s in SAMPLES), encoding="utf-8",
    )
    return run_dir


def _hisat2_tasks(samples=SAMPLES, *, qc=True, align=True, count=True, final=True):
    rows = []
    for index, sample in enumerate(samples):
        if qc:
            rows += [(f"a{index}/00000{n}", name, "COMPLETED") for n, name in enumerate((f"FASTQC_RAW (raw {sample})", f"FASTP_PREPARE ({sample})", f"FASTQC_PROCESSED (processed {sample})"))]
        if align:
            rows += [(f"b{index}/00000{n}", name, "COMPLETED") for n, name in enumerate((f"HISAT2_ALIGN ({sample})", f"SORT_LANE_BAM ({sample})", f"MERGE_AND_INDEX ({sample})"))]
        if count:
            rows += [(f"c{index}/00000{n}", name, "COMPLETED") for n, name in enumerate((f"PREPARE_COUNT_BAM ({sample})", f"FEATURECOUNTS ({sample})"))]
    if final:
        rows += [("d0/000000", "ASSEMBLE_COUNTS (canonical featureCounts matrix)", "COMPLETED"), ("d0/000001", "MULTIQC (MultiQC)", "COMPLETED")]
    return rows


def _downstream_tasks(enrichment_backends=()):
    rows = [("e0/000000", "L1_ANALYSIS (L1 expression QC)", "COMPLETED"), ("e0/000001", "L2_ANALYSIS (L2 DESeq2)", "COMPLETED")]
    rows += [(f"e1/00000{n}", f"ENRICHMENT_ANALYSIS ({module})", "COMPLETED") for n, module in enumerate(enrichment_backends)]
    rows.append(("e2/000000", "TECHNICAL_REPORT (HTML technical report)", "COMPLETED"))
    return rows


def _resources(run_dir: Path) -> None:
    (run_dir / "provenance" / "run_provenance.yaml").write_text(yaml.safe_dump({"runtime_resources": {
        "requested": {"cpus": 24, "memory_gib": 79}, "effective": {"cpus": 24, "memory_gib": 79}, "clamped": False,
        "policy": {"cpu_mode": "auto", "memory_mode": "auto", "detected": {"usable_cpus": 28, "usable_memory_gib": 94.0},
                   "os_reserve": {"cpus": 4, "memory_gib": 15}, "process_tuning": []},
    }}), encoding="utf-8")


def _render(run_dir: Path, now: datetime = NOW) -> str:
    return _norm(render_status(build_status(*select_run(run_dir.parents[2], run_dir.parent.name, run_dir.name), cache=_Cache()), now=now))


def _task_dir(work: Path, task_hash: str, *, begin: float | None, exited: bool = False) -> Path:
    prefix, rest = task_hash.split("/")
    task_dir = work / prefix / f"{rest}0123456789abcdef"
    task_dir.mkdir(parents=True)
    if begin is not None:
        (task_dir / ".command.begin").write_text("", encoding="utf-8")
        os.utime(task_dir / ".command.begin", (begin, begin))
    if exited:
        (task_dir / ".exitcode").write_text("0", encoding="utf-8")
    return task_dir


# ------------------------------------------------------------------ completed / failed / interrupted

def test_success_dashboard_for_the_hisat2_featurecounts_route(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="SUCCESS", phase="delivery", started_at="2026-09-26T17:43:38+08:00",
           completed_at="2026-09-26T19:04:12+08:00", delivery=str(run_dir / "delivery"))
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks())
    _trace(run_dir / "provenance" / "downstream.trace.txt", _downstream_tasks())
    _resources(run_dir)

    output = _status(tmp_path)

    assert output.startswith("PRJ / CASE-A\n━━━━")
    assert "\nSUCCESS 1h 20m 34s\n" in output
    assert (
        "Pipeline (HISAT2 → featureCounts)\n"
        "✓ Input / preflight\n"
        "✓ FASTQ QC / preprocessing 3 / 3 samples\n"
        "✓ HISAT2 alignment 3 / 3 samples\n"
        "✓ featureCounts 3 / 3 samples\n"
        "✓ MultiQC\n"
        "✓ L1 expression QC\n"
        "✓ L2 / DESeq2\n"
        "✓ Technical report\n"
        "✓ Delivery\n"
    ) in output
    assert "Tasks\n29 completed / 0 failed\n" in output
    assert f"Started 2026-09-26 17:43:38\nCompleted 2026-09-26 19:04:12\nRun {RUN_ID}\n" in output
    assert f"Delivery {run_dir / 'delivery'}" in output
    assert "Salmon" not in output and "tximport" not in output and "Enrichment" not in output
    assert "Running now" not in output and "Last update" not in output


def test_salmon_route_shows_salmon_and_tximport_stages_not_hisat2(tmp_path):
    run_dir = _route_run(tmp_path, "salmon")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    prefix = "NFCORE_RNASEQ:RNASEQ:"
    rows = [("f0/000000", "NFCORE_RNASEQ:PREPARE_GENOME:CUSTOM_GTFFILTER (genome.gtf)", "COMPLETED")]
    for index, sample in enumerate(SAMPLES):
        rows += [
            (f"f1/00000{index}", f"{prefix}FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE ({sample})", "COMPLETED"),
            (f"f2/00000{index}", f"{prefix}QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT ({sample})", "CACHED"),
        ]
    rows += [
        ("f3/000000", f"{prefix}QUANTIFY_PSEUDO_ALIGNMENT:QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT:TXIMETA_TXIMPORT (salmon)", "COMPLETED"),
        ("f4/000000", f"{prefix}MULTIQC_RNASEQ:MULTIQC", "COMPLETED"),
    ]
    _trace(_upstream_trace(run_dir), rows)

    output = _status(tmp_path)

    assert "Pipeline (nf-core/rnaseq Salmon → tximport)" in output
    for line in ("✓ Reference preparation", "✓ FASTQ QC / trimming 3 / 3 samples", "✓ Salmon quantification 3 / 3 samples",
                 "✓ tximport / gene summary", "✓ MultiQC"):
        assert line in output
    assert "HISAT2" not in output and "featureCounts" not in output
    assert "Tasks\n9 completed (3 cached) / 0 failed" in output


def test_enrichment_stage_is_shown_only_when_configured(tmp_path):
    configured = _route_run(tmp_path / "with", "salmon", enrichment=["go", "gsea"])
    _state(configured, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    _trace(configured / "provenance" / "downstream.trace.txt", _downstream_tasks(("go", "gsea-go", "gsea-kegg")))
    assert "✓ Enrichment (go, gsea)" in _status(tmp_path / "with")

    plain = _route_run(tmp_path / "without", "hisat2_featurecounts")
    _state(plain, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    _trace(plain / "provenance" / "downstream.trace.txt", _downstream_tasks())
    assert "Enrichment" not in _status(tmp_path / "without")

    l1_only = _route_run(tmp_path / "l1", "hisat2_featurecounts", preset="L1")
    _state(l1_only, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    output = _status(tmp_path / "l1")
    assert "✓ L1 expression QC" in output and "L2 / DESeq2" not in output


def test_enrichment_is_running_until_every_selected_backend_finished(tmp_path):
    run_dir = _route_run(tmp_path, "salmon", enrichment=["gsea"])
    _state(run_dir, status="RUNNING", phase="downstream", downstream_command=["nextflow", "-work-dir", str(tmp_path / "work")])
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    _trace(run_dir / "provenance" / "downstream.trace.txt", _downstream_tasks(("gsea-go",))[:-1])
    (run_dir / "logs" / "downstream.stdout.log").write_text("[e1/000009] Submitted process > ENRICHMENT_ANALYSIS (gsea-kegg)\n", encoding="utf-8")
    _task_dir(tmp_path / "work", "e1/000009", begin=NOW.timestamp() - 65)
    output = _render(run_dir)
    assert "▶ Enrichment (gsea) 1 done, 1 running" in output
    assert "○ Technical report waiting" in output
    assert "gsea-kegg ENRICHMENT_ANALYSIS 01:05" in output


def test_failed_run_marks_the_failing_stage_and_later_stages_not_started(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="FAILED", phase="upstream", completed_at="2026-09-25T23:00:00+08:00",
           error="upstream Nextflow failed with return code 1. Logs: x")
    rows = _hisat2_tasks(count=False, final=False) + [("c0/000000", "FEATURECOUNTS (S1)", "COMPLETED"), ("c1/000000", "FEATURECOUNTS (S2)", "FAILED")]
    _trace(run_dir / "provenance" / "upstream.trace.txt", rows)

    output = _status(tmp_path)

    assert "\nFAILED 2h 29m 28s\n" in output
    assert "✓ HISAT2 alignment 3 / 3 samples" in output
    assert "! featureCounts 1 failed task" in output
    assert "– MultiQC not started" in output
    assert "– L2 / DESeq2 not started" in output and "– Delivery" in output
    assert "Tasks\n19 completed / 1 failed" in output
    assert "Error: upstream Nextflow failed with return code 1." in output


def test_a_failed_attempt_that_succeeded_on_retry_is_not_a_failure(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [("aa/000001", "HISAT2_ALIGN (S1)", "FAILED"), *_hisat2_tasks()])
    output = _status(tmp_path)
    assert "0 failed\n1 failed attempt later succeeded on retry" in output


def test_completed_stage_never_pairs_a_check_mark_with_partial_sample_counts(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks(("S1", "S2"), count=False, final=False))
    output = _status(tmp_path)
    assert "✓ FASTQ QC / preprocessing\n" in output and "2 / 3" not in output


def test_recorded_interrupted_run_marks_the_stage_that_was_running(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="INTERRUPTED", phase="upstream", completed_at="2026-09-25T20:40:00+08:00",
           interrupted_by="SIGTERM", error="Interrupted by SIGTERM; the run did not finish.")
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks(align=False, count=False, final=False))
    (run_dir / "logs" / "upstream.stdout.log").write_text("[bb/000001] Submitted process > HISAT2_ALIGN (S1)\n", encoding="utf-8")

    output = _status(tmp_path)

    assert "\nINTERRUPTED 9m 28s\n" in output
    assert "✓ FASTQ QC / preprocessing 3 / 3 samples" in output
    assert "⏸ HISAT2 alignment 0 / 3 samples" in output
    assert "– featureCounts not started" in output
    assert "Error: Interrupted by SIGTERM" in output
    assert "Running now" not in output  # nothing is running in a finished run


def test_interrupted_without_task_evidence_is_reported_at_phase_level(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="INTERRUPTED", phase="upstream", completed_at="2026-09-25T20:40:00+08:00")
    output = _status(tmp_path)
    assert "⏸ Upstream no task trace recorded" in output and "– Delivery" in output


def test_failed_run_lists_the_failed_task_names(tmp_path):
    run_dir = _route_run(tmp_path, "salmon")
    _state(run_dir, status="FAILED", phase="upstream", completed_at="2026-09-25T23:00:00+08:00")
    name = "NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S3)"
    _trace(_upstream_trace(run_dir), [("aa/000001", name, "FAILED")])
    output = _status(tmp_path)
    assert f"Failed tasks\n{name}\n" in output


def test_tasks_aborted_as_a_side_effect_are_stopped_not_failed(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="FAILED", phase="upstream", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [
        ("aa/000001", "HISAT2_ALIGN (S2)", "FAILED"), ("aa/000002", "FASTQC_RAW (raw S3)", "ABORTED"), ("aa/000003", "FASTP_PREPARE (S3)", "ABORTED"),
    ])
    output = _status(tmp_path)
    assert "! HISAT2 alignment 1 failed task" in output
    assert "! FASTQ QC" not in output and "Tasks\n0 completed / 1 failed" in output


def test_interrupted_run_with_aborted_tasks_shows_interrupted_not_failed(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="INTERRUPTED", phase="upstream", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [("aa/000001", "HISAT2_ALIGN (S1)", "ABORTED")])
    output = _status(tmp_path)
    assert "⏸ HISAT2 alignment" in output and "!" not in output.split("Pipeline")[1].split("\n\n")[0]


def test_a_task_being_retried_while_running_is_not_shown_as_failed(tmp_path):
    work = tmp_path / "work"
    run_dir = _route_run(tmp_path, "salmon")
    _state(run_dir, status="RUNNING", phase="upstream", upstream_command=["nextflow", "-work-dir", str(work)])
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    name = "NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S3)"
    _trace(_upstream_trace(run_dir), [("aa/000001", name, "FAILED")])
    (run_dir / "logs" / "upstream.stdout.log").write_text(f"[aa/000001] Submitted process > {name}\n[bb/000002] Submitted process > {name}\n", encoding="utf-8")
    _task_dir(work, "bb/000002", begin=NOW.timestamp() - 10)
    output = _render(run_dir)
    assert "! Salmon" not in output and "▶ Salmon quantification 0 / 3 samples" in output
    assert "1 running / 0 queued / 0 failed\n" in output and "1 failed attempt later succeeded on retry" in output


def test_status_follows_a_work_dir_outside_the_execution_root(tmp_path):
    launch = tmp_path / "linux-cache" / "CASE-A" / RUN_ID / "launch"
    work = tmp_path / "large-volume" / "CASE-A" / RUN_ID / "work" / "upstream"
    run_dir = _route_run(tmp_path / "project", "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", phase="upstream", upstream_command=["nextflow", "run", "-work-dir", str(work)])
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    (run_dir / "provenance" / "run_provenance.yaml").write_text(yaml.safe_dump({
        "execution_root": str(launch.parent), "execution_launch_dir": str(launch), "execution_work_dir": str(work.parent),
    }), encoding="utf-8")
    launch.mkdir(parents=True)
    (launch / ".nextflow.log").write_text("heartbeat\n", encoding="utf-8")
    heartbeat = NOW.timestamp() - 5
    os.utime(launch / ".nextflow.log", (heartbeat, heartbeat))
    for path in (run_dir / "run_state.json", run_dir / "logs" / "rnaseq.process.json", run_dir / "provenance" / "run_provenance.yaml"):
        os.utime(path, (NOW.timestamp() - 600, NOW.timestamp() - 600))
    name = "HISAT2_ALIGN (S1)"
    (run_dir / "logs" / "upstream.stdout.log").write_text(f"[bb/000002] Submitted process > {name}\n", encoding="utf-8")
    os.utime(run_dir / "logs" / "upstream.stdout.log", (NOW.timestamp() - 600, NOW.timestamp() - 600))
    _task_dir(work, "bb/000002", begin=NOW.timestamp() - 65)
    output = _render(run_dir)
    assert "1 running / 0 queued / 0 failed\n" in output
    assert "S1 HISAT2_ALIGN 01:05" in output  # elapsed read from the relocated task directory
    assert "last write 0m 05s ago" in output  # heartbeat read from the launch directory under the execution root


def _multilane(run_dir: Path) -> None:
    (run_dir / "frozen" / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2,strandedness\n" + "".join(f"{s},{s}_L{l}_R1.fq.gz,{s}_L{l}_R2.fq.gz,unstranded\n" for s in SAMPLES for l in (1, 2)), encoding="utf-8",
    )


def test_multi_lane_qc_counts_a_sample_only_when_every_lane_finished(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _multilane(run_dir)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    rows = [(f"a{i}/00000{l}", f"FASTQC_PROCESSED (processed {s})", "COMPLETED") for i, s in enumerate(SAMPLES) for l in (1, 2) if l == 1 or s == "S1"]
    _trace(run_dir / "provenance" / "upstream.trace.txt", rows)
    output = _render(run_dir)
    assert "▶ FASTQ QC / preprocessing 1 / 3 samples" in output


def test_multi_lane_failure_is_not_hidden_by_the_other_lanes_success(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _multilane(run_dir)
    _state(run_dir, status="FAILED", phase="upstream", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [("aa/000001", "HISAT2_ALIGN (S1)", "COMPLETED"), ("aa/000002", "HISAT2_ALIGN (S1)", "FAILED")])
    output = _status(tmp_path)
    assert "! HISAT2 alignment 1 failed task" in output and "1 failed" in output.split("Tasks")[1] and "retry" not in output


def test_ascii_mode_is_pure_ascii(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T23:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks())
    text = render_status(build_status(*select_run(tmp_path), cache=_Cache()), ascii_only=True)
    text.encode("ascii")
    assert "Pipeline (HISAT2 -> featureCounts)" in text and "+ HISAT2 alignment" in text


def test_legacy_flat_run_names_the_project_and_a_legacy_case(tmp_path):
    project = tmp_path / "MyProject"
    run_dir = project / "runs" / RUN_ID
    run_dir.mkdir(parents=True)
    (run_dir / "run_state.json").write_text(json.dumps({"run_id": RUN_ID, "status": "SUCCESS", "phase": "delivery",
                                                         "started_at": "2026-09-25T20:30:32+08:00", "completed_at": "2026-09-25T21:00:00+08:00"}), encoding="utf-8")
    assert _status(project).startswith("MyProject / legacy\n")


# ------------------------------------------------------------------ live and stale RUNNING

def test_running_dashboard_shows_stage_progress_running_tasks_and_counts(tmp_path):
    work = tmp_path / "work" / "upstream"
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", phase="upstream", upstream_command=["nextflow", "run", "-work-dir", str(work)])
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    (run_dir / "logs" / "rnaseq.log").write_text("2026-09-25T20:30:32+08:00 state RUNNING phase=upstream\n", encoding="utf-8")
    _resources(run_dir)
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks(align=False, count=False, final=False) + _hisat2_tasks(("S1",), qc=False, count=False, final=False))
    (run_dir / "logs" / "upstream.stdout.log").write_text(
        "N E X T F L O W\n"
        "[b0/000000] Submitted process > HISAT2_ALIGN (S1)\n"
        "[bb/000002] Submitted process > HISAT2_ALIGN (S2)\n"
        "[bb/000003] Submitted process > HISAT2_ALIGN (S3)\n"
        "[bb/000004] Submitted process > SORT_LANE_BAM (S1)\n",
        encoding="utf-8",
    )
    _task_dir(work, "bb/000002", begin=NOW.timestamp() - 494)  # started 8m14s ago
    _task_dir(work, "bb/000003", begin=None)  # submitted, not started: queued
    _task_dir(work, "bb/000004", begin=NOW.timestamp() - 30, exited=True)  # finished, trace row pending

    output = _render(run_dir)

    assert "\nRUNNING 29m 28s elapsed\n" in output
    assert "✓ FASTQ QC / preprocessing 3 / 3 samples" in output
    assert "▶ HISAT2 alignment 1 / 3 samples" in output
    assert "○ featureCounts waiting" in output and "○ MultiQC waiting" in output
    assert "○ L2 / DESeq2 waiting" in output and "○ Delivery" in output
    assert "Running now\nS2 HISAT2_ALIGN 08:14\n" in output
    assert "S3 HISAT2_ALIGN" not in output
    assert "Tasks\n12 completed / 1 running / 2 queued / 0 failed\n(further tasks are created as Nextflow progresses; total not yet known)" in output
    assert "Last update 21:00:00" in output and "Activity last write" in output
    assert f"Log {run_dir / 'logs' / 'rnaseq.log'}" in output


def test_running_task_without_a_recorded_work_dir_has_no_invented_elapsed_time(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", phase="upstream")  # no upstream_command recorded yet
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    (run_dir / "logs" / "upstream.stdout.log").write_text("[aa/000001] Submitted process > FASTQC_RAW (raw S1)\n", encoding="utf-8")
    output = _render(run_dir)
    assert "○ FASTQ QC / preprocessing 1 queued" in output
    assert "Running now" not in output and "0 running / 1 queued" in output  # submitted, start time unknown


def test_resource_ceiling_is_labelled_as_configured_not_live_usage(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    _resources(run_dir)
    output = _render(run_dir)
    assert "Resources (configured Nextflow ceiling, not live usage)\nCPU 24 max\nMemory 79 GiB max\n" in output
    assert "Policy auto; usable 28 CPUs / 94.0 GiB; OS reserve 4 CPUs / 15 GiB" in output
    assert " / 24" not in output and " / 79" not in output  # never presented as used / allowed


def test_stale_running_state_whose_process_is_gone_is_interrupted_not_success(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", phase="downstream", downstream_command=["nextflow", "-work-dir", str(tmp_path / "work")])
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "pid": _dead_pid()}), encoding="utf-8")
    _trace(run_dir / "provenance" / "upstream.trace.txt", _hisat2_tasks())
    (run_dir / "logs" / "downstream.stdout.log").write_text("[e0/000000] Submitted process > L1_ANALYSIS (L1 expression QC)\n", encoding="utf-8")
    _task_dir(tmp_path / "work", "e0/000000", begin=NOW.timestamp() - 5)  # a leftover begin marker never makes it "running"
    # Even complete-looking outputs never make a dead RUNNING run a SUCCESS.
    (run_dir / "delivery" / "README.md").write_text("x", encoding="utf-8")

    output = _status(tmp_path)

    assert "\nINTERRUPTED" in output and "SUCCESS" not in output
    assert "Note: stale: recorded RUNNING, but rnaseq process" in output and "no longer exists" in output
    assert "✓ featureCounts" in output and "⏸ L1 expression QC" in output and "– Delivery" in output
    assert "Running now" not in output and "running" not in output.split("Tasks")[1].split("\n\n")[0]


def test_reused_pid_is_detected_by_process_start_time(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    identity = current_process_identity()
    if identity["start_ticks"] is None:
        pytest.skip("procfs start times unavailable")
    identity["start_ticks"] -= 1
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(identity), encoding="utf-8")
    assert "\nINTERRUPTED" in _status(tmp_path)


def test_legacy_running_state_without_process_record_is_unverified(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    output = _status(tmp_path)
    assert "\nRUNNING" in output
    assert "Note: unverified: no process identity was recorded" in output
    assert "▶ Upstream starting; no task submitted yet" in output
    assert "○ Downstream waiting" in output and "L2 / DESeq2" not in output  # no frozen plan: modules not guessed


def test_running_state_from_another_host_is_unverified(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({"pid": 1, "hostname": "some-other-host"}), encoding="utf-8")
    assert "unverified: launched on host some-other-host" in _status(tmp_path)


# ------------------------------------------------------------------ missing / partial artifacts, read-only, legacy runs

def test_status_with_missing_or_corrupt_optional_artifacts(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "provenance" / "run_provenance.yaml").write_text(": not [valid yaml", encoding="utf-8")
    (run_dir / "frozen" / "input_manifest.yaml").write_text(": not [valid yaml", encoding="utf-8")
    _upstream_trace(run_dir).parent.mkdir(parents=True)
    _upstream_trace(run_dir).write_text("garbage without header\n", encoding="utf-8")
    (run_dir / "logs" / "upstream.stdout.log").write_bytes(b"\xff\xfe binary \x00\n")
    # A malformed state of another run is skipped, not fatal.
    broken = _run_dir(tmp_path, case="OTHER")
    (broken / "run_state.json").write_text("{", encoding="utf-8")

    output = _status(tmp_path)
    assert " / grcm39-v130-salmon-r6\n" in output
    assert "▶ Upstream" in output


def test_status_with_only_a_run_state_and_no_runs(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="CREATED")
    output = _status(tmp_path)
    assert "○ Upstream" in output and "○ Delivery" in output

    empty = tmp_path / "empty"
    (empty / "planning").mkdir(parents=True)
    (empty / "planning" / "manifest.preview.yaml").write_text("{}", encoding="utf-8")
    assert "PLANNED: no case runs recorded yet" in _status(empty)
    nothing = tmp_path / "nothing-yet"
    nothing.mkdir()
    assert "No recorded case runs." in _status(nothing)


def test_raw_count_runs_report_upstream_not_applicable(tmp_path):
    run_dir = _run_dir(tmp_path, case="RAW")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T20:40:00+08:00")
    (run_dir / "frozen" / "input_manifest.yaml").write_text(yaml.safe_dump({"input": {"type": "raw_counts"}}), encoding="utf-8")
    output = _status(tmp_path)
    assert "Pipeline (raw counts → DESeq2)" in output and "– Upstream not applicable (raw-count input)" in output


def test_qc_preset_reports_downstream_skipped(tmp_path):
    run_dir = _route_run(tmp_path, "hisat2_featurecounts", preset="QC")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T20:40:00+08:00", downstream_skipped="technical_qc_only")
    output = _status(tmp_path)
    assert "– L1 / L2" in output and "skipped (technical QC only)" in output and "DESeq2" not in output


def test_released_v130_run_layout_remains_readable(tmp_path):
    """A run persisted by v1.3.0 (process record of a long-gone PID, policy resources) and a
    pre-9719374 run (no rnaseq.log, no process record, no resource policy) both render."""

    v130 = _route_run(tmp_path, "hisat2_featurecounts", case="V130-HISAT2", run_id="20260926-174337+0800")
    _state(v130, status="SUCCESS", phase="delivery", started_at="2026-09-26T17:43:38+08:00", completed_at="2026-09-26T19:04:12+08:00",
           delivery=str(v130 / "delivery"), command=["rnaseq", "run"], upstream_command=["nextflow"], downstream_command=["nextflow"])
    (v130 / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "pid": _dead_pid()}), encoding="utf-8")
    (v130 / "logs" / "rnaseq.log").write_text("x\n", encoding="utf-8")
    _trace(v130 / "provenance" / "upstream.trace.txt", _hisat2_tasks())
    _resources(v130)
    older = _route_run(tmp_path, "salmon", case="grcm39-v130-salmon-r6", run_id="20260925-203032+0800", enrichment=["go", "kegg"])
    _state(older, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    (older / "logs" / "upstream.stdout.log").write_text("", encoding="utf-8")
    (older / "provenance" / "run_provenance.yaml").write_text(yaml.safe_dump({"runtime_resources": {
        "requested": {"cpus": 22, "memory_gib": 20}, "effective": {"cpus": 22, "memory_gib": 20}}}), encoding="utf-8")

    latest = _status(tmp_path)
    assert "PRJ / V130-HISAT2" in latest and "\nSUCCESS 1h 20m 34s" in latest and "Note:" not in latest
    legacy = _status(tmp_path, "--case", "grcm39-v130-salmon-r6")
    assert "\nSUCCESS 4h 09m 04s" in legacy and "✓ Enrichment (go, kegg)" in legacy
    assert "Requested 22 CPUs / 20 GiB (pre-policy run)" in legacy
    assert "(no per-run rnaseq.log; run predates it)" in legacy


REAL_SJ = Path("/mnt/d/mrpan_rna/SJ_ctrl_vs_Gy22")


@pytest.mark.skipif(not (REAL_SJ / "runs" / "V130-HISAT2" / "20260926-174337+0800" / "run_state.json").is_file(), reason="real SJ project not present")
def test_real_completed_sj_v130_run_is_success_and_untouched():
    run_dir = REAL_SJ / "runs" / "V130-HISAT2" / "20260926-174337+0800"
    before = _tree_fingerprint(run_dir / "logs")
    output = _status(REAL_SJ, "--case", "V130-HISAT2", "--run", "20260926-174337+0800")
    assert "SJ / V130-HISAT2" in output and "\nSUCCESS 1h 20m 34s" in output
    assert "✓ HISAT2 alignment 7 / 7 samples" in output and "✓ featureCounts 7 / 7 samples" in output
    assert _tree_fingerprint(run_dir / "logs") == before


def test_status_never_modifies_the_project(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "pid": _dead_pid()}), encoding="utf-8")
    _trace(_upstream_trace(run_dir), [("aa/000001", "FASTQC", "COMPLETED")])
    before = _tree_fingerprint(tmp_path)
    _status(tmp_path)
    _status(tmp_path, "--all")
    runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "2"])  # final (stale) state: prints once
    assert _tree_fingerprint(tmp_path) == before


def test_status_never_contains_scientific_results(tmp_path):
    """Guard: status is execution monitoring; DE/enrichment/PCA results belong in the report."""

    run_dir = _route_run(tmp_path, "salmon", enrichment=["go", "kegg", "gsea"])
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00", delivery=str(run_dir / "delivery"))
    _trace(run_dir / "provenance" / "downstream.trace.txt", _downstream_tasks(("go", "kegg", "gsea-go", "gsea-kegg")))
    contrast = run_dir / "downstream" / "l2" / "contrasts" / "T_vs_C"
    contrast.mkdir(parents=True)
    (contrast / "backend_summary.json").write_text(json.dumps({"significant": 2165, "upregulated": 1562}), encoding="utf-8")
    (contrast / "significant.tsv").write_text("gene_id\tlog2FoldChange\tpadj\nG1\t2.0\t0.001\n", encoding="utf-8")
    (run_dir / "downstream" / "l1").mkdir(parents=True)
    (run_dir / "downstream" / "l1" / "pca_variance.tsv").write_text("component\tproportion_variance\nPC1\t0.62\n", encoding="utf-8")

    output = (_status(tmp_path) + _status(tmp_path, "--all")).lower()

    for forbidden in ("significant", "upregulated", "downregulated", "deg", "padj", "log2", "fold", "p-value", "pvalue", "pca", "pc1", "pathway", "2165", "1562"):
        assert forbidden not in output, forbidden
    source = (Path(__file__).parents[1] / "src" / "rnaseq" / "run_status.py").read_text(encoding="utf-8")
    assert "downstream/l2" not in source and '"l2"' not in source and "backend_summary" not in source and "all_genes" not in source


# ------------------------------------------------------------------ run selection and history

def _selection_project(tmp_path: Path) -> Path:
    older = _route_run(tmp_path, "hisat2_featurecounts", case="CASE-A", run_id="20260925-100000+0800")
    _state(older, status="FAILED", phase="upstream", started_at="2026-09-25T10:00:00+08:00", completed_at="2026-09-25T10:05:00+08:00", error="boom")
    retry = _route_run(tmp_path, "hisat2_featurecounts", case="CASE-A", run_id="20260925-120000+0800")
    _state(retry, status="SUCCESS", phase="delivery", started_at="2026-09-25T12:00:00+08:00", completed_at="2026-09-25T13:00:00+08:00",
           retry_of={"case_id": "CASE-A", "run_id": "20260925-100000+0800", "status": "FAILED"})
    other = _route_run(tmp_path, "salmon", case="CASE-B", run_id="20260925-110000+0800")
    _state(other, status="INTERRUPTED", phase="downstream", started_at="2026-09-25T11:00:00+08:00", completed_at="2026-09-25T11:30:00+08:00")
    return tmp_path


def test_default_selects_the_latest_run(tmp_path):
    output = _status(_selection_project(tmp_path))
    assert "PRJ / CASE-A" in output and "Run 20260925-120000+0800 (retry of CASE-A/20260925-100000+0800 (source FAILED))" in output


def test_case_selects_the_latest_run_of_that_case(tmp_path):
    project = _selection_project(tmp_path)
    assert "Run 20260925-110000+0800" in _status(project, "--case", "CASE-B")
    assert "Run 20260925-120000+0800" in _status(project, "--case", "CASE-A")


def test_run_selects_that_exact_run(tmp_path):
    project = _selection_project(tmp_path)
    output = _status(project, "--run", "20260925-100000+0800")
    assert "\nFAILED" in output and "Run 20260925-100000+0800" in output
    assert "Run 20260925-100000+0800" in _status(project, "--case", "CASE-A", "--run", "20260925-100000+0800")


def test_unknown_or_ambiguous_selection_is_an_error(tmp_path):
    project = _selection_project(tmp_path)
    missing = runner.invoke(app, ["status", str(project), "--case", "NOPE"])
    assert missing.exit_code == 1 and "no recorded run matches case NOPE" in missing.output
    wrong_pair = runner.invoke(app, ["status", str(project), "--case", "CASE-B", "--run", "20260925-100000+0800"])
    assert wrong_pair.exit_code == 1
    twin = _route_run(tmp_path, "hisat2_featurecounts", case="CASE-C", run_id="20260925-100000+0800")
    _state(twin, status="SUCCESS", phase="delivery", completed_at="2026-09-25T10:30:00+08:00")
    ambiguous = runner.invoke(app, ["status", str(project), "--run", "20260925-100000+0800"])
    assert ambiguous.exit_code == 2 and "exists in more than one case (CASE-A, CASE-C); add --case" in ambiguous.output


def test_all_is_a_compact_history_not_multiple_dashboards(tmp_path):
    project = _selection_project(tmp_path)
    output = _status(project, "--all")
    lines = output.splitlines()
    assert lines[0].split() == ["CASE", "RUN", "STATUS", "PHASE", "STARTED", "ELAPSED", "DELIVERY", "ATTEMPT"]
    assert lines[1].startswith("*") and "CASE-A" in lines[1] and "SUCCESS" in lines[1] and "retry of CASE-A/20260925-100000+0800 (source FAILED)" in lines[1]
    assert "CASE-B" in lines[2] and "INTERRUPTED" in lines[2] and "downstream" in lines[2]
    assert "CASE-A" in lines[3] and "FAILED" in lines[3] and "original" in lines[3]
    assert "━" not in output and "Pipeline" not in output
    assert "available" not in lines[1]  # no delivery directory exists in this fixture
    assert f"Run directories: {project / 'runs'}/CASE/RUN" in output
    only_b = _status(project, "--all", "--case", "CASE-B")
    assert "CASE-B" in only_b and "CASE-A" not in only_b
    assert runner.invoke(app, ["status", str(project), "--all", "--watch"]).exit_code == 2


def test_all_marks_stale_running_runs(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "pid": _dead_pid()}), encoding="utf-8")
    assert "INTERRUPTED (stale)" in _status(tmp_path, "--all")


# ------------------------------------------------------------------ --watch

def _live_run(tmp_path: Path, **values) -> Path:
    run_dir = _route_run(tmp_path, "hisat2_featurecounts")
    _state(run_dir, status="RUNNING", **values)
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    return run_dir


def test_watch_ctrl_c_stops_watching_without_touching_the_run(monkeypatch, tmp_path):
    _live_run(tmp_path, phase="upstream")
    signalled: list[tuple[int, int]] = []
    real_kill = os.kill

    def spy_kill(pid, sig):
        if sig != 0:
            signalled.append((pid, sig))
        return real_kill(pid, sig)

    def interrupt(_seconds):
        raise KeyboardInterrupt

    monkeypatch.setattr("os.kill", spy_kill)
    monkeypatch.setattr("rnaseq.cli.time.sleep", interrupt)
    before = _tree_fingerprint(tmp_path)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch"])
    assert result.exit_code == 0, result.output
    assert "\nRUNNING" in result.output
    assert "Refreshing every 10s; Ctrl-C stops watching (the run is not affected)." in result.output
    assert "Stopped watching; the run itself was not affected." in result.output
    assert "Traceback" not in result.output
    assert signalled == [] and _tree_fingerprint(tmp_path) == before


def test_watch_refreshes_until_the_run_reaches_a_final_state(monkeypatch, tmp_path):
    run_dir = _live_run(tmp_path, phase="downstream")
    sleeps: list[float] = []

    def finish(seconds):
        sleeps.append(seconds)
        _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T22:00:00+08:00")

    monkeypatch.setattr("rnaseq.cli.time.sleep", finish)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "5"])
    assert result.exit_code == 0, result.output
    assert sleeps == [5.0]
    assert result.output.count("PRJ / CASE-A") == 2
    assert result.output.split("PRJ / CASE-A")[-1].split("\n")[2].startswith("SUCCESS")


@pytest.mark.parametrize("final_state", ["FAILED", "INTERRUPTED"])
def test_watch_stops_at_every_terminal_state(monkeypatch, tmp_path, final_state):
    run_dir = _live_run(tmp_path, phase="upstream")
    ticks: list[float] = []

    def end(seconds):
        ticks.append(seconds)
        _state(run_dir, status=final_state, phase="upstream", completed_at="2026-09-25T21:00:00+08:00")

    monkeypatch.setattr("rnaseq.cli.time.sleep", end)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch"])
    assert result.exit_code == 0 and len(ticks) == 1
    assert result.output.split("PRJ / CASE-A")[-1].split("\n")[2].startswith(final_state)


def test_watch_stays_on_the_selected_run_when_a_newer_run_starts(monkeypatch, tmp_path):
    run_dir = _live_run(tmp_path, phase="upstream")

    def newer_run_then_finish(_seconds):
        newer = _route_run(tmp_path, "hisat2_featurecounts", case="NEWER", run_id="20260925-230000+0800")
        _state(newer, status="SUCCESS", phase="delivery", started_at="2026-09-25T23:00:00+08:00", completed_at="2026-09-25T23:10:00+08:00")
        _state(run_dir, status="FAILED", phase="upstream", completed_at="2026-09-25T21:00:00+08:00")

    monkeypatch.setattr("rnaseq.cli.time.sleep", newer_run_then_finish)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch"])
    assert result.exit_code == 0
    assert "NEWER" not in result.output and result.output.count("PRJ / CASE-A") == 2


def test_watch_interval_has_a_floor(tmp_path):
    result = runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "0.1"])
    assert result.exit_code != 0


def test_help_documents_only_the_implemented_options():
    result = runner.invoke(app, ["status", "--help"], env={"COLUMNS": "200", "TERM": "dumb", "NO_COLOR": "1"})
    assert result.exit_code == 0
    for option in ("--watch", "--interval", "--case", "--run", "--all"):
        assert option in result.output


# ------------------------------------------------------------------ per-run logs; same case ID twice

def test_two_runs_with_the_same_case_id_stay_distinguishable_and_logs_do_not_mix(monkeypatch, tmp_path):
    """The r6 scenario: a duplicate launch is SIGTERM'd while the official run succeeds."""

    from rnaseq.service import _Termination, execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    moments = iter([datetime(2026, 9, 25, 20, 30, 32), datetime(2026, 9, 25, 20, 31, 5)])
    from rnaseq import service
    original_create = service.create_case_run
    monkeypatch.setattr("rnaseq.service.create_case_run", lambda report, case_id: original_create(report, case_id, moment=next(moments)))

    def terminated(command, *, cwd, stdout_path, stderr_path):
        stdout_path.write_text("duplicate nextflow output\n", encoding="utf-8")
        raise _Termination("SIGTERM")

    monkeypatch.setattr("rnaseq.service._run_command", terminated)
    with pytest.raises(_Termination):
        execute_service_run(report, case_id="SAME-CASE")
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory(observed))
    official = execute_service_run(report, case_id="SAME-CASE")

    runs = sorted((tmp_path / "project" / "runs" / "SAME-CASE").iterdir())
    assert [run.name for run in runs] == ["20260925-203032+0800", "20260925-203105+0800"]
    duplicate, official_dir = runs
    assert official_dir == official.run_dir
    assert json.loads((duplicate / "run_state.json").read_text())["status"] == "INTERRUPTED"
    assert json.loads((duplicate / "run_state.json").read_text())["interrupted_by"] == "SIGTERM"
    assert json.loads((official_dir / "run_state.json").read_text())["status"] == "SUCCESS"

    duplicate_log = (duplicate / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    official_log = (official_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "SAME-CASE/20260925-203032+0800 started" in duplicate_log and "state INTERRUPTED" in duplicate_log
    assert "20260925-203105+0800" not in duplicate_log and "state SUCCESS" not in duplicate_log
    assert "SAME-CASE/20260925-203105+0800 started" in official_log and "run finished: SUCCESS" in official_log
    assert "20260925-203032+0800" not in official_log and "INTERRUPTED" not in official_log
    assert (duplicate / "logs" / "upstream.stdout.log").read_text() == "duplicate nextflow output\n"
    assert "duplicate" not in (official_dir / "logs" / "upstream.stdout.log").read_text()
    # Nothing nf-rna-owned is written at the project root.
    assert not list((tmp_path / "project").glob("*.log"))

    latest = _status(tmp_path / "project")
    assert "Run 20260925-203105+0800" in latest and "\nSUCCESS" in latest
    assert "✓ Salmon quantification" in latest and "Tasks\n1 completed / 0 failed" in latest
    older = _status(tmp_path / "project", "--case", "SAME-CASE", "--run", "20260925-203032+0800")
    assert "Run 20260925-203032+0800" in older and "\nINTERRUPTED" in older
    # Delivery still excludes internal logs.
    assert not any(path.name.startswith("rnaseq") for path in (official_dir / "delivery").rglob("*"))


def test_ctrl_c_during_a_run_is_recorded_as_interrupted(monkeypatch, tmp_path):
    report = salmon_case_project(tmp_path, monkeypatch)

    def interrupted(command, *, cwd, stdout_path, stderr_path):
        raise KeyboardInterrupt

    monkeypatch.setattr("rnaseq.service._run_command", interrupted)
    result = runner.invoke(app, ["run", str(tmp_path / "project"), "--case-id", "CASE-CTRL-C", "--yes"])
    assert result.exit_code == 130, result.output
    assert "INTERRUPTED: rnaseq received SIGINT" in result.output
    state_path = next((tmp_path / "project" / "runs" / "CASE-CTRL-C").glob("*/run_state.json"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    assert state["status"] == "INTERRUPTED" and state["interrupted_by"] == "SIGINT" and state["completed_at"]


def test_unexpected_error_is_recorded_as_failed_with_a_traceback_in_the_run_log(monkeypatch, tmp_path):
    from rnaseq.service import execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)

    def broken(command, *, cwd, stdout_path, stderr_path):
        raise KeyError("unexpected")

    monkeypatch.setattr("rnaseq.service._run_command", broken)
    with pytest.raises(KeyError):
        execute_service_run(report, case_id="CASE-BUG")
    run_dir = next((tmp_path / "project" / "runs" / "CASE-BUG").iterdir())
    state = json.loads((run_dir / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "FAILED" and "Unexpected KeyError" in state["error"]
    log = (run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "Traceback" in log and "KeyError: 'unexpected'" in log


def test_successful_run_log_records_lifecycle_and_process_identity(monkeypatch, tmp_path):
    from rnaseq.service import execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory([]))
    run = execute_service_run(report, case_id="CASE-LOG")
    log = (run.run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    for expected in (
        f"run CASE-LOG/{run.run_id} started by pid {os.getpid()}",
        "state RUNNING phase=upstream", "upstream Nextflow launched", "upstream Nextflow exited with code 0",
        "state RUNNING phase=downstream", "downstream Nextflow exited with code 0",
        "state RUNNING phase=delivery", "state SUCCESS phase=delivery", "run finished: SUCCESS",
    ):
        assert expected in log
    record = json.loads((run.run_dir / "logs" / "rnaseq.process.json").read_text(encoding="utf-8"))
    assert record["pid"] == os.getpid() and record["command"][:2] == ["rnaseq", "run"]
    # After the run, status reports SUCCESS from durable state even though this process is alive.
    status = build_status(*select_run(tmp_path / "project"), cache=_Cache())
    assert status.state == "SUCCESS" and status.phase == "completed"
    assert "CPU 14 max\nMemory 54 GiB max" in _norm(render_status(status))


def test_hisat2_route_records_an_upstream_trace_that_status_reads(tmp_path):
    from rnaseq.service import CaseRun, write_upstream_observer_config

    run_dir = _run_dir(tmp_path, case="H2-CASE")
    config = write_upstream_observer_config(CaseRun("H2-CASE", RUN_ID, run_dir, "2026-09-25T20:30:32+08:00"))
    text = config.read_text(encoding="utf-8")
    assert "trace {" in text and str((run_dir / "provenance" / "upstream.trace.txt").resolve()) in text
    assert "report" not in text and "timeline" not in text
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T21:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [("aa/000001", "HISAT2_ALIGN (S1)", "COMPLETED"), ("aa/000002", "FEATURECOUNTS (S1)", "COMPLETED")])
    output = _status(tmp_path)
    assert "✓ Upstream" in output and "Tasks\n2 completed / 0 failed" in output  # no frozen route: generic upstream stage


def test_real_sigterm_during_upstream_records_interrupted_and_stops_nextflow(monkeypatch, tmp_path):
    import signal
    import time

    from rnaseq.service import _Termination, execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    script = tmp_path / "fake-nextflow"
    script.write_text(
        "#!/bin/sh\n"
        f"trap 'echo terminated > {tmp_path}/terminated; exit 143' TERM\n"
        f"echo $$ > {tmp_path}/started\n"
        "while :; do sleep 0.1; done\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.setattr("rnaseq.service.build_nextflow_command", lambda *_args, **_kwargs: [str(script)])
    child = os.fork()
    if child == 0:
        try:
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            execute_service_run(report, case_id="CASE-SIGTERM")
            os._exit(0)
        except _Termination:
            os._exit(143)
        except BaseException:
            os._exit(99)
    deadline = time.monotonic() + 30
    while not (tmp_path / "started").exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    nextflow_pid = int((tmp_path / "started").read_text(encoding="utf-8"))
    os.kill(child, signal.SIGTERM)
    _pid, status = os.waitpid(child, 0)
    assert os.WEXITSTATUS(status) == 143
    assert (tmp_path / "terminated").read_text(encoding="utf-8").strip() == "terminated"
    with pytest.raises(ProcessLookupError):
        os.kill(nextflow_pid, 0)
    run_dir = next((tmp_path / "project" / "runs" / "CASE-SIGTERM").iterdir())
    state = json.loads((run_dir / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "INTERRUPTED" and state["interrupted_by"] == "SIGTERM" and state["phase"] == "upstream"
    assert "upstream Nextflow launched" in (run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    output = _status(tmp_path / "project")
    assert "\nINTERRUPTED" in output and "⏸ Upstream" in output
