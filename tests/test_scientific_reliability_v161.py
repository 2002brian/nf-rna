"""v1.6.1 scientific reliability: strandedness QC (F-1), actual thread provenance (F-2) and recorded
R warnings (F-3).  Expected values come from the Phase 3 validation (validation/phase3/report.md on
branch validation/v1.6.1-phase3, commit 64c6977): featureCounts -s 1 / -s 2 assignments of the
count-only BAMs of the D3, D4 and D5 runs."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from rnaseq.execution import TRACE_FIELDS, render_trace_fields_config, upstream_execution_threads
from rnaseq.hisat2_featurecounts import STRAND_QC_MIN_ASSIGNED, evaluate_strandedness, strandedness_verdict
from rnaseq.workflow_support import _r_warnings_html, _shrinkage_html, _upstream_qc_html

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------------ F-1 strandedness QC

@pytest.mark.parametrize(
    "case, declared, s1, s2, status",
    [
        ("D3 GEUVADIS, unstranded", "unstranded", 431639, 431694, "CONSISTENT"),
        ("D4 yeast, reverse (antisense overlaps: ratio 0.17)", "reverse", 6255, 35996, "CONSISTENT"),
        ("D5-F simulated forward", "forward", 190079, 18204, "CONSISTENT"),
        ("D5-R simulated reverse", "reverse", 18217, 189829, "CONSISTENT"),
        ("D5-SE simulated single-end reverse", "reverse", 10850, 188067, "CONSISTENT"),
        ("D5-NEG: reverse library declared forward", "forward", 18217, 189829, "CONTRADICTS_DECLARATION"),
        ("D4 declared forward", "forward", 6255, 35996, "CONTRADICTS_DECLARATION"),
        ("D5-U unstranded library declared reverse", "reverse", 104192, 103994, "APPEARS_UNSTRANDED"),
        ("D5-R reverse library declared unstranded", "unstranded", 18217, 189829, "APPEARS_STRANDED"),
    ],
)
def test_strandedness_verdict_on_phase3_evidence(case, declared, s1, s2, status):
    verdict = strandedness_verdict(declared, s1, s2)
    assert verdict["status"] == status, case
    assert verdict["declared"] == declared  # the declaration is reported, never replaced
    severity = {"CONSISTENT": "pass", "CONTRADICTS_DECLARATION": "warning", "APPEARS_UNSTRANDED": "warning", "APPEARS_STRANDED": "advisory"}[status]
    assert verdict["severity"] == severity


@pytest.mark.parametrize(
    "s1, s2, observed",
    [
        (STRAND_QC_MIN_ASSIGNED // 2, STRAND_QC_MIN_ASSIGNED // 2 - 1, "insufficient_evidence"),  # below the minimum evidence
        (0, 0, "insufficient_evidence"),
        (2000, 1000, "ambiguous"),   # ratio 2: neither stranded (>= 3x) nor balanced (<= 1.5x)
        (1000, 1800, "ambiguous"),
    ],
)
def test_unclear_evidence_is_inconclusive_not_forced(s1, s2, observed):
    for declared in ("unstranded", "forward", "reverse"):
        verdict = strandedness_verdict(declared, s1, s2)
        assert verdict["observed"] == observed and verdict["status"] == "INCONCLUSIVE" and verdict["severity"] == "inconclusive"


def test_one_sided_evidence_is_handled():
    assert strandedness_verdict("forward", 5000, 0)["status"] == "CONSISTENT"
    assert strandedness_verdict("reverse", 5000, 0)["status"] == "CONTRADICTS_DECLARATION"
    assert strandedness_verdict("forward", 5000, 0)["forward_over_reverse"] == "inf"
    with pytest.raises(ValueError):
        strandedness_verdict("auto", 10, 10)


def _summary(path: Path, assigned: int) -> Path:
    path.write_text(f"Status\t{path.stem}.bam\nAssigned\t{assigned}\nUnassigned_NoFeatures\t7\n", encoding="utf-8")
    return path


def test_run_level_status_is_the_most_severe_sample(tmp_path):
    evidence = {
        "ok": (_summary(tmp_path / "ok_s1.summary", 6255), _summary(tmp_path / "ok_s2.summary", 35996)),
        "few": (_summary(tmp_path / "few_s1.summary", 10), _summary(tmp_path / "few_s2.summary", 50)),
    }
    assert evaluate_strandedness("reverse", evidence)["status"] == "INCONCLUSIVE"
    evidence["wrong"] = (_summary(tmp_path / "w_s1.summary", 189829), _summary(tmp_path / "w_s2.summary", 18217))
    result = evaluate_strandedness("reverse", evidence)
    assert result["status"] == "WARNING" and result["samples"]["wrong"]["status"] == "CONTRADICTS_DECLARATION"
    assert result["policy"]["action"].startswith("report only")


def _process(text: str, name: str) -> str:
    match = re.search(rf"process {name} \{{(.*?)\n\}}\n", text, re.S)
    assert match, name
    return match.group(1)


def test_strand_evidence_reuses_the_count_only_bam_and_leaves_counting_unchanged():
    text = (ROOT / "workflow" / "hisat2_featurecounts.nf").read_text(encoding="utf-8")
    counting = _process(text, "FEATURECOUNTS")
    evidence = _process(text, "STRAND_EVIDENCE")
    # The counts used downstream are still made with the declared strandedness only.
    assert "-s ${strand}" in counting and "featureCountsStrand(params.strandedness)" in counting
    main = re.search(r"featureCounts (.*)", counting).group(1)
    for strand in ("1", "2"):
        line = next(item for item in re.findall(r"featureCounts (.*)", evidence) if f"-s {strand} " in item)
        # Same counting options as FEATURECOUNTS apart from the strand and output name.
        normalise = lambda command: re.sub(r"-o \S+|-s \S+", "", command).split()
        assert normalise(line) == normalise(main)
    workflow = text.split("workflow {", 1)[1]
    assert "STRAND_EVIDENCE(PREPARE_COUNT_BAM.out.countable" in workflow
    assert "STRAND_EVIDENCE.out" not in workflow  # never feeds counting, assembly or MultiQC


# ------------------------------------------------------------------ F-2 thread provenance

def test_trace_fields_add_cpus_without_changing_scheduling():
    config = render_trace_fields_config()
    assert config.strip() == f"trace.fields = '{TRACE_FIELDS}'"
    fields = TRACE_FIELDS.split(",")
    assert {"cpus", "process", "tag"} <= set(fields)
    # Nextflow's default columns stay (run_status parses traces by header).
    assert {"task_id", "hash", "native_id", "name", "status", "exit", "submit", "duration", "realtime"} <= set(fields)
    assert "resourceLimits" not in config and "cpus =" not in config


def _trace(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(fields) + "\n" + "".join("\t".join(row.get(name, "") for name in fields) + "\n" for row in rows), encoding="utf-8")


def test_hisat2_route_threads_come_from_the_runs_own_trace(tmp_path):
    fields = TRACE_FIELDS.split(",")
    rows = [
        {"process": "HISAT2_ALIGN", "name": "HISAT2_ALIGN (A)", "status": "COMPLETED", "cpus": "4"},
        {"process": "HISAT2_ALIGN", "name": "HISAT2_ALIGN (B)", "status": "COMPLETED", "cpus": "2"},  # e.g. lowered by resourceLimits
        {"process": "FEATURECOUNTS", "name": "FEATURECOUNTS (A)", "status": "CACHED", "cpus": "4"},
        {"process": "STRAND_EVIDENCE", "name": "STRAND_EVIDENCE (A)", "status": "COMPLETED", "cpus": "2"},
        {"process": "FASTQC_RAW", "name": "FASTQC_RAW (A)", "status": "COMPLETED", "cpus": "2"},
        {"process": "HISAT2_ALIGN", "name": "HISAT2_ALIGN (C)", "status": "FAILED", "cpus": "8"},
    ]
    _trace(tmp_path / "provenance" / "upstream.trace.txt", rows, fields)
    record = upstream_execution_threads(tmp_path, "hisat2_featurecounts")
    assert record["status"] == "RECORDED" and record["trace"] == "provenance/upstream.trace.txt"
    assert record["processes"]["HISAT2_ALIGN"] == {"process": "HISAT2_ALIGN", "tasks": 2, "cpus": [2, 4]}
    assert record["processes"]["FEATURECOUNTS"]["cpus"] == [4]
    assert "FASTQC_RAW" not in record["processes"]
    assert record["tool_versions"]["HISAT2_ALIGN"] == {"hisat2": "2.2.3"}
    assert "thread count" in record["note"]


def test_salmon_route_threads_and_versions_from_nfcore_records(tmp_path):
    fields = TRACE_FIELDS.split(",")
    info = tmp_path / "upstream" / "nfcore_rnaseq" / "pipeline_info"
    _trace(info / "execution_trace_2026-10-31_00-00-00.txt", [
        {"process": "NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT", "name": "NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (A)", "status": "COMPLETED", "cpus": "6"},
        {"process": "NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE", "name": "x:TRIMGALORE (A)", "status": "COMPLETED", "cpus": "4"},
    ], fields)
    (info / "nf_core_rnaseq_software_mqc_versions.yml").write_text("SALMON_QUANT:\n  salmon: 1.10.3\nTRIMGALORE:\n  trimgalore: 2.1.0\nMULTIQC:\n  multiqc: 1.33\n", encoding="utf-8")
    cmd = tmp_path / "upstream" / "nfcore_rnaseq" / "salmon" / "A" / "cmd_info.json"
    cmd.parent.mkdir(parents=True)
    cmd.write_text(json.dumps({"salmon_version": "1.10.3", "threads": "6", "libType": "ISR"}), encoding="utf-8")
    record = upstream_execution_threads(tmp_path, "salmon")
    assert record["status"] == "RECORDED"
    assert record["processes"]["SALMON_QUANT"]["cpus"] == [6] and record["processes"]["TRIMGALORE"]["cpus"] == [4]
    assert record["salmon_quant_threads_cmd_info"] == {"A": 6}
    assert record["tool_versions"] == {"SALMON_QUANT": {"salmon": "1.10.3"}, "TRIMGALORE": {"trimgalore": "2.1.0"}}


def test_runs_without_a_cpus_column_are_not_given_a_thread_count(tmp_path):
    _trace(tmp_path / "provenance" / "upstream.trace.txt", [{"name": "HISAT2_ALIGN (A)", "status": "COMPLETED"}],
           ["task_id", "hash", "native_id", "name", "status", "exit"])
    record = upstream_execution_threads(tmp_path, "hisat2_featurecounts")
    assert record["status"] == "NOT_RECORDED" and "no cpus column" in record["reason"] and record["processes"] == {}
    assert upstream_execution_threads(tmp_path / "other", "hisat2_featurecounts")["status"] == "NOT_RECORDED"


# ------------------------------------------------------------------ report rendering

def test_report_shows_strandedness_warning_and_threads():
    contract = {"upstream_qc": {
        "strandedness": {"status": "WARNING", "declared": "forward", "policy": {"stranded_if_ratio_at_least": 3.0, "unstranded_if_ratio_within": [0.6667, 1.5], "minimum_strand_assigned": 1000},
                         "samples": {"D5_R_A": strandedness_verdict("forward", 18217, 189829)}},
        "execution_threads": {"status": "RECORDED", "processes": {"HISAT2_ALIGN": {"cpus": [4], "tasks": 2}}, "note": "thread count note"},
    }}
    html = "".join(_upstream_qc_html(contract))
    assert "WARNING" in html and "D5_R_A" in html and "CONTRADICTS_DECLARATION" in html and "were not changed" in html
    assert "HISAT2_ALIGN" in html and "<td>4</td>" in html and "thread count note" in html
    unstated = "".join(_upstream_qc_html({"upstream_qc": {"strandedness": {"status": "NOT_EVALUATED", "reason": "Salmon route reason", "salmon_library_types": {"A": "ISR"}},
                                                        "execution_threads": {"status": "NOT_RECORDED", "reason": "no trace"}}}))
    assert "Not evaluated: Salmon route reason" in unstated and "ISR" in unstated and "Not recorded: no trace" in unstated
    assert _upstream_qc_html({"source": {"type": "raw_counts"}}) == []


def test_report_lists_r_warnings_verbatim_and_keeps_failure_separate():
    warning = {"step": "shrinkage: lfcShrink(type = apeglm)", "message": "the line search routine failed, unable to sufficiently decrease the function value",
               "call": "nbinomGLM(x = x, Y = YNZ)", "count": 6}
    text = _r_warnings_html({"status": "WARNINGS", "fit": [], "contrast": [warning]})
    assert "the line search routine failed" in text and "&times; 6" in text and "not errors" in text
    assert "none" in _r_warnings_html({"status": "NONE", "fit": [], "contrast": []})
    assert "before nf-rna v1.6.1" in _r_warnings_html(None)
    failed = _shrinkage_html({"deg_definition": {"rule": "r"}, "shrinkage": {"method": "apeglm", "status": "failed", "error": "boom"}, "r_warnings": {"fit": [], "contrast": [warning]}})
    assert "<strong>failed</strong>" in failed and "boom" in failed and "line search" in failed


def test_l2_script_records_warnings_without_muffling_them():
    script = (ROOT / "src" / "rnaseq" / "r" / "l2_analysis.R").read_text(encoding="utf-8")
    assert "withCallingHandlers" in script and "muffleWarning" not in script
    for step in ("fit: DESeq2::DESeq()", "results()", "shrinkage: lfcShrink(type = apeglm)", "shrinkage: lfcShrink(type = ashr)"):
        assert f'record_warnings("{step}"' in script
    # An error keeps its existing handling: shrinkage failure stays inside shrink_lfc's tryCatch.
    assert 'error = function(e) { info$error <<- conditionMessage(e); NULL }' in script
