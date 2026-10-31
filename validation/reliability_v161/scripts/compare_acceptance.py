"""v1.6.1 reliability acceptance (F-1, F-2, F-3) on small runs from nf-rna 9c37d84.

Compares every new run with the run of the same input made before the change (Phase 3 runs from
76c0568; FIX-D1 from 3cf0386) and checks the new QC/provenance records, the report, the run log
and the GUI data layer.  Writes one row per check.
Usage: python -I compare_acceptance.py <reliability root> <phase3 root> <v161 root> <featureCounts> <samtools> <out.tsv>
"""
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

REL, P3, V161, FC, SAMTOOLS, OUT = (Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], sys.argv[5], Path(sys.argv[6]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from rnaseq.gui import results  # noqa: E402
from rnaseq.hisat2_featurecounts import featurecounts_assigned, strandedness_verdict  # noqa: E402

rows = []


def add(check, item, observed, expected, ok):
    rows.append([check, item, str(observed), str(expected), "PASS" if ok else "FAIL"])


def run_dir(root: Path, project: str, case: str = "*") -> Path:
    found = sorted((root / project / "runs").glob(f"{case}/*/run_state.json"))
    assert len(found) == 1, (project, found)
    return found[0].parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_tables(new: Path, old: Path, pattern: str) -> tuple[int, list[str]]:
    files = sorted(p.relative_to(old) for p in old.glob(pattern) if p.is_file())
    differ = [str(f) for f in files if not (new / f).is_file() or sha(new / f) != sha(old / f)]
    return len(files), differ


# ---------------------------------------------------------------- F-1 on the HISAT2 route
for project, old_project, declared, expect_status, expect_sample in (
    ("rel_d5_neg", "p3_d5_neg", "forward", "WARNING", "CONTRADICTS_DECLARATION"),
    ("rel_d4_hf", "p3_d4_hf", "reverse", "CONSISTENT", "CONSISTENT"),
):
    new, old = run_dir(REL / "nfrna", project), run_dir(P3 / "nfrna", old_project)
    add("EXEC", project, json.loads((new / "run_state.json").read_text())["status"], "SUCCESS", json.loads((new / "run_state.json").read_text())["status"] == "SUCCESS")
    qc = json.loads((new / "upstream/hisat2_featurecounts/qc/strandedness/strandedness_qc.json").read_text())
    add("F-1", f"{project}: run-level strandedness status", qc["status"], expect_status, qc["status"] == expect_status)
    for sample, item in qc["samples"].items():
        add("F-1", f"{project}/{sample}: status (ratio -s1/-s2 {item['forward_over_reverse']:.3g})", item["status"], expect_sample, item["status"] == expect_sample)
        add("F-1", f"{project}/{sample}: declared kept", item["declared"], declared, item["declared"] == declared)
    contract = json.loads((new / "frozen/downstream_contract.json").read_text())
    add("F-1", f"{project}: frozen contract upstream_qc.strandedness", contract["upstream_qc"]["strandedness"]["status"], expect_status, contract["upstream_qc"]["strandedness"]["status"] == expect_status)
    log = (new / "logs" / "rnaseq.log").read_text(errors="replace") if (new / "logs" / "rnaseq.log").is_file() else ""
    logged = sum("strandedness QC" in line for line in log.splitlines())
    add("F-1", f"{project}: run log warnings", logged, len(qc["samples"]) if expect_status == "WARNING" else 0, logged == (len(qc["samples"]) if expect_status == "WARNING" else 0))
    report = next(new.glob("delivery/report*.html"), new / "downstream/report/report.html").read_text()
    add("F-1", f"{project}: report shows strandedness check", "Strandedness check" in report and expect_status in report, True, "Strandedness check" in report and expect_status in report)
    gui = results.qc_view(results.result_runs(REL / "nfrna" / project)[0]).strandedness
    add("F-1", f"{project}: GUI qc_view strandedness", gui.get("status"), expect_status, gui.get("status") == expect_status)
    n, differ = same_tables(new / "delivery/counts", old / "delivery/counts", "raw_counts.csv")
    add("UNCHANGED", f"{project}: raw_counts.csv vs {old_project}", f"{n - len(differ)}/{n} identical", "all", n and not differ)
    n, differ = same_tables(new / "upstream/hisat2_featurecounts/counts", old / "upstream/hisat2_featurecounts/counts", "**/*.counts.txt.summary")
    add("UNCHANGED", f"{project}: featureCounts summaries vs {old_project}", f"{n - len(differ)}/{n} identical", "all", n and not differ)
    n, differ = same_tables(new / "downstream/l1", old / "downstream/l1", "*.tsv")
    add("UNCHANGED", f"{project}: L1 tables vs {old_project}", f"{n - len(differ)}/{n} identical {differ}", "all", n and not differ)
    threads = yaml.safe_load((new / "provenance/run_provenance.yaml").read_text())["execution_threads"]
    hisat = threads["processes"].get("HISAT2_ALIGN", {})
    add("F-2", f"{project}: HISAT2_ALIGN cpus from trace", hisat.get("cpus"), [4], threads["status"] == "RECORDED" and hisat.get("cpus") == [4])
    add("F-2", f"{project}: FEATURECOUNTS / STRAND_EVIDENCE cpus", (threads["processes"]["FEATURECOUNTS"]["cpus"], threads["processes"]["STRAND_EVIDENCE"]["cpus"]), ([4], [2]),
        threads["processes"]["FEATURECOUNTS"]["cpus"] == [4] and threads["processes"]["STRAND_EVIDENCE"]["cpus"] == [2])
    trace = (new / "provenance/upstream.trace.txt").read_text().splitlines()
    cmd_cpus = {}
    work = Path(yaml.safe_load((new / "provenance/run_provenance.yaml").read_text())["execution_work_dir"]) / "upstream"
    for sh in work.glob("*/*/.command.sh"):
        for match in re.finditer(r"^\s*hisat2 -p (\d+)", sh.read_text(), re.M):
            cmd_cpus.setdefault("hisat2", set()).add(match.group(1))
    add("F-2", f"{project}: HISAT2 -p in executed .command.sh equals trace cpus", sorted(cmd_cpus.get("hisat2", [])), ["4"], sorted(cmd_cpus.get("hisat2", [])) == ["4"])
    add("F-2", f"{project}: trace header has cpus", "cpus" in trace[0].split("\t"), True, "cpus" in trace[0].split("\t"))
    gui_rows = dict(results.provenance(results.result_runs(REL / "nfrna" / project)[0]))
    add("F-2", f"{project}: GUI provenance thread rows", gui_rows.get("threads: HISAT2_ALIGN"), "4 CPU(s) per task, 2 task(s)", gui_rows.get("threads: HISAT2_ALIGN") == "4 CPU(s) per task, 2 task(s)")
    add("F-2", f"{project}: report thread table", "Upstream execution threads" in report and "HISAT2_ALIGN" in report, True, "Upstream execution threads" in report and "HISAT2_ALIGN" in report)

# ---------------------------------------------------------------- Salmon route (F-1 not evaluated; F-2 recorded)
new, old = run_dir(REL / "nfrna", "rel_d4_sal"), run_dir(P3 / "nfrna", "p3_d4_sal")
add("EXEC", "rel_d4_sal", json.loads((new / "run_state.json").read_text())["status"], "SUCCESS", json.loads((new / "run_state.json").read_text())["status"] == "SUCCESS")
contract = json.loads((new / "frozen/downstream_contract.json").read_text())
strand = contract["upstream_qc"]["strandedness"]
add("F-1", "rel_d4_sal: Salmon route strandedness", (strand["status"], strand.get("salmon_library_types")), ("NOT_EVALUATED", {"WT_REP2": "ISR", "RAP1_IAA_30M_REP1": "ISR"}),
    strand["status"] == "NOT_EVALUATED" and strand.get("salmon_library_types") == {"WT_REP2": "ISR", "RAP1_IAA_30M_REP1": "ISR"})
threads = yaml.safe_load((new / "provenance/run_provenance.yaml").read_text())["execution_threads"]
salmon_cpus = threads["processes"].get("SALMON_QUANT", {}).get("cpus")
cmd = threads.get("salmon_quant_threads_cmd_info")
add("F-2", "rel_d4_sal: SALMON_QUANT cpus (nf-core trace) and --threads (cmd_info.json)", (salmon_cpus, cmd), "trace cpus == cmd_info threads",
    threads["status"] == "RECORDED" and bool(cmd) and set(cmd.values()) <= set(salmon_cpus or []))
add("F-2", "rel_d4_sal: tool versions per process", threads.get("tool_versions"), "SALMON_QUANT salmon 1.10.3",
    (threads.get("tool_versions") or {}).get("SALMON_QUANT", {}).get("salmon") == "1.10.3")
n, differ = same_tables(new / "delivery/counts", old / "delivery/counts", "estimated_counts.csv")
add("UNCHANGED", "rel_d4_sal: estimated_counts.csv vs p3_d4_sal", f"{n - len(differ)}/{n} identical", "all", n and not differ)
n, differ = same_tables(new / "upstream/nfcore_rnaseq/salmon", old / "upstream/nfcore_rnaseq/salmon", "*/quant.sf")
add("UNCHANGED", "rel_d4_sal: quant.sf vs p3_d4_sal", f"{n - len(differ)}/{n} identical", "all", n and not differ)

# ---------------------------------------------------------------- F-3 on D3 counts (L2, apeglm warnings expected)
new = run_dir(REL / "nfrna", "rel_d3_counts")
old = run_dir(P3 / "nfrna", "p3_d3_hf")
add("EXEC", "rel_d3_counts", json.loads((new / "run_state.json").read_text())["status"], "SUCCESS", json.loads((new / "run_state.json").read_text())["status"] == "SUCCESS")
meta = json.loads((new / "downstream/l2/contrasts/YRI_vs_GBR/contrast_metadata.json").read_text())
warnings = [*meta["r_warnings"]["fit"], *meta["r_warnings"]["contrast"]]
line_search = sum(item["count"] for item in warnings if "line search" in item["message"])
add("F-3", "rel_d3_counts: apeglm line-search warnings recorded (Phase 3: 6 in the nf-rna work-dir stderr)", line_search, 6, line_search == 6)
add("F-3", "rel_d3_counts: warning status / shrinkage status", (meta["r_warnings"]["status"], meta["shrinkage"]["status"]), ("WARNINGS", "success"),
    meta["r_warnings"]["status"] == "WARNINGS" and meta["shrinkage"]["status"] == "success")
prov = json.loads((new / "downstream/l2/scientific_provenance.json").read_text())
add("F-3", "rel_d3_counts: L2 scientific_provenance.json has r_warnings", "r_warnings" in json.dumps(prov), True, "r_warnings" in json.dumps(prov))
delivered = list(new.glob("delivery/**/contrast_metadata*.json"))
add("F-3", "rel_d3_counts: contrast_metadata.json delivered with r_warnings", [str(p.relative_to(new)) for p in delivered],
    "in delivery", any("r_warnings" in p.read_text() for p in delivered))
report = next(new.glob("delivery/report*.html")).read_text()
add("F-3", "rel_d3_counts: report lists the warning verbatim", "the line search routine failed" in report, True, "the line search routine failed" in report)
deg = results.deg_view(results.result_runs(REL / "nfrna" / "rel_d3_counts")[0], "YRI_vs_GBR")
add("F-3", "rel_d3_counts: GUI deg_view r_warnings", len(deg.r_warnings or []), ">= 1", bool(deg.r_warnings))
for name in ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv"):
    a, b = new / "downstream/l2/contrasts/YRI_vs_GBR" / name, old / "downstream/l2/contrasts/YRI_vs_GBR" / name
    add("UNCHANGED", f"rel_d3_counts: {name} vs p3_d3_hf (same counts)", sha(a) == sha(b), True, sha(a) == sha(b))
add("UNCHANGED", "rel_d3_counts: L2 normalized_counts.tsv vs p3_d3_hf", sha(new / "downstream/l2/normalized_counts.tsv") == sha(old / "downstream/l2/normalized_counts.tsv"), True,
    sha(new / "downstream/l2/normalized_counts.tsv") == sha(old / "downstream/l2/normalized_counts.tsv"))

# ---------------------------------------------------------------- D1 with GO/KEGG/GSEA: every table vs FIX-D1 (3cf0386)
new, old = run_dir(REL / "nfrna", "rel_d1"), run_dir(V161 / "nfrna", "D1", "FIX-D1")
add("EXEC", "rel_d1", json.loads((new / "run_state.json").read_text())["status"], "SUCCESS", json.loads((new / "run_state.json").read_text())["status"] == "SUCCESS")
for part in ("l1", "l2"):
    n, differ = same_tables(new / "downstream" / part, old / "downstream" / part, "**/*.tsv")
    add("UNCHANGED", f"rel_d1: every downstream/{part} TSV (DE, ORA, GSEA) vs FIX-D1", f"{n - len(differ)}/{n} identical", "all", n > 0 and not differ)
    if differ:
        rows[-1][2] += " differ: " + ", ".join(differ[:10])
for path in sorted((new / "downstream/l2/contrasts").glob("*/contrast_metadata.json")):
    meta = json.loads(path.read_text())
    add("F-3", f"rel_d1: {path.parent.name} r_warnings status", meta["r_warnings"]["status"], "recorded", meta["r_warnings"]["status"] in {"NONE", "WARNINGS"})

# ---------------------------------------------------------------- F-1 INCONCLUSIVE on real reads (subsampled count-only BAM)
bam = run_dir(P3 / "nfrna", "p3_d5_hf_r") / "upstream/hisat2_featurecounts/bam/count_only/D5_R_A.countable.bam"
scratch = REL / "compare" / "inconclusive"
scratch.mkdir(parents=True, exist_ok=True)
sub = scratch / "D5_R_A.sub.bam"
subprocess.run([SAMTOOLS, "view", "-b", "-s", "7.002", "-o", str(sub), str(bam)], check=True)
gtf = P3 / "data/yeast_ref/genes.gtf.gz"
for s in (1, 2):
    subprocess.run([FC, "-T", "2", "-a", str(gtf), "-o", str(scratch / f"s{s}.txt"), "-t", "exon", "-g", "gene_id", "-s", str(s), "-Q", "0", "--primary", "-p", "--countReadPairs", "-B", "-C", str(sub)],
                   check=True, capture_output=True)
a1, a2 = featurecounts_assigned(scratch / "s1.txt.summary"), featurecounts_assigned(scratch / "s2.txt.summary")
verdict = strandedness_verdict("forward", a1, a2)
add("F-1", f"D5_R_A subsampled to 0.2% (-s1 {a1}, -s2 {a2}) declared forward", verdict["status"], "INCONCLUSIVE (insufficient evidence)", verdict["status"] == "INCONCLUSIVE")

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("check\titem\tobserved\texpected\tverdict\n")
    for r in rows:
        fh.write("\t".join(x.replace("\t", " ") for x in r) + "\n")
for r in rows:
    print(" | ".join(r))
print("FAIL rows:", sum(r[-1] == "FAIL" for r in rows), "of", len(rows))
