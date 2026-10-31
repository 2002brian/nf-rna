#!/usr/bin/env bash
# U-1 investigation: reproduce the GO_BP / GO_MF differences of D2 LP_vs_ML and locate the
# permutations involved.  Work files go to /home/brian/nf-rna-validation/v161/u1.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
P=/home/brian/.cache/nf-rna/runtime/prefixes/nf-rna-downstream-linux-64-0e721ff3cf1983f6-b0a7fe78cf4787e5
OUT="$HERE/results"; mkdir -p "$OUT"
cd /home/brian/nf-rna-validation/v161/u1
filter() { grep -v "ties in the preranked\|order of those tied\|^Warning\|^In \|preparePathwaysAndStats" || true; }
"$P/bin/Rscript" "$HERE/u1_repro.R" GO_BP BP 2>&1 | filter > "$OUT/repro_GO_BP.txt"
"$P/bin/Rscript" "$HERE/u1_repro.R" GO_MF MF 2>&1 | filter > "$OUT/repro_GO_MF.txt"
"$P/bin/Rscript" "$HERE/u1_counts.R" inputs_GO_BP.rds 2>&1 | filter > "$OUT/counts_GO_BP.txt"
"$P/bin/Rscript" "$HERE/u1_counts.R" inputs_GO_MF.rds 2>&1 | filter > "$OUT/counts_GO_MF.txt"
"$P/bin/Rscript" "$HERE/u1_tie.R" inputs_GO_BP.rds 2>&1 | filter > "$OUT/tie_GO_BP.txt"
"$P/bin/Rscript" "$HERE/u1_tie_es.R" inputs_GO_MF.rds 2>&1 | filter > "$OUT/tie_GO_MF.txt"
/home/brian/miniforge3/envs/nf-rna/bin/python -I "$HERE/u1_mc.py" > "$OUT/monte_carlo_error.txt"
