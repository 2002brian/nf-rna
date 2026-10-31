#!/bin/bash
# v1.6.1 reliability acceptance runs (small: D4/D5 FASTQ fixtures from Phase 3, raw-count L2).
R=/home/brian/nf-rna-validation/v161-reliability
echo $$ > $R/runner.pid
cd $R/nfrna
for p in "$@"; do
  echo "$(date -u +%FT%TZ) START $p"
  ( cd $p && $R/env-nfrna/bin/rnaseq run . --case-id ${p^^} --profile local --yes > $R/run_$p.log 2>&1 )
  rc=$?
  st=$(grep -ho '"status": "[A-Z_]*"' $p/runs/${p^^}/*/run_state.json 2>/dev/null | head -1)
  echo "$(date -u +%FT%TZ) END $p rc=$rc $st"
done
echo "$(date -u +%FT%TZ) ALL_DONE"
