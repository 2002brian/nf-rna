#!/bin/bash
W=/home/brian/nf-rna-validation/v161
cd $W/nfrna
for spec in D1:FIX-D1 D1:FIX-D1-REP D1s42:FIX-D1-S42 D1a:FIX-D1a D1b:FIX-D1b D2:FIX-D2; do
  d=${spec%%:*}; c=${spec##*:}
  echo "$(date -u +%FT%TZ) START $d $c"
  $W/env-nfrna/bin/rnaseq run $d --case-id $c --profile local --yes > $W/logs/run_$c.log 2>&1
  rc=$?
  st=$(grep -ho '"status": "[A-Z_]*"' $d/runs/$c/*/run_state.json 2>/dev/null | head -1)
  echo "$(date -u +%FT%TZ) END $d $c rc=$rc $st"
done
echo "$(date -u +%FT%TZ) ALL_DONE"
