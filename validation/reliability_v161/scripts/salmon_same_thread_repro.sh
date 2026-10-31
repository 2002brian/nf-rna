#!/usr/bin/env bash
# Explains the rel_d4_sal quant.sf difference: Salmon 1.10.3, identical trimmed reads, index and
# arguments, --threads 6 (as nf-core used), run three times.  NumReads never differs;
# EffectiveLength can, through the fragment-length distribution estimate.
set -euo pipefail
OUT=${1:?output tsv}
S=$(mktemp -d)
SAL=/home/brian/miniforge3/envs/nf-rna-reference-builder/bin/salmon
T=/home/brian/nf-rna-validation/v161-phase3/reference/scratch/D4/trim/WT_REP2
IDX=/home/brian/nf-rna-validation/v161-phase3/data/yeast_ref/salmon_build/salmon_index
for r in a b c; do
  "$SAL" quant --libType=ISR --index "$IDX" -1 "$T/WT_REP2_R1_val_1.fq.gz" -2 "$T/WT_REP2_R2_val_2.fq.gz" --threads 6 -o "$S/$r" > /dev/null 2>&1
done
echo -e "pair\teffective_length_rows_differing\tnumreads_rows_differing\tfrag_length_mean_1\tfrag_length_mean_2" > "$OUT"
for p in "a b" "a c" "b c"; do
  set -- $p
  d=$(paste <(tail -n+2 "$S/$1/quant.sf") <(tail -n+2 "$S/$2/quant.sf") | awk '{if($3!=$8)e++; if($5!=$10)n++} END{print e+0"\t"n+0}')
  m1=$(python3 -c "import json;print(json.load(open('$S/$1/aux_info/meta_info.json'))['frag_length_mean'])")
  m2=$(python3 -c "import json;print(json.load(open('$S/$2/aux_info/meta_info.json'))['frag_length_mean'])")
  echo -e "$1 vs $2\t$d\t$m1\t$m2" >> "$OUT"
done
rm -rf "$S"
