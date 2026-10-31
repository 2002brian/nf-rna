#!/bin/bash
set -u
REPO=/home/brian/nf-rna-validation/v161/d1_isolation/repo_view; mkdir -p $REPO/src/rnaseq; ln -sfn /home/brian/nf-rna-validation/v161/d1_isolation/r_snapshot $REPO/src/rnaseq/r
S=/home/brian/projects/nf-rna/validation/corrections_v161/scripts/d1_isolation.R
cp $S /home/brian/nf-rna-validation/v161/d1_isolation/d1_isolation.R
L1=/home/brian/nf-rna-validation/v1.6.0/nfrna/D1/runs/VAL-D1/20261031-003256+0800/downstream/l2
L2=/home/brian/nf-rna-validation/v1.6.0/nfrna/D2/runs/VAL-D2/20261031-005023+0800/downstream/l2
/home/brian/.cache/nf-rna/runtime/prefixes/nf-rna-downstream-linux-64-0e721ff3cf1983f6-b0a7fe78cf4787e5/bin/Rscript /home/brian/nf-rna-validation/v161/d1_isolation/d1_isolation.R $REPO $L1 trt_vs_untrt org.Hs.eg.db ~/.cache/nf-rna/genesets/kegg/hsa/hsa_kegg20261028_20261030T061908Z /home/brian/nf-rna-validation/v1.6.0/reference/outputs/D1/contrasts/trt_vs_untrt /home/brian/nf-rna-validation/v161/d1_isolation/D1_trt_vs_untrt.tsv
for c in LP_vs_Basal ML_vs_Basal LP_vs_ML; do
/home/brian/.cache/nf-rna/runtime/prefixes/nf-rna-downstream-linux-64-0e721ff3cf1983f6-b0a7fe78cf4787e5/bin/Rscript /home/brian/nf-rna-validation/v161/d1_isolation/d1_isolation.R $REPO $L2 $c org.Mm.eg.db ~/.cache/nf-rna/genesets/kegg/mmu/mmu_kegg20261028_20261030T065742Z /home/brian/nf-rna-validation/v1.6.0/reference/outputs/D2/contrasts/$c /home/brian/nf-rna-validation/v161/d1_isolation/D2_$c.tsv
done
echo ALL_DONE
