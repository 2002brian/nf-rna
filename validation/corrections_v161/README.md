# Evidence for the v1.6.1 scientific corrections

The corrections are in commit `3cf0386` on `fix/v1.6.1-scientific`. They are described in [`docs/dev/scientific_corrections_v161.md`](../../docs/dev/scientific_corrections_v161.md).

The frozen Phase 2 evidence is on branch `validation/v1.6.0` (commit `b886f11`) and is not modified here.

| Path | Content |
|---|---|
| `scripts/d1_isolation.R`, `scripts/run_d1_isolation.sh` | D-1 alone: recompute the v1.6.0 GSEA exactly, apply the correction, compare with v1.6.0 and the reference |
| `scripts/run_nfrna.sh` | The six corrected nf-rna runs (FIX-D1, FIX-D1-REP, FIX-D1-S42, FIX-D1a, FIX-D1b, FIX-D2) |
| `scripts/compare_v161.py` | Before/after and reference comparison (`python -I compare_v161.py <v1.6.0 root> <v161 root> <outdir>`) |
| `results/` | `d1_isolation.tsv`, `m1_deseq2_alpha.tsv`, `m1_ora.tsv`, `gsea_fixed_vs_reference_and_v160.tsv`, `reproducibility.tsv`, `seed_1_vs_42.tsv` |
| `evidence/` | Run identities (`fixed_runs.tsv`), install record (`direct_url.json`), the seed-42 project file and its diff from D1 |

Large data stays outside Git, in `/home/brian/nf-rna-validation/`:

- `v1.6.0/`: the Phase 2 runs and the reference.
- `v161/`: the corrected install (`env-nfrna`), the projects and runs, logs, and the D-1 isolation work directory.

Corrected runs used the locked downstream lock `0e721ff3…` and were installed from commit `3cf0386`.
