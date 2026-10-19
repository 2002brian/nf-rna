# Development installation

This document is for contributors. It is not required to install or run the official nf-rna release.

## Source checkout and editable environment

```bash
git clone https://github.com/2002brian/nf-rna.git
cd nf-rna
conda env create -f environment.yml
conda activate nf-rna
python -m pip install -e '.[dev]'
```

The shared Conda environment provides the development R stack for real-R tests. Contributors use editable installation so source changes are immediately visible; production users should use the non-editable, release-pinned installation in the [README](../README.md).

## Development runs

A development run executes from the source checkout and records its `HEAD` as
the source revision. nf-rna refuses to start a run while the checkout has
uncommitted changes, so commit (or stash) work before an end-to-end run. The
downstream wheel for a source checkout is built with `python -m build` from the
`dev` extra.

## Historical Docker image build (v1.2.1 and earlier)

Releases up to v1.2.1 ran downstream tasks in a first-party Docker image. The
`Dockerfile`, `environment.docker.yml`, and `.dockerignore` were removed from the
current tree after v1.3.0; nf-rna no longer builds, publishes, or runs a Docker
image. To reproduce the historical image, build it from the immutable `v1.2.1`
tag as documented in that release.

## Validation

```bash
python -m pip install --upgrade build
python -m pytest -q
python -m build
git diff --check
```

The release gate additionally installs each wheel and source distribution into a fresh environment, checks its bundled Nextflow assets from outside the checkout, and runs CLI/resource smoke checks. See `.github/workflows/release-gate.yml` for the authoritative CI sequence.

### R-dependent tests

Tests that run the R analysis scripts (L1/L2 backends, figures, enrichment) call `require_rscript()` or `require_r_packages()` and are **skipped** when `Rscript` or a required R package is missing. A plain `python -m pytest -q` therefore still passes without the R stack; CI installs base R plus the figure packages (ggplot2, pheatmap, jsonlite, colorspace, systemfonts, ggrepel) at the downstream lock's versions, so the figure tests run there; tests that also need Bioconductor packages (DESeq2, tximport, clusterProfiler, OrgDb) are skipped. `-rs` lists every skip and its reason.

For release validation, run these tests with the same locked R runtime the pipeline uses: the downstream prefix that `rnaseq run` provisions from `workflow/envs/locks/nf-rna-downstream-<platform>.lock.yml`. It lives under `${RNASEQ_RUNTIME_ROOT:-${XDG_CACHE_HOME:-$HOME/.cache}/nf-rna/runtime}/prefixes/` and is named `nf-rna-downstream-<platform>-<first 16 hex of the lock SHA-256>-<first 16 hex of the wheel SHA-256>`; the lock SHA-256 is in `workflow/envs/locks/SHA256SUMS`. Put its `bin` directory on `PATH` after the development environment, so pytest keeps the development Python and resolves `Rscript` to the locked R, and set `NF_RNA_REQUIRE_R=1` so a missing runtime or package fails the test instead of skipping it:

```bash
LOCK=$(awk '/linux-64/ {print substr($1, 1, 16)}' workflow/envs/locks/SHA256SUMS)
R_PREFIX=$(ls -d "${XDG_CACHE_HOME:-$HOME/.cache}"/nf-rna/runtime/prefixes/nf-rna-downstream-linux-64-"$LOCK"-* | tail -n 1)
PATH="$CONDA_PREFIX/bin:$R_PREFIX/bin:$PATH" NF_RNA_REQUIRE_R=1 python -m pytest -q -rs
```

The prefix exists only after a downstream run (for example the `examples/nfcore_smoke_test` smoke run) has provisioned it. Check that `Rscript` resolves into `$R_PREFIX` before trusting a pass. With `NF_RNA_REQUIRE_R=1`, a pass with no R-related skips means every R-dependent test actually ran against that runtime.
