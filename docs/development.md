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
`Dockerfile` remains for that historical path only; v1.3.0 does not use it:

```bash
docker build --build-arg NF_RNA_SOURCE_REVISION="$(git rev-parse HEAD)" -t nf-rna:dev .
docker run --rm nf-rna:dev rnaseq --version
```

## Validation

```bash
python -m pip install --upgrade build
python -m pytest -q
python -m build
git diff --check
```

The release gate additionally installs each wheel and source distribution into a fresh environment, checks its bundled Nextflow assets from outside the checkout, and runs CLI/resource smoke checks. See `.github/workflows/release-gate.yml` for the authoritative CI sequence.
