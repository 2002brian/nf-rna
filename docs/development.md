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

The shared Conda environment provides the development and image-build R stack. Contributors use editable installation so source changes are immediately visible; production users should use the non-editable, release-pinned installation in the [README](../README.md).

## Local execution-image build

Only build locally when developing or validating the image itself:

```bash
docker build --build-arg NF_RNA_SOURCE_REVISION="$(git rev-parse HEAD)" -t nf-rna:dev .
docker run --rm nf-rna:dev rnaseq --version
```

Set `runtime.execution_image: nf-rna:dev` only in a development project. It is not an official image and is not valid for production acceptance.

## Validation

```bash
python -m pip install --upgrade build
python -m pytest -q
python -m build
git diff --check
```

The release gate additionally installs each wheel and source distribution into a fresh environment, checks its bundled Nextflow assets from outside the checkout, and runs CLI/resource smoke checks. See `.github/workflows/release-gate.yml` for the authoritative CI sequence.
