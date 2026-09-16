English | [繁體中文](README_zh-TW.md)

# nf-rna

`nf-rna` is a reproducible bulk RNA-seq workflow for FASTQ and raw-count projects. `rnaseq` is the control plane, Nextflow is the execution plane, and Docker is the process runtime. Production users enter through `rnaseq`; they do not build an execution image or invoke internal Nextflow modules.

The currently published release is `v1.0.0`. The installation UX changes in this branch are intended for the next patch release; no new release is implied by this documentation.

## Canonical installation for the next published patch release

After a patch release containing these changes is published, substitute its actual published version below. Do not use this block to imply that `v1.0.1` already exists.

```bash
RELEASE_VERSION=<published-version>
python3.11 -m venv ~/.venvs/nf-rna-${RELEASE_VERSION}
source ~/.venvs/nf-rna-${RELEASE_VERSION}/bin/activate
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
docker pull ghcr.io/2002brian/nf-rna:${RELEASE_VERSION}
rnaseq --version
rnaseq doctor
```

This is a non-editable Git installation; no PyPI distribution exists. The CLI derives a new project's default image from its own package version, so a released CLI `X.Y.Z` creates `ghcr.io/2002brian/nf-rna:X.Y.Z`. The same applies to a prerelease such as `1.0.1rc1`, which selects the explicit `:1.0.1rc1` tag if that image is published.

Host prerequisites are Python 3.11+ with `venv`, Git, Docker with a running daemon, Bash, Java 17+, and Nextflow on `PATH`. Java and Nextflow run on the host. R, DESeq2, HISAT2, featureCounts, SAMtools, and Salmon are container concerns. See the [Nextflow installation guide](https://docs.seqera.io/nextflow/install).

## Canonical Quick Start for that release

```bash
mkdir -p ~/projects/rnaseq-projects
cd ~/projects/rnaseq-projects
rnaseq new
cd <new-project>
# Add FASTQs to input/fastq/, or a count matrix to input/counts.csv.
# Complete metadata.csv and contrasts.csv.
rnaseq validate .
rnaseq plan .
rnaseq doctor .
rnaseq run . --case-id CASE-001 --profile local --yes
```

No YAML replacement is required in that path. FASTQ projects additionally need an execution-ready reference; register an existing checksum-bound managed reference with `rnaseq reference register /absolute/reference-root`, or configure a supported reference route. Details are in the [quick start](docs/quickstart.md).

## Historical v1.0.0 compatibility

`v1.0.0` is immutable and its wizard writes `execution_image: nf-rna:latest`. It remains installable with its matching qualified runtime:

```bash
python3.11 -m venv ~/.venvs/nf-rna-1.0.0
source ~/.venvs/nf-rna-1.0.0/bin/activate
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.0.0"
docker pull ghcr.io/2002brian/nf-rna:1.0.0
```

For a v1.0.0 project, replace the historic local image reference before planning or running:

```bash
sed -i 's|execution_image: nf-rna:latest|execution_image: ghcr.io/2002brian/nf-rna:1.0.0|' project.yaml
```

This is a historical compatibility step, not the future happy path.

## Runtime identity and updates

`rnaseq doctor PROJECT` reports the requested image and Docker-observed image ID, repository digest, and OCI revision label. Formal analyses should use an explicit release tag or the qualified immutable v1.0.0 digest:

```text
ghcr.io/2002brian/nf-rna@sha256:ee60405181783ff075a1f4a9f452990c651e91152f5a9837c2a5df44a838decb
```

`:latest` is a convenience tag only. Do not use it, `git pull`, or unpinned image updates inside an active production analysis.

## Documentation

- [Detailed Quick Start and reference setup](docs/quickstart.md)
- [Runtime and reproducibility contract](docs/runtime.md)
- [Scientific contract](docs/scientific_contract.md)
- [Architecture](docs/architecture.md)
- [Development installation](docs/development.md)

Cloning the source, editable installation, local image builds, and tests are contributor activities and are intentionally kept out of the production path.
