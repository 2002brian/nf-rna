# nf-rna Quick Start

This guide is for production users on Linux or WSL. Install Python 3.11+, Git, Docker with a running daemon, Bash, Java 17+, and Nextflow before starting. Nextflow runs on the host; Docker runs analysis tasks. See the [official Nextflow installation guide](https://docs.seqera.io/nextflow/install) for Java and Nextflow.

## Install a released nf-rna version

nf-rna is installed from Git tags; it is not published to PyPI. Replace `X.Y.Z` with the released version you intend to use:

```bash
RELEASE_VERSION=X.Y.Z
python3.11 -m venv ~/.venvs/nf-rna-${RELEASE_VERSION}
source ~/.venvs/nf-rna-${RELEASE_VERSION}/bin/activate
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
docker pull ghcr.io/2002brian/nf-rna:${RELEASE_VERSION}
rnaseq --version
rnaseq doctor
```

The installed package contains required workflow assets. `rnaseq doctor` is read-only: it checks the local environment and image without pulling or building anything.

## Create and run a project

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
rnaseq status .
```

The wizard writes `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`; it does not import data unless explicit import flags are supplied. FASTQ projects require an execution-ready reference.

## References and inputs

Register an existing checksum-bound managed reference once:

```bash
rnaseq reference register /absolute/reference-root
```

Registration validates and records a local pointer to the existing manifest; it does not download, copy, or build assets. If you need to create an index rather than adopt one, use the optional builder workflow in [Runtime](runtime.md#managed-reference-builder-host-native).

For `input/counts.csv`, the first column is `gene_id` and remaining columns are sample IDs. For FASTQ, accurately declare `input.layout`, `input.preprocessing`, `upstream.strandedness`, the selected quantification backend, and a matching reference. HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` strandedness.

Only `--profile local` is implemented. HPC, SLURM, and Apptainer are not currently supported.

## Reproducible runtime

The release contract is:

```text
CLI X.Y.Z  ↔  Git tag vX.Y.Z  ↔  ghcr.io/2002brian/nf-rna:X.Y.Z
```

New projects use the image that matches the installed CLI version. A prerelease CLI uses its explicit matching prerelease image tag; publish and qualify that image before use. For formal analyses, use a versioned tag or an immutable digest rather than `:latest`.

`rnaseq doctor PROJECT` reports requested and Docker-observed image identity. Run provenance records the requested image, observed image ID/repository digest, OCI revision label when available, CLI version, and executed workflow hashes.

## v1.0.0 compatibility

The immutable `v1.0.0` release predates the version-matched project default. Its `rnaseq new` wizard writes `execution_image: nf-rna:latest`. Users running that historical release with the official GHCR image should set this field before planning or running:

```yaml
runtime:
  execution_image: ghcr.io/2002brian/nf-rna:1.0.0
```

For example, update the generated `project.yaml` in an editor after creating the project. This is a v1.0.0 compatibility note only; it is not part of the normal Quick Start.

The qualified immutable v1.0.0 image identity is:

```text
ghcr.io/2002brian/nf-rna@sha256:ee60405181783ff075a1f4a9f452990c651e91152f5a9837c2a5df44a838decb
```
