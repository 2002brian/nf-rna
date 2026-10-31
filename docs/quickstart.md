# nf-rna Quick Start

This guide is for production users on Linux x86-64 or Windows WSL2. Install Python 3.11+, Git, Conda (for example Miniforge) with the `conda-forge` and `bioconda` channels configured in that order (`conda config --add channels bioconda && conda config --add channels conda-forge && conda config --set channel_priority strict`), Bash, Java 17+, and Nextflow before starting. Nextflow runs on the host and runs every analysis task in a Conda environment; Docker is not required. See the [official Nextflow installation guide](https://docs.seqera.io/nextflow/install) for Java and Nextflow.

## Install a released nf-rna version

The latest stable release is `v1.6.1`, and this guide describes it (see the [1.6.1 changes](../CHANGELOG.md#161--2026-10-31)). **v1.6.1 changes results:** GSEA leading edges and Monte Carlo results, DESeq2 padj, DEG sets and ORA can differ from `v1.6.0`, so recompute projects to obtain corrected results; it also adds a strandedness check (HISAT2 route) and records upstream thread counts and R warnings without changing any analysis value. `v1.6.0` added the optional browser GUI (see [Browser GUI](#browser-gui-optional)) and includes the v1.5.1 test fixes. **If you upgrade from v1.4.x, read the [1.5.0 changes](../CHANGELOG.md#150--2026-10-30) first:** in v1.5.0 GO/KEGG enrichment calls a term significant when the Benjamini-Hochberg adjusted p-value is ≤ 0.05, and KEGG gene sets come from a local snapshot: for results you will publish, pin the snapshot ID that `rnaseq validate` reports with `annotation.enrichment.kegg.snapshot` in `project.yaml`. For an earlier release, follow the documentation at that release's tag; to use the development version, see [Development installation](development.md).

nf-rna is installed from Git tags; it is not published to PyPI. Installing from
Git lets pip record the exact source commit, which every run records; nf-rna
refuses to run from an install without a commit. Replace `X.Y.Z` with the
released version you intend to use:

```bash
RELEASE_VERSION=X.Y.Z
conda create --yes --name nf-rna python=3.11
conda activate nf-rna
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
rnaseq --version
rnaseq doctor
```

The installed package contains required workflow assets. `rnaseq doctor` is read-only: it checks Java, Nextflow, Conda and its channel order, the nf-core/rnaseq pin, the downstream Conda lock, and the nf-rna source revision without creating environments or downloading anything, and ends with `Overall: READY` or `Overall: NOT READY`. The first run creates the locked downstream Conda environment; Nextflow creates the nf-core/rnaseq and HISAT2/featureCounts process environments in a shared Conda cache.

Releases `v1.1.2` through `v1.2.1` used a Docker execution image published as
`ghcr.io/2002brian/nf-rna:X.Y.Z`; `v1.2.1` is the last Docker-based release.
`v1.0.0` predates GHCR publication, and `v1.1.0`/`v1.1.1` have no official image.

## Create and run a project

Both supported workflows are template-first:

```text
Staged L1 -> L2: rnaseq new (L1) -> add inputs, fill metadata.csv -> validate -> plan -> run L1 -> review QC -> rnaseq design -> run L2
Direct L2:       rnaseq new (L2) -> add inputs, fill metadata.csv -> rnaseq design -> validate -> run L2 (includes the L1 QC)
```

Staged L1 -> L2:

```bash
mkdir -p ~/projects/rnaseq-projects
cd ~/projects/rnaseq-projects
rnaseq new                 # choose L1
cd <new-project>
# Add FASTQs to input/fastq/, or a count matrix to input/counts.csv, and complete metadata.csv.
rnaseq validate .
rnaseq plan .
rnaseq doctor .
rnaseq run . --case-id L1-001 --profile local --yes   # L1: QC and exploratory analysis
rnaseq status . --all     # RUN column: the run ID of the SUCCESS run of L1-001
# Review the L1 QC (sample boxplots, PCA) in runs/L1-001/<run-id>/delivery/, then confirm the model and contrasts:
rnaseq design .
rnaseq run . --case-id L2-001 --profile local --reuse-upstream L1-001/<run-id> --yes   # FASTQ projects
rnaseq run . --case-id L2-001 --profile local --yes                                    # raw-count projects
rnaseq clean . --case-id L1-001 --run <run-id> --dry-run   # after SUCCESS: show what would be removed
rnaseq clean . --case-id L1-001 --run <run-id>             # remove the disposable Nextflow work
```

Direct L2:

```bash
rnaseq new                 # choose L2
cd <new-project>
# Add the inputs and complete metadata.csv, then configure the contrasts before the first run:
rnaseq design .            # writes contrasts.csv and refreshes the planning artifacts
rnaseq validate .
rnaseq doctor .
rnaseq run . --case-id L2-001 --profile local --yes   # upstream, L1 QC, DESeq2 and enrichment in one run
```

An L2 project fails `rnaseq validate` and `rnaseq plan` until `rnaseq design` has written its contrasts; a separate L1 run is never needed for L2. The wizard writes `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`; it does not import data unless explicit import flags are supplied. FASTQ projects require an execution-ready reference.

### FASTQ names and sample IDs

FASTQs must be gzipped (`.fastq.gz` or `.fq.gz`) and named `<sample>[_L<nnn>]_R1|_R2[_<nnn>]`. The sample ID is everything before the optional lane and the read marker, so `Ctrl1_S1_L001_R1_001.fastq.gz` belongs to sample `Ctrl1_S1` (the Illumina `_S1` is kept) and `Ctrl1_R1.fastq.gz` to `Ctrl1`. Lanes of one sample are combined, paired-end data needs a matching `_R2` for each `_R1`, and `_1`/`_2` or `.R1` names are rejected. Each `sample_id` in `metadata.csv` must match exactly; `rnaseq validate` lists the samples found on only one side.

### Case ID, run ID, and monitoring

`--case-id` is a name you choose (for example `L1-001`). nf-rna creates a run ID for every execution, its Asia/Taipei start time `YYYYMMDD-HHMMSS+0800`, prints `Case run: SUCCESS (CASE-ID/RUN-ID)` at the end, and stores the run in `runs/<case-id>/<run-id>/`; final results are in its `delivery/` directory. `rnaseq status . --all` lists both IDs. `--reuse-upstream` and `rnaseq retry --retry-of` take `CASE-ID/RUN-ID`; `rnaseq clean` takes `--case-id` and `--run`.

`rnaseq run` runs in the foreground. Watch it from a second terminal with `rnaseq status . --watch` (Ctrl-C there stops only the watching). Closing the run's terminal or SSH session, logging out, or Ctrl-C stops the run and records it as INTERRUPTED; shutdown, restart, sleep, and `wsl --shutdown` stop or suspend it as well, so keep the machine awake until it finishes. To survive a closed terminal, start the run inside `tmux` or with `nohup rnaseq run . --case-id L1-001 --profile local --yes > L1-001.out 2>&1 &`; neither survives shutdown or sleep.

### Metadata templates

The wizard asks for the analysis scope, defaulting to `L1` (FASTQ projects can also choose `qc`; non-interactively, `--preset`); `rnaseq design` later moves an L1 project to L2. It then asks for the design (two-group, paired two-group, multi-group), the pairing column for paired designs, and one question: is there a known technical batch? It does not ask for condition levels, batch levels, reference levels, or contrasts; those are discovered from the completed metadata. The standard column names are `sample_id`, `condition`, and `batch`:

| Design | Known batch | `metadata.csv` header |
|---|---|---|
| two-group / multi-group | no | `sample_id,condition` |
| two-group / multi-group | yes | `sample_id,condition,batch` |
| paired | no | `sample_id,<pair_id>,condition` |
| paired | yes | `sample_id,<pair_id>,condition,batch` |

`contrasts.csv` starts header-only (`contrast_id,factor,numerator,denominator`). The provisional model written to `project.yaml` is `~ condition`, `~ batch + condition`, `~ <pair_id> + condition`, or `~ <pair_id> + batch + condition`.

### L1 is QC, not batch correction

L1 runs without contrasts. It uses the provisional model only as the DESeq2 object its normalization and QC need. A `batch` term in that model does not mean the plotted values are batch-corrected: no `removeBatchEffect` and no corrected expression matrix is created, and the raw-count boxplot, VST boxplot and PCA are uncorrected. Formal differential inference happens in L2.

The two L1 boxplots compare sample distributions: `boxplot_log2_raw_counts` shows log2(raw count + 1) before normalization (depth and composition differences), and `boxplot_vst` shows the blind VST values used by PCA and correlation (whether normalization aligned the samples). `sample_distribution_qc.tsv` records each sample's Tukey box statistics for both.

### Configuring L2 with `rnaseq design`

`rnaseq design .` works before or after an L1 run. It validates the model (full rank, no perfect confounding of batch and condition, and residual degrees of freedom after categorical variables are expanded), prints the batch-by-condition sample counts, and proposes contrasts from the condition levels in the order they first appear in `metadata.csv`:

- two-group: `Treatment_vs_Control` (later vs earlier; log2FC > 0 means higher in Treatment); answer `swap` to reverse it (only the contrast changes; `design.condition_order` keeps the first-appearance order);
- multi-group: all pairwise contrasts, later vs earlier (`DrugA_vs_Control`, `DrugB_vs_Control`, `DrugB_vs_DrugA`).

Nothing is written until you confirm. Confirmation writes `contrasts.csv`, freezes `design.condition_order` in `project.yaml`, sets the scope to L2, and refreshes the planning artifacts (so `rnaseq run` does not report a stale plan; you do not need to re-run `rnaseq plan`). An existing non-empty `contrasts.csv` is kept unless you pass `--force`; `--condition-order A,B,C` sets the order explicitly and `--yes` accepts non-interactively. For custom metadata, covariates, contrasts, or formulas use the explicit import flags (`--metadata`, `--contrasts`, `--condition-column`, `--covariate`, `--continuous-covariate`, `--pair-id`), which are never silently ignored.

### Reusing L1 upstream results (FASTQ projects)

For a FASTQ project (Salmon/tximport or HISAT2 + featureCounts), the L2 run can reuse the upstream quantification of a successful L1 run instead of repeating it. The value is exactly `CASE-ID/RUN-ID`:

```bash
rnaseq run . --case-id L2-001 --profile local --reuse-upstream L1-001/20261029-101010+0800 --yes
```

Take the run ID from the SUCCESS row of `rnaseq status . --all`. Compatibility checks must pass:

- the source must be a SUCCESS run of the same project with a validated upstream handoff;
- its frozen nf-core parameters, Nextflow config, and samplesheet must match the new run byte for byte;
- its input, upstream, and reference identity must match the new run.

The frozen Nextflow config contains the local CPU/memory ceiling, so changing `execution.max_cpus` or `execution.max_memory_gb`, or using `auto` on a machine with a different capacity, makes the runs incompatible; keep the same limits and machine, or run without `--reuse-upstream`. Otherwise the new run fails with an incompatibility error and no upstream results are copied. Scope, design, metadata, and contrasts may differ, so the changes made by `rnaseq design` do not block reuse. Results are copied from the durable L1 run record, so `rnaseq clean` on the L1 run does not prevent reuse. Raw-count projects have no upstream stage; do not use this option for them.

## Browser GUI (optional)

The same workflow is available in a local browser GUI, from project creation to results and downloads:

```bash
python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.1"
tmux new -s nf-rna-gui
rnaseq gui                                   # serves http://127.0.0.1:8501 on the server
ssh -N -L 8501:127.0.0.1:8501 user@server    # on your laptop, then open http://127.0.0.1:8501
```

The GUI only offers folders listed in `~/.config/nf-rna/gui.yaml`; see [Browser GUI](gui.md) for the pages, settings and the optional `systemd --user` service.

## References and inputs

Register an existing checksum-bound managed reference once:

```bash
rnaseq reference register /absolute/reference-root
```

Registration validates and records a local pointer to the existing manifest; it does not download, copy, or build assets. If you need to create an index rather than adopt one, use the optional builder workflow in [Runtime](runtime.md#managed-reference-builder-host-native).

How `rnaseq new` picks the reference for a FASTQ project:

- **Interactive wizard:** it offers the compatible registered references for the species and backend; declining one leads to the iGenomes, local, or custom choice.
- **Non-interactive (flags):** the registry is never consulted. Name a registered or other managed reference explicitly:

  ```bash
  rnaseq new ... --method salmon --reference-source local --reference-root /absolute/reference-root
  ```

  `--reference-root` or `--reference-manifest` alone also select `local`, and the custom flags (`--reference-fasta`, `--reference-gtf`, ...) select `custom`. With no reference flag at all, Salmon uses iGenomes (`GRCh38` for human, `GRCm39` for mouse) and HISAT2 + featureCounts defaults to `local`, so it needs `--reference-root` (or the custom flags). Reference flags that do not belong to the selected source are refused rather than ignored.

The project review printed by `rnaseq new` states the selected reference and whether it is the default, and `planning/analysis_plan.md` records the reference source (for iGenomes, the genome).

For `input/counts.csv`, the first column is `gene_id` and remaining columns are sample IDs. For FASTQ, accurately declare `input.layout`, `input.preprocessing`, `upstream.strandedness`, the selected quantification backend, and a matching reference. HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` strandedness.

Only `--profile local` is implemented. HPC, SLURM, and Apptainer are not currently supported.

## Reproducible runtime

The release contract is:

```text
CLI X.Y.Z  ↔  Git tag vX.Y.Z  ↔  one source commit, recorded in every run
```

Install an explicit release tag for formal analyses. `rnaseq doctor PROJECT`
reports the source revision, the downstream Conda lock, and any provisioned
runtime. Run provenance records the nf-rna version and source commit, the
downstream Conda lock filename and SHA-256, the installed wheel SHA-256, the
bundled R-script inventory, the nf-core/rnaseq version and revision, the
reference manifest and file checksums, and executed workflow hashes.

## Historical v1.0.0 compatibility

The immutable `v1.0.0` release predates both GHCR publication and the
version-matched project default. It has no official GHCR image. Its `rnaseq
new` wizard writes `execution_image: nf-rna:latest`; users working from that
historical source release should build that local development image before
planning or running:

```yaml
runtime:
  execution_image: nf-rna:latest
```

For example, update the generated `project.yaml` in an editor after creating the project. This is a v1.0.0 compatibility note only; it is not part of the normal Quick Start.
