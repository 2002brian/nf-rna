English | [繁體中文](README_zh-TW.md)

# nf-rna

`nf-rna` is a reproducible bulk RNA-seq workflow for FASTQ and raw-count projects. It validates declared inputs and design, performs the selected QC and analysis stages, and produces figures, tables, a technical report, and frozen provenance.

`rnaseq` is the control plane, Nextflow is the execution plane, and Conda provides every process environment. The supported runtime for v1.4.1 is **Linux x86-64 or Windows WSL2, with Nextflow and Conda**; Docker is not required. Production users enter through `rnaseq`; they do not invoke internal Nextflow modules.

## Features

- FASTQ projects using nf-core/rnaseq with Salmon/tximport or first-party HISAT2 + featureCounts.
- Raw-count projects with validated metadata, explicit contrasts, QC, differential expression, and optional GO/KEGG ORA and preranked GSEA.
- Immutable runs with frozen configuration, the exact nf-rna source commit, the locked Conda runtime identity, workflow hashes, and curated delivery artifacts.

## Architecture

nf-rna provides one validated, count-based downstream workflow for two input routes. FASTQ projects first perform upstream quantification; count-matrix projects enter the downstream path after their supplied counts are validated. Both routes converge on the same expression analysis.

```mermaid
flowchart LR
    A[FASTQ] --> B["Quantification<br/>Salmon / tximport<br/>or HISAT2 / featureCounts"]
    C[Count matrix] --> D[Validated counts]
    B --> D
    D --> E["L1<br/>Expression QC, PCA,<br/>sample correlation"]
    E --> F{L2 selected?}
    F -- Yes --> G["L2: DESeq2<br/>optional GO/KEGG ORA<br/>and preranked GSEA"]
    F -- No --> H["Outputs<br/>figures and tables<br/>technical report<br/>curated delivery<br/>frozen provenance"]
    G --> H
```

L1 provides expression quality control and exploration. L2 is explicitly selected: it performs DESeq2 for declared contrasts and may be followed independently by GO ORA, KEGG ORA, and preranked GO/KEGG GSEA; it is never inferred from metadata or contrasts. Each route produces figures, tables, a technical report, curated delivery artifacts, and frozen provenance.

## Inputs and outputs

### Inputs

- **FASTQ:** place reads in `input/fastq/`. FASTQ projects quantify with the selected Salmon/tximport or HISAT2 + featureCounts route and require an execution-ready reference compatible with that route.
- **Count matrix:** place an imported gene-count matrix in `input/counts.csv` when quantification was completed elsewhere and nf-rna should perform the validated downstream analysis.
- **Design files:** `metadata.csv` provides sample IDs and design variables; `contrasts.csv` declares the directional contrasts used by L2. `rnaseq new` writes both as templates (a header-only `metadata.csv` with standard column names, and a header-only `contrasts.csv`); `rnaseq design` fills in `contrasts.csv` from your completed metadata.
- **References:** FASTQ projects declare the selected reference route explicitly. Reference identity and checksums are validated and frozen with the run rather than inferred from local files.

### Outputs

Every authorized execution creates an immutable run under `runs/<case-id>/<run-id>/`. The run keeps frozen inputs and configuration, the nf-rna source commit, Conda runtime and workflow identity, and provenance alongside its analysis artifacts.

- **Upstream count data:** FASTQ projects retain upstream QC and the Salmon/tximport or featureCounts count handoff; count-matrix projects retain their validated imported count source.
- **L1 QC:** filtering and normalization records, transformed expression for visualization, PCA, sample correlation, and QC tables and figures. Two per-sample boxplots compare sample distributions: `boxplot_log2_raw_counts` (log2(raw count + 1), before normalization, showing depth and composition differences) and `boxplot_vst` (the blind VST values used by PCA and correlation, showing whether normalization aligned the samples). Their Tukey box statistics per sample are in `sample_distribution_qc.tsv`. Plotted values are never batch-corrected.
- **L2 results:** when selected, DESeq2 results for each declared contrast, with associated tables and figures.
- **Enrichment:** independently enabled GO ORA, KEGG ORA, and preranked GO (BP, MF, CC)/KEGG GSEA results. From v1.5 a term is significant when p.adjust (Benjamini-Hochberg) ≤ 0.05, and the ORA universe is the genes with a non-NA DESeq2 padj, mapped to Entrez and intersected with the gene-set database. KEGG gene sets come from a local snapshot of KEGG REST (`~/.cache/nf-rna/genesets/kegg/`); without a pin, the latest local snapshot is used and `rnaseq validate` reports it as not pinned. **For results you will publish, pin the snapshot** in `project.yaml` with `annotation.enrichment.kegg.snapshot: <snapshot ID>` (the ID is shown by `rnaseq validate` and recorded in the KEGG provenance).
- **Delivery:** a technical HTML report and a curated delivery package of figures, tables, methods/version records, provenance, and a delivery manifest.

## Requirements

Use a Linux x86-64 workstation or Windows with WSL2 (native Windows and macOS are not supported by v1.4.1) with:

- Python 3.11+ and Git;
- Conda (for example Miniforge) on `PATH`, with the `conda-forge` and `bioconda` channels configured in that order, as nf-core requires (`conda config --add channels bioconda && conda config --add channels conda-forge && conda config --set channel_priority strict`);
- Bash, Java 17+, and Nextflow on `PATH`; and
- sufficient local storage for references, Conda environments, and Nextflow work data.

Docker is not required. Java, Nextflow, and Conda run on the host. nf-core/rnaseq 3.26.0 runs with `-profile conda`; HISAT2/featureCounts use an exactly pinned Conda environment; R, DESeq2, tximport, and the enrichment packages run in a locked linux-64 Conda environment that nf-rna creates on first use. Normal users do not install these tools by hand. See the [Nextflow installation guide](https://docs.seqera.io/nextflow/install).

## Installation

> **Latest stable release: `v1.4.1`.** This README describes the patch release `v1.4.1` (see the [1.4.1 changes](CHANGELOG.md#141--2026-10-29)): corrected wizard instructions, clearer validation and upstream-reuse errors, and the beginner workflow documentation below. The analysis itself is unchanged from `v1.4.0`. For an earlier release, use the README at that release's tag. To use the development version, follow [Development installation](docs/development.md).

nf-rna has no PyPI distribution. Install the released CLI from its Git tag.
Installing from Git lets pip record the exact source commit, which every run
must record; nf-rna refuses to run from an install that has no commit, such as
a plain directory, sdist, or wheel install. Replace `X.Y.Z` with the released
version you intend to use:

```bash
RELEASE_VERSION=X.Y.Z
conda create --yes --name nf-rna python=3.11
conda activate nf-rna
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
rnaseq --version
rnaseq doctor
```

The first run creates the locked downstream Conda environment under the nf-rna
runtime cache (one environment per lock and installed nf-rna build); later runs
reuse it. nf-core/rnaseq and HISAT2/featureCounts process environments are
created by Nextflow in a shared upstream Conda cache.

The release contract is deliberate:

```text
CLI X.Y.Z  ↔  Git tag vX.Y.Z  ↔  one source commit, recorded in every run
```

`rnaseq doctor` is read-only: it checks Java, Nextflow, Conda and its channel
configuration, the nf-core/rnaseq pin, the downstream Conda lock, and the nf-rna
source revision, then prints `Overall: READY` or `Overall: NOT READY` (exit code 1).
It does not create environments or download pipelines.

### Historical Docker releases

Releases `v1.1.2` through `v1.2.1` ran their processes in Docker and published
official images as `ghcr.io/2002brian/nf-rna:X.Y.Z`, qualified for
`linux/amd64`. `v1.2.1` is the last Docker-based release; v1.3.0 does not
require or qualify a Docker image. `v1.0.0` predates GHCR publication, and
`v1.1.0`/`v1.1.1` remain valid source releases without official images because
their publication workflows stopped before image build or push.

## Quick Start

There are two supported workflows. Both start from the template written by `rnaseq new`:

- **Staged L1 → L2** (default): run L1 for QC, review it, then let `rnaseq design` configure differential expression and run L2.
- **Direct L2**: choose `L2` in the wizard, run `rnaseq design` before the first run, and run L2 once. An L2 run always performs the L1 QC as well, so no separate L1 run is needed.

`rnaseq design` is what creates the contrasts: it reads the completed `metadata.csv`, proposes the contrasts, writes `contrasts.csv`, sets the scope to L2, and refreshes the planning artifacts (you do not need to re-run `rnaseq plan` afterwards). It does not need a previous run. An L2 project with an empty `contrasts.csv` fails `rnaseq validate` and `rnaseq plan` until `rnaseq design` has been run.

### FASTQ file names and sample IDs

Put gzipped reads (`.fastq.gz` or `.fq.gz`) in `input/fastq/`. nf-rna reads the sample ID from each file name: the name must end in `_R1` or `_R2` (optionally followed by `_001`), optionally preceded by a three-digit lane `_L001`, and **the sample ID is everything before the lane/read marker**:

| File name | Sample ID | Lane |
|---|---|---|
| `Ctrl1_S1_L001_R1_001.fastq.gz` | `Ctrl1_S1` | `L001` |
| `Ctrl1_S1_L002_R1_001.fastq.gz` | `Ctrl1_S1` | `L002` (lanes of one sample are combined) |
| `Ctrl1_R1.fastq.gz` | `Ctrl1` | – |

The Illumina sample number (`_S1`) is part of the sample ID. Each `sample_id` in `metadata.csv` must match exactly, so for `Ctrl1_S1_L001_R1_001.fastq.gz` write `Ctrl1_S1` (or rename the files to `Ctrl1_L001_R1_001.fastq.gz` etc. before validating). Paired-end data needs an `_R2` file for every `_R1` file of the same sample and lane. Names such as `Ctrl1_1.fastq.gz`, `Ctrl1.R1.fastq.gz`, and uncompressed `.fastq` files are rejected. `rnaseq validate` lists every sample found only in the FASTQ files or only in `metadata.csv`.

### Workflow A: L1, then L2

```bash
mkdir -p ~/projects/rnaseq-projects
cd ~/projects/rnaseq-projects
rnaseq new                 # choose L1, the design, and whether there is a known technical batch
cd <new-project>
# Add FASTQs to input/fastq/ (or a count matrix to input/counts.csv) and complete metadata.csv.
rnaseq validate .
rnaseq plan .
rnaseq doctor .
rnaseq run . --case-id L1-001 --profile local --yes     # L1: QC and exploratory analysis (foreground; see below)
rnaseq status . --all      # the RUN column shows the run ID of the SUCCESS run of case L1-001
# Review the L1 QC in runs/L1-001/<run-id>/delivery/: sample boxplots, PCA, library sizes.
rnaseq design .            # proposes contrasts and asks for confirmation; also refreshes the plan
rnaseq run . --case-id L2-001 --profile local --reuse-upstream L1-001/<run-id> --yes   # FASTQ projects
rnaseq run . --case-id L2-001 --profile local --yes     # raw-count projects (no upstream to reuse)
```

### Workflow B: direct L2

```bash
rnaseq new                 # choose L2
cd <new-project>
# Add FASTQs to input/fastq/ (or a count matrix to input/counts.csv) and complete metadata.csv.
rnaseq design .            # required first: writes contrasts.csv and refreshes the plan
rnaseq validate .          # expect "Validation result: PASS" and "Execution readiness: READY"
rnaseq doctor .
rnaseq run . --case-id L2-001 --profile local --yes     # one run: upstream, L1 QC, DESeq2, enrichment
```

### Case ID and run ID

- The **case ID** is the name you choose with `--case-id` (letters, digits, `.`, `_`, `-`; for example `L1-001`). Running the same case ID again creates another run under that case; it never overwrites one.
- The **run ID** is created by nf-rna for every execution: its start time in Asia/Taipei time, `YYYYMMDD-HHMMSS+0800` (for example `20261029-101010+0800`, with a `-01` suffix if two runs start in the same second). `rnaseq run` prints both at the end (`Case run: SUCCESS (L1-001/20261029-101010+0800)`), the RUN column of `rnaseq status . --all` lists them, and every run lives in `runs/<case-id>/<run-id>/`.
- Commands that select one run take both: `rnaseq status . --case L1-001 --run <run-id>`, `rnaseq clean . --case-id L1-001 --run <run-id>`, and the single value `CASE-ID/RUN-ID` for `--reuse-upstream` and `rnaseq retry --retry-of`.

### Running, monitoring, and interruptions

`rnaseq run` runs in the foreground: the terminal stays busy until the run finishes, which can take hours for FASTQ projects. Follow it from a **second terminal**:

```bash
rnaseq status .            # execution dashboard of the latest run
rnaseq status . --watch    # redraw until SUCCESS, FAILED or INTERRUPTED (Ctrl-C stops watching only)
rnaseq status . --all      # one line per recorded run
```

Keep the machine awake and the run's terminal open until the run ends:

- Ctrl-C, closing the terminal or SSH session, or logging out stops Nextflow and records the run as INTERRUPTED.
- Shutdown, restart, Windows updates that restart the machine, and `wsl --shutdown` kill the run without recording an end state; `rnaseq status` then shows it as INTERRUPTED (stale). Sleep or hibernation suspends it, and long pauses can make tasks fail. On Windows, set the power plan so the PC does not sleep during a run.
- A recorded FAILED or INTERRUPTED run can be retried with `rnaseq retry . --retry-of CASE-ID/RUN-ID`; for a stale run, start a new run.

To keep a run going when the terminal is closed, start it inside `tmux` (install it separately; detach with Ctrl-b d, return with `tmux attach`) or with `nohup` (requires `--yes`, because no confirmation prompt can be answered):

```bash
nohup rnaseq run . --case-id L1-001 --profile local --yes > L1-001.out 2>&1 &
```

Neither protects against shutdown, sleep, or `wsl --shutdown`.

### Results and cleanup

Final results are delivered to `runs/<case-id>/<run-id>/delivery/`. Each run's disposable Nextflow work directory (often hundreds of GB for FASTQ runs) is kept until you remove it. Check what would be removed first, then clean:

```bash
rnaseq clean . --case-id L1-001 --run <run-id> --dry-run   # shows the path, state and size; deletes nothing
rnaseq clean . --case-id L1-001 --run <run-id>             # asks before deleting
```

Cleanup removes only the work directory; `runs/<case-id>/<run-id>/` (including `delivery/`) is never removed, and a cleaned L1 run can still be reused by L2. See [Runtime](docs/runtime.md#nextflow-work-directories-and-rnaseq-clean).

### What the wizard asks, and what the templates contain

Besides the species, input, and reference choices, `rnaseq new` asks for the analysis scope. The default is `L1`, the start of the normal journey; FASTQ projects can also choose `qc`, and `rnaseq design` later moves an L1 project to L2. For the design, it asks only what cannot be read from your metadata later: the design (two-group, paired two-group, or multi-group), the pairing column for paired designs, and one yes/no question: *is there a known technical batch?* It never asks for condition levels, batch levels, a reference level, or contrasts. The generated `metadata.csv` header depends on those answers:

| Design | Known batch | `metadata.csv` header | Provisional model |
|---|---|---|---|
| two-group / multi-group | no | `sample_id,condition` | `~ condition` |
| two-group / multi-group | yes | `sample_id,condition,batch` | `~ batch + condition` |
| paired | no | `sample_id,<pair_id>,condition` | `~ <pair_id> + condition` |
| paired | yes | `sample_id,<pair_id>,condition,batch` | `~ <pair_id> + batch + condition` |

`contrasts.csv` starts as the header only: `contrast_id,factor,numerator,denominator`. Both files are ordinary CSV files you can open and edit. The same choices are available non-interactively as `--preset`, `--design-type`, `--pair-id`, and `--batch/--no-batch`.

### L1: QC and exploratory analysis

L1 needs no contrasts. It produces the raw-count and VST sample boxplots, PCA, sample correlations, and library-size QC. L1 is QC/exploratory only:

- A `batch` term in the provisional model does **not** mean the plotted values are batch-corrected. No `removeBatchEffect` or corrected expression matrix is created; boxplots and PCA are always uncorrected, and a `batch` column is shown only as an annotation.
- nf-rna never uses these plots to decide whether batch belongs in the model; you declare known batch at project creation.
- Formal differential inference happens in L2.

### `rnaseq design`: configuring L2

`rnaseq design PROJECT` reads the completed `metadata.csv` and prepares L2, after reviewing L1 (Workflow A) or before the first run (Workflow B):

- It validates the design before anything is written: the model matrix must be full rank (no perfectly confounded batch and condition, and for a paired design batch must vary within pairs) and must leave residual degrees of freedom after categorical variables are expanded. Partial confounding and single-sample batches are reported as warnings, and the batch-by-condition sample counts are shown.
- **Two-group:** it reads the two condition levels in the order they first appear in `metadata.csv` and proposes *later vs earlier*, for example `Treatment_vs_Control` (log2FC > 0 = higher in Treatment). Answer `y` to accept, `swap` to reverse the direction, or `n` to change nothing (`--yes`, `--swap` for non-interactive use). `swap` changes only the contrast in `contrasts.csv` (the authority on direction); `design.condition_order` stays the first-appearance (or `--condition-order`) order.
- **Multi-group:** it generates every pairwise contrast, later level vs earlier level, in first-appearance order. For `Control, DrugA, DrugB` it writes `DrugA_vs_Control`, `DrugB_vs_Control`, `DrugB_vs_DrugA`. k levels give k(k-1)/2 contrasts. There is no per-comparison selection; edit `contrasts.csv` to remove comparisons. Per-contrast adjusted p-values are not corrected across contrasts.
- Use `--condition-order Control,DrugA,DrugB` to set the order explicitly (reference level first).
- Only after you confirm does it write `contrasts.csv`, freeze the order as `design.condition_order` in `project.yaml` (so re-sorting metadata rows later cannot flip a contrast), set the scope to `L2`, and refresh the planning artifacts. The project's `project.yaml` formula stays as derived: `~ batch + condition`, or `~ <pair_id> + batch + condition` for paired designs.
- An existing non-empty `contrasts.csv` is never overwritten unless you pass `--force`. QC-only projects and projects with invalid metadata/design are refused. A project already set to L2 with empty contrasts is accepted.

### Reusing L1 upstream results in L2 (FASTQ projects)

For a FASTQ project, the L2 run can reuse the upstream quantification of a successful L1 run instead of repeating it. `--reuse-upstream` takes exactly `CASE-ID/RUN-ID` of that run, for example:

```bash
rnaseq run . --case-id L2-001 --profile local --reuse-upstream L1-001/20261029-101010+0800 --yes
```

Take the run ID from the RUN column of the SUCCESS row of `rnaseq status . --all` (a case can hold several runs, including FAILED ones); `rnaseq design` also prints the newest successful run in its `Next:` line. Both Salmon/tximport and HISAT2 + featureCounts projects are supported. Compatibility checks must pass:

- the source must be a SUCCESS run of the same project with a validated upstream handoff;
- its frozen nf-core parameters, Nextflow config, and samplesheet must match the new run byte for byte;
- its input, upstream, and reference identity must match the new run.

The frozen Nextflow config records the local CPU/memory ceiling. Changing `execution.max_cpus` or `execution.max_memory_gb` between the two runs, or running with `auto` on a machine whose capacity differs from the one that ran L1, therefore makes reuse incompatible: keep the same limits and machine, or run L2 without `--reuse-upstream`. If a check fails, the new run fails with an incompatibility error and no upstream results are copied. The analysis scope, design, metadata, and contrasts may differ, so the changes made by `rnaseq design` do not block reuse. The upstream results are copied from the durable L1 run record, so removing the L1 run's Nextflow work with `rnaseq clean` does not prevent reuse. Raw-count projects have no upstream stage; do not use this option for them.

New projects use hardware-aware automatic CPU/memory limits (`execution.max_cpus: auto`, `execution.max_memory_gb: auto`); explicit integer limits in `project.yaml` are preserved. See [Local resource policy](docs/runtime.md#local-resource-policy).

The wizard creates a project scaffold; it does not infer scientific inputs. FASTQ projects also need an execution-ready reference. Register an existing checksum-bound managed reference with `rnaseq reference register /absolute/reference-root`, or configure a supported reference route. See the [detailed Quick Start](docs/quickstart.md).

## Covariate-aware DESeq2 designs

New schema 1.3 projects declare the type of every additive formula variable. Categorical variables are fitted as R factors; continuous variables are finite numeric adjustment covariates. Existing schema 1.0–1.2 projects remain readable with their historical categorical/clearly numeric behavior.

```yaml
design:
  type: multi_group
  formula: "~ batch + age + condition"
  variables:
    batch: categorical
    age: continuous
    condition: categorical
```

Supported additive examples are `~ condition`, `~ age + condition`, `~ batch + condition`, and `~ batch + age + condition`. Contrasts remain categorical and directional—for example `Treatment_vs_Control,condition,Treatment,Control`; a continuous variable adjusts the estimate and cannot be used as a numerator/denominator contrast.

### Advanced: imported metadata, custom covariates, contrasts, and formulas

The simple wizard covers `condition` and an optional `batch`. For anything else, import your own files with explicit flags (flags-only; the wizard never imports): `--metadata`, `--contrasts`, `--condition-column`, `--covariate batch` for categorical adjustments, `--continuous-covariate age` for numeric adjustments, and `--pair-id`. `--batch` with `--metadata` is shorthand for `--covariate batch`. These options are never ignored: without `--metadata` they are rejected with an explanation, and a metadata file with no `condition` column requires `--condition-column`. Imported contrasts are preserved, and `rnaseq design` will not replace them unless you pass `--force`; a formula that is not the standard derivation is kept as written.

Known batch belongs in the DESeq2 design (`~ batch + condition`), not in a batch-corrected raw count matrix. nf-rna preserves original counts for inference, records the fitted typed design in provenance and reports, and adds batch metadata (plus a batch-colored PCA when `batch` is declared categorical) to exploratory VST QC. Continuous-effect hypothesis tests are outside this milestone.

## Reproducibility

Every run records the nf-rna source commit, the downstream Conda lock and its SHA-256, the installed nf-rna wheel SHA-256, the bundled R-script inventory, the nf-core/rnaseq version and revision, and the reference manifest and file checksums. nf-rna refuses to start a run when the source commit is unknown or the source checkout has uncommitted changes.

For formal analyses, install an explicit release tag. `rnaseq doctor PROJECT` reports the source revision, the Conda lock, and the provisioned runtime before you run.

## Documentation

- [Detailed Quick Start and v1.0.0 compatibility](docs/quickstart.md)
- [Runtime and reproducibility contract](docs/runtime.md)
- [Scientific contract](docs/scientific_contract.md)
- [Architecture](docs/architecture.md)

## Development

Cloning the source, editable installation, and tests are contributor activities. See [Development installation](docs/development.md).
