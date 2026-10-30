# nf-rna v1.6.0

> **Upgrading from v1.4.x? The enrichment significance semantics changed in v1.5.0.**
> GO/KEGG ORA and GSEA now call a term significant when the Benjamini–Hochberg
> adjusted p-value (p.adjust) is ≤ 0.05 (v1.4.x: ORA used pvalue ≤ 0.05 and
> qvalue ≤ 0.2; GSEA used p.adjust ≤ 0.05 and pvalue ≤ 0.05). The ORA universe is
> now the genes with a non-NA DESeq2 padj, mapped to Entrez and intersected with
> the gene-set database, and KEGG gene sets come from a local, SHA-256-verified
> snapshot — **pin `annotation.enrichment.kegg.snapshot` for results you will
> publish.** Results produced with v1.4.x are not recomputed and keep their
> earlier semantics; do not compare enrichment results across that boundary
> without re-running. Details: the
> [1.5.0 CHANGELOG entry](https://github.com/2002brian/nf-rna/blob/main/CHANGELOG.md#150--2026-10-30).

v1.6.0 adds an optional local browser GUI and includes the v1.5.1 test fixes
(v1.5.1 was not released separately). DESeq2, GO/KEGG ORA, GSEA, the "v1.5"
significance semantics, the downstream Conda lock, the Nextflow workflows, the R
scripts and every result table and figure are unchanged from v1.5.0; v1.5.x
results do not need to be recomputed.

## Browser GUI

`rnaseq gui` serves a Streamlit GUI on `127.0.0.1` (port 8501 by default) and
uses the same Python services as the CLI. It keeps no database: projects are
`rnaseq new` projects and runs are ordinary `rnaseq run` runs.

- **Dashboard**: projects and recent runs under the configured data folders, the
  mount state of those folders, `rnaseq doctor`, registered references.
- **New Project**: the `rnaseq new` wizard (species, FASTQ or raw counts,
  Salmon or HISAT2 + featureCounts, reference, scope, design, enrichment,
  resources).
- **Input & Metadata**: import a count matrix or a FASTQ folder (files are
  recognised by the validator's naming rule and copied into the project); a
  metadata grid that accepts paste from Excel and can be saved incomplete.
- **Contrasts**: numerator/denominator contrasts with their log2FC direction,
  written by the `rnaseq design` transaction; KEGG snapshot pinning.
- **Run**: validate and plan, then start the run in the background (one active
  GUI run at a time). Runs survive closing the browser and restarting the GUI;
  status and logs are read from the run files.
- **Results**: QC and PCA; per contrast the DEG counts and thresholds, volcano,
  MA plot, heatmap and a searchable DEG table with raw and shrunken log2 fold
  changes and the shrinkage method; GO/KEGG ORA and GSEA tables with the member
  genes of a selected pathway, dotplots, simplified GO dotplots and
  running-score plots; significance semantics, ORA universe, GSEA ranking and
  seed, KEGG snapshot and pin; provenance. Pre-v1.5 enrichment output is
  labelled legacy. Missing outputs are shown as unavailable with their reason.
- **Reports & Downloads**: the run's self-contained `report.html`, the DEG,
  count, QC and enrichment TSVs, PNG/PDF figures, and the delivery package as a
  zip checked against its manifest.

The GUI shows only what the pipeline wrote: figures are the R pipeline's PNG
files and tables are displayed exactly as written.

`report.html` gains Samples, Reference and quantification, LFC shrinkage, and
Software versions and provenance sections, for CLI users too.

## Install

Supported runtime: Linux x86-64 or Windows WSL2 with Conda, Java 17+ and
Nextflow (see `docs/runtime.md`). nf-rna is installed from the Git tag:

```bash
conda create --yes --name nf-rna python=3.11
conda activate nf-rna
python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.0"   # without the GUI: "git+https://github.com/2002brian/nf-rna.git@v1.6.0"
rnaseq --version      # 1.6.0
rnaseq doctor
```

Start the GUI on the server and open it through an SSH port forward:

```bash
tmux new -s nf-rna-gui
rnaseq gui
# on your laptop:
ssh -N -L 8501:127.0.0.1:8501 user@server     # then open http://127.0.0.1:8501
```

The first run of v1.6.0 provisions a new downstream environment from the
unchanged lock (several minutes). See `docs/gui.md` for the GUI settings
(`~/.config/nf-rna/gui.yaml`) and an optional `systemd --user` service.
