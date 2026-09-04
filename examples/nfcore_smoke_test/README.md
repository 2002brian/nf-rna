English | [繁體中文](README_zh-TW.md)

# nf-core/rnaseq 3.26.0 local smoke fixture

This fixture is for Milestone 2 upstream-execution acceptance only. It is not a biological study and must not be used for downstream interpretation.

Source: the public `nf-core/rnaseq` release `3.26.0` (`e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`) [`conf/test.config`](https://github.com/nf-core/rnaseq/blob/3.26.0/conf/test.config), using its official `nf-core/test-datasets` inputs.

- FASTQs: `GSE110004/SRR6357072_{1,2}.fastq.gz` and `GSE110004/SRR6357076_{1,2}.fastq.gz` from the public `nf-core/test-datasets/rnaseq/testdata` dataset. They are renamed only to the control-plane's conventional `sample_R1/R2` form.
- Reference: the fixed test-dataset revision `626c8fab639062eade4b10747e919341cbf9b41a`: `reference/genome.fasta`, `reference/genes.gtf.gz`, and `reference/transcriptome.fasta`. A prebuilt `salmon.tar.gz` is neither configured nor published: its feature identifiers are incompatible with the GTF-derived tx2gene mapping in the pinned workflow's merged SummarizedExperiment step.

`build_reference.py` derives `reference/transcriptome.gtf_matched.fasta` from that official transcript FASTA, retaining exactly the 124 transcript IDs in the official GTF and excluding the unannotated `Gfp_transgene_gene` record. The smoke project builds its Salmon index from this matched derivative. This is the minimal fixture correction after real runs showed that the full FASTA's extra unannotated record made nf-core's merged transcript SummarizedExperiment fail.

The four FASTQs and small reference files in this directory are deliberately tracked public test fixtures. They are not clinical or client data. Immutable runs, workflow work directories, delivery artifacts, and the prebuilt Salmon archive are ignored by Git. For an explicit local smoke run:

```bash
rnaseq plan examples/nfcore_smoke_test
rnaseq run examples/nfcore_smoke_test --profile local --yes
```
