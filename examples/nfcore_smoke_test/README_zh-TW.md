[English](README.md) | 繁體中文

# nf-core/rnaseq 3.26.0 local smoke fixture

> 本文件為繁體中文翻譯；若與英文版內容有差異，以英文版 README.md 為準。

此 fixture 僅用於 Milestone 2 的 upstream-execution acceptance；它不是 biological study，也不得用於 downstream interpretation。

來源：公開的 `nf-core/rnaseq` release `3.26.0`（`e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`）之 [`conf/test.config`](https://github.com/nf-core/rnaseq/blob/3.26.0/conf/test.config)，並使用其官方 `nf-core/test-datasets` input。

- FASTQ：來自公開 `nf-core/test-datasets/rnaseq/testdata` dataset 的 `GSE110004/SRR6357072_{1,2}.fastq.gz` 與 `GSE110004/SRR6357076_{1,2}.fastq.gz`。僅重新命名為 control plane 慣用的 `sample_R1/R2` 格式。
- Reference：固定 test-dataset revision `626c8fab639062eade4b10747e919341cbf9b41a` 中的 `reference/genome.fasta`、`reference/genes.gtf.gz` 與 `reference/transcriptome.fasta`。prebuilt `salmon.tar.gz` 既未設定也不會發布：其 feature identifier 與已固定 workflow 的 merged SummarizedExperiment step 所用、GTF-derived tx2gene mapping 不相容。

`build_reference.py` 由該官方 transcript FASTA 衍生 `reference/transcriptome.gtf_matched.fasta`，精確保留 official GTF 中的 124 個 transcript ID，並排除未註解的 `Gfp_transgene_gene` record。smoke project 由此 matched derivative 建立 Salmon index。這是在真實執行顯示完整 FASTA 的額外未註解 record 會使 nf-core 的 merged transcript SummarizedExperiment 失敗後，所做的最小 fixture 修正。

此 directory 中的四個 FASTQ 與小型 reference file 是刻意追蹤的公開 test fixture，並非 clinical 或 client data。immutable run、workflow work directory、delivery artifact 與 prebuilt Salmon archive 都由 Git 忽略。若要明確執行 local smoke run：

```bash
rnaseq plan examples/nfcore_smoke_test
rnaseq run examples/nfcore_smoke_test --profile local --yes
```
