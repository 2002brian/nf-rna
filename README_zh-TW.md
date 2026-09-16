[English](README.md) | 繁體中文

# nf-rna

`nf-rna` 是可重現的 bulk RNA-seq workflow，支援 FASTQ 與 raw-count 專案。它會驗證已宣告的輸入與實驗設計，執行選定的 QC 與分析階段，並產生圖表、表格、technical report 與凍結的 provenance。

`rnaseq` 是 control plane，Nextflow 是 execution plane，Docker 是 process runtime。production 使用者透過 `rnaseq` 執行，不需要自行 build execution image 或呼叫內部 Nextflow module。

## 功能

- FASTQ 專案可使用 nf-core/rnaseq 的 Salmon/tximport，或 first-party HISAT2 + featureCounts。
- Raw-count 專案提供已驗證 metadata、明確 contrast、QC、differential expression 與選用的 preranked GSEA。
- Immutable run 會凍結設定、execution-image identity、workflow hash，並建立 curated delivery artifact。

## 系統需求

請使用支援的 Linux 或 WSL workstation，並準備：

- Python 3.11+（含 `venv`）與 Git；
- 已啟動的 Docker daemon；
- 位於 `PATH` 的 Bash、Java 17+ 與 Nextflow；以及
- 足夠儲存 reference 與 Nextflow work data 的本機空間。

Java 與 Nextflow 在 host 執行。R、DESeq2、HISAT2、featureCounts、SAMtools 與 Salmon 由 execution image 提供，一般使用者不需要在 host 安裝。請參考 [Nextflow 官方文件](https://docs.seqera.io/nextflow/install)。

## 安裝

nf-rna 沒有 PyPI distribution。請從 Git tag 安裝已發布的 CLI，並 pull 對應的 official image。將 `X.Y.Z` 換成想使用的已發布版本：

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

release contract 為：

```text
CLI X.Y.Z  ↔  Git tag vX.Y.Z  ↔  ghcr.io/2002brian/nf-rna:X.Y.Z
```

`rnaseq doctor` 是 read-only 檢查；它會驗證 host 與本機 image，但不會 pull 或 build image。

## 快速開始

```bash
mkdir -p ~/projects/rnaseq-projects
cd ~/projects/rnaseq-projects
rnaseq new
cd <new-project>
# 將 FASTQ 放入 input/fastq/，或將 count matrix 放入 input/counts.csv。
# 完成 metadata.csv 與 contrasts.csv。
rnaseq validate .
rnaseq plan .
rnaseq doctor .
rnaseq run . --case-id CASE-001 --profile local --yes
```

wizard 只建立 project scaffold，不會推測科學輸入。FASTQ 專案還需要 execution-ready reference；可用 `rnaseq reference register /absolute/reference-root` 註冊既有的 checksum-bound managed reference，或設定支援的 reference route。詳細內容請見 [Quick Start](docs/quickstart.md)。

## 可重現性

新專案使用與已安裝 CLI 版本相符的 execution image。prerelease CLI 也使用明確相符的 prerelease tag；使用 prerelease 前必須先發布並驗證該 image。

正式分析應使用明確的 release tag 或 immutable digest。`ghcr.io/2002brian/nf-rna:latest` 只是 convenience tag，不應用於進行中的 production analysis。`rnaseq doctor PROJECT` 會報告 requested 與 Docker 實際觀察到的 image identity。

## 文件

- [Detailed Quick Start and v1.0.0 compatibility](docs/quickstart.md)
- [Runtime and reproducibility contract](docs/runtime.md)
- [Scientific contract](docs/scientific_contract.md)
- [Architecture](docs/architecture.md)

## 開發

git clone、editable install、local image build 與 tests 都是 contributor activity。請見 [Development installation](docs/development.md)。
