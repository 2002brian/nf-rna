[English](README.md) | 繁體中文

# nf-rna

`nf-rna` 是可重現的 bulk RNA-seq workflow。`rnaseq` 是 control plane；Nextflow 是 execution plane；Docker 是 process runtime。一般 production 使用者應透過 `rnaseq` 執行，不需要自行 build execution image 或呼叫內部 Nextflow module。

目前已發布的 release 是 `v1.0.0`。本 branch 的 installation UX 變更是為下一個 patch release 準備；本文件不宣稱已有新的 release。

## 下一個已發布 patch release 的 canonical 安裝路徑

下一個包含此變更的 patch release 發布後，請將下列 `RELEASE_VERSION` 換成實際已發布的版本；不要據此假設 `v1.0.1` 已存在。

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

此為 non-editable Git 安裝；目前沒有 PyPI distribution。CLI 會將自己的 package version 映射為新專案的預設 image：CLI `X.Y.Z` 對應 `ghcr.io/2002brian/nf-rna:X.Y.Z`。例如 prerelease `1.0.1rc1` 需要明確對應已發布的 `:1.0.1rc1` image，絕不以 `latest` 取代。

Linux/WSL host 需要 Python 3.11+（含 `venv`）、Git、已啟動的 Docker daemon、Bash、Java 17+ 與 Nextflow。Java/Nextflow 在 host 執行；R、DESeq2、HISAT2、featureCounts、SAMtools 與 Salmon 在 execution image 中。請參考 [Nextflow 官方文件](https://docs.seqera.io/nextflow/install)。

## 下一個 release 的 canonical Quick Start

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

這個路徑不需要手動替換 YAML。FASTQ 專案仍需要 execution-ready reference；詳細設定請見 [quick start](docs/quickstart.md)。

## 歷史 v1.0.0 相容性

不可變的 `v1.0.0` wizard 仍會寫入 `execution_image: nf-rna:latest`。它可使用已驗證的 image，但在 plan/run 前需要一次相容性替換：

```bash
python3.11 -m venv ~/.venvs/nf-rna-1.0.0
source ~/.venvs/nf-rna-1.0.0/bin/activate
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.0.0"
docker pull ghcr.io/2002brian/nf-rna:1.0.0
sed -i 's|execution_image: nf-rna:latest|execution_image: ghcr.io/2002brian/nf-rna:1.0.0|' project.yaml
```

這只屬於 v1.0.0 historical compatibility，不是下一個 patch release 的正常 Quick Start。

## Runtime identity 與文件

`rnaseq doctor PROJECT` 會報告 requested image 與 Docker 實際觀察到的 image identity。正式分析應使用 version tag 或 immutable digest；`:latest` 僅是 convenience tag。

- [Detailed Quick Start](docs/quickstart.md)
- [Runtime and reproducibility](docs/runtime.md)
- [Scientific contract](docs/scientific_contract.md)
- [Architecture](docs/architecture.md)
- [Development installation](docs/development.md)

git clone、editable install、local Docker build 與 tests 都是 contributor workflow，已與一般 production installation 分離。
