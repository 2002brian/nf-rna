# v1.7.0 HISAT2 中間產物生命週期與 upstream 重用修正

- 分支：`fix/v1.7.0-hisat2-artifact-reuse`，基底為 `feature/v1.7.0` 的 `5cba156`
- 日期：2026-11-04
- 觸發案例：NGS1130901_hisat2
  - L1 run `20261104-180729+0800`
  - L2 run `20261104-222119+0800`，啟動參數為 `--reuse-upstream`

## 1. 現象

L1 的 `upstream/` 共 180 個檔案、323.71 GB（以 `find -printf %s` 唯讀量測）。組成如下：

| 類別 | 檔案數 | 大小 |
|---|---|---|
| `alignment/lane_summaries/*.lane.sam` | 8 | 247.34 GB |
| `bam/*.bam`、`*.bai` | 16 | 26.16 GB |
| `bam/count_only/*.countable.bam`、`*.bai` | 16 | 23.07 GB |
| `qc/fastp/prepared/*.fastq.gz` | 16 | 26.14 GB |
| 其餘（counts、summaries、QC、MultiQC） | 124 | 約 0.99 GB |

L2 freeze 期間觀察到的情形：

- 單一 Python 行程同時讀取 L1 的 SAM 並寫出副本。
- 38 分鐘以上，實體讀寫各約 198 GB。
- 沒有任何 HISAT2、Nextflow 或 R 子行程。
- 這個 run 之後於 23:18:41 因 SIGTERM 變成 `INTERRUPTED`（phase=freeze），留下 284 GB 不完整的 `upstream/`。

## 2. 根本原因

### A. 過大的已發佈檔案（workflow）

- `HISAT2_ALIGN` 的輸出 `path("*.lane.sam")` 是傳給 `SORT_LANE_BAM` 的 channel 中介檔。但它的 `publishDir` 沒有 `pattern`，Nextflow 因此把**所有**輸出都複製到 `alignment/lane_summaries/`。
- `FASTP_PREPARE` 也有同樣的問題：`prepared/*.fastq.gz` 一併被發佈到 `qc/fastp/prepared/`。pretrimmed 模式下，這些檔案是原始 FASTQ 的逐位元組副本。
- SAM 在產生 BAM 之後就不再需要。
  - 只有 `SORT_LANE_BAM` 會讀它（經由 channel，在 work dir 內）。
  - 在 `src/`、`workflow/`、R 程式、delivery、status 與 GUI 中，沒有任何程式碼引用 `*.lane.sam`、`lane_summaries` 或 `qc/fastp/prepared`。
- `publishDir` 是唯一的發佈來源：沒有任何 Python 程式碼把這些檔案複製進 run 目錄。

### B. L2 重複複製（Python）

- `rnaseq.service.reuse_upstream_if_compatible`（`--reuse-upstream`）與 `_execute_retry_service_run`（重試一個曾經重用 upstream 的 run）都呼叫 `shutil.copytree(source/"upstream", run/"upstream", dirs_exist_ok=True)`，複製整棵樹。
- 原本的實作沒有任何選擇性處理，只會忽略 `._*` 檔。

L2 實際讀取的內容如下：

| 讀取者 | 檔案 |
|---|---|
| `generate_hisat2_featurecounts_handoff` | `counts/canonical_counts.csv`、`counts/sample_map.csv`、`multiqc/multiqc_report.html`、`qc/strandedness/*.summary` |
| staging 與 delivery | canonical matrix |
| delivery | MultiQC |

handoff 中的 `alignment.bam_directory` 只是一個路徑記錄，沒有任何消費者。

為了可重現性，L2 需要的是「自己持有分析所用的輸入」。這一點只對上述小檔成立。複製 BAM、SAM 與 FASTQ 只是實作上的選擇，並非可重現性的要求：L2 的分析不讀這些檔案，它們的 provenance 屬於產生它們的 L1。

## 3. 修正

### A. Workflow（`workflow/hisat2_featurecounts.nf`）

- `HISAT2_ALIGN`：`publishDir ..., pattern: '*.hisat2.summary'`
- `FASTP_PREPARE`：`publishDir ..., pattern: 'fastp_*'`，只發佈 fastp 的 html 與 json 報告

channel 輸出宣告與所有 script 都沒有改變，因此：

- SAM 與 prepared FASTQ 仍然傳給下游 task，並留在 work dir，由 `rnaseq clean` 清除。
- 沒有新增任何刪除指令。

### B. 選擇性重用（新模組 `src/rnaseq/upstream_reuse.py`）

`copy_reusable_upstream()` 取代上述兩處 `copytree`，行為如下：

- 只在 `upstream/hisat2_featurecounts/` 底下處理以下檔案：
  - `*.sam`、`*.bam`、`*.bai`、`*.cram`、`*.crai`、`*.csi`
  - `qc/fastp/prepared/*.fastq(.gz)`、`*.fq(.gz)`
- 上述檔案**不複製**，改寫入 `referenced_alignment_artifacts.json`。每個項目記錄 path、kind、size_bytes、source_run 與 source_path，並附上 `reused_from` 與 `source_pipeline`（產生這些輸出的 workflow 身分）。
- 其餘檔案以 `shutil.copy2` 複製，並逐檔驗證大小。
- 先在同一層的 `.upstream.partial-<uuid>` 中暫存，完成後才 `os.replace` 到目標位置。任何失敗都會移除暫存目錄，讓 `upstream/` 保持空白；`_run_outcome` 會記錄 `FAILED`，或在收到 signal 時記錄 `INTERRUPTED`。
- 目的目錄非空時拒絕執行，不會合併。
- 鏈式重用（重試一個重用過的 run）時，保留原本的引用與擁有者，不會把擁有者改成中間的 run。
- Salmon（`nfcore_rnaseq`）的樹完全照舊複製。

`generate_hisat2_featurecounts_handoff` 的變更：

- 偵測到引用記錄時，`alignment` 增加 `storage: referenced_in_source_run`、`reference_record`、`reused_from`、`referenced_files` 與 `referenced_bytes`。
- 增加頂層的 `reuse.source_pipeline`。
- 既有欄位全部保留。

`reuse_upstream_if_compatible` 的身分比對完全不變，而且在任何複製之前執行。比對項目為：

- `nfcore.params.json`、`nfcore.local.config` 與 `samplesheet.csv` 逐位元組相同
- input、upstream 與 reference 的身分，包含 FASTQ 的 SHA-256
- 來源 run 必須為 SUCCESS，且有 handoff

### 評估過但未採用的方案

| 方案 | 不採用的原因 |
|---|---|
| hardlink | 跨裝置會失敗；共用 inode 時，任一路徑的寫入都會改變另一個 run 的證據 |
| symlink | 來源 run 封存或刪除後會失效；在 run 目錄內留下指向外部的路徑 |
| reflink | 依檔案系統而定，行為不可預期 |

## 4. 相容性

| 情境 | 行為 |
|---|---|
| 修正前的 L1（含 lane SAM 與 prepared FASTQ） | 可重用；SAM、BAM 與 FASTQ 被引用，不複製。`hisat2_featurecounts.nf` 的 SHA-256 不在重用比對中 |
| 修正後的 L1 | 不再發佈 SAM 與 prepared FASTQ |
| 修正前完整複製的 L2 作為重試來源 | 引用其中的副本，擁有者記為該 L2 |
| 修正後的 L2 作為重試來源 | 保留原本 L1 的引用 |
| Salmon 重用 | 不變 |
| 來源 run 被封存或刪除 | L2 的分析與 delivery 不受影響（從未讀取這些檔案）；引用記錄的 source_path 會失效 |
| `upstream_handoff_manifest.yaml` 的讀取端 | 只增加欄位，沒有移除 |

## 5. 驗證

### 新增的測試：`tests/test_hisat2_artifact_reuse.py`

Workflow 契約：

- 兩個 publishDir 的 `pattern`
- 沒有未過濾的 SAM 發佈
- 必要輸出全部保留
- 科學指令列逐字不變；與 `5cba156` 的差異只在 publishDir 與註解

真實 Nextflow 執行：

- 使用替身工具，分 pretrimmed 與 fastp 兩種模式。
- 確認沒有發佈任何 `.sam` 或 `.fastq.gz`，所有必要輸出都存在，且 SAM 仍然傳到 sort。

選擇性複製：

- 修正前與修正後的 L1
- Salmon 照舊完整複製
- 鏈式重用
- 修正前的完整副本作為來源
- 跨裝置（EXDEV）、磁碟滿（ENOSPC）與短寫入時沒有殘留
- 非空目的目錄、缺少來源與壞掉的記錄都會被拒絕

服務層（只模擬 Nextflow）：

- L1 → design → L2 重用成功，包含修正前與修正後的 L1
- L2 的 delivery `raw_counts.csv` 與 staged matrix 和 L1 逐位元組相同
- 身分不符時，在複製前拒絕，且沒有產生任何檔案
- 來源缺少 canonical matrix 時為 `FAILED`
- 複製途中 ENOSPC 時為 `FAILED`，且 `upstream/` 為空
- 重試一個失敗的重用 L2 後為 `SUCCESS`，引用鏈指回 L1

### 測試結果

| 測試 | 結果 |
|---|---|
| 全套測試（第一個 commit，使用鎖定的 downstream R prefix 與 Nextflow 25） | 1362 passed、3 skipped |
| 加入 FASTP 修正後的全套測試 | 1364 passed、3 skipped |
| `test_featurecounts_real_semantics.py` 與 pinned conda 測試（使用真實 SAMtools 1.21 與 featureCounts 2.0.6） | 3 passed |

## 6. 效能

### 實測

這是一個合成的 8 樣本 fixture，各類檔案的比例依照真實 run，總大小 0.96 GB。腳本為 `nf-rna-validation/hisat2-artifact-reuse/measure_reuse.py`，跑了 3 次。

| 方法 | 複製檔案數 | 複製位元組 | 放大倍數 | 時間 | rchar / wchar |
|---|---|---|---|---|---|
| 原本的 `copytree` | 115 | 958,541,533 | 1.000 | 0.349–0.357 s | 958.5 MB / 958.5 MB |
| 選擇性複製 | 60 | 23,678,023 | 0.025 | 0.019–0.021 s | 23.7 MB / 23.7 MB |

時間數字受 page cache 影響，只能當作相對比較。

### 推估（非實測）

NGS1130901 的 L1 若以修正後的程式碼重用：

- 不複製的部分：322.72 GB（SAM、BAM、FASTQ）
- 實際複製的部分：約 0.99 GB
- 複製量減少約 99.7%

以觀察到的約 87 MB/s 估算：

- 原本的 L2 freeze 約 62 分鐘
- 修正後約 10–15 秒

修正後的新 L1 不再發佈 SAM 與 prepared FASTQ，`upstream/` 預期約 50 GB，主要是 BAM。

## 7. 是否阻擋 v1.7.0 發佈

建議**列為發佈前必修**。理由：

- 這不是科學錯誤，所有數值輸出都不變。
- 但它讓每次 HISAT2 L1 多佔約 270 GB 磁碟。
- 每次 L2 重用都多一次等量的讀寫，而且會長時間停在 D 狀態，實際上讓 L2 很難完成。這次觀察到的 run 已被中斷。
- 修正範圍小，只有發佈過濾與複製選擇，而且已有回歸測試覆蓋。

## 8. 剩餘風險與遷移

- 既有 run 中的 `*.lane.sam` 與 `qc/fastp/prepared/` 不會被自動刪除。確認不需要後，可以手動刪除；建議先確認 BAM 完好。
- 中斷的 L2 `20261104-222119+0800` 留有 284 GB 不完整的副本。不建議 `rnaseq retry` 它：它的 contract 中 `upstream_handoff` 為 null，重試會被視為「不重用」而重跑整個 HISAT2 upstream。應改用修正後的程式碼，重新啟動一個 `--reuse-upstream` 指向 L1 的新 L2。
- 引用記錄只有大小，沒有 SHA-256。對數百 GB 計算雜湊等於重新讀取所有資料，與修正的目的相違。L2 的科學結果並不依賴這些檔案。
- 重用的 L2 的 handoff `pipeline` 仍然是目前安裝的 workflow 身分（既有行為）。產生輸出的 workflow 身分改記在 `reuse.source_pipeline`。
- 自訂（非受管）參考的 HISAT2 index 是目錄，這類 run 無法 `rnaseq retry`（既有限制，與本修正無關）。
