# nf-rna 獨立雲端程式碼稽核報告（2026-11-04）

- 稽核對象：`2002brian/nf-rna`，基準 commit `c75e962`（`develop`，v1.4.0 之後）
- 稽核方式：靜態閱讀 + 以合成資料（synthetic fixtures）實際重現；未存取任何醫院伺服器、原始 FASTQ 或私有資料
- 執行環境限制：雲端容器有 Python 3.13，但**沒有** R/DESeq2、Nextflow、Conda、HISAT2、featureCounts。因此 R 與 Nextflow 端的結論來自程式碼閱讀，並在下文標示「靜態確認」；Python 端的結論都有可重跑的測試佐證
- 本報告附帶的變更**只有**文件與 `tests/audit/` 下的獨立測試，未修改任何正式程式、GUI、繪圖程式或科學預設值

## 0. 重要前提：本 repo 沒有 Streamlit GUI

在 `c75e962` 中搜尋 `streamlit`，零結果。GUI 程式碼目前應該只存在於本地開發環境，尚未提交。因此第 5 節的「GUI 架構風險」只能從 GUI 必然會呼叫的 service/CLI 契約推論，不是對 GUI 程式碼本身的稽核。

## 1. 測試基準

```bash
python -m pytest -q -p no:cacheprovider
```

| 狀態 | 數量 |
|---|---|
| 稽核前 | 838 passed, 8 failed, 186 skipped |
| 加入 `tests/audit/` 後 | 843 passed, 8 failed, 186 skipped, 14 xfailed |

8 個失敗都是**環境耦合**造成，不是產品邏輯錯誤（詳見 T-01）：

- 7 個失敗是因為宿主沒有 `conda` 或 `nextflow` 可執行檔，涵蓋 `test_execution.py` 4 個、`test_retry.py` 1 個、`test_runtime_dispatch.py` 2 個
- `test_cleanup.py::test_symlink_loop_is_refused_without_traceback` 在 Python 3.13 失敗。程式仍正確拒絕了 symlink，只是錯誤訊息改走 `src/rnaseq/cleanup.py:273` 的分支；CI 只跑 3.11，所以沒被發現

`tests/audit/` 的慣例：已確認的缺陷寫成 `xfail(strict=True)`，斷言內容是「安全的預期行為」。缺陷修好後該測試會 XPASS，strict 模式會把它判為失敗，迫使有人刻意移除 marker，作為修復的驗收點。

---

## 2. 已確認缺陷（Confirmed defects）

### D-01【高】L1 → L2 `--reuse-upstream` 在常見路徑上必然被拒絕

- **位置**：
  - `src/rnaseq/service.py:268` 把 `sha256_source`（`computed` 或 `cache`）寫進凍結的 `input.files[*]`
  - `src/rnaseq/service.py:2131` 的 `_upstream_identity` 比較的是整個 `input` 區段
  - `src/rnaseq/service.py:2167` 的不相等判斷因此拋出 `incompatible`
- **機制**：
  1. 第一次 run（典型情況是 L1）會計算 SHA256，記錄為 `sha256_source: computed`
  2. 第二次 run（`rnaseq design` 之後的 L2）從 checksum cache 取值，記錄為 `cache`
  3. 兩次的路徑、大小、mtime、SHA256 完全相同，只有這個簿記欄位不同，就被判定為上游不相容
- **影響**：CHANGELOG 1.4.0 主打的「L1 → design → L2 重用上游」在實務上的第一次重用就會失敗，使用者只能重跑整段 Salmon 或 HISAT2，對大型 FASTQ 代價很高。科學結果不會出錯，屬於 fail-closed。
- **為何既有測試沒抓到**：`tests/test_service.py:946-993` 使用手寫的 `_UPSTREAM_MANIFEST`，裡面沒有 `sha256_source`、`resolved_path` 和 `identity_contract`。
- **證據**：`tests/audit/test_audit_reuse_checksum_source.py`。以 `--runxfail` 執行時，差異正好落在 `input.files`。
- **建議修補**（最小、無科學影響）：

  ```python
  # service.py _upstream_identity
  def _identity_input(section):
      if not isinstance(section, dict):
          return section
      files = [{k: v for k, v in f.items() if k != "sha256_source"} if isinstance(f, dict) else f
               for f in section.get("files", [])]
      return {**section, "files": files}
  ...
  "input": _identity_input(manifest.get("input")),
  ```

  另一種做法是把 `sha256_source` 移出 `input`，改放到 provenance。

### D-02【中】`contrast_id` 未受限制，直接變成 R 輸出路徑的一部分

- **位置**：
  - `src/rnaseq/validators.py:577` 起只檢查非空與唯一
  - `src/rnaseq/r/l2_analysis.R:106` 使用 `file.path(cfg$output_dir, "contrasts", contrast_id)` 並 `dir.create(recursive=TRUE)`
  - `src/rnaseq/workflow_support.py:739` 以同樣方式讀回
- **來源**：只有 `rnaseq design` 產生的 ID 會經過 `CONTRAST_ID_PATTERN` 檢查（`src/rnaseq/design_setup.py:35`）。手寫的 contrasts.csv 或 `rnaseq new --contrasts` 匯入的檔案都不會檢查。
- **重現**：`contrast_id = ../../escape`、`a/b`、`with space` 都能通過 `validate`。
- **影響**：輸出可能寫到 `contrasts/` 之外，報告和 delivery 找不到檔案或放錯位置，破壞 run 目錄的封閉性。這是本機自傷型問題，不是遠端攻擊面。
- **證據**：`tests/audit/test_audit_validation_gaps.py::test_unsafe_contrast_id_is_rejected`
- **建議**：在 `validate_contrasts` 對 `contrast_id` 套用與 `CONTRAST_ID_PATTERN` 相同的規則，產生 `unsafe_contrast_id` error。

### D-03【中】設計變數值為字面 `NA` 時，Python 視為合法 level，R 視為缺值

- **位置**：
  - `src/rnaseq/validators.py:514` 只拒絕空白
  - `src/rnaseq/r/design_metadata.R:16` 的 `read.csv(..., colClasses="character")` 預設 `na.strings="NA"`，檔案第 4 行的註解也明說只有字面 NA 會被當成缺值
- **重現**：把 `Treatment` 改成 `NA`，並設定 contrast `NA_vs_Control`，結果 `validate` 仍然 PASS。
- **影響**（靜態確認）：R 端 factor 會出現 NA，`DESeqDataSetFromMatrix` 在 L1 會因設計變數含 NA 而失敗。FASTQ 專案要到上游跑完數小時後才會報錯。
- **證據**：`test_literal_na_design_level_is_flagged`
- **建議**：Python 端對設計變數與 `pair_id` 的字面 `NA` 回報 error。或者 R 端改用 `na.strings=character()`，但這會改變既有語意，需要評估。

### D-04【中】`sample_id = "NA"` 在 R 端變成 NA rowname

- **位置**：`src/rnaseq/r/design_metadata.R:16-18` 的 `rownames(metadata) <- metadata$sample_id`
- **重現**：counts 和 metadata 都把 `T3` 改成 `NA`，結果 PASS。
- **影響**（靜態確認）：R 會拋出 `missing values in 'row.names' are not allowed`，L1 失敗。
- **證據**：`test_literal_na_sample_id_is_flagged`
- **建議**：與 D-03 一起處理。

### D-05【中】Python 數值語法比 R 寬鬆，驗證通過的 counts 在 R 端變成 NA

- **位置**：
  - `src/rnaseq/validators.py:348` 的 `Decimal(raw_value)`
  - `src/rnaseq/validators.py:522` 的 `float(...)`，用於連續變數
  - `src/rnaseq/r/l1_analysis.R:16` 與 `src/rnaseq/r/l2_analysis.R:21` 的 `read.csv`
- **機制**：Python 接受 PEP 515 底線寫法（`1_000`）和非 ASCII 數字（`١٢`）。R 的 `read.csv` 會把該欄讀成字元，`storage.mode <- "numeric"` 之後變成 NA，DESeq2 再拒絕含 NA 的 counts。連續設計變數同理，例如 `as.numeric("1_000")` 會得到 NA。
- **證據**：`test_count_values_r_cannot_parse_are_rejected[1_000 / ١٢]`。對照測試 `4e1` 仍然合法。
- **建議**：數值欄位先用 ASCII regex `^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$` 檢查格式，再轉成 Decimal 或 float。

### D-06【中，HISAT2 E2E 相關】FASTQ 檔名衍生的 sample ID 未加引號插入 Nextflow shell script

- **位置**：
  - `src/rnaseq/validators.py:38` 的 `(?P<sample>.+?)` 允許任意字元
  - `workflow/hisat2_featurecounts.nf` 的 `MERGE_AND_INDEX`（第 169 行起）、`PREPARE_COUNT_BAM`、`FEATURECOUNTS`（第 193 行）都直接使用 `${sample}`，沒有加引號
- **重現**：`a b_R1.fastq.gz`、`a$(id)_R1.fastq.gz`、`a;b_…`、`a'b_…` 全部被接受為 sample ID。
- **影響**：
  - 含空白的 ID 會讓 `samtools merge -o a b.merged.bam` 失敗。這發生在該樣本完成比對之後，等於白跑
  - 含 `$()` 或 `;` 的 ID 會在 task shell 中執行命令
  - Salmon 路徑由 nf-core 的 samplesheet schema（不允許空白）擋下一部分；HISAT2 第一方路徑沒有對應的防線
- **證據**：`test_shell_unsafe_fastq_sample_ids_are_rejected[...]`
- **建議**：
  - 在 `validate_fastq` 與 metadata 驗證中限制 `sample_id` 為 `^[A-Za-z0-9][A-Za-z0-9._-]*$`。這會讓含空白的既有專案驗證失敗，屬於刻意的破壞性收緊，需要在 CHANGELOG 中說明
  - nf script 改用 `"${sample}"` 加引號作為縱深防禦

### D-07【低到中，provenance】delivery 內的 `run_state_*.json` 永遠是 `RUNNING`

- **位置**：
  - `src/rnaseq/service.py:1655` 在 `assemble_delivery` 中複製 `run.state_path`
  - 呼叫順序在 `src/rnaseq/service.py:2288-2290`、`2309-2311`、`2418-2420`、`2440-2442`，四處都是先寫 `RUNNING/phase=delivery`，再 assemble，最後才寫 `SUCCESS`
- **影響**：客戶收到的交付包記錄的是未完成狀態，`completed_at` 為 null，與 run 目錄內的最終狀態不一致。這不影響 `delivery_manifest.yaml` 的雜湊自洽。
- **證據**：`tests/audit/test_audit_delivery_state.py::test_delivered_run_state_is_final`（`'RUNNING' == 'SUCCESS'` 失敗）。對照測試確認 run 本身是 SUCCESS，且 `verify_delivery_manifest` 通過。
- **建議**：兩種做法擇一：
  - 不交付 run_state，只保留 manifest 與 provenance
  - delivery 中寫入一份明確標為 `delivery_snapshot` 的狀態，記錄 `status_at_packaging: RUNNING` 與說明，避免與最終狀態混淆

  不建議在 manifest 寫完之後再改動 delivery 內容。

### D-08【低到中，QC 報告】L1 的 PCA、箱型圖、相關圖分組變數與 Python 的 condition factor 不一致

- **位置**：
  - `src/rnaseq/r/l1_analysis.R:60` 的 `color_name <- tail(all.vars(formula), 1)`
  - Python 端是 `src/rnaseq/validators.py:679`：formula 含 `condition` 時優先用它
- **機制**：自動產生的 formula 都把 `condition` 放在最後（`src/rnaseq/project.py:154`、`src/rnaseq/cli.py:782`），所以預設專案不受影響。但使用者自訂 `~ condition + batch` 時（`rnaseq design` 會保留自訂 formula），L1 會用 batch 著色並標為 group；`~ condition + age` 則會把連續變數 age 當成 factor 著色。
- **證據**：靜態確認。容器內沒有 R，所以沒有寫可執行測試。
- **建議**：由 control plane 在 downstream contract 中凍結 `primary_factor`，R 端直接讀取，不再自行推斷。這只影響圖的分組，不影響任何統計量。由於涉及繪圖程式，請交給本地 session 處理。

---

## 3. 高優先風險（High-priority risks，尚非已證實缺陷）

| ID | 風險 | 位置 | 說明與建議 |
|---|---|---|---|
| R-01 | 狀態檔與合約檔的寫入**非原子** | `src/rnaseq/service.py:180-182`（`write_text`）、`src/rnaseq/service.py:675-699`（`_write_state`）、`src/rnaseq/execution.py:1111` | SIGKILL、OOM 或斷電可能留下截斷的 `run_state.json`。<br>`_write_state` 遇到 `JSONDecodeError` 會**靜默重置為 `{}`**（`service.py:680-681`），遺失 `started_at`、`command` 等欄位。<br>`load_run_states` 會略過壞檔（`execution.py:1568` 起），run 因此從清單中消失。<br>`status --watch` 和未來的 GUI 輪詢可能讀到半寫的檔案。<br>建議：同目錄 tmp → `fsync` → `os.replace`。checksum cache（`service.py:243-246`）已經這樣做，可直接沿用 |
| R-02 | 凍結輸入的 TOCTOU | `src/rnaseq/service.py:749-756` 複製當下的專案檔，`:758` 再對當下的檔案重算雜湊 | `execute_service_run` 會再執行一次 `require_fresh_plan`，所以互動確認期間的修改會被擋下；剩餘時間窗只有毫秒級。<br>但驗證用的 report（設計型別、pairing）與實際凍結的檔案之間，沒有「凍結副本雜湊 = 計畫雜湊」的直接斷言。<br>建議：複製後以凍結副本計算雜湊，並與 `planning/manifest.preview.yaml` 比對 |
| R-03 | 同一樣本同時有「合併 FASTQ」與「各 lane FASTQ」時會靜默重複計數 | `src/rnaseq/validators.py:418` 把沒有 lane 的檔案標為 `L000` 並與 `001` 並存 | 實驗室常見情況：資料夾裡同時放了 `S1_R1.fastq.gz`（合併檔）和 `S1_L001_R1_001.fastq.gz`。目前會兩份都比對、合併，**reads 加倍**。這是有科學意義的錯誤，因為 library size 與離散度估計都會受影響。<br>證據：`test_merged_and_lane_fastqs_for_one_sample_are_flagged`（xfail）。<br>建議：同一 sample 同時出現 `L000` 與編號 lane 時回報 error 或 warning。lane token 也應統一，目前混用 `001` 和 `L000` |
| R-04 | FASTQ checksum cache 以 `size + mtime_ns` 為鍵 | `src/rnaseq/service.py:211-248` | 已寫入文件的設計取捨。`cp -p`、`rsync -t` 或保留時間戳的原地修復會讓舊的 SHA256 被沿用。建議在 production acceptance 模式提供 `--rehash-fastq`，或定期全量重算 |
| R-05 | **HISAT2 acceptance 證據已過期** | `docs/runtime.md:263`、`:273` | 驗收紀錄的 workflow SHA 是 `e5408f…` 與 `ead01d…`（2026-09-05）。目前 `workflow/hisat2_featurecounts.nf` 的 SHA256 是 `b62c0e6c…`，最後一次修改是 `d4bc48e`（2026-09-23，WSL2 Conda runtime）。現行 Conda runtime 版本的 workflow 沒有任何 E2E 驗收紀錄，見第 7 節 |
| R-06 | Groovy slashy regex `/\\.fastq$/` 很可能永遠不會匹配 | `workflow/hisat2_featurecounts.nf:73-74`、`:109` | 在 slashy string 中，`\\.` 是「字面反斜線 + 任意字元」。檔名會因此變成 `X_R1.fastq.fastq.gz`、`X_R1.fastq.fastq.lane.sam`，只影響外觀，不影響計數。**需要用 Nextflow 實跑確認**。建議改為 `/\.f(ast)?q$/` |
| R-07 | `ASSEMBLE_COUNTS` 用 `.replace('.counts','')` 從檔名取回 sample ID | `workflow/hisat2_featurecounts.nf:233` | Groovy 的 `replace` 會替換所有出現的子字串，sample `rna.counts1` 會變成 `rna1`。之後的 `_stage_featurecounts_matrix`（`service.py:398-403`）會因樣本集合不一致而大聲失敗，不會靜默出錯。建議改用 `it.name.replaceFirst(/\.counts\.txt$/, '')`，或者直接從 tuple 傳入 sample |
| R-08 | 自訂 reference 依 basename 凍結 | `src/rnaseq/service.py:787` | 兩個不同目錄下同名的檔案（例如兩個 `genome.fa`）會互相覆蓋。建議以 key 作為子目錄 |
| R-09 | 中斷處理期間第二次 Ctrl-C | `src/rnaseq/service.py:1268-1278` | 等待 300 秒期間若再收到 `KeyboardInterrupt`，會跳過 `kill()`。由於使用了 `start_new_session=True`，Nextflow 會成為孤兒程序。建議在 `finally` 中保證 `kill` 並記錄 PID |
| R-10 | 沒有專案層級的鎖 | `rnaseq design` 與 `rnaseq run` 之間 | `create_case_run` 以 `mkdir` 保證 run ID 唯一（`service.py:607-629`），但兩個 CLI 或 GUI 行程可能同時改寫 `planning/` 和 `contrasts.csv`。GUI 上線後風險會升高（見第 5 節） |

## 4. 科學可重現性與方法一致性觀察

以下項目已逐一檢查並判定為**正確或已有文件記載**，不屬於缺陷：

- **Strandedness**：HISAT2 路徑由 `models.py:542-545` 拒絕 `auto`。`hisat2_strand_option` 與 `featurecounts_arguments` 的對應正確：forward 對應 FR/F 與 `-s 1`，reverse 對應 RF/R 與 `-s 2`。`MERGE_AND_INDEX` 會檢查同一樣本的各 lane strandedness 一致。
- **featureCounts 計數政策**：`-t exon -g gene_id -Q 0 --primary`；雙端資料加 `-p --countReadPairs -B -C`；另外以 `-F 0x900` 產生只用於計數的 BAM，保留原始 BAM 作為診斷用。與 `scientific_contract.md` 一致。
- **矩陣組裝**：`hisat2_featurecounts.py` 要求每個檔案恰好一個計數欄、基因集合一致、計數為非負整數。staging（`service.py:379-415`）依凍結 metadata 的順序重排，並驗證樣本集合完全一致。
- **Salmon / tximport**：`sample_inputs.R` 以 sample 名稱而非位置選取 `quant.sf`，並驗證 `colnames(txi$counts)` 與 metadata 一致。
- **配對設計**：每個 pair 必須在 numerator 和 denominator 各恰有一個樣本，L2 至少需要 2 個完整 pair。rank 檢查（`validators.py:1002-1095`）使用與 R 相同的加法 model matrix；level 排序不同不影響 rank。
- **L2**：只 fit 一次 DESeq2，各 contrast 明確以 `c(factor, num, den)` 取出；不做 shrinkage，保留 independent filtering。significance 判準為 `padj < α` 且 `|log2FC| ≥ t`，與文件一致。
  - `l2_analysis.R:17` 會以字母順序重新 `factor()` contrast factor，這只影響 provenance 中記錄的係數名稱，不影響任何結果。
- **L1/L2 基因宇集**：L2 以相同規則**重新計算**過濾結果（`l2_analysis.R:36-39`），而不是讀取 L1 保留的基因清單。只要輸入相同，兩者必然相等，但目前**沒有斷言**。建議 L2 讀取 L1 的 `backend_summary.json` 或基因清單並比對，失敗時停止。
- **VST / PCA**：
  - 少於 1000 個基因時改用精確的 `varianceStabilizingTransformation`，已記錄在 summary
  - PCA 使用**所有保留基因**，而 DESeq2 的 `plotPCA` 預設只用前 500 個高變異基因。這是有文件記載的方法選擇，但與其他工具的結果比較時，必須說明這個差異
  - 相關係數是 blind VST 的 Pearson 相關
- **Salmon 的 size factor**：tximport 路徑的 `sizeFactors` 為 NULL，報表中記為 NA，並以 `scalar_size_factors_available` 標示，處理正確。

## 5. Streamlit GUI 架構風險（從契約推論）

由於 GUI 程式碼不在 repo 中，以下是本地 GUI 實作時應避免的陷阱：

1. **輪詢讀取非原子寫入的檔案**（R-01）：GUI 讀取 `run_state.json`、`execution.log` 或 trace 時，要容忍 `JSONDecodeError` 並重試，不能把它當成「沒有 run」。最好的做法是先修 R-01。
2. **在 Streamlit 行程內呼叫 `execute_service_run`**：
   - `_terminate_on_signals` 只在 main thread 安裝 signal handler（`service.py:1234`），Streamlit 的 script thread 不是 main thread，因此 SIGTERM 和 SIGHUP **不會**被轉成 INTERRUPTED
   - Streamlit 每次 rerun 都會重新執行整個腳本，有重複啟動 run 的風險

   建議 GUI 一律用 `subprocess` 啟動 `rnaseq run --yes`，並以 `run_status` 的 process record 判斷 run 是否存活。
3. **驗證結果快取**：`st.cache_data` 若以檔案路徑為鍵，metadata 修改後會回傳舊的 `ValidationReport`。鍵應包含 `size + mtime_ns`，或直接使用 `render_manifest` 的雜湊。
4. **在 GUI 中編輯 metadata 或 contrasts 時**，要經過與 CLI 相同的 `validate_project` 與 `rnaseq design` 交易流程，並使用原子寫入。必須拒絕 D-02 到 D-06 列出的值；在正式程式修好之前，GUI 端也可以先擋。
5. **GUI 與 CLI 一致性**：GUI 不應自行推斷 primary factor（見 D-08），應直接使用 `condition_factor()` 或凍結合約中的值。

## 6. 缺少的測試（Missing tests）

| ID | 缺口 | 本次是否補上 |
|---|---|---|
| T-01 | 測試不夠 hermetic：7 個測試依賴宿主是否有 conda/nextflow，1 個依賴 Python 版本的錯誤訊息。CI 只在 3.11 跑，而且只在 `main` 與 PR 上觸發 | 未補，屬於正式測試修改，建議以 monkeypatch `check_upstream_conda` 和 `check_nextflow` 解決 |
| T-02 | reuse 與真實凍結 manifest 的整合測試（D-01） | ✅ `tests/audit/test_audit_reuse_checksum_source.py` |
| T-03 | Python 與 R 解析語法的對等性（D-03/04/05） | ✅ `tests/audit/test_audit_validation_gaps.py` |
| T-04 | 路徑與 shell 安全的識別字（D-02/06） | ✅ 同上 |
| T-05 | 合併檔與各 lane 檔重複（R-03） | ✅ 同上，含正常多 lane 的對照測試 |
| T-06 | 交付快照與最終狀態一致（D-07） | ✅ `tests/audit/test_audit_delivery_state.py` |
| T-07 | 截斷的 `run_state.json` 復原行為（R-01） | 未補：目前行為是靜默重置，需要先決定預期語意 |
| T-08 | L1 primary factor 與 Python 一致（D-08） | 未補：需要 R，且涉及繪圖程式 |
| T-09 | L2 宇集等於 L1 宇集的斷言 | 未補：需要 R |
| T-10 | HISAT2 workflow 的 Nextflow stub-run（`-stub`）或 nf-test，涵蓋檔名、單雙端、多 lane、`.fq.gz` | 未補：需要 Nextflow |

## 7. HISAT2 + featureCounts E2E 準備度評估

**結論：條件式就緒（Conditionally ready）。科學計數政策正確，但目前的 workflow 版本沒有 E2E 證據，且有兩個會在比對之後才觸發的失敗模式。**

| 面向 | 評估 |
|---|---|
| 計數政策與 strandedness | ✅ 正確，並有 `test_featurecounts_*semantics` 測試 |
| 參考資料契約（checksum、runtime splice sites、`requires_smoke_validation` 閘門） | ✅ 嚴謹；未驗證的版本組合不視為 execution-ready |
| 多 lane 合併與樣本身分 | ✅ 依 sample ID 合併並檢查 strandedness 一致；⚠️ R-03 會造成重複計數 |
| 識別字安全 | ❌ D-06：含空白或 shell 字元的檔名會在比對後失敗或被注入 |
| 上游重用（L1 → L2） | ❌ D-01：常見路徑上被拒絕 |
| 現行 workflow 的驗收紀錄 | ❌ R-05：最後一次驗收早於 Conda runtime 的修改 |
| 檔名與樣本對應 | ⚠️ R-06 / R-07，失敗時會明確報錯 |

**放行前的建議順序**：

1. 修 D-01 與 D-06。
2. 合成資料煙霧測試：使用 `tests/fixtures/hisat2_paired_raw_v2`，在現行 Conda runtime 上重跑一次 paired QC 與一次 raw L2，並把新的 workflow SHA 寫進 `docs/runtime.md`。
3. 跑一次 L1 → `rnaseq design` → L2 的 `--reuse-upstream` 實測，確認 D-01 修好。
4. 加入 `nextflow run -stub` 的 CI 工作，至少檢查 R-06 / R-07 的檔名行為。

## 8. 優先行動計畫

| 優先 | 項目 | 範圍 | 與本地 GUI 或繪圖工作衝突？ |
|---|---|---|---|
| P0 | D-01 reuse 身分排除 `sha256_source` | `service.py` 約 10 行 | 否 |
| P0 | D-06 sample ID 字元限制，nf script 加引號 | `validators.py`、`hisat2_featurecounts.nf` | 否 |
| P1 | D-02 contrast_id 限制；D-03/04/05 改為 R 等價的解析規則 | `validators.py` | 否 |
| P1 | R-01 原子寫入與壞檔處理 | `service.py`、`execution.py` | 否，但 GUI 依賴它 |
| P1 | R-03 合併檔與 lane 檔重複偵測 | `validators.py` | 否 |
| P1 | R-05 現行 workflow 的 HISAT2 E2E 驗收 | 執行加文件 | 否 |
| P2 | D-07 交付快照語意 | `service.py` | 否 |
| P2 | D-08 凍結 `primary_factor` | contract 加 `l1_analysis.R` | **是（繪圖）**，請交給本地 session |
| P2 | T-01 測試 hermetic 化；CI 加 3.13；`develop` push 觸發 CI | `tests/`、`.github/` | 否 |
| P3 | R-02、R-04、R-06 到 R-10 | 各處 | 否 |

## 9. 重現方式

```bash
pip install -e '.[dev]'
```

```bash
python -m pytest -q -p no:cacheprovider tests/audit -rxX
```

```bash
python -m pytest -q -p no:cacheprovider tests/audit --runxfail
```

第二個指令看到 14 個 xfailed 就是預期結果。第三個指令會顯示每個缺陷的實際失敗斷言。

## 10. 稽核範圍限制

- 未在 R 或 Nextflow 中實際執行。D-03、D-04、D-05 的 R 端後果，以及 D-08、R-06，都是靜態確認。
- 依「衝突避免」要求，未深入審查 `figures.R` 與 figure policy，也沒有 GUI 原始碼可供審查。
- 富集分析的 R 模組（GO、KEGG、GSEA）只檢查了它們和契約的介面，沒有逐行審查統計邏輯。
