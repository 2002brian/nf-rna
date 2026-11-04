# nf-rna v1.7.0 雲端稽核修補紀錄（2026-11-04）

- 原稽核報告：`docs/dev/audit/2026-11-04-cloud-audit.zh-TW.md`（基準 `c75e962`，v1.4.0 之後，原樣保留作為歷史證據）
- 本次基準：`feature/v1.7.0` @ `e1edcc6`
- 修補分支：`fix/v1.7.0-audit-remediation`
- 原則：先在 v1.7.0 上重現，再修補；只做聚焦修改，不改統計方法、Salmon/HISAT2 科學預設、normalization、VST、PCA、相關係數或 DESeq2 計算

## 1. 重新評估結果

原稽核分支的 14 個 strict-xfail 測試，移植到 v1.7.0 後以 `--runxfail` 執行，**14 個全部重現**。

| ID | v1.7.0 狀態 | 處理 |
|---|---|---|
| D-01 L1→L2 reuse 因 `sha256_source` 被拒 | 仍存在；以真實 L1 → `rnaseq design` → L2 流程重現 | 已修 |
| D-02 `contrast_id` 路徑穿越 | 仍存在 | 已修 |
| D-03 設計變數字面 `NA` | 仍存在 | 已修 |
| D-04 `sample_id = NA` | 仍存在 | 已修 |
| D-05 `1_000`、非 ASCII 數字 | 仍存在 | 已修 |
| D-06 sample ID 未受限、插入 shell | 仍存在 | 已修（Python 驗證 + workflow 防線） |
| D-07 交付包 run_state 為 RUNNING | 仍存在 | 已修 |
| D-08 L1 圖分組變數 | 仍存在（v1.7.0 圖形修改未處理 `color_name`） | 已修（僅靜態測試，雲端無 R） |
| R-01 非原子寫入 | 仍存在 | 已修 |
| R-03 合併檔與 lane 檔重複計數 | 仍存在 | 已修 |
| GUI 在行程內執行 run（稽核 §5-2） | 已不存在：`gui/jobs.py` 以 `subprocess.Popen(start_new_session=True)` 啟動 | 無需處理 |
| R-05 HISAT2 驗收證據過期 | 仍存在，且本次 workflow 再變更 | 需本地 E2E |
| R-06 slashy regex `/\\.fastq$/` | 仍存在（外觀問題） | 未修，需 Nextflow 實跑 |
| R-07 `.replace('.counts','')` | 仍存在（會大聲失敗） | 未修 |
| R-02、R-04、R-08、R-09、R-10、T-01 | 未變 | 未修，列為剩餘風險 |

## 2. 測試

雲端容器：Python 3.13，無 R、Nextflow、Conda、HISAT2、featureCounts。

| 狀態 | 結果 |
|---|---|
| 修補前（`e1edcc6`） | 995 passed, 8 failed, 194 skipped |
| 修補後 | 1135 passed, 8 failed, 194 skipped |
| `tests/audit/` | 140 passed |

8 個失敗在修補前後的錯誤輸出逐字相同，全部為環境因素：7 個需要 `conda`／`nextflow` 執行檔，1 個是 Python 3.13 的 symlink 錯誤訊息差異（稽核 T-01）。194 個 skip 主要是 Rscript 不存在（167 個）。R 端與 Nextflow 端的行為必須在本地 Ubuntu Server 驗證。
