# v1.8.0：PPI 開發準備

狀態：實驗性前處理，尚未接入正式 `rnaseq run`、Nextflow、GUI 或交付報告。
目前套件版本仍為 1.7.0；分支名稱不代表已發布 1.8.0。

第二階段已加入 `ppi_network.py`：離線 snapshot 匯入、明確 ID namespace 映射、
一對多排除與稽核、保留孤立節點的 induced network 擷取、證據分數保留及輸出
checksum。操作與格式見 [PPI 文件](../ppi.md)。目前仍需操作者提供正規化映射表，
尚未提供原始 aliases／API 自動轉換器；版本與來源真實性仍屬操作者宣告。

## 本階段已實作

`src/rnaseq/ppi.py` 使用 Python 標準函式庫，讀取同一個 L2 contrast 的
`all_genes.tsv`、`significant.tsv`、`contrast_metadata.json`，驗證 significant
集合與記錄的 DEG 門檻一致，並驗證 raw LFC、padj 與完整結果表一致。
此檢查不建立新的 DEG 定義，不使用 shrunken LFC 選基因。

產生 `all.tsv`、`up.tsv`、`down.tsv`；方向以 raw LFC 為準，正值代表 numerator
較高。若原設定允許 LFC=0 的 DEG，該基因只保留在 all 並記錄數量。
保留原始 gene_id；只有明確宣告 ENSEMBL 時，才為 query_id 移除尾端數字版本號。
query_id 碰撞完整列出，不任意挑選或合併。symbol 僅為註解。

同一批讀入的來源 bytes 會複製至輸出，並記錄 SHA-256、大小、物種 taxon ID、
contrast 方向與 DEG 定義。輸出 checksum 記於 summary.json，summary 最後寫入作為
完成標記。既有目錄一律拒絕覆寫；這個 helper 尚未提供多程序發布／崩潰復原協定。
沒有 DEG 時輸出 EMPTY；其餘成功候選集為 READY_FOR_MAPPING，mapping_status 永遠
為 NOT_RUN，不能宣稱已完成 PPI。此 helper 記錄來源內容，尚不驗證 immutable run
的來源真實性。

開發者使用方式（從 checkout，指定新的輸出目錄）：

```bash
PYTHONPATH=src python -m rnaseq.ppi \
  --contrast-dir /path/to/run/downstream/l2/contrasts/Treatment_vs_Control \
  --outdir /path/to/new/ppi_candidates \
  --organism 'Homo sapiens' --input-id-type ENSEMBL
```

## 整合路線（1、2 的離線核心已實作）

1. **STRING 資源與映射**：先準備指定版本、物種的本機 snapshot；manifest 記錄
   provider、版本、taxon ID、原始 URL、取得時間與每個檔案 checksum。
   映射表保留 gene_id、query_id、STRING protein ID 與所有一對多結果；另列
   unmapped／ambiguous、mapping rate 與每個 protein 的來源 gene_id。
   不將 symbol 的第一筆結果當作可靠映射。映射門檻在接入前另行訂定。
2. **網路擷取**：明確選擇 functional association 或 physical network，固定
   confidence 門檻並保留各 evidence channel。預設只使用輸入蛋白，不自動加入
   鄰居；保留孤立節點、輸出節點與邊表。STRING association 不能一律標為實體結合。
   全網路擷取與圖面顯示上限分開；圖面篩選不能改變完整結果表。
3. **Nextflow 整合**：新增候選整理 → ID 映射 → 網路擷取 → 報告的獨立 tasks。
   Python 負責驗證、凍結與交付，Nextflow 負責排程。第一版僅開放 L2 DEG 來源；
   WGCNA 來源另案開發。analysis 設定保持明確 opt-in；舊 schema 遷移規則需測試。
4. **來源重用與交付**：整合時須驗證 SUCCESS run、delivery manifest、物種、
   ID type、contrast 與來源 checksum；將來源 run、pipeline source revision、
   workflow/script hash、runtime lock 與 STRING snapshot 一起凍結。
   GUI、CLI 與 delivery 使用同一份結果，不各自重算。PPI 統計富集檢定另訂背景
   與多重比較契約，本階段不啟用。

## 驗收

- 離線測試涵蓋 DEG 集合一致性、raw/shrunken LFC 方向、NA padj、嚴格 padj 邊界、
  空結果、零 LFC、重複 ID、版本號碰撞、固定輸出、checksum 及拒絕覆寫。
- 映射階段必須另外涵蓋錯誤物種、缺少 snapshot、一對多、未映射與 checksum 不符。
- Nextflow 階段須驗證未選 PPI 時不建立 PPI tasks、來源不符時在執行前拒絕、
  無邊與無映射結果可明確區分；既有 DEG／ORA／GSEA 數值維持一致。

參考：[STRING API](https://string-db.org/help/api/)；正式使用須固定版本位址。
