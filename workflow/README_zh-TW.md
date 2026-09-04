[English](README.md) | 繁體中文

# Upstream execution boundary

> 本文件為繁體中文翻譯；若與英文版內容有差異，以英文版 README.md 為準。

local Docker profile 會呼叫外部且已固定版本的 `nf-core/rnaseq 3.26.0` pipeline，接著執行此 first-party DSL2 downstream workflow：

```text
rnaseq CLI → nf-core/rnaseq → standardized upstream outputs → downstream Nextflow workflow → R modules
```

Python control plane 會凍結 version-pinned input contract，並記錄穩定的 nf-core handoff boundary。凍結的 `analysis_level` 是 graph selector：`L1` 專案只執行 L1 與 L1 technical report；`L2` 專案執行 L1、L2、選用的 GO preranked GSEA（BP/MF/CC）與 KEGG preranked GSEA，接著產生 technical HTML report。graph 絕不會從 contrast、metadata 或可用 workflow module 推斷 L2。production 中不會執行 GO 或 KEGG ORA。server executor setting 仍刻意延後決定。

預設 downstream workflow 在任一時間最多只執行一個 `ENRICHMENT_ANALYSIS` task（`maxForks 1`），以避免 concurrent R enrichment job 耗盡 local Docker host。這不會改變 module selection 或 calculation；未來可由 server-specific Nextflow configuration 覆寫此 process directive。

每個 enrichment task 都會在 immutable L2 output 下發布各自的 backend directory：`downstream/l2/enrichment/gsea_go/` 與 `downstream/l2/enrichment/gsea_kegg/`。process 不會將共用的 `enrichment/` directory 作為 task output 發布。相反地，每個 task 會相對於 L2 publish root 發布 module-specific `enrichment/*` output，因此在強制 `overwrite: false` 時，一個成功的 backend 不會遮蔽另一個 backend。

啟用 GSEA 時，final report 會接收兩個 GSEA backend 所收集的 output channel。因此，只有 GO GSEA 與 KEGG GSEA 都完成後才會啟動；任何已啟用 backend 缺失或失敗，都會阻止 report 與 delivery finalization。

## Local resource policy

checked-in local configuration 使用明確且非科學性的 runtime class：SMALL（1 CPU、2 GiB、2 h）、MEDIUM（4 CPUs、8 GiB、8 h）與 LARGE（6 CPUs、12 GiB、12 h）。L1 與 technical report 使用 SMALL；L2 與每個 GSEA backend 使用 MEDIUM。`ENRICHMENT_ANALYSIS` 保持 `maxForks 1`，凍結的 upstream local configuration 也將 nf-core `SALMON_QUANT` 限制為同一時間一個 MEDIUM task。global local resource ceiling 為 LARGE。

這些限制用於避免 local Docker memory oversubscription，不會變更任何 input、model、filtering、threshold、ranking、enrichment calculation 或已發布 artifact。未來 server profile 可以明確覆寫 resource directive。local run 前請使用 `rnaseq doctor [PROJECT]` 檢查 host/Docker capacity 與 architecture；dynamic upstream nf-core image architecture 必須在凍結的 Nextflow trace 中檢查。
