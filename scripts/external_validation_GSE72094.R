#!/usr/bin/env Rscript
# =====================================================================
# external_validation_GSE72094.R
# 外部验证：在独立 LUAD bulk 队列 GSE72094 上重复 TCGA 的核心结果
#   - ECM communication score 是否仍是不良预后因子（Cox, 校正分期）
#   - ECM score 与免疫抑制标志物的相关性
# 目的：为合并版加一个独立验证队列，提升高分期刊命中率。
#
# 运行环境：conda activate seg_env（需 GEOquery / GSVA / survival 等）
# 运行：    Rscript external_validation_GSE72094.R
# 说明：GSE72094 ≈ 442 例 LUAD，平台 GPL15048，含生存与分期信息。
# =====================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(GSVA)
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
})

set.seed(1)
proj    <- "~/caf_spatial_project"
out_dir <- file.path(proj, "results/merged_external_GSE72094")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 与论文完全一致的两套基因集
ECM_GENES    <- c("COL1A1","COL1A2","COL4A1","COL6A1","COL6A2","COL6A3","FN1")
CXCLMK_GENES <- c("CXCL12","MDK","CXCL2")
IMMUNE_MARKERS <- c("CD8A","FOXP3","PDCD1","CD274")

# ---------------------------------------------------------------------
# 1. 下载 GSE72094
# ---------------------------------------------------------------------
cat("== 下载 GSE72094（首次较慢）==\n")
gse <- getGEO("GSE72094", GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]
expr <- exprs(eset)
pdat <- pData(eset)
fdat <- fData(eset)

cat("[诊断] 表达矩阵:", paste(dim(expr), collapse=" x "), "\n")
cat("[诊断] pData 列名（找生存/分期字段）:\n"); print(colnames(pdat))
cat("[诊断] fData 列名（找基因符号字段）:\n"); print(colnames(fdat))

# ---------------------------------------------------------------------
# 2. 探针 -> 基因符号（GPL15048 注释里通常有 GeneSymbol 列；按实际列名改）
# ---------------------------------------------------------------------
sym_col <- colnames(fdat)[grepl("symbol|gene", colnames(fdat), ignore.case = TRUE)][1]
cat("[诊断] 使用基因符号列:", sym_col, "\n")
sym <- as.character(fdat[[sym_col]])
keep <- !is.na(sym) & sym != "" & sym != "---"
expr <- expr[keep, ]; sym <- sym[keep]

# 同名探针取均值最高的一个
ord  <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[ord, ]; sym <- sym[ord]
expr <- expr[!duplicated(sym), ]
rownames(expr) <- sym[!duplicated(sym)]
cat("[诊断] 去重后基因数:", nrow(expr), "\n")

# 缺失基因提示
miss <- setdiff(c(ECM_GENES, CXCLMK_GENES), rownames(expr))
if (length(miss)) cat("!! 缺失基因（将被跳过）:", paste(miss, collapse=", "), "\n")

# ---------------------------------------------------------------------
# 3. ssGSEA 打分（与论文一致）
# ---------------------------------------------------------------------
gs <- list(ECM = intersect(ECM_GENES, rownames(expr)),
           CXCLMK = intersect(CXCLMK_GENES, rownames(expr)))
par <- ssgseaParam(as.matrix(expr), gs)         # GSVA >= 1.50 新接口
scores <- t(gsva(par))
scores <- as.data.frame(scores)
scores$sample <- rownames(scores)
cat("[诊断] 打分完成，样本数:", nrow(scores), "\n")

# ---------------------------------------------------------------------
# 4. 整理生存与分期（GSE72094 字段名因版本而异，下面做模糊匹配）
# ---------------------------------------------------------------------
getcol <- function(pat) {
  hit <- colnames(pdat)[grepl(pat, colnames(pdat), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA
}
col_status <- getcol("vital|status|death")
col_time   <- getcol("survival.*day|time.*day|os.*time|follow")
col_stage  <- getcol("stage")
cat("[诊断] 生存状态列:", col_status, "| 时间列:", col_time, "| 分期列:", col_stage, "\n")

surv <- data.frame(
  sample = rownames(pdat),
  status_raw = if (!is.na(col_status)) pdat[[col_status]] else NA,
  time_raw   = if (!is.na(col_time))   pdat[[col_time]]   else NA,
  stage_raw  = if (!is.na(col_stage))  pdat[[col_stage]]  else NA,
  stringsAsFactors = FALSE
)
# 状态编码：Dead/dead/1 -> 1，Alive -> 0
surv$status <- ifelse(grepl("dead|decease|^1$", surv$status_raw, ignore.case = TRUE), 1, 0)
surv$time   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", surv$time_raw)))
# 分期简化为 I/II/III/IV -> 1..4
st <- toupper(gsub("[^IV]", "", toupper(surv$stage_raw)))
surv$stage_num <- match(st, c("I","II","III","IV"))

df <- merge(scores, surv, by = "sample")
df <- df[!is.na(df$time) & df$time > 0, ]
cat("[诊断] 可用于生存分析样本:", nrow(df), "\n")

# ---------------------------------------------------------------------
# 5. Cox（校正分期）+ KM
# ---------------------------------------------------------------------
df$ECM_z    <- scale(df$ECM)[,1]
df$CXCLMK_z <- scale(df$CXCLMK)[,1]

cat("\n================ Cox（ECM + CXCL/MK + stage）================\n")
cox_fit <- coxph(Surv(time, status) ~ ECM_z + CXCLMK_z + stage_num, data = df)
print(summary(cox_fit))

cat("\n================ 仅 ECM（校正分期）================\n")
print(summary(coxph(Surv(time, status) ~ ECM_z + stage_num, data = df)))

# KM：按 ECM 中位数分组
df$ECM_grp <- ifelse(df$ECM > median(df$ECM), "ECM-high", "ECM-low")
km <- survfit(Surv(time, status) ~ ECM_grp, data = df)
p_km <- ggsurvplot(km, data = df, pval = TRUE, risk.table = TRUE,
                   palette = c("#d73027","#4575b4"),
                   title = "GSE72094: OS by ECM communication score")
pdf(file.path(out_dir, "GSE72094_KM_ECM.pdf"), width = 7, height = 6)
print(p_km); dev.off()

# ---------------------------------------------------------------------
# 6. ECM score 与免疫标志物相关性
# ---------------------------------------------------------------------
cat("\n================ ECM score vs 免疫标志物（Spearman）================\n")
imm_res <- lapply(IMMUNE_MARKERS, function(g) {
  if (!g %in% rownames(expr)) return(data.frame(gene=g, rho=NA, p=NA))
  v <- expr[g, df$sample]
  ct <- suppressWarnings(cor.test(df$ECM, as.numeric(v), method = "spearman"))
  data.frame(gene = g, rho = round(unname(ct$estimate),3), p = signif(ct$p.value,3))
})
imm_res <- do.call(rbind, imm_res)
print(imm_res, row.names = FALSE)

# ---------------------------------------------------------------------
# 7. 保存
# ---------------------------------------------------------------------
write.csv(df,      file.path(out_dir, "GSE72094_scores_survival.csv"), row.names = FALSE)
write.csv(imm_res, file.path(out_dir, "GSE72094_ECM_immune_corr.csv"), row.names = FALSE)
saveRDS(cox_fit,   file.path(out_dir, "GSE72094_cox_fit.rds"))

cat("\n完成。把 Cox 的 ECM HR/95%CI/p 和 KM 的 p 发我，\n")
cat("我加进合并版作为独立外部验证（新增一段 Results + 一个补充图）。\n")
