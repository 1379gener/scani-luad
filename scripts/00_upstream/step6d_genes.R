suppressPackageStartupMessages({
  library(SummarizedExperiment); library(survival)
  library(survminer); library(ggplot2); library(dplyr)
})

DATA_DIR   <- path.expand("~/caf_spatial_project/data/TCGA_LUAD")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/05_tcga_survival")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/05_tcga_survival")

se       <- readRDS(file.path(DATA_DIR, "TCGA_LUAD_SE.rds"))
expr_mat <- log2(assay(se, "tpm_unstrand") + 1)
rownames(expr_mat) <- as.data.frame(rowData(se))$gene_name
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

# 构建 patient_id -> 列名 映射（取第一个样本）
sample_patient <- substr(colnames(expr_mat), 1, 12)
pid_to_col <- tapply(colnames(expr_mat), sample_patient, `[`, 1)

df_surv <- read.csv(file.path(RESULT_DIR, "TCGA_LUAD_niche_survival_v2.csv"),
                    stringsAsFactors=FALSE)
df_surv$patient_id <- substr(df_surv$sample_id, 1, 12)
df_surv$stage_simple <- case_when(
  df_surv$ajcc_pathologic_stage %in% c("Stage I","Stage IA","Stage IB") ~ "I",
  df_surv$ajcc_pathologic_stage %in% c("Stage IIA","Stage IIB")         ~ "II",
  df_surv$ajcc_pathologic_stage %in% c("Stage IIIA","Stage IIIB")       ~ "III",
  df_surv$ajcc_pathologic_stage == "Stage IV"                            ~ "IV",
  TRUE ~ NA_character_
)

df2 <- df_surv[!is.na(df_surv$stage_simple) & df_surv$stage_simple != "I" &
               !is.na(df_surv$OS_time) & df_surv$OS_time > 0, ]
cat(sprintf("Stage II-IV: n=%d, 死亡=%d\n", nrow(df2), sum(df2$OS_event)))

key_genes <- c("ACTA2","COL1A1","COL1A2","MMP2","TAGLN","FAP",
               "POSTN","PDPN","S100A4","THBS2","FN1","VIM",
               "SPARC","TIMP1","CXCL12","TGFb1","TGFB1")
key_genes <- intersect(key_genes, rownames(expr_mat))

cat("\n单基因预后分析（Stage II-IV）:\n")
gene_results <- data.frame()
for (g in key_genes) {
  col_ids  <- pid_to_col[df2$patient_id]
  valid    <- !is.na(col_ids)
  expr_g   <- rep(NA_real_, nrow(df2))
  expr_g[valid] <- as.numeric(expr_mat[g, col_ids[valid]])
  df2$gene_expr <- expr_g

  df2$grp <- ifelse(df2$gene_expr >= median(df2$gene_expr, na.rm=TRUE),
                    "High","Low")
  lr  <- tryCatch(survdiff(Surv(OS_time, OS_event) ~ grp, data=df2),
                  error=function(e) NULL)
  if (is.null(lr)) next
  p   <- 1 - pchisq(lr$chisq, df=1)
  cox <- coxph(Surv(OS_time, OS_event) ~ I(scale(gene_expr)), data=df2)
  hr  <- exp(coef(cox))
  hr_p <- summary(cox)$coefficients[,"Pr(>|z|)"]
  cat(sprintf("  %-10s KM_p=%.4f  HR=%.3f  p=%.4f\n", g, p, hr, hr_p))
  gene_results <- rbind(gene_results,
    data.frame(gene=g, KM_p=p, HR=round(hr,3), Cox_p=round(hr_p,4)))
}

gene_results <- gene_results[order(gene_results$KM_p), ]
write.csv(gene_results,
          file.path(RESULT_DIR, "single_gene_survival_StageII_IV.csv"),
          row.names=FALSE)

# 最显著基因做KM图
best_gene <- gene_results$gene[1]
cat(sprintf("\n最佳单基因: %s (KM_p=%.4f)\n",
            best_gene, gene_results$KM_p[1]))

col_ids       <- pid_to_col[df2$patient_id]
valid         <- !is.na(col_ids)
df2$best_expr <- NA_real_
df2$best_expr[valid] <- as.numeric(expr_mat[best_gene, col_ids[valid]])
df2$grp <- ifelse(df2$best_expr >= median(df2$best_expr, na.rm=TRUE),
                  "High","Low")

fit <- survfit(Surv(OS_time/30, OS_event) ~ grp, data=df2)
pdf(file.path(FIG_DIR, paste0("KM_", best_gene, "_StageII_IV.pdf")),
    width=7, height=7)
print(ggsurvplot(fit, data=df2,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.25,
  palette=c("#d73027","#4575b4"), legend.labs=c("High","Low"),
  xlab="Time (months)", ylab="Overall Survival",
  title=paste0(best_gene, " Expression vs OS (TCGA-LUAD Stage II-IV)"),
  ggtheme=theme_classic()))
dev.off()

cat("✅ 完成\n")
