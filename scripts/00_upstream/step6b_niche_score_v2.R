suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GSVA)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- path.expand("~/caf_spatial_project/data/TCGA_LUAD")
ST_RES     <- path.expand("~/caf_spatial_project/results/03_spatial_GSE307534")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/05_tcga_survival")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/05_tcga_survival")

cat("\n[1] 读取数据...\n")
se         <- readRDS(file.path(DATA_DIR, "TCGA_LUAD_SE.rds"))
core_genes <- readLines(file.path(ST_RES, "CAF_Niche_final_geneset.txt"))
clin       <- read.csv(file.path(DATA_DIR, "TCGA_LUAD_clinical.csv"),
                       stringsAsFactors=FALSE)

cat("\n[2] 提取表达矩阵...\n")
expr_mat <- log2(assay(se, "tpm_unstrand") + 1)
gene_info <- as.data.frame(rowData(se))
rownames(expr_mat) <- gene_info$gene_name
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

valid_genes <- intersect(core_genes, rownames(expr_mat))
cat(sprintf("  有效基因: %d\n", length(valid_genes)))

cat("\n[3] 计算三种评分方法...\n")

# 方法A: 基因集平均表达（最简单可靠）
score_mean <- colMeans(expr_mat[valid_genes, ], na.rm=TRUE)

# 方法B: ssGSEA 不normalize
gsva_raw <- gsva(expr_mat,
                 list(CAF_Niche=valid_genes),
                 method="ssgsea",
                 ssgsea.norm=FALSE,
                 verbose=FALSE)
score_ssgsea <- as.numeric(gsva_raw["CAF_Niche", ])
names(score_ssgsea) <- colnames(gsva_raw)

# 方法C: GSVA（不同算法）
gsva_score <- gsva(expr_mat,
                   list(CAF_Niche=valid_genes),
                   method="gsva",
                   verbose=FALSE)
score_gsva <- as.numeric(gsva_score["CAF_Niche", ])
names(score_gsva) <- colnames(gsva_score)

cat(sprintf("  mean_expr 范围: %.3f ~ %.3f\n",
            min(score_mean), max(score_mean)))
cat(sprintf("  ssGSEA(raw) 范围: %.3f ~ %.3f\n",
            min(score_ssgsea), max(score_ssgsea)))
cat(sprintf("  GSVA 范围: %.3f ~ %.3f\n",
            min(score_gsva), max(score_gsva)))

cat("\n[4] 整合临床数据...\n")
surv_df <- data.frame(
  sample_id    = colnames(expr_mat),
  score_mean   = score_mean,
  score_ssgsea = score_ssgsea,
  score_gsva   = score_gsva,
  patient_id   = substr(colnames(expr_mat), 1, 12),
  stringsAsFactors=FALSE
)

clin$patient_id <- substr(clin$patient, 1, 12)
surv_df <- merge(surv_df,
  clin[, c("patient_id","vital_status","days_to_death",
           "days_to_last_follow_up","ajcc_pathologic_stage",
           "age_at_index")],
  by="patient_id", all.x=TRUE)

surv_df$OS_time  <- ifelse(surv_df$vital_status=="Dead",
                            surv_df$days_to_death,
                            surv_df$days_to_last_follow_up)
surv_df$OS_event <- ifelse(surv_df$vital_status=="Dead", 1, 0)
surv_df <- surv_df[!is.na(surv_df$OS_time) & surv_df$OS_time > 0, ]

cat("\n[5] 比较三种评分方法的生存分析结果...\n")
for (score_col in c("score_mean","score_ssgsea","score_gsva")) {
  df_tmp <- surv_df[!is.na(surv_df[[score_col]]), ]
  df_tmp$grp <- ifelse(df_tmp[[score_col]] >= median(df_tmp[[score_col]]),
                       "High", "Low")
  lr <- survdiff(Surv(OS_time, OS_event) ~ grp, data=df_tmp)
  p  <- 1 - pchisq(lr$chisq, df=1)
  cat(sprintf("  %-15s: Log-rank p = %.4f\n", score_col, p))
}

# 用最好的方法继续分析
cat("\n[6] 用 score_mean 做完整分析...\n")
surv_df$niche_score <- surv_df$score_mean
surv_df$niche_group <- ifelse(
  surv_df$niche_score >= median(surv_df$niche_score),
  "Niche_High", "Niche_Low")

# 全队列KM
fit <- survfit(Surv(OS_time/30, OS_event) ~ niche_group, data=surv_df)
lr  <- survdiff(Surv(OS_time/30, OS_event) ~ niche_group, data=surv_df)
p   <- 1 - pchisq(lr$chisq, df=1)
cat(sprintf("  全队列 OS p = %.4f\n", p))

pdf(file.path(FIG_DIR, "KM_OS_NicheScore_v2.pdf"), width=7, height=7)
print(ggsurvplot(fit, data=surv_df,
  pval=TRUE, pval.method=TRUE, conf.int=TRUE,
  risk.table=TRUE, risk.table.height=0.25,
  palette=c("#d73027","#4575b4"),
  legend.labs=c("Niche_High","Niche_Low"),
  legend.title="CAF Niche Score",
  xlab="Time (months)", ylab="Overall Survival",
  title="TCGA-LUAD: CAF Niche Score vs OS",
  ggtheme=theme_classic()))
dev.off()

# Stage II-IV 亚组
surv_late <- surv_df[grepl("Stage II|Stage III|Stage IV",
                            surv_df$ajcc_pathologic_stage), ]
cat(sprintf("  晚期样本(II-IV): %d, 死亡: %d\n",
            nrow(surv_late), sum(surv_late$OS_event)))

surv_late$niche_group <- ifelse(
  surv_late$niche_score >= median(surv_late$niche_score),
  "Niche_High", "Niche_Low")

fit_late <- survfit(Surv(OS_time/30, OS_event) ~ niche_group,
                    data=surv_late)
lr_late  <- survdiff(Surv(OS_time/30, OS_event) ~ niche_group,
                     data=surv_late)
p_late   <- 1 - pchisq(lr_late$chisq, df=1)
cat(sprintf("  晚期亚组 OS p = %.4f\n", p_late))

pdf(file.path(FIG_DIR, "KM_OS_NicheScore_StageII_IV.pdf"), width=7, height=7)
print(ggsurvplot(fit_late, data=surv_late,
  pval=TRUE, pval.method=TRUE, conf.int=TRUE,
  risk.table=TRUE, risk.table.height=0.25,
  palette=c("#d73027","#4575b4"),
  legend.labs=c("Niche_High","Niche_Low"),
  legend.title="CAF Niche Score",
  xlab="Time (months)", ylab="Overall Survival",
  title="TCGA-LUAD Stage II-IV: CAF Niche Score vs OS",
  ggtheme=theme_classic()))
dev.off()

write.csv(surv_df,
          file.path(RESULT_DIR, "TCGA_LUAD_niche_survival_v2.csv"),
          row.names=FALSE)

cat("\n✅ 完成！\n")
