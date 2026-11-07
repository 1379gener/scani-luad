suppressPackageStartupMessages({
  library(SummarizedExperiment); library(GSVA)
  library(survival); library(survminer)
  library(ggplot2); library(dplyr); library(patchwork)
})

DATA_DIR   <- path.expand("~/caf_spatial_project/data/TCGA_LUAD")
SC_RES     <- path.expand("~/caf_spatial_project/results/01_scRNA")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/05_tcga_survival")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/05_tcga_survival")

cat("[1] 读取数据...\n")
se        <- readRDS(file.path(DATA_DIR, "TCGA_LUAD_SE.rds"))
CAF_SIGS  <- readRDS(file.path(SC_RES, "CAF_subtype_signatures.rds"))
df_surv   <- read.csv(
  file.path(RESULT_DIR, "TCGA_LUAD_niche_survival_v2.csv"),
  stringsAsFactors=FALSE)

expr_mat <- log2(assay(se, "tpm_unstrand") + 1)
rownames(expr_mat) <- as.data.frame(rowData(se))$gene_name
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
pid_to_col <- tapply(colnames(expr_mat),
                     substr(colnames(expr_mat),1,12), `[`, 1)

cat("[2] 计算所有评分...\n")
# ssGSEA亚型评分
sig_list <- lapply(CAF_SIGS, function(g) intersect(g, rownames(expr_mat)))
gsva_res <- gsva(expr_mat, sig_list,
                 method="ssgsea", ssgsea.norm=FALSE, verbose=FALSE)
score_df <- as.data.frame(t(gsva_res))
score_df$patient_id <- substr(rownames(score_df), 1, 12)
score_df$myCAF_iCAF_ratio <- score_df$myCAF / (score_df$iCAF + 1e-6)

# 关键单基因
key_genes <- c("TIMP1","PDPN","CXCL12","TAGLN","FAP")
key_genes <- intersect(key_genes, rownames(expr_mat))
for (g in key_genes) {
  col_ids <- pid_to_col[score_df$patient_id]
  valid   <- !is.na(col_ids)
  vals    <- rep(NA_real_, nrow(score_df))
  vals[valid] <- as.numeric(expr_mat[g, col_ids[valid]])
  score_df[[g]] <- vals
}

cat("[3] 整合临床数据...\n")
df_surv$patient_id <- substr(df_surv$sample_id, 1, 12)
df_surv$stage_simple <- case_when(
  df_surv$ajcc_pathologic_stage %in%
    c("Stage I","Stage IA","Stage IB")    ~ "I",
  df_surv$ajcc_pathologic_stage %in%
    c("Stage IIA","Stage IIB")            ~ "II",
  df_surv$ajcc_pathologic_stage %in%
    c("Stage IIIA","Stage IIIB")          ~ "III",
  df_surv$ajcc_pathologic_stage == "Stage IV" ~ "IV",
  TRUE ~ NA_character_
)

df <- merge(score_df, df_surv[, c("patient_id","OS_time","OS_event",
                                   "stage_simple","age_at_index")],
            by="patient_id")
df2 <- df[!is.na(df$stage_simple) & df$stage_simple != "I" &
            !is.na(df$OS_time) & df$OS_time > 0, ]
cat(sprintf("  Stage II-IV: n=%d, 死亡=%d\n", nrow(df2), sum(df2$OS_event)))

cat("[4] 构建整合评分（CAF_Index）...\n")
# 标准化后整合：myCAF_score + TIMP1 - PDPN（根据方向）
df2$CAF_Index <- scale(df2$myCAF)[,1] +
                 scale(df2$TIMP1)[,1] -
                 scale(df2$PDPN)[,1]

df2$grp_index <- ifelse(df2$CAF_Index >= median(df2$CAF_Index),
                        "CAF_High","CAF_Low")
lr_idx <- survdiff(Surv(OS_time/30, OS_event) ~ grp_index, data=df2)
p_idx  <- 1 - pchisq(lr_idx$chisq, df=1)
cat(sprintf("  CAF_Index KM p=%.4f\n", p_idx))

cat("[5] 生成图表...\n")

# 图1: myCAF Score KM
df2$grp_myCAF <- ifelse(df2$myCAF >= median(df2$myCAF),
                        "myCAF_High","myCAF_Low")
fit1 <- survfit(Surv(OS_time/30, OS_event) ~ grp_myCAF, data=df2)
lr1  <- survdiff(Surv(OS_time/30, OS_event) ~ grp_myCAF, data=df2)
p1   <- 1 - pchisq(lr1$chisq, df=1)

# 图2: TIMP1 KM
df2$grp_timp1 <- ifelse(df2$TIMP1 >= median(df2$TIMP1, na.rm=TRUE),
                        "TIMP1_High","TIMP1_Low")
fit2 <- survfit(Surv(OS_time/30, OS_event) ~ grp_timp1, data=df2)
lr2  <- survdiff(Surv(OS_time/30, OS_event) ~ grp_timp1, data=df2)
p2   <- 1 - pchisq(lr2$chisq, df=1)

# 图3: CAF_Index KM
fit3 <- survfit(Surv(OS_time/30, OS_event) ~ grp_index, data=df2)

pdf(file.path(FIG_DIR, "KM_CAF_combined_StageII_IV.pdf"), width=18, height=7)
km1 <- ggsurvplot(fit1, data=df2,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.28,
  palette=c("#d73027","#4575b4"),
  legend.labs=c("myCAF_High","myCAF_Low"),
  xlab="Time (months)", ylab="Overall Survival",
  title=sprintf("myCAF Score (p=%.3f)\nTCGA-LUAD Stage II-IV", p1),
  ggtheme=theme_classic())
km2 <- ggsurvplot(fit2, data=df2,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.28,
  palette=c("#d73027","#4575b4"),
  legend.labs=c("TIMP1_High","TIMP1_Low"),
  xlab="Time (months)", ylab="Overall Survival",
  title=sprintf("TIMP1 Expression (p=%.3f)\nTCGA-LUAD Stage II-IV", p2),
  ggtheme=theme_classic())
km3 <- ggsurvplot(fit3, data=df2,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.28,
  palette=c("#d73027","#4575b4"),
  legend.labs=c("CAF_High","CAF_Low"),
  xlab="Time (months)", ylab="Overall Survival",
  title=sprintf("CAF Index (p=%.3f)\nTCGA-LUAD Stage II-IV", p_idx),
  ggtheme=theme_classic())
arrange_ggsurvplots(list(km1, km2, km3), ncol=3, nrow=1)
dev.off()

cat("[6] 多变量Cox回归...\n")
df_cox <- df2[!is.na(df2$age_at_index) & !is.na(df2$TIMP1), ]
df_cox$stage_f <- factor(df_cox$stage_simple, levels=c("II","III","IV"))

# 模型1: myCAF + 临床
cox1 <- coxph(Surv(OS_time,OS_event) ~
                I(scale(myCAF)) + age_at_index + stage_f,
              data=df_cox)
# 模型2: TIMP1 + 临床
cox2 <- coxph(Surv(OS_time,OS_event) ~
                I(scale(TIMP1)) + age_at_index + stage_f,
              data=df_cox)
# 模型3: CAF_Index + 临床
cox3 <- coxph(Surv(OS_time,OS_event) ~
                I(scale(CAF_Index)) + age_at_index + stage_f,
              data=df_cox)

cat("\n  模型1 (myCAF):\n")
print_cox <- function(cox_obj, label) {
  r <- summary(cox_obj)$coefficients
  cat(sprintf("  %-20s HR=%.3f CI=[%.3f-%.3f] p=%.4f\n",
    label,
    exp(r[1,"coef"]),
    exp(r[1,"coef"]-1.96*r[1,"se(coef)"]),
    exp(r[1,"coef"]+1.96*r[1,"se(coef)"]),
    r[1,"Pr(>|z|)"]))
}
print_cox(cox1, "myCAF Score")
print_cox(cox2, "TIMP1")
print_cox(cox3, "CAF_Index")

# Forest plot (最优模型)
pdf(file.path(FIG_DIR, "Cox_forest_TIMP1_StageII_IV.pdf"), width=7, height=5)
ggforest(cox2, data=df_cox,
         main="Multivariate Cox: TIMP1 + Clinical\nTCGA-LUAD Stage II-IV")
dev.off()

cat("[7] AUC比较（3年生存预测）...\n")
library(timeROC)
t_eval <- 3 * 365
for (score_col in c("myCAF","TIMP1","CAF_Index")) {
  vals <- df_cox[[score_col]]
  roc_res <- tryCatch(
    timeROC(T=df_cox$OS_time, delta=df_cox$OS_event,
            marker=vals, cause=1,
            times=t_eval, iid=FALSE),
    error=function(e) NULL)
  if (!is.null(roc_res)) {
    auc <- roc_res$AUC[2]
    cat(sprintf("  %-20s 3-year AUC=%.3f\n", score_col, auc))
  }
}

# 保存整合评分
write.csv(df2[, c("patient_id","myCAF","iCAF","apCAF",
                   "myCAF_iCAF_ratio","TIMP1","PDPN",
                   "CAF_Index","OS_time","OS_event","stage_simple")],
          file.path(RESULT_DIR, "CAF_final_scores.csv"),
          row.names=FALSE)

cat("\n✅ Step 6e 完成！\n")
cat(sprintf("  myCAF KM p=%.4f\n", p1))
cat(sprintf("  TIMP1 KM p=%.4f\n", p2))
cat(sprintf("  CAF_Index KM p=%.4f\n", p_idx))
