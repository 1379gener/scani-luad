suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GSVA)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(pROC)
})

DATA_DIR   <- path.expand("~/caf_spatial_project/data/TCGA_LUAD")
SC_RES     <- path.expand("~/caf_spatial_project/results/01_scRNA")
ST_RES     <- path.expand("~/caf_spatial_project/results/03_spatial_GSE307534")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/05_tcga_survival")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/05_tcga_survival")
dir.create(RESULT_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG_DIR,    showWarnings=FALSE, recursive=TRUE)

# ── 1. 读取数据 ──────────────────────────────────────────
cat("[1] 读取TCGA-LUAD数据...\n")
se       <- readRDS(file.path(DATA_DIR, "TCGA_LUAD_SE.rds"))
CAF_SIGS <- readRDS(file.path(SC_RES,  "CAF_subtype_signatures.rds"))

expr_mat <- log2(assay(se, "tpm_unstrand") + 1)
gene_info <- as.data.frame(rowData(se))
rownames(expr_mat) <- gene_info$gene_name
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
cat("  表达矩阵:", nrow(expr_mat), "基因 x", ncol(expr_mat), "样本\n")

# ── 2. 临床数据整理 ──────────────────────────────────────
cat("[2] 整合临床数据...\n")
clin <- as.data.frame(colData(se))
clin$patient_id <- substr(rownames(clin), 1, 12)
clin$OS_time    <- ifelse(clin$vital_status == "Dead",
                           clin$days_to_death,
                           clin$days_to_last_follow_up)
clin$OS_event   <- ifelse(clin$vital_status == "Dead", 1, 0)
clin$stage_simple <- case_when(
  clin$ajcc_pathologic_stage %in%
    c("Stage I","Stage IA","Stage IB")    ~ "I",
  clin$ajcc_pathologic_stage %in%
    c("Stage IIA","Stage IIB")            ~ "II",
  clin$ajcc_pathologic_stage %in%
    c("Stage IIIA","Stage IIIB")          ~ "III",
  clin$ajcc_pathologic_stage == "Stage IV" ~ "IV",
  TRUE ~ NA_character_
)
clin <- clin[!is.na(clin$OS_time) & clin$OS_time > 0, ]
cat("  有效样本:", nrow(clin),
    " 死亡:", sum(clin$OS_event, na.rm=TRUE), "\n")
cat("  分期分布:\n")
print(table(clin$stage_simple, useNA="ifany"))

# ── 3. ssGSEA 计算 CAF 亚型评分 ─────────────────────────
cat("[3] 计算ssGSEA评分 (myCAF/iCAF/apCAF)...\n")
sig_filt <- lapply(CAF_SIGS, function(g) intersect(g, rownames(expr_mat)))
cat("  Signature基因覆盖:\n")
for (nm in names(sig_filt))
  cat("   ", nm, ":", length(sig_filt[[nm]]), "/", length(CAF_SIGS[[nm]]), "\n")

scores_raw <- gsva(expr_mat, sig_filt,
                   method="ssgsea", ssgsea.norm=FALSE, verbose=FALSE)
scores_df            <- as.data.frame(t(scores_raw))
scores_df$sample_id  <- rownames(scores_df)
scores_df$patient_id <- substr(rownames(scores_df), 1, 12)

# ── 4. 提取 TIMP1 / PDPN 表达 ───────────────────────────
cat("[4] 提取TIMP1/PDPN表达...\n")
pid_to_col <- tapply(colnames(expr_mat),
                     substr(colnames(expr_mat),1,12), `[`, 1)

scores_df$TIMP1 <- as.numeric(
  expr_mat["TIMP1", pid_to_col[scores_df$patient_id]])
scores_df$PDPN  <- as.numeric(
  expr_mat["PDPN",  pid_to_col[scores_df$patient_id]])

# ── 5. 合并临床数据 ──────────────────────────────────────
cat("[5] 合并临床数据...\n")
merged <- merge(
  scores_df,
  clin[, c("patient_id","OS_time","OS_event",
           "stage_simple","age_at_index")],
  by = "patient_id"
)
merged <- merged[!duplicated(merged$patient_id), ]
cat("  合并后总样本:", nrow(merged),
    " 死亡:", sum(merged$OS_event), "\n")

# Stage II-IV 亚组
df2 <- merged[!is.na(merged$stage_simple) &
              merged$stage_simple != "I", ]
cat("  Stage II-IV: n=", nrow(df2),
    " 死亡:", sum(df2$OS_event), "\n")

# ── 6. 构建 CAF Index ────────────────────────────────────
cat("[6] 构建CAF Index = Z(myCAF) + Z(TIMP1) - Z(PDPN)...\n")
df2$myCAF_Z   <- scale(df2$myCAF)[,1]
df2$TIMP1_Z   <- scale(df2$TIMP1)[,1]
df2$PDPN_Z    <- scale(df2$PDPN)[,1]
df2$CAF_Index <- df2$myCAF_Z + df2$TIMP1_Z - df2$PDPN_Z
df2$CAF_group <- ifelse(df2$CAF_Index >= median(df2$CAF_Index),
                         "CAF_High", "CAF_Low")

cat("  CAF Index范围:",
    round(min(df2$CAF_Index),3), "~",
    round(max(df2$CAF_Index),3), "\n")
cat("  CAF_High n=", sum(df2$CAF_group=="CAF_High"),
    " CAF_Low n=", sum(df2$CAF_group=="CAF_Low"), "\n")

# ── 7. 各评分汇总 ────────────────────────────────────────
cat("[7] 各评分预后价值汇总 (Stage II-IV)...\n")
cat(sprintf("%-20s %8s %8s %10s %8s\n",
            "Score","KM_p","HR","95%CI","Cox_p"))

summary_tbl <- data.frame()
for (sc in c("myCAF","iCAF","apCAF","TIMP1","CAF_Index")) {
  vals     <- df2[[sc]]
  df2$grp  <- ifelse(vals >= median(vals,na.rm=TRUE),"H","L")
  lr       <- survdiff(Surv(OS_time,OS_event) ~ grp, data=df2)
  km_p     <- 1 - pchisq(lr$chisq, df=1)
  cox      <- coxph(Surv(OS_time,OS_event) ~ I(scale(vals)), data=df2)
  res      <- summary(cox)$coefficients
  hr       <- exp(res[,"coef"])
  lo       <- exp(res[,"coef"] - 1.96*res[,"se(coef)"])
  hi       <- exp(res[,"coef"] + 1.96*res[,"se(coef)"])
  cox_p    <- res[,"Pr(>|z|)"]
  cat(sprintf("%-20s %8.4f %8.3f %5.3f-%-5.3f %8.4f\n",
              sc, km_p, hr, lo, hi, cox_p))
  summary_tbl <- rbind(summary_tbl,
    data.frame(Score=sc, KM_p=round(km_p,4),
               HR=round(hr,3), CI_lo=round(lo,3),
               CI_hi=round(hi,3), Cox_p=round(cox_p,4)))
}
write.csv(summary_tbl,
          file.path(RESULT_DIR,"score_summary_StageII_IV.csv"),
          row.names=FALSE)

# ── 8. CAF Index KM 图 ──────────────────────────────────
cat("[8] 生成KM图...\n")
surv_obj   <- Surv(df2$OS_time/30, df2$OS_event)
fit_idx    <- survfit(surv_obj ~ CAF_group, data=df2)
lr_idx     <- survdiff(surv_obj ~ CAF_group, data=df2)
p_idx      <- 1 - pchisq(lr_idx$chisq, df=1)
cat("  ★ CAF Index KM Log-rank p =", round(p_idx,4), "\n")

p_km <- ggsurvplot(fit_idx, data=df2,
  pval=TRUE, pval.method=TRUE, conf.int=TRUE,
  risk.table=TRUE, risk.table.height=0.28,
  palette=c("#C0392B","#2980B9"),
  legend.labs=c(paste0("CAF_High (n=",
                        sum(df2$CAF_group=="CAF_High"),")"),
                paste0("CAF_Low (n=",
                        sum(df2$CAF_group=="CAF_Low"),")")),
  legend.title="CAF Index",
  xlab="Time (months)", ylab="Overall Survival Probability",
  title=sprintf("CAF Index vs OS\nTCGA-LUAD Stage II-IV (n=%d, p=%.4f)",
                nrow(df2), p_idx),
  ggtheme=theme_bw(base_size=12))
pdf(file.path(FIG_DIR,"KM_CAFIndex_final.pdf"), width=8, height=7)
print(p_km); dev.off()

# myCAF KM
df2$grp_myCAF <- ifelse(df2$myCAF >= median(df2$myCAF),
                         "myCAF_High","myCAF_Low")
fit_my <- survfit(surv_obj ~ grp_myCAF, data=df2)
lr_my  <- survdiff(surv_obj ~ grp_myCAF, data=df2)
p_my   <- 1 - pchisq(lr_my$chisq, df=1)

p_km2 <- ggsurvplot(fit_my, data=df2,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE,
  risk.table.height=0.28,
  palette=c("#E74C3C","#3498DB"),
  legend.labs=c("myCAF_High","myCAF_Low"),
  xlab="Time (months)", ylab="Overall Survival Probability",
  title=sprintf("myCAF Score vs OS\nTCGA-LUAD Stage II-IV (p=%.4f)", p_my),
  ggtheme=theme_bw(base_size=12))
pdf(file.path(FIG_DIR,"KM_myCAF_final.pdf"), width=8, height=7)
print(p_km2); dev.off()

# ── 9. 多变量 Cox 回归 ───────────────────────────────────
cat("[9] 多变量Cox回归 (CAF_Index + Stage)...\n")
df_cox      <- df2[!is.na(df2$stage_simple) &
                   !is.na(df2$age_at_index), ]
df_cox$stage_f <- factor(df_cox$stage_simple, levels=c("II","III","IV"))
cox_mv      <- coxph(Surv(OS_time,OS_event) ~ CAF_Index + stage_f,
                     data=df_cox)
res_mv      <- summary(cox_mv)$coefficients
cat("  多变量Cox结果:\n")
for (i in 1:nrow(res_mv)) {
  cat(sprintf("  %-18s HR=%.3f (%.3f-%.3f) p=%.4f\n",
    rownames(res_mv)[i],
    exp(res_mv[i,"coef"]),
    exp(res_mv[i,"coef"]-1.96*res_mv[i,"se(coef)"]),
    exp(res_mv[i,"coef"]+1.96*res_mv[i,"se(coef)"]),
    res_mv[i,"Pr(>|z|)"]))
}

pdf(file.path(FIG_DIR,"Forest_CAFIndex_multivariate.pdf"),
    width=8, height=5)
ggforest(cox_mv, data=df_cox,
         main="Multivariate Cox: CAF Index + Stage\nTCGA-LUAD Stage II-IV")
dev.off()

# ── 10. ROC / AUC (3年OS) ───────────────────────────────
cat("[10] ROC/AUC分析 (3年OS)...\n")
df2$OS_3yr <- ifelse(df2$OS_time <= 1095 & df2$OS_event==1, 1, 0)
df_roc     <- df2[df2$OS_time >= 1095 | df2$OS_event==1, ]
cat("  有效样本(3年AUC):", nrow(df_roc),
    " 3年内死亡:", sum(df_roc$OS_3yr), "\n")

auc_tbl <- data.frame()
for (sc in c("myCAF","TIMP1","CAF_Index")) {
  roc_obj <- roc(df_roc$OS_3yr, df_roc[[sc]],
                 quiet=TRUE, direction="auto")
  auc_val <- round(as.numeric(auc(roc_obj)), 3)
  cat(sprintf("  %-15s AUC = %.3f\n", sc, auc_val))
  auc_tbl <- rbind(auc_tbl,
                   data.frame(Score=sc, AUC=auc_val))
}
write.csv(auc_tbl, file.path(RESULT_DIR,"AUC_3year_OS.csv"),
          row.names=FALSE)

# ROC曲线图
pdf(file.path(FIG_DIR,"ROC_3year_OS.pdf"), width=6, height=6)
cols <- c("#e74c3c","#2980b9","#27ae60")
plot(NULL, xlim=c(1,0), ylim=c(0,1),
     xlab="1 - Specificity", ylab="Sensitivity",
     main="ROC Curves: 3-year OS\nTCGA-LUAD Stage II-IV")
abline(a=1, b=-1, col="grey70", lty=2)
for (i in seq_along(c("myCAF","TIMP1","CAF_Index"))) {
  sc      <- c("myCAF","TIMP1","CAF_Index")[i]
  roc_obj <- roc(df_roc$OS_3yr, df_roc[[sc]],
                 quiet=TRUE, direction="auto")
  lines(roc_obj, col=cols[i], lwd=2)
  text(0.25, 0.22-(i-1)*0.08,
       sprintf("%s  AUC=%.3f", sc,
               round(as.numeric(auc(roc_obj)),3)),
       col=cols[i], adj=0, cex=0.9)
}
dev.off()

# ── 11. CAF Index vs Stage 相关性 ───────────────────────
cat("[11] CAF Index与分期相关性...\n")
kw <- kruskal.test(CAF_Index ~ stage_simple,
                   data=merged[!is.na(merged$stage_simple),])
cat("  Kruskal-Wallis p =", round(kw$p.value,4), "\n")
for (s in c("I","II","III","IV")) {
  sub <- merged[!is.na(merged$stage_simple) &
                merged$stage_simple==s, ]
  cat(sprintf("  Stage %s: n=%d  median CAF_Index=%.3f\n",
              s, nrow(sub),
              median(sub$CAF_Index, na.rm=TRUE)))
}

# ── 12. 保存最终结果 ─────────────────────────────────────
cat("[12] 保存最终结果...\n")
write.csv(
  df2[, c("patient_id","myCAF","iCAF","apCAF",
          "TIMP1","PDPN","CAF_Index","CAF_group",
          "OS_time","OS_event","stage_simple")],
  file.path(RESULT_DIR,"CAF_final_scores.csv"),
  row.names=FALSE
)

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║     TCGA-LUAD CAF Index 分析结果汇总          ║\n")
cat("╠══════════════════════════════════════════════╣\n")
cat(sprintf("║  队列: Stage II-IV  n=%-4d  死亡=%-3d         ║\n",
            nrow(df2), sum(df2$OS_event)))
cat(sprintf("║  KM Log-rank p = %-8.4f                  ║\n", p_idx))
cat(sprintf("║  单变量 HR = %-6.3f  (p<0.0001)           ║\n",
            exp(coxph(Surv(OS_time,OS_event)~CAF_Index,
                      data=df2)$coefficients)))
cat(sprintf("║  多变量 HR = %-6.3f  p=%-7.4f            ║\n",
            exp(res_mv["CAF_Index","coef"]),
            res_mv["CAF_Index","Pr(>|z|)"]))
cat(sprintf("║  3年OS AUC = %-5.3f                       ║\n",
            auc_tbl$AUC[auc_tbl$Score=="CAF_Index"]))
cat(sprintf("║  分期Kruskal p = %-6.4f (分期独立性)      ║\n",
            kw$p.value))
cat("╚══════════════════════════════════════════════╝\n")
cat("\n✅ Step6完成！图像保存至:", FIG_DIR, "\n")
