#!/usr/bin/env Rscript
# verify_scani_components.R
# 两个核实任务：
# 1. SCANI 公式去掉 PDPN 后，TCGA 多变量 Cox 表现如何？
#    （如果去掉后 HR 也类似，说明 PDPN 只是噪声补丁；如果去掉后 HR 显著下降，说明 PDPN 残差信号有用）
# 2. 找到原始 50-gene signature 来源，确认是哪个文件 + 用什么 cutoff
#
# 输出: ~/caf_spatial_project/results/suppl_tables_v90/

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

OUT_DIR <- "~/caf_spatial_project/results/suppl_tables_v90"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================
# 任务 1: SCANI 去掉 PDPN 的对比
# =============================================================
cat("\n========================================================\n")
cat("  任务 1: SCANI 公式中 PDPN 的边际贡献\n")
cat("========================================================\n")

scani <- read.csv("~/caf_spatial_project/results/05_tcga_survival/CAF_final_scores.csv",
                  stringsAsFactors = FALSE)
clin  <- read.csv("~/caf_spatial_project/results/05_tcga_survival/TCGA_LUAD_niche_survival.csv",
                  stringsAsFactors = FALSE)

clin_uniq <- clin %>%
  select(patient_id, age_at_index, gender) %>%
  distinct(patient_id, .keep_all = TRUE)

merged <- scani %>%
  inner_join(clin_uniq, by = "patient_id") %>%
  mutate(
    time = as.numeric(OS_time),
    event = as.numeric(OS_event),
    age = as.numeric(age_at_index),
    sex = factor(gender),
    stage_group = factor(stage_simple, levels = c("II","III","IV"))
  ) %>%
  filter(!is.na(stage_group), !is.na(time), !is.na(event), !is.na(age))

# 重新 Z-score 三组分（在 TCGA 内）
zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)

merged$Z_myCAF <- zsc(merged$myCAF)
merged$Z_TIMP1 <- zsc(merged$TIMP1)
merged$Z_PDPN  <- zsc(merged$PDPN)

# 三种 SCANI 变体
merged$SCANI_full   <- merged$Z_myCAF + merged$Z_TIMP1 - merged$Z_PDPN  # v89 公式
merged$SCANI_noPDPN <- merged$Z_myCAF + merged$Z_TIMP1                  # 去 PDPN
merged$SCANI_TIMP1only <- merged$Z_TIMP1                                # 仅 TIMP1

# 对比 v89 原版与重算版
cor_check <- cor(merged$SCANI_full, merged$CAF_Index, use = "pairwise.complete.obs")
cat("Sanity check: SCANI_full (重算) vs CAF_Index (v89 原版) correlation =",
    round(cor_check, 4), "\n")
if (cor_check < 0.99) {
  cat("⚠️ 重算 SCANI 与 CAF_Index 不完全一致！可能 v89 用了不同 Z 参数\n")
} else {
  cat("✅ 重算正确\n")
}

cat("\n--- 三种 SCANI 变体的多变量 Cox（含 stage + age + sex） ---\n")

run_cox <- function(score_col, label) {
  formula <- as.formula(paste("Surv(time, event) ~", score_col, "+ stage_group + age + sex"))
  model <- coxph(formula, data = merged)
  ci <- summary(model)$conf.int
  cf <- summary(model)$coefficients
  cat(sprintf("\n  %s:\n", label))
  cat(sprintf("    %s HR = %.3f, 95%%CI [%.3f, %.3f], p = %.3g\n",
              score_col,
              ci[score_col, "exp(coef)"],
              ci[score_col, "lower .95"],
              ci[score_col, "upper .95"],
              cf[score_col, "Pr(>|z|)"]))
  cat(sprintf("    Concordance = %.3f, AIC = %.1f\n",
              summary(model)$concordance[1],
              AIC(model)))
  return(list(model=model, ci=ci, cf=cf))
}

m_full   <- run_cox("SCANI_full",     "M2_full (Z(myCAF) + Z(TIMP1) - Z(PDPN))   [v89 公式]")
m_noPDPN <- run_cox("SCANI_noPDPN",   "M2_noPDPN (Z(myCAF) + Z(TIMP1))           [去 PDPN]")
m_only   <- run_cox("SCANI_TIMP1only","M2_TIMP1only (Z(TIMP1) alone)             [仅 TIMP1]")

# AIC 比较
cat("\n--- AIC 比较 (lower = better) ---\n")
cat("  SCANI_full:      AIC =", round(AIC(m_full$model), 2), "\n")
cat("  SCANI_noPDPN:    AIC =", round(AIC(m_noPDPN$model), 2), "\n")
cat("  SCANI_TIMP1only: AIC =", round(AIC(m_only$model), 2), "\n")

# 解读
cat("\n--- 解读 ---\n")
hr_full <- summary(m_full$model)$conf.int["SCANI_full", "exp(coef)"]
hr_no   <- summary(m_noPDPN$model)$conf.int["SCANI_noPDPN", "exp(coef)"]
delta_hr <- hr_full - hr_no
cat(sprintf("  HR change when removing PDPN: %.3f → %.3f  (Δ = %+.3f)\n",
            hr_full, hr_no, delta_hr))

if (abs(delta_hr) < 0.05) {
  cat("  → PDPN 对 SCANI 的边际贡献很小；论文应当承认 PDPN 是冗余\n")
} else if (delta_hr > 0) {
  cat("  → 加 PDPN（负号）使 SCANI 的 HR 更大；PDPN 残差信号有意义\n")
} else {
  cat("  → 加 PDPN 反而降低 HR；PDPN 项可能是噪声\n")
}

# =============================================================
# 任务 2: 找原始 50-gene signature
# =============================================================
cat("\n\n========================================================\n")
cat("  任务 2: 原始 50-gene signature 来源核查\n")
cat("========================================================\n")

GENESET_FILE <- "~/caf_spatial_project/results/03_spatial_GSE307534/CAF_Niche_final_geneset.txt"
if (file.exists(GENESET_FILE)) {
  geneset <- readLines(GENESET_FILE)
  cat(sprintf("Found %s\n  Number of genes: %d\n  First 10: %s\n",
              GENESET_FILE, length(geneset),
              paste(head(geneset, 10), collapse=", ")))
  cat(sprintf("  Last 5: %s\n", paste(tail(geneset, 5), collapse=", ")))
}

# 看 gene_frequency 在 50-gene signature 上
freq <- read.csv("~/caf_spatial_project/results/03_spatial_GSE307534/CAF_Niche_gene_frequency.csv",
                 stringsAsFactors = FALSE)
cat(sprintf("\nFull gene frequency table: nrow = %d\n", nrow(freq)))
cat("Columns:", paste(colnames(freq), collapse=", "), "\n")

# 50-gene 与 frequency table 的对应
in_freq <- geneset[geneset %in% freq$gene]
cat(sprintf("50-gene signature ∩ frequency table: %d / %d\n",
            length(in_freq), length(geneset)))

# 查看 50-gene signature 在 freq 中的实际 n_samples 分布
sub_freq <- freq[freq$gene %in% geneset, ]
sub_freq <- sub_freq[match(geneset, sub_freq$gene), ]
cat("\n50-gene signature 在 frequency table 中的 n_samples 分布:\n")
print(summary(sub_freq$n_samples))
cat("\nn_samples 分布表:\n")
print(table(sub_freq$n_samples))

# 真实 cutoff 验证：如果用不同 n_samples 阈值，会得到多少基因
cat("\n--- 真实 cutoff 验证 (基于 n_samples) ---\n")
cat("如果改 cutoff，得到的基因数（必须再加上 logFC 等其他筛选条件）：\n")
for (ct in c(6, 7, 8, 9, 10, 12, 15)) {
  n_above <- sum(freq$n_samples >= ct)
  cat(sprintf("  n_samples >= %d:  %d 个候选基因\n", ct, n_above))
}

cat("\n注意: 真正的 50-gene signature 是 n_samples >= 8 + 其他筛选(如 mean_logFC, p_adj)\n")
cat("      要做严格的 cutoff 验证，需要重跑 step3 spatial analysis 用不同阈值\n")

# 保存 v90 公式验证结果
cox_summary <- data.frame(
  model = c("SCANI_full (v89: Z(myCAF) + Z(TIMP1) - Z(PDPN))",
            "SCANI_noPDPN (Z(myCAF) + Z(TIMP1))",
            "SCANI_TIMP1only (Z(TIMP1))"),
  HR = c(round(summary(m_full$model)$conf.int["SCANI_full", "exp(coef)"], 3),
         round(summary(m_noPDPN$model)$conf.int["SCANI_noPDPN", "exp(coef)"], 3),
         round(summary(m_only$model)$conf.int["SCANI_TIMP1only", "exp(coef)"], 3)),
  CI_low = c(round(summary(m_full$model)$conf.int["SCANI_full", "lower .95"], 3),
             round(summary(m_noPDPN$model)$conf.int["SCANI_noPDPN", "lower .95"], 3),
             round(summary(m_only$model)$conf.int["SCANI_TIMP1only", "lower .95"], 3)),
  CI_high = c(round(summary(m_full$model)$conf.int["SCANI_full", "upper .95"], 3),
              round(summary(m_noPDPN$model)$conf.int["SCANI_noPDPN", "upper .95"], 3),
              round(summary(m_only$model)$conf.int["SCANI_TIMP1only", "upper .95"], 3)),
  p_value = c(format(summary(m_full$model)$coefficients["SCANI_full", "Pr(>|z|)"], scientific = TRUE, digits = 3),
              format(summary(m_noPDPN$model)$coefficients["SCANI_noPDPN", "Pr(>|z|)"], scientific = TRUE, digits = 3),
              format(summary(m_only$model)$coefficients["SCANI_TIMP1only", "Pr(>|z|)"], scientific = TRUE, digits = 3)),
  AIC = c(round(AIC(m_full$model), 2),
          round(AIC(m_noPDPN$model), 2),
          round(AIC(m_only$model), 2)),
  stringsAsFactors = FALSE
)
write.csv(cox_summary,
          file.path(OUT_DIR, "Suppl_Table_SZ_SCANI_formula_validation.csv"),
          row.names = FALSE)

cat("\n✅ Saved:", file.path(OUT_DIR, "Suppl_Table_SZ_SCANI_formula_validation.csv"), "\n")
print(cox_summary)
