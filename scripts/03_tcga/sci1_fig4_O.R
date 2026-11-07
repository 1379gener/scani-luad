#!/usr/bin/env Rscript
# ==============================================================================
# sci1_fig4_O.R —— SCI-1 Fig 4（O 口径：先去重、再在 166 例上 Z 标准化）+ TCGA 头条数字表
#
#   读入（只读）：CAF_final_scores.csv（192 行）、TCGA_LUAD_niche_survival_v2.csv（取 age）
#   写出（只写新目录 ~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/）：
#     tables/SCANI_TCGA_O_166.csv          166 例 O 口径 SCANI（供 verify 与 repo 使用）
#     tables/headline_numbers_TCGA_O.csv   单变量 / 多变量 / KM / AUC
#     figures/Fig4A.{tiff,png,pdf}         中位分组 KM
#     figures/Fig4B.{tiff,png,pdf}         3 年时间依赖 ROC
#   出图只读上面写出的 csv，脚本内无字面量数据向量。
#
# 用法：Rscript ~/sci1_fig4_O.R
# ==============================================================================

suppressPackageStartupMessages({ library(survival); library(ggplot2) })
H    <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS")
dir.create(file.path(OUTD, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTD, "tables"),  recursive = TRUE, showWarnings = FALSE)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_fig4_O_%Y%m%d_%H%M.txt"))
zz  <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
cat("SCI-1 Fig 4（O 口径）\n时间: ", format(Sys.time()), "\n日志: ", LOG, "\n", sep = "")

# ---------------------------------------------------------------- 1. 数据（先去重、断言）
sc <- read.csv(file.path(C, "results/05_tcga_survival/CAF_final_scores.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(sc) == 192)
d <- sc[!duplicated(sc$patient_id), ]
stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)

cl <- read.csv(file.path(C, "results/05_tcga_survival/TCGA_LUAD_niche_survival_v2.csv"), stringsAsFactors = FALSE)
cl <- cl[!duplicated(cl$patient_id), c("patient_id", "age_at_index")]
d  <- merge(d, cl, by = "patient_id", all.x = TRUE)
stopifnot(nrow(d) == 166, !anyNA(d$age_at_index))

zs <- function(x) as.numeric(scale(x))
d$SCANI <- zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)          # O：在 166 例上标准化
out <- d[, c("patient_id", "myCAF", "TIMP1", "PDPN", "SCANI", "OS_time", "OS_event", "stage_simple", "age_at_index")]
write.csv(out, file.path(OUTD, "tables/SCANI_TCGA_O_166.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 2. 头条数字（从刚写出的 csv 读回）
d <- read.csv(file.path(OUTD, "tables/SCANI_TCGA_O_166.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(d) == 166)
d$z     <- zs(d$SCANI)
d$stage <- factor(d$stage_simple, levels = c("II", "III", "IV"))
d$grp   <- factor(ifelse(d$SCANI >= median(d$SCANI), "High", "Low"), levels = c("Low", "High"))

u  <- summary(coxph(Surv(OS_time, OS_event) ~ z, d))
m  <- summary(coxph(Surv(OS_time, OS_event) ~ z + age_at_index + stage, d))
km <- survdiff(Surv(OS_time, OS_event) ~ grp, d); km_p <- 1 - pchisq(km$chisq, 1)
roc <- timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = d$SCANI,
                        cause = 1, times = 3 * 365, iid = FALSE)
auc <- unname(roc$AUC[2])

hn <- data.frame(
  analysis = c("univariable_perSD", "multivariable_age_stage_perSD", "km_median_split", "auc_3yr"),
  estimate = c(u$conf.int[1, 1], m$conf.int["z", 1], NA, auc),
  lower    = c(u$conf.int[1, 3], m$conf.int["z", 3], NA, NA),
  upper    = c(u$conf.int[1, 4], m$conf.int["z", 4], NA, NA),
  p        = c(u$coefficients[1, 5], m$coefficients["z", 5], km_p, NA),
  n        = c(u$n, m$n, nrow(d), nrow(d)),
  events   = c(u$nevent, m$nevent, sum(d$OS_event), sum(d$OS_event)))
write.csv(hn, file.path(OUTD, "tables/headline_numbers_TCGA_O.csv"), row.names = FALSE)
cat("\n头条数字（O 口径）:\n"); print(format(hn, digits = 4), row.names = FALSE)
cat("\n分组人数: "); print(table(d$grp))

# ---------------------------------------------------------------- 3. Fig 4A：KM（中位分组）
fit <- survfit(Surv(OS_time / 30.44, OS_event) ~ grp, d)
kmdf <- do.call(rbind, lapply(seq_along(fit$strata), function(i) {
  idx <- if (i == 1) seq_len(fit$strata[1]) else sum(fit$strata[1:(i - 1)]) + seq_len(fit$strata[i])
  data.frame(time = c(0, fit$time[idx]), surv = c(1, fit$surv[idx]),
             group = sub("grp=", "", names(fit$strata)[i]))
}))
ns <- table(d$grp)
kmdf$group <- factor(kmdf$group, levels = names(ns), labels = paste0(names(ns), " (n=", ns, ")"))
pA <- ggplot(kmdf, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = setNames(c("#4575B4", "#D73027"), levels(kmdf$group)), name = "SCANI") +
  scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x = max(kmdf$time) * 0.62, y = 0.95, size = 3.3,
           label = sprintf("Log-rank p = %.2f", km_p)) +
  labs(x = "Time (months)", y = "Overall survival probability") +
  theme_classic(base_size = 11) + theme(legend.position = "top")

# ---------------------------------------------------------------- 4. Fig 4B：3 年 ROC
rdf <- data.frame(FPR = roc$FP[, 2], TPR = roc$TP[, 2])
rdf <- rdf[order(rdf$FPR, rdf$TPR), ]
pB <- ggplot(rdf, aes(FPR, TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_step(colour = "#B2182B", linewidth = 0.8) +
  annotate("text", x = 0.62, y = 0.12, size = 3.3, label = sprintf("3-year AUC = %.3f", auc)) +
  coord_equal() + labs(x = "1 - specificity", y = "Sensitivity") +
  theme_classic(base_size = 11)

for (nm in c("Fig4A", "Fig4B")) {
  p <- if (nm == "Fig4A") pA else pB
  ggsave(file.path(OUTD, "figures", paste0(nm, ".tiff")), p, width = 5, height = 4.5, dpi = 600, compression = "lzw")
  ggsave(file.path(OUTD, "figures", paste0(nm, ".png")),  p, width = 5, height = 4.5, dpi = 600)
  ggsave(file.path(OUTD, "figures", paste0(nm, ".pdf")),  p, width = 5, height = 4.5)
}
cat("\n已写出 Fig4A / Fig4B 与两张表，目录: ", OUTD, "\n", sep = "")
cat("把 ", LOG, " 贴回 chat，并发 Fig4A.png / Fig4B.png。\n", sep = "")
