#!/usr/bin/env Rscript
# ==============================================================================
# sci1_remark_extras.R —— SCI-1 改按 REMARK 报告后需要补的几项（只读计算，结果入补充表 S4）
#   REMARK 条目 6  中位随访（反向 KM）
#   REMARK 条目 13 队列基本特征：年龄、性别、分期、事件，含缺失数
#   REMARK 条目 14 SCANI 与标准预后变量的关系（年龄 Spearman；分期 Kruskal–Wallis；性别 Wilcoxon）
#   REMARK 条目 16/17 多变量模型中所有协变量的 HR 与 95% CI（不只报 SCANI）
#   REMARK 条目 18 比例风险假设检验（cox.zph，单变量与多变量模型）
#   不新增任何预测性能分析；不改动已报告的头条数字（脚本会断言头条数字可复现）
#
#   只读：tables/SCANI_TCGA_O_166_v2.csv、tables/SCANI_GSE72094_both_versions.csv、
#         caf_spatial_project/results/05_tcga_survival/TCGA_LUAD_niche_survival_v2.csv（性别）、
#         01_bases/A_bulk_OS/derived/GSE72094/GSE72094_clinical_clean.csv（性别）
#   只写（已存在即停止）：tables/remark_S4_characteristics_v2.csv、tables/remark_S4_models_v2.csv
#   日志：/tmp/sci1_remark_<时间>.txt
# 用法：Rscript ~/sci1_remark_extras.R
# ==============================================================================
suppressPackageStartupMessages(library(survival))
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
F1 <- file.path(T, "remark_S4_characteristics_v2.csv"); F2 <- file.path(T, "remark_S4_models_v2.csv")
for (f in c(F1, F2)) if (file.exists(f)) stop("已存在，不覆盖：", f)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_remark_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
zs <- function(x) as.numeric(scale(x))
revkm <- function(t, e) { s <- summary(survfit(Surv(t, 1 - e) ~ 1))$table; unname(c(s["median"], s["0.95LCL"], s["0.95UCL"])) }
fmt <- function(m) { s <- summary(m); data.frame(term = rownames(s$conf.int), HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
                                                 p = s$coefficients[, ncol(s$coefficients)], row.names = NULL) }
CH <- list(); MO <- list()
addch <- function(cohort, var, level, value) CH[[length(CH) + 1]] <<- data.frame(cohort, var, level, value)

# ============================== TCGA ==========================================
hr("TCGA-LUAD（II–IV 期，v2 同 aliquot）")
d <- read.csv(file.path(T, "SCANI_TCGA_O_166_v2.csv"), stringsAsFactors = FALSE)
d <- d[!duplicated(d$patient_id), ]; stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)
# 性别：建模临床表无此列，取自 TCGA-CDR 补充表 S1（Liu 2018；md5 a4591b2dcee39591f59e5e25a6ce75fa）
fcdr <- file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/raw/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx")
stopifnot(unname(tools::md5sum(fcdr)) == "a4591b2dcee39591f59e5e25a6ce75fa")
cdr <- as.data.frame(readxl::read_excel(fcdr, sheet = 1, guess_max = 20000)); cdr <- cdr[cdr$type %in% "LUAD", ]
stopifnot(!anyDuplicated(cdr$bcr_patient_barcode))
d$sex <- tolower(cdr$gender[match(d$patient_id, cdr$bcr_patient_barcode)])
cat("性别来自 TCGA-CDR；匹配到 ", sum(!is.na(d$sex)), "/166\n", sep = "")
d$stage <- factor(d$stage_simple, levels = c("II", "III", "IV"))
# 头条数字复现断言
u <- summary(coxph(Surv(OS_time, OS_event) ~ zs(SCANI), d))$conf.int[1, 1]
stopifnot(abs(u - 1.4629) < 1e-3)
m <- coxph(Surv(OS_time, OS_event) ~ zs(SCANI) + age_at_index + stage, d)
stopifnot(abs(summary(m)$conf.int[1, 1] - 1.4866) < 1e-3)
cat("头条数字复现：单变量 HR ", round(u, 4), "；多变量 HR ", round(summary(m)$conf.int[1, 1], 4), "\n", sep = "")
f <- revkm(d$OS_time, d$OS_event)
cat(sprintf("中位随访（反向 KM）：%.0f 天（%.0f–%.0f）\n", f[1], f[2], f[3]))
addch("TCGA", "patients/deaths", "", sprintf("%d/%d", nrow(d), sum(d$OS_event)))
addch("TCGA", "median follow-up, months (reverse KM, 95% CI)", "", sprintf("%.1f (%.1f–%.1f)", f[1] / 30.44, f[2] / 30.44, f[3] / 30.44))
addch("TCGA", "age, years", "median (IQR); missing", sprintf("%.0f (%.0f–%.0f); %d", median(d$age_at_index, na.rm = TRUE),
      quantile(d$age_at_index, .25, na.rm = TRUE), quantile(d$age_at_index, .75, na.rm = TRUE), sum(is.na(d$age_at_index))))
for (s in c(names(table(d$sex)), NA)) addch("TCGA", "sex", ifelse(is.na(s), "missing", s), sum(if (is.na(s)) is.na(d$sex) else d$sex %in% s))
for (s in levels(d$stage)) addch("TCGA", "stage", s, sum(d$stage %in% s))
print(do.call(rbind, CH), row.names = FALSE)
cat("\nSCANI 与标准预后变量：\n")
ra <- cor.test(d$SCANI, d$age_at_index, method = "spearman", exact = FALSE); kw <- kruskal.test(SCANI ~ stage, d)
cat(sprintf("  年龄 Spearman rho=%.3f P=%.3g；分期 Kruskal–Wallis P=%.3g\n", ra$estimate, ra$p.value, kw$p.value))
print(tapply(d$SCANI, d$stage, function(x) round(median(x), 2)))
if (length(unique(na.omit(d$sex))) == 2) { w <- wilcox.test(SCANI ~ sex, d); cat(sprintf("  性别 Wilcoxon P=%.3g\n", w$p.value))
  addch("TCGA", "SCANI vs sex", "Wilcoxon P", signif(w$p.value, 3)) }
addch("TCGA", "SCANI vs age", "Spearman rho (P)", sprintf("%.2f (%.2g)", ra$estimate, ra$p.value))
addch("TCGA", "SCANI vs stage", "Kruskal–Wallis P", signif(kw$p.value, 3))
cat("\n多变量模型（全部协变量）：\n"); x <- fmt(m); print(x, digits = 3); MO[[1]] <- cbind(cohort = "TCGA", model = "multivariable", x)
mu <- coxph(Surv(OS_time, OS_event) ~ zs(SCANI), d); MO[[2]] <- cbind(cohort = "TCGA", model = "univariable", fmt(mu))
cat("\n比例风险检验（cox.zph）：\n"); z1 <- cox.zph(mu); z2 <- cox.zph(m); print(z1$table); print(z2$table)
MO[[3]] <- data.frame(cohort = "TCGA", model = c("PH univariable SCANI", paste0("PH multivariable ", rownames(z2$table))),
                      term = c("zph", rep("zph", nrow(z2$table))), HR = NA, lo = NA, hi = NA, p = c(z1$table[1, "p"], z2$table[, "p"]))

# ============================== GSE72094 =====================================
hr("GSE72094（全基因组版 SCANI）")
g <- read.csv(file.path(T, "SCANI_GSE72094_both_versions.csv"), stringsAsFactors = FALSE)
cat("列名：", paste(names(g), collapse = ", "), "\n")
idc <- names(g)[1]; g <- g[!duplicated(g[[idc]]), ]; stopifnot(nrow(g) == 398)
sc <- "SCANI_full"; tt <- "OS_days"; ev <- "OS_event"; st <- "Stage"; ag <- "age"   # 列名按 2026-11-05 实际表头写死
stopifnot(all(c(sc, tt, ev, st, ag, "stage_num") %in% names(g)))
cat("使用列：SCANI=", sc, " time=", tt, " event=", ev, " stage=", st, " age=", ag, "\n", sep = "")
g$S <- g[[sc]]; g$t <- g[[tt]]; g$e <- g[[ev]]; stopifnot(sum(g$e) == 113)
# 分期：用表中原有的 stage_num（与已报告的多变量模型同一编码）；与 Stage 原值交叉核对
print(table(Stage = g[[st]], stage_num = g$stage_num, useNA = "ifany"))
g$age <- suppressWarnings(as.numeric(g[[ag]]))
cg <- read.csv(file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/GSE72094/GSE72094_clinical_clean.csv"), stringsAsFactors = FALSE)
cg <- cg[!duplicated(cg$gsm), ]; g$sex <- tolower(cg$sex[match(g[[idc]], cg$gsm)])
u <- summary(coxph(Surv(t, e) ~ zs(S), g))$conf.int[1, 1]; stopifnot(abs(u - 1.2334) < 1e-3)
cc <- complete.cases(g[, c("age", "stage_num")])
m2 <- coxph(Surv(t, e) ~ zs(S) + age + stage_num, g[cc, ])
cat(sprintf("头条数字复现：单变量 HR %.4f；多变量 n=%d/%d HR %.4f（已报告 1.20，n=393/111）\n", u, sum(cc), sum(g$e[cc]), summary(m2)$conf.int[1, 1]))
if (sum(cc) != 393 || abs(summary(m2)$conf.int[1, 1] - 1.2043) > 1e-3) cat("★ 多变量复现不一致：分期编码与原脚本不同，先别用下面的多变量结果，把这段贴回\n")
f <- revkm(g$t, g$e)
addch("GSE72094", "patients/deaths", "", sprintf("%d/%d", nrow(g), sum(g$e)))
addch("GSE72094", "median follow-up, months (reverse KM, 95% CI)", "", sprintf("%.1f (%.1f–%.1f)", f[1] / 30.44, f[2] / 30.44, f[3] / 30.44))
addch("GSE72094", "age, years", "median (IQR); missing", sprintf("%.0f (%.0f–%.0f); %d", median(g$age, na.rm = TRUE),
      quantile(g$age, .25, na.rm = TRUE), quantile(g$age, .75, na.rm = TRUE), sum(is.na(g$age))))
for (s in c(names(table(g$sex)), NA)) addch("GSE72094", "sex", ifelse(is.na(s), "missing", s), sum(if (is.na(s)) is.na(g$sex) else g$sex %in% s))
for (s in c(sort(unique(na.omit(g$stage_num))), NA)) addch("GSE72094", "stage", ifelse(is.na(s), "missing", paste0("stage ", s)), sum(if (is.na(s)) is.na(g$stage_num) else g$stage_num %in% s))
ra <- cor.test(g$S, g$age, method = "spearman", exact = FALSE); rs <- cor.test(g$S, g$stage_num, method = "spearman", exact = FALSE)
addch("GSE72094", "SCANI vs age", "Spearman rho (P)", sprintf("%.2f (%.2g)", ra$estimate, ra$p.value))
addch("GSE72094", "SCANI vs stage", "Spearman rho (P)", sprintf("%.2f (%.2g)", rs$estimate, rs$p.value))
if (length(unique(na.omit(g$sex))) == 2) addch("GSE72094", "SCANI vs sex", "Wilcoxon P", signif(wilcox.test(S ~ sex, g)$p.value, 3))
print(do.call(rbind, CH)[sapply(CH, function(r) r$cohort) == "GSE72094", ], row.names = FALSE)
cat("\n多变量模型（全部协变量）：\n"); x <- fmt(m2); print(x, digits = 3)
MO[[4]] <- cbind(cohort = "GSE72094", model = "multivariable", x)
MO[[5]] <- cbind(cohort = "GSE72094", model = "univariable", fmt(coxph(Surv(t, e) ~ zs(S), g)))
z1 <- cox.zph(coxph(Surv(t, e) ~ zs(S), g)); z2 <- cox.zph(m2)
cat("\n比例风险检验：\n"); print(z1$table); print(z2$table)
MO[[6]] <- data.frame(cohort = "GSE72094", model = c("PH univariable SCANI", paste0("PH multivariable ", rownames(z2$table))),
                      term = "zph", HR = NA, lo = NA, hi = NA, p = c(z1$table[1, "p"], z2$table[, "p"]))

write.csv(do.call(rbind, CH), F1, row.names = FALSE); write.csv(do.call(rbind, MO), F2, row.names = FALSE)
cat("\n写出：", F1, "\n      ", F2, "\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
