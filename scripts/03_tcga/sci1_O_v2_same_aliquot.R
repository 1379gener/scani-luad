#!/usr/bin/env Rscript
# ==============================================================================
# sci1_O_v2_same_aliquot.R —— SCI-1 TCGA 主分析 v2：每例病人的三个组分取自同一个 aliquot
#
#   规则（与 Methods 声明一致、事先固定，不看结果）：每例病人只用一个 aliquot，
#   即 step6d/step6e 的 pid_to_col（每病人第一列）；myCAF 也取该列的 ssGSEA 分数。
#   v1（SCANI_TCGA_O_166.csv）中只有 TCGA-44-2665 的 myCAF 取自另一个 aliquot（01A-01R-0946-07），
#   其余 165 例三个组分已同列。本脚本从表达矩阵现算该列的 myCAF（不手填），替换后在 166 例上 Z 标准化。
#
#   只读：TCGA_LUAD_SE.rds（经 step6d 第 1–17 行）、CAF_final_scores.csv、TCGA_LUAD_niche_survival_v2.csv、CAF_signatures.csv
#   只写（已存在即停止，不覆盖 v1）：
#     tables/SCANI_TCGA_O_166_v2.csv、tables/headline_numbers_TCGA_O_v2.csv
#   日志：/tmp/sci1_O_v2_<时间>.txt
# 用法：Rscript ~/sci1_O_v2_same_aliquot.R
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(GSVA) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTT <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
F_O2 <- file.path(OUTT, "SCANI_TCGA_O_166_v2.csv"); F_HN2 <- file.path(OUTT, "headline_numbers_TCGA_O_v2.csv")
for (f in c(F_O2, F_HN2)) if (file.exists(f)) stop("已存在，不覆盖：", f)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_O_v2_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
stopifnot(packageVersion("GSVA") == "1.48.3")
cat("SCI-1 O 口径 v2（同 aliquot）  ", format(Sys.time()), "\n", sep = "")

# ---- 表达矩阵与 aliquot 规则（照 step6d 第 1–17 行，写操作屏蔽）
e6 <- new.env(parent = globalenv()); for (f in c("write.csv", "saveRDS", "pdf", "png", "sink")) assign(f, function(...) invisible(NULL), envir = e6)
eval(parse(text = readLines(file.path(C, "scripts/step6d_genes.R"), warn = FALSE)[1:17]), envir = e6)
em <- get("expr_mat", e6); p2c <- get("pid_to_col", e6)
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)

# ---- 数据：先去重、断言
sc <- read.csv(file.path(C, "results/05_tcga_survival/CAF_final_scores.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(sc) == 192)
d <- sc[!duplicated(sc$patient_id), ]; stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)
cl <- read.csv(file.path(C, "results/05_tcga_survival/TCGA_LUAD_niche_survival_v2.csv"), stringsAsFactors = FALSE)
cl <- cl[!duplicated(cl$patient_id), c("patient_id", "age_at_index")]
d <- merge(d, cl, by = "patient_id", all.x = TRUE); stopifnot(nrow(d) == 166, !anyNA(d$age_at_index))

# ---- 每例病人的 myCAF 都按其 pid_to_col 列现算，与 v1 逐例比较（预期只有 TCGA-44-2665 不同）
cols <- unname(unlist(p2c[d$patient_id])); stopifnot(length(cols) == 166, !anyNA(cols))
my_same <- as.numeric(gsva(em[, cols], list(myCAF = sig), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)["myCAF", ])
dif <- abs(my_same - d$myCAF); chg <- d$patient_id[dif > 1e-6]
cat("myCAF 与 v1 不同的病人：", if (length(chg)) paste(chg, collapse = ", ") else "无", "\n", sep = "")
stopifnot(identical(chg, "TCGA-44-2665"))
cat("TCGA-44-2665：v1 ", round(d$myCAF[d$patient_id == chg], 2), " → v2 ", round(my_same[d$patient_id == chg], 2),
    "（aliquot ", p2c[[chg]], "）\n", sep = "")
d$myCAF <- my_same

zs <- function(x) as.numeric(scale(x))
d$SCANI <- zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)
out <- d[, c("patient_id", "myCAF", "TIMP1", "PDPN", "SCANI", "OS_time", "OS_event", "stage_simple", "age_at_index")]
out$aliquot <- cols
write.csv(out, F_O2, row.names = FALSE)

# ---- 头条数字（从刚写出的文件读回）
d <- read.csv(F_O2, stringsAsFactors = FALSE); stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)
d$stage <- factor(d$stage_simple, levels = c("II", "III", "IV"))
uni <- function(x) { s <- summary(coxph(Surv(OS_time, OS_event) ~ zs(x), d)); c(s$conf.int[1, c(1, 3, 4)], s$coefficients[1, 5]) }
auc <- function(x) unname(timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = x, cause = 1, times = 3 * 365, iid = FALSE)$AUC[2])
m <- summary(coxph(Surv(OS_time, OS_event) ~ zs(SCANI) + age_at_index + stage, d))
g <- d$SCANI >= median(d$SCANI); kp <- 1 - pchisq(survdiff(Surv(OS_time, OS_event) ~ g, d)$chisq, 1)
hn <- rbind(
  data.frame(analysis = "SCANI_univariable_perSD", t(setNames(uni(d$SCANI), c("estimate", "lower", "upper", "p")))),
  data.frame(analysis = "SCANI_multivariable_age_stage_perSD", estimate = m$conf.int[1, 1], lower = m$conf.int[1, 3], upper = m$conf.int[1, 4], p = m$coefficients[1, 5]),
  data.frame(analysis = "SCANI_km_median_split", estimate = NA, lower = NA, upper = NA, p = kp),
  data.frame(analysis = "SCANI_auc_3yr", estimate = auc(d$SCANI), lower = NA, upper = NA, p = NA),
  data.frame(analysis = "myCAF_univariable_perSD", t(setNames(uni(d$myCAF), c("estimate", "lower", "upper", "p")))),
  data.frame(analysis = "myCAF_auc_3yr", estimate = auc(d$myCAF), lower = NA, upper = NA, p = NA),
  data.frame(analysis = "myCAF_PDPN_spearman", estimate = cor(d$myCAF, d$PDPN, method = "spearman"), lower = NA, upper = NA, p = NA))
write.csv(hn, F_HN2, row.names = FALSE)
print(format(hn, digits = 4), row.names = FALSE)
cat("\n（TIMP1、PDPN 单基因结果不变：两者在 v1 中已取自同一 aliquot。）\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
