#!/usr/bin/env Rscript
# ==============================================================================
# sci1_aliquot_2665_check.R —— 只读：TCGA-44-2665 的 3 个 aliquot 各自的 myCAF ssGSEA 分数，
#   确认现行 O 口径保留的 myCAF（28749.32）来自哪个 aliquot，是否与单基因所用 aliquot
#   （TCGA-44-2665-01B-06R-A277-07）一致；并报告三种取法下的头条统计量（敏感性，不挑）。
#   ssgsea.norm=FALSE 时每个样本的分数只取决于该样本自身的基因排序，因此只对需要的列打分即可。
#   不写任何结果文件；日志 /tmp/sci1_aliquot_<时间>.txt
# 用法：Rscript ~/sci1_aliquot_2665_check.R
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(GSVA) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
S1T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
LOG <- file.path("/tmp", format(Sys.time(), "sci1_aliquot_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
cat("TCGA-44-2665 aliquot 核对  ", format(Sys.time()), "  GSVA ", as.character(packageVersion("GSVA")), "\n", sep = "")
stopifnot(packageVersion("GSVA") == "1.48.3")

sandbox <- function() {
  e <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
  for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "ggsave", "dev.off", "sink"))
    assign(f, noop, envir = e)
  e
}
e6 <- sandbox(); eval(parse(text = readLines(file.path(C, "scripts/step6d_genes.R"), warn = FALSE)[1:17]), envir = e6)
em <- get("expr_mat", e6); p2c <- get("pid_to_col", e6)
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]
stopifnot(length(sig) == 50)
cat("表达矩阵：", nrow(em), " 基因 × ", ncol(em), " 列；签名匹配 ", sum(sig %in% rownames(em)), "/50\n", sep = "")
ss <- function(cols) as.numeric(gsva(em[, cols, drop = FALSE], list(myCAF = sig), method = "ssgsea",
                                     ssgsea.norm = FALSE, verbose = FALSE)["myCAF", ])

hr("1. 预处理一致性：非重复病人的 myCAF 是否与 CAF_final_scores.csv 逐值一致")
sc <- read.csv(file.path(C, "results/05_tcga_survival/CAF_final_scores.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(sc) == 192)
s166 <- sc[!duplicated(sc$patient_id), ]; stopifnot(nrow(s166) == 166)
set.seed(1); chk <- sample(setdiff(s166$patient_id, "TCGA-44-2665"), 8)
v <- ss(unname(unlist(p2c[chk])))
print(data.frame(patient = chk, recomputed = round(v, 3), in_file = round(s166$myCAF[match(chk, s166$patient_id)], 3)), row.names = FALSE)
cat("最大绝对差：", signif(max(abs(v - s166$myCAF[match(chk, s166$patient_id)])), 3), "\n", sep = "")

hr("2. TCGA-44-2665 三个 aliquot 的 myCAF")
cols <- grep("^TCGA-44-2665", colnames(em), value = TRUE)
a <- data.frame(aliquot = cols, myCAF = round(ss(cols), 2))
a$used_for_TIMP1_PDPN <- a$aliquot == p2c[["TCGA-44-2665"]]
a$kept_in_O <- abs(a$myCAF - 28749.32) < 0.05
print(a, row.names = FALSE)

hr("3. 敏感性：该病人分别取三个 aliquot 的 myCAF（其余 165 例不变，Z 在 166 例上做）")
o <- read.csv(file.path(S1T, "SCANI_TCGA_O_166.csv"), stringsAsFactors = FALSE); stopifnot(nrow(o) == 166, sum(o$OS_event) == 80)
zs <- function(x) as.numeric(scale(x)); i <- which(o$patient_id == "TCGA-44-2665"); stopifnot(length(i) == 1)
o$stage <- factor(o$stage_simple, levels = c("II", "III", "IV"))
res <- do.call(rbind, lapply(seq_len(nrow(a)), function(k) {
  d <- o; d$myCAF[i] <- a$myCAF[k]; d$SCANI <- zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN); d$z <- zs(d$SCANI)
  u <- summary(coxph(Surv(OS_time, OS_event) ~ z, d)); m <- summary(coxph(Surv(OS_time, OS_event) ~ z + age_at_index + stage, d))
  g <- d$SCANI >= median(d$SCANI); kp <- 1 - pchisq(survdiff(Surv(OS_time, OS_event) ~ g, d)$chisq, 1)
  auc <- unname(timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = d$SCANI, cause = 1, times = 3 * 365, iid = FALSE)$AUC[2])
  data.frame(aliquot = a$aliquot[k], uni_HR = u$conf.int[1, 1], lo = u$conf.int[1, 3], hi = u$conf.int[1, 4], p = u$coefficients[1, 5],
             mv_HR = m$conf.int["z", 1], mv_p = m$coefficients["z", 5], KM_p = kp, AUC3 = auc)
}))
print(format(res, digits = 4), row.names = FALSE)
cat("\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
