#!/usr/bin/env Rscript
# ==============================================================================
# sci1_estimate_corr.R —— SCANI 及各组分与 ESTIMATE 基质 / 免疫评分（及 GSE72094 纯度）的 Spearman 相关，写成 csv 供稿件读取
#   TCGA-LUAD：CDR 226 例（tables/cdr226/SCANI_TCGA_CDR226.csv），ESTIMATE 用 tidyestimate 在 step6d 的 expr_mat 上重算（RNA-seq，无纯度）
#   GSE72094：398 例（tables/SCANI_GSE72094_both_versions.csv 的全基因组版），ESTIMATE 评分取 reviewer_prep/reviewer_D_estimate_scores.csv；
#             纯度按 ESTIMATE 芯片公式 cos(0.6049872018 + 0.0001467884 × ESTIMATE score) 由 stromal + immune 计算（Yoshihara 2013）
#   只读上述文件与 step6d_genes.R 第 1–17 行（写操作屏蔽）
#   只写（已存在即停止）：tables/cdr226/estimate_correlations.csv、tables/cdr226/estimate_scores_TCGA_CDR226.csv
#   日志：/tmp/sci1_estimate_corr_<时间>.txt
# 用法：Rscript ~/sci1_estimate_corr.R
# ==============================================================================
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables"); OD <- file.path(T, Sys.getenv("CDR_SUB", "cdr226"))
F1 <- file.path(OD, "estimate_correlations.csv"); F2 <- file.path(OD, "estimate_scores_TCGA_CDR226.csv")
for (f in c(F1, F2)) if (file.exists(f)) stop("已存在，不覆盖：", f)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_estimate_corr_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
stopifnot(requireNamespace("tidyestimate", quietly = TRUE))
zs <- function(x) as.numeric(scale(x))

# ---- TCGA 226
d <- read.csv(file.path(OD, "SCANI_TCGA_CDR226.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(d) == 226, !anyDuplicated(d$patient_id), sum(d$OS_event) == 115)
stopifnot(max(abs(d$SCANI - (zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)))) < 1e-8)
e6 <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "dev.off", "sink")) assign(f, noop, envir = e6)
eval(parse(text = readLines(file.path(C, "scripts/step6d_genes.R"), warn = FALSE)[1:17]), envir = e6)
EM <- get("expr_mat", e6); stopifnot(all(d$aliquot %in% colnames(EM)))
stopifnot(max(abs(as.numeric(EM["TIMP1", d$aliquot]) - d$TIMP1)) < 1e-8)
f <- tidyestimate::filter_common_genes(as.data.frame(EM[, d$aliquot]), id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE)
es <- tidyestimate::estimate_score(f, is_affymetrix = FALSE); print(head(es))
gcol <- function(df, p) df[[grep(p, names(df), ignore.case = TRUE, value = TRUE)[1]]]
m <- match(d$aliquot, as.character(gcol(es, "^sample"))); stopifnot(!anyNA(m))
d$stromal <- gcol(es, "^stromal")[m]; d$immune <- gcol(es, "^immune")[m]
write.csv(d[, c("patient_id", "aliquot", "stromal", "immune")], F2, row.names = FALSE); cat("写出 ", F2, "\n", sep = "")

# ---- GSE72094
g <- read.csv(file.path(T, "SCANI_GSE72094_both_versions.csv"), stringsAsFactors = FALSE); g <- g[!duplicated(g$sample_id), ]
stopifnot(nrow(g) == 398, sum(g$OS_event) == 113)
g$myCAF <- g$myCAF_full; g$SCANI <- g$SCANI_full
stopifnot(max(abs(g$SCANI - (zs(g$myCAF) + zs(g$TIMP1) - zs(g$PDPN)))) < 1e-8)
ge <- read.csv(file.path(T, "reviewer_prep/reviewer_D_estimate_scores.csv"), stringsAsFactors = FALSE); ge <- ge[ge$cohort == "GSE72094", ]
m <- match(g$sample_id, ge$id); stopifnot(!anyNA(m))
g$stromal <- ge$stromal[m]; g$immune <- ge$immune[m]
g$purity <- cos(0.6049872018 + 0.0001467884 * (g$stromal + g$immune))

# ---- 相关
rows <- list()
for (cn in c("TCGA-LUAD", "GSE72094")) { x <- if (cn == "TCGA-LUAD") d else g
  refs <- c("stromal", "immune", if (cn == "GSE72094") "purity")
  for (v in c("SCANI", "myCAF", "TIMP1", "PDPN")) for (r in refs) {
    ct <- cor.test(x[[v]], x[[r]], method = "spearman", exact = FALSE)
    rows[[length(rows) + 1]] <- data.frame(cohort = cn, n = nrow(x), score = v, estimate = r, rho = unname(ct$estimate), p = ct$p.value) } }
R <- do.call(rbind, rows); print(format(R, digits = 3), row.names = FALSE)
write.csv(R, F1, row.names = FALSE); cat("写出 ", F1, "\n完成。把 ", F1, " 发回（或 cat 打出来贴回）。\n", sep = "")
