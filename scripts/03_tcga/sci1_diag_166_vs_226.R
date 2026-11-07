#!/usr/bin/env Rscript
# ==============================================================================
# sci1_diag_166_vs_226.R —— 只读诊断：为什么 TCGA II–IV 期按 CDR 有 226 例（115 死亡），分析表只有 166 例（80 死亡）
#   1) 找出 226 中不在 166 里的病人，看他们在原始临床表 TCGA_LUAD_niche_survival_v2.csv / CAF_final_scores.csv 中是否存在、分期和生存是什么
#   2) 对比 166 例：CDR 分期 vs 分析表分期
#   3) 搜索哪个脚本写出了 TCGA_LUAD_niche_survival_v2.csv，打印其中含 filter / subset / stage / OS 的行
#   4) 附带：GSE72094 / TCGA 中 myCAF 单独 HR 在随机 50 基因集单独 HR 分布中的位置
# 只读；不写任何文件；日志写 /tmp/sci1_diag_166_vs_226_<时间>.txt
# 用法：Rscript ~/sci1_diag_166_vs_226.R
# ==============================================================================
suppressPackageStartupMessages(library(survival))
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables"); RP <- file.path(T, "reviewer_prep")
LOG <- file.path("/tmp", format(Sys.time(), "sci1_diag_166_vs_226_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
stage_roman <- function(s) { s <- toupper(gsub("STAGE|\\s", "", toupper(as.character(s)))); s <- sub("[ABC][0-9]?$", "", s)
  ifelse(s %in% c("I", "II", "III", "IV"), s, NA) }

tc <- read.csv(file.path(T, "SCANI_TCGA_O_166_v2.csv"), stringsAsFactors = FALSE); tc <- tc[!duplicated(tc$patient_id), ]
stopifnot(nrow(tc) == 166, sum(tc$OS_event) == 80)
al <- read.csv(file.path(RP, "reviewer_C_TCGA_allstage_scores.csv"), stringsAsFactors = FALSE); stopifnot(!anyDuplicated(al$patient_id))
cat("全期别表：", nrow(al), " 例；分期：\n", sep = ""); print(table(al$stage4, useNA = "ifany"))
a24 <- al[al$stage4 %in% c("II", "III", "IV"), ]; stopifnot(nrow(a24) == 226)
extra <- a24[!a24$patient_id %in% tc$patient_id, ]; miss <- tc[!tc$patient_id %in% a24$patient_id, ]
cat("CDR II–IV 期 226 例中不在 166 里的：", nrow(extra), " 例（死亡 ", sum(extra$e), "）；166 中不在 226 里的：", nrow(miss), "\n", sep = "")

hr("1. 这些病人在原始临床表 / 打分表中的情况")
nf <- file.path(C, "results/05_tcga_survival/TCGA_LUAD_niche_survival_v2.csv"); ns <- read.csv(nf, stringsAsFactors = FALSE, check.names = FALSE)
cat("TCGA_LUAD_niche_survival_v2.csv：", nrow(ns), " 行；列名：\n  ", paste(names(ns), collapse = ", "), "\n", sep = "")
idc <- intersect(c("patient_id", "bcr_patient_barcode", "submitter_id", "case_submitter_id", "sample"), names(ns))[1]
cat("病人 id 列：", idc, "；唯一病人数：", length(unique(substr(ns[[idc]], 1, 12))), "\n", sep = "")
ns$pid12 <- substr(ns[[idc]], 1, 12); nsu <- ns[!duplicated(ns$pid12), ]
for (k in grep("stage|OS|os_|vital|days|time|event|niche|group|cluster|purity|sample_type|subset", names(nsu), ignore.case = TRUE, value = TRUE)) {
  v <- nsu[[k]]; if (length(unique(v)) <= 15) { cat("\n列 ", k, " 分布：\n", sep = ""); print(table(v, useNA = "ifany")) } else cat("\n列 ", k, "：连续/多值，缺失 ", sum(is.na(v)), "\n", sep = "") }
cat("\n226 中多出的 ", nrow(extra), " 例：在 niche_survival_v2 中 ", sum(extra$patient_id %in% nsu$pid12), " 例\n", sep = "")
cf <- read.csv(file.path(C, "results/05_tcga_survival/CAF_final_scores.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cid <- intersect(c("patient_id", "bcr_patient_barcode", "submitter_id", "sample"), names(cf))[1]
cat("CAF_final_scores.csv：", nrow(cf), " 行，唯一病人 ", length(unique(substr(cf[[cid]], 1, 12))), "；多出的病人中在打分表里的：",
    sum(extra$patient_id %in% substr(cf[[cid]], 1, 12)), "\n", sep = "")
inn <- extra[extra$patient_id %in% nsu$pid12, ]
if (nrow(inn)) { cat("\n多出病人在 niche_survival_v2 中的记录（前 30 例）：\n")
  print(head(merge(inn[, c("patient_id", "stage4", "t", "e")], nsu, by.x = "patient_id", by.y = "pid12"), 30)) }
cat("\n多出病人的 CDR 分期 / 生存概况：\n"); print(table(extra$stage4)); print(summary(extra$t)); cat("死亡：", sum(extra$e), "\n")
cat("\n多出病人列表：\n"); print(extra$patient_id)

hr("2. 166 例：CDR 分期 vs 分析表分期")
m <- match(tc$patient_id, al$patient_id); print(table(CDR = al$stage4[m], analysis = tc$stage_simple, useNA = "ifany"))

hr("3. 哪个脚本写出了 TCGA_LUAD_niche_survival_v2.csv / 被哪些脚本读入")
sc <- list.files(file.path(C, "scripts"), pattern = "\\.(R|r|py|sh)$", full.names = TRUE, recursive = FALSE)
for (f in sc) { L <- readLines(f, warn = FALSE); h <- grep("niche_survival_v2", L)
  if (length(h)) { cat("\n--- ", basename(f), "（第 ", paste(h, collapse = ","), " 行提到）\n", sep = "")
    k <- grep("filter|subset|stage|Stage|OS|os_time|vital|days_to|niche|write|>=|<=|%in%|!is.na|complete", L)
    for (i in k[k <= max(h) + 5]) cat(sprintf("%4d: %s\n", i, substr(L[i], 1, 180))) } }

hr("4. myCAF 单独 HR 在随机 50 基因集单独 HR 分布中的位置")
A <- read.csv(file.path(RP, "reviewer_A_component_variants.csv")); R <- read.csv(file.path(RP, "reviewer_E_random_draws.csv"))
for (cn in c("TCGA", "GSE72094")) { o <- A$HR[A$cohort == cn & A$variant == "myCAF_only"]; r <- R[R$cohort == cn, ]
  cat(sprintf("%s：myCAF 单独 HR %.3f；随机集单独 HR ≥ 该值的比例 %.3f；|log HR| ≥ 该值的比例 %.3f（n=%d）\n", cn, o,
              mean(r$rand_alone_HR >= o), mean(abs(log(r$rand_alone_HR)) >= abs(log(o))), nrow(r))) }
cat("\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
