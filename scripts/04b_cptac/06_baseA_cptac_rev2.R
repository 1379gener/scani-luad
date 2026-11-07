#!/usr/bin/env Rscript
# ==============================================================================
# 06_baseA_cptac.R (rev2) —— 底座 A · CPTAC-LUAD【部分整备】
#
# 诊断结论（06a_cptac_diagnose.R）：
#   - os_data.csv 是 clinical_raw 长表的自笛卡尔积，不含任何生存时间；
#     OS_event 列唯一值只有 0（零事件）。该文件弃用。
#   - cell2020_OS.xlsx 的 Annotions_S1A（225 x 67）不含 OS/survival/vital 任何一列。
#   - 结论：CPTAC 目前没有可用的生存数据。
#
# 因此本脚本产出【部分整备】结果：
#   表达矩阵 + 分期 + 纯度/基质评分 + 人口学，OS_time / OS_event 明确留空为 NA。
#   manifest 中标记 has_survival = FALSE，供下游据此排除本队列的生存分析。
#
# 同时验证 Paper 1 的 CPTAC 分期独立性 n=72 是否来自本 xlsx。
#
# 输出：~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/CPTAC_LUAD/
# ==============================================================================

suppressPackageStartupMessages({ library(data.table) })
if (!requireNamespace("readxl", quietly = TRUE))
  stop("需要 readxl：install.packages('readxl')")

SRC <- path.expand("~/caf_spatial_project/data/CPTAC_LUAD")
D   <- path.expand("~/thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/CPTAC_LUAD")
dir.create(D, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(D, format(Sys.time(), "build_%Y%m%d_%H%M.log"))
sink(LOG, split = TRUE); on.exit(sink(), add = TRUE)

hr <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n", strrep("=", 68), "\n", sep="")
sh <- function(t) cat("\n  --- ", t, "\n", sep="")
kv <- function(k, v) cat(sprintf("  %-34s %s\n", paste0(k,":"), paste(v, collapse=", ")))

cat("底座 A · CPTAC-LUAD 部分整备\n时间: ", format(Sys.time()), "\n日志: ", LOG, "\n", sep="")
cat("\n⚠️  本队列无可用生存数据，仅整备表达 / 分期 / 纯度。见脚本头部说明。\n")

# ============================================================ 1. 表达矩阵
hr("1. 表达矩阵")
ex <- fread(file.path(SRC, "cptac_rna.cct"), sep = "\t", header = TRUE, showProgress = FALSE)
g  <- as.character(ex[[1]])
M  <- as.matrix(ex[, -1, with = FALSE]); if (!is.numeric(M)) storage.mode(M) <- "numeric"
rownames(M) <- g
kv("原始维度", paste(dim(M), collapse=" x "))
kv("重复基因名", sum(duplicated(g)))
if (anyDuplicated(g)) {                       # 统一规则：同名取均表达最大
  rm2 <- rowMeans(M, na.rm = TRUE)
  o <- order(rownames(M), -rm2); M <- M[o, , drop=FALSE]
  M <- M[!duplicated(rownames(M)), , drop=FALSE]
}
M <- M[order(rownames(M)), , drop=FALSE]
kv("去重后维度", paste(dim(M), collapse=" x "))
smp <- M[seq_len(min(2000,nrow(M))), seq_len(min(20,ncol(M))), drop=FALSE]
expr_med <- median(smp, na.rm = TRUE)
kv("取值范围（抽样）", sprintf("[%.3g, %.3g]  中位 %.3g",
    min(smp,na.rm=TRUE), max(smp,na.rm=TRUE), expr_med))
kv("NA 比例", sprintf("%.3f%%", 100*mean(is.na(M))))
cat("    -> 取值含负数且中位约 2.4，为已中心化的 log2 表达；不再做 log 变换\n")
ids <- colnames(M)
kv("样本 ID 形态", head(ids, 4))

# ============================================================ 2. xlsx 元数据
hr("2. cell2020_OS.xlsx / Annotions_S1A")
x <- suppressMessages(readxl::read_excel(file.path(SRC, "cell2020_OS.xlsx"),
                                          sheet = "Annotions_S1A"))
x <- as.data.table(x)
kv("维度", paste(dim(x), collapse=" x "))

sh("2.1 生存相关列复查（确认确实没有）")
sv_cols <- grep("surv|os[._ ]|overall|follow|vital|death|censor|status$|days|month|event",
                colnames(x), ignore.case = TRUE, value = TRUE)
kv("命中的列", ifelse(length(sv_cols), paste(sv_cols, collapse=" | "), "（无）"))
cat("    -> 上列若全为 *.mutation.status / QC.status / Smoking.Status 之类，\n")
cat("       则确认本表不含生存信息。\n")

sh("2.2 Type 分布（Tumor / NAT）")
tb <- table(x$Type, useNA = "ifany")
for (i in seq_along(tb)) cat(sprintf("    %-10s %4d\n", names(tb)[i], tb[i]))

sh("2.3 ID 形态与匹配")
kv("Sample.ID 示例", head(x$Sample.ID, 4))
kv("Participant 示例", head(x$Participant, 4))
kv("表达矩阵 ∩ Sample.ID", length(intersect(ids, x$Sample.ID)))
kv("表达矩阵 ∩ Participant", length(intersect(ids, x$Participant)))
cat("    -> Sample.ID 与表达列名同为带点格式，应可直接对齐，无需转换\n")

# ============================================================ 3. 验证 n=72
hr("3. 验证 Paper 1 的 CPTAC 分期独立性 n=72")
tum <- x[Type == "Tumor"]
kv("Type == Tumor", nrow(tum))
kv("其中 Sample.ID 在表达矩阵中", sum(tum$Sample.ID %in% ids))

sh("Stage 取值分布（Tumor 且在表达矩阵中）")
sel <- tum[Sample.ID %in% ids]
tb <- table(sel$Stage, useNA = "ifany")
for (i in seq_along(tb)) cat(sprintf("    %-8s %4d\n",
    ifelse(is.na(names(tb)[i]), "NA", names(tb)[i]), tb[i]))

n_stage <- sum(!is.na(sel$Stage) & nzchar(as.character(sel$Stage)))
sh("与 n=72 的核对")
cat(sprintf("    Tumor ∩ 表达矩阵 ∩ Stage 非空  =  %d\n", n_stage))
cat(sprintf("    Paper 1 报告                    =  72\n"))
cat(sprintf("    判定: %s\n", ifelse(n_stage == 72, "✅ 命中 —— Paper 1 的分期来自本 xlsx",
    sprintf("⚠️ 不符（差 %d）—— 分期来源需另查", n_stage - 72))))
if (n_stage != 72) {
  cat("\n    其他可能的口径：\n")
  cat(sprintf("      全部 Tumor 且 Stage 非空（不限表达矩阵）   = %d\n",
      sum(!is.na(tum$Stage) & nzchar(as.character(tum$Stage)))))
  cat(sprintf("      clinical_raw 的 STAGE 属性                  = 48（见诊断脚本）\n"))
  cat(sprintf("      表达矩阵样本数                              = %d\n", length(ids)))
}

# ============================================================ 4. 组装临床表
hr("4. 组装临床表（OS 明确留空）")
norm_stage <- function(s) {
  s <- toupper(trimws(as.character(s)))
  s[s %in% c("NA","","--",".")] <- NA
  s <- sub("^STAGE\\s*", "", s); s <- sub("[ABC]$", "", s)
  s <- unname(c("1"="I","2"="II","3"="III","4"="IV",
                "I"="I","II"="II","III"="III","IV"="IV")[s])
  s[!s %in% c("I","II","III","IV")] <- NA
  s
}

m <- match(ids, x$Sample.ID)
pick <- function(col) if (col %in% colnames(x)) x[[col]][m] else NA
cl <- data.table(
  sample_id   = ids,
  participant = as.character(pick("Participant")),
  type        = as.character(pick("Type")),
  stage_raw   = as.character(pick("Stage")),
  age         = suppressWarnings(as.numeric(pick("Age"))),
  sex         = tolower(as.character(pick("Gender"))),
  smoking     = as.character(pick("Smoking.Status")),
  purity_estimate_rna = suppressWarnings(as.numeric(pick("Tumor.Purity.byESTIMATE.RNAseq"))),
  purity_tsnet        = suppressWarnings(as.numeric(pick("TSNet Purity"))),
  estimate_score      = suppressWarnings(as.numeric(pick("ESTIMATEScore"))),
  estimate_immune     = suppressWarnings(as.numeric(pick("ESTIMATE ImmuneScore"))),
  estimate_stromal    = suppressWarnings(as.numeric(pick("ESTIMATE StromalScore"))),
  OS_time  = NA_real_,        # ★ 无可用来源
  OS_event = NA_integer_
)
cl[, stage := norm_stage(stage_raw)]
kv("匹配到 xlsx 的样本", sum(!is.na(m)))
kv("Type 分布", paste(sprintf("%s=%d", names(table(cl$type, useNA="ifany")),
                              table(cl$type, useNA="ifany")), collapse="  "))

sh("各字段的非缺失数")
for (cn in setdiff(colnames(cl), "sample_id"))
  cat(sprintf("    %-24s %4d / %d\n", cn, sum(!is.na(cl[[cn]])), nrow(cl)))

sh("分期归一后分布")
tb <- table(cl$stage, useNA = "ifany")
for (i in seq_along(tb)) cat(sprintf("    %-6s %4d\n",
    ifelse(is.na(names(tb)[i]),"NA",names(tb)[i]), tb[i]))

sh("★ 纯度字段可用性（Paper 2 实验二需要）")
for (p in c("purity_tsnet","purity_estimate_rna","estimate_stromal")) {
  v <- cl[[p]]; v <- v[!is.na(v)]
  if (length(v)) {
    cat(sprintf("    %-22s n=%3d  范围 [%.3f, %.3f]  中位 %.3f\n",
        p, length(v), min(v), max(v), median(v)))
  } else {
    cat(sprintf("    %-22s 全缺失\n", p))
  }
}
cat("\n    注：ESTIMATE StromalScore 不得作为主分析的纯度协变量\n")
cat("        （与被测 CAF signature 基因重叠，等同用自身校正自身）。\n")
cat("        主分析优先 TSNet Purity；ESTIMATE 系列仅作敏感性分析。\n")

# ============================================================ 5. 落盘
hr("5. 输出")
stopifnot(identical(colnames(M), cl$sample_id))
f_e <- file.path(D, "CPTAC_LUAD_expr_symbol.rds")
f_c <- file.path(D, "CPTAC_LUAD_clinical_clean.csv")
f_m <- file.path(D, "cohort_manifest_CPTAC.csv")
saveRDS(M, f_e)
write.csv(cl, f_c, row.names = FALSE)
write.csv(data.frame(
  cohort = "CPTAC-LUAD", platform = "RNA-seq (LinkedOmics .cct)",
  n_samples = ncol(M), n_genes = nrow(M),
  has_survival = FALSE,                       # ★ 下游据此排除本队列的生存分析
  survival_note = "os_data.csv 为长表自笛卡尔积，无生存时间且事件全为0；cell2020_OS.xlsx 不含 OS 列",
  n_with_stage = sum(!is.na(cl$stage)),
  n_with_purity_tsnet = sum(!is.na(cl$purity_tsnet)),
  clinical_source = "cell2020_OS.xlsx / Annotions_S1A",
  expr_scale = "log2, 已中心化（含负值）", need_log = FALSE,
  expr_median_sample = round(expr_med, 3),
  gene_collapse = "同名取均表达最大",
  source_md5_expr = tools::md5sum(file.path(SRC, "cptac_rna.cct")),
  source_md5_xlsx = tools::md5sum(file.path(SRC, "cell2020_OS.xlsx")),
  built = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE),
  f_m, row.names = FALSE)
for (f in c(f_e, f_c, f_m)) kv(basename(f), sprintf("%.2f MB", file.info(f)$size/1048576))

hr("6. 若要补齐 OS，可用的重下途径")
cat("  A. cBioPortal（研究 ID: luad_cptac_2020）\n")
cat("     https://www.cbioportal.org/study/summary?id=luad_cptac_2020\n")
cat("     下载 Clinical Data (TSV)，含 OS_STATUS / OS_MONTHS 两列（宽表，非长表）\n")
cat("     ★ 注意：现有的 clinical_raw.csv 就是从这里导出的长表，但缺 OS 两个属性\n\n")
cat("  B. Cell 2020 原文补充材料 Table S1B 及之后\n")
cat("     现有 xlsx 只含 S1A（样本元数据）。生存数据可能在 S1B/S1C\n\n")
cat("  C. GDC / PDC（Proteomic Data Commons）的 CPTAC-3 临床表\n\n")
cat("  在补齐之前，本队列不进入 Paper 2 的任何生存分析。\n")

hr("sessionInfo"); print(sessionInfo())
cat("\n日志: ", LOG, "\n", sep="")
