#!/usr/bin/env Rscript
# ==============================================================================
# sci1_reviewer_prep.R —— 投稿前预演审稿意见（只读计算；结果先不进稿件）
#   对应模拟审稿意见：
#     A  M1/M2  组分单独与组分变体（SCANI、TIMP1、myCAF、PDPN、去 PDPN、PDPN 正号……）
#               两个队列：单变量 HR/SD、校正年龄+分期 HR、C-index、3 年 AUC
#     B  M1     嵌套似然比检验：TIMP1 单独 vs TIMP1 + [Z(myCAF) − Z(PDPN)]（两个队列）
#     C  M4     分期人群：TCGA 全期别（I–IV，用 TCGA-CDR 生存）与 I 期、II–IV 期；GSE72094 I 期与 II–IV 期
#     D  M5     混杂：加入性别；ESTIMATE 基质/免疫评分（tidyestimate，若已安装）；
#               GSE72094 的 KRAS/STK11/TP53/吸烟（若临床表里有）
#     E  M1     随机基因集基准（Venet 2011）：随机 50 基因替换 myCAF 组分，及随机 50 基因单独
#     F  M3     函数形态：自然样条（df=3）与线性比较的似然比检验（只看形态，不做分时段 PH 分析）
#     G  minor  3 年时仍在风险中的人数；TCGA-CDR 生存与分析表生存的一致性（166 例）
#
#   规矩：只读源文件；输出已存在即停止、不覆盖；日志写 /tmp；读入队列表先去重并 stopifnot 断言样本量；
#         先复现已报告的头条数字，复现不了就停。各节用 tryCatch 隔开，一节失败不影响其他节。
#   只读：tables/SCANI_TCGA_O_166_v2.csv、tables/SCANI_GSE72094_both_versions.csv、
#         caf_spatial_project/scripts/step6d_genes.R（第 1–17 行，写操作屏蔽）、
#         caf_spatial_project/scripts/step22_GSE72094_validation.R（第 1–145 行，写操作屏蔽）、
#         results/01_scRNA/CAF_signatures.csv、TCGA-CDR-SupplementalTableS1.xlsx、
#         01_bases/A_bulk_OS/derived/GSE72094/GSE72094_clinical_clean.csv、data/GSE72094_parsed/clinical.csv
#   只写（已存在即停止）：tables/reviewer_prep/ 目录下 reviewer_*.csv
#   日志：/tmp/sci1_reviewer_prep_<时间>.txt
# 用法：Rscript ~/sci1_reviewer_prep.R          （随机基因集默认 1000 组；N_RAND=200 Rscript ... 可先试跑）
# 预计耗时：随机基因集一节最久（两个队列各 1000 组 ssGSEA），用 3 个进程
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(splines) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
N_RAND <- as.integer(Sys.getenv("N_RAND", "1000"))
OD <- file.path(T, if (N_RAND == 1000) "reviewer_prep" else paste0("reviewer_prep_test_N", N_RAND))   # 试跑写到另一个目录
if (dir.exists(OD) && length(list.files(OD))) stop("已存在且非空，不覆盖：", OD)
dir.create(OD, showWarnings = FALSE)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_reviewer_prep_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
NCORE <- 3L; set.seed(20261106)
hr <- function(t) cat("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74), "\n", sep = "")
zs <- function(x) as.numeric(scale(x))
wout <- function(x, name) { f <- file.path(OD, name); if (file.exists(f)) stop("已存在，不覆盖：", f)
  write.csv(x, f, row.names = FALSE); cat("  → 写出 ", f, "\n", sep = "") }
section <- function(name, expr) { hr(name); t0 <- Sys.time()
  r <- tryCatch({ force(expr); "OK" }, error = function(e) { cat("★ 本节出错，已跳过：", conditionMessage(e), "\n"); "ERROR" })
  cat(sprintf("[%s] %s，用时 %.1f 分钟\n", name, r, as.numeric(difftime(Sys.time(), t0, units = "mins")))); invisible(r) }
cat("SCI-1 审稿预演  ", format(Sys.time()), "\nR ", R.version.string, "；survival ", as.character(packageVersion("survival")), "\n", sep = "")
stopifnot(requireNamespace("GSVA", quietly = TRUE), packageVersion("GSVA") == "1.48.3", requireNamespace("timeROC", quietly = TRUE))

# ---------------------------------------------------------------- 工具函数
D3 <- 3 * 365
cox1 <- function(t, e, x) { s <- summary(coxph(Surv(t, e) ~ zs(x)))
  c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
coxadj <- function(dd, x, covs) {           # dd 已是完整病例；x 在 dd 内重新标准化
  dd$x. <- zs(x); s <- summary(coxph(as.formula(paste("Surv(t, e) ~ x. +", paste(covs, collapse = " + "))), data = dd))
  c(HR = s$conf.int["x.", 1], lo = s$conf.int["x.", 3], hi = s$conf.int["x.", 4], p = s$coefficients["x.", 5], n = s$n, ev = s$nevent) }
cidx <- function(t, e, x) { cc <- concordance(Surv(t, e) ~ x, reverse = TRUE); c(C = unname(cc$concordance), C_se = sqrt(unname(cc$var))) }
auc3 <- function(t, e, x) unname(timeROC::timeROC(T = t, delta = e, marker = x, cause = 1, times = D3, iid = FALSE)$AUC[2])
ssg <- function(mat, sets) GSVA::gsva(mat, sets, method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE, parallel.sz = NCORE)
sandbox <- function() { e <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
  for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "tiff", "jpeg", "ggsave",
              "dev.off", "fwrite", "write_csv", "sink")) assign(f, noop, envir = e)
  assign("dir.create", function(...) invisible(TRUE), envir = e); e }
run_lines <- function(file, upto, env) eval(parse(text = readLines(file, warn = FALSE)[seq_len(upto)]), envir = env)
stage_roman <- function(s) { s <- toupper(gsub("STAGE|\\s", "", toupper(as.character(s)))); s <- sub("[ABC][0-9]?$", "", s)
  ifelse(s %in% c("I", "II", "III", "IV"), s, NA) }
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)

# ---------------------------------------------------------------- 0. 读入与头条数字复现
hr("0. 读入两个分析表，复现已报告的头条数字")
tc <- read.csv(file.path(T, "SCANI_TCGA_O_166_v2.csv"), stringsAsFactors = FALSE)
tc <- tc[!duplicated(tc$patient_id), ]; stopifnot(nrow(tc) == 166, sum(tc$OS_event) == 80)
tc$t <- tc$OS_time; tc$e <- tc$OS_event; tc$stage <- factor(tc$stage_simple, levels = c("II", "III", "IV")); tc$age <- tc$age_at_index
stopifnot(max(abs(tc$SCANI - (zs(tc$myCAF) + zs(tc$TIMP1) - zs(tc$PDPN)))) < 1e-8)
stopifnot(abs(cox1(tc$t, tc$e, tc$SCANI)[["HR"]] - 1.4629) < 1e-3,
          abs(coxadj(tc, tc$SCANI, c("age", "stage"))[["HR"]] - 1.4866) < 1e-3)
gs <- read.csv(file.path(T, "SCANI_GSE72094_both_versions.csv"), stringsAsFactors = FALSE)
gs <- gs[!duplicated(gs$sample_id), ]; stopifnot(nrow(gs) == 398, sum(gs$OS_event) == 113)
gs$t <- gs$OS_days; gs$e <- gs$OS_event; gs$myCAF <- gs$myCAF_full; gs$SCANI <- gs$SCANI_full
gs$age <- suppressWarnings(as.numeric(gs$age))
stopifnot(max(abs(gs$SCANI - (zs(gs$myCAF) + zs(gs$TIMP1) - zs(gs$PDPN)))) < 1e-8)
gcc <- complete.cases(gs[, c("age", "stage_num")]); stopifnot(sum(gcc) == 393, sum(gs$e[gcc]) == 111)
stopifnot(abs(cox1(gs$t, gs$e, gs$SCANI)[["HR"]] - 1.2334) < 1e-3,
          abs(coxadj(gs[gcc, ], gs$SCANI[gcc], c("age", "stage_num"))[["HR"]] - 1.2043) < 1e-3)
cat("头条数字复现通过：TCGA 1.4629 / 1.4866；GSE72094 1.2334 / 1.2043（n=393/111）\n")
cat("GSE72094 分期分布（stage_num）：\n"); print(table(gs$stage_num, useNA = "ifany"))

# 组分变体：每个组分先在队列内 Z 标准化，再按符号相加
VAR <- list(SCANI = c(myCAF = 1, TIMP1 = 1, PDPN = -1), TIMP1_only = c(TIMP1 = 1), myCAF_only = c(myCAF = 1),
            PDPN_only_neg = c(PDPN = -1), myCAF_plus_TIMP1 = c(myCAF = 1, TIMP1 = 1), TIMP1_minus_PDPN = c(TIMP1 = 1, PDPN = -1),
            myCAF_minus_PDPN = c(myCAF = 1, PDPN = -1), all_three_PDPN_positive = c(myCAF = 1, TIMP1 = 1, PDPN = 1))
score <- function(d, w) Reduce(`+`, lapply(names(w), function(k) w[[k]] * zs(d[[k]])))
COH <- list(TCGA = list(d = tc, covs = c("age", "stage")), GSE72094 = list(d = gs, covs = c("age", "stage_num")))

# ---------------------------------------------------------------- A
section("A. 组分单独与组分变体（M1、M2）", {
  A <- do.call(rbind, lapply(names(COH), function(cn) { d <- COH[[cn]]$d; covs <- COH[[cn]]$covs
    cc <- complete.cases(d[, covs])
    do.call(rbind, lapply(names(VAR), function(v) { x <- score(d, VAR[[v]]); u <- cox1(d$t, d$e, x)
      a <- coxadj(d[cc, ], x[cc], covs); ci <- cidx(d$t, d$e, x)
      data.frame(cohort = cn, variant = v, n = nrow(d), events = sum(d$e), HR = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], p = u[["p"]],
                 adj_n = a[["n"]], adj_HR = a[["HR"]], adj_lo = a[["lo"]], adj_hi = a[["hi"]], adj_p = a[["p"]],
                 C = ci[["C"]], C_se = ci[["C_se"]], AUC3y = auc3(d$t, d$e, x)) })) }))
  print(format(A, digits = 3), row.names = FALSE)
  cat("核对：TCGA 的 SCANI / TIMP1_only / myCAF_only / PDPN 行应与 Table 1 一致（PDPN_only_neg 的 HR 为 1/0.84）\n")
  wout(A, "reviewer_A_component_variants.csv") })

# ---------------------------------------------------------------- B
section("B. 嵌套似然比检验：TIMP1 单独 vs TIMP1 + [Z(myCAF) − Z(PDPN)]（M1）", {
  B <- do.call(rbind, lapply(names(COH), function(cn) { d <- COH[[cn]]$d; covs <- COH[[cn]]$covs
    d$zT <- zs(d$TIMP1); d$zMP <- zs(d$myCAF) - zs(d$PDPN); d$zM <- zs(d$myCAF); d$zP <- zs(d$PDPN)
    m0 <- coxph(Surv(t, e) ~ zT, d); m1 <- coxph(Surv(t, e) ~ zT + zMP, d); m2 <- coxph(Surv(t, e) ~ zT + zM + zP, d)
    l1 <- anova(m0, m1); l2 <- anova(m0, m2); s1 <- summary(m1)$conf.int; s2 <- summary(m2)$conf.int
    cat(cn, "：TIMP1 + (myCAF − PDPN) 模型\n"); print(round(s1, 3)); cat(cn, "：TIMP1 + myCAF + PDPN 各自估计\n"); print(round(s2, 3))
    data.frame(cohort = cn, comparison = c("TIMP1 vs TIMP1 + [Z(myCAF)-Z(PDPN)] (1 df)", "TIMP1 vs TIMP1 + myCAF + PDPN (2 df)"),
               LR_chisq = c(l1$Chisq[2], l2$Chisq[2]), df = c(l1$Df[2], l2$Df[2]), p = c(l1[[ncol(l1)]][2], l2[[ncol(l2)]][2]),
               added_term_HR = c(s1["zMP", 1], NA), added_term_lo = c(s1["zMP", 3], NA), added_term_hi = c(s1["zMP", 4], NA),
               C_TIMP1 = unname(concordance(m0)$concordance), C_full = c(unname(concordance(m1)$concordance), unname(concordance(m2)$concordance))) }))
  print(format(B, digits = 3), row.names = FALSE)
  cat("注：TCGA 中的检验是乐观的（组分在此队列选出）；GSE72094 的结果才是独立证据\n")
  wout(B, "reviewer_B_nested_LRT.csv") })

# ---------------------------------------------------------------- 表达矩阵（C、D、E 共用）
EM <- NULL; P2C <- NULL; FULL <- NULL
section("表达矩阵：TCGA（step6d 第 1–17 行）与 GSE72094 全基因组（step22 第 1–145 行 + 同 symbol 取均值）", {
  e6 <- sandbox(); run_lines(file.path(C, "scripts/step6d_genes.R"), 17, e6)
  EM <<- get("expr_mat", e6); P2C <<- get("pid_to_col", e6)
  cat("TCGA expr_mat：", nrow(EM), " 基因 x ", ncol(EM), " 列；pid_to_col 病人数 ", length(P2C), "\n", sep = "")
  cat("列名样本类型（barcode 第 14–15 位）：\n"); print(table(substr(colnames(EM), 14, 15)))
  cat("expr_mat 列名示例：", paste(head(colnames(EM), 3), collapse = " | "), "\n", sep = "")
  p2n <- function(v) if (is.numeric(v)) colnames(EM)[v] else as.character(v)
  cols <- p2n(unname(unlist(P2C[tc$patient_id]))); stopifnot(length(cols) == 166, !anyNA(cols), all(cols %in% colnames(EM)))
  stopifnot(all(as.character(unname(unlist(P2C[tc$patient_id]))) == as.character(tc$aliquot)))
  tc$col <- cols
  stopifnot(max(abs(as.numeric(EM["TIMP1", cols]) - tc$TIMP1)) < 1e-8, max(abs(as.numeric(EM["PDPN", cols]) - tc$PDPN)) < 1e-8)
  my <- as.numeric(ssg(EM[, cols], list(myCAF = sig))["myCAF", ]); stopifnot(max(abs(my - tc$myCAF)) < 1e-6)
  cat("TCGA 166 例：aliquot、TIMP1、PDPN、myCAF 与 v2 表完全一致\n")
  e22 <- sandbox(); run_lines(file.path(C, "scripts/step22_GSE72094_validation.R"), 145, e22)
  annot <- get("annot", e22); lines <- get("lines", e22); su <- get("samples_use", e22)
  b <- grep("^!series_matrix_table_begin", lines); en <- grep("^!series_matrix_table_end", lines)
  mat <- read.delim(text = paste(lines[(b + 1):(en - 1)], collapse = "\n"), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(mat) <- mat[[1]]; mat <- as.matrix(mat[, -1]); storage.mode(mat) <- "numeric"
  sym <- annot$GeneSymbol[match(rownames(mat), annot$ID)]; keep <- !is.na(sym) & nzchar(sym)
  mat <- mat[keep, su, drop = FALSE]; sym <- sym[keep]
  rs <- rowsum(mat, sym); f <- rs / as.vector(table(sym)[rownames(rs)]); f <- f[rowSums(is.na(f)) == 0, ]
  f <- f[, gs$sample_id]; stopifnot(ncol(f) == 398)
  stopifnot(max(abs(f["TIMP1", ] - gs$TIMP1)) < 1e-8, max(abs(f["PDPN", ] - gs$PDPN)) < 1e-8)
  my2 <- as.numeric(ssg(f, list(myCAF = sig))["myCAF", ]); stopifnot(max(abs(my2 - gs$myCAF)) < 1e-6)
  FULL <<- f; cat("GSE72094：", nrow(f), " 基因 x 398 样本；TIMP1、PDPN、myCAF 与分析表完全一致\n", sep = "") })

# ---------------------------------------------------------------- C
CDR <- NULL
section("C. 分期人群（M4）：TCGA 全期别 I–IV（TCGA-CDR 生存）；GSE72094 I 期与 II–IV 期", {
  fcdr <- file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/raw/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx")
  stopifnot(unname(tools::md5sum(fcdr)) == "a4591b2dcee39591f59e5e25a6ce75fa")
  cdr <- as.data.frame(readxl::read_excel(fcdr, sheet = 1, guess_max = 20000)); cdr <- cdr[cdr$type %in% "LUAD", ]
  stopifnot(!anyDuplicated(cdr$bcr_patient_barcode))
  cdr$OS.time <- suppressWarnings(as.numeric(cdr$OS.time)); cdr$OS <- suppressWarnings(as.numeric(cdr$OS))
  cdr$age <- suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis)); cdr$sex <- tolower(cdr$gender)
  cdr$stage4 <- stage_roman(cdr$ajcc_pathologic_tumor_stage); CDR <<- cdr
  cat("TCGA-CDR LUAD：", nrow(cdr), " 例；分期：\n", sep = ""); print(table(cdr$stage4, useNA = "ifany"))
  # 每例病人取 pid_to_col 的列；不在 pid_to_col 中的病人取其第一个原发肿瘤（01）列（同一规则）
  p2n <- function(v) if (is.numeric(v)) colnames(EM)[v] else as.character(v)
  pts <- cdr$bcr_patient_barcode; col <- vapply(pts, function(p) {
    if (p %in% names(P2C)) return(p2n(P2C[[p]])[1])
    cs <- colnames(EM)[substr(colnames(EM), 1, 12) == p & substr(colnames(EM), 14, 15) == "01"]; if (length(cs)) cs[1] else NA_character_ }, "")
  cat("有表达列的病人：", sum(!is.na(col)), "；其中来自 pid_to_col ", sum(pts %in% names(P2C)), "\n", sep = "")
  cat("所用列的样本类型：\n"); print(table(substr(na.omit(col), 14, 15)))
  if (sum(!is.na(col)) < 400) cat("★ 有表达列的病人不足 400：step6d 的 expr_mat 可能只含 II–IV 期，全期别结果不可用，请把本节输出贴回\n")
  a <- data.frame(patient_id = pts, aliquot = col, t = cdr$OS.time, e = cdr$OS, age = cdr$age, sex = cdr$sex, stage4 = cdr$stage4, stringsAsFactors = FALSE)
  a <- a[!is.na(a$aliquot) & substr(a$aliquot, 14, 15) == "01" & !is.na(a$t) & a$t > 0 & !is.na(a$e), ]
  cat("全期别分析人群（有原发肿瘤表达、CDR 生存时间 >0）：", nrow(a), " 例，", sum(a$e), " 例死亡\n", sep = "")
  a$myCAF <- as.numeric(ssg(EM[, a$aliquot], list(myCAF = sig))["myCAF", ])
  a$TIMP1 <- as.numeric(EM["TIMP1", a$aliquot]); a$PDPN <- as.numeric(EM["PDPN", a$aliquot])
  a$SCANI <- score(a, VAR$SCANI)                     # 在全期别人群内 Z 标准化
  wout(a, "reviewer_C_TCGA_allstage_scores.csv")
  m <- match(tc$patient_id, a$patient_id)
  cat("166 例中能在全期别人群找到：", sum(!is.na(m)), "；myCAF 一致（最大差）：", signif(max(abs(a$myCAF[m] - tc$myCAF), na.rm = TRUE), 3), "\n", sep = "")
  cat("166 例 CDR 分期 vs 分析表分期：\n"); print(table(CDR = a$stage4[m], analysis = tc$stage_simple, useNA = "ifany"))
  rows <- list(); add <- function(lbl, d, x, covs) { u <- cox1(d$t, d$e, x); cc <- complete.cases(d[, covs])
    ad <- if (length(covs)) coxadj(d[cc, ], x[cc], covs) else c(HR = NA, lo = NA, hi = NA, p = NA, n = NA, ev = NA)
    rows[[length(rows) + 1]] <<- data.frame(analysis = lbl, n = nrow(d), events = sum(d$e), HR = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], p = u[["p"]],
      adj_covariates = paste(covs, collapse = "+"), adj_n = ad[["n"]], adj_HR = ad[["HR"]], adj_lo = ad[["lo"]], adj_hi = ad[["hi"]], adj_p = ad[["p"]]) }
  a$stageF <- factor(a$stage4, levels = c("I", "II", "III", "IV"))
  add("TCGA all stages I-IV, SCANI (Z within all-stage cohort)", a, a$SCANI, c("age", "stageF"))
  add("TCGA all stages I-IV, TIMP1 only", a, a$TIMP1, c("age", "stageF"))
  for (s in list(c("I"), c("II", "III", "IV"))) { d <- a[a$stage4 %in% s, ]; d$stageF <- droplevels(factor(d$stage4, levels = c("I", "II", "III", "IV")))
    covs <- if (length(s) > 1) c("age", "stageF") else "age"; lb <- paste0("TCGA stage ", paste(s, collapse = "/"))
    add(paste0(lb, ", SCANI re-standardized in subset"), d, score(d, VAR$SCANI), covs)
    add(paste0(lb, ", TIMP1 only"), d, d$TIMP1, covs) }
  for (s in list(1, 2:4)) { d <- gs[gs$stage_num %in% s, ]; covs <- if (length(s) > 1) c("age", "stage_num") else "age"
    lb <- paste0("GSE72094 stage ", paste(as.roman(s), collapse = "/"))
    add(paste0(lb, ", SCANI re-standardized in subset"), d, score(d, VAR$SCANI), covs)
    add(paste0(lb, ", SCANI from full-cohort standardization"), d, d$SCANI, covs)
    add(paste0(lb, ", TIMP1 only"), d, d$TIMP1, covs) }
  Cc <- do.call(rbind, rows); print(format(Cc, digits = 3), row.names = FALSE)
  wout(Cc, "reviewer_C_stage_subsets.csv") })

# ---------------------------------------------------------------- D
section("D. 混杂（M5）：性别、ESTIMATE 基质/免疫评分、GSE72094 突变与吸烟", {
  rows <- list(); add <- function(cn, lbl, d, covs) { cc <- complete.cases(d[, covs]); a <- coxadj(d[cc, ], d$SCANI[cc], covs)
    rows[[length(rows) + 1]] <<- data.frame(cohort = cn, model = lbl, covariates = paste(covs, collapse = "+"), n = a[["n"]], events = a[["ev"]],
                                             SCANI_HR = a[["HR"]], lo = a[["lo"]], hi = a[["hi"]], p = a[["p"]]) }
  # 性别
  fcdr <- file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/raw/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx")
  cdr <- if (is.null(CDR)) { x <- as.data.frame(readxl::read_excel(fcdr, sheet = 1, guess_max = 20000)); x <- x[x$type %in% "LUAD", ]; x$sex <- tolower(x$gender); x } else CDR
  tc$sex <- cdr$sex[match(tc$patient_id, cdr$bcr_patient_barcode)]; stopifnot(sum(!is.na(tc$sex)) == 166)
  cg <- read.csv(file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/GSE72094/GSE72094_clinical_clean.csv"), stringsAsFactors = FALSE)
  cg <- cg[!duplicated(cg$gsm), ]; gs$sex <- tolower(cg$sex[match(gs$sample_id, cg$gsm)])
  cat("GSE72094 性别：\n"); print(table(gs$sex, useNA = "ifany"))
  add("TCGA", "age + stage (reported)", tc, c("age", "stage")); add("TCGA", "+ sex", tc, c("age", "stage", "sex"))
  add("GSE72094", "age + stage (reported)", gs, c("age", "stage_num")); add("GSE72094", "+ sex", gs, c("age", "stage_num", "sex"))
  # ESTIMATE（tidyestimate，CRAN）
  if (requireNamespace("tidyestimate", quietly = TRUE) && !is.null(EM) && !is.null(FULL)) tryCatch({
    est <- function(mat, affy) { df <- as.data.frame(mat); f <- tidyestimate::filter_common_genes(df, id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE)
      tidyestimate::estimate_score(f, is_affymetrix = affy) }
    et <- est(EM[, tc$col], FALSE); eg <- est(FULL, TRUE); print(head(et)); print(head(eg))
    gc <- function(df, pat) { h <- grep(pat, names(df), ignore.case = TRUE, value = TRUE); if (length(h)) df[[h[1]]] else NULL }
    getv <- function(df, pat, ids) { v <- gc(df, pat); if (is.null(v)) return(NULL); v[match(ids, as.character(gc(df, "^sample")))] }
    tc$stromal <- getv(et, "^stromal", tc$col); tc$immune <- getv(et, "^immune", tc$col)
    gs$stromal <- getv(eg, "^stromal", gs$sample_id); gs$immune <- getv(eg, "^immune", gs$sample_id)
    pu <- getv(eg, "^purity", gs$sample_id); if (!is.null(pu)) gs$purity <- pu
    stopifnot(sum(!is.na(tc$stromal)) == 166, sum(!is.na(gs$stromal)) == 398)
    for (cn in c("TCGA", "GSE72094")) { d <- if (cn == "TCGA") tc else gs
      cat(sprintf("%s：SCANI 与 ESTIMATE 基质评分 Spearman %.3f；与免疫评分 %.3f；myCAF 与基质 %.3f；TIMP1 与基质 %.3f\n", cn,
                  cor(d$SCANI, d$stromal, method = "spearman"), cor(d$SCANI, d$immune, method = "spearman"),
                  cor(d$myCAF, d$stromal, method = "spearman"), cor(d$TIMP1, d$stromal, method = "spearman"))) }
    if (!is.null(gs$purity)) cat(sprintf("GSE72094：SCANI 与 ESTIMATE 纯度 Spearman %.3f（纯度公式只适用于芯片）\n", cor(gs$SCANI, gs$purity, method = "spearman")))
    tc$zstromal <- zs(tc$stromal); tc$zimmune <- zs(tc$immune); gs$zstromal <- zs(gs$stromal); gs$zimmune <- zs(gs$immune)
    add("TCGA", "+ ESTIMATE stromal + immune", tc, c("age", "stage", "zstromal", "zimmune"))
    add("GSE72094", "+ ESTIMATE stromal + immune", gs, c("age", "stage_num", "zstromal", "zimmune"))
    wout(rbind(data.frame(cohort = "TCGA", id = tc$patient_id, stromal = tc$stromal, immune = tc$immune),
               data.frame(cohort = "GSE72094", id = gs$sample_id, stromal = gs$stromal, immune = gs$immune)), "reviewer_D_estimate_scores.csv")
  }, error = function(e) cat("★ ESTIMATE 部分出错，已跳过（性别、突变部分照常）：", conditionMessage(e), "\n")) else cat("★ 未安装 tidyestimate（或表达矩阵一节失败），跳过 ESTIMATE。安装：install.packages('tidyestimate')，装好后需要把 reviewer_prep 目录改名再重跑\n")
  # GSE72094 突变 / 吸烟：在两个临床表里按列名搜
  cand <- list(cg = cg); f2 <- file.path(C, "data/GSE72094_parsed/clinical.csv")
  if (file.exists(f2)) { x <- read.csv(f2, stringsAsFactors = FALSE, check.names = FALSE); cand$parsed <- x }
  for (nm in names(cand)) cat("临床表 ", nm, " 列名：", paste(names(cand[[nm]]), collapse = ", "), "\n", sep = "")
  pick <- function(pat) { for (nm in names(cand)) { x <- cand[[nm]]; hit <- grep(pat, names(x), ignore.case = TRUE, value = TRUE)
    if (length(hit)) { idc <- intersect(c("gsm", "geo_accession", "sample_id", "GSM", "Sample"), names(x))[1]
      if (is.na(idc)) { cat("  ", nm, " 有 ", hit[1], " 但无样本 id 列，跳过\n", sep = ""); next }
      x <- x[!duplicated(x[[idc]]), ]; v <- x[[hit[1]]][match(gs$sample_id, x[[idc]])]
      cat("  ", pat, " → ", nm, "$", hit[1], "：\n", sep = ""); print(table(v, useNA = "ifany")); return(v) } }; NULL }
  for (g in c("kras", "stk11", "tp53", "smok")) { v <- pick(g)
    if (!is.null(v)) { v <- tolower(trimws(as.character(v))); v[v %in% c("", "na", "n/a", "unknown", "not available")] <- NA
      if (length(unique(na.omit(v))) >= 2 && sum(!is.na(v)) >= 200) { gs[[paste0("cv_", g)]] <- factor(v) } else cat("  ", g, "：非缺失 <200 或只有一个水平，不入模型\n") } }
  extra <- grep("^cv_", names(gs), value = TRUE)
  if (length(extra)) { for (x in extra) add("GSE72094", paste("+", x), gs, c("age", "stage_num", x))
    add("GSE72094", "+ all of: sex, mutation/smoking found", gs, c("age", "stage_num", "sex", extra)) }
  Dd <- do.call(rbind, rows); print(format(Dd, digits = 3), row.names = FALSE); wout(Dd, "reviewer_D_confounding.csv") })

# ---------------------------------------------------------------- E
section(paste0("E. 随机基因集基准（M1；Venet 2011），", N_RAND, " 组随机 50 基因"), {
  stopifnot(!is.null(EM), !is.null(FULL))
  uni <- setdiff(intersect(rownames(EM), rownames(FULL)), c(sig, "TIMP1", "PDPN"))
  cat("随机基因池（两平台共有、去掉签名与 TIMP1/PDPN）：", length(uni), " 个基因\n", sep = "")
  sets <- lapply(seq_len(N_RAND), function(i) sample(uni, 50)); names(sets) <- paste0("R", seq_len(N_RAND))
  res <- list()
  for (cn in c("TCGA", "GSE72094")) { d <- COH[[cn]]$d
    mat <- if (cn == "TCGA") EM[, tc$col] else FULL[, gs$sample_id]
    sc <- ssg(mat, sets); stopifnot(nrow(sc) == N_RAND)
    rr <- do.call(rbind, lapply(seq_len(N_RAND), function(i) { r <- as.numeric(sc[i, ])
      a <- cox1(d$t, d$e, r); b <- cox1(d$t, d$e, zs(r) + zs(d$TIMP1) - zs(d$PDPN))
      data.frame(cohort = cn, set = rownames(sc)[i], rand_alone_HR = a[["HR"]], rand_alone_p = a[["p"]], rand_SCANI_HR = b[["HR"]], rand_SCANI_p = b[["p"]]) }))
    obs <- cox1(d$t, d$e, d$SCANI); obm <- cox1(d$t, d$e, d$myCAF)
    cat(sprintf("\n%s：真实 myCAF 单独 HR %.3f (P=%.3g)；随机 50 基因单独 P<0.05 的比例 %.1f%%（HR>1 且 P<0.05：%.1f%%）\n", cn, obm[["HR"]], obm[["p"]],
                100 * mean(rr$rand_alone_p < 0.05), 100 * mean(rr$rand_alone_p < 0.05 & rr$rand_alone_HR > 1)))
    cat(sprintf("%s：真实 SCANI HR %.3f；用随机基因集替换 myCAF 后的 HR 中位数 %.3f（2.5–97.5%%：%.3f–%.3f）；真实值所处百分位 %.1f%%；随机版 HR ≥ 真实值的比例 %.3f\n",
                cn, obs[["HR"]], median(rr$rand_SCANI_HR), quantile(rr$rand_SCANI_HR, .025), quantile(rr$rand_SCANI_HR, .975),
                100 * mean(rr$rand_SCANI_HR < obs[["HR"]]), mean(rr$rand_SCANI_HR >= obs[["HR"]])))
    res[[cn]] <- rr }
  R <- do.call(rbind, res); wout(R, "reviewer_E_random_draws.csv")
  S <- do.call(rbind, lapply(names(res), function(cn) { rr <- res[[cn]]; d <- COH[[cn]]$d; obs <- cox1(d$t, d$e, d$SCANI)
    data.frame(cohort = cn, n_sets = nrow(rr), obs_SCANI_HR = obs[["HR"]], rand_SCANI_HR_median = median(rr$rand_SCANI_HR),
               rand_SCANI_HR_q025 = quantile(rr$rand_SCANI_HR, .025), rand_SCANI_HR_q975 = quantile(rr$rand_SCANI_HR, .975),
               frac_rand_SCANI_HR_ge_obs = mean(rr$rand_SCANI_HR >= obs[["HR"]]),
               frac_rand_alone_p_lt_05 = mean(rr$rand_alone_p < 0.05), frac_rand_alone_p_lt_05_HRgt1 = mean(rr$rand_alone_p < 0.05 & rr$rand_alone_HR > 1)) }))
  print(format(S, digits = 3), row.names = FALSE); wout(S, "reviewer_E_random_summary.csv") })

# ---------------------------------------------------------------- F
section("F. 函数形态（M3）：自然样条 df=3 vs 线性（只看形态；不做分时段 PH 分析）", {
  Fr <- do.call(rbind, lapply(names(COH), function(cn) { d <- COH[[cn]]$d; d$z <- zs(d$SCANI)
    m0 <- coxph(Surv(t, e) ~ z, d); m1 <- coxph(Surv(t, e) ~ ns(z, df = 3), d); l <- anova(m0, m1)
    grid <- seq(quantile(d$z, .05), quantile(d$z, .95), length.out = 7)
    lp <- predict(m1, newdata = data.frame(z = grid), type = "lp") - predict(m1, newdata = data.frame(z = 0), type = "lp")
    cat(cn, "：样条 vs 线性 LR P = ", signif(l[[ncol(l)]][2], 3), "；相对 z=0 的 HR（5%–95% 分位范围）：\n", sep = "")
    print(data.frame(z = round(grid, 2), HR_vs_z0 = round(exp(lp), 3)), row.names = FALSE)
    data.frame(cohort = cn, z = grid, HR_vs_z0 = exp(lp), nonlinearity_LR_p = l[[ncol(l)]][2]) }))
  wout(Fr, "reviewer_F_spline.csv") })

# ---------------------------------------------------------------- G
section("G. 3 年时仍在风险中的人数；TCGA-CDR 生存与分析表一致性（minor）", {
  G <- do.call(rbind, lapply(names(COH), function(cn) { d <- COH[[cn]]$d; hi <- d$SCANI >= median(d$SCANI)
    data.frame(cohort = cn, group = c("all", "high", "low"), at_risk_3y = c(sum(d$t >= D3), sum(d$t[hi] >= D3), sum(d$t[!hi] >= D3)),
               events_by_3y = c(sum(d$e == 1 & d$t <= D3), sum(d$e[hi] == 1 & d$t[hi] <= D3), sum(d$e[!hi] == 1 & d$t[!hi] <= D3))) }))
  print(G, row.names = FALSE); wout(G, "reviewer_G_at_risk_3y.csv")
  if (!is.null(CDR)) { m <- match(tc$patient_id, CDR$bcr_patient_barcode)
    dt <- tc$t - CDR$OS.time[m]; de <- tc$e != CDR$OS[m]
    cat(sprintf("166 例：CDR 生存时间缺失 %d；时间不同 %d 例（|差|中位 %.0f 天）；事件状态不同 %d 例\n", sum(is.na(dt)), sum(abs(dt) > 0, na.rm = TRUE),
                median(abs(dt[abs(dt) > 0]), na.rm = TRUE), sum(de, na.rm = TRUE)))
    d2 <- tc[!is.na(m), ]; d2$t <- CDR$OS.time[m][!is.na(m)]; d2$e <- CDR$OS[m][!is.na(m)]; d2 <- d2[!is.na(d2$t) & d2$t > 0, ]
    u <- cox1(d2$t, d2$e, d2$SCANI)
    cat(sprintf("用 CDR 生存重算 TCGA（n=%d/%d）：SCANI HR %.3f (%.3f–%.3f) P=%.3g\n", nrow(d2), sum(d2$e), u[["HR"]], u[["lo"]], u[["hi"]], u[["p"]]))
    wout(data.frame(analysis = "TCGA SCANI with CDR OS", n = nrow(d2), events = sum(d2$e), HR = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], p = u[["p"]],
                    n_time_differs = sum(abs(dt) > 0, na.rm = TRUE), n_event_differs = sum(de, na.rm = TRUE)), "reviewer_G_CDR_OS_check.csv")
  } else cat("C 节未读到 CDR，跳过一致性检查\n") })

cat("\n全部完成。输出目录：", OD, "\n把 ", LOG, " 全文贴回。\n", sep = "")
