#!/usr/bin/env Rscript
# ==============================================================================
# sci1_cdr226_main.R —— TCGA 发现队列改为 TCGA-CDR 口径（II–IV 期、OS>0、有原发肿瘤 RNA-seq）：n=226
#   原因：原分析表（step6b 由 GDC 临床字段生成）漏掉 60 例——死亡病人缺分期、存活病人缺随访，缺失与生死状态相关。
#   本脚本按稿件写明的纳入标准，用 TCGA-CDR（Liu 2018）的分期、OS、年龄、性别重建队列，重算全部 TCGA 数字。
#   SCANI 公式不变（组分不重选）；每个组分在 226 例内重新 Z 标准化。
#   同 aliquot 规则：每例病人三个组分取同一列——pid_to_col 中已有的病人沿用其列（与 v2 相同）；
#   其余病人取 expr_mat 中该病人第一个原发肿瘤（样本类型 01）列。脚本会报告这两种规则在 166 例中是否一致。
#
#   输出目录 tables/cdr226/（已存在且非空即停止，不覆盖）：
#     SCANI_TCGA_CDR226.csv                主分析表（226 行）
#     headline_numbers_TCGA_CDR226.csv     头条数字（单变量、校正、AUC、KM、随访、PH、C-index、相关）
#     supp_export_cdr226.txt               ### S1a（16 候选基因）与 ### S1b（SCANI vs 组分）两段，格式同 supp_export_for_S1_S2.txt
#     candidate16_univariable_TCGA_CDR226.csv、single_gene_vs_SCANI_TCGA_CDR226.csv
#     remark_S4_characteristics_TCGA_CDR226.csv、remark_S4_models_TCGA_CDR226.csv（只含 TCGA 行，格式同 remark_S4_*_v2.csv）
#     cohort_flow_TCGA_CDR226.csv          患者流程各步人数（Figure 1 / REMARK 12）
#     sensitivity_166_vs_226.csv           原 166 例、226 例、226 中的 166 例、60 例新增病人
#     reviewer_TCGA_CDR226.csv             审稿预演 A/B/D/F/G 各项在 226 例上的结果
#     reviewer_E_random_TCGA_CDR226.csv    随机 50 基因集基准（TCGA 226）
#   图：figures/Figure4_cdr226.{tiff,png,pdf}（不覆盖已有文件）
#   日志：/tmp/sci1_cdr226_<时间>.txt
# 用法：Rscript ~/sci1_cdr226_main.R        （随机基因集默认 1000 组，约 5–10 分钟；N_RAND=50 可试跑，输出到 cdr226_test_N50/）
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(splines); library(ggplot2); library(cowplot) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
B <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS"); T <- file.path(B, "tables")
N_RAND <- as.integer(Sys.getenv("N_RAND", "1000")); NCORE <- 3L; set.seed(20261106)
OD <- file.path(T, if (N_RAND == 1000) "cdr226" else paste0("cdr226_test_N", N_RAND))
FIGN <- if (N_RAND == 1000) "Figure4_cdr226" else paste0("Figure4_cdr226_test_N", N_RAND)
if (dir.exists(OD) && length(list.files(OD))) stop("已存在且非空，不覆盖：", OD)
for (ext in c(".tiff", ".png", ".pdf")) if (file.exists(file.path(B, "figures", paste0(FIGN, ext)))) stop("图已存在，不覆盖：", FIGN, ext)
dir.create(OD, showWarnings = FALSE)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_cdr226_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74), "\n", sep = "")
zs <- function(x) as.numeric(scale(x))
wout <- function(x, name) { f <- file.path(OD, name); if (file.exists(f)) stop("已存在，不覆盖：", f)
  write.csv(x, f, row.names = FALSE); cat("  → 写出 ", f, "\n", sep = ""); invisible(f) }
section <- function(name, expr) { hr(name); t0 <- Sys.time()
  r <- tryCatch({ force(expr); "OK" }, error = function(e) { cat("★ 本节出错，已跳过：", conditionMessage(e), "\n"); "ERROR" })
  cat(sprintf("[%s] %s，用时 %.1f 分钟\n", name, r, as.numeric(difftime(Sys.time(), t0, units = "mins")))); invisible(r) }
cat("SCI-1 CDR 226 主分析  ", format(Sys.time()), "\n", R.version.string, "\n", sep = "")
stopifnot(requireNamespace("GSVA", quietly = TRUE), packageVersion("GSVA") == "1.48.3", requireNamespace("timeROC", quietly = TRUE))
D3 <- 3 * 365
ssg <- function(mat, sets, par = 1L) GSVA::gsva(mat, sets, method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE, parallel.sz = par)
cox1 <- function(t, e, x) { s <- summary(coxph(Surv(t, e) ~ zs(x)))
  c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
coxadj <- function(dd, x, covs) { dd$x. <- zs(x); s <- summary(coxph(as.formula(paste("Surv(t, e) ~ x. +", paste(covs, collapse = " + "))), data = dd))
  c(HR = s$conf.int["x.", 1], lo = s$conf.int["x.", 3], hi = s$conf.int["x.", 4], p = s$coefficients["x.", 5], n = s$n, ev = s$nevent) }
auc3 <- function(t, e, x) unname(timeROC::timeROC(T = t, delta = e, marker = x, cause = 1, times = D3, iid = FALSE)$AUC[2])
cidx <- function(t, e, x) { cc <- concordance(Surv(t, e) ~ x, reverse = TRUE); c(C = unname(cc$concordance), C_se = sqrt(unname(cc$var))) }
kmp <- function(t, e, x) { g <- x >= median(x); 1 - pchisq(survdiff(Surv(t, e) ~ g)$chisq, 1) }
revkm <- function(t, e) { s <- summary(survfit(Surv(t, 1 - e) ~ 1))$table; unname(c(s["median"], s["0.95LCL"], s["0.95UCL"])) }
fmt <- function(m) { s <- summary(m); data.frame(term = rownames(s$conf.int), HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
                                                 p = s$coefficients[, ncol(s$coefficients)], row.names = NULL) }
sandbox <- function() { e <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
  for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "tiff", "jpeg", "ggsave", "dev.off", "fwrite", "write_csv", "sink"))
    assign(f, noop, envir = e); assign("dir.create", function(...) invisible(TRUE), envir = e); e }
stage_roman <- function(s) { s <- toupper(gsub("STAGE|\\s", "", toupper(as.character(s)))); s <- sub("[ABC][0-9]?$", "", s)
  ifelse(s %in% c("I", "II", "III", "IV"), s, NA) }
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)
CAND <- c("ACTA2", "COL1A1", "COL1A2", "MMP2", "TAGLN", "FAP", "POSTN", "PDPN", "S100A4", "THBS2", "FN1", "VIM", "SPARC", "TIMP1", "CXCL12", "TGFB1")

# =============================================================== 1. 建队列
hr("1. 用 TCGA-CDR 建队列，逐步计数")
o166 <- read.csv(file.path(T, "SCANI_TCGA_O_166_v2.csv"), stringsAsFactors = FALSE); o166 <- o166[!duplicated(o166$patient_id), ]
stopifnot(nrow(o166) == 166, sum(o166$OS_event) == 80)
stopifnot(abs(cox1(o166$OS_time, o166$OS_event, o166$SCANI)[["HR"]] - 1.4629) < 1e-3)
fcdr <- file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/raw/TCGA_CDR/TCGA-CDR-SupplementalTableS1.xlsx")
stopifnot(unname(tools::md5sum(fcdr)) == "a4591b2dcee39591f59e5e25a6ce75fa")
cdr <- as.data.frame(readxl::read_excel(fcdr, sheet = 1, guess_max = 20000)); cdr <- cdr[cdr$type %in% "LUAD", ]
stopifnot(!anyDuplicated(cdr$bcr_patient_barcode))
cdr$t <- suppressWarnings(as.numeric(cdr$OS.time)); cdr$e <- suppressWarnings(as.numeric(cdr$OS))
cdr$age <- suppressWarnings(as.numeric(cdr$age_at_initial_pathologic_diagnosis)); cdr$sex <- tolower(cdr$gender)
cdr$stage4 <- stage_roman(cdr$ajcc_pathologic_tumor_stage)
e6 <- sandbox(); eval(parse(text = readLines(file.path(C, "scripts/step6d_genes.R"), warn = FALSE)[1:17]), envir = e6)
EM <- get("expr_mat", e6); P2C <- get("pid_to_col", e6)
p2n <- function(v) if (is.numeric(v)) colnames(EM)[v] else as.character(v)
first01 <- function(p) { cs <- colnames(EM)[substr(colnames(EM), 1, 12) == p & substr(colnames(EM), 14, 15) == "01"]; if (length(cs)) cs[1] else NA_character_ }
cdr$col_p2c <- vapply(cdr$bcr_patient_barcode, function(p) if (p %in% names(P2C)) p2n(P2C[[p]])[1] else NA_character_, "")
cdr$col_01 <- vapply(cdr$bcr_patient_barcode, first01, "")
cdr$col <- ifelse(!is.na(cdr$col_p2c), cdr$col_p2c, cdr$col_01)
k166 <- cdr$bcr_patient_barcode %in% o166$patient_id
cat("166 例中：pid_to_col 的列 == 第一个原发肿瘤列：", sum(cdr$col_p2c[k166] == cdr$col_01[k166], na.rm = TRUE), "/166\n", sep = "")
FL <- list(); step <- function(lbl, keep) { FL[[length(FL) + 1]] <<- data.frame(step = lbl, n = sum(keep), deaths = sum(cdr$e[keep] %in% 1)); keep }
k <- step("TCGA-CDR LUAD patients", rep(TRUE, nrow(cdr)))
k <- step("with primary-tumour RNA-seq", k & !is.na(cdr$col) & substr(cdr$col, 14, 15) == "01")
k <- step("with OS time > 0 and OS status", k & !is.na(cdr$t) & cdr$t > 0 & !is.na(cdr$e)); k3 <- k
invisible(step("  excluded: stage unknown", k3 & is.na(cdr$stage4))); invisible(step("  excluded: stage I", k3 & cdr$stage4 %in% "I"))
k <- step("stage II–IV (analysis cohort)", k3 & cdr$stage4 %in% c("II", "III", "IV"))
FLd <- do.call(rbind, FL)
print(FLd, row.names = FALSE); wout(FLd, "cohort_flow_TCGA_CDR226.csv")
d <- data.frame(patient_id = cdr$bcr_patient_barcode[k], aliquot = cdr$col[k], OS_time = cdr$t[k], OS_event = cdr$e[k],
                age = cdr$age[k], sex = cdr$sex[k], stage = cdr$stage4[k], in_original_166 = k166[k], stringsAsFactors = FALSE)
stopifnot(nrow(d) == 226, sum(d$OS_event) == 115, !anyDuplicated(d$patient_id), all(o166$patient_id %in% d$patient_id))
d$myCAF <- as.numeric(ssg(EM[, d$aliquot], list(myCAF = sig))["myCAF", ])
d$TIMP1 <- as.numeric(EM["TIMP1", d$aliquot]); d$PDPN <- as.numeric(EM["PDPN", d$aliquot])
m <- match(o166$patient_id, d$patient_id)
stopifnot(all(d$aliquot[m] == o166$aliquot), max(abs(d$TIMP1[m] - o166$TIMP1)) < 1e-8, max(abs(d$PDPN[m] - o166$PDPN)) < 1e-8)
cat("166 例 myCAF 与 v2 表最大差：", signif(max(abs(d$myCAF[m] - o166$myCAF)), 3), "（ssGSEA 每样本独立，应 <1e-6）\n", sep = "")
stopifnot(max(abs(d$myCAF[m] - o166$myCAF)) < 1e-6)
cat("166 例 OS 与 CDR 一致：时间不同 ", sum(d$OS_time[m] != o166$OS_time), "，事件不同 ", sum(d$OS_event[m] != o166$OS_event),
    "；年龄（CDR vs 原表）差异 >0 的 ", sum(abs(d$age[m] - o166$age_at_index) > 0, na.rm = TRUE), " 例\n", sep = "")
d$SCANI <- zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)
d$stage <- factor(d$stage, levels = c("II", "III", "IV")); d$t <- d$OS_time; d$e <- d$OS_event
cat("年龄缺失：", sum(is.na(d$age)), "；性别：\n", sep = ""); print(table(d$sex, useNA = "ifany")); print(table(d$stage))
wout(d[, c("patient_id", "aliquot", "myCAF", "TIMP1", "PDPN", "SCANI", "OS_time", "OS_event", "stage", "age", "sex", "in_original_166")], "SCANI_TCGA_CDR226.csv")
cc <- complete.cases(d[, c("age", "stage")])

# =============================================================== 2. 头条数字
section("2. 头条数字（226）", {
  u <- cox1(d$t, d$e, d$SCANI); mv <- coxph(Surv(t, e) ~ zs(SCANI) + age + stage, d[cc, ]); s <- summary(mv)
  ci <- cidx(d$t, d$e, d$SCANI); fu <- revkm(d$t, d$e); z1 <- cox.zph(coxph(Surv(t, e) ~ zs(SCANI), d)); z2 <- cox.zph(mv)
  ra <- cor.test(d$SCANI, d$age, method = "spearman", exact = FALSE); kw <- kruskal.test(SCANI ~ stage, d); ws <- wilcox.test(SCANI ~ sex, d)
  rmp <- cor.test(d$myCAF, d$PDPN, method = "spearman", exact = FALSE)
  Hd <- data.frame(metric = c("n", "deaths", "uni_HR", "uni_lo", "uni_hi", "uni_p", "mv_n", "mv_deaths", "mv_HR", "mv_lo", "mv_hi", "mv_p",
                              "AUC3y", "C", "C_se", "KM_median_split_p", "followup_months", "followup_lo", "followup_hi",
                              "PH_uni_p", "PH_mv_SCANI_p", "PH_mv_global_p", "rho_age", "p_age", "KW_stage_p", "Wilcoxon_sex_p", "rho_myCAF_PDPN", "at_risk_3y"),
                   value = c(nrow(d), sum(d$e), u, s$n, s$nevent, s$conf.int[1, 1], s$conf.int[1, 3], s$conf.int[1, 4], s$coefficients[1, 5],
                             auc3(d$t, d$e, d$SCANI), ci, kmp(d$t, d$e, d$SCANI), fu / 30.44,
                             z1$table[1, "p"], z2$table[1, "p"], z2$table["GLOBAL", "p"], ra$estimate, ra$p.value, kw$p.value, ws$p.value, rmp$estimate, sum(d$t >= D3)))
  print(format(Hd, digits = 4), row.names = FALSE); wout(Hd, "headline_numbers_TCGA_CDR226.csv") })

# =============================================================== 3. Table S2
section("3. Table S2：16 个候选基因（A）与 SCANI vs 组分（B）", {
  miss <- setdiff(CAND, rownames(EM)); if (length(miss)) stop("expr_mat 缺候选基因：", paste(miss, collapse = ","))
  A <- do.call(rbind, lapply(CAND, function(g) { x <- as.numeric(EM[g, d$aliquot]); u <- cox1(d$t, d$e, x)
    data.frame(gene = g, HR_perSD = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], Cox_p = u[["p"]], KM_p = kmp(d$t, d$e, x)) }))
  A$Cox_p_bonf16 <- pmin(1, A$Cox_p * 16); A <- A[order(A$Cox_p), ]
  Bt <- do.call(rbind, lapply(c("SCANI", "TIMP1", "PDPN", "myCAF"), function(nm) { u <- cox1(d$t, d$e, d[[nm]])
    data.frame(score = nm, HR_perSD = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], Cox_p = u[["p"]], AUC3y = auc3(d$t, d$e, d[[nm]])) }))
  print(format(A, digits = 3), row.names = FALSE); print(format(Bt, digits = 3), row.names = FALSE)
  fa <- wout(A, "candidate16_univariable_TCGA_CDR226.csv"); fb <- wout(Bt, "single_gene_vs_SCANI_TCGA_CDR226.csv")
  txt <- c("### S1a candidate genes (n=226; TCGA-CDR)", paste0("# file: ", fa, "  md5 ", unname(tools::md5sum(fa))), readLines(fa), "",
           "### S1b SCANI vs components (n=226; TCGA-CDR)", paste0("# file: ", fb, "  md5 ", unname(tools::md5sum(fb))), readLines(fb), "")
  f <- file.path(OD, "supp_export_cdr226.txt"); writeLines(txt, f); cat("  → 写出 ", f, "\n", sep = "") })

# =============================================================== 4. Table S4（TCGA 行）
section("4. Table S4：队列特征、SCANI 与临床变量、多变量全部协变量、PH 检验（TCGA 226）", {
  CH <- list(); addch <- function(var, level, value) CH[[length(CH) + 1]] <<- data.frame(cohort = "TCGA", var = var, level = level, value = as.character(value))
  fu <- revkm(d$t, d$e)
  addch("patients/deaths", "", sprintf("%d/%d", nrow(d), sum(d$e)))
  addch("median follow-up, months (reverse KM, 95% CI)", "", sprintf("%.1f (%.1f–%.1f)", fu[1] / 30.44, fu[2] / 30.44, fu[3] / 30.44))
  addch("age, years", "median (IQR); missing", sprintf("%.0f (%.0f–%.0f); %d", median(d$age, na.rm = TRUE), quantile(d$age, .25, na.rm = TRUE),
        quantile(d$age, .75, na.rm = TRUE), sum(is.na(d$age))))
  for (s in c(names(table(d$sex)), NA)) addch("sex", ifelse(is.na(s), "missing", s), sum(if (is.na(s)) is.na(d$sex) else d$sex %in% s))
  for (s in levels(d$stage)) addch("stage", s, sum(d$stage %in% s))
  ws <- wilcox.test(SCANI ~ sex, d); ra <- cor.test(d$SCANI, d$age, method = "spearman", exact = FALSE); kw <- kruskal.test(SCANI ~ stage, d)
  addch("SCANI vs sex", "Wilcoxon P", signif(ws$p.value, 3)); addch("SCANI vs age", "Spearman rho (P)", sprintf("%.2f (%.2g)", ra$estimate, ra$p.value))
  addch("SCANI vs stage", "Kruskal–Wallis P", signif(kw$p.value, 3))
  CHd <- do.call(rbind, CH); print(CHd, row.names = FALSE)
  dd <- d[cc, ]; mv <- coxph(Surv(t, e) ~ zs(SCANI) + age + stage, dd); mu <- coxph(Surv(t, e) ~ zs(SCANI), d)
  z1 <- cox.zph(mu); z2 <- cox.zph(mv); print(fmt(mv), digits = 3); print(z1$table); print(z2$table)
  rn <- rownames(z2$table); rn[rn == "age"] <- "age_at_index"   # 与 v2 文件的标签一致，方便构建脚本复用
  MO <- rbind(cbind(cohort = "TCGA", model = "multivariable", within(fmt(mv), term[term == "age"] <- "age_at_index")),
              cbind(cohort = "TCGA", model = "univariable", fmt(mu)),
              data.frame(cohort = "TCGA", model = c("PH univariable SCANI", paste0("PH multivariable ", rn)), term = "zph", HR = NA, lo = NA, hi = NA,
                         p = c(z1$table[1, "p"], z2$table[, "p"])))
  cat("多变量完整病例：", nrow(dd), " 例，", sum(dd$e), " 例死亡\n", sep = "")
  wout(CHd, "remark_S4_characteristics_TCGA_CDR226.csv"); wout(MO, "remark_S4_models_TCGA_CDR226.csv") })

# =============================================================== 5. 敏感性：166 vs 226
section("5. 敏感性分析：原 166 例 / 226 例 / 226 中的 166 例 / 60 例新增", {
  rows <- list(); add <- function(lbl, dd, x, adj = TRUE) { u <- cox1(dd$t, dd$e, x); c2 <- complete.cases(dd[, c("age", "stage")])
    a <- if (adj) coxadj(dd[c2, ], x[c2], c("age", "stage")) else c(HR = NA, lo = NA, hi = NA, p = NA, n = NA, ev = NA)
    rows[[length(rows) + 1]] <<- data.frame(analysis = lbl, n = nrow(dd), deaths = sum(dd$e), HR = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], p = u[["p"]],
      adj_n = a[["n"]], adj_HR = a[["HR"]], adj_lo = a[["lo"]], adj_hi = a[["hi"]], adj_p = a[["p"]], AUC3y = auc3(dd$t, dd$e, x), KM_p = kmp(dd$t, dd$e, x)) }
  oo <- o166; oo$t <- oo$OS_time; oo$e <- oo$OS_event; oo$age <- oo$age_at_index; oo$stage <- factor(oo$stage_simple, levels = c("II", "III", "IV"))
  add("Original analysis table (GDC clinical; n=166), SCANI standardized in 166", oo, oo$SCANI)
  add("TCGA-CDR cohort (n=226), SCANI standardized in 226 [primary]", d, d$SCANI)
  s1 <- d[d$in_original_166, ]; add("Original 166 patients, SCANI standardized in 226, CDR clinical data", s1, s1$SCANI)
  s2 <- d[!d$in_original_166, ]; add("60 patients not in the original table, SCANI standardized in 226", s2, s2$SCANI)
  S <- do.call(rbind, rows); print(format(S, digits = 3), row.names = FALSE); wout(S, "sensitivity_166_vs_226.csv")
  cat("60 例新增病人的分期：\n"); print(table(s2$stage)); cat("死亡：", sum(s2$e), "/", nrow(s2), "\n") })

# =============================================================== 6. 审稿预演各项（226）
section("6. 审稿预演 A/B/D/F/G 在 226 例上", {
  R <- list(); put <- function(item, lbl, v) R[[length(R) + 1]] <<- data.frame(item = item, analysis = lbl, t(v))
  VAR <- list(SCANI = c(myCAF = 1, TIMP1 = 1, PDPN = -1), TIMP1_only = c(TIMP1 = 1), myCAF_only = c(myCAF = 1), PDPN_only_neg = c(PDPN = -1),
              myCAF_plus_TIMP1 = c(myCAF = 1, TIMP1 = 1), TIMP1_minus_PDPN = c(TIMP1 = 1, PDPN = -1), myCAF_minus_PDPN = c(myCAF = 1, PDPN = -1),
              all_three_PDPN_positive = c(myCAF = 1, TIMP1 = 1, PDPN = 1))
  score <- function(w) Reduce(`+`, lapply(names(w), function(k) w[[k]] * zs(d[[k]])))
  for (v in names(VAR)) { x <- score(VAR[[v]]); u <- cox1(d$t, d$e, x); a <- coxadj(d[cc, ], x[cc], c("age", "stage")); ci <- cidx(d$t, d$e, x)
    put("A_variant", v, c(HR = u[["HR"]], lo = u[["lo"]], hi = u[["hi"]], p = u[["p"]], adj_HR = a[["HR"]], adj_lo = a[["lo"]], adj_hi = a[["hi"]], adj_p = a[["p"]],
                          C = ci[["C"]], AUC3y = auc3(d$t, d$e, x))) }
  d$zT <- zs(d$TIMP1); d$zMP <- zs(d$myCAF) - zs(d$PDPN); d$zM <- zs(d$myCAF); d$zP <- zs(d$PDPN)
  m0 <- coxph(Surv(t, e) ~ zT, d); m1 <- coxph(Surv(t, e) ~ zT + zMP, d); m2 <- coxph(Surv(t, e) ~ zT + zM + zP, d)
  l1 <- anova(m0, m1); l2 <- anova(m0, m2)
  put("B_LRT", "TIMP1 vs TIMP1 + [Z(myCAF)-Z(PDPN)] (1 df)", c(chisq = l1$Chisq[2], p = l1[[ncol(l1)]][2], C_TIMP1 = unname(concordance(m0)$concordance), C_full = unname(concordance(m1)$concordance)))
  put("B_LRT", "TIMP1 vs TIMP1 + myCAF + PDPN (2 df)", c(chisq = l2$Chisq[2], p = l2[[ncol(l2)]][2], C_TIMP1 = unname(concordance(m0)$concordance), C_full = unname(concordance(m2)$concordance)))
  a <- coxadj(d[cc, ], d$SCANI[cc], c("age", "stage", "sex")); put("D_confounding", "+ sex", c(HR = a[["HR"]], lo = a[["lo"]], hi = a[["hi"]], p = a[["p"]], n = a[["n"]]))
  if (requireNamespace("tidyestimate", quietly = TRUE)) tryCatch({
    f <- tidyestimate::filter_common_genes(as.data.frame(EM[, d$aliquot]), id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE)
    es <- tidyestimate::estimate_score(f, is_affymetrix = FALSE); gc <- function(p) es[[grep(p, names(es), ignore.case = TRUE, value = TRUE)[1]]]
    smp <- as.character(gc("^sample")); d$zstromal <- zs(gc("^stromal")[match(d$aliquot, smp)]); d$zimmune <- zs(gc("^immune")[match(d$aliquot, smp)])
    stopifnot(!anyNA(d$zstromal))
    cat(sprintf("SCANI 与 ESTIMATE 基质 Spearman %.3f；免疫 %.3f\n", cor(d$SCANI, d$zstromal, method = "spearman"), cor(d$SCANI, d$zimmune, method = "spearman")))
    a <- coxadj(d[cc, ], d$SCANI[cc], c("age", "stage", "zstromal", "zimmune"))
    put("D_confounding", "+ ESTIMATE stromal + immune", c(HR = a[["HR"]], lo = a[["lo"]], hi = a[["hi"]], p = a[["p"]], n = a[["n"]],
        rho_stromal = cor(d$SCANI, d$zstromal, method = "spearman"), rho_immune = cor(d$SCANI, d$zimmune, method = "spearman")))
  }, error = function(e) cat("★ ESTIMATE 出错，跳过：", conditionMessage(e), "\n")) else cat("★ 未安装 tidyestimate，跳过 ESTIMATE\n")
  d$z <- zs(d$SCANI); f0 <- coxph(Surv(t, e) ~ z, d); f1 <- coxph(Surv(t, e) ~ ns(z, df = 3), d); l <- anova(f0, f1)
  grid <- seq(quantile(d$z, .05), quantile(d$z, .95), length.out = 7)
  lp <- predict(f1, newdata = data.frame(z = grid), type = "lp") - predict(f1, newdata = data.frame(z = 0), type = "lp")
  put("F_spline", "nonlinearity LR test", c(p = l[[ncol(l)]][2]))
  for (i in seq_along(grid)) put("F_spline", sprintf("HR at z=%.2f vs z=0", grid[i]), c(HR = unname(exp(lp[i]))))
  hi <- d$SCANI >= median(d$SCANI)
  put("G_at_risk_3y", "all / high / low", c(all = sum(d$t >= D3), high = sum(d$t[hi] >= D3), low = sum(d$t[!hi] >= D3)))
  COLS <- c("item", "analysis", "HR", "lo", "hi", "p", "adj_HR", "adj_lo", "adj_hi", "adj_p", "C", "AUC3y", "chisq", "C_TIMP1", "C_full", "n",
            "rho_stromal", "rho_immune", "all", "high", "low")
  RR <- do.call(rbind, lapply(R, function(x) { bad <- setdiff(names(x), COLS); if (length(bad)) stop("未预期的列：", paste(bad, collapse = ","))
    x[setdiff(COLS, names(x))] <- NA; x[, COLS] }))
  print(format(RR, digits = 3), row.names = FALSE); wout(RR, "reviewer_TCGA_CDR226.csv") })

# =============================================================== 7. 随机基因集（226）
section(paste0("7. 随机 50 基因集基准（TCGA 226；", N_RAND, " 组）"), {
  uni <- setdiff(rownames(EM), c(sig, "TIMP1", "PDPN")); fg <- file.path(T, "reviewer_prep/reviewer_E_random_draws.csv")
  cat("随机基因池：", length(uni), " 个（TCGA 全部基因去掉签名与 TIMP1/PDPN；与 reviewer_prep 的两平台交集池不同）\n", sep = "")
  sets <- lapply(seq_len(N_RAND), function(i) sample(uni, 50)); names(sets) <- paste0("R", seq_len(N_RAND))
  sc <- ssg(EM[, d$aliquot], sets, NCORE)
  rr <- do.call(rbind, lapply(seq_len(N_RAND), function(i) { r <- as.numeric(sc[i, ]); a <- cox1(d$t, d$e, r); b <- cox1(d$t, d$e, zs(r) + zs(d$TIMP1) - zs(d$PDPN))
    data.frame(set = rownames(sc)[i], rand_alone_HR = a[["HR"]], rand_alone_p = a[["p"]], rand_SCANI_HR = b[["HR"]], rand_SCANI_p = b[["p"]]) }))
  obs <- cox1(d$t, d$e, d$SCANI); obm <- cox1(d$t, d$e, d$myCAF)
  cat(sprintf("真实 SCANI HR %.3f；随机替换 myCAF 后 HR 中位数 %.3f（2.5–97.5%%：%.3f–%.3f）；随机版 HR ≥ 真实值比例 %.3f\n", obs[["HR"]],
              median(rr$rand_SCANI_HR), quantile(rr$rand_SCANI_HR, .025), quantile(rr$rand_SCANI_HR, .975), mean(rr$rand_SCANI_HR >= obs[["HR"]])))
  cat(sprintf("myCAF 单独 HR %.3f；随机集单独 HR ≥ 该值比例 %.3f；随机集单独 P<0.05 比例 %.3f\n", obm[["HR"]], mean(rr$rand_alone_HR >= obm[["HR"]]), mean(rr$rand_alone_p < 0.05)))
  wout(rbind(data.frame(set = "SUMMARY", rand_alone_HR = obm[["HR"]], rand_alone_p = mean(rr$rand_alone_HR >= obm[["HR"]]),
                        rand_SCANI_HR = obs[["HR"]], rand_SCANI_p = mean(rr$rand_SCANI_HR >= obs[["HR"]])), rr), "reviewer_E_random_TCGA_CDR226.csv")
  cat("（首行 SUMMARY：rand_alone_HR=真实 myCAF HR，rand_alone_p=随机≥真实的比例；rand_SCANI_HR=真实 SCANI HR，rand_SCANI_p=随机≥真实的比例）\n") })

# =============================================================== 8. Figure 4
section("8. Figure 4（226）：KM + number at risk + 3 年 ROC", {
  XMAX <- 120; BRK <- seq(0, XMAX, 24)
  d$grp <- factor(ifelse(d$SCANI >= median(d$SCANI), "High", "Low"), levels = c("Low", "High"))
  km_p <- kmp(d$t, d$e, d$SCANI); fit <- survfit(Surv(OS_time / 30.44, OS_event) ~ grp, d); ns <- table(d$grp)
  labs_g <- setNames(paste0(names(ns), " (n=", ns, ")"), names(ns)); gcol <- setNames(c("#4575B4", "#D73027"), labs_g)
  kmdf <- do.call(rbind, lapply(seq_along(fit$strata), function(i) {
    idx <- if (i == 1) seq_len(fit$strata[1]) else sum(fit$strata[1:(i - 1)]) + seq_len(fit$strata[i]); g <- sub("grp=", "", names(fit$strata)[i])
    data.frame(time = c(0, fit$time[idx]), surv = c(1, fit$surv[idx]), cens = c(0, fit$n.censor[idx]), group = unname(labs_g[g])) }))
  kmdf$group <- factor(kmdf$group, levels = labs_g)
  p4A <- ggplot(kmdf, aes(time, surv, colour = group)) + geom_step(linewidth = 0.8) +
    geom_point(data = subset(kmdf, cens > 0), shape = 3, size = 1.6, stroke = 0.6, show.legend = FALSE) +
    scale_colour_manual(values = gcol, name = "SCANI") + scale_x_continuous(breaks = BRK) + scale_y_continuous(limits = c(0, 1)) +
    coord_cartesian(xlim = c(0, XMAX)) + annotate("text", x = XMAX * 0.7, y = 0.95, size = 3.3, label = sprintf("Log-rank p = %.3f", km_p)) +
    labs(x = NULL, y = "Overall survival probability") + theme_classic(base_size = 11) +
    theme(legend.position = "top", legend.text = element_text(margin = margin(l = 4, r = 14)))
  s <- summary(fit, times = BRK, extend = TRUE)
  rt <- data.frame(time = s$time, n = s$n.risk, group = unname(labs_g[sub("grp=", "", as.character(s$strata))])); rt$group <- factor(rt$group, levels = rev(labs_g))
  cat("Number at risk:\n"); print(xtabs(n ~ group + time, rt))
  pRisk <- ggplot(rt, aes(time, group, label = n, colour = group)) + geom_text(size = 3.2, show.legend = FALSE) + scale_colour_manual(values = gcol) +
    scale_x_continuous(breaks = BRK) + coord_cartesian(xlim = c(0, XMAX)) + scale_y_discrete(labels = function(x) sub(" \\(.*", "", x)) +
    labs(x = "Time (months)", y = NULL, title = "Number at risk") + theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 10), axis.line.y = element_blank(), axis.ticks.y = element_blank())
  km_block <- plot_grid(p4A, pRisk, ncol = 1, rel_heights = c(4, 1.25), align = "v", axis = "lr")
  roc <- timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = d$SCANI, cause = 1, times = D3, iid = FALSE); auc <- unname(roc$AUC[2])
  rdf <- data.frame(FPR = roc$FP[, 2], TPR = roc$TP[, 2]); rdf <- rdf[order(rdf$FPR, rdf$TPR), ]
  p4B <- ggplot(rdf, aes(FPR, TPR)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") + geom_step(colour = "#B2182B", linewidth = 0.8) +
    annotate("text", x = 0.65, y = 0.12, size = 3.3, label = sprintf("3-year AUC = %.3f", auc)) + coord_equal() +
    labs(x = "1 - specificity", y = "Sensitivity") + theme_classic(base_size = 11)
  fig4 <- plot_grid(km_block, p4B, nrow = 1, labels = c("A", "B"), label_size = 14, rel_widths = c(1.2, 1))
  ff <- function(ext) file.path(B, "figures", paste0(FIGN, ext))
  ggsave(ff(".tiff"), fig4, width = 10, height = 5, dpi = 600, compression = "lzw", bg = "white")
  ggsave(ff(".png"), fig4, width = 10, height = 5, dpi = 600, bg = "white"); ggsave(ff(".pdf"), fig4, width = 10, height = 5, bg = "white")
  cat(sprintf("Figure 4：log-rank p = %.3f；3 年 AUC = %.3f；最长随访 %.1f 月；>%d 月仍随访 %d 例\n", km_p, auc, max(d$OS_time) / 30.44, XMAX, sum(d$OS_time / 30.44 > XMAX)))
  cat("  → 写出 ", ff(".png"), "\n", sep = "") })

cat("\n全部完成。输出：", OD, "\n把 ", LOG, " 全文、", OD, "/supp_export_cdr226.txt，以及 ", FIGN, ".png 发回。\n", sep = "")
