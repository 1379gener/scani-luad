#!/usr/bin/env Rscript
# ==============================================================================
# sci1_recheck2.R —— 只读复算原始数据；只写 ~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables/
#   A. 16 个候选基因：照 step6d 读入，按病人去重到 166 后逐个单变量 Cox（per SD）+ KM
#   A2. TIMP1 / PDPN 单基因对照（O 口径 166 例）：Cox per SD + 3 年 AUC，与 SCANI 并列
#   B. GSE72094：step22 实际只抽取了 55 个探针（21 个基因），myCAF ssGSEA 只用 12/50 基因、
#      在 21 基因的宇宙里打分。这里按 Methods 所写的做法（全基因组塌缩、同一 50 基因、同一 GSVA 设置）
#      重算一次，并与原分数并列报告。两个结果都照实输出，不挑。
# 用法：Rscript ~/sci1_recheck2.R
# ==============================================================================
suppressPackageStartupMessages({ library(survival) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTT <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
LOG <- file.path("/tmp", format(Sys.time(), "sci1_recheck2_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
cat("SCI-1 复核 2  ", format(Sys.time()), "\n日志: ", LOG, "\n", sep = "")
cat("GSVA ", as.character(packageVersion("GSVA")), "；survival ", as.character(packageVersion("survival")), "\n", sep = "")

sandbox <- function() {
  e <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
  for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "tiff", "jpeg",
              "ggsave", "dev.off", "fwrite", "write_csv")) assign(f, noop, envir = e)
  assign("dir.create", function(...) invisible(TRUE), envir = e); e
}
run_lines <- function(file, upto, env) eval(parse(text = readLines(file, warn = FALSE)[seq_len(upto)]), envir = env)
coxsd <- function(t, e, x) { s <- summary(coxph(Surv(t, e) ~ as.numeric(scale(x))))
  c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
kmp <- function(t, e, x) { g <- x >= median(x); 1 - pchisq(survdiff(Surv(t, e) ~ g)$chisq, 1) }
auc3 <- function(t, e, x) unname(timeROC::timeROC(T = t, delta = e, marker = x, cause = 1, times = 3 * 365, iid = FALSE)$AUC[2])

# ---------------------------------------------------------------- A
hr("A. 16 个候选基因（step6d 读入，按病人去重）")
e6 <- sandbox(); run_lines(file.path(C, "scripts/step6d_genes.R"), 37, e6)
d  <- get("df2", e6); em <- get("expr_mat", e6); p2c <- get("pid_to_col", e6); kg <- get("key_genes", e6)
cat("step6d 的 df2：", nrow(d), " 行；重复病人行：", sum(duplicated(d$patient_id)), "\n", sep = "")
d <- d[!duplicated(d$patient_id), ]
stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)
cols <- p2c[d$patient_id]; stopifnot(!anyNA(cols))
resA <- do.call(rbind, lapply(kg, function(g) {
  x <- as.numeric(em[g, cols]); r <- coxsd(d$OS_time, d$OS_event, x)
  data.frame(gene = g, HR_perSD = r[["HR"]], lo = r[["lo"]], hi = r[["hi"]], Cox_p = r[["p"]], KM_p = kmp(d$OS_time, d$OS_event, x))
}))
resA$Cox_p_bonf16 <- pmin(1, resA$Cox_p * nrow(resA))
resA <- resA[order(resA$Cox_p), ]
cat("去重后 166 例、80 事件；候选基因 ", nrow(resA), " 个\n", sep = "")
print(format(resA, digits = 3), row.names = FALSE)
cat("Cox p<0.05 的基因：", paste(resA$gene[resA$Cox_p < 0.05], collapse = ", "), "\n", sep = "")
write.csv(resA, file.path(OUTT, "candidate16_univariable_TCGA_O166.csv"), row.names = FALSE)

hr("A2. TIMP1 / PDPN 单基因对照（O 口径表，166 例）")
o <- read.csv(file.path(OUTT, "SCANI_TCGA_O_166.csv"))
stopifnot(nrow(o) == 166, !anyDuplicated(o$patient_id), sum(o$OS_event) == 80)
m <- match(o$patient_id, d$patient_id)
cat("O 表 TIMP1 与 step6d 表达的相关（应≈1，确认同一 aliquot / 同一变换）：",
    round(cor(o$TIMP1, as.numeric(em["TIMP1", cols[m]])), 4), "\n", sep = "")
resA2 <- do.call(rbind, lapply(c("SCANI", "TIMP1", "PDPN", "myCAF"), function(v) {
  r <- coxsd(o$OS_time, o$OS_event, o[[v]])
  data.frame(score = v, HR_perSD = r[["HR"]], lo = r[["lo"]], hi = r[["hi"]], Cox_p = r[["p"]],
             AUC3y = auc3(o$OS_time, o$OS_event, o[[v]]))
}))
print(format(resA2, digits = 3), row.names = FALSE)
cat("注：PDPN 的 HR<1 时 AUC<0.5 属方向问题，比较时看 |AUC−0.5|\n")
write.csv(resA2, file.path(OUTT, "single_gene_vs_SCANI_TCGA_O166.csv"), row.names = FALSE)

# ---------------------------------------------------------------- B
hr("B1. GSE72094：照 step22 读入（第 1–145 行，写操作已屏蔽）")
e22 <- sandbox(); run_lines(file.path(C, "scripts/step22_GSE72094_validation.R"), 145, e22)
annot <- get("annot", e22); lines <- get("lines", e22); su <- get("samples_use", e22); gsub_old <- get("gene_sub", e22)
cat("as-run gene_sub：", nrow(gsub_old), " 基因 x ", ncol(gsub_old), " 样本\n", sep = "")

hr("B2. 全基因组矩阵（同一注释、同一塌缩规则：同 symbol 取均值）")
b <- grep("^!series_matrix_table_begin", lines); e <- grep("^!series_matrix_table_end", lines)
mat <- read.delim(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(mat) <- mat[[1]]; mat <- as.matrix(mat[, -1]); storage.mode(mat) <- "numeric"
cat("探针矩阵：", nrow(mat), " x ", ncol(mat), "\n", sep = "")
sym <- annot$GeneSymbol[match(rownames(mat), annot$ID)]
keep <- !is.na(sym) & nzchar(sym); mat <- mat[keep, su, drop = FALSE]; sym <- sym[keep]
rs <- rowsum(mat, sym); full <- rs / as.vector(table(sym)[rownames(rs)])   # 同 symbol 取均值
full <- full[rowSums(is.na(full)) == 0, ]
cat("全基因组基因矩阵：", nrow(full), " 基因 x ", ncol(full), " 样本\n", sep = "")
chk <- intersect(rownames(gsub_old), rownames(full))
cat("与 as-run gene_sub 共有基因 ", length(chk), " 个，数值最大差：",
    signif(max(abs(full[chk, su] - gsub_old[chk, su])), 3), "（应为 0，确认塌缩一致）\n", sep = "")

sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]
stopifnot(length(sig) == 50)
use <- intersect(sig, rownames(full))
cat("全基因组下 myCAF 匹配：", length(use), "/50；缺失：", paste(setdiff(sig, use), collapse = ", "), "\n", sep = "")

suppressPackageStartupMessages(library(GSVA))
gs_full <- as.numeric(gsva(full, list(myCAF = use), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)["myCAF", su])
zs <- function(x) as.numeric(scale(x))

hr("B3. 两版分数并列的生存分析")
old <- read.csv(file.path(C, "results/09_GSE72094_validation/SCANI_GSE72094_scores.csv"), stringsAsFactors = FALSE)
old <- old[!duplicated(old$sample_id), ]
df <- data.frame(sample_id = su, myCAF_full = gs_full, TIMP1 = full["TIMP1", su], PDPN = full["PDPN", su])
df <- merge(df, old[, c("sample_id", "SCANI", "OS_days", "OS_event", "Stage", "age")], by = "sample_id")
df <- df[!is.na(df$OS_days) & df$OS_days > 0, ]
stopifnot(nrow(df) == 398, sum(df$OS_event) == 113)
df$SCANI_full <- zs(df$myCAF_full) + zs(df$TIMP1) - zs(df$PDPN)
names(df)[names(df) == "SCANI"] <- "SCANI_asrun"
cat("as-run SCANI 与全基因组 SCANI 的 Spearman：", round(cor(df$SCANI_asrun, df$SCANI_full, method = "spearman"), 3), "\n", sep = "")
df$stage_num <- suppressWarnings(as.numeric(gsub("[A-Za-z]", "", df$Stage)))   # 照抄 step22
cc <- !is.na(df$stage_num) & !is.na(df$age)
resB <- do.call(rbind, lapply(c("SCANI_asrun", "SCANI_full"), function(v) {
  u <- coxsd(df$OS_days, df$OS_event, df[[v]])
  dd <- df[cc, ]; dd$x <- zs(dd[[v]])                      # 与 precheck 一致：在完整病例上标准化
  mm <- summary(coxph(Surv(OS_days, OS_event) ~ x + age + stage_num, data = dd))
  data.frame(score = v, n = nrow(df), uni_HR = u[["HR"]], uni_lo = u[["lo"]], uni_hi = u[["hi"]], uni_p = u[["p"]],
             KM_p = kmp(df$OS_days, df$OS_event, df[[v]]),
             mv_n = mm$n, mv_HR = mm$conf.int[1, 1], mv_lo = mm$conf.int[1, 3], mv_hi = mm$conf.int[1, 4], mv_p = mm$coefficients[1, 5])
}))
print(format(resB, digits = 3), row.names = FALSE)
cat("核对：as-run 行应重现 HR 1.21 (1.00–1.46) p 0.044；KM 0.031；多变量 1.19 (0.99–1.43) p 0.066\n")
write.csv(resB, file.path(OUTT, "GSE72094_asrun_vs_fullgenome.csv"), row.names = FALSE)
write.csv(df, file.path(OUTT, "SCANI_GSE72094_both_versions.csv"), row.names = FALSE)
cat("\n完成。把 ", LOG, " 全文贴回 chat。\n", sep = "")
