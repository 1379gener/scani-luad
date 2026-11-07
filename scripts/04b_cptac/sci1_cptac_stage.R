#!/usr/bin/env Rscript
# ==============================================================================
# sci1_cptac_stage.R —— CPTAC-LUAD：SCANI 与病理分期的 Spearman 相关（稿件 Results / Table 1），按现行方法重算
#   背景：稿件现写 rho=0.15, P=0.21（n=72），来自 2026-03-15 的早期结果 CPTAC_CAF_stage_validation.csv（rho 0.148, P 0.214）；
#         服务器上另有 0.157/0.187（06_cptac_validation.R 目标值）与 -0.056/0.71（生存子集）。本脚本：
#   A. 核对旧结果：读 CPTAC_CAF_stage_validation.csv，复算 rho/P，检查其 CAF_Index 是否 = Z(myCAF)+Z(TIMP1)-Z(PDPN)
#   B. 按现行方法重算：表达 = 底座 A 的 CPTAC_LUAD_expr_symbol.rds（全基因组），分期 = CPTAC_LUAD_clinical_clean.csv；
#      myCAF = 50 基因 ssGSEA（GSVA 1.48.3，ssgsea.norm=FALSE，全基因组排序）；TIMP1、PDPN 取 log2(RPKM+1)；
#      三组分在「有表达且有分期」的肿瘤样本内 Z 标准化；分期 I–IV 编码 1–4；Spearman（exact=FALSE）
#   C. A 与 B 的逐样本比较（共同样本上 SCANI 的 Spearman）
#   只读；只写（已存在即停止）：tables/cptac/SCANI_CPTAC_stage.csv、tables/cptac/cptac_stage_spearman.csv
#   日志：/tmp/sci1_cptac_stage_<时间>.txt。结构不符合预期时打印信息后停止，不猜。
# 用法：Rscript ~/sci1_cptac_stage.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables/cptac")
F1 <- file.path(OD, "SCANI_CPTAC_stage.csv"); F2 <- file.path(OD, "cptac_stage_spearman.csv")
for (f in c(F1, F2)) if (file.exists(f)) stop("已存在，不覆盖：", f)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_cptac_stage_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
zs <- function(x) as.numeric(scale(x))
stopifnot(requireNamespace("GSVA", quietly = TRUE), as.character(packageVersion("GSVA")) == "1.48.3")
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)
stage_num <- function(s) { s <- toupper(trimws(as.character(s))); s <- sub("^STAGE\\s*", "", s); s <- sub("[ABC][0-9]?$", "", s)
  unname(c(I = 1, II = 2, III = 3, IV = 4)[s]) }

# ---------------------------------------------------------------- A. 旧结果
hr("A. 旧结果 CPTAC_CAF_stage_validation.csv（2026-03-15）")
fo <- file.path(C, "results/05_tcga_survival/CPTAC_CAF_stage_validation.csv")
old <- NULL
if (file.exists(fo)) {
  old <- fread(fo); cat("行数 ", nrow(old), "；列名：", paste(names(old), collapse = ", "), "\n", sep = ""); print(head(old, 3))
  ic <- grep("CAF_Index|SCANI", names(old), value = TRUE)[1]; sc <- grep("stage", names(old), ignore.case = TRUE, value = TRUE)
  cat("指数列：", ic, "；分期相关列：", paste(sc, collapse = ", "), "\n")
  for (s in sc) { v <- old[[s]]; vn <- if (is.numeric(v)) v else stage_num(v); ok <- !is.na(vn) & !is.na(old[[ic]])
    if (sum(ok) > 10) { ct <- cor.test(old[[ic]][ok], vn[ok], method = "spearman", exact = FALSE)
      cat(sprintf("  %s：n=%d，rho=%.4f，P=%.4f\n", s, sum(ok), ct$estimate, ct$p.value)) } }
  if (all(c("myCAF", "TIMP1", "PDPN") %in% names(old))) {
    re <- zs(old$myCAF) + zs(old$TIMP1) - zs(old$PDPN)
    cat(sprintf("旧 %s 与 Z(myCAF)+Z(TIMP1)-Z(PDPN) 的最大差：%.3g（≈0 说明公式一致、在这些样本内标准化）\n", ic, max(abs(re - old[[ic]]))))
  }
} else cat("未找到 ", fo, "\n")

# ---------------------------------------------------------------- B. 现行方法重算
hr("B. 按现行方法重算（底座 A 的 CPTAC 表达与临床）")
DA <- file.path(H, "thesis_LUAD_CAF/01_bases/A_bulk_OS/derived/CPTAC_LUAD")
mf <- file.path(DA, "cohort_manifest_CPTAC.csv"); if (file.exists(mf)) { cat("cohort_manifest：\n"); print(fread(mf)) }
ex <- readRDS(file.path(DA, "CPTAC_LUAD_expr_symbol.rds"))
cat("表达对象：class=", paste(class(ex), collapse = "/"), "；dim=", paste(dim(ex), collapse = " x "), "\n", sep = "")
if (is.data.frame(ex)) { rn <- if (!is.null(rownames(ex)) && !all(rownames(ex) == seq_len(nrow(ex)))) rownames(ex) else ex[[1]]
  num <- vapply(ex, is.numeric, TRUE); ex <- as.matrix(ex[, num, drop = FALSE]); rownames(ex) <- rn }
stopifnot(is.matrix(ex), is.numeric(ex), !is.null(rownames(ex)), !anyDuplicated(rownames(ex)))
cat("行名示例：", paste(head(rownames(ex), 3), collapse = ", "), "；列名示例：", paste(head(colnames(ex), 3), collapse = ", "), "\n")
cat(sprintf("数值范围：%.3g – %.3g；中位数 %.3g；NA 比例 %.3g\n", min(ex, na.rm = TRUE), max(ex, na.rm = TRUE), median(ex, na.rm = TRUE), mean(is.na(ex))))
if (max(ex, na.rm = TRUE) > 50) { cat("★ 最大值 >50，判定为线性 RPKM，按 Methods 取 log2(x+1)\n"); ex <- log2(ex + 1) } else cat("判定为已是 log 尺度，不再变换\n")
miss <- setdiff(c(sig, "TIMP1", "PDPN"), rownames(ex)); cat("50 个签名基因匹配：", sum(sig %in% rownames(ex)), "/50；缺：", paste(miss, collapse = ","), "\n")
stopifnot(all(c("TIMP1", "PDPN") %in% rownames(ex)))
cl <- fread(file.path(DA, "CPTAC_LUAD_clinical_clean.csv")); cat("临床表：", nrow(cl), " 行；列名：", paste(names(cl), collapse = ", "), "\n", sep = "")
idc <- intersect(c("sample_id", "sample", "patient_id", "case_id", "id"), names(cl)); stc <- intersect(c("stage", "Stage", "stage_simple"), names(cl))
if (!length(idc) || !length(stc)) stop("临床表里找不到样本 id 或 stage 列，请把上面的列名贴回")
cl <- cl[!duplicated(cl[[idc[1]]])]; cl$stage_n <- stage_num(cl[[stc[1]]])
cat("分期原值分布：\n"); print(table(cl[[stc[1]]], useNA = "ifany"))
common <- intersect(colnames(ex), cl[[idc[1]]][!is.na(cl$stage_n)])
cat("表达列 ", ncol(ex), "；有分期 ", sum(!is.na(cl$stage_n)), "；两者共同 ", length(common), "\n", sep = "")
if (length(common) < 30) stop("样本 id 对不上（共同样本 <30）。表达列名与临床 id 示例见上，请贴回")
X <- ex[, common, drop = FALSE]; X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
cat("去掉含 NA 的基因后：", nrow(X), " 个基因；签名基因 ", sum(sig %in% rownames(X)), "/50\n", sep = "")
my <- as.numeric(GSVA::gsva(X, list(myCAF = intersect(sig, rownames(X))), method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)["myCAF", ])
d <- data.frame(sample_id = common, myCAF = my, TIMP1 = as.numeric(X["TIMP1", ]), PDPN = as.numeric(X["PDPN", ]),
                stage = cl[[stc[1]]][match(common, cl[[idc[1]]])], stage_n = cl$stage_n[match(common, cl[[idc[1]]])])
d$SCANI <- zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)
ct <- cor.test(d$SCANI, d$stage_n, method = "spearman", exact = FALSE)
cat(sprintf("\n★ 现行方法：n=%d，SCANI 与分期 Spearman rho=%.4f，P=%.4f\n", nrow(d), ct$estimate, ct$p.value))
print(table(d$stage_n)); print(round(tapply(d$SCANI, d$stage_n, median), 3))
comp <- lapply(c("myCAF", "TIMP1", "PDPN"), function(v) { x <- cor.test(d[[v]], d$stage_n, method = "spearman", exact = FALSE)
  data.frame(score = v, n = nrow(d), rho = unname(x$estimate), p = x$p.value) })
R <- rbind(data.frame(score = "SCANI", n = nrow(d), rho = unname(ct$estimate), p = ct$p.value), do.call(rbind, comp))
R$signature_genes_matched <- sum(sig %in% rownames(X)); R$genes_in_universe <- nrow(X)
print(R)

# ---------------------------------------------------------------- C. 新旧逐样本比较
hr("C. 新旧比较")
if (!is.null(old)) { oid <- intersect(c("sample_id", "sample", "Sample", "patient_id", "id"), names(old))[1]; ic <- grep("CAF_Index|SCANI", names(old), value = TRUE)[1]
  if (!is.na(oid)) { m <- match(d$sample_id, old[[oid]]); cat("共同样本：", sum(!is.na(m)), "\n")
    if (sum(!is.na(m)) > 10) { cat(sprintf("新 SCANI 与旧 %s 的 Spearman：%.3f\n", ic, cor(d$SCANI[!is.na(m)], old[[ic]][m[!is.na(m)]], method = "spearman")))
      for (v in c("myCAF", "TIMP1", "PDPN")) if (v %in% names(old)) cat(sprintf("  %s 新旧 Spearman：%.3f\n", v, cor(d[[v]][!is.na(m)], old[[v]][m[!is.na(m)]], method = "spearman"))) }
  } else cat("旧表无可识别的样本 id 列\n") }

dir.create(OD, showWarnings = FALSE, recursive = TRUE)
write.csv(d, F1, row.names = FALSE); write.csv(R, F2, row.names = FALSE)
cat("\n写出：", F1, "\n      ", F2, "\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
