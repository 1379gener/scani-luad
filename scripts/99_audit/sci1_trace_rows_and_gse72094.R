#!/usr/bin/env Rscript
# ==============================================================================
# sci1_trace_rows_and_gse72094.R —— 只读溯源：174 / 192 / 166 的来历；GSE72094 的 55 探针 / 21 基因 / 12 签名基因
#   不写任何结果文件；只写日志 /tmp/sci1_trace_<时间>.txt
#   A. step6d（组分筛选）的 df2：行数、病人数、哪些病人重复、各重复几行
#   B. step6e 产物 CAF_final_scores.csv（打分表）：行数、重复结构；TCGA-44-2665 的 27 行里
#      不同的 (myCAF, TIMP1, PDPN) 取值有几组，按 patient_id 去重时保留的是哪一组（=现行 O 口径所用）
#   C. 两表的 166 例是否为同一批病人；O 口径表与 B 的去重结果是否逐值一致
#   D. step22（GSE72094 早期实现）：55 个探针 → 21 个基因的完整列表；其中哪些是 myCAF 签名基因、哪些是 TIMP1/PDPN、其余是什么
# 用法：Rscript ~/sci1_trace_rows_and_gse72094.R
# ==============================================================================
suppressPackageStartupMessages(library(survival))
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
S1T <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
LOG <- file.path("/tmp", format(Sys.time(), "sci1_trace_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
sandbox <- function() {
  e <- new.env(parent = globalenv()); noop <- function(...) invisible(NULL)
  for (f in c("write.csv", "write.table", "saveRDS", "save", "writeLines", "pdf", "png", "tiff", "jpeg",
              "ggsave", "dev.off", "fwrite", "write_csv", "sink")) assign(f, noop, envir = e)
  assign("dir.create", function(...) invisible(TRUE), envir = e); e
}
run_lines <- function(file, upto, env) eval(parse(text = readLines(file, warn = FALSE)[seq_len(upto)]), envir = env)
dupinfo <- function(ids) { t <- table(ids); t <- t[t > 1]; if (length(t)) print(sort(t, decreasing = TRUE)) else cat("无重复\n") }
cat("SCI-1 溯源  ", format(Sys.time()), "\n日志：", LOG, "\n", sep = "")

hr("A. step6d_genes.R（TIMP1/PDPN 的筛选输入）")
f6d <- file.path(C, "scripts/step6d_genes.R"); cat("md5：", unname(tools::md5sum(f6d)), "\n", sep = "")
L <- readLines(f6d, warn = FALSE)
cat("第 1–37 行中读文件 / 合并 / 筛选的语句：\n")
idx <- grep("read|merge|subset|filter|stage|OS_time|pid_to_col|df2", L[1:37]); cat(sprintf("  L%-3d %s\n", idx, L[idx]), sep = "")
e6 <- sandbox(); run_lines(f6d, 37, e6)
d6 <- get("df2", e6); p2c <- get("pid_to_col", e6)
cat("\ndf2：", nrow(d6), " 行，", length(unique(d6$patient_id)), " 个病人，事件（按行）", sum(d6$OS_event), "\n", sep = "")
cat("重复病人及行数：\n"); dupinfo(d6$patient_id)
cat("单基因取表达列的规则：pid_to_col 每病人对应列数 = ", paste(unique(lengths(p2c)), collapse = "/"), "\n", sep = "")
if ("TCGA-44-2665" %in% names(p2c)) cat("TCGA-44-2665 → 列 ", paste(p2c[["TCGA-44-2665"]], collapse = ", "), "\n", sep = "")

hr("B. CAF_final_scores.csv（step6e 打分表）")
fs <- file.path(C, "results/05_tcga_survival/CAF_final_scores.csv"); cat("md5：", unname(tools::md5sum(fs)), "\n", sep = "")
sc <- read.csv(fs, stringsAsFactors = FALSE)
cat(nrow(sc), " 行，", length(unique(sc$patient_id)), " 个病人；列：", paste(names(sc), collapse = ", "), "\n", sep = "")
cat("重复病人及行数：\n"); dupinfo(sc$patient_id)
x <- sc[sc$patient_id == "TCGA-44-2665", intersect(c("patient_id", "sample_id", "aliquot", "myCAF", "TIMP1", "PDPN", "OS_time", "OS_event", "stage_simple"), names(sc))]
cat("\nTCGA-44-2665 共 ", nrow(x), " 行；不同的 (myCAF, TIMP1, PDPN) 组合：\n", sep = "")
print(unique(x[, intersect(c("myCAF", "TIMP1", "PDPN"), names(x))]))
cat("不同的 (OS_time, OS_event, stage_simple) 组合：\n"); print(unique(x[, intersect(c("OS_time", "OS_event", "stage_simple"), names(x))]))
cat("按 patient_id 去重保留的第一行（现行 O 口径所用）：\n"); print(x[1, ])

hr("C. 两表与 O 口径表的一致性")
s166 <- sc[!duplicated(sc$patient_id), ]; d166 <- d6[!duplicated(d6$patient_id), ]
cat("打分表去重后 ", nrow(s166), " 例，事件 ", sum(s166$OS_event), "；筛选表去重后 ", nrow(d166), " 例，事件 ", sum(d166$OS_event), "\n", sep = "")
cat("两表病人集合相同：", setequal(s166$patient_id, d166$patient_id), "\n", sep = "")
o <- read.csv(file.path(S1T, "SCANI_TCGA_O_166.csv"), stringsAsFactors = FALSE)
m <- merge(o, s166, by = "patient_id", suffixes = c(".O", ".src"))
for (v in c("myCAF", "TIMP1", "PDPN")) cat(v, " 与 O 表最大绝对差：", signif(max(abs(m[[paste0(v, ".O")]] - m[[paste0(v, ".src")]])), 3), "\n", sep = "")
em <- get("expr_mat", e6)
tim <- as.numeric(em["TIMP1", unlist(lapply(p2c[o$patient_id], `[`, 1))])
cat("O 表 TIMP1 与 step6d 每病人第一列表达的 Spearman：", round(cor(o$TIMP1, tim, method = "spearman"), 5), "\n", sep = "")

hr("D. GSE72094 早期实现（step22）的 55 探针 / 21 基因")
f22 <- file.path(C, "scripts/step22_GSE72094_validation.R"); cat("md5：", unname(tools::md5sum(f22)), "\n", sep = "")
L <- readLines(f22, warn = FALSE)
idx <- grep("gene_sub|genes_needed|key_genes|%in%|probe|myCAF|TIMP1|PDPN|gsva", L[1:145]); cat(sprintf("  L%-3d %s\n", idx, L[idx]), sep = "")
e22 <- sandbox(); run_lines(f22, 145, e22)
gsub_old <- get("gene_sub", e22); annot <- get("annot", e22)
sig <- read.csv(file.path(C, "results/01_scRNA/CAF_signatures.csv")); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)
g21 <- rownames(gsub_old)
cat("\ngene_sub：", length(g21), " 基因 × ", ncol(gsub_old), " 样本\n", sep = "")
cat("其中 myCAF 签名基因（", sum(g21 %in% sig), "）：", paste(sort(g21[g21 %in% sig]), collapse = ", "), "\n", sep = "")
cat("TIMP1 / PDPN：", paste(intersect(c("TIMP1", "PDPN"), g21), collapse = ", "), "\n", sep = "")
cat("其余：", paste(sort(setdiff(g21, c(sig, "TIMP1", "PDPN"))), collapse = ", "), "\n", sep = "")
pr <- annot$ID[annot$GeneSymbol %in% g21]
cat("注释表中对应这 ", length(g21), " 个基因的探针数：", length(pr), "\n", sep = "")
print(table(annot$GeneSymbol[annot$GeneSymbol %in% g21]))
cat("\n完成。把 ", LOG, " 全文贴回。\n", sep = "")
