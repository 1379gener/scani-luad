#!/usr/bin/env Rscript
# ==============================================================================
# sci1_supp_export.R —— 导出 Table S1 / S2 所需数据（只读；打印为 CSV 文本，便于整段贴回）
#   S1a  16 个候选基因去重后单变量结果        tables/candidate16_univariable_TCGA_O166.csv
#   S1b  SCANI 与单组分对照                    tables/single_gene_vs_SCANI_TCGA_O166.csv
#   S1c  myCAF 50 基因签名 + FindAllMarkers 统计  results/01_scRNA/CAF_signatures.csv + CAF_DEGs_all_subtypes.csv
#   S2   GSE72094 两版结果                     tables/GSE72094_asrun_vs_fullgenome.csv
#   同时把全部内容写到 tables/supp_export_for_S1_S2.txt
# 用法：Rscript ~/sci1_supp_export.R
# ==============================================================================
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
TB <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/tables")
OUT <- file.path(TB, "supp_export_for_S1_S2.txt")
zz <- file(OUT, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
sec <- function(t, f) cat("\n### ", t, "\n# file: ", sub(H, "~", f), "  md5 ", tools::md5sum(f), "\n", sep = "")
dump <- function(d) cat(capture.output(write.csv(d, row.names = FALSE)), sep = "\n")

f1 <- file.path(TB, "candidate16_univariable_TCGA_O166.csv"); sec("S1a candidate genes (n=166)", f1)
d1 <- read.csv(f1); stopifnot(nrow(d1) == 16); dump(d1)

f2 <- file.path(TB, "single_gene_vs_SCANI_TCGA_O166.csv"); sec("S1b SCANI vs components (n=166)", f2)
d2 <- read.csv(f2); stopifnot(nrow(d2) == 4); dump(d2)

f3 <- file.path(C, "results/01_scRNA/CAF_signatures.csv"); f4 <- file.path(C, "results/01_scRNA/CAF_DEGs_all_subtypes.csv")
sec("S1c myCAF signature with FindAllMarkers statistics", f4)
sig <- read.csv(f3); sig <- sig$Gene[sig$Subtype == "myCAF"]; stopifnot(length(sig) == 50)
de <- read.csv(f4, stringsAsFactors = FALSE)
cat("# DEG 表列名：", paste(names(de), collapse = ", "), "\n", sep = "")
gcol <- intersect(c("gene", "Gene", "X"), names(de))[1]; ccol <- intersect(c("cluster", "Subtype", "subtype"), names(de))[1]
cat("# 分组取值：", paste(unique(de[[ccol]]), collapse = ", "), "\n", sep = "")
m <- de[de[[ccol]] == "myCAF" & de[[gcol]] %in% sig, ]
m <- m[match(sig, m[[gcol]]), ]
stopifnot(nrow(m) == 50, !anyNA(m[[gcol]]))
m$rank <- seq_len(50)
dump(m[, c("rank", gcol, intersect(c("avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"), names(m)))])

f5 <- file.path(TB, "GSE72094_asrun_vs_fullgenome.csv"); sec("S2 GSE72094 genome-wide (primary) vs 21-gene (sensitivity)", f5)
d5 <- read.csv(f5); stopifnot(nrow(d5) == 2); dump(d5)

cat("\n# 完成。整段贴回 chat；同内容已写入 ", sub(H, "~", OUT), "\n", sep = "")
