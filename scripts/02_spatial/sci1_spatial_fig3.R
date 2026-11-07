#!/usr/bin/env Rscript
# ==============================================================================
# sci1_spatial_fig3.R —— 空间 50 基因：补全代码出处 + 从逐切片 DE 表复算 + 重跑富集 + Fig 3
#
#   A. 打印 step3_spatial_GSE307534.R 的模块打分（约 L60–185）与基因整合（L274–315）原文
#   B. 列出所有做过 enrichGO / enrichKEGG 的脚本行（旧富集结果的出处）
#   C. 从 CAF_Niche_DEGs_all_LUAD.csv 按代码规则复算 778 → 50，断言与 geneset.txt 一致
#   D. 用当前 50 基因重新跑 GO-BP 与 KEGG（旧 GO/KEGG csv 早于当前 geneset，不再使用）
#   E. Fig 3：A 基因在切片中的复现次数分布；B GO-BP 前 10；C KEGG 前 10（只读 C/D 写出的 csv）
#
#   只读 caf_spatial_project；只写 ~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/{tables,figures}
# 用法：Rscript ~/sci1_spatial_fig3.R
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
H   <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
SP  <- file.path(C, "results/03_spatial_GSE307534")
S3  <- file.path(C, "scripts/step3_spatial_GSE307534.R")
OUTD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS")
TB <- function(f) file.path(OUTD, "tables", f)
LOG <- file.path("/tmp", format(Sys.time(), "sci1_spatial_fig3_%Y%m%d_%H%M.txt"))
zz <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70), "\n", sep = "")
cat("空间 Fig 3  ", format(Sys.time()), "\n日志: ", LOG, "\n", sep = "")

# ---------------------------------------------------------------- A. 代码原文
hr("A. step3_spatial_GSE307534.R 原文")
x <- readLines(S3, warn = FALSE)
cat("md5 ", tools::md5sum(S3), "；mtime ", format(file.mtime(S3)), "\n", sep = "")
for (r in list(60:185, 274:315)) {
  r <- r[r <= length(x)]
  cat("\n---- L", min(r), "–", max(r), "\n", sep = ""); cat(sprintf("%5d  %s", r, x[r]), sep = "\n")
}

# ---------------------------------------------------------------- B. 旧富集出处
hr("B. 含 enrichGO / enrichKEGG 的脚本")
scr <- list.files(file.path(C, "scripts"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
for (f in scr) {
  y <- readLines(f, warn = FALSE); i <- grep("enrichGO|enrichKEGG|universe", y)
  i <- i[!grepl("^\\s*#", y[i])]
  if (length(i)) cat(sprintf("%s  L%d: %s\n", sub(H, "~", f), i, trimws(substr(y[i], 1, 140))), sep = "")
}
old <- file.path(SP, c("GO_BP_enrichment.csv", "KEGG_enrichment.csv", "CAF_Niche_final_geneset.txt"))
cat("\n旧文件时间：\n"); print(data.frame(file = basename(old), mtime = format(file.mtime(old))), row.names = FALSE)

# ---------------------------------------------------------------- C. 复算 50 基因
hr("C. 从逐切片 DE 表复算")
de <- read.csv(file.path(SP, "CAF_Niche_DEGs_all_LUAD.csv"), stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(paste(de$sample, de$gene)))
nsec <- length(unique(de$sample)); stopifnot(nsec == 26)
min_s <- max(2, floor(nsec / 3))                       # 与 step3 L284 相同
sig <- de[de$p_val_adj < 0.05 & de$avg_log2FC > 0.25, ]
agg <- data.frame(gene = names(table(sig$gene)), n_sections = as.integer(table(sig$gene)))
agg$mean_log2FC <- as.numeric(tapply(sig$avg_log2FC, sig$gene, mean)[agg$gene])
agg <- agg[order(-agg$mean_log2FC), ]
rec <- agg[agg$n_sections >= min_s, ]
top <- head(rec, 50)
gs  <- trimws(readLines(file.path(SP, "CAF_Niche_final_geneset.txt"), warn = FALSE)); gs <- gs[nzchar(gs)]
cat("切片数 ", nsec, "；min_samples = ", min_s, "；≥", min_s, " 切片复现 ", nrow(rec),
    " 个；前 50 与 geneset 集合相等 = ", setequal(top$gene, gs), "\n", sep = "")
stopifnot(setequal(top$gene, gs))
write.csv(agg, TB("spatial_gene_recurrence_all.csv"), row.names = FALSE)
write.csv(transform(top, rank = seq_len(50)), TB("spatial_geneset50_derivation.csv"), row.names = FALSE)
myc <- readRDS(file.path(C, "results/01_scRNA/CAF_subtype_signatures.rds"))$myCAF
cat("与 SCANI myCAF 50 基因重叠：", length(intersect(gs, myc)), "（", paste(intersect(gs, myc), collapse = ", "), "）\n", sep = "")

# ---------------------------------------------------------------- D. 富集
hr("D. GO-BP / KEGG 富集（当前 50 基因）")
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
cat("clusterProfiler ", as.character(packageVersion("clusterProfiler")),
    "；org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")), "\n", sep = "")
id <- bitr(gs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("SYMBOL→ENTREZ 映射 ", length(unique(id$SYMBOL)), "/50；未映射：",
    paste(setdiff(gs, id$SYMBOL), collapse = ", "), "\n", sep = "")
go <- enrichGO(id$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_df <- as.data.frame(go); write.csv(go_df, TB("spatial_GO_BP_enrichment.csv"), row.names = FALSE)
cat("GO-BP 显著条目 ", nrow(go_df), "\n", sep = ""); print(head(go_df[, c("Description", "GeneRatio", "p.adjust")], 10), row.names = FALSE)

kg <- tryCatch(enrichKEGG(id$ENTREZID, organism = "hsa", pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2),
               error = function(e) { cat("KEGG 失败（多为无法联网）：", conditionMessage(e), "\n"); NULL })
if (!is.null(kg)) {
  kg <- setReadable(kg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kg_df <- as.data.frame(kg); write.csv(kg_df, TB("spatial_KEGG_enrichment.csv"), row.names = FALSE)
  cat("KEGG 查询日期 ", format(Sys.Date()), "；显著条目 ", nrow(kg_df), "\n", sep = "")
  print(head(kg_df[, c("Description", "GeneRatio", "p.adjust")], 10), row.names = FALSE)
}
cat("\n（仅供对照，不进稿）旧 GO_BP_enrichment.csv 前 3 行：\n")
ob <- read.csv(old[1], stringsAsFactors = FALSE)
print(head(ob[, intersect(c("Description", "GeneRatio", "p.adjust"), names(ob))], 3), row.names = FALSE)

# ---------------------------------------------------------------- E. Fig 3
hr("E. Fig 3")
wrap <- function(s, w = 34) vapply(s, function(z) paste(strwrap(z, w), collapse = "\n"), "")
ra <- read.csv(TB("spatial_gene_recurrence_all.csv"))
pA <- ggplot(ra, aes(n_sections)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "white", boundary = 0.5) +
  geom_vline(xintercept = min_s - 0.5, linetype = 2, colour = "#B2182B") +
  annotate("text", x = min_s, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
           label = sprintf("Recurrent in at least %d sections: %d genes", min_s, sum(ra$n_sections >= min_s))) +
  scale_x_continuous(breaks = seq(0, nsec, 2)) +
  labs(x = sprintf("Number of invasive LUAD sections (of %d) with significant up-regulation", nsec),
       y = "Genes") + theme_classic(base_size = 11) + theme(plot.margin = margin(20, 10, 5, 10))

enr_panel <- function(f, ttl) {
  if (!file.exists(f)) return(ggplot() + annotate("text", 0, 0, label = paste(ttl, ": not available")) + theme_void())
  e <- read.csv(f, stringsAsFactors = FALSE)
  if (!nrow(e)) return(ggplot() + annotate("text", 0, 0, label = paste(ttl, ": no significant terms")) + theme_void())
  e <- head(e[order(e$p.adjust), ], 10)
  e$Description <- factor(wrap(e$Description), levels = rev(wrap(e$Description)))
  ggplot(e, aes(-log10(p.adjust), Description)) +
    geom_col(fill = "#4393C3", width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.3, size = 2.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "-log10 (adjusted p)", y = NULL, title = ttl) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.margin = margin(5, 10, 5, 25))
}
pB <- enr_panel(TB("spatial_GO_BP_enrichment.csv"), "GO biological process")
pC <- enr_panel(TB("spatial_KEGG_enrichment.csv"), "KEGG pathway")
fig3 <- plot_grid(pA, plot_grid(pB, pC, nrow = 1, labels = c("B", "C"), label_size = 14),
                  ncol = 1, labels = c("A", ""), label_size = 14, rel_heights = c(0.7, 1))
for (ext in c(".tiff", ".png", ".pdf")) {
  f <- file.path(OUTD, "figures", paste0("Figure3", ext))
  if (ext == ".tiff") ggsave(f, fig3, width = 12, height = 8.5, dpi = 600, compression = "lzw", bg = "white")
  else if (ext == ".png") ggsave(f, fig3, width = 12, height = 8.5, dpi = 600, bg = "white")
  else ggsave(f, fig3, width = 12, height = 8.5, bg = "white")
}
cat("Figure 3 已写出。把 ", LOG, " 全文贴回，并发 Figure3.png。\n", sep = "")
