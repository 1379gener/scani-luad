#!/usr/bin/env Rscript
# ==============================================================================
# sci1_fig2_and_zorder.R —— SCI-1：新版 Fig 2 + 补充表 S1 + Z 标准化顺序敏感性分析
#
#   第 1 部分（会写文件，只写入新目录 ~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/）
#     Fig 2A  2,992 个成纤维细胞按组织来源 × 亚群的构成（堆叠柱，数据全部来自 meta.data）
#     Fig 2B  三个亚群的 marker 点图（基因取自各亚群签名文件的前 5 个，不手写）
#     S1 Table  SCANI 使用的 myCAF 50 基因及其 FindAllMarkers 统计量
#     fig2_counts_origin_by_subtype.csv  供 verify 核对
#
#   第 2 部分（只读）：TCGA 的 Z 标准化在去重前（现行，变体 N）还是去重后（变体 O）
#     并排输出两种口径下的：单变量 per-SD HR、中位分组 KM p、3 年 AUC、多变量 HR
#
# 用法：Rscript ~/sci1_fig2_and_zorder.R
# 输出：/tmp/sci1_fig2_zorder_<时间>.txt（贴回 chat）+ 上述图表文件
# ==============================================================================

suppressPackageStartupMessages({ library(survival) })
H    <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS")
dir.create(file.path(OUTD, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTD, "tables"),  recursive = TRUE, showWarnings = FALSE)
LOG  <- file.path("/tmp", format(Sys.time(), "sci1_fig2_zorder_%Y%m%d_%H%M.txt"))
zz   <- file(LOG, open = "wt"); sink(zz, split = TRUE); on.exit({ sink(); close(zz) }, add = TRUE)
hr   <- function(t) cat("\n", strrep("=", 72), "\n", t, "\n", strrep("=", 72), "\n", sep = "")
safe <- function(expr, label) tryCatch(expr, error = function(e) {
  cat("  !! ", label, " 失败: ", conditionMessage(e), "\n", sep = ""); invisible(NULL) })
cat("SCI-1 Fig 2 + Z 顺序敏感性\n时间: ", format(Sys.time()), "\n日志: ", LOG, "\n输出目录: ", OUTD, "\n", sep = "")

# ---- 标签映射（文字标签，不是数据） ------------------------------------------
SUBTYPE_MAP <- c("Myofibroblasts"      = "Myofibroblasts (myCAF)",
                 "COL13A1+ matrix FBs" = "Matrix FBs (iCAF-like)",
                 "COL14A1+ matrix FBs" = "Matrix FBs (iCAF-like)",
                 "FB-like cells"       = "FB-like cells")
ORIGIN_MAP  <- c("tLung" = "Primary tumour", "tL/B" = "Primary tumour",
                 "nLung" = "Normal lung", "mBrain" = "Brain metastasis", "mLN" = "Lymph node metastasis")
SUB_LEVELS  <- c("Myofibroblasts (myCAF)", "Matrix FBs (iCAF-like)", "FB-like cells")
ORI_LEVELS  <- c("Primary tumour", "Normal lung", "Lymph node metastasis", "Brain metastasis")

# ==============================================================================
hr("1. Fig 2 与 S1 Table")
# ==============================================================================
safe({
  suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
  obj <- readRDS(file.path(C, "results/01_scRNA/seu_CAF_annotated.rds"))
  stopifnot(ncol(obj) == 2992)
  md <- obj@meta.data
  unk_s <- setdiff(unique(md$Cell_subtype), names(SUBTYPE_MAP))
  unk_o <- setdiff(unique(md$Sample_Origin), names(ORIGIN_MAP))
  if (length(unk_s) || length(unk_o)) stop("未映射的标签: ", paste(c(unk_s, unk_o), collapse = ", "))
  obj$FB_group <- factor(unname(SUBTYPE_MAP[md$Cell_subtype]), levels = SUB_LEVELS)
  obj$Origin   <- factor(unname(ORIGIN_MAP[md$Sample_Origin]), levels = ORI_LEVELS)

  tab <- as.data.frame(table(Origin = obj$Origin, Group = obj$FB_group))
  write.csv(tab, file.path(OUTD, "tables/fig2_counts_origin_by_subtype.csv"), row.names = FALSE)
  cat("  来源 × 亚群计数:\n"); print(addmargins(table(obj$Origin, obj$FB_group)))
  stopifnot(sum(tab$Freq) == 2992)

  # --- Fig 2A：堆叠柱（每个来源一根，柱内为亚群构成），柱顶标注 n
  tot <- aggregate(Freq ~ Origin, tab, sum)
  pA <- ggplot(tab, aes(x = Origin, y = Freq, fill = Group)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_text(data = tot, aes(x = Origin, y = Freq, label = paste0("n=", Freq)),
              inherit.aes = FALSE, vjust = -0.4, size = 3.2) +
    scale_fill_manual(values = setNames(c("#B2182B", "#4393C3", "#7F7F7F"), SUB_LEVELS), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Number of fibroblasts") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")

  # --- Fig 2B：marker 点图；基因 = 三个签名文件各自的前 5 个（数据驱动）
  sigs <- readRDS(file.path(C, "results/01_scRNA/CAF_subtype_signatures.rds"))
  feats <- unique(c(head(sigs$myCAF, 5), head(sigs$iCAF, 5), head(sigs$apCAF, 5)))
  feats <- intersect(feats, rownames(obj))
  cat("\n  点图基因（各签名前 5）: ", paste(feats, collapse = ", "), "\n", sep = "")
  Idents(obj) <- "FB_group"
  pB <- DotPlot(obj, features = feats, group.by = "FB_group") + RotatedAxis() +
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(size = 8))

  for (nm in c("Fig2A", "Fig2B")) {
    p <- if (nm == "Fig2A") pA else pB
    w <- if (nm == "Fig2A") 5.5 else 8; h <- 4.5
    ggsave(file.path(OUTD, "figures", paste0(nm, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw")
    ggsave(file.path(OUTD, "figures", paste0(nm, ".png")),  p, width = w, height = h, dpi = 600)
    ggsave(file.path(OUTD, "figures", paste0(nm, ".pdf")),  p, width = w, height = h)
  }
  cat("  已写出 Fig2A / Fig2B（tiff 600 dpi + png + pdf）\n")

  # --- S1 Table：SCANI 实际使用的 myCAF 50 基因
  deg <- read.csv(file.path(C, "results/01_scRNA/CAF_DEGs_all_subtypes.csv"), check.names = FALSE)
  s1  <- deg[deg$cluster == "myCAF" & deg$gene %in% sigs$myCAF, ]
  s1  <- s1[match(sigs$myCAF, s1$gene), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
  s1$rank <- seq_len(nrow(s1))
  stopifnot(nrow(s1) == 50, !anyNA(s1$gene))
  write.csv(s1, file.path(OUTD, "tables/S1_Table_myCAF_signature_50genes.csv"), row.names = FALSE)
  cat("  已写出 S1 Table（50 基因；PDPN 排第 ", which(s1$gene == "PDPN"), " 位）\n", sep = "")
}, "Fig 2 / S1")

# ==============================================================================
hr("2. Z 标准化顺序敏感性（TCGA，只读）")
# ==============================================================================
safe({
  f <- file.path(C, "results/05_tcga_survival/CAF_final_scores.csv")
  raw <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  cat("  列名: ", paste(names(raw), collapse = ", "), "\n", sep = "")
  cat("  原始行数: ", nrow(raw), "\n", sep = "")
  need <- c("patient_id", "myCAF", "TIMP1", "PDPN", "CAF_Index", "OS_time", "OS_event")
  miss <- setdiff(need, names(raw)); if (length(miss)) stop("缺列: ", paste(miss, collapse = ", "))

  # 变体 N（现行）：Z 在 192 行上做 → 去重（保留首条，与金标准复现脚本一致）→ 166
  N <- raw[!duplicated(raw$patient_id), ]
  # 变体 O（正确顺序）：先去重（同样保留首条）→ 在 166 行上 Z
  O <- N
  O$CAF_Index <- as.numeric(scale(O$myCAF) + scale(O$TIMP1) - scale(O$PDPN))
  stopifnot(nrow(N) == 166, nrow(O) == 166, sum(N$OS_event) == 80)
  cat("  N 中 CAF_Index 与文件原值一致性（应为 TRUE）: ",
      isTRUE(all.equal(N$CAF_Index, raw$CAF_Index[!duplicated(raw$patient_id)])), "\n", sep = "")
  cat("  N 与 O 的 SCANI Spearman: ", sprintf("%.5f", cor(N$CAF_Index, O$CAF_Index, method = "spearman")), "\n", sep = "")

  stc <- grep("stage", names(raw), ignore.case = TRUE, value = TRUE)
  agc <- grep("^age", names(raw), ignore.case = TRUE, value = TRUE)
  sxc <- grep("^(sex|gender)", names(raw), ignore.case = TRUE, value = TRUE)
  cat("  协变量列: stage=", paste(stc, collapse = "/"), " age=", paste(agc, collapse = "/"),
      " sex=", paste(sxc, collapse = "/"), "\n", sep = "")

  res <- function(d, lab) {
    d$z <- as.numeric(scale(d$CAF_Index))
    u   <- summary(coxph(Surv(OS_time, OS_event) ~ z, d))
    km  <- survdiff(Surv(OS_time, OS_event) ~ I(CAF_Index >= median(CAF_Index)), d)
    auc <- NA
    if (requireNamespace("timeROC", quietly = TRUE)) {
      r <- timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = d$CAF_Index,
                            cause = 1, times = 3 * 365, iid = FALSE)
      auc <- unname(r$AUC[2])
    }
    mvs <- "（无协变量列）"
    if (length(stc) && length(agc) && length(sxc)) {
      dd <- d; dd$.st <- factor(dd[[stc[1]]]); dd$.ag <- as.numeric(dd[[agc[1]]]); dd$.sx <- factor(dd[[sxc[1]]])
      m  <- summary(coxph(Surv(OS_time, OS_event) ~ z + .st + .ag + .sx, dd))
      mvs <- sprintf("HR %.3f (%.3f–%.3f), p=%.4f, n=%d", m$conf.int["z", 1], m$conf.int["z", 3],
                     m$conf.int["z", 4], m$coefficients["z", 5], m$n)
    }
    cat(sprintf("\n  [%s]\n    单变量 per SD : HR %.3f (%.3f–%.3f), p=%.4f\n    中位分组 KM   : p=%.4f\n    3 年 AUC      : %s\n    多变量 per SD : %s\n",
                lab, u$conf.int[1, 1], u$conf.int[1, 3], u$conf.int[1, 4], u$coefficients[1, 5],
                1 - pchisq(km$chisq, 1), ifelse(is.na(auc), "（timeROC 未装）", sprintf("%.3f", auc)), mvs))
  }
  res(N, "N 现行：Z 在去重前（192 行）")
  res(O, "O 正确顺序：先去重再 Z（166 行）")
  cat("\n  论文现值：HR 1.46 (1.13–1.89), p=0.004；KM p=0.23；AUC 0.635；多变量 HR 1.45, p=0.005\n")
  cat("  注：多变量里 stage 按原始取值作分类变量；若与论文口径不同，以 N 行是否复现 1.45/0.005 为准判断。\n")
}, "Z 顺序敏感性")

hr("完成"); cat("  把 ", LOG, " 的全部内容贴回 chat；Fig2A.png / Fig2B.png 也发给我看一眼。\n", sep = "")
