#!/usr/bin/env Rscript
# ==============================================================================
# sci1_figs_final.R —— SCI-1 投稿用合成图 Figure 2 / Figure 4（IJMS：一图一号）
#
#   Figure 2 = A 各来源构成（横向 100% 堆叠）+ B 签名基因点图
#   Figure 4 = A KM（加删失标记 + number-at-risk 表，横轴显示到 120 个月）+ B 3 年 ROC
#
#   只读：tables/fig2_counts_origin_by_subtype.csv、tables/SCANI_TCGA_O_166.csv、
#         caf_spatial_project/results/01_scRNA/{seu_CAF_annotated.rds, CAF_subtype_signatures.rds}
#   只写：~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/figures/Figure2.* 与 Figure4.*
#   脚本内无字面量数据向量；统计量全部现算。
#
# 用法：Rscript ~/sci1_figs_final.R
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(survival); library(cowplot) })
H    <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS")
SUB_LEVELS <- c("Myofibroblasts (myCAF)", "Matrix FBs (iCAF-like)", "FB-like cells")
ORI_LEVELS <- c("Primary tumour", "Normal lung", "Lymph node metastasis", "Brain metastasis")
COLS <- setNames(c("#B2182B", "#4393C3", "#7F7F7F"), SUB_LEVELS)
XMAX <- 120; BRK <- seq(0, XMAX, 24)          # KM 横轴显示范围（月），仅影响显示，不影响检验

save3 <- function(p, nm, w, h) {
  f <- function(ext) file.path(OUTD, "figures", paste0(nm, ext))
  ggsave(f(".tiff"), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(f(".png"),  p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(f(".pdf"),  p, width = w, height = h, bg = "white")
}

# ============================================================== Figure 2A
tab <- read.csv(file.path(OUTD, "tables/fig2_counts_origin_by_subtype.csv"), stringsAsFactors = FALSE)
stopifnot(sum(tab$Freq) == 2992)
tab$Origin <- factor(tab$Origin, levels = rev(ORI_LEVELS))       # 横向：原发瘤在最上
tab$Group  <- factor(tab$Group,  levels = rev(SUB_LEVELS))       # myCAF 在最左
tot <- aggregate(Freq ~ Origin, tab, sum)
tab <- merge(tab, setNames(tot, c("Origin", "Total")), by = "Origin")
tab$Pct <- 100 * tab$Freq / tab$Total

p2A <- ggplot(tab, aes(Pct, Origin, fill = Group)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(data = tot, aes(101, Origin, label = paste0("n = ", Freq)),
            inherit.aes = FALSE, hjust = 0, size = 3.2) +
  scale_fill_manual(values = COLS, breaks = SUB_LEVELS, name = NULL) +
  scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) sub("Lymph node metastasis", "Lymph node\nmetastasis", x)) +
  labs(x = "Fibroblasts within tissue of origin (%)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(margin = margin(l = 4, r = 14)),
        plot.margin = margin(5, 10, 5, 5))

# ============================================================== Figure 2B
if (Sys.getenv("SCI1_SKIP_2B") == "1") {                          # 仅本地无 Seurat 时测试排版
  p2B <- ggplot() + annotate("text", 0, 0, label = "2B placeholder") + theme_void()
} else {
  suppressPackageStartupMessages(library(Seurat))
  obj <- readRDS(file.path(C, "results/01_scRNA/seu_CAF_annotated.rds"))
  stopifnot(ncol(obj) == 2992)
  MAP <- c("Myofibroblasts" = SUB_LEVELS[1], "COL13A1+ matrix FBs" = SUB_LEVELS[2],
           "COL14A1+ matrix FBs" = SUB_LEVELS[2], "FB-like cells" = SUB_LEVELS[3])
  obj$FB_group <- factor(unname(MAP[obj@meta.data$Cell_subtype]), levels = rev(SUB_LEVELS))
  stopifnot(!anyNA(obj$FB_group))
  sigs  <- readRDS(file.path(C, "results/01_scRNA/CAF_subtype_signatures.rds"))
  feats <- list(myCAF = head(sigs$myCAF, 5), `Matrix FBs` = head(sigs$iCAF, 5), `FB-like` = head(sigs$apCAF, 5))
  feats <- lapply(feats, intersect, rownames(obj))
  gset  <- data.frame(gene = unlist(feats), set = rep(names(feats), lengths(feats)))
  gset  <- gset[!duplicated(gset$gene), ]
  dd <- DotPlot(obj, features = gset$gene, group.by = "FB_group")$data
  dd$set <- factor(gset$set[match(as.character(dd$features.plot), gset$gene)], levels = names(feats))
  dd$id  <- factor(dd$id, levels = rev(SUB_LEVELS))
  p2B <- ggplot(dd, aes(features.plot, id)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_colour_gradient(low = "lightgrey", high = "#B2182B", name = "Scaled average\nexpression") +
    scale_size_area(max_size = 6, name = "Percent\nexpressed") +
    facet_grid(cols = vars(set), scales = "free_x", space = "free_x") +
    scale_y_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
          panel.grid = element_blank(), strip.background = element_blank(),
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.title = element_text(size = 9), plot.margin = margin(5, 10, 5, 25))
}
fig2 <- plot_grid(p2A, p2B, ncol = 1, labels = c("A", "B"), label_size = 14, rel_heights = c(0.75, 1),
                  align = "v", axis = "l")
save3(fig2, "Figure2", 9.5, 7)
cat("Figure 2 已写出\n")

# ============================================================== Figure 4 数据
d <- read.csv(file.path(OUTD, "tables/SCANI_TCGA_O_166.csv"), stringsAsFactors = FALSE)
d <- d[!duplicated(d$patient_id), ]
stopifnot(nrow(d) == 166, sum(d$OS_event) == 80)
d$grp <- factor(ifelse(d$SCANI >= median(d$SCANI), "High", "Low"), levels = c("Low", "High"))
km   <- survdiff(Surv(OS_time, OS_event) ~ grp, d); km_p <- 1 - pchisq(km$chisq, 1)
fit  <- survfit(Surv(OS_time / 30.44, OS_event) ~ grp, d)
ns   <- table(d$grp)
labs_g <- setNames(paste0(names(ns), " (n=", ns, ")"), names(ns))
gcol   <- setNames(c("#4575B4", "#D73027"), labs_g)

# ---- KM 曲线 + 删失标记
kmdf <- do.call(rbind, lapply(seq_along(fit$strata), function(i) {
  idx <- if (i == 1) seq_len(fit$strata[1]) else sum(fit$strata[1:(i - 1)]) + seq_len(fit$strata[i])
  g <- sub("grp=", "", names(fit$strata)[i])
  data.frame(time = c(0, fit$time[idx]), surv = c(1, fit$surv[idx]),
             cens = c(0, fit$n.censor[idx]), group = unname(labs_g[g]))
}))
kmdf$group <- factor(kmdf$group, levels = labs_g)
p4A <- ggplot(kmdf, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.8) +
  geom_point(data = subset(kmdf, cens > 0), shape = 3, size = 1.6, stroke = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = gcol, name = "SCANI") +
  scale_x_continuous(breaks = BRK) + scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, XMAX)) +
  annotate("text", x = XMAX * 0.7, y = 0.95, size = 3.3, label = sprintf("Log-rank p = %.2f", km_p)) +
  labs(x = NULL, y = "Overall survival probability") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(margin = margin(l = 4, r = 14)))

# ---- number at risk
s  <- summary(fit, times = BRK, extend = TRUE)
rt <- data.frame(time = s$time, n = s$n.risk, group = unname(labs_g[sub("grp=", "", as.character(s$strata))]))
rt$group <- factor(rt$group, levels = rev(labs_g))
cat("\nNumber at risk:\n"); print(xtabs(n ~ group + time, rt))
cat("最长随访（月）：", round(max(d$OS_time) / 30.44, 1), "；", XMAX, " 个月后仍在随访的病人数：",
    sum(d$OS_time / 30.44 > XMAX), "\n", sep = "")
pRisk <- ggplot(rt, aes(time, group, label = n, colour = group)) +
  geom_text(size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = gcol) +
  scale_x_continuous(breaks = BRK) + coord_cartesian(xlim = c(0, XMAX)) +
  scale_y_discrete(labels = function(x) sub(" \\(.*", "", x)) +
  labs(x = "Time (months)", y = NULL, title = "Number at risk") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 10), axis.line.y = element_blank(), axis.ticks.y = element_blank())
km_block <- plot_grid(p4A, pRisk, ncol = 1, rel_heights = c(4, 1.25), align = "v", axis = "lr")

# ---- 3 年 ROC
roc <- timeROC::timeROC(T = d$OS_time, delta = d$OS_event, marker = d$SCANI,
                        cause = 1, times = 3 * 365, iid = FALSE)
auc <- unname(roc$AUC[2])
rdf <- data.frame(FPR = roc$FP[, 2], TPR = roc$TP[, 2]); rdf <- rdf[order(rdf$FPR, rdf$TPR), ]
p4B <- ggplot(rdf, aes(FPR, TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_step(colour = "#B2182B", linewidth = 0.8) +
  annotate("text", x = 0.65, y = 0.12, size = 3.3, label = sprintf("3-year AUC = %.3f", auc)) +
  coord_equal() + labs(x = "1 - specificity", y = "Sensitivity") +
  theme_classic(base_size = 11)

fig4 <- plot_grid(km_block, p4B, nrow = 1, labels = c("A", "B"), label_size = 14, rel_widths = c(1.2, 1))
save3(fig4, "Figure4", 10, 5)
cat(sprintf("\n核对：log-rank p = %.3f；3 年 AUC = %.3f\n", km_p, auc))
cat("Figure 4 已写出。请发回 Figure2.png、Figure4.png 与上面的 number-at-risk 表。\n")
