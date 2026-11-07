#!/usr/bin/env Rscript
# ==============================================================================
# sci1_fig4_p3.R —— 重画 Figure 4（TCGA-CDR 226），唯一改动：KM 图上 log-rank P 显示 3 位小数（0.016），与正文、图注一致
#   画图代码与 sci1_cdr226_main.R 第 8 节逐行相同；输入只读 tables/cdr226/SCANI_TCGA_CDR226.csv
#   只写（已存在即停止）：figures/Figure4_cdr226_p3.{tiff,png,pdf}
#   用法：Rscript ~/sci1_fig4_p3.R
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(ggplot2); library(cowplot) })
H <- path.expand("~"); B <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS"); FIGN <- "Figure4_cdr226_p3"
for (ext in c(".tiff", ".png", ".pdf")) if (file.exists(file.path(B, "figures", paste0(FIGN, ext)))) stop("图已存在，不覆盖：", FIGN, ext)
D3 <- 3 * 365
kmp <- function(t, e, x) { g <- x >= median(x); 1 - pchisq(survdiff(Surv(t, e) ~ g)$chisq, 1) }
d <- read.csv(file.path(B, "tables/cdr226/SCANI_TCGA_CDR226.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(d) == 226, !anyDuplicated(d$patient_id), sum(d$OS_event) == 115)
d$t <- d$OS_time; d$e <- d$OS_event

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
  cat("  → 写出 ", ff(".png"), "\n", sep = "") 
stopifnot(abs(km_p - 0.016) < 0.0005, abs(auc - 0.621) < 0.0005)
cat("md5：", tools::md5sum(ff(".png")), "\n完成。把这几行输出贴回，并把 ", ff(".png"), " 下载发给我。\n", sep = "")
