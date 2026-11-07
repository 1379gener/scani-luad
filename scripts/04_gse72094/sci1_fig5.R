#!/usr/bin/env Rscript
# ==============================================================================
# sci1_fig5.R —— SCI-1 Figure 5：GSE72094 外部验证 KM（全基因组版 SCANI；删失标记 + number at risk）
#   只读 tables/SCANI_GSE72094_both_versions.csv（列 SCANI_full）；脚本内无字面量数据
#   只写 ~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/figures/Figure5.*
# 用法：Rscript ~/sci1_fig5.R
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(survival); library(cowplot) })
H <- path.expand("~"); C <- file.path(H, "caf_spatial_project")
OUTD <- file.path(H, "thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS")

# 主分析：全基因组版 SCANI（sci1_recheck2.R 写出；50/50 签名基因、全基因组排序）
d <- read.csv(file.path(OUTD, "tables/SCANI_GSE72094_both_versions.csv"), stringsAsFactors = FALSE)
d <- d[!duplicated(d$sample_id), ]
d$SCANI <- d$SCANI_full
d <- d[!is.na(d$OS_days) & d$OS_days > 0 & !is.na(d$OS_event) & !is.na(d$SCANI), ]
stopifnot(nrow(d) == 398, sum(d$OS_event) == 113)

d$grp <- factor(ifelse(d$SCANI >= median(d$SCANI), "High", "Low"), levels = c("Low", "High"))
km_p  <- 1 - pchisq(survdiff(Surv(OS_days, OS_event) ~ grp, d)$chisq, 1)
p_gt  <- 1 - pchisq(survdiff(Surv(OS_days, OS_event) ~ I(SCANI > median(SCANI)), d)$chisq, 1)
u     <- summary(coxph(Surv(OS_days, OS_event) ~ scale(SCANI), d))
cat(sprintf("n = %d，事件 %d；分组 %s\n", nrow(d), sum(d$OS_event), paste(names(table(d$grp)), table(d$grp), collapse = " / ")))
cat(sprintf("log-rank p（>= 中位数）= %.4f；（> 中位数）= %.4f；预期 0.0038（全基因组版）\n", km_p, p_gt))
cat(sprintf("单变量 per SD：HR %.2f (%.2f–%.2f), p = %.3f\n", u$conf.int[1, 1], u$conf.int[1, 3], u$conf.int[1, 4], u$coefficients[1, 5]))

d$mo <- d$OS_days / 30.44
XMAX <- 12 * ceiling(max(d$mo) / 12); BRK <- seq(0, XMAX, 12)
fit  <- survfit(Surv(mo, OS_event) ~ grp, d)
ns   <- table(d$grp)
labs_g <- setNames(paste0(names(ns), " (n=", ns, ")"), names(ns))
gcol   <- setNames(c("#4575B4", "#D73027"), labs_g)

kmdf <- do.call(rbind, lapply(seq_along(fit$strata), function(i) {
  idx <- if (i == 1) seq_len(fit$strata[1]) else sum(fit$strata[1:(i - 1)]) + seq_len(fit$strata[i])
  g <- sub("grp=", "", names(fit$strata)[i])
  data.frame(time = c(0, fit$time[idx]), surv = c(1, fit$surv[idx]),
             cens = c(0, fit$n.censor[idx]), group = unname(labs_g[g]))
}))
kmdf$group <- factor(kmdf$group, levels = labs_g)
pKM <- ggplot(kmdf, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.8) +
  geom_point(data = subset(kmdf, cens > 0), shape = 3, size = 1.4, stroke = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = gcol, name = "SCANI") +
  scale_x_continuous(breaks = BRK) + scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, XMAX)) +
  annotate("text", x = XMAX * 0.72, y = 0.95, size = 3.3, label = sprintf("Log-rank p = %.3f", km_p)) +
  labs(x = NULL, y = "Overall survival probability") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(margin = margin(l = 4, r = 14)))

s  <- summary(fit, times = BRK, extend = TRUE)
rt <- data.frame(time = s$time, n = s$n.risk, group = unname(labs_g[sub("grp=", "", as.character(s$strata))]))
rt$group <- factor(rt$group, levels = rev(labs_g))
cat("\nNumber at risk:\n"); print(xtabs(n ~ group + time, rt))
pRisk <- ggplot(rt, aes(time, group, label = n, colour = group)) +
  geom_text(size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = gcol) +
  scale_x_continuous(breaks = BRK) + coord_cartesian(xlim = c(0, XMAX)) +
  scale_y_discrete(labels = function(x) sub(" \\(.*", "", x)) +
  labs(x = "Time (months)", y = NULL, title = "Number at risk") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 10), axis.line.y = element_blank(), axis.ticks.y = element_blank())

fig5 <- plot_grid(pKM, pRisk, ncol = 1, rel_heights = c(4, 1.25), align = "v", axis = "lr")
for (ext in c(".tiff", ".png", ".pdf")) {
  f <- file.path(OUTD, "figures", paste0("Figure5", ext))
  if (ext == ".tiff") ggsave(f, fig5, width = 6, height = 5.5, dpi = 600, compression = "lzw", bg = "white")
  else if (ext == ".png") ggsave(f, fig5, width = 6, height = 5.5, dpi = 600, bg = "white")
  else ggsave(f, fig5, width = 6, height = 5.5, bg = "white")
}
cat("\nFigure 5 已写出。贴回以上输出，并发 Figure5.png。\n")
