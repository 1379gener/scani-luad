suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  # library(survminer)
  library(timeROC)
  library(GSVA)
  library(ggplot2)
})

BASE   <- path.expand("~/caf_spatial_project")
OUTDIR <- file.path(BASE, "results/09_GSE72094_validation")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ── Step 1: 读取probe→gene映射 ─────────────────────────────────
cat("=== Step 1: 读取注释文件 ===\n")
annot_raw <- readLines(file.path(BASE, "data/GPL15048.txt"))
# 找到表格起始行（ID\tGB_LIST...）
header_line <- grep("^ID\tGB_LIST", annot_raw)
cat("Table header at line:", header_line, "\n")

annot <- read.table(
  file.path(BASE, "data/GPL15048.txt"),
  skip       = header_line - 1,
  header     = TRUE,
  sep        = "\t",
  quote      = "",
  fill       = TRUE,
  comment.char = "^",
  stringsAsFactors = FALSE
)
cat(sprintf("Annotation: %d probes, columns: %s\n",
            nrow(annot), paste(colnames(annot), collapse=", ")))

# 构建probe→symbol映射（去重：同基因取第一个probe）
probe2gene <- annot[annot$GeneSymbol != "" & !is.na(annot$GeneSymbol),
                    c("ID","GeneSymbol")]
colnames(probe2gene) <- c("probe_id","symbol")
cat(sprintf("Probes with gene symbol: %d\n", nrow(probe2gene)))

# 检查关键基因
key_genes <- c("TIMP1","PDPN","ACTA2","TAGLN","MYH11","POSTN","COL11A1")
for(g in key_genes) {
  n <- sum(probe2gene$symbol == g)
  cat(sprintf("  %s: %d probes\n", g, n))
}

# ── Step 2: 解析表达矩阵 ──────────────────────────────────────
cat("\n=== Step 2: 解析表达矩阵 ===\n")
fname <- file.path(BASE, "data/GSE72094_series_matrix.txt.gz")
lines <- readLines(gzcon(file(fname, "rb")))
close(file(fname, "rb"))

# 找矩阵边界
tbl_start <- grep("!series_matrix_table_begin", lines) + 1
tbl_end   <- grep("!series_matrix_table_end",   lines) - 1
cat(sprintf("Matrix: rows %d to %d (%d data rows)\n",
            tbl_start, tbl_end, tbl_end - tbl_start))

# 解析header（sample IDs）
header  <- strsplit(lines[tbl_start], "\t")[[1]]
samples <- gsub('"', '', header[-1])
cat(sprintf("Samples: %d\n", length(samples)))

# 只提取我们需要的probes（目标基因对应的probe）
target_genes <- unique(c(
  # SCANI组件
  "TIMP1", "PDPN",
  # myCAF markers（用于ssGSEA）
  probe2gene$symbol[probe2gene$symbol %in% c(
    "ACTA2","TAGLN","MYH11","POSTN","COL11A1","CTHRC1",
    "COL10A1","MMP11","FAP","INHBA","LRRC15","SFRP2",
    "IGFBP7","CALD1","ADAMTS2","DCN","LUM","COL1A1","COL3A1"
  )]
))

target_probes <- probe2gene$probe_id[probe2gene$symbol %in% target_genes]
cat(sprintf("Target probes to extract: %d\n", length(target_probes)))

# 读取表达矩阵（只保留目标probe行）
expr_list <- list()
for(i in (tbl_start+1):tbl_end) {
  parts    <- strsplit(lines[i], "\t")[[1]]
  probe_id <- gsub('"', '', parts[1])
  if(probe_id %in% target_probes) {
    vals <- suppressWarnings(as.numeric(gsub('"','',parts[-1])))
    expr_list[[probe_id]] <- vals
  }
}
cat(sprintf("Extracted %d probe rows\n", length(expr_list)))

# 转为矩阵
expr_mat <- do.call(rbind, expr_list)
colnames(expr_mat) <- samples
cat(sprintf("Expression matrix: %d x %d\n", nrow(expr_mat), ncol(expr_mat)))

# 同基因多probe取均值，转为gene×sample矩阵
probe_ids  <- rownames(expr_mat)
gene_syms  <- probe2gene$symbol[match(probe_ids, probe2gene$probe_id)]
df_expr    <- as.data.frame(expr_mat)
df_expr$symbol <- gene_syms
gene_mat <- df_expr %>%
  group_by(symbol) %>%
  summarise(across(everything(), mean), .groups="drop") %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()
cat(sprintf("Gene matrix: %d genes x %d samples\n",
            nrow(gene_mat), ncol(gene_mat)))

# ── Step 3: 解析临床数据 ──────────────────────────────────────
cat("\n=== Step 3: 解析临床数据 ===\n")
clin_df <- read.csv(file.path(BASE, "data/GSE72094_parsed/clinical.csv"),
                    stringsAsFactors=FALSE)

parse_val <- function(x, field) {
  m <- regmatches(x, regexpr(paste0(field, ": ([^|]+)"), x))
  if(length(m)==0) return(NA)
  sub(paste0(field, ": "), "", m)
}

clin <- data.frame(
  sample_id    = clin_df$sample_id,
  vital_status = clin_df$vital_status,
  OS_days      = as.numeric(clin_df$survival_time_in_days),
  Stage        = clin_df$Stage,
  age          = as.numeric(clin_df$age_at_diagnosis),
  gender       = clin_df$gender,
  stringsAsFactors = FALSE
)
clin$OS_event  <- ifelse(clin$vital_status == "Dead", 1,
                  ifelse(clin$vital_status == "Alive", 0, NA))
clin$OS_months <- clin$OS_days / 30.44
clin <- clin[!is.na(clin$OS_event) & !is.na(clin$OS_months) &
             clin$OS_months > 0, ]
cat(sprintf("Clinical data: %d samples with OS\n", nrow(clin)))
cat(sprintf("Events: %d (%.0f%%)\n",
            sum(clin$OS_event), 100*mean(clin$OS_event)))

# ── Step 4: 计算SCANI ─────────────────────────────────────────
cat("\n=== Step 4: 计算SCANI ===\n")
samples_use <- intersect(clin$sample_id, colnames(gene_mat))
cat(sprintf("Samples with both expr and clinical: %d\n", length(samples_use)))

gene_sub <- gene_mat[, samples_use]
clin_sub <- clin[clin$sample_id %in% samples_use, ]

# myCAF ssGSEA
sig <- read.csv(file.path(BASE, "results/01_scRNA/CAF_signatures.csv"))
mycaf_genes <- sig$Gene[sig$Subtype == "myCAF"]
mycaf_use   <- mycaf_genes[mycaf_genes %in% rownames(gene_sub)]
cat(sprintf("myCAF genes in matrix: %d/%d\n",
            length(mycaf_use), length(mycaf_genes)))

gsva_res <- gsva(gene_sub, list(myCAF=mycaf_use),
                 method="ssgsea", ssgsea.norm=FALSE, verbose=FALSE)
mycaf_scores <- as.numeric(gsva_res["myCAF", samples_use])

# TIMP1 and PDPN
timp1_expr <- as.numeric(gene_sub["TIMP1", samples_use])
pdpn_expr  <- as.numeric(gene_sub["PDPN",  samples_use])

# Z-score
zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
scores_df <- data.frame(
  sample_id = samples_use,
  myCAF_z   = zsc(mycaf_scores),
  TIMP1_z   = zsc(timp1_expr),
  PDPN_z    = zsc(pdpn_expr)
)
scores_df$SCANI <- scores_df$myCAF_z + scores_df$TIMP1_z - scores_df$PDPN_z
cat(sprintf("SCANI range: %.2f to %.2f\n",
            min(scores_df$SCANI), max(scores_df$SCANI)))

# 合并
df_final <- merge(clin_sub, scores_df, by="sample_id")
df_final$SCANI_group <- ifelse(df_final$SCANI >= median(df_final$SCANI),
                                "SCANI_High", "SCANI_Low")
df_final$SCANI_group <- factor(df_final$SCANI_group,
                                levels=c("SCANI_Low","SCANI_High"))

# ── Step 5: KM ────────────────────────────────────────────────
cat("\n=== Step 5: KM分析 ===\n")
fit_km <- survfit(Surv(OS_months, OS_event) ~ SCANI_group, data=df_final)
km_test <- survdiff(Surv(OS_months, OS_event) ~ SCANI_group, data=df_final)
km_p    <- 1 - pchisq(km_test$chisq, df=1)
cat(sprintf("KM log-rank p = %.4f\n", km_p))

p_km <- ggplot(
  data.frame(time=fit_km$time, surv=fit_km$surv,
             group=rep(names(fit_km$strata), fit_km$strata)) %>%
    dplyr::mutate(group=gsub("SCANI_group=","",group)),
  aes(x=time, y=surv, color=group)) +
  geom_step(linewidth=1) +
  scale_color_manual(values=c("SCANI_Low"="#3498DB","SCANI_High"="#E74C3C")) +
  annotate("text", x=max(fit_km$time)*0.6, y=0.95,
           label=sprintf("Log-rank p=%.4f", km_p), size=4) +
  labs(title=sprintf("SCANI Validation GSE72094 (n=%d)", nrow(df_final)),
       x="Time (months)", y="Overall Survival") +
  theme_classic(base_size=12)

ggsave(file.path(OUTDIR,"KM_SCANI_GSE72094.pdf"), p_km, width=7, height=6)
ggsave(file.path(OUTDIR,"KM_SCANI_GSE72094.png"), p_km, width=7, height=6, dpi=300)
cat("KM saved\n")

# ── Step 6: Cox ───────────────────────────────────────────────
cat("\n=== Step 6: Cox回归 ===\n")
cox_uni <- coxph(Surv(OS_months, OS_event) ~ SCANI, data=df_final)
s1 <- summary(cox_uni)
cat(sprintf("Univariable: HR=%.3f (95%%CI: %.3f-%.3f), p=%.4f\n",
            s1$coefficients[,"exp(coef)"],
            s1$conf.int[,"lower .95"],
            s1$conf.int[,"upper .95"],
            s1$coefficients[,"Pr(>|z|)"]))

df_mv <- df_final[!is.na(df_final$age) & !is.na(df_final$Stage), ]
df_mv$Stage_num <- suppressWarnings(
  as.numeric(gsub("[A-Za-z]","", df_mv$Stage)))
df_mv <- df_mv[!is.na(df_mv$Stage_num), ]
cat(sprintf("Multivariable cohort: %d samples\n", nrow(df_mv)))

cox_mv <- coxph(Surv(OS_months, OS_event) ~ SCANI + Stage_num + age,
                data=df_mv)
s2 <- summary(cox_mv)
cat(sprintf("Multivariable: HR=%.3f (95%%CI: %.3f-%.3f), p=%.4f\n",
            s2$coefficients["SCANI","exp(coef)"],
            s2$conf.int["SCANI","lower .95"],
            s2$conf.int["SCANI","upper .95"],
            s2$coefficients["SCANI","Pr(>|z|)"]))

# ── Step 7: AUC ───────────────────────────────────────────────
cat("\n=== Step 7: 时间依赖AUC ===\n")
roc <- timeROC(T=df_final$OS_months, delta=df_final$OS_event,
               marker=df_final$SCANI, cause=1,
               times=c(12,36,60), iid=TRUE)
cat(sprintf("AUC 1yr=%.3f, 3yr=%.3f, 5yr=%.3f\n",
            roc$AUC[1], roc$AUC[2], roc$AUC[3]))

# ── Step 8: 保存 ──────────────────────────────────────────────
write.csv(df_final, file.path(OUTDIR,"SCANI_GSE72094_scores.csv"),
          row.names=FALSE)

summary_df <- data.frame(
  Cohort  = "GSE72094", N=nrow(df_final),
  Events  = sum(df_final$OS_event),
  KM_p    = round(km_p,4),
  HR_uni  = round(s1$coefficients[,"exp(coef)"],3),
  p_uni   = round(s1$coefficients[,"Pr(>|z|)"],4),
  HR_mv   = round(s2$coefficients["SCANI","exp(coef)"],3),
  p_mv    = round(s2$coefficients["SCANI","Pr(>|z|)"],4),
  AUC_3yr = round(roc$AUC[2],3)
)
write.csv(summary_df, file.path(OUTDIR,"validation_summary.csv"),
          row.names=FALSE)

cat("\n=== 完成！===\n")
print(summary_df)
cat("输出目录:", OUTDIR, "\n")
