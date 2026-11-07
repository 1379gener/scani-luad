# =============================================================================
# step3_spatial_GSE307534.R
# 空间转录组 CAF 生态位映射 - GSE307534（56个Visium样本）
# 数据: GSE307534 (Normal/AAH/AIS/MIA/LUAD 各阶段)
# 输入: CAF_subtype_signatures.rds (来自 step2)
# 运行: nohup Rscript step3_spatial_GSE307534.R > logs/step3_GSE307534.log 2>&1 &
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
})

set.seed(42)

DATA_DIR   <- path.expand("~/caf_spatial_project/data/GSE307534")
SC_RES     <- path.expand("~/caf_spatial_project/results/01_scRNA")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/03_spatial_GSE307534")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/03_spatial_GSE307534")
dir.create(RESULT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG_DIR,    recursive=TRUE, showWarnings=FALSE)

plot_spatial <- function(seu, color_by, title="", discrete=FALSE, pt_size=0.3) {
  df <- data.frame(
    x     = seu$pxl_col,
    y     = -seu$pxl_row,
    value = seu@meta.data[[color_by]]
  )
  p <- ggplot(df, aes(x=x, y=y, color=value)) +
    geom_point(size=pt_size) + coord_fixed() +
    theme_void() + ggtitle(title) +
    theme(legend.title=element_text(size=8),
          legend.text=element_text(size=7),
          plot.title=element_text(size=9))
  if (!discrete) {
    p <- p + scale_color_gradientn(
      colors=c("lightgrey","yellow","red"), name="Score")
  }
  return(p)
}

load_visium_sample <- function(sample_dir, sample_name) {
  mtx_dir    <- file.path(sample_dir, "filtered_feature_bc_matrix")
  coord_file <- file.path(sample_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(coord_file)) {
    coord_file <- file.path(sample_dir, "spatial", "tissue_positions_list.csv")
  }

  raw <- Read10X(data.dir=mtx_dir, gene.column=2)
  counts <- if (is.list(raw)) raw[["Gene Expression"]] else raw

  coords <- read.csv(coord_file, header=TRUE, stringsAsFactors=FALSE)
  if ("pxl_row_in_fullres" %in% colnames(coords)) {
    coords <- coords %>%
      rename(pxl_row = pxl_row_in_fullres,
             pxl_col = pxl_col_in_fullres)
  }
  rownames(coords) <- coords$barcode
  coords_tissue <- coords[coords$in_tissue == 1, ]

  common <- intersect(colnames(counts), rownames(coords_tissue))
  if (length(common) == 0) return(NULL)

  counts        <- counts[, common]
  coords_tissue <- coords_tissue[common, ]

  seu <- CreateSeuratObject(
    counts  = counts,
    assay   = "Spatial",
    project = sample_name
  )

  seu$array_row <- coords_tissue[colnames(seu), "array_row"]
  seu$array_col <- coords_tissue[colnames(seu), "array_col"]
  seu$pxl_row   <- coords_tissue[colnames(seu), "pxl_row"]
  seu$pxl_col   <- coords_tissue[colnames(seu), "pxl_col"]
  seu$sample_id <- sample_name

  stage <- gsub("P[0-9]+_", "", sample_name)
  stage <- gsub("[0-9]+$", "", stage)
  seu$stage <- stage

  return(seu)
}

# =============================================================================
# STEP 1: 读取 CAF Signatures
# =============================================================================
cat("\n[STEP 1] 读取 CAF Signatures...\n")
CAF_SIGS <- readRDS(file.path(SC_RES, "CAF_subtype_signatures.rds"))
for (nm in names(CAF_SIGS)) {
  cat(sprintf("  %s: %d 个基因\n", nm, length(CAF_SIGS[[nm]])))
}

# =============================================================================
# STEP 2: 读取所有56个样本
# =============================================================================
cat("\n[STEP 2] 读取 GSE307534 样本...\n")

sample_dirs  <- list.dirs(DATA_DIR, recursive=FALSE)
sample_dirs  <- sample_dirs[!grepl("\\.tar", sample_dirs)]
sample_names <- basename(sample_dirs)

cat(sprintf("  发现 %d 个样本目录\n", length(sample_dirs)))

seurat_list <- list()
for (i in seq_along(sample_dirs)) {
  sname <- sample_names[i]
  cat(sprintf("  [%d/%d] 读取: %s\n", i, length(sample_dirs), sname))
  seu <- tryCatch(
    load_visium_sample(sample_dirs[i], sname),
    error = function(e) {
      cat(sprintf("    ⚠ 失败: %s\n", e$message))
      return(NULL)
    }
  )
  if (!is.null(seu)) seurat_list[[sname]] <- seu
}

cat(sprintf("\n  成功读取: %d / %d 个样本\n",
            length(seurat_list), length(sample_dirs)))
print(table(sapply(seurat_list, function(x) unique(x$stage))))

# =============================================================================
# STEP 3: 逐样本 CAF Signature 打分
# =============================================================================
cat("\n[STEP 3] CAF Signature 打分（逐样本）...\n")

score_results <- list()

for (sname in names(seurat_list)) {
  seu <- seurat_list[[sname]]

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-",
                                               assay="Spatial")
  keep <- colnames(seu)[
    seu$nCount_Spatial   > 100 &
    seu$nFeature_Spatial > 50  &
    seu$percent.mt       < 30
  ]
  if (length(keep) < 10) next
  seu <- seu[, keep]
  seu <- NormalizeData(seu, assay="Spatial", verbose=FALSE)

  for (subtype in names(CAF_SIGS)) {
    valid <- intersect(CAF_SIGS[[subtype]], rownames(seu))
    if (length(valid) < 5) next
    seu <- AddModuleScore(
      seu,
      features = list(valid),
      name     = paste0(subtype, "_score"),
      nbin     = 10,
      assay    = "Spatial",
      seed     = 42
    )
  }

  seurat_list[[sname]] <- seu

  score_row <- data.frame(
    sample_id  = sname,
    stage      = unique(seu$stage),
    n_spots    = ncol(seu),
    myCAF_mean = ifelse("myCAF_score1" %in% colnames(seu@meta.data),
                        mean(seu$myCAF_score1), NA),
    iCAF_mean  = ifelse("iCAF_score1"  %in% colnames(seu@meta.data),
                        mean(seu$iCAF_score1),  NA),
    apCAF_mean = ifelse("apCAF_score1" %in% colnames(seu@meta.data),
                        mean(seu$apCAF_score1), NA)
  )
  score_results[[sname]] <- score_row
  cat(sprintf("  ✓ %s: %d spots, myCAF=%.3f, iCAF=%.3f\n",
              sname, ncol(seu),
              ifelse(is.na(score_row$myCAF_mean), 0, score_row$myCAF_mean),
              ifelse(is.na(score_row$iCAF_mean),  0, score_row$iCAF_mean)))
}

score_df <- do.call(rbind, score_results)
write.csv(score_df, file.path(RESULT_DIR, "CAF_scores_per_sample.csv"),
          row.names=FALSE)

# =============================================================================
# STEP 4: 可视化 CAF 分数随肿瘤进展变化
# =============================================================================
cat("\n[STEP 4] 可视化 CAF 分数随肿瘤进展变化...\n")

stage_order <- c("Normal","AAH","AIS","MIA","LUAD")
score_df$stage <- factor(score_df$stage, levels=stage_order)

pdf(file.path(FIG_DIR, "CAF_scores_by_stage.pdf"), width=14, height=5)
p1 <- ggplot(score_df, aes(x=stage, y=myCAF_mean, fill=stage)) +
  geom_boxplot(outlier.size=0.5) +
  geom_jitter(width=0.15, size=1.5, alpha=0.6) +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  theme_classic() +
  labs(title="myCAF Score across LUAD Progression",
       x="Stage", y="Mean myCAF Score") +
  theme(legend.position="none")

p2 <- ggplot(score_df, aes(x=stage, y=iCAF_mean, fill=stage)) +
  geom_boxplot(outlier.size=0.5) +
  geom_jitter(width=0.15, size=1.5, alpha=0.6) +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  theme_classic() +
  labs(title="iCAF Score across LUAD Progression",
       x="Stage", y="Mean iCAF Score") +
  theme(legend.position="none")

p3 <- ggplot(score_df, aes(x=stage, y=apCAF_mean, fill=stage)) +
  geom_boxplot(outlier.size=0.5) +
  geom_jitter(width=0.15, size=1.5, alpha=0.6) +
  scale_fill_brewer(palette="RdYlBu", direction=-1) +
  theme_classic() +
  labs(title="apCAF Score across LUAD Progression",
       x="Stage", y="Mean apCAF Score") +
  theme(legend.position="none")

print(p1 + p2 + p3)
dev.off()

# =============================================================================
# STEP 5: 定义 CAF Niche（仅 LUAD 样本）
# =============================================================================
cat("\n[STEP 5] 定义 CAF Niche（LUAD 样本）...\n")

luad_samples <- names(seurat_list)[
  sapply(seurat_list, function(x) unique(x$stage) == "LUAD")
]
cat(sprintf("  LUAD 样本数: %d\n", length(luad_samples)))

niche_markers_list <- list()

for (sname in luad_samples) {
  seu <- seurat_list[[sname]]
  if (!all(c("myCAF_score1","iCAF_score1") %in% colnames(seu@meta.data))) next

  for (q in c(0.75, 0.65, 0.60)) {
    thr_my <- quantile(seu$myCAF_score1, q)
    thr_ic <- quantile(seu$iCAF_score1,  q)
    seu$CAF_niche <- "Other"
    seu$CAF_niche[seu$myCAF_score1 >= thr_my] <- "myCAF_high"
    seu$CAF_niche[seu$iCAF_score1  >= thr_ic] <- "iCAF_high"
    seu$CAF_niche[
      seu$myCAF_score1 >= thr_my &
      seu$iCAF_score1  >= thr_ic
    ] <- "CAF_Niche"
    if (sum(seu$CAF_niche == "CAF_Niche") >= 10) break
  }

  n_niche <- sum(seu$CAF_niche == "CAF_Niche")
  if (n_niche < 10) next
  seurat_list[[sname]] <- seu

  Idents(seu) <- "CAF_niche"
  markers <- tryCatch(
    FindMarkers(seu, ident.1="CAF_Niche", ident.2="Other",
                test.use="wilcox", min.pct=0.1,
                logfc.threshold=0.25, only.pos=TRUE,
                assay="Spatial", verbose=FALSE),
    error = function(e) NULL
  )
  if (!is.null(markers) && nrow(markers) > 0) {
    markers$gene   <- rownames(markers)
    markers$sample <- sname
    niche_markers_list[[sname]] <- markers
    cat(sprintf("  ✓ %s: %d Niche spots, %d DEGs\n",
                sname, n_niche, nrow(markers)))
  }
}

# =============================================================================
# STEP 6: 整合核心基因集
# =============================================================================
cat("\n[STEP 6] 整合 Niche 核心基因集...\n")

all_markers <- do.call(rbind, niche_markers_list)
write.csv(all_markers,
          file.path(RESULT_DIR, "CAF_Niche_DEGs_all_LUAD.csv"),
          row.names=FALSE)

min_samples <- max(2, floor(length(niche_markers_list) / 3))
gene_freq <- all_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  group_by(gene) %>%
  summarise(
    n_samples  = n(),
    mean_logFC = mean(avg_log2FC),
    mean_pct1  = mean(pct.1),
    .groups    = "drop"
  ) %>%
  filter(n_samples >= min_samples) %>%
  arrange(desc(mean_logFC))

cat(sprintf("  在 >= %d 个样本中重复的基因: %d 个\n",
            min_samples, nrow(gene_freq)))

core_genes <- head(gene_freq$gene, 50)
cat(sprintf("  核心基因集: %d 个\n", length(core_genes)))
cat(sprintf("  Top 10: %s\n", paste(head(core_genes,10), collapse=", ")))

write.csv(gene_freq,
          file.path(RESULT_DIR, "CAF_Niche_gene_frequency.csv"),
          row.names=FALSE)
saveRDS(core_genes,
        file.path(RESULT_DIR, "CAF_Niche_final_geneset.rds"))
writeLines(core_genes,
           file.path(RESULT_DIR, "CAF_Niche_final_geneset.txt"))

# =============================================================================
# STEP 7: 空间可视化
# =============================================================================
cat("\n[STEP 7] 空间可视化...\n")

rep_samples <- c("P4_Normal","P1_AAH","P8_AIS","P10_MIA","P1_LUAD")
rep_samples <- rep_samples[rep_samples %in% names(seurat_list)]

pdf(file.path(FIG_DIR, "ST_myCAF_score_by_stage.pdf"), width=20, height=4)
plots <- lapply(rep_samples, function(sname) {
  seu <- seurat_list[[sname]]
  if (!"myCAF_score1" %in% colnames(seu@meta.data)) return(NULL)
  stage <- unique(seu$stage)
  plot_spatial(seu, "myCAF_score1",
               title=paste0(stage, " (", sname, ")"))
})
plots <- Filter(Negate(is.null), plots)
if (length(plots) > 0) print(wrap_plots(plots, nrow=1))
dev.off()

NICHE_COLS <- c("CAF_Niche"="#d62728","myCAF_high"="#ff7f0e",
                "iCAF_high"="#1f77b4","Other"="#e8e8e8")

luad_with_niche <- luad_samples[
  sapply(luad_samples, function(s)
    "CAF_niche" %in% colnames(seurat_list[[s]]@meta.data))
]

if (length(luad_with_niche) > 0) {
  n_show <- min(6, length(luad_with_niche))
  pdf(file.path(FIG_DIR, "ST_CAF_Niche_LUAD_samples.pdf"),
      width=18, height=ceiling(n_show/3)*6)
  plots <- lapply(luad_with_niche[1:n_show], function(sname) {
    seu <- seurat_list[[sname]]
    df  <- data.frame(
      x     = seu$pxl_col,
      y     = -seu$pxl_row,
      niche = seu$CAF_niche
    )
    ggplot(df, aes(x=x, y=y, color=niche)) +
      geom_point(size=0.3) + coord_fixed() +
      scale_color_manual(values=NICHE_COLS) +
      theme_void() + ggtitle(sname) +
      theme(plot.title=element_text(size=9),
            legend.text=element_text(size=7))
  })
  print(wrap_plots(plots, ncol=3))
  dev.off()
}

# =============================================================================
# STEP 8: 保存
# =============================================================================
cat("\n[STEP 8] 保存...\n")
saveRDS(score_df,     file.path(RESULT_DIR, "CAF_score_summary.rds"))
saveRDS(seurat_list,  file.path(RESULT_DIR, "seurat_list_GSE307534.rds"))

cat("\n✅ Step 3 (GSE307534) 完成！\n")
cat(sprintf("  → %s/CAF_Niche_final_geneset.txt\n", RESULT_DIR))
cat(sprintf("  → %s/CAF_scores_per_sample.csv\n",   RESULT_DIR))
cat(sprintf("  → %s/\n",                            FIG_DIR))

# 补充：Top 6 Niche 基因空间表达图（代表性LUAD样本）
cat("\n[补充] Top 6 基因空间表达图...\n")

plot_spatial_gene <- function(seu, gene, title=gene) {
  expr <- FetchData(seu, vars=gene, layer="data")[, 1]
  df <- data.frame(
    x     = seu$pxl_col,
    y     = -seu$pxl_row,
    value = expr
  )
  ggplot(df, aes(x=x, y=y, color=value)) +
    geom_point(size=0.3) + coord_fixed() + theme_void() +
    scale_color_gradientn(colors=c("lightgrey","yellow","red"), name="Expr") +
    ggtitle(title) +
    theme(legend.title=element_text(size=8), plot.title=element_text(size=9))
}

core_genes <- readLines(
  "~/caf_spatial_project/results/03_spatial_GSE307534/CAF_Niche_final_geneset.txt"
)
seurat_list <- readRDS(
  "~/caf_spatial_project/results/03_spatial_GSE307534/seurat_list_GSE307534.rds"
)

# 用 P5_LUAD（myCAF分数最高的样本）
rep_luad <- "P5_LUAD"
seu_rep  <- seurat_list[[rep_luad]]
top6     <- intersect(head(core_genes, 6), rownames(seu_rep))

pdf("~/caf_spatial_project/figures/03_spatial_GSE307534/ST_Niche_top_genes.pdf",
    width=18, height=12)
plots <- lapply(top6, function(g) plot_spatial_gene(seu_rep, g))
print(patchwork::wrap_plots(plots, ncol=3))
dev.off()
cat("  ✓ ST_Niche_top_genes.pdf 生成完成\n")
