# =============================================================================
# step2_scrna_caf.R
# scRNA-seq CAF 亚群鉴定
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

DATA_DIR   <- path.expand("~/caf_spatial_project/data/GSE131907")
RESULT_DIR <- path.expand("~/caf_spatial_project/results/01_scRNA")
FIG_DIR    <- path.expand("~/caf_spatial_project/figures/01_scRNA")
dir.create(RESULT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(FIG_DIR,    recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# STEP 1: 读取数据
# =============================================================================
cat("\n[STEP 1] 读取数据 (约3-5分钟)...\n")

mat <- readRDS(file.path(DATA_DIR,
        "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds"))
meta <- read.table(
  file.path(DATA_DIR, "GSE131907_Lung_Cancer_cell_annotation.txt"),
  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE
)

cat(sprintf("  矩阵: %d genes x %d cells\n", nrow(mat), ncol(mat)))
cat(sprintf("  注释: %d cells\n", nrow(meta)))

common <- intersect(colnames(mat), rownames(meta))
cat(sprintf("  对齐后: %d cells\n", length(common)))
mat  <- mat[, common]
meta <- meta[common, ]

# =============================================================================
# STEP 2: 创建 Seurat 对象
# =============================================================================
cat("\n[STEP 2] 创建 Seurat 对象...\n")

seu <- CreateSeuratObject(
  counts    = mat,
  meta.data = meta,
  project   = "GSE131907_LUAD"
)
seu <- SetAssayData(seu, slot="data",
                    new.data=as(as.matrix(mat), "dgCMatrix"))

cat(sprintf("  Seurat 对象: %d cells\n", ncol(seu)))

# =============================================================================
# STEP 3: 降维
# =============================================================================
cat("\n[STEP 3] 降维聚类...\n")

seu <- FindVariableFeatures(seu, selection.method="vst",
                             nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE)
seu <- RunPCA(seu, npcs=50, verbose=FALSE)

pdf(file.path(FIG_DIR, "PCA_elbow.pdf"), width=8, height=5)
ElbowPlot(seu, ndims=50)
dev.off()

cat("  运行 Harmony...\n")
seu <- RunHarmony(seu, group.by.vars="Sample",
                  reduction="pca", max.iter.harmony=30)

seu <- RunUMAP(seu, reduction="harmony", dims=1:30, verbose=FALSE)
seu <- FindNeighbors(seu, reduction="harmony", dims=1:30, verbose=FALSE)
seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)

pdf(file.path(FIG_DIR, "UMAP_cell_types.pdf"), width=10, height=8)
DimPlot(seu, reduction="umap", group.by="Cell_type",
        label=TRUE, repel=TRUE, pt.size=0.1) +
  ggtitle("GSE131907 - Cell Types")
dev.off()

pdf(file.path(FIG_DIR, "UMAP_sample_origin.pdf"), width=10, height=8)
DimPlot(seu, reduction="umap", group.by="Sample_Origin",
        pt.size=0.1) +
  ggtitle("Sample Origin")
dev.off()

cat("  UMAP 完成 ✓\n")

# =============================================================================
# STEP 4: 提取 CAF 亚群
# =============================================================================
cat("\n[STEP 4] 提取 CAF 亚群...\n")

caf_subtypes <- c("COL14A1+ matrix FBs", "COL13A1+ matrix FBs",
                  "Myofibroblasts", "FB-like cells")

seu_CAF <- subset(seu, subset = Cell_subtype %in% caf_subtypes)
cat(sprintf("  CAF 总数: %d cells\n", ncol(seu_CAF)))
print(table(seu_CAF$Cell_subtype))

seu_CAF$CAF_subtype <- dplyr::recode(
  seu_CAF$Cell_subtype,
  "Myofibroblasts"      = "myCAF",
  "COL14A1+ matrix FBs" = "iCAF",
  "COL13A1+ matrix FBs" = "iCAF",
  "FB-like cells"       = "apCAF"
)

cat("\n  CAF 亚型分布:\n")
print(table(seu_CAF$CAF_subtype))

seu_CAF <- FindVariableFeatures(seu_CAF, nfeatures=1500, verbose=FALSE)
seu_CAF <- ScaleData(seu_CAF, verbose=FALSE)
seu_CAF <- RunPCA(seu_CAF, npcs=30, verbose=FALSE)
seu_CAF <- RunUMAP(seu_CAF, dims=1:20, verbose=FALSE)

pdf(file.path(FIG_DIR, "CAF_UMAP_subtypes.pdf"), width=8, height=7)
DimPlot(seu_CAF, reduction="umap", group.by="CAF_subtype",
        cols=c("myCAF"="#d62728","iCAF"="#1f77b4","apCAF"="#2ca02c"),
        label=TRUE, pt.size=0.5) +
  ggtitle("CAF Subtypes")
dev.off()

caf_markers_plot <- intersect(
  c("ACTA2","TAGLN","COL1A1","IL6","CXCL12","COL14A1","CD74","HLA-DRA","CXCL13"),
  rownames(seu_CAF)
)

pdf(file.path(FIG_DIR, "CAF_DotPlot_markers.pdf"), width=12, height=5)
DotPlot(seu_CAF, features=caf_markers_plot, group.by="CAF_subtype") +
  RotatedAxis() +
  ggtitle("CAF Subtype Marker Expression")
dev.off()

# =============================================================================
# STEP 5: DEG 分析
# =============================================================================
cat("\n[STEP 5] 差异表达分析 (约10-20分钟)...\n")

Idents(seu_CAF) <- "CAF_subtype"

all_markers <- FindAllMarkers(
  seu_CAF,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5,
  test.use        = "wilcox",
  verbose         = FALSE
)

cat(sprintf("  总 DEG 数: %d\n", nrow(all_markers)))
write.csv(all_markers,
          file.path(RESULT_DIR, "CAF_DEGs_all_subtypes.csv"),
          row.names=FALSE)

get_signature <- function(markers, subtype, n=50) {
  markers %>%
    filter(cluster == subtype,
           p_val_adj < 0.05,
           avg_log2FC > 0.5) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n) %>%
    pull(gene)
}

sig_myCAF <- get_signature(all_markers, "myCAF")
sig_iCAF  <- get_signature(all_markers, "iCAF")
sig_apCAF <- get_signature(all_markers, "apCAF")

cat(sprintf("\n  myCAF (%d genes): %s ...\n",
            length(sig_myCAF), paste(sig_myCAF[1:5], collapse=", ")))
cat(sprintf("  iCAF  (%d genes): %s ...\n",
            length(sig_iCAF),  paste(sig_iCAF[1:5],  collapse=", ")))
cat(sprintf("  apCAF (%d genes): %s ...\n",
            length(sig_apCAF), paste(sig_apCAF[1:5], collapse=", ")))

CAF_SIGNATURES <- list(
  myCAF = sig_myCAF,
  iCAF  = sig_iCAF,
  apCAF = sig_apCAF
)

top10 <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  pull(gene) %>%
  unique()

pdf(file.path(FIG_DIR, "CAF_DEG_heatmap.pdf"), width=14, height=10)
DoHeatmap(seu_CAF, features=top10, group.by="CAF_subtype", size=3) +
  scale_fill_gradientn(colors=c("#3182bd","white","#d62728"))
dev.off()

# =============================================================================
# STEP 6: 保存
# =============================================================================
cat("\n[STEP 6] 保存结果...\n")

saveRDS(seu_CAF,         file.path(RESULT_DIR, "seu_CAF_annotated.rds"))
saveRDS(CAF_SIGNATURES,  file.path(RESULT_DIR, "CAF_subtype_signatures.rds"))

sig_df <- data.frame(
  Subtype = rep(names(CAF_SIGNATURES), sapply(CAF_SIGNATURES, length)),
  Gene    = unlist(CAF_SIGNATURES)
)
write.csv(sig_df, file.path(RESULT_DIR, "CAF_signatures.csv"), row.names=FALSE)

cat("\n✅ Step 2 完成！\n")
cat(sprintf("  → %s/seu_CAF_annotated.rds\n",        RESULT_DIR))
cat(sprintf("  → %s/CAF_subtype_signatures.rds\n",   RESULT_DIR))
cat(sprintf("  → %s/CAF_DEGs_all_subtypes.csv\n",    RESULT_DIR))
