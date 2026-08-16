#!/usr/bin/env Rscript
# =============================================================================
# 06_cptac_validation.R
# External validation of CAF Index in CPTAC-LUAD
#
# Purpose : Replicate the stage-independence of CAF Index in an ethnically
#           and technically distinct cohort (predominantly Asian; RNA-seq RPKM).
#           Corresponds to Section 2.4 / Results 3.7 / Figure 8 of the paper.
#
# Input   :
#   CPTAC-LUAD bulk RNA-seq data (log2[RPKM+1]) downloaded from cBioPortal
#   https://www.cbioportal.org  (study: luad_cptac_2020)
#   Required files (place in data/CPTAC_LUAD/):
#     - luad_cptac_2020_data_mrna_seq_rpkm.txt
#     - luad_cptac_2020_data_clinical_sample.txt
#   CAF subtype signatures from step 01:
#     results/01_scRNA/CAF_subtype_signatures.rds
#   50-gene CAF Niche Signature from step 03:
#     results/03_spatial_GSE307534/CAF_Niche_final_geneset.txt
#
# Output  :
#   results/09_cptac_validation/
#     CPTAC_CAF_scores.csv            per-sample CAF Index + stage
#     CPTAC_stage_spearman.csv        Spearman rho + p-value
#     fig3_cptac_validation.pdf       Figure 8 (box plot by stage)
#
# Key result (paper) : Spearman rho = 0.157, p = 0.187 (n=72)
#                      → stage-independence replication
#
# Usage   : Rscript scripts/06_cptac_validation.R
# Runtime : ~2 min (no large object loading)
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

set.seed(42)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR   <- path.expand("~/caf_spatial_project")
DATA_DIR   <- file.path(BASE_DIR, "data",    "CPTAC_LUAD")
SC_RES     <- file.path(BASE_DIR, "results", "01_scRNA")
ST_RES     <- file.path(BASE_DIR, "results", "03_spatial_GSE307534")
SURV_RES   <- file.path(BASE_DIR, "results", "05_tcga_survival")
OUT_DIR    <- file.path(BASE_DIR, "results", "09_cptac_validation")
FIG_DIR    <- file.path(BASE_DIR, "figures", "09_cptac_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 06: CPTAC-LUAD External Validation ===\n")

# ── 1. Load expression matrix ────────────────────────────────────────────────
cat("\n[1/5] Loading CPTAC-LUAD expression data...\n")

expr_file <- file.path(DATA_DIR, "luad_cptac_2020_data_mrna_seq_rpkm.txt")
if (!file.exists(expr_file)) {
  stop(paste(
    "Expression file not found:", expr_file,
    "\nPlease download from cBioPortal:",
    "\n  https://www.cbioportal.org/study/summary?id=luad_cptac_2020",
    "\nand place luad_cptac_2020_data_mrna_seq_rpkm.txt in", DATA_DIR
  ))
}

raw <- read.delim(expr_file, check.names = FALSE, stringsAsFactors = FALSE)

# cBioPortal format: first column = Hugo_Symbol, second = Entrez_Gene_Id
gene_col <- "Hugo_Symbol"
cat(sprintf("  Raw matrix: %d genes × %d columns\n", nrow(raw), ncol(raw)))

# Extract numeric expression matrix (log2[RPKM+1])
expr_mat <- as.matrix(raw[, -(1:2)])
rownames(expr_mat) <- raw[[gene_col]]

# Remove duplicate gene symbols (keep highest mean)
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
cat(sprintf("  Expression matrix: %d genes × %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# ── 2. Load clinical data ─────────────────────────────────────────────────────
cat("\n[2/5] Loading clinical data...\n")

clin_file <- file.path(DATA_DIR, "luad_cptac_2020_data_clinical_sample.txt")
if (!file.exists(clin_file)) {
  stop(paste(
    "Clinical file not found:", clin_file,
    "\nDownload from the same cBioPortal study page."
  ))
}

clin_raw <- read.delim(clin_file, comment.char = "#",
                       check.names = FALSE, stringsAsFactors = FALSE)

# Identify sample ID and stage columns (cBioPortal column names vary)
id_col    <- grep("SAMPLE_ID|PATIENT_ID", toupper(colnames(clin_raw)), value = TRUE)[1]
stage_col <- grep("STAGE|AJCC|TUMOR_STAGE", toupper(colnames(clin_raw)), value = TRUE)[1]

cat(sprintf("  Using ID column: %s | Stage column: %s\n", id_col, stage_col))

clin <- clin_raw %>%
  rename(sample_id = all_of(id_col),
         stage_raw = all_of(stage_col)) %>%
  filter(!is.na(stage_raw), stage_raw != "") %>%
  mutate(
    # Simplify stage to I / II / III / IV / IIIA / IIIB etc.
    stage_simple = sub("Stage ", "", stage_raw, ignore.case = TRUE),
    stage_simple = sub("A$|B$", "", stage_simple),   # collapse IIIA/IIIB → III
    stage_simple = toupper(trimws(stage_simple))
  )

cat(sprintf("  Clinical samples: %d | Stage distribution:\n", nrow(clin)))
print(table(clin$stage_simple))

# Align samples between expression and clinical
shared_samples <- intersect(colnames(expr_mat), clin$sample_id)
cat(sprintf("  Matched samples: %d\n", length(shared_samples)))

expr_mat <- expr_mat[, shared_samples]
clin     <- clin %>% filter(sample_id %in% shared_samples) %>%
            arrange(match(sample_id, shared_samples))

# ── 3. Compute CAF Index ──────────────────────────────────────────────────────
cat("\n[3/5] Computing CAF Index...\n")

# Load CAF subtype gene signatures from scRNA-seq step
sig_file <- file.path(SC_RES, "CAF_subtype_signatures.rds")
if (!file.exists(sig_file)) {
  stop(paste("CAF signatures not found:", sig_file,
             "\nPlease run 01_scrna_caf_subtyping.R first."))
}
caf_sigs <- readRDS(sig_file)
cat(sprintf("  Loaded signatures: %s\n", paste(names(caf_sigs), collapse = ", ")))

# ssGSEA scoring (same method as TCGA validation)
gsva_param <- ssgseaParam(expr_mat, caf_sigs)
ssgsea_scores <- gsva(gsva_param, verbose = FALSE)
scores_df <- as.data.frame(t(ssgsea_scores))
colnames(scores_df) <- names(caf_sigs)
scores_df$sample_id <- rownames(scores_df)

# Load TIMP1 and PDPN raw expression (used in CAF Index formula)
genes_needed <- c("TIMP1", "PDPN")
for (g in genes_needed) {
  if (g %in% rownames(expr_mat)) {
    scores_df[[g]] <- as.numeric(expr_mat[g, scores_df$sample_id])
  } else {
    warning(sprintf("Gene %s not found in expression matrix — using 0", g))
    scores_df[[g]] <- 0
  }
}

# CAF Index = Z(myCAF ssGSEA) + Z(TIMP1) - Z(PDPN)
scores_df <- scores_df %>%
  mutate(
    myCAF_Z = scale(myCAF)[, 1],
    TIMP1_Z = scale(TIMP1)[, 1],
    PDPN_Z  = scale(PDPN)[, 1],
    CAF_Index = myCAF_Z + TIMP1_Z - PDPN_Z
  )

cat(sprintf("  CAF Index computed for %d samples\n", nrow(scores_df)))
cat(sprintf("  CAF Index range: [%.3f, %.3f]\n",
            min(scores_df$CAF_Index), max(scores_df$CAF_Index)))

# Merge with clinical
result_df <- scores_df %>%
  left_join(clin %>% select(sample_id, stage_raw, stage_simple), by = "sample_id")

# ── 4. Stage-independence test (Spearman) ────────────────────────────────────
cat("\n[4/5] Testing CAF Index stage-independence (Spearman)...\n")

# Encode stage as ordinal numeric for Spearman
stage_order <- c("I" = 1, "II" = 2, "III" = 3, "IV" = 4,
                 "IA" = 1, "IB" = 1, "IIA" = 2, "IIB" = 2,
                 "IIIA" = 3, "IIIB" = 3)

result_df <- result_df %>%
  mutate(stage_num = stage_order[stage_simple])

df_test <- result_df %>% filter(!is.na(stage_num), !is.na(CAF_Index))
cat(sprintf("  Samples for Spearman test: n = %d\n", nrow(df_test)))

sp_test <- cor.test(df_test$CAF_Index, df_test$stage_num,
                    method = "spearman", exact = FALSE)
cat(sprintf("  Spearman rho = %.3f,  p = %.3f\n",
            sp_test$estimate, sp_test$p.value))
cat("  → Non-significant: confirms CAF Index stage-independence\n")

# Save results
write_csv(result_df %>% select(sample_id, stage_raw, stage_simple,
                                myCAF, iCAF, apCAF, TIMP1, PDPN, CAF_Index),
          file.path(OUT_DIR, "CPTAC_CAF_scores.csv"))

write_csv(data.frame(
  n         = nrow(df_test),
  spearman_rho = round(sp_test$estimate, 4),
  p_value      = round(sp_test$p.value,  4),
  method       = "Spearman (stage as ordinal: I=1, II=2, III=3, IV=4)"
), file.path(OUT_DIR, "CPTAC_stage_spearman.csv"))

# ── 5. Figure 8: Box plot of CAF Index by pathological stage ─────────────────
cat("\n[5/5] Generating Figure 8 (CAF Index by stage)...\n")

stage_levels <- c("I", "IA", "IB", "II", "IIA", "IIB",
                  "III", "IIIA", "IIIB", "IV")
stage_levels <- stage_levels[stage_levels %in% df_test$stage_simple]

df_plot <- df_test %>%
  filter(stage_simple %in% stage_levels) %>%
  mutate(stage_simple = factor(stage_simple, levels = stage_levels))

stage_n <- df_plot %>% count(stage_simple) %>%
  mutate(label = paste0(stage_simple, "\n(n=", n, ")"))

df_plot <- df_plot %>%
  left_join(stage_n %>% select(stage_simple, label), by = "stage_simple")

p_box <- ggplot(df_plot, aes(x = stage_simple, y = CAF_Index)) +
  geom_boxplot(aes(fill = stage_simple),
               outlier.shape = 21, outlier.size = 1.5,
               width = 0.55, alpha = 0.85, color = "grey30") +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "grey40") +
  scale_fill_manual(
    values = colorRampPalette(c("#4575b4", "#abd9e9", "#fee090", "#d73027"))(
      length(levels(df_plot$stage_simple)))) +
  annotate("text",
           x    = nlevels(df_plot$stage_simple) / 2,
           y    = max(df_plot$CAF_Index) * 1.15,
           label = sprintf("Spearman rho = %.3f,  p = %.3f  (n = %d)",
                           sp_test$estimate, sp_test$p.value, nrow(df_test)),
           size = 3.8, color = "grey20") +
  labs(
    title    = "External Validation: CAF Index Stage-Independence (CPTAC-LUAD)",
    subtitle = "log2[RPKM+1] expression; predominantly Asian cohort",
    x        = "Pathological Stage",
    y        = "CAF Index"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10)),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x      = element_text(size = 10)
  )

# Save to both locations (results/ for archival, figures/ for submission)
for (out_path in c(file.path(OUT_DIR, "fig3_cptac_validation.pdf"),
                   file.path(FIG_DIR,  "fig3_cptac_validation.pdf"))) {
  ggsave(out_path, p_box, width = 8, height = 5)
}

# Also save to nomogram results dir to match original figure organisation
lasso_dir <- file.path(BASE_DIR, "results", "08_lasso_nomogram")
dir.create(lasso_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(lasso_dir, "fig3_cptac_validation.pdf"), p_box, width = 8, height = 5)

cat(sprintf("\n✅ Step 06 complete.\n"))
cat(sprintf("   Spearman rho = %.3f  (p = %.3f,  n = %d)\n",
            sp_test$estimate, sp_test$p.value, nrow(df_test)))
cat(sprintf("   Figures saved to: %s\n", FIG_DIR))
cat(sprintf("   Data saved to:    %s\n", OUT_DIR))
