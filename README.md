# scani-luad

Code and derived data for:

> **A cancer-associated fibroblast–derived transcriptomic index (SCANI) and overall survival in lung adenocarcinoma: derivation and evaluation in an independent cohort.**
> Wang Y, Xu W, Yin J, Yin H. Submitted to *Translational Cancer Research*.

SCANI (Standardized CAF Activation Net Index) = Z(myCAF ssGSEA score) + Z(TIMP1) − Z(PDPN), with each component
Z-standardized within the analysed cohort.
The myCAF score is a single-sample gene set enrichment score (GSVA 1.48.3, `method = "ssgsea"`, `ssgsea.norm = FALSE`)
of a 50-gene myofibroblastic CAF signature derived from the GSE131907 single-cell atlas.

**Version.** Tag `v9.1` corresponds to the submitted manuscript (version 9, revised 7 November 2026; Figure 1 wording,
Figure 4 log-rank P shown to three decimals). Tag `v9` is the first release of version 9. The TCGA-LUAD discovery cohort is defined with the
TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR): stage II–IV, positive overall-survival time and primary-tumour RNA sequencing,
**n = 226 (115 deaths)**. Earlier versions used a 166-patient table built from Genomic Data Commons clinical fields; see *Analysis history*.

---

## Quick check

```bash
Rscript verify_numbers.R
```

`verify_numbers.R` needs only `data_derived/` and the R packages `survival` and (for the AUC checks) `timeROC`. It

* **recomputes** from the patient-level derived tables the Cox, Kaplan–Meier, time-dependent AUC, C-index and reverse
  Kaplan–Meier follow-up results for TCGA-LUAD (n = 226), its components, the earlier 166-patient table, the 60 patients
  missing from that table, GSE72094 (genome-wide and 21-gene versions, and components) and the CPTAC-LUAD stage correlation;
* **checks** the additional analyses that need the expression matrices (cohort flow, stage subgroups, further adjustment,
  random gene sets, ESTIMATE correlations, spline tests) against the analysis output files;
* prints PASS / FAIL / SKIP for each reported value and writes `verify_output/verify_numbers_result.csv`.

Single-cell counts and the spatial program are not covered by `verify_numbers.R`; they are
asserted inside the scripts that produce them.

---

## Data sources (all public)

| Dataset | Use | Source |
|---|---|---|
| GSE131907 | LUAD single-cell atlas; fibroblast annotation by the original authors; myCAF signature | GEO |
| GSE307534 | Visium spatial transcriptomics, 56 sections / 25 patients (descriptive only) | GEO |
| TCGA-LUAD | RNA-seq (Genomic Data Commons); clinical data from TCGA-CDR (Liu et al., *Cell* 2018, Table S1) | GDC; TCGA-CDR |
| GSE72094 | Independent validation cohort (GPL15048 microarray) | GEO |
| CPTAC-LUAD | SCANI vs pathological stage only (no survival) | LinkedOmics; Gillette et al. 2020 |

Raw data are not redistributed. `data_derived/` contains the signature, per-patient scores and analysis outputs.

## Environment

R 4.3.1; GSVA 1.48.3 (old `gsva()` interface); survival 3.5.7; timeROC 0.4; Seurat 4.3; ggplot2; cowplot; readxl;
tidyestimate 1.1.1 (ESTIMATE scores); Python 3 with matplotlib (Figure 1 only).

**Paths.** The analysis scripts were run on the authors' server and use absolute paths under `~/caf_spatial_project/` (upstream
pipeline and raw data) and `~/thesis_LUAD_CAF/05_manuscripts/SCI-1_IJMS/` (outputs). To rerun them, recreate this layout or edit
the path variables at the top of each script. All scripts are read-only with respect to their inputs and stop rather than
overwrite an existing output.

---

## Repository layout and run order

```
scripts/
  00_upstream/        project pipeline (only the parts listed below are used by the SCANI analysis)
  01_single_cell/     myCAF signature table, Figure 2
  02_spatial/         spatial program, Figure 3
  03_tcga/            TCGA-LUAD analyses, Figure 4
  04_gse72094/        GSE72094 scores, Figure 5
  04b_cptac/          CPTAC-LUAD stage analysis
  05_supplementary/   Tables S4 and S5 inputs
  06_figure1/         workflow diagram
  99_audit/           read-only provenance checks
data_derived/         derived tables (see below)
figures/              submitted figures
verify_numbers.R
```

Run in this order (each script lists its inputs and outputs in its header):

| # | Script | Main outputs | Used for |
|---|---|---|---|
| 0a | `00_upstream/step2_scrna_caf.R` | `CAF_signatures.csv`, `seu_CAF_annotated.rds`, `CAF_subtype_signatures.rds` | myCAF signature |
| 0b | `00_upstream/step3_spatial_GSE307534.R` | per-section differential expression, `CAF_Niche_final_geneset.txt` | spatial program |
| 0c | `00_upstream/step6b_niche_score_v2.R` | `TCGA_LUAD_niche_survival_v2.csv` (GDC clinical fields) | earlier 166-patient table only |
| 0d | `00_upstream/step6d_genes.R` (lines 1–17) | TCGA expression matrix `expr_mat`, `pid_to_col` (first aliquot per patient) | all TCGA scores |
| 0e | `00_upstream/step6e_final_model.R` | `CAF_final_scores.csv` | earlier 166-patient table only |
| 0f | `00_upstream/step22_GSE72094_validation.R` (lines 1–145) | GSE72094 series matrix and probe annotation in memory | GSE72094 scores |
| 1 | `01_single_cell/sci1_fig2_and_zorder.R` | `S1_Table_myCAF_signature_50genes.csv`, `fig2_counts_origin_by_subtype.csv` | Table S1, Figure 2A |
| 2 | `03_tcga/sci1_fig4_O.R` | `SCANI_TCGA_O_166.csv` (intermediate input for steps 3 and 6; its figures are superseded) | — |
| 3 | `01_single_cell/sci1_figs_final.R` | Figure 2 (its Figure 4 section is superseded) | Figure 2 |
| 4 | `02_spatial/sci1_spatial_fig3.R` | recurrence table, GO/KEGG enrichment, Figure 3 | Figure 3 |
| 5 | `03_tcga/sci1_O_v2_same_aliquot.R` | `SCANI_TCGA_O_166_v2.csv` | earlier 166-patient table (Table S5 A) |
| 6 | `04_gse72094/sci1_recheck2.R` | `SCANI_GSE72094_both_versions.csv`, `GSE72094_asrun_vs_fullgenome.csv` | GSE72094 (Table 1, Table S3) |
| 7 | `01_single_cell/sci1_supp_export.R` | text export of the signature statistics and GSE72094 versions | Tables S1, S3 |
| 8 | `04_gse72094/sci1_fig5.R` | Figure 5 | Figure 5 |
| 8b | `04b_cptac/06_baseA_cptac_rev2.R` | cleaned CPTAC-LUAD expression (`CPTAC_LUAD_expr_symbol.rds`) and clinical table | CPTAC inputs |
| 8c | `04b_cptac/sci1_cptac_stage.R` | `cptac/SCANI_CPTAC_stage.csv`, `cptac/cptac_stage_spearman.csv` | CPTAC-LUAD SCANI vs stage (Results, Table 1) |
| 9 | `05_supplementary/sci1_remark_extras.R` | `remark_S4_*_v2.csv` | Table S4, GSE72094 rows |
| 10 | `05_supplementary/sci1_reviewer_prep.R` | `reviewer_prep/` | Table S5, GSE72094 parts and TCGA stage I–IV |
| 11 | `03_tcga/sci1_diag_166_vs_226.R` | log only | reason for the 60 missing patients |
| 12 | `03_tcga/sci1_cdr226_main.R` | `cdr226/`, `SCANI_TCGA_CDR226.csv`, Figure 4 | **primary TCGA analysis**: Table 1, Table S2, Table S4 (TCGA rows), Table S5, Figure 4 |
| 12b | `03_tcga/sci1_fig4_p3.R` | `Figure4_cdr226_p3.*` | submitted Figure 4 (same plotting code as step 12, log-rank P to three decimals) |
| 13 | `05_supplementary/sci1_estimate_corr.R` | `cdr226/estimate_correlations.csv` | Table S5 E |
| 14 | `06_figure1/make_fig1.py`, `06_figure1/make_fig1_pptx.py` | Figure 1 (image and editable PowerPoint) | Figure 1 |
| 15 | `verify_numbers.R` | `verify_output/` | check |

CPTAC-LUAD: gene-level log2 RNA-seq data from LinkedOmics (`cptac_rna.cct`) and pathological stage from the sample annotation
of Gillette et al. (*Cell* 2020, Table S1A); 76 of 110 tumour samples have recorded stage (I–III). Earlier, unaudited CPTAC
results on the server (n = 72; rho 0.148 / 0.157) used a different scoring implementation and are superseded by `sci1_cptac_stage.R`.

## Derived data

| File | Content |
|---|---|
| `CAF_signatures.csv`, `S1_Table_myCAF_signature_50genes.csv` | the 50-gene myCAF signature (and the other states) |
| `fig2_counts_origin_by_subtype.csv` | fibroblast counts by state and tissue of origin |
| `SCANI_TCGA_CDR226.csv` | **primary** TCGA cohort: per-patient myCAF, TIMP1, PDPN, SCANI, OS, stage, age, sex, aliquot, membership of the earlier table |
| `cdr226/` | TCGA outputs: cohort flow, headline numbers, Table S2 (candidate genes, components), Table S4 (TCGA rows), sensitivity 166 vs 226, additional analyses, random gene sets, ESTIMATE scores and correlations |
| `SCANI_GSE72094_both_versions.csv` | GSE72094 per-sample scores: genome-wide (`SCANI_full`, primary) and 21-gene (`SCANI_asrun`, sensitivity) |
| `remark_S4_*_v2.csv` | Table S4 inputs (GSE72094 rows are used; TCGA rows are superseded by `cdr226/`) |
| `reviewer_prep/` | additional analyses: component combinations, likelihood-ratio tests, stage subgroups, further adjustment, random gene sets, spline tests |
| `SCANI_TCGA_O_166_v2.csv` | earlier 166-patient table (sensitivity analysis) |
| `intermediate/SCANI_TCGA_O_166.csv` | first version of the 166-patient table; input to steps 3 and 6 only; not reported |
| `spatial/` | spatial gene recurrence, the 50-gene program, enrichment results |
| `cptac/` | CPTAC-LUAD per-sample myCAF, TIMP1, PDPN, SCANI and stage; Spearman results |

---

## Analysis history (why there are several TCGA tables)

1. **Duplicated records.** The earlier scoring table had 192 rows for 166 patients (one patient, TCGA-44-2665, had three aliquots and
   nine identical clinical records); the component-screening table had 174 rows. All analyses de-duplicate to one record per patient
   before standardization. TIMP1 and PDPN were selected from 16 candidate genes on the 174-row table; the formula was not changed afterwards.
2. **Same aliquot (v1 → v2).** In v1 the myCAF component of TCGA-44-2665 came from a different aliquot than TIMP1/PDPN. v2 takes all
   three components from the first-listed aliquot (`sci1_O_v2_same_aliquot.R`); any of the three aliquots changes the HR by < 0.01.
3. **Discovery cohort (166 → 226).** The 166-patient table was built from GDC clinical fields in which pathological stage was missing
   for some deceased patients and follow-up time for some living patients. Defining the cohort with TCGA-CDR as stated in the Methods
   gives 226 patients (115 deaths); 60 eligible patients (35 deaths) had been missing (`sci1_diag_166_vs_226.R`). The 226-patient cohort
   is the primary analysis; the 166-patient table is reported as a sensitivity analysis (HR per SD 1.46 vs 1.45).
4. **GSE72094 scoring.** The original implementation (`step22`) scored the myCAF signature within a 55-probe (21-gene) subset. The
   primary analysis uses the method stated in the Methods (genome-wide probe collapsing, all 50 signature genes); the 21-gene version
   is a sensitivity analysis (Table S3).
5. **Additional analyses** in Table S5 were not prespecified and were added after internal review.

Removed from earlier releases: `external_validation_GSE72094.R` (it computed an ECM/CXCL score, not SCANI) and the previous
`verify_numbers.R` / README (166-patient table without the same-aliquot rule).

## Licence

To be added by the authors.
