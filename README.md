# SCANI: a CAF-derived transcriptomic prognostic index for lung adenocarcinoma

Analysis code and derived data for the manuscript:

> **SCANI, a cancer-associated-fibroblast–derived transcriptomic index, provides modest but externally validated prognostic information beyond TNM stage in lung adenocarcinoma.** Wang Y, Xu W, Yin J, Yin H.

SCANI (SpatioNiche CAF Activation Index) is computable from bulk RNA-seq:

```
SCANI = Z(myCAF ssGSEA) + Z(TIMP1) - Z(PDPN)
```

each component Z-scored within each cohort.

## Key results (all independently recomputed from source data)

| Analysis | Cohort | Estimate | p |
|---|---|---|---|
| Univariable Cox (per SD) | TCGA-LUAD (n=166) | HR 1.46 (1.13-1.89) | 0.004 |
| Multivariable Cox (+stage, age, sex) | TCGA-LUAD (n=166) | HR 1.45 | 0.005 |
| Time-dependent ROC (3-year) | TCGA-LUAD (n=166) | AUC 0.635 | - |
| Median-split Kaplan-Meier | TCGA-LUAD (n=166) | not significant | 0.23 |
| External validation, log-rank | GSE72094 (n=398) | High vs low | 0.031 |
| External validation, Cox (per SD) | GSE72094 (n=398) | HR 1.21 | 0.044 |
| Stage independence, Spearman | CPTAC-LUAD (n=72) | rho 0.15 | 0.21 |

SCANI is significant as a continuous prognostic variable in two independent cohorts but does not support single-cutpoint stratification (non-significant median split). It is reported as a graded stromal correlate of outcome, not a stand-alone classifier.

## Data sources (public)

| Dataset | Accession | Use |
|---|---|---|
| scRNA-seq (LUAD) | GEO GSE131907 | CAF subtyping, signature |
| Visium spatial | GEO GSE307534 | 50-gene signature derivation |
| Bulk RNA-seq (discovery) | TCGA-LUAD | SCANI construction, survival |
| Proteogenomics | CPTAC-LUAD (cBioPortal) | stage independence |
| Microarray survival | GEO GSE72094 | external validation |

Raw data are obtained from the above repositories and are not redistributed here.

## Reproduce

```bash
Rscript scripts/verify_numbers.R
```

`verify_numbers.R` is the reproducibility gate: it recomputes all reported statistics (HRs, AUC, log-rank p, enrichment p, sample sizes) directly from the derived-data files and errors out if any value disagrees with the manuscript. All checks currently PASS.

## Repository layout

```
data_derived/   verified derived data (SCANI scores, signature, enrichment)
scripts/        analysis scripts + verify_numbers.R
figures/        manuscript figures (regenerated from real data)
environment/    sessionInfo.txt
```

## Notes

The TCGA discovery cohort comprises 166 unique stage II-IV patients (80 deaths); a duplicated patient record in an intermediate table was collapsed prior to analysis. The 50-gene CAF signature (`data_derived/CAF_signature_50genes.csv`) was derived by differential expression of high-CAF-activity spatial regions.

## Use of AI tools

Large language models assisted with code drafting and manuscript language. All numerical results were independently recomputed from source data by the authors, who take full responsibility for the integrity of the data and analyses.

## License

MIT (code). Please cite the original data repositories when using derived data.
