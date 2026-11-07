#!/usr/bin/env Rscript
# ==============================================================================
# verify_numbers.R —— recompute the reported statistics of the SCANI manuscript (TCR v9) from the derived
# patient-level tables in data_derived/, and check every statistic that is read from an analysis output file.
#
#   Recomputed from patient-level data (Cox, Kaplan–Meier, time-dependent AUC, C-index, reverse-KM follow-up):
#     TCGA-LUAD stage II–IV (TCGA-CDR cohort, n=226), its components, the earlier 166-patient table, the 60 added patients;
#     GSE72094 genome-wide SCANI (primary) and the 21-gene implementation (sensitivity), and the components;
#     CPTAC-LUAD SCANI vs pathological stage (Spearman).
#   Checked against the analysis output files (not recomputed here, because they need the expression matrices):
#     cohort flow, stage subgroups, further adjustment, random gene sets, ESTIMATE correlations, spline tests.
#
#   Each check compares the value, rounded as in the manuscript, with the reported text. Result: PASS / FAIL / SKIP.
#   Usage (from the repository root):  Rscript verify_numbers.R
#   Requires: R >= 4.3, survival; timeROC for the AUC checks (otherwise those checks are SKIP).
# ==============================================================================
suppressPackageStartupMessages(library(survival))
DD <- "data_derived"
stopifnot(dir.exists(DD))
zs <- function(x) as.numeric(scale(x))
HAS_TROC <- requireNamespace("timeROC", quietly = TRUE)
pf <- function(x) { if (x < 0.001) return("<0.001"); if (x > 0.99) return(">0.99"); if (x < 0.10) sprintf("%.3f", x) else sprintf("%.2f", x) }
f2 <- function(x) sprintf("%.2f", x)
RES <- list()
chk <- function(id, got, expected) {
  st <- if (is.null(got) || (length(got) == 1 && is.na(got))) "SKIP" else if (identical(got, expected)) "PASS" else "FAIL"
  RES[[length(RES) + 1]] <<- data.frame(check = id, reported = expected, recomputed = if (is.null(got)) NA else got, status = st)
}
rd <- function(f) { p <- file.path(DD, f); if (!file.exists(p)) { message("missing: ", p); return(NULL) }; read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) }
cox_s <- function(t, e, x) { s <- summary(coxph(Surv(t, e) ~ zs(x)))
  sprintf("%s (%s–%s), P%s", f2(s$conf.int[1, 1]), f2(s$conf.int[1, 3]), f2(s$conf.int[1, 4]),
          { p <- pf(s$coefficients[1, 5]); if (substr(p, 1, 1) %in% c("<", ">")) p else paste0("=", p) }) }
coxadj_s <- function(dd, x, covs) { dd$x. <- zs(x); s <- summary(coxph(as.formula(paste("Surv(t, e) ~ x. +", paste(covs, collapse = "+"))), dd))
  sprintf("%s (%s–%s), P%s; n=%d/%d", f2(s$conf.int["x.", 1]), f2(s$conf.int["x.", 3]), f2(s$conf.int["x.", 4]),
          { p <- pf(s$coefficients["x.", 5]); if (substr(p, 1, 1) %in% c("<", ">")) p else paste0("=", p) }, s$n, s$nevent) }
kmp <- function(t, e, x) pf(1 - pchisq(survdiff(Surv(t, e) ~ I(x >= median(x)))$chisq, 1))
auc3 <- function(t, e, x) if (HAS_TROC) f2(timeROC::timeROC(T = t, delta = e, marker = x, cause = 1, times = 3 * 365, iid = FALSE)$AUC[2]) else NA
cidx <- function(t, e, x) f2(concordance(Surv(t, e) ~ x, reverse = TRUE)$concordance)
revkm <- function(t, e) { s <- summary(survfit(Surv(t, 1 - e) ~ 1))$table / 30.44; sprintf("%.1f (%.1f–%.1f)", s["median"], s["0.95LCL"], s["0.95UCL"]) }

# ---------------------------------------------------------------- TCGA-LUAD, TCGA-CDR cohort (primary)
d <- rd("SCANI_TCGA_CDR226.csv")
if (!is.null(d)) {
  stopifnot(nrow(d) == 226, !anyDuplicated(d$patient_id), sum(d$OS_event) == 115)
  stopifnot(max(abs(d$SCANI - (zs(d$myCAF) + zs(d$TIMP1) - zs(d$PDPN)))) < 1e-8)
  d$t <- d$OS_time; d$e <- d$OS_event; d$stage <- factor(d$stage, levels = c("II", "III", "IV"))
  cc <- complete.cases(d[, c("age", "stage")])
  chk("TCGA n/deaths", sprintf("%d/%d", nrow(d), sum(d$e)), "226/115")
  chk("TCGA SCANI univariable HR per SD", cox_s(d$t, d$e, d$SCANI), "1.45 (1.17–1.80), P<0.001")
  chk("TCGA SCANI adjusted age+stage", coxadj_s(d[cc, ], d$SCANI[cc], c("age", "stage")), "1.47 (1.18–1.83), P<0.001; n=223/112")
  chk("TCGA median-split log-rank P", kmp(d$t, d$e, d$SCANI), "0.016")
  chk("TCGA 3-year AUC", auc3(d$t, d$e, d$SCANI), "0.62")
  chk("TCGA C-index", cidx(d$t, d$e, d$SCANI), "0.60")
  chk("TCGA median follow-up, months", revkm(d$t, d$e), "30.0 (26.5–37.1)")
  chk("TCGA at risk at 3 years", as.character(sum(d$t >= 3 * 365)), "45")
  chk("TCGA TIMP1 alone", cox_s(d$t, d$e, d$TIMP1), "1.30 (1.07–1.58), P=0.007")
  chk("TCGA PDPN alone", cox_s(d$t, d$e, d$PDPN), "0.86 (0.71–1.04), P=0.13")
  chk("TCGA myCAF alone", cox_s(d$t, d$e, d$myCAF), "1.07 (0.89–1.28), P=0.47")
  chk("TCGA TIMP1 3-year AUC", auc3(d$t, d$e, d$TIMP1), "0.58")
  chk("TCGA myCAF 3-year AUC", auc3(d$t, d$e, d$myCAF), "0.48")
  chk("TCGA myCAF–PDPN Spearman rho", f2(cor(d$myCAF, d$PDPN, method = "spearman")), "0.65")
  ct <- cor.test(d$SCANI, d$age, method = "spearman", exact = FALSE)
  chk("TCGA SCANI vs age rho, P", sprintf("%s, P=%s", f2(ct$estimate), pf(ct$p.value)), "0.13, P=0.057")
  chk("TCGA SCANI vs stage Kruskal–Wallis P", pf(kruskal.test(SCANI ~ stage, d)$p.value), "0.35")
  chk("TCGA SCANI vs sex Wilcoxon P", pf(wilcox.test(SCANI ~ sex, d)$p.value), "0.33")
  chk("TCGA stage II/III/IV", paste(table(d$stage), collapse = "/"), "120/81/25")
  chk("TCGA female/male", paste(table(d$sex)[c("female", "male")], collapse = "/"), "110/116")
  mv <- coxph(Surv(t, e) ~ zs(SCANI) + age + stage, d[cc, ])
  chk("TCGA proportional hazards, global P", pf(cox.zph(mv)$table["GLOBAL", "p"]), "0.30")
  n60 <- d[!d$in_original_166, ]
  chk("TCGA 60 added patients", cox_s(n60$t, n60$e, n60$SCANI), "1.44 (0.97–2.14), P=0.071")
  # likelihood-ratio test: TIMP1 vs TIMP1 + [Z(myCAF) - Z(PDPN)]
  d$zT <- zs(d$TIMP1); d$zMP <- zs(d$myCAF) - zs(d$PDPN)
  l <- anova(coxph(Surv(t, e) ~ zT, d), coxph(Surv(t, e) ~ zT + zMP, d))
  chk("TCGA LRT adding myCAF-PDPN to TIMP1", pf(l[[ncol(l)]][2]), "0.018")
}

# ---------------------------------------------------------------- TCGA-LUAD, earlier 166-patient table (sensitivity)
o <- rd("SCANI_TCGA_O_166_v2.csv")
if (!is.null(o)) {
  o <- o[!duplicated(o$patient_id), ]; stopifnot(nrow(o) == 166, sum(o$OS_event) == 80)
  o$t <- o$OS_time; o$e <- o$OS_event; o$age <- o$age_at_index; o$stage <- factor(o$stage_simple, levels = c("II", "III", "IV"))
  chk("Earlier table (166) univariable", cox_s(o$t, o$e, o$SCANI), "1.46 (1.13–1.89), P=0.004")
  chk("Earlier table (166) adjusted", coxadj_s(o, o$SCANI, c("age", "stage")), "1.49 (1.14–1.93), P=0.003; n=166/80")
}

# ---------------------------------------------------------------- GSE72094
g <- rd("SCANI_GSE72094_both_versions.csv")
if (!is.null(g)) {
  g <- g[!duplicated(g$sample_id), ]; stopifnot(nrow(g) == 398, sum(g$OS_event) == 113)
  stopifnot(max(abs(g$SCANI_full - (zs(g$myCAF_full) + zs(g$TIMP1) - zs(g$PDPN)))) < 1e-8)
  g$t <- g$OS_days; g$e <- g$OS_event; g$age <- suppressWarnings(as.numeric(g$age)); gc <- complete.cases(g[, c("age", "stage_num")])
  chk("GSE72094 SCANI univariable", cox_s(g$t, g$e, g$SCANI_full), "1.23 (1.00–1.51), P=0.045")
  chk("GSE72094 SCANI adjusted age+stage", coxadj_s(g[gc, ], g$SCANI_full[gc], c("age", "stage_num")), "1.20 (0.99–1.47), P=0.067; n=393/111")
  chk("GSE72094 median-split log-rank P", kmp(g$t, g$e, g$SCANI_full), "0.004")
  chk("GSE72094 3-year AUC", auc3(g$t, g$e, g$SCANI_full), "0.59")
  chk("GSE72094 C-index", cidx(g$t, g$e, g$SCANI_full), "0.60")
  chk("GSE72094 median follow-up, months", revkm(g$t, g$e), "29.9 (29.0–31.2)")
  m2 <- coxph(Surv(t, e) ~ zs(SCANI_full) + age + stage_num, g[gc, ])
  chk("GSE72094 PH test for SCANI (univariable)", pf(cox.zph(coxph(Surv(t, e) ~ zs(SCANI_full), g))$table[1, "p"]), "0.011")
  chk("GSE72094 21-gene univariable", cox_s(g$t, g$e, g$SCANI_asrun), "1.21 (1.00–1.46), P=0.044")
  chk("GSE72094 21-gene log-rank P", kmp(g$t, g$e, g$SCANI_asrun), "0.031")
  chk("GSE72094 21-gene adjusted", coxadj_s(g[gc, ], g$SCANI_asrun[gc], c("age", "stage_num")), "1.19 (0.99–1.42), P=0.066; n=393/111")
  chk("GSE72094 TIMP1 alone", cox_s(g$t, g$e, g$TIMP1), "1.11 (0.92–1.35), P=0.28")
  chk("GSE72094 myCAF alone", cox_s(g$t, g$e, g$myCAF_full), "1.32 (1.08–1.60), P=0.006")
  chk("GSE72094 PDPN alone", cox_s(g$t, g$e, g$PDPN), "1.09 (0.89–1.32), P=0.41")
  g$zT <- zs(g$TIMP1); g$zMP <- zs(g$myCAF_full) - zs(g$PDPN)
  l <- anova(coxph(Surv(t, e) ~ zT, g), coxph(Surv(t, e) ~ zT + zMP, g))
  chk("GSE72094 LRT adding myCAF-PDPN to TIMP1", pf(l[[ncol(l)]][2]), "0.040")
  st <- cor.test(g$SCANI_full, g$stage_num, method = "spearman", exact = FALSE)
  chk("GSE72094 SCANI vs stage rho, P", sprintf("%s, P=%s", f2(st$estimate), pf(st$p.value)), "0.12, P=0.014")
}

# ---------------------------------------------------------------- CPTAC-LUAD (stage only)
cp <- rd("cptac/SCANI_CPTAC_stage.csv")
if (!is.null(cp)) {
  stopifnot(nrow(cp) == 76, !anyDuplicated(cp$sample_id))
  stopifnot(max(abs(cp$SCANI - (zs(cp$myCAF) + zs(cp$TIMP1) - zs(cp$PDPN)))) < 1e-8)
  ct <- cor.test(cp$SCANI, cp$stage_n, method = "spearman", exact = FALSE)
  chk("CPTAC SCANI vs stage rho, P (n)", sprintf("%s, P=%s (n=%d)", f2(ct$estimate), pf(ct$p.value), nrow(cp)), "0.14, P=0.22 (n=76)")
  chk("CPTAC stage I/II/III", paste(table(cp$stage_n), collapse = "/"), "46/17/13")
}

# ---------------------------------------------------------------- values read from analysis output files
fl <- rd("cdr226/cohort_flow_TCGA_CDR226.csv")
if (!is.null(fl)) chk("TCGA cohort flow (n)", paste(fl$n, collapse = " > "), "522 > 517 > 504 > 8 > 270 > 226")
sc <- rd("reviewer_prep/reviewer_C_stage_subsets.csv")
if (!is.null(sc)) { r <- function(a) sc[sc$analysis == a, ]
  x <- r("TCGA all stages I-IV, SCANI (Z within all-stage cohort)"); chk("TCGA all stages", sprintf("%d: %s, adj %s", x$n, f2(x$HR), f2(x$adj_HR)), "504: 1.34, adj 1.29")
  x <- r("TCGA stage I, SCANI re-standardized in subset"); chk("TCGA stage I", sprintf("%d: %s, P=%s", x$n, f2(x$HR), pf(x$p)), "270: 1.06, P=0.62")
  x <- r("GSE72094 stage I, SCANI re-standardized in subset"); chk("GSE72094 stage I", sprintf("%d: %s, P=%s", x$n, f2(x$HR), pf(x$p)), "254: 1.10, P=0.51")
  x <- r("GSE72094 stage II/III/IV, SCANI re-standardized in subset"); chk("GSE72094 stage II-IV", sprintf("%d: %s, P=%s; adj %s, P=%s", x$n, f2(x$HR), pf(x$p), f2(x$adj_HR), pf(x$adj_p)), "139: 1.16, P=0.27; adj 1.25, P=0.13") }
dc <- rd("reviewer_prep/reviewer_D_confounding.csv")
if (!is.null(dc)) { r <- function(m) dc[dc$cohort == "GSE72094" & dc$model == m, ]
  chk("GSE72094 + sex", f2(r("+ sex")$SCANI_HR), "1.16"); chk("GSE72094 + ESTIMATE", f2(r("+ ESTIMATE stromal + immune")$SCANI_HR), "1.18")
  x <- r("+ all of: sex, mutation/smoking found"); chk("GSE72094 + all covariates", sprintf("%s, P=%s", f2(x$SCANI_HR), pf(x$p)), "1.14, P=0.19") }
rv <- rd("cdr226/reviewer_TCGA_CDR226.csv")
if (!is.null(rv)) { r <- function(i, a) rv[rv$item == i & rv$analysis == a, ]
  chk("TCGA + sex", f2(r("D_confounding", "+ sex")$HR), "1.47"); chk("TCGA + ESTIMATE", f2(r("D_confounding", "+ ESTIMATE stromal + immune")$HR), "1.60")
  chk("TCGA spline non-linearity P", pf(r("F_spline", "nonlinearity LR test")$p), "0.53")
  chk("TCGA TIMP1-PDPN 3-year AUC", f2(r("A_variant", "TIMP1_minus_PDPN")$AUC3y), "0.66") }
sp <- rd("reviewer_prep/reviewer_F_spline.csv"); if (!is.null(sp)) chk("GSE72094 spline non-linearity P", pf(sp$nonlinearity_LR_p[sp$cohort == "GSE72094"][1]), "0.030")
te <- rd("cdr226/reviewer_E_random_TCGA_CDR226.csv")
if (!is.null(te)) { s <- te[te$set == "SUMMARY", ]; rr <- te[te$set != "SUMMARY", ]; stopifnot(nrow(rr) == 1000)
  chk("TCGA random sets: HR >= SCANI", sprintf("%.1f%%", 100 * mean(rr$rand_SCANI_HR >= s$rand_SCANI_HR)), "3.8%")
  chk("TCGA random sets: HR >= myCAF alone", sprintf("%.1f%%", 100 * mean(rr$rand_alone_HR >= s$rand_alone_HR)), "22.6%")
  chk("TCGA random sets: P<0.05", sprintf("%.1f%%", 100 * mean(rr$rand_alone_p < 0.05)), "9.8%") }
ge <- rd("reviewer_prep/reviewer_E_random_draws.csv"); ga <- rd("reviewer_prep/reviewer_A_component_variants.csv")
if (!is.null(ge) && !is.null(ga)) { rr <- ge[ge$cohort == "GSE72094", ]; stopifnot(nrow(rr) == 1000)
  obs <- ga$HR[ga$cohort == "GSE72094" & ga$variant == "SCANI"]; my <- ga$HR[ga$cohort == "GSE72094" & ga$variant == "myCAF_only"]
  chk("GSE72094 random sets: HR >= SCANI", sprintf("%.1f%%", 100 * mean(rr$rand_SCANI_HR >= obs)), "5.7%")
  chk("GSE72094 random sets: HR >= myCAF alone", sprintf("%.1f%%", 100 * mean(rr$rand_alone_HR >= my)), "6.2%")
  chk("GSE72094 random sets: P<0.05", sprintf("%.1f%%", 100 * mean(rr$rand_alone_p < 0.05)), "34.2%") }
ec <- rd("cdr226/estimate_correlations.csv")
if (!is.null(ec)) { r <- function(c, s, e) f2(ec$rho[ec$cohort == c & ec$score == s & ec$estimate == e])
  chk("ESTIMATE stromal rho, TCGA myCAF/PDPN/SCANI", paste(r("TCGA-LUAD", "myCAF", "stromal"), r("TCGA-LUAD", "PDPN", "stromal"), r("TCGA-LUAD", "SCANI", "stromal"), sep = "/"), "0.73/0.72/0.31")
  chk("ESTIMATE stromal rho, GSE72094 myCAF/PDPN/SCANI", paste(r("GSE72094", "myCAF", "stromal"), r("GSE72094", "PDPN", "stromal"), r("GSE72094", "SCANI", "stromal"), sep = "/"), "0.62/0.64/0.29")
  chk("ESTIMATE purity rho, GSE72094 SCANI", r("GSE72094", "SCANI", "purity"), "-0.13") }

R <- do.call(rbind, RES)
for (i in seq_len(nrow(R))) cat(sprintf("[%s] %-48s reported: %-40s recomputed: %s\n", R$status[i], R$check[i], R$reported[i], R$recomputed[i]))
dir.create("verify_output", showWarnings = FALSE)
write.csv(R, file.path("verify_output", "verify_numbers_result.csv"), row.names = FALSE)
cat(sprintf("\nPASS %d   FAIL %d   SKIP %d\n", sum(R$status == "PASS"), sum(R$status == "FAIL"), sum(R$status == "SKIP")))
if (any(R$status == "FAIL")) quit(status = 1)
