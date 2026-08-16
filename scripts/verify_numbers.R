#!/usr/bin/env Rscript
# verify_numbers.R — 复现门闸:从派生数据复算论文每个关键数字,不符即报错
suppressMessages({ library(survival); library(timeROC) })
DD <- "data_derived"; pass <- TRUE
chk <- function(name, computed, claimed, tol) {
  ok <- is.finite(computed) && abs(computed - claimed) <= tol
  cat(sprintf("[%s] %-30s computed=%.4f claimed=%.4f\n", ifelse(ok,"PASS","FAIL"), name, computed, claimed))
  if (!ok) pass <<- FALSE
}
# TCGA (n=166, 去重)
s <- read.csv(file.path(DD,"SCANI_TCGA_scores.csv")); s <- s[!duplicated(s$patient_id),]
chk("TCGA n", nrow(s), 166, 0); chk("TCGA events", sum(s$OS_event==1), 80, 0)
s$z <- as.numeric(scale(s$CAF_Index))
chk("TCGA univariable HR", exp(coef(coxph(Surv(OS_time,OS_event)~z,s))), 1.46, 0.02)
km <- survdiff(Surv(OS_time,OS_event)~I(s$z>=median(s$z)),s)
chk("TCGA KM log-rank p", pchisq(km$chisq,1,lower.tail=FALSE), 0.23, 0.03)
tr <- timeROC(T=s$OS_time,delta=s$OS_event,marker=s$z,cause=1,times=c(365,1095,1825),iid=FALSE)
chk("TCGA 3yr AUC", tr$AUC[2], 0.635, 0.02)
# GSE72094 (n=398)
g <- read.csv(file.path(DD,"SCANI_GSE72094_scores.csv")); g <- g[!is.na(g$OS_days)&!is.na(g$OS_event),]
chk("GSE72094 n", nrow(g), 398, 0)
kmg <- survdiff(Surv(OS_days,OS_event)~I(g$SCANI>=median(g$SCANI)),g)
chk("GSE72094 log-rank p", pchisq(kmg$chisq,1,lower.tail=FALSE), 0.031, 0.005)
chk("GSE72094 HR/SD", exp(coef(coxph(Surv(OS_days,OS_event)~scale(SCANI),g))), 1.21, 0.03)
# 富集 top term
go <- read.csv(file.path(DD,"GO_BP_enrichment.csv")); ke <- read.csv(file.path(DD,"KEGG_enrichment.csv"))
chk("GO top p.adj (log10)", log10(min(go$p.adjust)), log10(4.90e-10), 0.2)
chk("KEGG top p.adj (log10)", log10(min(ke$p.adjust)), log10(1.40e-19), 0.2)
cat("\n", strrep("=",44), "\n", sep="")
if (pass) { cat("ALL CHECKS PASSED\n"); quit(status=0) } else { cat("SOME CHECKS FAILED\n"); quit(status=1) }
