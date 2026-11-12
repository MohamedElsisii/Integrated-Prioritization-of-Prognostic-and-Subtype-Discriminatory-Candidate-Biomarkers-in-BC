# =============================================================================
# Code S2 - Kaplan-Meier screening of subtype-specific DEGs
#
# For each unique DEG of one PAM50 subtype, patients are split at the optimal
# expression cut-point (surv_cutpoint, minimum group proportion 10%) and overall
# survival is compared with the log-rank test. For every gene the script also
# reports:
#   - the number of patients and deaths in each expression group
#   - the hazard ratio (high vs low expression) with 95% CI
#   - a log-rank p-value corrected for cut-point selection (maxstat, Lausen 1992)
#   - a cut-point-free Cox hazard ratio per 1-SD increase in expression
# P-values are adjusted with the Benjamini-Hochberg method. For the final
# candidates, a conditional Monte Carlo permutation p-value is also computed.
#
# Run once per subtype by setting SUBTYPE below.
#
# Input (working directory)
#   <SUBTYPE>_Significant.rda   character vector of unique DEGs (from Code S1)
#   <SUBTYPE>_Samples.rda       sample IDs of the subtype
#   BRCA.rds                    raw counts (matrix or SummarizedExperiment)
#   Clinical_Data.csv           submitter_id, vital_status, days_to_death,
#                               days_to_last_follow_up
# Output
#   <SUBTYPE>_survival_results.csv   (Tables S7-S9)
#   survival_plots_<SUBTYPE>/        Kaplan-Meier plots (Files S1-S3)
# =============================================================================

suppressPackageStartupMessages({
  library(survminer)
  library(survival)
  library(tidyverse)
  library(DESeq2)
  library(maxstat)
})

SUBTYPE <- "LumA"   # "LumA", "LumB" or "Basal"
final_candidates <- list(
  LumA  = "F2RL2",
  LumB  = c("GLYATL2", "ST8SIA2", "PID1", "LRIG3"),
  Basal = c("TRIM59", "LINC01615", "DLL3", "CIP2A")
)[[SUBTYPE]]

genes_file    <- paste0(SUBTYPE, "_Significant.rda")
samples_file  <- paste0(SUBTYPE, "_Samples.rda")
counts_file   <- "BRCA.rds"
clinical_file <- "Clinical_Data.csv"
out_dir       <- paste0("survival_plots_", SUBTYPE)
dir.create(out_dir, showWarnings = FALSE)
set.seed(42)

trim_submitter <- function(x) sub("^([^-]+-[^-]+-[^-]+).*$", "\\1", x)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
Genes   <- unique(as.character(readRDS(genes_file)))
samples <- trim_submitter(gsub("\\.", "-", unique(as.character(readRDS(samples_file)))))
message(SUBTYPE, ": ", length(Genes), " genes, ", length(samples), " samples")

counts_obj <- readRDS(counts_file)
counts_mat <- if (inherits(counts_obj, "SummarizedExperiment")) assay(counts_obj) else as.matrix(counts_obj)
counts_sub <- counts_mat[, trim_submitter(colnames(counts_mat)) %in% samples, drop = FALSE]

clinical <- read.csv(clinical_file, stringsAsFactors = FALSE) %>%
  mutate(submitter_id_trim = trim_submitter(submitter_id),
         deceased = case_when(tolower(vital_status) == "alive" ~ 0,
                              tolower(vital_status) == "dead"  ~ 1,
                              TRUE ~ NA_real_),               # "Not Reported" is excluded
         overall_survival = as.numeric(coalesce(days_to_death, days_to_last_follow_up))) %>%
  filter(submitter_id_trim %in% trim_submitter(colnames(counts_sub))) %>%
  select(submitter_id_trim, deceased, overall_survival) %>%
  drop_na() %>%
  filter(overall_survival > 0) %>%
  distinct(submitter_id_trim, .keep_all = TRUE)
message(nrow(clinical), " patients with follow-up, ", sum(clinical$deceased), " deaths")

# -----------------------------------------------------------------------------
# 2. Variance-stabilizing transformation
# -----------------------------------------------------------------------------
counts_sub <- counts_sub[, match(clinical$submitter_id_trim, trim_submitter(colnames(counts_sub)))]
colnames(counts_sub) <- clinical$submitter_id_trim
counts_sub <- counts_sub[rowSums(counts_sub) >= 10, ]
dds <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                              colData = data.frame(row.names = clinical$submitter_id_trim),
                              design = ~ 1)
counts_vst <- assay(vst(dds, blind = FALSE))

gene_data <- function(gene) {
  data.frame(submitter_id_trim = colnames(counts_vst), expr = counts_vst[gene, ]) %>%
    left_join(clinical, by = "submitter_id_trim") %>%
    drop_na()
}

# -----------------------------------------------------------------------------
# 3. Survival analysis per gene
# -----------------------------------------------------------------------------
results <- list()
for (gene in Genes) {
  if (!gene %in% rownames(counts_vst)) { message("Not in expression data: ", gene); next }
  df <- gene_data(gene)
  if (nrow(df) < 10) next

  cut_result <- surv_cutpoint(df, time = "overall_survival", event = "deceased",
                              variables = "expr", minprop = 0.1)
  cp <- cut_result$cutpoint$cutpoint
  df$strata <- factor(ifelse(df$expr > cp, "HIGH", "LOW"), levels = c("LOW", "HIGH"))   # LOW = reference
  if (length(unique(df$strata)) < 2) next

  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = df)
  km  <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                    title = paste("Overall survival:", gene),
                    subtitle = paste0("Cut-point: ", round(cp, 2)))

  cox_cat <- summary(coxph(Surv(overall_survival, deceased) ~ strata, data = df))$conf.int
  df$expr_sd <- as.numeric(scale(df$expr))
  cox_cont <- summary(coxph(Surv(overall_survival, deceased) ~ expr_sd, data = df))
  ms <- maxstat.test(Surv(overall_survival, deceased) ~ expr, data = df,
                     smethod = "LogRank", pmethod = "Lau92", minprop = 0.1, maxprop = 0.9)

  results[[gene]] <- data.frame(
    gene            = gene,
    n_patients      = nrow(df),
    n_deaths        = sum(df$deceased),
    n_high          = sum(df$strata == "HIGH"),
    deaths_high     = sum(df$deceased[df$strata == "HIGH"]),
    n_low           = sum(df$strata == "LOW"),
    deaths_low      = sum(df$deceased[df$strata == "LOW"]),
    cutoff          = cp,
    HR_high_vs_low  = cox_cat["strataHIGH", "exp(coef)"],
    HR_lower95      = cox_cat["strataHIGH", "lower .95"],
    HR_upper95      = cox_cat["strataHIGH", "upper .95"],
    prognosis       = ifelse(cox_cat["strataHIGH", "exp(coef)"] < 1, "Good", "Poor"),
    logrank_p_value = surv_pvalue(fit)$pval,
    maxstat_p_value = ms$p.value,
    HR_per_SD       = cox_cont$conf.int["expr_sd", "exp(coef)"],
    HR_per_SD_lower = cox_cont$conf.int["expr_sd", "lower .95"],
    HR_per_SD_upper = cox_cont$conf.int["expr_sd", "upper .95"],
    HR_per_SD_p     = cox_cont$coefficients["expr_sd", "Pr(>|z|)"]
  )
  ggsave(file.path(out_dir, paste0("Survival_KM_", gene, ".png")), plot = km$plot, width = 7, height = 6)
}

# -----------------------------------------------------------------------------
# 4. Multiple-testing correction and output
# -----------------------------------------------------------------------------
results_df <- bind_rows(results) %>%
  mutate(logrank_adj_p_value = p.adjust(logrank_p_value, method = "BH"),
         maxstat_adj_p_value = p.adjust(maxstat_p_value, method = "BH"),
         HR_per_SD_adj_p     = p.adjust(HR_per_SD_p, method = "BH"))
write.csv(results_df, paste0(SUBTYPE, "_survival_results.csv"), row.names = FALSE)
message("Significant by log-rank (BH < 0.05): ", sum(results_df$logrank_adj_p_value < 0.05))
message("Significant after cut-point correction (BH < 0.05): ", sum(results_df$maxstat_adj_p_value < 0.05))

# -----------------------------------------------------------------------------
# 5. Final candidates: cut-point-corrected and permutation p-values
# -----------------------------------------------------------------------------
for (g in intersect(final_candidates, rownames(counts_vst))) {
  df <- gene_data(g)
  p_lau  <- maxstat.test(Surv(overall_survival, deceased) ~ expr, data = df, smethod = "LogRank",
                         pmethod = "Lau92", minprop = 0.1, maxprop = 0.9)$p.value
  p_perm <- maxstat.test(Surv(overall_survival, deceased) ~ expr, data = df, smethod = "LogRank",
                         pmethod = "condMC", B = 10000, minprop = 0.1, maxprop = 0.9)$p.value
  p_cont <- summary(coxph(Surv(overall_survival, deceased) ~ scale(expr), data = df))$coefficients[1, 5]
  cat(g, ": Lau92 p =", signif(p_lau, 3), "| permutation p =", signif(p_perm, 3),
      "| continuous Cox p =", signif(p_cont, 3), "\n")
}
