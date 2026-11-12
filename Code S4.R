# =============================================================================
# Code S4 - Final Cox models and figures for the prognostic candidates
#
# For each candidate (F2RL2 in Luminal A; GLYATL2, ST8SIA2, PID1, LRIG3 in Luminal B;
# TRIM59, LINC01615, DLL3, CIP2A in Basal-like):
#   - association of expression with age (Spearman) and stage (Kruskal-Wallis)
#   - optimal cut-point (surv_cutpoint, minimum group proportion 10%); low expression = reference
#   - Cox models on the same patients: unadjusted, + age, + stage, age + stage
#   - proportional hazards (PH) tested with Schoenfeld residuals; pre-specified remedy:
#       age violates -> age as time-varying effect (age x log time)
#       stage violates -> stratification by stage
#   - final model, and the same model with expression as a continuous variable (HR per 1 SD)
#   - number of patients and deaths per expression group; BH adjustment within subtype
#   - ROC curve of continuous expression against vital status (apparent AUC)
# Stage is modeled as I-II vs III-IV.
#
# Input (working directory)
#   BRCA.rds, LumA_Samples.rda, LumB_Samples.rda, Basal_Samples.rda, Clinical_Data.csv
#   Stage_Data.csv (downloaded from GDC with TCGAbiolinks on the first run)
# Output
#   S4_results_stage_binary.csv               all models (Table S38)
#   Table_prognostic_models_stage_binary.csv  final and continuous models
#   S4_plots/Figure_S4 ... Figure_S9 (.png, .tiff)
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(survminer); library(DESeq2); library(tidyverse); library(pROC)
})

counts_file   <- "BRCA.rds"
clinical_file <- "Clinical_Data.csv"
stage_file    <- "Stage_Data.csv"
subtypes <- list(     
  LumA  = list(samples = "LumA_Samples.rda",  genes = "F2RL2"),
  LumB  = list(samples = "LumB_Samples.rda",  genes = c("GLYATL2", "ST8SIA2", "PID1", "LRIG3")),
  Basal = list(samples = "Basal_Samples.rda", genes = c("TRIM59", "LINC01615", "DLL3", "CIP2A"))
)
plot_dir <- "S4_plots"; dir.create(plot_dir, showWarnings = FALSE)
# "binary" = Stage I-II vs III-IV (some four-level stage groups contain no deaths);
# "four"   = Stage I / II / III / IV
stage_mode <- "binary"
trim_submitter <- function(x) sub("^([^-]+-[^-]+-[^-]+).*$", "\\1", x)

# -----------------------------------------------------------------------------
# 1. Clinical data (survival and age from Clinical_Data.csv; stage from GDC)
# -----------------------------------------------------------------------------
if (!file.exists(stage_file)) {
  library(TCGAbiolinks)
  gdc <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
  write.csv(gdc %>% select(submitter_id, ajcc_pathologic_stage), stage_file, row.names = FALSE)
}
stage_raw <- read.csv(stage_file, stringsAsFactors = FALSE)

# substages collapsed: "Stage IIA" -> "Stage II"
stage_roman <- toupper(gsub("STAGE|\\s", "", toupper(stage_raw$ajcc_pathologic_stage)))
stage_roman <- sub("[ABC]$", "", stage_roman)
stage_raw$stage <- factor(ifelse(stage_roman %in% c("I", "II", "III", "IV"),
                                 paste("Stage", stage_roman), NA),
                          levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

clinical <- read.csv(clinical_file, stringsAsFactors = FALSE) %>%
  mutate(submitter_id_trim = trim_submitter(submitter_id),
         deceased = case_when(tolower(vital_status) == "alive" ~ 0,
                              tolower(vital_status) == "dead"  ~ 1, TRUE ~ NA_real_),
         overall_survival = as.numeric(coalesce(days_to_death, days_to_last_follow_up))) %>%
  left_join(stage_raw %>% transmute(submitter_id_trim = trim_submitter(submitter_id), stage),
            by = "submitter_id_trim") %>%
  select(submitter_id_trim, deceased, overall_survival, age, stage) %>%
  distinct(submitter_id_trim, .keep_all = TRUE)

counts_obj <- readRDS(counts_file)
counts_mat <- if (inherits(counts_obj, "SummarizedExperiment")) assay(counts_obj) else as.matrix(counts_obj)

# -----------------------------------------------------------------------------
# 2. Models for one gene
# -----------------------------------------------------------------------------
hr_of <- function(fit, term) {
  ci <- summary(fit)$conf.int; co <- summary(fit)$coefficients
  i <- grep(term, rownames(ci))[1]
  tibble(HR = ci[i, "exp(coef)"], lower95 = ci[i, "lower .95"], upper95 = ci[i, "upper .95"],
         p = co[i, ncol(co)])
}
tt_age <- function(x, t, ...) x * log(t + 1)

validate_gene <- function(gene, expr_vst, clin, subtype) {
  df <- data.frame(submitter_id_trim = colnames(expr_vst), expr = expr_vst[gene, ]) %>%
    inner_join(clin, by = "submitter_id_trim") %>%
    rename(OS_time = overall_survival, OS_status = deceased) %>%
    filter(!is.na(expr), !is.na(age), !is.na(stage), !is.na(OS_time), !is.na(OS_status), OS_time > 0) %>%
    mutate(stage = droplevels(stage))
  if (stage_mode == "binary")
    df$stage <- factor(ifelse(df$stage %in% c("Stage I", "Stage II"), "I-II", "III-IV"),
                       levels = c("I-II", "III-IV"))
  
  rho_age <- cor(df$expr, df$age, method = "spearman")
  p_age   <- cor.test(df$expr, df$age, method = "spearman", exact = FALSE)$p.value
  p_stage <- kruskal.test(expr ~ stage, data = df)$p.value
  
  cut <- surv_cutpoint(df, time = "OS_time", event = "OS_status", variables = "expr", minprop = 0.1)
  df$Expression <- relevel(factor(surv_categorize(cut)$expr), ref = "low")
  
  m_full <- coxph(Surv(OS_time, OS_status) ~ Expression + age + stage, data = df)
  zph <- cox.zph(m_full)$table
  age_viol <- zph["age", "p"] < 0.05; stage_viol <- zph["stage", "p"] < 0.05
  
  f_final <- if (age_viol && stage_viol) Surv(OS_time, OS_status) ~ Expression + age + tt(age) + strata(stage) else
    if (age_viol)               Surv(OS_time, OS_status) ~ Expression + age + tt(age) + stage else
      if (stage_viol)             Surv(OS_time, OS_status) ~ Expression + age + strata(stage) else
        Surv(OS_time, OS_status) ~ Expression + age + stage
  label <- c("conventional", "tt(age)", "strata(stage)", "tt(age) + strata(stage)")[1 + age_viol + 2 * stage_viol]
  fit_with <- function(f) if (age_viol) coxph(f, data = df, tt = tt_age) else coxph(f, data = df)
  m_final <- fit_with(f_final)
  
  m_unadj <- coxph(Surv(OS_time, OS_status) ~ Expression, data = df)
  m_age   <- coxph(Surv(OS_time, OS_status) ~ Expression + age, data = df)
  m_stage <- coxph(Surv(OS_time, OS_status) ~ Expression + stage, data = df)
  b0 <- coef(m_unadj)["Expressionhigh"]
  pct <- function(m) 100 * (coef(m)["Expressionhigh"] - b0) / b0
  
  df$expr_sd <- as.numeric(scale(df$expr))
  m_cont <- fit_with(update(f_final, . ~ . - Expression + expr_sd))
  
  roc_obj <- roc(df$OS_status, df$expr, quiet = TRUE)
  auc <- as.numeric(ci(roc_obj))
  fig_store[[gene]] <<- list(m_final = m_final, label = label, roc = roc_obj, auc = auc,
                             n = nrow(df), deaths = sum(df$OS_status))
  
  # PH test of the final model (not available for models with a tt() term)
  zph_final <- tryCatch(cox.zph(m_final)$table, error = function(e) NULL)
  ph_expr_final <- if (!is.null(zph_final) && "Expression" %in% rownames(zph_final))
    zph_final["Expression", "p"] else NA_real_
  if (zph["Expression", "p"] < 0.05)
    message("  PH assumption violated for ", gene, " expression (p = ", signif(zph["Expression", "p"], 2), ")")
  
  bind_rows(
    hr_of(m_unadj, "Expressionhigh") %>% mutate(model = "1 unadjusted (same patients)"),
    hr_of(m_age,   "Expressionhigh") %>% mutate(model = "2 + age only",   pct_change_logHR = pct(m_age)),
    hr_of(m_stage, "Expressionhigh") %>% mutate(model = "3 + stage only", pct_change_logHR = pct(m_stage)),
    hr_of(m_full,  "Expressionhigh") %>% mutate(model = "4 age + stage",  pct_change_logHR = pct(m_full)),
    hr_of(m_final, "Expressionhigh") %>% mutate(model = paste("5 FINAL:", label)),
    hr_of(m_cont,  "expr_sd")        %>% mutate(model = paste("6 continuous per SD:", label))
  ) %>%
    mutate(subtype = subtype, gene = gene, n_patients = nrow(df), n_deaths = sum(df$OS_status),
           n_high = sum(df$Expression == "high"), deaths_high = sum(df$OS_status[df$Expression == "high"]),
           n_low  = sum(df$Expression == "low"),  deaths_low  = sum(df$OS_status[df$Expression == "low"]),
           PH_p_age = zph["age", "p"], PH_p_stage = zph["stage", "p"], PH_p_expr = zph["Expression", "p"],
           PH_p_expr_final_model = ph_expr_final, PH_expr_violated = zph["Expression", "p"] < 0.05,
           rho_age = rho_age, p_age = p_age, p_stage = p_stage,
           AUC = auc[2], AUC_lower = auc[1], AUC_upper = auc[3], .before = 1)
}

# -----------------------------------------------------------------------------
# 3. Loop over subtypes and genes (VST as in Code S2)
# -----------------------------------------------------------------------------
all_res <- list(); fig_store <- list()
for (st in names(subtypes)) {
  samp <- trim_submitter(gsub("\\.", "-", unique(as.character(readRDS(subtypes[[st]]$samples)))))
  clin_st <- clinical %>% filter(submitter_id_trim %in% samp, !is.na(deceased), !is.na(overall_survival),
                                 overall_survival > 0)
  cm <- counts_mat[, trim_submitter(colnames(counts_mat)) %in% clin_st$submitter_id_trim, drop = FALSE]
  cm <- cm[, !duplicated(trim_submitter(colnames(cm))), drop = FALSE]
  colnames(cm) <- trim_submitter(colnames(cm))
  cm <- cm[rowSums(cm) >= 10, ]
  vsd <- vst(DESeqDataSetFromMatrix(round(cm), data.frame(row.names = colnames(cm), id = colnames(cm)), ~ 1),
             blind = FALSE)
  expr_vst <- assay(vsd)
  message(st, ": ", ncol(expr_vst), " patients, ", sum(clin_st$deceased), " deaths, ",
          sum(is.na(clin_st$stage)), " without stage (excluded from the Cox models)")
  print(table(Stage = clin_st$stage, Died = clin_st$deceased))   # deaths per stage group
  for (g in subtypes[[st]]$genes) {
    if (!g %in% rownames(expr_vst)) { message("  ", g, " not found"); next }
    all_res[[paste(st, g)]] <- validate_gene(g, expr_vst, clin_st, st)
    message("  done: ", g)
  }
}
res <- bind_rows(all_res) %>%
  mutate(model_type = sub(":.*", "", model)) %>%
  group_by(subtype, model_type) %>% mutate(p_BH_within_subtype = p.adjust(p, method = "BH")) %>%
  ungroup() %>% select(-model_type)
write_csv(res, paste0("S4_results_stage_", stage_mode, ".csv"))

# final and continuous models only
table_s <- res %>% filter(grepl("^5 FINAL|^6 continuous", model)) %>%
  transmute(Subtype = subtype, Gene = gene, Model = sub("^[0-9] ", "", model),
            `Patients (deaths)` = paste0(n_patients, " (", n_deaths, ")"),
            `High: patients (deaths)` = paste0(n_high, " (", deaths_high, ")"),
            `Low: patients (deaths)`  = paste0(n_low, " (", deaths_low, ")"),
            `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, lower95, upper95),
            p = signif(p, 3), `BH-adjusted p (within subtype)` = signif(p_BH_within_subtype, 3),
            `PH p (expression)` = signif(PH_p_expr, 3))
write_csv(table_s, paste0("Table_prognostic_models_stage_", stage_mode, ".csv"))

# console summary
print(res %>% select(gene, model, HR, lower95, upper95, p, p_BH_within_subtype, n_patients, n_deaths,
                     deaths_high, n_high, deaths_low, n_low), n = Inf, width = Inf)
print(res %>% filter(grepl("^5 FINAL", model)) %>%
        select(gene, PH_p_expr, PH_p_expr_final_model, PH_expr_violated))
# -----------------------------------------------------------------------------
# 4. Figures S4-S9: forest plots of the final models and ROC curves
# -----------------------------------------------------------------------------
library(ggpubr)
term_label <- function(x) dplyr::case_when(
  grepl("^Expression", x) ~ "Expression: high vs low",
  grepl("^tt\\(age", x)    ~ "Age \u00d7 log(time)",
  x == "age"              ~ "Age (per year)",
  grepl("^stage", x)      ~ "Stage III\u2013IV vs I\u2013II",
  TRUE ~ x)
forest_plot <- function(gene) {
  f <- fig_store[[gene]]; s <- summary(f$m_final)
  d <- tibble(term = rownames(s$conf.int), HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
              p = s$coefficients[, ncol(s$coefficients)]) %>%
    mutate(label = term_label(term),
           txt = sprintf("%.2f (%.2f\u2013%.2f)  p %s", HR, lo, hi,
                         ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, format = "f", digits = 3)))),
           label = factor(label, levels = rev(label)))
  model_txt <- c(conventional = "adjusted for age and stage", `tt(age)` = "age as time-varying effect; adjusted for stage",
                 `strata(stage)` = "adjusted for age; stratified by stage",
                 `tt(age) + strata(stage)` = "age as time-varying effect; stratified by stage")[f$label]
  rng <- c(min(d$lo), max(d$hi)); rng <- c(min(rng[1], 0.8), max(rng[2], 1.25))
  ggplot(d, aes(x = HR, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_point(size = 2.6, shape = 15) +
    geom_text(aes(label = txt), x = Inf, hjust = -0.08, size = 3.1) +
    scale_x_log10(limits = c(rng[1] / 1.3, rng[2] * 1.3)) +
    coord_cartesian(clip = "off") +
    labs(title = gene, subtitle = sprintf("%s; %d patients, %d deaths", model_txt, f$n, f$deaths),
         x = "Hazard ratio (log scale)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold.italic"), plot.subtitle = element_text(size = 8.5),
          plot.margin = margin(5.5, 175, 5.5, 5.5))
}
roc_plot <- function(gene) {
  f <- fig_store[[gene]]
  ggroc(f$roc, legacy.axes = TRUE, linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    annotate("text", x = 0.98, y = 0.06, hjust = 1, size = 3.4,
             label = sprintf("AUC = %.2f (95%% CI %.2f\u2013%.2f)", f$auc[2], f$auc[1], f$auc[3])) +
    coord_equal() + labs(title = gene, x = "False positive rate (1 \u2212 specificity)", y = "True positive rate (sensitivity)") +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold.italic"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 300, bg = "white",
         compression = "lzw")
}
panel <- function(genes, fun, ncol = 2) ggarrange(plotlist = lapply(genes, fun), ncol = ncol,
                                                  nrow = ceiling(length(genes) / ncol),
                                                  labels = letters[seq_along(genes)])
save_fig(forest_plot("F2RL2"), "Figure_S4", 8, 3.2)
save_fig(roc_plot("F2RL2"), "Figure_S5", 5, 5)
save_fig(panel(c("GLYATL2", "ST8SIA2", "PID1", "LRIG3"), forest_plot), "Figure_S6", 15, 7)
save_fig(panel(c("GLYATL2", "LRIG3", "PID1", "ST8SIA2"), roc_plot), "Figure_S7", 10, 10)
save_fig(panel(c("TRIM59", "LINC01615", "DLL3", "CIP2A"), forest_plot), "Figure_S8", 15, 7)
save_fig(panel(c("DLL3", "LINC01615", "TRIM59", "CIP2A"), roc_plot), "Figure_S9", 10, 10)
message("Figures S4-S9 saved in ", plot_dir)