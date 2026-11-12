# =============================================================================
# Code S11 - Sensitivity of the final Cox models to treatment and receptor status
#
# Each final model of Code S4 (same cut-point and PH remedy) is refitted with one
# additional covariate at a time: radiotherapy, systemic drug therapy, both treatments,
# ER, PR, HER2 status, or triple-negative status. Patients with a missing covariate are
# excluded, and the model without the covariate is refitted on the same patients for
# comparison. A covariate is not fitted if its smaller group has fewer than 10 patients.
#
# Input (working directory): BRCA.rds, LumA/LumB/Basal_Samples.rda, Clinical_Data.csv,
#   Stage_Data.csv (from Code S4). Treatment (GDC clinical) and receptor status (BCR Biotab)
#   are downloaded with TCGAbiolinks on the first run and saved as Treatment_Receptor_Data.csv.
# Output: Treatment_receptor_sensitivity.csv (Table S42)
# =============================================================================
suppressPackageStartupMessages({ library(survival); library(survminer); library(DESeq2); library(tidyverse) })

counts_file <- "BRCA.rds"; clinical_file <- "Clinical_Data.csv"; stage_file <- "Stage_Data.csv"
tr_file <- "Treatment_Receptor_Data.csv"
subtypes <- list(
  LumA  = list(samples = "LumA_Samples.rda",  genes = "F2RL2"),
  LumB  = list(samples = "LumB_Samples.rda",  genes = c("GLYATL2", "ST8SIA2", "PID1", "LRIG3")),
  Basal = list(samples = "Basal_Samples.rda", genes = c("TRIM59", "LINC01615", "DLL3", "CIP2A"))
)
GDC_DIR <- if (.Platform$OS.type == "windows") "C:/gdc" else "GDCdata"  # short path (Windows path-length limit)
MIN_GROUP <- 10   # a binary covariate is used only if its smaller group has >= 10 patients
trim_submitter <- function(x) sub("^([^-]+-[^-]+-[^-]+).*$", "\\1", x)

# -----------------------------------------------------------------------------
# 1. Treatment and receptor status
# -----------------------------------------------------------------------------
if (!file.exists(tr_file)) {
  library(TCGAbiolinks)
  yn <- function(x) { x <- tolower(as.character(x))
  ifelse(grepl("yes", x), 1, ifelse(grepl("(^|[^a-z])no([^a-z]|$)", x), 0, NA)) }
  gdc <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
  pc <- grep("pharmaceutical_treatment_or_therapy", names(gdc), value = TRUE)[1]
  rc <- grep("radiation_treatment_or_therapy", names(gdc), value = TRUE)[1]
  cat("Treatment columns used:", pc, "|", rc, "\n")
  tr <- tibble(submitter_id = gdc$submitter_id,
               systemic_therapy = if (is.na(pc)) NA else yn(gdc[[pc]]),
               radiotherapy     = if (is.na(rc)) NA else yn(gdc[[rc]]))
  
  q <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical",
                data.type = "Clinical Supplement", data.format = "BCR Biotab")
  dir.create(GDC_DIR, showWarnings = FALSE, recursive = TRUE)
  GDCdownload(q, directory = GDC_DIR, method = "api", files.per.chunk = 5)
  pat <- tryCatch(GDCprepare(q, directory = GDC_DIR)$clinical_patient_brca, error = function(e) {
    # fallback: read the downloaded patient table directly
    f <- list.files(GDC_DIR, pattern = "clinical_patient_brca\\.txt$", recursive = TRUE, full.names = TRUE)
    if (!length(f)) stop("nationwidechildrens.org_clinical_patient_brca.txt not found in ", GDC_DIR)
    read.delim(f[1], stringsAsFactors = FALSE, check.names = FALSE)
  })
  pat <- pat %>% filter(grepl("^TCGA-", bcr_patient_barcode))
  cat("Receptor status found for", nrow(pat), "patients\n")
  pn <- function(x) { x <- tolower(as.character(x)); ifelse(x == "positive", 1, ifelse(x == "negative", 0, NA)) }
  ihc  <- tolower(pat$her2_status_by_ihc)
  fish <- if ("her2_fish_status" %in% names(pat)) tolower(pat$her2_fish_status) else rep(NA, nrow(pat))
  her2 <- ifelse(ihc %in% "positive" | fish %in% "positive", 1,
                 ifelse(ihc %in% "negative" | fish %in% "negative", 0, NA))   # IHC positive or FISH positive
  rec <- tibble(submitter_id = pat$bcr_patient_barcode,
                ER = pn(pat$er_status_by_ihc), PR = pn(pat$pr_status_by_ihc), HER2 = her2)
  write.csv(full_join(tr, rec, by = "submitter_id"), tr_file, row.names = FALSE)
}
tr <- read.csv(tr_file) %>%
  transmute(submitter_id_trim = trim_submitter(submitter_id), radiotherapy, systemic_therapy, ER, PR, HER2,
            TNBC = ifelse(ER == 0 & PR == 0 & HER2 == 0, 1,
                          ifelse(ER == 1 | PR == 1 | HER2 == 1, 0, NA))) %>%
  distinct(submitter_id_trim, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# 2. Survival, age and stage (as in Code S4)
# -----------------------------------------------------------------------------
stage_raw <- read.csv(stage_file, stringsAsFactors = FALSE)
stage_roman <- sub("[ABC]$", "", toupper(gsub("STAGE|\\s", "", toupper(stage_raw$ajcc_pathologic_stage))))
stage_raw$stage <- factor(ifelse(stage_roman %in% c("I", "II", "III", "IV"), paste("Stage", stage_roman), NA),
                          levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
clinical <- read.csv(clinical_file, stringsAsFactors = FALSE) %>%
  mutate(submitter_id_trim = trim_submitter(submitter_id),
         deceased = case_when(tolower(vital_status) == "alive" ~ 0, tolower(vital_status) == "dead" ~ 1),
         overall_survival = as.numeric(coalesce(days_to_death, days_to_last_follow_up))) %>%
  left_join(stage_raw %>% transmute(submitter_id_trim = trim_submitter(submitter_id), stage),
            by = "submitter_id_trim") %>%
  select(submitter_id_trim, deceased, overall_survival, age, stage) %>%
  distinct(submitter_id_trim, .keep_all = TRUE) %>%
  left_join(tr, by = "submitter_id_trim")

counts_obj <- readRDS(counts_file)
counts_mat <- if (inherits(counts_obj, "SummarizedExperiment")) assay(counts_obj) else as.matrix(counts_obj)
tt_age <- function(x, t, ...) x * log(t + 1)
hr_of <- function(fit, term) {
  ci <- summary(fit)$conf.int; co <- summary(fit)$coefficients; i <- grep(term, rownames(ci))[1]
  c(HR = ci[i, "exp(coef)"], lower95 = ci[i, "lower .95"], upper95 = ci[i, "upper .95"], p = co[i, ncol(co)])
}

# -----------------------------------------------------------------------------
# 3. Final model + one covariate at a time
# -----------------------------------------------------------------------------
sens_gene <- function(gene, expr_vst, clin, subtype) {
  df <- data.frame(submitter_id_trim = colnames(expr_vst), expr = expr_vst[gene, ]) %>%
    inner_join(clin, by = "submitter_id_trim") %>%
    rename(OS_time = overall_survival, OS_status = deceased) %>%
    filter(!is.na(expr), !is.na(age), !is.na(stage), !is.na(OS_time), !is.na(OS_status), OS_time > 0)
  df$stage <- factor(ifelse(df$stage %in% c("Stage I", "Stage II"), "I-II", "III-IV"), levels = c("I-II", "III-IV"))
  
  cut <- surv_cutpoint(df, time = "OS_time", event = "OS_status", variables = "expr", minprop = 0.1)
  df$Expression <- relevel(factor(surv_categorize(cut)$expr), ref = "low")
  zph <- cox.zph(coxph(Surv(OS_time, OS_status) ~ Expression + age + stage, data = df))$table
  age_viol <- zph["age", "p"] < 0.05; stage_viol <- zph["stage", "p"] < 0.05
  f_final <- if (age_viol && stage_viol) Surv(OS_time, OS_status) ~ Expression + age + tt(age) + strata(stage) else
    if (age_viol)               Surv(OS_time, OS_status) ~ Expression + age + tt(age) + stage else
      if (stage_viol)             Surv(OS_time, OS_status) ~ Expression + age + strata(stage) else
        Surv(OS_time, OS_status) ~ Expression + age + stage
  label <- c("conventional", "tt(age)", "strata(stage)", "tt(age) + strata(stage)")[1 + age_viol + 2 * stage_viol]
  fit_with <- function(f, d) if (age_viol) coxph(f, data = d, tt = tt_age) else coxph(f, data = d)
  
  sets <- list("Radiotherapy" = "radiotherapy", "Systemic therapy" = "systemic_therapy",
               "Radiotherapy + systemic therapy" = c("radiotherapy", "systemic_therapy"),
               "ER status" = "ER", "PR status" = "PR", "HER2 status" = "HER2", "Triple-negative status" = "TNBC")
  out <- list(as_tibble(as.list(hr_of(fit_with(f_final, df), "Expressionhigh"))) %>%
                mutate(covariate = "None (final model, Code S4)", n = nrow(df), deaths = sum(df$OS_status),
                       note = ""))
  for (nm in names(sets)) {
    cs <- sets[[nm]]
    d2 <- df[complete.cases(df[, cs, drop = FALSE]), ]
    small <- sapply(cs, function(v) min(table(factor(d2[[v]], levels = 0:1))))
    if (nrow(d2) == 0 || any(small < MIN_GROUP)) {
      out[[nm]] <- tibble(covariate = nm, n = nrow(d2), deaths = sum(d2$OS_status),
                          note = paste0("not fitted: smaller group < ", MIN_GROUP, " patients"))
      next
    }
    base <- tryCatch(hr_of(fit_with(f_final, d2), "Expressionhigh"), error = function(e) NULL)
    adj  <- tryCatch(hr_of(fit_with(update(f_final, as.formula(paste(". ~ . +", paste(cs, collapse = " + ")))), d2),
                           "Expressionhigh"), error = function(e) NULL)
    if (is.null(base) || is.null(adj)) { out[[nm]] <- tibble(covariate = nm, note = "model failed"); next }
    out[[nm]] <- as_tibble(as.list(adj)) %>%
      mutate(covariate = nm, n = nrow(d2), deaths = sum(d2$OS_status),
             HR_base_same_patients = base[["HR"]], p_base_same_patients = base[["p"]],
             pct_change_logHR = 100 * (log(adj[["HR"]]) - log(base[["HR"]])) / log(base[["HR"]]),
             note = ifelse(adj[["upper95"]] > 1e3 | adj[["lower95"]] < 1e-3, "unstable estimate (sparse data)", ""))
  }
  bind_rows(out) %>% mutate(subtype = subtype, gene = gene, final_model = label, .before = 1)
}

# -----------------------------------------------------------------------------
# 4. Loop over subtypes and genes (VST as in Codes S2 and S4)
# -----------------------------------------------------------------------------
res <- list()
for (st in names(subtypes)) {
  samp <- trim_submitter(gsub("\\.", "-", unique(as.character(readRDS(subtypes[[st]]$samples)))))
  clin_st <- clinical %>% filter(submitter_id_trim %in% samp, !is.na(deceased), !is.na(overall_survival),
                                 overall_survival > 0)
  cat("\n====", st, "- covariate availability (0 / 1 / missing) ====\n")
  for (v in c("radiotherapy", "systemic_therapy", "ER", "PR", "HER2", "TNBC"))
    cat(sprintf("  %-17s %s\n", v, paste(table(factor(clin_st[[v]], levels = 0:1), useNA = "always"), collapse = " / ")))
  cm <- counts_mat[, trim_submitter(colnames(counts_mat)) %in% clin_st$submitter_id_trim, drop = FALSE]
  cm <- cm[, !duplicated(trim_submitter(colnames(cm))), drop = FALSE]
  colnames(cm) <- trim_submitter(colnames(cm))
  cm <- cm[rowSums(cm) >= 10, ]
  expr_vst <- assay(vst(DESeqDataSetFromMatrix(round(cm), data.frame(row.names = colnames(cm)), ~ 1), blind = FALSE))
  for (g in subtypes[[st]]$genes) if (g %in% rownames(expr_vst)) res[[paste(st, g)]] <- sens_gene(g, expr_vst, clin_st, st)
}
res <- bind_rows(res)
write_csv(res, "Treatment_receptor_sensitivity.csv")
print(res %>% transmute(gene, covariate, n, deaths, HR = round(HR, 2), CI = sprintf("%.2f-%.2f", lower95, upper95),
                        p = signif(p, 2), HR_base = round(HR_base_same_patients, 2),
                        pct_change = round(pct_change_logHR), note), n = Inf, width = Inf)