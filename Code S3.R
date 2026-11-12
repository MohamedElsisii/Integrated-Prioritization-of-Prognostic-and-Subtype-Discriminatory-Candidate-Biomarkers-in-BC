# ==================================
# Script: Code S3 - Univariate & Multivariate Cox.R
#
# Purpose:
# This script performs survival analysis on a list of candidate genes
# (identified in 'Code S2') for a specific breast cancer subtype.
# It runs two types of Cox Proportional Hazards regression:
# 1. Univariate (Uni-Cox): Tests the prognostic power of each gene's
#    expression ALONE.
# 2. Multivariate (Multi-Cox): Tests the gene's prognostic power
#    WHILE CONTROLLING for key clinical covariates (Age and Stage).
#
# Output:
# - A summary CSV comparing Uni- and Multi-Cox results for all genes.
# - Individual forest plots for both models for each gene.
# - A full summary table for each gene's multivariate model.
# ==================================

# Suppress startup messages for a cleaner console
suppressPackageStartupMessages({
  library(survival)    # Core package for survival analysis (Surv, coxph)
  library(survminer)   # For plotting survival models (ggforest)
  library(DESeq2)      # For VST normalization
  library(tidyverse)   # For data manipulation (dplyr, tidyr, etc.)
  library(broom)       # For tidying model outputs (tidy)
  library(ggplot2)     # For saving plots (ggsave)
})

# ----------------------------
# User options
# ----------------------------
# --- Input Files ---
# File from 'Code S2', containing significant genes from Kaplan-Meier
gene_file     <- "LumB_Survival_Sig.rda"    
# File from 'Code S1', containing sample IDs for this subtype
samples_file  <- "LumB_Samples.rda"         
# The raw TCGA RNA-seq counts data (RDS format)
counts_file   <- "BRCA.rds"                 
# The master clinical data file (CSV format)
clinical_file <- "Clinical_Data.csv"      

# --- Output Settings ---
# Directory to save all plots and result tables
out_dir       <- "Cox_results_LumB"         
# Subtype name used for file naming
subtype_name  <- "LumB"                     

# Create the output directory if it doesn't already exist
dir.create(out_dir, showWarnings = FALSE)

# ----------------------------
# 1. Helper Functions
# ----------------------------

#' @title Trim TCGA Submitter ID
#' @description Standardizes TCGA barcodes (e.g., 'TCGA.A1.A0SD.01A' -> 'TCGA-A1-A0SD')
#' @param x A character vector of TCGA barcodes
#' @return A character vector of trimmed barcodes
trim_submitter <- function(x) {
  sub("^([^-]+-[^-]+-[^-]+).*$","\\1", gsub("\\.", "-", x))
}

#' @title Make Expression DataFrame for Modeling
#' @description Merges a single gene's VST expression with clinical data and
#' filters for complete survival data and specified covariates.
#' @param gene The gene symbol (e.g., "F2RL2")
#' @param expr_vst The VST-normalized expression matrix
#' @param clinical_df The clinical data frame
#' @param covariates A vector of covariate names (e.g., c("age_at_diagnosis"))
#'                   If NULL, only filters for survival data (for univariate).
#'                   If provided, filters for complete cases (no NAs)
#'                   for all specified covariates (for multivariate).
#' @return A data frame ready for modeling, or NULL if data is insufficient.
make_expr_df <- function(gene, expr_vst, clinical_df, covariates = NULL) {
  # Check if gene exists in the expression matrix
  if (!gene %in% rownames(expr_vst)) return(NULL)
  
  # 1. Create a base data frame with gene expression
  df <- data.frame(
    submitter_id_trim = colnames(expr_vst),
    counts = as.numeric(expr_vst[gene, ]) # 'counts' is now our VST expression
  ) %>%
    # 2. Join with all clinical data
    left_join(clinical_df, by = "submitter_id_trim") %>%
    # 3. Filter for valid survival data (required for all models)
    filter(!is.na(overall_survival), !is.na(deceased), overall_survival > 0)
  
  # 4. If covariates are specified (for multivariate model),
  #    drop all rows that have an NA in ANY of those columns.
  if (!is.null(covariates)) {
    df <- df %>% drop_na(all_of(covariates))
  }
  
  # 5. Check if enough data remains to build a model
  if (nrow(df) < 10) return(NULL) # Not enough samples
  return(df)
}

#' @title Plot Univariate Forest Plot
#' @description Creates and saves a ggforest plot for a univariate model
#' @param fit A coxph model object
#' @param dat The data frame used to fit the model
#' @param gene_name The name of the gene
#' @param out_dir The output directory
plot_univariate_with_ggforest <- function(fit, dat, gene_name, out_dir) {
  # Skip if the model failed
  if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
  
  # Define file name
  fname <- file.path(out_dir, paste0("Univariate_Cox_Regression_Plot_", gene_name, ".png"))
  
  # Use tryCatch to prevent a single plot failure from stopping the loop
  tryCatch({
    gf <- ggforest(fit, data = dat, main = paste("Univariate Cox Regression for", gene_name))
    ggsave(filename = fname, plot = gf, width = 7, height = 5, dpi = 300)
    return(TRUE)
  }, error = function(e) {
    message("Failed to plot ", gene_name, ": ", e$message)
    return(NULL)
  })
}

#' @title Extract Roman Numeral Stage
#' @description Cleans messy TCGA stage data ("Stage IIIA", "STAGE IV")
#' into standardized Roman numerals ("I", "II", "III", "IV").
#' @param x A character vector of stage information
#' @return A character vector of standardized Roman numerals (I, II, III, IV) or NA
extract_stage_roman <- function(x) {
  s <- toupper(as.character(x)) # To uppercase
  s[s == ""] <- NA_character_   # Handle empty strings
  
  # Remove "STAGE", "STG", dots, colons, and spaces
  s_clean <- gsub("STAGE|STG\\.?|\\.|\\:|\\s+", "", s)
  
  # Pattern to find Roman numerals IV, III, II, or I
  pat <- "(IV|III|II|I)"
  m <- regmatches(s_clean, regexpr(pat, s_clean, perl = TRUE))
  
  m[m == ""] <- NA_character_                 # If no match, set to NA
  m[! m %in% c("I","II","III","IV")] <- NA_character_ # Ensure only valid stages
  
  # Safety check to ensure output length matches input length
  if (length(m) != length(x)) {
    m2 <- rep(NA_character_, length(x))
    m2[seq_along(m)] <- m
    return(m2)
  }
  return(m)
}

# ----------------------------
# 2. Data Loading & Preprocessing
# ----------------------------

message("Loading and preparing data...")

# Load the list of significant genes from the Kaplan-Meier analysis (Code S2)
gene_names_obj <- readRDS(gene_file)

# Handle if the saved object is a data frame or a simple vector
if (is.data.frame(gene_names_obj) || is.matrix(gene_names_obj)) {
  gene_names <- unique(as.character(gene_names_obj[,1]))
} else {
  gene_names <- unique(as.character(gene_names_obj))
}
message("Total genes to test: ", length(gene_names))

# Load the sample IDs for the specific subtype (e.g., LumB)
Subtype_samples <- readRDS(samples_file) %>% unique() %>% as.character()
# Trim sample IDs to match clinical data format
Subtype_samples_trim <- trim_submitter(Subtype_samples)

# Load the full raw count matrix
counts_obj <- readRDS(counts_file)
# Check if it's a SummarizedExperiment object or a plain matrix
if (inherits(counts_obj, "SummarizedExperiment")) {
  counts_mat <- assay(counts_obj, "counts") # Extract counts
} else {
  counts_mat <- as.matrix(counts_obj)
}

# Load and prepare master clinical data
clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
clinical$submitter_id_trim <- trim_submitter(clinical$submitter_id)

# Create binary survival status (0 = Alive, 1 = Deceased)
clinical$deceased <- ifelse(tolower(clinical$vital_status) == "alive", 0, 1)

# Create overall survival time (in days)
# coalesce finds the first non-NA value (uses days_to_death if available,
# otherwise falls back to days_to_last_follow_up)
clinical$overall_survival <- coalesce(clinical$days_to_death, clinical$days_to_last_follow_up)

# --- Aligning Samples (CRITICAL STEP) ---
# We need to ensure the counts matrix and clinical data
# have the exact same samples, in the exact same order.

# 1. Filter clinical data to ONLY our subtype's samples
clinical_subtype <- clinical %>% filter(submitter_id_trim %in% Subtype_samples_trim)

# 2. Find samples that are in BOTH the counts matrix AND our subtype clinical file
common_samples <- intersect(trim_submitter(colnames(counts_mat)), clinical_subtype$submitter_id_trim)

if (length(common_samples) == 0) {
  stop("No common samples between clinical data and counts matrix.")
}

# 3. Filter and re-order clinical data to match the common list
clinical_subtype <- clinical_subtype %>%  
  filter(submitter_id_trim %in% common_samples) %>%
  arrange(submitter_id_trim) # Sort alphabetically

# 4. Filter and re-order counts matrix to match the clinical data
#    This `match` operation is key to aligning them perfectly.
counts_subtype <- counts_mat[, match(clinical_subtype$submitter_id_trim, trim_submitter(colnames(counts_mat))), drop = FALSE]

# 5. Create the 'coldata' (column metadata) required by DESeq2
coldata <- data.frame(row.names = clinical_subtype$submitter_id_trim, sample = clinical_subtype$submitter_id_trim)

# 6. Standardize the column names of the final counts matrix
colnames(counts_subtype) <- trim_submitter(colnames(counts_subtype))

message("Aligned samples for ", subtype_name, ": ", ncol(counts_subtype))

# ----------------------------
# 3. VST Transformation
# ----------------------------
# We normalize the raw counts to make expression values
# more comparable across samples and genes for the Cox model.
# VST (Variance Stabilizing Transformation) is good for this.

# Create a DESeqDataSet object (design = ~1 means no differential analysis)
dds <- DESeqDataSetFromMatrix(countData = round(counts_subtype), colData = coldata, design = ~1)
# Pre-filter for rows with at least 10 counts total
dds <- dds[rowSums(counts(dds)) >= 10, ]
# Apply the VST
vsd <- vst(dds, blind = FALSE)
# Extract the normalized expression matrix
expr_vst <- assay(vsd)

# ===================================================================
# 4. DEFINE AND PREPARE COVARIATES (Age and Stage Only)
# ===================================================================
# Here we prepare the clinical variables for the multivariate model

# --- Define the user-specified core clinical covariates ---
clinical_covariates <- c("age_at_diagnosis", "ajcc_pathologic_stage")  

# compute roman vector outside mutate to avoid length/evaluation issues
ajcc_stage_roman_vec <- extract_stage_roman(clinical_subtype$ajcc_pathologic_stage)

# now mutate using that precomputed vector
clinical_subtype_multi <- clinical_subtype %>%
  mutate(
    # Convert age from days to years
    age_at_diagnosis = suppressWarnings(as.numeric(age_at_diagnosis) / 365.25),
    
    # Use the pre-cleaned roman numeral vector
    ajcc_stage_roman = ajcc_stage_roman_vec,        
    
    # Convert "I", "II" etc. into clean factors ("Stage I", "Stage II")
    ajcc_pathologic_stage = case_when(
      ajcc_stage_roman == "I"   ~ "Stage I",
      ajcc_stage_roman == "II"  ~ "Stage II",
      ajcc_stage_roman == "III" ~ "Stage III",
      ajcc_stage_roman == "IV"  ~ "Stage IV",
      TRUE ~ NA_character_ # All others become NA
    ),
    
    # Set the factor levels to ensure "Stage I" is the reference
    ajcc_pathologic_stage = factor(ajcc_pathologic_stage, levels = c("Stage I","Stage II","Stage III","Stage IV"))
  )

# Define the base formula for the multivariable model
# 'counts' will be the gene expression variable
multi_formula_base <- paste("Surv(overall_survival, deceased) ~ counts +", paste(clinical_covariates, collapse = " + "))

# ===================================================================
# 5. RUN UNIVARIATE AND MULTIVARIATE COX REGRESSION
# ===================================================================

message("Running Uni- and Multivariable Cox models on ALL genes with Age and Stage adjustment...")

# Create empty lists to store the model results for each gene
uni_list <- list()
multi_list <- list()

# Loop through every gene in our significant list
for (g in gene_names) {
  
  # --- UNIVARIATE MODEL ---
  # Create data frame for this gene (no covariate filtering needed)
  dfg_uni <- make_expr_df(g, expr_vst, clinical_subtype)
  if (is.null(dfg_uni)) next # Skip if gene had < 10 samples
  
  # Run the model: Survival ~ Gene Expression
  fit_uni <- tryCatch(
    coxph(Surv(overall_survival, deceased) ~ counts, data = dfg_uni),
    error = function(e) { return(NULL) } # Return NULL if model fails
  )
  
  # If the model was successful...
  if (!is.null(fit_uni)) {
    # Tidy the output (get HR, p-value) and store it
    uni_list[[g]] <- tidy(fit_uni, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(gene = g, n_uni = nrow(dfg_uni))
    
    # Plot the univariate forest plot
    plot_univariate_with_ggforest(fit_uni, dfg_uni, g, out_dir)  
  }
  
  # --- MULTIVARIATE MODEL ---
  # Create data frame, this time filtering for complete cases of Age and Stage
  dfg_multi <- make_expr_df(g, expr_vst, clinical_subtype_multi, covariates = clinical_covariates)
  
  # Check sample size again (NAs in covariates might reduce sample size)
  if (is.null(dfg_multi) || nrow(dfg_multi) < 20) { # Use a larger threshold
    message("Skipping multi-Cox for ", g, " (data incomplete or N<20 after covariate filter)")
    next # Skip this gene
  }
  
  # Run the model: Survival ~ Gene Expression + Age + Stage
  fit_multi <- tryCatch(
    coxph(as.formula(multi_formula_base), data = dfg_multi),
    error = function(e) {
      message("Error with multivariable fit for ", g, ": ", e$message)
      return(NULL) # Return NULL if model fails
    }
  )
  
  # If the multivariate model was successful...
  if (!is.null(fit_multi)) {
    # Tidy the output...
    gene_result <- tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "counts") %>%  # ...but ONLY keep the row for our gene ("counts")
      mutate(gene = g, n_multi = nrow(dfg_multi))
    
    # Store the gene-specific result
    multi_list[[g]] <- gene_result
    
    # Save the FULL multivariable model summary (including Age and Stage)
    # This is useful for checking the effect of covariates.
    write.csv(tidy(fit_multi, exponentiate = TRUE, conf.int = TRUE),
              file = file.path(out_dir, paste0("Multivariate_Model_Summary_Full_", g, ".csv")),  
              row.names = FALSE)
  }
  
  # --- Plot Multivariate Forest Plot ---
  multivariate_plot_fname <- file.path(out_dir, paste0("Multivariate_Cox_Regression_Plot_", g, ".png"))
  
  tryCatch({
    # Create the forest plot for the full model
    multi_plot <- ggforest(
      fit_multi,  
      data = dfg_multi, # Use the model's data frame
      main = paste("Multivariate Cox Regression for", g),
      # Adjust column positions to prevent text overlap
      cpositions = c(0.03, 0.20, 0.38), 
      fontsize = 1.0 
    )
    
    # Save the plot
    ggsave(
      filename = multivariate_plot_fname,  
      plot = multi_plot,  
      width = 12, # Wider plot to accommodate all variables
      height = 6, 
      dpi = 300
    )
    message("Saved Multivariate Cox Regression Plot for ", g)
  }, error = function(e) {
    message("Failed to plot multivariate forest plot for ", g, ": ", e$message)
  })
} # End of gene loop

# ----------------------------
# 6. Final Save & Merge Results
# ----------------------------

message("Combining and saving final results...")

# Combine all univariate results into one data frame
if (length(uni_list) > 0) {
  uni_df_final <- bind_rows(uni_list) %>%
    # Calculate FDR-adjusted p-values
    mutate(p.adjusted_uni = p.adjust(p.value, method = "fdr")) %>%
    # Select and rename columns for clarity
    select(gene, n_uni, HR_uni = estimate, p.value_uni = p.value, p.adjusted_uni)
} else {
  stop("No valid univariate results to save.")
}

# Combine all multivariate results into one data frame
if (length(multi_list) > 0) {
  multi_df_final <- bind_rows(multi_list) %>%
    # Calculate FDR-adjusted p-values
    mutate(p.adjusted_multi = p.adjust(p.value, method = "fdr")) %>%
    # Select and rename columns
    select(gene, n_multi, HR_multi = estimate, p.value_multi = p.value, p.adjusted_multi)
} else {
  message("No valid multivariable results to save.")
  multi_df_final <- NULL # Set to NULL if it's empty
}

# Merge the univariate and multivariate results for a final comparison table
if (!is.null(multi_df_final)) {
  final_results <- full_join(uni_df_final, multi_df_final, by = "gene")
} else {
  final_results <- uni_df_final # If multi-cox failed, just save uni-cox
}

# Save the final summary CSV
write.csv(final_results, file = file.path(out_dir, paste0(subtype_name, "_final_cox_comparison_age_stage.csv")), row.names = FALSE)
message("Saved final comparison results.")

message("Full Cox analysis complete.")