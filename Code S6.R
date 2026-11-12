# ==================================
# Script: Code S6 - Correlation with PAM50 Genes.R
#
# Purpose:
# This script calculates the correlation between candidate biomarkers for each
# breast cancer subtype and the standard PAM50 signature genes.
# It automatically selects between Pearson (parametric) and Spearman (non-parametric)
# correlation based on the normality of the candidate gene's expression data.
#
# Key Outputs:
# 1. Per-subtype CSV files containing all correlation results with normality tests.
# 2. Per-subtype CSV files containing only significant correlations after FDR correction.
# ==================================

# Step 1: Install and load necessary R packages
# These packages are used for data manipulation, reading/writing files, and accessing PAM50 data.
suppressPackageStartupMessages({
  library(dplyr)     # For data manipulation (filter, arrange, etc.)
  library(pheatmap)  # (Optional) For heatmap visualization if needed later
  library(readr)     # For efficient reading/writing of CSV files
  library(tidyr)     # For data tidying
  library(genefu)    # Contains the PAM50 gene signature data
  library(stats)     # For statistical tests (cor.test, shapiro.test)
})

# Load the standard PAM50 gene list from the 'genefu' package
data(pam50)
pam50_genes = pam50$centroids.map$probe

# Step 2: Load your master gene expression data
# This file should contain the expression matrix (genes x samples)
BRCA <- readRDS("BRCA.rda")

# Step 3: Define your candidate genes and subtypes as a NAMED LIST
# You can manually define small lists or load larger lists from files.
luma_candidates <- c("F2RL2")
lumb_candidates <- c("GLYATL2","ST8SIA2","PID1","LRIG3","FOXC1","STAC")
basal_candidates <- readRDS("Basal_Candidates.rda")
her2_candidates <- c("HMGN2P15")

# Combine all candidate lists into one master list, named by subtype
candidate_gene_list <- list(
  "LumA" = luma_candidates,
  "LumB" = lumb_candidates,
  "Basal" = basal_candidates,
  "Her2" = her2_candidates
)

# Step 4: Loop through each subtype and perform the correlation analysis
# Define thresholds for what is considered "significant"
alpha_adj <- 0.05      # Significance threshold for FDR-adjusted p-value
min_abs_r <- 0.3       # Minimum absolute correlation coefficient (|r|) to be considered relevant

# ---- Main Analysis Loop ----
# Iterate through each subtype defined in candidate_gene_list
for (subtype in names(candidate_gene_list)) {
  message("\nProcessing subtype: ", subtype)
  
  # 4a. Load sample IDs for the current subtype
  subtype_samples_file <- paste0(subtype, "_Samples.rda")
  if (!file.exists(subtype_samples_file)) {
    message("Warning: Sample file for ", subtype, " not found. Skipping.")
    next
  }
  
  subtype_samples <- readRDS(subtype_samples_file)
  subtype_samples <- gsub("\\.", "-", subtype_samples) # Standardize sample ID format
  
  current_candidates <- candidate_gene_list[[subtype]]
  
  # 4b. Identify genes present in your expression matrix
  # It's important to only use genes that actually have data.
  common_candidate_genes <- intersect(current_candidates, rownames(BRCA))
  common_pam50_genes    <- intersect(pam50_genes, rownames(BRCA))
  
  if (length(common_candidate_genes) == 0 || length(common_pam50_genes) == 0) {
    message("Warning: One or more genes for ", subtype, " not found in BRCA. Skipping.")
    next
  }
  
  # 4c. Identify samples present in your expression matrix
  samples_in_data <- intersect(colnames(BRCA), subtype_samples)
  # Correlation requires at least 3 paired observations
  if (length(samples_in_data) < 3) {
    message("Warning: Not enough samples (need >=3) for ", subtype, " to compute correlation. Skipping.")
    next
  }
  
  # 4d. Subset the expression matrix for just this subtype's samples and relevant genes
  # Transpose (t()) so that rows are samples and columns are genes, which is easier for cor.test
  sub_mat <- BRCA[c(common_candidate_genes, common_pam50_genes), samples_in_data, drop = FALSE]
  filtered_expr <- t(sub_mat)
  
  # ------------------------
  # 4e. Precompute Normality (Shapiro-Wilk test)
  # We test if the expression of each candidate gene follows a normal distribution.
  # This determines which correlation method to use later.
  # ------------------------
  normal_alpha <- 0.05 # Threshold for Shapiro-Wilk test (p < 0.05 means NOT normal)
  
  cand_shapiro_p <- sapply(common_candidate_genes, function(g) {
    vals <- as.numeric(filtered_expr[, g])
    vals <- vals[is.finite(vals)] # Remove any infinite/NA values
    if (length(vals) < 3) return(NA_real_) # Need >= 3 samples for test
    # Perform test, return NA if it fails for any reason
    p <- tryCatch(shapiro.test(vals)$p.value, error = function(e) NA_real_)
    return(p)
  }, USE.NAMES = TRUE)
  
  # (Optional) Compute normality for PAM50 genes just for reporting purposes
  pam50_shapiro_p <- sapply(common_pam50_genes, function(g) {
    vals <- as.numeric(filtered_expr[, g])
    vals <- vals[is.finite(vals)]
    if (length(vals) < 3) return(NA_real_)
    p <- tryCatch(shapiro.test(vals)$p.value, error = function(e) NA_real_)
    return(p)
  }, USE.NAMES = TRUE)
  
  # 4f. Perform Pairwise Correlations
  # We will test every candidate gene against every PAM50 gene.
  results_list <- vector("list", length(common_candidate_genes) * length(common_pam50_genes))
  idx <- 1
  
  # Counters for summary statistics
  cnt_skipped_pairs <- 0
  cnt_pearson <- 0
  cnt_spearman <- 0
  
  for (cand_gene in common_candidate_genes) {
    # Decide correlation method based on candidate gene's normality
    p_sh_cand <- cand_shapiro_p[cand_gene]
    # If p >= 0.05, data is normal -> use Pearson. Otherwise -> use Spearman.
    method_for_candidate <- if (!is.na(p_sh_cand) && p_sh_cand >= normal_alpha) "pearson" else "spearman"
    
    # Update counters
    if (method_for_candidate == "pearson") {
      cnt_pearson <- cnt_pearson + length(common_pam50_genes)
    } else {
      cnt_spearman <- cnt_spearman + length(common_pam50_genes)
    }
    
    for (pam50_gene in common_pam50_genes) {
      
      # Extract expression values for the gene pair
      x <- as.numeric(filtered_expr[, cand_gene])
      y <- as.numeric(filtered_expr[, pam50_gene])
      # Keep only samples where BOTH genes have valid data
      good_idx <- which(is.finite(x) & is.finite(y))
      
      # Skip if too few paired observations
      if (length(good_idx) < 3) {
        cnt_skipped_pairs <- cnt_skipped_pairs + 1
        # Record as NA result
        results_list[[idx]] <- list(
          Candidate_Gene = cand_gene, PAM50_Gene = pam50_gene, Correlation_r = NA_real_,
          P_value = NA_real_, Method = NA_character_, N = length(good_idx),
          Shapiro_p_x = p_sh_cand, Shapiro_p_y = pam50_shapiro_p[pam50_gene]
        )
        idx <- idx + 1
        next
      }
      
      xf <- x[good_idx]
      yf <- y[good_idx]
      method_used <- method_for_candidate
      
      # Perform the correlation test
      cor_result <- tryCatch(
        cor.test(xf, yf, method = method_used, alternative = "two.sided"),
        error = function(e) NULL
      )
      
      # Store the results
      if (is.null(cor_result)) {
        results_list[[idx]] <- list(
          Candidate_Gene = cand_gene, PAM50_Gene = pam50_gene, Correlation_r = NA_real_,
          P_value = NA_real_, Method = method_used, N = length(good_idx),
          Shapiro_p_x = p_sh_cand, Shapiro_p_y = pam50_shapiro_p[pam50_gene]
        )
      } else {
        results_list[[idx]] <- list(
          Candidate_Gene = cand_gene, PAM50_Gene = pam50_gene,
          Correlation_r = as.numeric(cor_result$estimate), # The correlation coefficient (r or rho)
          P_value = cor_result$p.value,                    # The raw p-value
          Method = method_used,                            # Method used (pearson/spearman)
          N = length(good_idx),                            # Number of samples used
          Shapiro_p_x = p_sh_cand, Shapiro_p_y = pam50_shapiro_p[pam50_gene]
        )
      }
      idx <- idx + 1
    }
  }
  
  # 4g. Compile and Process Results
  # Convert list of lists to a data frame
  results_list <- results_list[seq_len(idx-1)]
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  # Ensure columns are correct types
  results_df[] <- lapply(results_df, function(x) if (is.factor(x)) as.character(x) else x)
  results_df$Correlation_r <- as.numeric(as.character(results_df$Correlation_r))
  results_df$P_value <- as.numeric(as.character(results_df$P_value))
  results_df$N <- as.integer(as.character(results_df$N))
  results_df$Shapiro_p_x <- as.numeric(as.character(results_df$Shapiro_p_x))
  results_df$Shapiro_p_y <- as.numeric(as.character(results_df$Shapiro_p_y))
  
  # Calculate FDR-adjusted p-values (Benjamini-Hochberg method)
  results_df$P_adj <- NA_real_
  non_na_p <- !is.na(results_df$P_value)
  if (any(non_na_p)) results_df$P_adj[non_na_p] <- p.adjust(results_df$P_value[non_na_p], method = "BH")
  
  # Handle extremely small p-values for display (avoid true 0)
  zero_idx <- which(!is.na(results_df$P_adj) & results_df$P_adj == 0)
  if (length(zero_idx) > 0) {
    results_df$P_adj[zero_idx] <- .Machine$double.xmin
  }
  
  # Add helpful columns for filtering and interpretation
  results_df$Direction <- ifelse(is.na(results_df$Correlation_r), NA, 
                                 ifelse(results_df$Correlation_r > 0, "positive", "negative"))
  results_df$Abs_r <- abs(results_df$Correlation_r)
  
  # Create a filtered data frame with only significant results
  significant_by_adj <- results_df %>%
    filter(!is.na(P_adj) & P_adj < alpha_adj & Abs_r >= min_abs_r) %>%
    arrange(P_adj)
  
  # Print summary to console
  message("Summary for ", subtype, ": ",
          nrow(significant_by_adj),
          " significant pairs (P_adj <", alpha_adj, " and |r| >=", min_abs_r, ").")
  message("Pairs skipped due to too few observations: ", cnt_skipped_pairs)
  message("Pearson tests run: ", cnt_pearson, " ; Spearman tests run: ", cnt_spearman)
  
  # 4h. Save Results to Files
  # Save full results (all pairs)
  write_csv(results_df, paste0(subtype, "_all_correlations_with_normality.csv"))
  saveRDS(results_df, paste0(subtype, "_all_correlations_with_normality.rds"))
  
  # Save only significant results
  write_csv(significant_by_adj, paste0(subtype, "_significant_correlations_adjP.csv"))
  saveRDS(significant_by_adj, paste0(subtype, "_significant_correlations_adjP.rds"))
  
  # Preview top significant results
  if (nrow(significant_by_adj) > 0) {
    print(head(significant_by_adj, 10))
  } else {
    message("No significant correlations after FDR correction for ", subtype, ".")
  }
  
} # End of subtype loop