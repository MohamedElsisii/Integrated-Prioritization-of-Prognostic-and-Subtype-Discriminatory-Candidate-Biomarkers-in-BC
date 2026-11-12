# ==================================
# Script: Code S8 - Functional Rewiring of lncRNAs.R
#
# Purpose:
# This script performs a functional rewiring analysis to understand how
# long non-coding RNAs (lncRNAs) change their regulatory relationships
# with biological pathways in cancer. It calculates the correlation
# between lncRNA expression and pathway activity (GSVA scores) in both
# normal and tumor tissues, then tests for significant changes in these
# correlations.
#
# Key Outputs:
# 1. A CSV summary of gained/lost pathways for each lncRNA.
# 2. A CSV summary of which lncRNAs are gained/lost for each pathway.
# 3. Plots visualizing the overall "Functional Recoding Index" (FRI) and
#    top rewired pathways.
# ==================================

# --- 1. SETUP ---
# Load necessary libraries, suppressing startup messages for cleaner output
suppressPackageStartupMessages({
  library(GSVA)      # For Gene Set Variation Analysis (pathway activity scores)
  library(msigdbr)   # Contains MSigDB gene sets (e.g., Hallmark pathways)
  library(dplyr)     # For data manipulation
  library(purrr)     # For functional programming (map functions)
  library(ggplot2)   # For visualization
  library(pheatmap)  # (Optional) For heatmap visualization
  library(tidyr)     # For data tidying (pivot_wider, etc.)
  library(tibble)    # For working with tibbles (modern data frames)
})

# Load expression data and sample IDs
# These files should contain your master expression matrices and specific sample lists.
brca_expr <- readRDS("BRCA.rda")
normal_expr <- readRDS("Normal.rds")
# Standardize sample IDs to match the expression matrix column names
basal_samples <- gsub("\\.", "-", readRDS("Basal_Samples.rda"))

# Get Hallmark gene sets for humans from MSigDB
# We use Hallmark pathways as they represent well-defined biological states.
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# --- 2. CORE ANALYSIS FUNCTION ---

#' @title Calculate Functional Recoding for One lncRNA
#' @description Computes the correlation between a single lncRNA and ALL Hallmark
#' pathways in both normal and tumor samples, then tests for a significant
#' difference in these correlations.
#' @param lnc Gene symbol of the lncRNA
#' @param normal_expr Expression matrix for normal samples
#' @param basal_expr Expression matrix for tumor (e.g., Basal-like) samples
#' @param gsva_normal Pre-calculated GSVA scores for normal samples
#' @param gsva_basal Pre-calculated GSVA scores for tumor samples
#' @return A tibble with correlation stats and p-values for all pathways for this lncRNA.
calculate_recoding <- function(lnc, normal_expr, basal_expr, gsva_normal, gsva_basal) {
  # Check if lncRNA exists in both expression matrices
  if (!lnc %in% rownames(normal_expr) || !lnc %in% rownames(basal_expr)) {
    warning(paste("Skipping", lnc, "- not found in expression data."))
    return(NULL)
  }
  
  # Helper function: Fisher's z-test for difference between two independent correlations
  fisher_p <- function(r1, r2, n1, n2) {
    z1 <- 0.5 * log((1 + r1) / (1 - r1))
    z2 <- 0.5 * log((1 + r2) / (1 - r2))
    z_diff <- (z1 - z2) / sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
    2 * pnorm(abs(z_diff), lower.tail = FALSE) # Two-tailed p-value
  }
  
  # Calculate Spearman correlations between lncRNA expression and ALL pathway scores
  # t() is used because cor() expects rows as samples when comparing a vector to a matrix.
  cor_norm <- cor(t(gsva_normal), as.numeric(normal_expr[lnc, ]), method = "spearman")
  cor_basal <- cor(t(gsva_basal), as.numeric(basal_expr[lnc, ]), method = "spearman")
  
  # Combine results into a tidy data frame
  tibble(
    lncRNA = lnc,
    Pathway = rownames(cor_norm),
    r_normal = as.vector(cor_norm),
    r_tumor = as.vector(cor_basal)
  ) %>%
    mutate(
      # delta_r: Positive = Gained correlation in tumor, Negative = Lost correlation
      delta_r = r_tumor - r_normal,
      # Calculate p-value for each pathway's correlation change
      p_value = map2_dbl(r_normal, r_tumor, ~fisher_p(.x, .y, ncol(normal_expr), ncol(basal_expr))),
      # Adjust p-values for multiple testing (Benjamini-Hochberg)
      adj_p = p.adjust(p_value, method = "BH")
    )
}

# --- 3. EXECUTE MAIN ANALYSIS ---

# Pre-calculate GSVA scores for all samples.
# This is done ONCE outside the loop for efficiency.
message("Calculating GSVA scores for Normal samples...")
gsva_param <- gsvaParam(exprData = as.matrix(normal_expr), geneSets = hallmarks)
gsva_normal <- gsva(gsva_param, verbose = FALSE)

message("Calculating GSVA scores for Tumor samples...")
gsva_param_basal <- gsvaParam(exprData = as.matrix(brca_expr[, basal_samples]), geneSets = hallmarks)
gsva_basal <- gsva(gsva_param_basal, verbose = FALSE)

# Define the list of lncRNAs to analyze
lncRNAs <- c(
  "LINC02343","PCAT18","LINC00504","KLHDC7B-DT","MAGEB17-AS1","KIRREL3-AS1",
  "LINC00993","SLIT3-AS2","PPP1R3B-DT","GASK1B-AS1","LINC01411","CT62",
  "LINC02568","LINC01843","KIRREL3-AS4","LINC01016","NRIP3-DT","NTM-AS3",
  "LINC02732","MAPT-IT1","LINC01087","VIPR1-AS1","LINC01956","LINC00511",
  "DEPDC1-AS1","LINC01615"
)

message("Running rewiring analysis for ", length(lncRNAs), " lncRNAs...")
# Use map_dfr to loop through all lncRNAs and combine results into one big data frame
full_results <- map_dfr(lncRNAs, ~calculate_recoding(
  .x, normal_expr, brca_expr[, basal_samples], gsva_normal, gsva_basal
))

# --- 4a. CREATE LNC-CENTRIC SUMMARY FILE ---
# This summary shows which PATHWAYS are changed for each LNC.

# Start with a complete list of lncRNAs to ensure all are included in the final summary
all_lncRNAs_list <- distinct(full_results, lncRNA)

# Summarize GAINED pathways (significant increase in correlation)
gained_summary <- full_results %>%
  filter(adj_p < 0.05, delta_r > 0) %>%
  group_by(lncRNA) %>%
  summarise(
    Gained_Pathways = paste(Pathway, collapse = ", "),
    Num_Gained = n(),
    Mean_Gained_Delta_r = mean(delta_r, na.rm = TRUE),
    .groups = "drop"
  )

# Summarize LOST pathways (significant decrease in correlation)
lost_summary <- full_results %>%
  filter(adj_p < 0.05, delta_r < 0) %>%
  group_by(lncRNA) %>%
  summarise(
    Lost_Pathways = paste(Pathway, collapse = ", "),
    Num_Lost = n(),
    Mean_Lost_Delta_r = mean(delta_r, na.rm = TRUE),
    .groups = "drop"
  )

# Merge everything into one final lncRNA summary table
lnc_pathway_summary <- all_lncRNAs_list %>%
  left_join(gained_summary, by = "lncRNA") %>%
  left_join(lost_summary, by = "lncRNA") %>%
  # Replace NAs with 0s for numeric columns where appropriate
  mutate(across(starts_with("Num_"), ~ifelse(is.na(.), 0, .)),
         across(starts_with("Mean_"), ~ifelse(is.na(.), 0, .)))

write.csv(lnc_pathway_summary, "lncRNA_pathway_summary_with_stats.csv", row.names = FALSE)
cat("Successfully created 'lncRNA_pathway_summary_with_stats.csv'\n")

# --- 4b. CREATE PATHWAY-CENTRIC SUMMARY FILE ---
# This summary shows which LNCs are changed for each PATHWAY.

all_pathways_list <- distinct(full_results, Pathway)

# Summarize which lncRNAs GAINED correlation with each pathway
gained_lnc_summary <- full_results %>%
  filter(adj_p < 0.05, delta_r > 0) %>%
  group_by(Pathway) %>%
  summarise(
    Gained_lncRNAs = paste(lncRNA, collapse = ", "),
    Num_Gained_lncRNAs = n(), 
    .groups = "drop"
  )

# Summarize which lncRNAs LOST correlation with each pathway
lost_lnc_summary <- full_results %>%
  filter(adj_p < 0.05, delta_r < 0) %>%
  group_by(Pathway) %>%
  summarise(
    Lost_lncRNAs = paste(lncRNA, collapse = ", "),
    Num_Lost_lncRNAs = n(), 
    .groups = "drop"
  )

# Merge into one final pathway summary table
pathway_lncRNA_summary <- all_pathways_list %>%
  left_join(gained_lnc_summary, by = "Pathway") %>%
  left_join(lost_lnc_summary, by = "Pathway") %>%
  mutate(
    Num_Gained_lncRNAs = ifelse(is.na(Num_Gained_lncRNAs), 0, Num_Gained_lncRNAs),
    Num_Lost_lncRNAs = ifelse(is.na(Num_Lost_lncRNAs), 0, Num_Lost_lncRNAs)
  )

write.csv(pathway_lncRNA_summary, "pathway_lncRNA_summary.csv", row.names = FALSE)
cat("Successfully created 'pathway_lncRNA_summary.csv' with lncRNA counts.\n\n")

# --- 5. POST-PROCESSING & VISUALIZATION FUNCTIONS ---

#' @title Analyze and Visualize Recoding Results
#' @description Generates summary plots for Functional Recoding Index (FRI)
#' and top gained/lost pathways.
#' @param results The full results data frame from calculate_recoding
#' @return A list containing ggplot objects
analyze_and_visualize_recoding <- function(results) {
  
  # --- 5a. Functional Recoding Index (FRI) ---
  # FRI is the mean absolute change in correlation across ALL pathways for a gene.
  # It represents how much a gene's functional "wiring" has changed overall.
  summary_fri <- results %>%
    group_by(lncRNA) %>%
    summarise(
      FRI = mean(abs(delta_r), na.rm = TRUE),
      # Check if it had ANY significant rewiring event
      min_adjP = min(adj_p, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_fri <- ggplot(summary_fri, aes(x = reorder(lncRNA, FRI), y = FRI, fill = min_adjP < 0.05)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#8b5cf6", "FALSE" = "grey80"), name = "Significant\n(Any Pathway)") +
    labs(title = "Functional Recoding Index (FRI) of lncRNAs", y = "FRI (Mean |Δr|)", x = NULL) +
    theme_minimal()
  
  # --- 5b. Analyze commonly gained/lost pathways ---
  # Identify global trends: which pathways are most frequently rewired across all lncRNAs?
  sig_results <- results %>% filter(adj_p < 0.05)
  
  pathway_trends <- sig_results %>%
    mutate(direction = ifelse(delta_r > 0, "Gain", "Loss")) %>%
    group_by(Pathway, direction) %>%
    summarise(n_lnc = n(), .groups = "drop") %>%
    pivot_wider(names_from = direction, values_from = n_lnc, values_fill = 0)
  
  # Plot Top 10 Gained Pathways
  plot_top_gained <- pathway_trends %>%
    arrange(desc(Gain)) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Pathway, Gain), y = Gain)) +
    geom_col(fill = "#1f77b4") +
    coord_flip() +
    labs(title = "Top 10 Pathways Gained Across lncRNAs", y = "Number of lncRNAs", x = NULL) +
    theme_minimal()
  
  # Plot Top 10 Lost Pathways
  plot_top_lost <- pathway_trends %>%
    arrange(desc(Loss)) %>%
    head(10) %>%
    ggplot(aes(x = reorder(Pathway, Loss), y = Loss)) +
    geom_col(fill = "#d62728") +
    coord_flip() +
    labs(title = "Top 10 Pathways Lost Across lncRNAs", y = "Number of lncRNAs", x = NULL) +
    theme_minimal()
  
  # --- Save intermediate summary tables ---
  write.csv(summary_fri, "summary_fri.csv", row.names = FALSE)
  write.csv(pathway_trends, "pathway_trends.csv", row.names = FALSE)
  
  # Return plots as a list
  list(
    FRI_Plot = plot_fri,
    Top_Gained_Plot = plot_top_gained,
    Top_Lost_Plot = plot_top_lost
  )
}

# --- 6. RUN POST-PROCESSING & VIEW/SAVE RESULTS ---
message("Generating summary plots...")
analysis_results <- analyze_and_visualize_recoding(full_results)

# You can view these plots interactively in RStudio by running these lines:
# analysis_results$FRI_Plot
# analysis_results$Top_Gained_Plot
# analysis_results$Top_Lost_Plot

# Save the main FRI plot to a file
ggsave("FRI_plot.png", analysis_results$FRI_Plot, width = 8, height = 6)
message("Analysis complete. FRI plot saved to FRI_plot.png")