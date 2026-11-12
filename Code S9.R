# ==================================
# Script: Code S9 - Functional Rewiring of Protein-Coding Genes.R
#
# Purpose:
# This script performs a functional rewiring analysis for protein-coding genes
# identified as either prognostic or diagnostic biomarkers in specific subtypes.
# It groups genes by Subtype and Marker Type (e.g., "Basal-like_Diagnostic")
# and calculates the change in correlation (rewiring) between each gene and
# all Hallmark pathways from normal to tumor tissue.
#
# Key Outputs:
# 1. Per-group CSV summaries of rewired pathways for each gene.
# 2. Per-group CSV summaries of which genes are rewired for each pathway.
# ==================================

# --- 1. SETUP ---
# Load necessary libraries
suppressPackageStartupMessages({
  library(GSVA)      # For pathway activity scoring
  library(msigdbr)   # For Hallmark gene sets
  library(dplyr)     # Data manipulation
  library(purrr)     # Functional programming (map functions)
  library(ggplot2)   # Visualization
  library(pheatmap)  # Heatmaps
  library(tidyr)     # Data tidying
  library(tibble)    # Modern data frames
})

# Load master expression data
# brca_expr: Tumor samples (genes x samples)
# normal_expr: Normal tissue samples (genes x samples)
brca_expr <- readRDS("BRCA.rda")
normal_expr <- readRDS("Normal.rds")

# Load sample lists for each relevant subtype
# Standardize sample names by replacing dots with dashes to match expression matrix
lumA_samples <- gsub("\\.", "-", readRDS("LumA_Samples.rda"))
lumB_samples <- gsub("\\.", "-", readRDS("LumB_Samples.rda"))
basal_samples <- gsub("\\.", "-", readRDS("Basal_Samples.rda"))

# Get Hallmark gene sets for humans from MSigDB
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# --- 2. DEFINE GENES AND ANALYSIS GROUPS ---
# Manually define which genes are prognostic for specific subtypes based on prior analysis
prognostic_lumB <- c("GLYATL2", "ST8SIA2", "PID1", "LRIG3")
prognostic_basal <- c("TRIM59", "DLL3", "CIP2A")

# Create a master table of all protein-coding genes to be analyzed.
# We assign each gene to its Subtype and determine its 'Marker_Type' (Prognostic vs. Diagnostic).
protein_coding_genes <- tibble(
  Subtype = c(rep("Basal-like", 146), "Luminal A", rep("Luminal B", 6)),
  Gene = c(
    # ... (Long list of Basal-like diagnostic genes) ...
    "POU4F1", "DACH1", "TFF3", "GRPR", "PTPRT", "INPP5J", "ANKRD30B", "SLC7A13", "CYP4F8", "ARHGEF38",
    "RAET1L", "BPIFB2", "ACER2", "THSD4", "A2ML1", "DNAH5", "CFAP70", "ANKRD30A", "SYT9", "TTC6",
    "CFAP46", "OPRK1", "FAM241A", "C9orf152", "B3GNT4", "DDN", "AR", "ZIC1", "HCAR1", "SMIM31",
    "CLSTN2", "ABLIM3", "PPP1R3C", "NCBP2L", "ZIC4", "MSX2", "DNAH7", "FER1L5", "ACRV1", "AFF3",
    "LONRF2", "ABCC8", "CYP4Z1", "SCGB1D2", "SCGB2A2", "AUNIP", "INPP4B", "DHRS2", "GJB6", "PSAT1",
    "KLHDC7A", "HPDL", "SRARP", "RAB6C", "ERBB4", "TNS2", "ACADSB", "INA", "IGDCC3", "IGSF9B", "PIF1",
    "CA12", "FZD9", "TICRR", "LIN7A", "SPMIP5", "APOF", "HASPIN", "TTC36", "RAI2", "NXNL2", "TPSG1",
    "CHEK1", "ATP1A4", "KCNK5", "CENPW", "BCAS1", "DYNLRB2", "SLC7A5", "CYBRD1", "GABBR2", "GAL",
    "ACSM1", "GP2", "PLIN5", "SPDEF", "FAM47E", "SYTL5", "CDCA7", "IQUB", "KRT222", "IGF2BP3",
    "ACOX2", "RGS22", "SLC26A9", "PTGER3", "REEP6", "SYBU", "PHGR1", "KCND3", "GJB3", "SLC7A2",
    "CDC7", "SERPINA11", "LAMP3", "ADAMTS15", "TPRG1", "CHAD", "SCUBE2", "SIDT1", "PLCD4", "ART3",
    "KLHDC1", "CAMKV", "CFAP69", "TMEM232", "UGT2B15", "LYPD6", "IL6ST", "DNAI7", "FMO5", "TBC1D9",
    "NOSTRIN", "RCOR2", "PGLYRP2", "CCDC170", "SLC5A7", "GREB1", "RAPGEF3", "LRP8", "MLPH",
    "IL12RB2", "MYRIP", "PRICKLE2", "MAPT", "UBXN10", "CAPN8", "DPYSL5", "CYP4X1", "SLC40A1",
    "PGBD5", "FYB2", "DNALI1", "TRIM59", "DLL3", "CIP2A",
    # Luminal A gene
    "F2RL2",
    # Luminal B genes
    "FOXC1", "GLYATL2", "STAC", "ST8SIA2", "PID1", "LRIG3"
  )
) %>%
  # Automatically classify genes based on the subtype and predefined prognostic lists
  mutate(
    Marker_Type = case_when(
      Subtype == "Luminal A" ~ "Prognostic", # Only one gene, known to be prognostic
      Subtype == "Luminal B" & Gene %in% prognostic_lumB ~ "Prognostic",
      Subtype == "Luminal B" ~ "Diagnostic", # The rest are diagnostic
      Subtype == "Basal-like" & Gene %in% prognostic_basal ~ "Prognostic",
      Subtype == "Basal-like" ~ "Diagnostic",
      TRUE ~ "Unknown" # Safety catch
    )
  )

# Identify unique analysis groups (e.g., Basal-like_Prognostic, Luminal B_Diagnostic)
# We will loop through these groups to perform targeted analyses.
analysis_groups <- protein_coding_genes %>%
  distinct(Subtype, Marker_Type)

# --- 3. CORE ANALYSIS FUNCTION (GENERALIZED) ---

#' @title Calculate Functional Recoding for One Gene in One Subtype
#' @description Computes correlation between a gene and all pathways in normal
#' vs. a specific tumor subtype, then tests for significant change.
#' @param gene_symbol Gene to analyze
#' @param subtype_name Name of the subtype (for logging warnings)
#' @param normal_expr Normal expression matrix
#' @param subtype_expr Tumor subtype expression matrix
#' @param gsva_normal GSVA scores for normal samples
#' @param gsva_subtype GSVA scores for tumor subtype samples
#' @return A tibble with rewiring stats for all pathways for this gene.
calculate_recoding <- function(gene_symbol, subtype_name, normal_expr, subtype_expr, gsva_normal, gsva_subtype) {
  # Verify gene exists in both datasets
  if (!gene_symbol %in% rownames(normal_expr) || !gene_symbol %in% rownames(subtype_expr)) {
    warning(paste("Skipping", gene_symbol, "for", subtype_name, "- not found in expression data."))
    return(NULL)
  }
  
  # Fisher's z-test for difference in correlations
  fisher_p <- function(r1, r2, n1, n2) {
    z1 <- 0.5 * log((1 + r1) / (1 - r1)); z2 <- 0.5 * log((1 + r2) / (1 - r2))
    z_diff <- (z1 - z2) / sqrt((1/(n1 - 3)) + (1/(n2 - 3)))
    2 * pnorm(abs(z_diff), lower.tail = FALSE)
  }
  
  # Calculate Spearman correlations
  cor_norm <- cor(t(gsva_normal), as.numeric(normal_expr[gene_symbol, ]), method = "spearman")
  cor_subtype <- cor(t(gsva_subtype), as.numeric(subtype_expr[gene_symbol, ]), method = "spearman")
  
  # Combine into a results table
  tibble(
    Gene = gene_symbol,
    Pathway = rownames(cor_norm),
    r_normal = as.vector(cor_norm),
    r_tumor = as.vector(cor_subtype)
  ) %>%
    mutate(
      delta_r = r_tumor - r_normal,
      # Calculate p-value for correlation change
      p_value = map2_dbl(r_normal, r_tumor, ~fisher_p(.x, .y, ncol(normal_expr), ncol(subtype_expr))),
      # FDR adjustment
      adj_p = p.adjust(p_value, method = "BH")
    )
}

# --- 4. EXECUTE MAIN ANALYSIS (LOOPING OVER ANALYSIS GROUPS) ---

# Pre-calculate GSVA for normal samples ONCE, as it's the common baseline.
message("Calculating GSVA scores for Normal samples...")
gsva_param_normal <- gsvaParam(exprData = as.matrix(normal_expr), geneSets = hallmarks)
gsva_normal <- gsva(gsva_param_normal, verbose = FALSE)

# Store sample lists in a named list for easy access in the loop
subtype_samples_list <- list(`Luminal A` = lumA_samples, `Luminal B` = lumB_samples, `Basal-like` = basal_samples)

# MASTER LOOP: Iterate through each unique (Subtype, Marker_Type) group
for (i in 1:nrow(analysis_groups)) {
  current_subtype <- analysis_groups$Subtype[i]
  current_marker_type <- analysis_groups$Marker_Type[i]
  # Create a clean name for output files (e.g., "Basal-like_Diagnostic")
  group_name <- paste(current_subtype, current_marker_type, sep = "_")
  
  cat(sprintf("\n--- Processing Group: %s ---\n", group_name))
  
  # 4a. Prepare Subtype-Specific Data
  # Get the specific samples for the current subtype
  current_samples <- subtype_samples_list[[current_subtype]]
  subtype_expr <- brca_expr[, current_samples]
  
  # Calculate GSVA scores SPECIFICALLY for this tumor subtype
  message("  Calculating GSVA scores for ", current_subtype, "...")
  gsva_param_subtype <- gsvaParam(exprData = as.matrix(subtype_expr), geneSets = hallmarks)
  gsva_subtype <- gsva(gsva_param_subtype, verbose = FALSE)
  
  # 4b. Identify Genes for this Group
  genes_to_process <- protein_coding_genes %>%
    filter(Subtype == current_subtype, Marker_Type == current_marker_type) %>%
    pull(Gene)
  
  # 4c. Run Rewiring Analysis
  message("  Running rewiring analysis for ", length(genes_to_process), " genes...")
  group_results <- map_dfr(genes_to_process, ~calculate_recoding(
    .x,                  # gene_symbol
    current_subtype,     # subtype_name
    normal_expr,
    subtype_expr,
    gsva_normal,
    gsva_subtype
  ))
  
  # --- 5. CREATE AND SAVE GROUP-SPECIFIC SUMMARY FILES ---
  if (nrow(group_results) > 0) {
    # A. GENE-CENTRIC Summary: For each gene, what pathways changed?
    gene_centric_summary <- group_results %>%
      group_by(Gene) %>%
      summarise(
        Num_Gained = sum(adj_p < 0.05 & delta_r > 0, na.rm = TRUE),
        Num_Lost = sum(adj_p < 0.05 & delta_r < 0, na.rm = TRUE),
        # Calculate mean Delta_r only for significant pathways
        Mean_Gained_Delta_r = mean(delta_r[adj_p < 0.05 & delta_r > 0], na.rm = TRUE),
        Mean_Lost_Delta_r = mean(delta_r[adj_p < 0.05 & delta_r < 0], na.rm = TRUE),
        # Create comma-separated lists of changed pathways
        Gained_Pathways = paste(Pathway[adj_p < 0.05 & delta_r > 0], collapse = ", "),
        Lost_Pathways = paste(Pathway[adj_p < 0.05 & delta_r < 0], collapse = ", "),
        .groups = "drop"
      ) %>%
      # Add metadata columns and clean up NaNs (from mean of empty sets)
      mutate(
        Subtype = current_subtype, Marker_Type = current_marker_type, .before = 1,
        Mean_Gained_Delta_r = ifelse(is.nan(Mean_Gained_Delta_r), 0, Mean_Gained_Delta_r),
        Mean_Lost_Delta_r = ifelse(is.nan(Mean_Lost_Delta_r), 0, Mean_Lost_Delta_r)
      )
    
    write.csv(gene_centric_summary, paste0(group_name, "_gene_summary.csv"), row.names = FALSE)
    cat(sprintf("  Successfully created '%s_gene_summary.csv'\n", group_name))
    
    # B. PATHWAY-CENTRIC Summary: For each pathway, which genes changed?
    pathway_centric_summary <- group_results %>%
      filter(adj_p < 0.05) %>% # Only consider significant changes
      group_by(Pathway) %>%
      summarise(
        Num_Gained_Genes = sum(delta_r > 0, na.rm = TRUE),
        Num_Lost_Genes = sum(delta_r < 0, na.rm = TRUE),
        Gained_Genes = paste(Gene[delta_r > 0], collapse = ", "),
        Lost_Genes = paste(Gene[delta_r < 0], collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(Subtype = current_subtype, Marker_Type = current_marker_type, .before = 1)
    
    write.csv(pathway_centric_summary, paste0(group_name, "_pathway_summary.csv"), row.names = FALSE)
    cat(sprintf("  Successfully created '%s_pathway_summary.csv'\n", group_name))
    
  } else {
    cat("  No results generated for this group. Skipping file creation.\n")
  }
} # End of master loop

cat("\n--- All Groups Processed ---\n")