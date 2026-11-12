# =============================================================================
# Code S5 - Molecular subtype-discriminatory performance (one-versus-rest ROC)
#
# For each candidate gene of each subtype, the subtype of interest is the positive
# class and the other three subtypes are the negative class. Discrimination is
# quantified as the mean AUC over repeated stratified 10-fold cross-validation
# (3 repeats). As a sensitivity analysis, the AUC is also computed on the same
# folds with the direction of the comparison learned from the training samples
# of each split (AUC_CV_fixeddir).
# Selection criteria: BH-adjusted Wilcoxon p < 0.05, mean CV AUC >= 0.8, CV AUC SD < 0.05.
#
# Input (working directory)
#   BRCA.rda                   normalized expression matrix (genes x samples)
#   <Subtype>_Candidates.rda   candidate genes for each subtype
#   <Subtype>_Samples.rda      sample IDs for each subtype
# Output
#   ROC_<Subtype>_summary.csv, ROC_summary_by_subtype_with_fixeddir.csv (Tables S19-S22)
#   roc_cv10_plots/            ROC curve of each gene
#   File_S4.pdf                ROC curves of the 181 Basal-like candidates (File S4)
# =============================================================================

library(pROC)
library(dplyr)
library(readr)
library(ggplot2)
library(png)

# ---- Parameters ----
plot_dir <- "roc_cv10_plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
thresholds <- list(
  LumA  = list(min_pos = 15, min_neg = 15),
  LumB  = list(min_pos = 12, min_neg = 12),
  Basal = list(min_pos = 10, min_neg = 10),
  Her2  = list(min_pos = 7,  min_neg = 7)
)

k_folds <- 10
cv_repeats <- 3
set.seed(42)  
log_transform <- FALSE
# ---- Helpers ----
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

load_rds_get <- function(path){
  if (!file.exists(path)) stop("File not found: ", path)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  env <- new.env()
  nm <- tryCatch(load(path, envir = env), error = function(e) character(0))
  if(length(nm) == 0) stop("No readable object found in ", path)
  if(length(nm) > 1) warning("Multiple objects in ", path, " — returning first: ", nm[1])
  return(get(nm[1], envir = env))
}

norm_samples <- function(x) gsub("\\.", "-", x)

# ROC plot: x-axis = 1 - specificity, dashed line = chance
plot_roc_curve <- function(labels, scores, out_file, title = NULL, subtitle = NULL){
  roc_obj <- tryCatch(pROC::roc(labels, scores, quiet = TRUE), error = function(e) NULL)
  if(is.null(roc_obj)) return(FALSE)
  p <- ggroc(roc_obj, legacy.axes = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    coord_equal() +
    labs(title = title, subtitle = subtitle,
         x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
    theme_classic() +
    theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9))
  ggsave(out_file, plot = p, width = 6, height = 6, dpi = 300)
  return(TRUE)
}

# ---- Load BRCA ----
BRCA_obj <- load_rds_get("BRCA.rda")
if (is.list(BRCA_obj) && length(BRCA_obj)==1) BRCA <- BRCA_obj[[1]] else BRCA <- BRCA_obj
if (!is.matrix(BRCA)) BRCA <- as.matrix(BRCA)
mode(BRCA) <- "numeric"

# ---- Load candidates ----
cand_files <- list.files(pattern = "_Candidates\\.rda$", ignore.case = TRUE)
candidate_gene_list <- list()
for (f in cand_files) {
  subname <- sub("_Candidates\\.rda$", "", f, ignore.case = TRUE)
  obj <- load_rds_get(f)
  candidate_gene_list[[subname]] <- unique(as.character(obj))
}

# ---- Load samples ----
sample_files <- list.files(pattern = "_Samples\\.rda$", ignore.case = TRUE)
meta_list <- list()
for (f in sample_files) {
  subname <- sub("_Samples\\.rda$", "", f, ignore.case = TRUE)
  samp <- load_rds_get(f)
  samp <- norm_samples(as.character(samp))
  meta_list[[subname]] <- data.frame(SampleID = samp, Subtype = subname, stringsAsFactors = FALSE)
}
metadata <- bind_rows(meta_list)
metadata <- metadata %>% filter(SampleID %in% colnames(BRCA))
metadata$Subtype <- factor(metadata$Subtype, levels = unique(metadata$Subtype))
BRCA <- BRCA[, unique(metadata$SampleID), drop = FALSE]

# ---- Cross-validation function (stratified, repeated) ----
cv_auc <- function(scores, labels, k = 10, repeats = 3){
  pos_idx <- which(labels==1)
  neg_idx <- which(labels==0)
  if (length(pos_idx) < 2 || length(neg_idx) < 2) return(list(mean = NA_real_, sd = NA_real_, n_folds = 0, aucs = numeric(0)))
  all_aucs <- c()
  all_aucs_fixed <- c()
  for(rep_i in seq_len(repeats)){
    folds_pos <- split(sample(pos_idx), rep(1:k, length.out=length(pos_idx)))
    folds_neg <- split(sample(neg_idx), rep(1:k, length.out=length(neg_idx)))
    for(i in seq_len(k)){
      test_idx <- c(folds_pos[[i]], folds_neg[[i]])
      if(length(test_idx) < 4) next
      if(length(unique(labels[test_idx])) < 2) next
      roc_obj <- tryCatch(pROC::roc(labels[test_idx], scores[test_idx], quiet = TRUE), error = function(e) NULL)
      if(!is.null(roc_obj)) all_aucs <- c(all_aucs, as.numeric(pROC::auc(roc_obj)))
      # direction learned on the training samples of the split, applied to the test fold
      train_pos <- setdiff(pos_idx, test_idx); train_neg <- setdiff(neg_idx, test_idx)
      dir <- if (median(scores[train_pos]) >= median(scores[train_neg])) "<" else ">"
      roc_fix <- tryCatch(pROC::roc(labels[test_idx], scores[test_idx], levels = c(0, 1),
                                    direction = dir, quiet = TRUE), error = function(e) NULL)
      if(!is.null(roc_fix)) all_aucs_fixed <- c(all_aucs_fixed, as.numeric(pROC::auc(roc_fix)))
    }
  }
  if(length(all_aucs) == 0) return(list(mean = NA_real_, sd = NA_real_, n_folds = 0, aucs = numeric(0),
                                        mean_fixed = NA_real_, sd_fixed = NA_real_))
  return(list(mean = mean(all_aucs), sd = ifelse(length(all_aucs)>1, sd(all_aucs), NA_real_), n_folds = length(all_aucs), aucs = all_aucs,
              mean_fixed = mean(all_aucs_fixed),
              sd_fixed = ifelse(length(all_aucs_fixed)>1, sd(all_aucs_fixed), NA_real_)))
}

# ---- Loop subtypes -> genes ----
results <- list()
for (sub in names(candidate_gene_list)) {
  cat("\nProcessing subtype:", sub, "with", length(candidate_gene_list[[sub]]), "candidates\n")
  genes <- candidate_gene_list[[sub]]
  subtype_results <- list()
  pos_samps <- metadata$SampleID[metadata$Subtype == sub]
  neg_samps <- setdiff(colnames(BRCA), pos_samps)
  if(length(pos_samps) == 0){
    warning("No positive samples found for subtype: ", sub, " — skipping.")
    next
  }
  for (g in genes){
    if (!(g %in% rownames(BRCA))){
      subtype_results[[length(subtype_results)+1]] <- data.frame(Subtype=sub, Gene=g, InMatrix=FALSE, stringsAsFactors=FALSE)
      next
    }
    expr_pos <- as.numeric(BRCA[g, pos_samps])
    expr_neg <- as.numeric(BRCA[g, neg_samps])
    if(log_transform){
      expr_pos <- log2(expr_pos + 1)
      expr_neg <- log2(expr_neg + 1)
    }
    npos <- sum(!is.na(expr_pos)); nneg <- sum(!is.na(expr_neg))
    mean_pos <- mean(expr_pos, na.rm=TRUE)
    mean_neg <- mean(expr_neg, na.rm=TRUE)
    mean_diff <- mean_pos - mean_neg
    w_p <- tryCatch(wilcox.test(expr_pos, expr_neg, alternative="two.sided")$p.value, error=function(e) NA_real_)
    labels <- c(rep(1,length(expr_pos)), rep(0,length(expr_neg)))
    scores <- c(expr_pos, expr_neg)
    ok <- !is.na(scores) & !is.na(labels)
    auc_cv_mean <- NA_real_; auc_cv_sd <- NA_real_; auc_nfolds <- 0
    auc_fix_mean <- NA_real_; auc_fix_sd <- NA_real_
    thr <- thresholds[[sub]]
    if(sum(ok & labels==1) >= thr$min_pos && sum(ok & labels==0) >= thr$min_neg){
      cv <- cv_auc(scores[ok], labels[ok], k=k_folds, repeats=cv_repeats)
      auc_cv_mean <- cv$mean; auc_cv_sd <- cv$sd; auc_nfolds <- cv$n_folds
      auc_fix_mean <- cv$mean_fixed; auc_fix_sd <- cv$sd_fixed
      fname <- file.path(plot_dir, paste0("ROC_", safe_name(sub), "_", safe_name(g), ".png"))
      subtitle <- paste0("AUC=", ifelse(is.na(auc_cv_mean), "NA", formatC(auc_cv_mean, digits=3)),
                         ", sd=", ifelse(is.na(auc_cv_sd), "NA", formatC(auc_cv_sd, digits=3)),
                         ", folds=", auc_nfolds)
      if(auc_nfolds > 0){
        tryCatch(plot_roc_curve(labels[ok], scores[ok], fname,
                                title = paste0(g, " | ", sub),
                                subtitle = subtitle), error=function(e) NULL)
      }
    } else {
      message("Skipping CV for ", g, " in ", sub, " due to small sample counts: pos=", sum(ok & labels==1), " neg=", sum(ok & labels==0))
    }
    subtype_results[[length(subtype_results)+1]] <- data.frame(
      Subtype = sub, Gene = g, InMatrix = TRUE, n_pos = npos, n_neg = nneg,
      Mean_Pos = mean_pos, Mean_Neg = mean_neg, Mean_Diff = mean_diff,
      Wilcox_p = w_p, AUC_CV = auc_cv_mean, AUC_CV_SD = auc_cv_sd, AUC_nfolds = auc_nfolds,
      AUC_CV_fixeddir = auc_fix_mean, AUC_CV_SD_fixeddir = auc_fix_sd,
      stringsAsFactors = FALSE
    )
  }
  res_df_sub <- bind_rows(subtype_results)
  # adjust p-values (FDR) per subtype
  res_df_sub$Wilcox_p_num <- as.numeric(res_df_sub$Wilcox_p)
  res_df_sub$Wilcox_p_num[is.na(res_df_sub$Wilcox_p_num)] <- 1
  res_df_sub$Wilcox_p_adj <- p.adjust(res_df_sub$Wilcox_p_num, method = "BH")
  res_df_sub <- res_df_sub %>% arrange(desc(AUC_CV))
  file_name <- paste0("ROC_", safe_name(sub), "_summary.csv")
  readr::write_csv(res_df_sub, file_name)
  message("Wrote ", file_name)
  results[[sub]] <- res_df_sub
}

# Combined summary
if(length(results) > 0){
  res_df <- bind_rows(results)
  res_df <- res_df %>% arrange(Subtype, desc(AUC_CV))
  readr::write_csv(res_df, "ROC_summary_by_subtype.csv")
  message("Wrote ROC_summary_by_subtype.csv")
  
  # does fixing the direction change the selected candidates?
  sel <- function(auc, sd) res_df$Wilcox_p_adj < 0.05 & auc >= 0.8 & sd < 0.05
  res_df$selected_original <- sel(res_df$AUC_CV, res_df$AUC_CV_SD)
  res_df$selected_fixeddir <- sel(res_df$AUC_CV_fixeddir, res_df$AUC_CV_SD_fixeddir)
  print(res_df %>% group_by(Subtype) %>%
          summarise(selected_original = sum(selected_original, na.rm = TRUE),
                    selected_fixeddir = sum(selected_fixeddir, na.rm = TRUE),
                    lost = sum(selected_original & !selected_fixeddir, na.rm = TRUE),
                    max_AUC_change = if (any(selected_original, na.rm = TRUE))
                      max(abs(AUC_CV - AUC_CV_fixeddir)[selected_original], na.rm = TRUE) else NA_real_,
                    .groups = "drop"))
  print(res_df %>% filter(Gene %in% c("FOXC1", "STAC", "HMGN2P15")) %>%
          select(Subtype, Gene, AUC_CV, AUC_CV_SD, AUC_CV_fixeddir, AUC_CV_SD_fixeddir))
  readr::write_csv(res_df, "ROC_summary_by_subtype_with_fixeddir.csv")
} else {
  message("No results to save.")
}

sel_genes <- res_df %>%
  filter(Subtype == "Basal", Wilcox_p_adj < 0.05, AUC_CV >= 0.8, AUC_CV_SD < 0.05) %>%
  arrange(Gene) %>% pull(Gene)
cat("Basal-like candidates:", length(sel_genes), "\n")

# File S4: ROC curves of the Basal-like candidates
dir.create("File_S4_plots", showWarnings = FALSE)
files <- file.path(plot_dir, paste0("ROC_Basal_", safe_name(sel_genes), ".png"))
cat("Missing plots:", sum(!file.exists(files)), "\n")
file.copy(files, "File_S4_plots", overwrite = TRUE)

# one plot per page, alphabetical
pdf("File_S4.pdf", width = 6, height = 6)
for (f in files) {
  grid::grid.newpage()
  grid::grid.raster(png::readPNG(f))
}
dev.off()
