# DEGs_pipeline_cleaned.R
# Cleaned, documented, and optimized differential expression pipeline
# - Input: gene-level count CSV files matching FILE_PATTERN (Ensembl IDs in first column)
# - Outputs: per-contrast TREAT results (CSV + RDS) and optional MA plots
# Notes:
#  - This script expects precomputed sample lists (LumA_Samples.rda etc.) and two RDS
#    expression matrices: BRCA_RNA_DEGS.rds and Normal_RNA_DEGS.rds (genes x samples).
#  - Mapping Ensembl->Symbol attempts biomaRt once (if available) then falls back to
#    org.Hs.eg.db (offline). By default we keep only mapped genes; set keep_entrez_fallback
#    to TRUE to keep ENTrez fallback names like ENTREZ_12345.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(biomaRt)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(data.table)
  library(dplyr)
  library(tibble)
})

# ----------------------------
# User options
# ----------------------------
FILE_PATTERN        <- "_RNA\\.csv$"    # pattern for per-sample/count CSVs (first col = Ensembl ID)
SAVE_SUFFIX         <- ".rds"            # suffix for saved processed matrices
set.seed(123)

out_dir <- "DE_results"
dir.create(out_dir, showWarnings = FALSE)
save_plots <- TRUE
use_voom_quality_weights <- TRUE
min_count_filter <- 10   # currently unused but left for future extension
fdr_cut <- 0.05
lfc_cut <- 2
keep_entrez_fallback_global <- FALSE  # keep ENTREZ fallback names when mapping?

# Pre-saved sample lists (expect these files to exist)
LumA  <- readRDS("LumA_Samples.rda")
LumB  <- readRDS("LumB_Samples.rda")
Her2  <- readRDS("Her2_Samples.rda")
Basal <- readRDS("Basal_Samples.rda")
Normal <- readRDS("Normal_Samples.rda")

# Normalize sample ID format to match matrix colnames (if dots vs dashes)
normalize_ids <- function(x) gsub("\\.", "-", x)
LumA <- normalize_ids(LumA); LumB <- normalize_ids(LumB)
Basal <- normalize_ids(Basal); Her2 <- normalize_ids(Her2); Normal <- normalize_ids(Normal)

# ----------------------------
# Helpers
# ----------------------------
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

# Single-shot mapping using biomaRt (if available) or org.Hs.eg.db fallback.
# Returns a data.table with columns ENSEMBL and SYMBOL (SYMBOL may contain ENTREZ_ fallback)
map_ensembl_batch <- function(ensembl_ids, mart = NULL, keep_entrez_fallback = FALSE) {
  ens_clean <- unique(sub("\\..*$", "", ensembl_ids))
  # Try biomaRt if mart provided
  if (!is.null(mart)) {
    log_msg("Attempting biomaRt mapping for ", length(ens_clean), " IDs...")
    mm <- tryCatch(getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                         filters = "ensembl_gene_id",
                         values = ens_clean,
                         mart = mart), error = function(e) { warning("biomaRt failed: ", e$message); NULL })
    if (!is.null(mm) && nrow(mm) > 0) {
      dt <- as.data.table(mm)
      setnames(dt, c("ensembl_gene_id", "external_gene_name"), c("ENSEMBL", "SYMBOL"))
      dt[, ENSEMBL := sub("\\..*$", "", ENSEMBL)]
      dt <- dt[!(is.na(SYMBOL) | SYMBOL == "")]
      return(dt[, .(ENSEMBL, SYMBOL)])
    }
  }
  # Fallback to org.Hs.eg.db
  log_msg("Using org.Hs.eg.db fallback mapping...")
  sel <- AnnotationDbi::select(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
                               columns = c("SYMBOL", "ENTREZID"))
  dt <- as.data.table(sel)
  dt[, ENSEMBL := sub("\\..*$", "", ENSEMBL)]
  if (keep_entrez_fallback) {
    dt[, SYMBOL := ifelse(is.na(SYMBOL) | SYMBOL == "", paste0("ENTREZ_", ENTREZID), SYMBOL)]
    dt <- dt[, .(ENSEMBL, SYMBOL)]
  } else {
    dt <- dt[!(is.na(SYMBOL) | SYMBOL == ""), .(ENSEMBL, SYMBOL)]
  }
  unique(dt, by = "ENSEMBL")
}

# Fast pre-processing for a single CSV counts file (first column = Ensembl ID)
preprocess_count_matrix_fast <- function(file_path, mart = NULL, use_sparse = TRUE, keep_entrez_fallback = FALSE) {
  log_msg("Reading:", file_path)
  dt <- fread(file_path, header = TRUE, data.table = TRUE, showProgress = FALSE)
  if (ncol(dt) < 2) stop("Expected at least 2 columns: EnsemblID + at least 1 sample column")
  setnames(dt, 1, "EnsemblID")
  dt[, EnsemblID := sub("\\..*$", "", EnsemblID)]
  
  mapping <- map_ensembl_batch(dt$EnsemblID, mart = mart, keep_entrez_fallback = keep_entrez_fallback)
  setkey(mapping, ENSEMBL); setkey(dt, EnsemblID)
  # Join keeping only mapped genes (nomatch = 0)
  joined <- mapping[dt, nomatch = 0]
  # 'joined' will have ENSEMBL, SYMBOL, and the sample columns from dt
  # Drop ENSEMBL column and rename SYMBOL -> Gene
  joined[, ENSEMBL := NULL]
  setnames(joined, "SYMBOL", "Gene")
  
  # Identify sample columns and coerce to integer (by reference)
  sample_cols <- setdiff(colnames(joined), "Gene")
  for (cname in sample_cols) set(joined, j = cname, value = as.integer(joined[[cname]]))
  
  # Aggregate duplicate gene names by summing counts (fast data.table approach)
  agg <- joined[, lapply(.SD, sum, na.rm = TRUE), by = Gene, .SDcols = sample_cols]
  
  # Convert to matrix (genes x samples)
  mat <- as.matrix(agg[, ..sample_cols])
  rownames(mat) <- agg$Gene
  storage.mode(mat) <- "integer"
  
  # Optionally convert to sparse matrix if many zeros
  if (use_sparse) {
    nz <- sum(mat != 0)
    tot <- length(mat)
    sparsity <- 1 - nz / tot
    log_msg(sprintf("Matrix: %d genes x %d samples; sparsity=%.2f", nrow(mat), ncol(mat), sparsity))
    if (sparsity > 0.5) {
      sp <- Matrix::Matrix(mat, sparse = TRUE)
      return(sp)
    }
  }
  mat
}

# ----------------------------
# Batch processing: convert *_RNA.csv -> .rds matrices
# ----------------------------
mart <- NULL
try({ mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") }, silent = TRUE)

files <- list.files(pattern = FILE_PATTERN)
if (length(files) > 0) {
  for (f in files) {
    out <- sub("\\.csv$", SAVE_SUFFIX, f, ignore.case = TRUE)
    if (file.exists(out)) next
    res <- preprocess_count_matrix_fast(f, mart = mart, use_sparse = TRUE, keep_entrez_fallback = keep_entrez_fallback_global)
    saveRDS(res, out)
    rm(res)
  }
} else {
  log_msg("No raw CSV files found for batch processing. Skipping CSV -> RDS step.")
}

# ----------------------------
# Differential expression analysis
# ----------------------------
log_msg("Loading expression matrices...")
brca <- readRDS("BRCA_RNA_DEGS.rds")
normal_expr <- readRDS("Normal_RNA_DEGS.rds")
if (is.data.frame(brca)) brca <- as.matrix(brca)
if (is.data.frame(normal_expr)) normal_expr <- as.matrix(normal_expr)

# Subset sample lists to available columns in matrices
lumA_samps <- intersect(LumA, colnames(brca))
lumB_samps <- intersect(LumB, colnames(brca))
basal_samps <- intersect(Basal, colnames(brca))
her2_samps <- intersect(Her2, colnames(brca))
normal_samps <- intersect(Normal, colnames(normal_expr))

log_msg(sprintf("Samples kept: LumA=%d LumB=%d Basal=%d Her2=%d Normal=%d",
                length(lumA_samps), length(lumB_samps), length(basal_samps), length(her2_samps), length(normal_samps)))

# Subset matrices and keep only common genes between tumor & normal
brca_sub <- brca[, unique(c(lumA_samps, lumB_samps, basal_samps, her2_samps)), drop = FALSE]
normal_sub <- normal_expr[, normal_samps, drop = FALSE]
common_genes <- intersect(rownames(brca_sub), rownames(normal_sub))
log_msg("Common genes between tumor & normal: ", length(common_genes))
brca_sub <- brca_sub[common_genes, , drop = FALSE]
normal_sub <- normal_sub[common_genes, , drop = FALSE]
counts_all <- cbind(brca_sub, normal_sub)

# Metadata: assign each column to a condition
samples_all <- colnames(counts_all)
condition <- ifelse(samples_all %in% normal_samps, "Normal",
                    ifelse(samples_all %in% lumA_samps, "LumA",
                           ifelse(samples_all %in% lumB_samps, "LumB",
                                  ifelse(samples_all %in% basal_samps, "Basal",
                                         ifelse(samples_all %in% her2_samps, "Her2", NA)))))
if (any(is.na(condition))) stop("Some samples could not be assigned to a condition. Check sample lists vs matrix colnames.")
meta <- data.frame(sample = samples_all,
                   condition = factor(condition, levels = c("Normal","LumA","LumB","Basal","Her2")),
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$sample

# DGEList / filtering / normalization
log_msg("Constructing DGEList and filtering low-expression genes...")
dge <- DGEList(counts = counts_all)
storage.mode(dge$counts) <- "integer"
keep <- filterByExpr(dge, group = meta$condition)
log_msg(sprintf("Keeping %d / %d genes after filterByExpr.", sum(keep), nrow(dge)))
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Design + voom
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- sub("^condition", "", colnames(design))
log_msg("Design matrix columns: ", paste(colnames(design), collapse = ", "))
if (use_voom_quality_weights) {
  log_msg("Running voomWithQualityWeights...")
  v <- voomWithQualityWeights(dge, design = design, plot = FALSE)
} else {
  log_msg("Running voom...")
  v <- voom(dge, design = design, plot = FALSE)
}

# Fit linear model, contrasts, and eBayes
fit <- lmFit(v, design)
cont_matrix <- makeContrasts(
  LumA_vs_Normal = LumA - Normal,
  LumB_vs_Normal = LumB - Normal,
  Basal_vs_Normal = Basal - Normal,
  Her2_vs_Normal = Her2 - Normal,
  levels = design
)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# Apply TREAT and save results
lfc_thresholds <- unique(lfc_cut)
for (lfc_thr in lfc_thresholds) {
  log_msg("Applying TREAT with LFC threshold = ", lfc_thr)
  fit_treat <- treat(fit2, lfc = lfc_thr)
  for (cn in colnames(cont_matrix)) {
    res_t <- topTreat(fit_treat, coef = cn, number = Inf, sort.by = "P")
    res_t$Significant <- (res_t$adj.P.Val < fdr_cut)
    out_csv <- file.path(out_dir, paste0("DE_", cn, "_TREAT_lfc", lfc_thr, ".csv"))
    out_rds <- file.path(out_dir, paste0("DE_", cn, "_TREAT_lfc", lfc_thr, ".rds"))
    data.table::fwrite(as.data.table(res_t, keep.rownames = "Gene"), out_csv)
    saveRDS(res_t, out_rds)
    if (save_plots) {
      png(file.path(out_dir, paste0("MAplot_", cn, "_TREAT_lfc", lfc_thr, ".png")), width = 1600, height = 1200, res = 200)
      plotMD(fit_treat, column = which(colnames(cont_matrix) == cn), main = paste0(cn, " (TREAT lfc=", lfc_thr, ")"))
      abline(h = c(-lfc_thr, lfc_thr), col = "blue", lty = 2)
      dev.off()
    }
  }
}

log_msg("DE analysis complete. Results saved to:", normalizePath(out_dir))
