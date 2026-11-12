# =============================================================================
# Code S7 - Pseudogene-parent assignment and decoupling analysis
#
# Part A - parent assignment
#   Query: longest mature transcript of each candidate pseudogene (GENCODE v49).
#   Database: GENCODE v49 protein-coding transcripts.
#   BLASTN with DUST filtering and no limit on the number of hits; hits retained at
#   identity >= 80%, query coverage >= 50%, E-value <= 1e-10. Transcripts are collapsed
#   to genes and the highest-scoring gene is the putative parent (flagged as ambiguous
#   if a second gene scores >= 95% of the best).
# Part B - decoupling statistics
#   Spearman correlation of each pseudogene-parent pair in normal and tumor samples;
#   delta r = r_tumor - r_normal tested with a permutation test (10,000 label shuffles)
#   and Fisher's z-test with the Bonett-Wright variance; bootstrap 95% CI (2,000
#   resamples); BH correction. Repeated after excluding normal samples from patients
#   who also contributed a tumor sample.
#
# Input (working directory)
#   gencode.v49.transcripts.fa.gz   GENCODE v49 transcript sequences
#   BRCA.rda, Normal.rds            tumor and normal expression matrices
#   Her2_Samples.rda, Basal_Samples.rda
# Software: BLAST+ (blastn, makeblastdb); R packages dplyr, ggplot2, biomaRt
# Output
#   Pseudogene_BLAST_hits_filtered.csv, Pseudogene_parent_assignment.csv (Table S28)
#   Pseudogene_Decoupling_Results.csv (Table S29), Pseudogene_Decoupling_Plot.png (Figure 5)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

# =============================== PART A: PARENT ASSIGNMENT ===============================

PSEUDOGENES <- c("HMGN2P15", "TPTE2P2", "TMEM161BP1", "VN1R53P", "TNRC18P1", "AIDAP2", "CPHL1P",
                 "RAB6C-AS1", "GTF2IP7", "RAD17P1", "POTEKP", "AZGP1P1", "PLAC9P1", "CYP4Z2P", "RN7SL314P")
TRANSCRIPTS_FILE <- "gencode.v49.transcripts.fa.gz"  # GENCODE v49 transcripts (.fa or .fa.gz)
BLAST_DIR <- ""   # folder containing blastn/makeblastdb; "" = search PATH and the working directory
BLAST_WORK <- if (.Platform$OS.type == "windows") "C:/pseudoblast" else file.path(tempdir(), "pseudoblast")
# BLAST database paths must not contain spaces
MIN_IDENT <- 80; MIN_QCOV <- 50; MAX_EVALUE <- 1e-10; AMBIG <- 0.95

# ---------------------------------------------------------------- 1. GENCODE transcripts
lines <- readLines(gzfile(TRANSCRIPTS_FILE))
is_h <- startsWith(lines, ">")
grp <- cumsum(is_h)
hdr <- sub("^>", "", lines[is_h])
seqs <- vapply(split(lines[!is_h], grp[!is_h]), paste, character(1), collapse = "")
rm(lines); gc()
if (grepl("|", hdr[1], fixed = TRUE)) {
  # full GENCODE headers: >ENST|ENSG|OTTHUMG|OTTHUMT|transcript_name|gene_name|length|biotype|
  f <- strsplit(hdr, "|", fixed = TRUE)
  tx <- tibble(tx_id = vapply(f, `[`, "", 1), gene_id = vapply(f, `[`, "", 2),
               gene_name = vapply(f, `[`, "", 6), biotype = vapply(f, function(x) x[8], ""),
               len = nchar(seqs), seq = unname(seqs))
} else {
  # headers with transcript IDs only: gene name and biotype from Ensembl 115 (= GENCODE v49)
  cat("Headers contain only transcript IDs - retrieving gene names and biotypes from Ensembl 115...\n")
  mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl", version = 115)
  ann <- biomaRt::getBM(c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name",
                          "transcript_biotype"), mart = mart)
  tx <- tibble(tx_id = sub("\\s.*$", "", hdr), len = nchar(seqs), seq = unname(seqs)) %>%
    mutate(key = sub("\\..*$", "", tx_id)) %>%
    left_join(ann, by = c("key" = "ensembl_transcript_id")) %>%
    transmute(tx_id, gene_id = ensembl_gene_id, gene_name = external_gene_name,
              biotype = transcript_biotype, len, seq)
  cat(sum(is.na(tx$gene_name)), "of", nrow(tx), "transcripts could not be annotated (ignored)\n")
  tx <- tx %>% filter(!is.na(gene_name), gene_name != "")
}

# ---------------------------------------------------------------- 2. queries: longest transcript of each pseudogene
query <- tx %>% filter(gene_name %in% PSEUDOGENES) %>% group_by(gene_name) %>%
  slice_max(len, n = 1, with_ties = FALSE) %>% ungroup()
missing <- setdiff(PSEUDOGENES, query$gene_name)
if (length(missing)) cat("Not found in GENCODE v49 transcripts:", paste(missing, collapse = ", "), "\n")
dir.create(BLAST_WORK, showWarnings = FALSE, recursive = TRUE)
if (grepl(" ", BLAST_WORK)) stop("BLAST_WORK must not contain spaces: ", BLAST_WORK)
writeLines(paste0(">", query$gene_name, "\n", query$seq), file.path(BLAST_WORK, "pseudogene_queries.fa"))

# ---------------------------------------------------------------- 3. database: protein-coding transcripts
db <- tx %>% filter(biotype == "protein_coding", !gene_name %in% PSEUDOGENES)
writeLines(paste0(">", db$tx_id, "\n", db$seq), file.path(BLAST_WORK, "pc_transcripts.fa"))
if (!nzchar(BLAST_DIR)) {                       # look for a BLAST+ installation in the working directory
  cand <- list.files(".", pattern = "^blastn(\\.exe)?$", recursive = TRUE, full.names = TRUE)
  if (length(cand)) BLAST_DIR <- normalizePath(dirname(cand[1]), winslash = "/")
}
bin <- function(x) {
  exe <- if (nzchar(BLAST_DIR)) file.path(BLAST_DIR, x) else unname(Sys.which(x))
  if (.Platform$OS.type == "windows" && !grepl("\\.exe$", exe)) exe <- paste0(exe, ".exe")
  if (!nzchar(exe) || !file.exists(exe))
    stop(x, " not found. Install BLAST+ or set BLAST_DIR at the top of the script to its 'bin' folder.")
  exe
}
MAKEBLASTDB <- bin("makeblastdb"); BLASTN <- bin("blastn")
cat("Using BLAST in:", dirname(BLASTN), "\n")
old_wd <- setwd(BLAST_WORK)
st <- system2(MAKEBLASTDB, c("-in", "pc_transcripts.fa", "-dbtype", "nucl", "-out", "pc_db"))
if (st != 0) { setwd(old_wd); stop("makeblastdb failed (exit code ", st, ")") }

# ---------------------------------------------------------------- 4. BLASTN
st <- system2(BLASTN, c("-task", "blastn", "-query", "pseudogene_queries.fa", "-db", "pc_db",
                        "-dust", "yes", "-evalue", "1e-5", "-max_target_seqs", "5000",
                        "-num_threads", "4",
                        "-outfmt", shQuote("6 qseqid sseqid pident length qlen slen evalue bitscore qcovs", type = "cmd"),
                        "-out", "pseudogene_vs_pc.tsv"))
if (st != 0 || !file.exists("pseudogene_vs_pc.tsv")) { setwd(old_wd); stop("blastn failed (exit code ", st, ")") }
file.copy("pseudogene_vs_pc.tsv", file.path(old_wd, "pseudogene_vs_pc.tsv"), overwrite = TRUE)
setwd(old_wd)
hits <- read.delim("pseudogene_vs_pc.tsv", header = FALSE,
                   col.names = c("pseudogene", "tx_id", "pident", "length", "qlen", "slen",
                                 "evalue", "bitscore", "qcovs")) %>%
  left_join(db %>% select(tx_id, gene_id, gene_name), by = "tx_id")

# best HSP per pseudogene-transcript pair, then thresholds
filt <- hits %>% group_by(pseudogene, tx_id) %>% slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(pident >= MIN_IDENT, qcovs >= MIN_QCOV, evalue <= MAX_EVALUE)

# transcripts collapsed to genes
genes <- filt %>% group_by(pseudogene, gene_name) %>%
  slice_max(bitscore, n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(pseudogene, desc(bitscore))

# ---------------------------------------------------------------- 5. genes overlapping the pseudogene locus (biomaRt)
genes$overlaps_locus <- NA
coords <- tryCatch({
  mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
  biomaRt::getBM(c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                 filters = "external_gene_name",
                 values = unique(c(PSEUDOGENES, genes$gene_name)), mart = mart)
}, error = function(e) { cat("biomaRt failed - overlap check skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(coords)) {
  coords <- coords %>% filter(chromosome_name %in% c(1:22, "X", "Y"))
  ov <- function(a, b) {
    A <- coords[coords$external_gene_name == a, ]; B <- coords[coords$external_gene_name == b, ]
    if (!nrow(A) || !nrow(B)) return(NA)
    any(outer(seq_len(nrow(A)), seq_len(nrow(B)), Vectorize(function(i, j)
      A$chromosome_name[i] == B$chromosome_name[j] &&
        A$start_position[i] <= B$end_position[j] && B$start_position[j] <= A$end_position[i])))
  }
  genes$overlaps_locus <- mapply(ov, genes$pseudogene, genes$gene_name)
  cat("Overlap check done:", sum(genes$overlaps_locus, na.rm = TRUE), "overlapping gene(s) excluded;",
      sum(is.na(genes$overlaps_locus)), "gene(s) without coordinates\n")
}
write.csv(genes, "Pseudogene_BLAST_hits_filtered.csv", row.names = FALSE)

# ---------------------------------------------------------------- 6. parent assignment
assign <- lapply(PSEUDOGENES, function(p) {
  g <- genes %>% filter(pseudogene == p, is.na(overlaps_locus) | !overlaps_locus)
  if (!p %in% query$gene_name) return(tibble(Pseudogene = p, Parent = NA, Status = "not in GENCODE v49"))
  if (!nrow(g)) return(tibble(Pseudogene = p, Parent = NA,
                              Status = "no protein-coding gene passed the thresholds"))
  top <- g[1, ]
  second <- if (nrow(g) > 1) g[2, ] else NULL
  amb <- !is.null(second) && second$bitscore >= AMBIG * top$bitscore
  tibble(Pseudogene = p, Parent = top$gene_name, Identity = top$pident, Query_coverage = top$qcovs,
         Alignment_length = top$length, Query_length = top$qlen, Evalue = top$evalue, Bitscore = top$bitscore,
         Second_gene = if (is.null(second)) NA else second$gene_name,
         Second_bitscore = if (is.null(second)) NA else second$bitscore,
         Status = if (amb) "assigned (ambiguous: paralog scores >= 95% of best)" else "assigned")
}) %>% bind_rows()
# parent used in the original submission, for comparison
original_parent <- c("HMGN2", "TPTE2", "TMEM161B", "VN1R53", "TNRC18", "AIDA", "CPHL1", "RAB6C",
                     "GTF2I", "RAD17", "POTEK", "AZGP1", "PLAC9", "CYP4Z1", "7SL")
assign$Original_parent <- original_parent[match(assign$Pseudogene, PSEUDOGENES)]
assign$Agrees_with_original <- assign$Parent == assign$Original_parent
write.csv(assign, "Pseudogene_parent_assignment.csv", row.names = FALSE)
print(assign %>% select(Pseudogene, Parent, Identity, Query_coverage, Second_gene, Status, Original_parent), width = Inf)

rm(tx, seqs, db, hits); invisible(gc())

# =============================== PART B: DECOUPLING ANALYSIS ===============================
set.seed(42)
B_PERM <- 10000; B_BOOT <- 2000

# ---------------------------------------------------------------- 7. expression data
brca_expr   <- readRDS("BRCA.rda")
normal_expr <- readRDS("Normal.rds")
her2_samples  <- gsub("\\.", "-", readRDS("Her2_Samples.rda"))
basal_samples <- gsub("\\.", "-", readRDS("Basal_Samples.rda"))
colnames(brca_expr)   <- gsub("\\.", "-", colnames(brca_expr))
colnames(normal_expr) <- gsub("\\.", "-", colnames(normal_expr))

# pairs from Part A; HMGN2P15 is the HER2-enriched candidate, the others are Basal-like
pairs <- assign %>% filter(!is.na(Parent)) %>%
  transmute(Pseudogene, Parent, Status,
            Subtype = ifelse(Pseudogene == "HMGN2P15", "HER2-enriched", "Basal-like"))
print(pairs)

# ---------------------------------------------------------------- 8. tests
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
fisher_bw <- function(r1, n1, r2, n2) {   # Bonett-Wright variance for Spearman r
  se <- sqrt((1 + r1^2 / 2) / (n1 - 3) + (1 + r2^2 / 2) / (n2 - 3))
  2 * pnorm(abs(atanh(r2) - atanh(r1)) / se, lower.tail = FALSE)
}
test_pair <- function(ps, pa, norm_cols, tum_expr) {
  xn <- as.numeric(normal_expr[ps, norm_cols]); yn <- as.numeric(normal_expr[pa, norm_cols])
  xt <- as.numeric(tum_expr[ps, ]);             yt <- as.numeric(tum_expr[pa, ])
  rn <- sp(xn, yn); rt <- sp(xt, yt); d <- rt - rn
  n1 <- length(xn); n2 <- length(xt)
  # permutation of normal/tumor labels
  X <- c(xn, xt); Y <- c(yn, yt); N <- n1 + n2
  dperm <- replicate(B_PERM, { i <- sample.int(N, n1); sp(X[-i], Y[-i]) - sp(X[i], Y[i]) })
  p_perm <- (1 + sum(abs(dperm) >= abs(d))) / (B_PERM + 1)
  # bootstrap CI
  dboot <- replicate(B_BOOT, { a <- sample.int(n1, replace = TRUE); b <- sample.int(n2, replace = TRUE)
  sp(xt[b], yt[b]) - sp(xn[a], yn[a]) })
  ci <- quantile(dboot, c(.025, .975), na.rm = TRUE)
  c(n_normal = n1, n_tumor = n2, R_normal = rn, R_tumor = rt, Delta_r = d,
    CI_low = ci[[1]], CI_high = ci[[2]],
    P_fisher_BW = fisher_bw(rn, n1, rt, n2), P_permutation = p_perm)
}

run <- function(exclude_matched) {
  res <- lapply(seq_len(nrow(pairs)), function(i) {
    ps <- pairs$Pseudogene[i]; pa <- pairs$Parent[i]
    samp <- if (pairs$Subtype[i] == "HER2-enriched") her2_samples else basal_samples
    samp <- intersect(samp, colnames(brca_expr))
    if (!all(c(ps, pa) %in% rownames(brca_expr)) || !all(c(ps, pa) %in% rownames(normal_expr))) {
      warning("Skipping ", ps, "-", pa, ": not in expression matrix"); return(NULL) }
    norm_cols <- colnames(normal_expr)
    if (exclude_matched)   # exclude normals from the same patients
      norm_cols <- norm_cols[!substr(norm_cols, 1, 12) %in% substr(samp, 1, 12)]
    cat(ps, "-", pa, if (exclude_matched) "(no shared patients)" else "", "\n")
    data.frame(Subtype = pairs$Subtype[i], Pseudogene = ps, Parent = pa,
               Normals = if (exclude_matched) "unmatched only" else "all",
               t(test_pair(ps, pa, norm_cols, brca_expr[, samp, drop = FALSE])))
  })
  out <- bind_rows(res)
  out$FDR_fisher_BW   <- p.adjust(out$P_fisher_BW, "BH")
  out$FDR_permutation <- p.adjust(out$P_permutation, "BH")
  out
}
cat("Matched normals shared with Basal-like tumors:",
    sum(substr(colnames(normal_expr), 1, 12) %in% substr(basal_samples, 1, 12)), "\n")
cat("Matched normals shared with HER2-enriched tumors:",
    sum(substr(colnames(normal_expr), 1, 12) %in% substr(her2_samples, 1, 12)), "\n")

all_res <- bind_rows(run(FALSE), run(TRUE))
write.csv(all_res, "Pseudogene_Decoupling_Results.csv", row.names = FALSE)
print(all_res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ---------------------------------------------------------------- 9. Figure 5 (all normal samples)
main <- all_res %>% filter(Normals == "all") %>%
  mutate(label = paste0(Pseudogene, " / ", Parent),
         star = ifelse(FDR_permutation < 0.05, "*", ""))
p <- ggplot(main, aes(x = reorder(label, Delta_r), y = Delta_r, fill = Subtype)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.25) +
  geom_text(aes(label = star, y = ifelse(Delta_r >= 0, CI_high, CI_low)),
            hjust = ifelse(main$Delta_r >= 0, -0.6, 1.6), vjust = 0.75, size = 6, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Pseudogene / parent gene",
       y = expression(Delta * r ~ "(" * r[tumor] - r[normal] * ")"), fill = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave("Pseudogene_Decoupling_Plot.png", p, width = 8, height = 6, dpi = 300)