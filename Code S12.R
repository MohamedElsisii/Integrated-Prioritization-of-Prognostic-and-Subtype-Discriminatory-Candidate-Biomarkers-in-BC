# =============================================================================
# Code S12 - External validation in SCAN-B (GEO: GSE96058)
#
# GSE96058: 3,273 primary breast tumors with RNA-seq, PAM50 subtype and overall survival;
# technical replicates and Normal-like tumors are excluded.
#   (A) Prognostic candidates: Cox models within the PAM50 subtype with expression as a
#       continuous variable (HR per 1 SD), stratified by sequencing platform, unadjusted
#       and adjusted for age, tumor size and lymph node status; median split for comparison.
#       Replicated = same direction as in TCGA-BRCA and p < 0.05.
#   (B) Molecular subtype-discriminatory candidates: one-versus-rest AUC (DeLong 95% CI)
#       against SCAN-B PAM50 labels, with the direction of the comparison fixed from TCGA-BRCA.
# Genes absent under their current symbol are matched through previous symbols (org.Hs.eg.db).
#
# Input (working directory)
#   Table S19.csv - Table S22.csv (TCGA ROC results), Table S27.csv (candidates and biotypes)
#   GSE96058 clinical data and expression file are downloaded on the first run.
# Packages: data.table, dplyr, survival, pROC, org.Hs.eg.db, GEOquery
# Output: GSE96058_prognostic_validation.csv, GSE96058_discrimination_validation.csv (Table S43)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(survival); library(pROC) })

GSE_DIR <- if (.Platform$OS.type == "windows") "C:/gse96058" else "gse96058"
dir.create(GSE_DIR, showWarnings = FALSE, recursive = TRUE)
expr_file <- file.path(GSE_DIR, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
clin_file <- "GSE96058_clinical.csv"
options(timeout = 3600)

# TCGA-BRCA final-model hazard ratios (Table S38)
prognostic <- tibble::tribble(
  ~subtype, ~gene,       ~tcga_HR,
  "LumA",   "F2RL2",     0.42,
  "LumB",   "GLYATL2",   10.1,
  "LumB",   "ST8SIA2",   9.42,
  "LumB",   "PID1",      4.20,
  "LumB",   "LRIG3",     4.17,
  "Basal",  "TRIM59",    0.17,
  "Basal",  "LINC01615", 15.5,
  "Basal",  "DLL3",      2.89,
  "Basal",  "CIP2A",     0.25)

# -----------------------------------------------------------------------------
# 1. Clinical data (GEO series matrix)
# -----------------------------------------------------------------------------
if (!file.exists(clin_file)) {
  library(GEOquery)
  gse <- getGEO("GSE96058", GSEMatrix = TRUE, destdir = GSE_DIR, getGPL = FALSE)
  ph  <- bind_rows(lapply(gse, function(e) {
    p <- pData(e); p$platform_id <- as.character(p$platform_id); mutate(p, across(everything(), as.character)) }))
  write.csv(ph, clin_file, row.names = FALSE)
}
ph <- read.csv(clin_file, check.names = FALSE, stringsAsFactors = FALSE)
col <- function(pattern) { m <- grep(pattern, names(ph), ignore.case = TRUE, value = TRUE); if (length(m)) m[1] else NA }
cols <- c(title = "^title$", platform = "^platform_id$", pam50 = "pam50.*ch1", age = "^age.*ch1",
          size = "tumor size.*ch1", node = "lymph node status.*ch1", os_days = "overall survival days.*ch1",
          os_event = "overall survival event.*ch1")
found <- sapply(cols, col); print(found)
if (any(is.na(found[c("title", "pam50", "os_days", "os_event")]))) stop("Required clinical columns not found in GSE96058_clinical.csv")

clin <- tibble(sample   = ph[[found["title"]]],
               platform = ph[[found["platform"]]],
               pam50    = ph[[found["pam50"]]],
               age      = suppressWarnings(as.numeric(ph[[found["age"]]])),
               size     = if (is.na(found["size"])) NA_real_ else suppressWarnings(as.numeric(ph[[found["size"]]])),
               node_raw = if (is.na(found["node"])) NA_character_ else ph[[found["node"]]],
               os_time  = suppressWarnings(as.numeric(ph[[found["os_days"]]])),
               os_event = suppressWarnings(as.numeric(ph[[found["os_event"]]]))) %>%
  filter(!grepl("repl", sample, ignore.case = TRUE)) %>%   # technical replicates
  mutate(node = case_when(grepl("NodeNegative|negative|^N0|^0$", node_raw, ignore.case = TRUE) ~ 0,
                          grepl("NodePositive|positive|N[1-3]|SubMicro", node_raw, ignore.case = TRUE) ~ 1),
         pam50 = recode(pam50, Her2 = "Her2", HER2 = "Her2", Basal = "Basal", LumA = "LumA", LumB = "LumB",
                        Normal = "Normal"))
cat("\nPAM50 in GSE96058:\n"); print(table(clin$pam50, useNA = "ifany"))
cat("Lymph node coding:\n"); print(table(clin$node_raw, clin$node, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 2. Expression (log2-transformed FPKM as provided); candidate genes only
# -----------------------------------------------------------------------------
if (!file.exists(expr_file))
  download.file(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/", basename(expr_file)),
                expr_file, mode = "wb")
s27 <- read.csv("Table S27.csv", stringsAsFactors = FALSE)
disc <- s27 %>% filter(Biomarker != "Prognostic") %>%
  mutate(subtype = recode(Subtype, "Luminal A" = "LumA", "Luminal B" = "LumB", "Basal-like" = "Basal",
                          "HER2-enriched" = "Her2")) %>% select(subtype, Gene, Biotype)
need <- unique(c(prognostic$gene, disc$Gene))

ex <- fread(expr_file, data.table = FALSE)
genes <- ex[[1]]
# SCAN-B uses a 2014 annotation; genes missing under their current symbol are matched through
# unambiguous previous symbols (e.g. CIP2A = KIAA1524, SRARP = C1orf64)
alias_map <- setNames(need[need %in% genes], need[need %in% genes])
miss <- setdiff(need, genes)
if (length(miss) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  ok <- miss[miss %in% AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")]
  eg <- if (length(ok)) suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ok,
                                                               keytype = "SYMBOL", columns = "ENTREZID")) else data.frame(ENTREZID = character())
  eg <- eg[!is.na(eg$ENTREZID), ]
  if (nrow(eg)) {
    al <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(eg$ENTREZID),
                                                 keytype = "ENTREZID", columns = "ALIAS"))
    al <- merge(eg, al, by = "ENTREZID")
    official <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
    al <- al[al$ALIAS %in% genes & !al$ALIAS %in% need &   # alias present in SCAN-B, not another candidate,
               !al$ALIAS %in% official, ]  # and not the current symbol of another gene
    al <- al[!duplicated(al$SYMBOL) & !duplicated(al$ALIAS), ]
    alias_map <- c(alias_map, setNames(al$ALIAS, al$SYMBOL))
    if (nrow(al)) { cat("Genes found under an older symbol:\n"); print(al[, c("SYMBOL", "ALIAS")], row.names = FALSE) }
  }
}
ex <- ex[match(alias_map, genes), , drop = FALSE]
rownames(ex) <- names(alias_map); ex <- as.matrix(ex[, -1])
ex <- ex[, colnames(ex) %in% clin$sample, drop = FALSE]
clin <- clin %>% filter(sample %in% colnames(ex), pam50 %in% c("LumA", "LumB", "Her2", "Basal"))
ex <- ex[, clin$sample, drop = FALSE]
cat("\nPatients analysed (Normal-like and replicates excluded):", ncol(ex), "\n")
cat("Candidate genes found in GSE96058:", sum(need %in% rownames(ex)), "of", length(need), "\n")
cat("Not found:", paste(setdiff(need, rownames(ex)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 3. (A) Prognostic validation
# -----------------------------------------------------------------------------
hr <- function(fit, term) { s <- summary(fit); i <- grep(term, rownames(s$conf.int))[1]
c(HR = s$conf.int[i, 1], lower = s$conf.int[i, 3], upper = s$conf.int[i, 4], p = s$coefficients[i, ncol(s$coefficients)]) }
prog_res <- lapply(seq_len(nrow(prognostic)), function(i) {
  g <- prognostic$gene[i]; st <- prognostic$subtype[i]
  if (!g %in% rownames(ex)) return(tibble(subtype = st, gene = g, note = "gene not in GSE96058"))
  d <- clin %>% mutate(expr = ex[g, ]) %>% filter(pam50 == st, !is.na(os_time), !is.na(os_event), os_time > 0)
  d$expr_sd <- as.numeric(scale(d$expr))
  d$high    <- as.numeric(d$expr > median(d$expr))
  m_un  <- coxph(Surv(os_time, os_event) ~ expr_sd + strata(platform), data = d)
  da    <- d %>% filter(!is.na(age), !is.na(size), !is.na(node))
  m_adj <- coxph(Surv(os_time, os_event) ~ expr_sd + age + size + node + strata(platform), data = da)
  m_med <- coxph(Surv(os_time, os_event) ~ high + strata(platform), data = d)
  u <- hr(m_un, "expr_sd"); a <- hr(m_adj, "expr_sd"); md <- hr(m_med, "high")
  tibble(subtype = st, gene = g, tcga_HR = prognostic$tcga_HR[i],
         n = nrow(d), deaths = sum(d$os_event), n_adj = nrow(da), deaths_adj = sum(da$os_event),
         HR_per_SD = u[["HR"]], lower = u[["lower"]], upper = u[["upper"]], p = u[["p"]],
         HR_per_SD_adj = a[["HR"]], lower_adj = a[["lower"]], upper_adj = a[["upper"]], p_adj = a[["p"]],
         HR_median_split = md[["HR"]], p_median_split = md[["p"]],
         same_direction = sign(log(u[["HR"]])) == sign(log(prognostic$tcga_HR[i])),
         replicated = same_direction & u[["p"]] < 0.05, note = "")
}) %>% bind_rows()
write.csv(prog_res, "GSE96058_prognostic_validation.csv", row.names = FALSE)
cat("\n==== (A) PROGNOSTIC VALIDATION (HR per 1 SD, stratified by platform) ====\n")
print(prog_res %>% transmute(subtype, gene, n, deaths, tcga_HR,
                             `HR/SD (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR_per_SD, lower, upper), p = signif(p, 2),
                             `adj HR/SD` = sprintf("%.2f (%.2f-%.2f)", HR_per_SD_adj, lower_adj, upper_adj),
                             p_adj = signif(p_adj, 2), HR_median = round(HR_median_split, 2),
                             p_median = signif(p_median_split, 2), same_direction, replicated, note), width = Inf)

# -----------------------------------------------------------------------------
# 4. (B) Discrimination validation (direction fixed from TCGA-BRCA)
# -----------------------------------------------------------------------------
tcga <- bind_rows(lapply(c("Table S19.csv", "Table S20.csv", "Table S21.csv", "Table S22.csv"), function(f)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)[, c("Subtype", "Gene", "Mean_Diff", "AUC_CV")])) %>%
  filter(!is.na(Gene)) %>% rename(subtype = Subtype)
disc <- disc %>% left_join(tcga, by = c("subtype", "Gene"))
disc_res <- lapply(seq_len(nrow(disc)), function(i) {
  g <- disc$Gene[i]; st <- disc$subtype[i]
  if (!g %in% rownames(ex) || is.na(disc$Mean_Diff[i]))
    return(tibble(subtype = st, gene = g, biotype = disc$Biotype[i], tcga_AUC = disc$AUC_CV[i], note = "gene not in GSE96058"))
  y <- as.numeric(clin$pam50 == st)
  r <- roc(y, ex[g, ], levels = c(0, 1), direction = ifelse(disc$Mean_Diff[i] > 0, "<", ">"), quiet = TRUE)
  ci <- as.numeric(ci.auc(r))
  tibble(subtype = st, gene = g, biotype = disc$Biotype[i], tcga_AUC = disc$AUC_CV[i],
         n_pos = sum(y), n_neg = sum(1 - y), AUC = ci[2], lower = ci[1], upper = ci[3],
         mean_diff_gse = mean(ex[g, y == 1]) - mean(ex[g, y == 0]),
         same_direction = sign(mean_diff_gse) == sign(disc$Mean_Diff[i]), note = "")
}) %>% bind_rows()
write.csv(disc_res, "GSE96058_discrimination_validation.csv", row.names = FALSE)

cat("\n==== (B) DISCRIMINATION VALIDATION (one-vs-rest AUC, direction fixed from TCGA) ====\n")
print(disc_res %>% filter(subtype != "Basal") %>%
        transmute(subtype, gene, tcga_AUC = round(tcga_AUC, 3), AUC = round(AUC, 3),
                  CI = sprintf("%.3f-%.3f", lower, upper), same_direction, note))
b <- disc_res %>% filter(subtype == "Basal", note == "")
cat(sprintf("\nBasal-like: %d of %d candidates found in GSE96058\n", nrow(b), sum(disc_res$subtype == "Basal")))
cat(sprintf("  same direction of difference as TCGA: %d (%.1f%%)\n", sum(b$same_direction), 100 * mean(b$same_direction)))
cat(sprintf("  AUC >= 0.80: %d (%.1f%%); AUC >= 0.70: %d (%.1f%%)\n", sum(b$AUC >= 0.8), 100 * mean(b$AUC >= 0.8),
            sum(b$AUC >= 0.7), 100 * mean(b$AUC >= 0.7)))
cat(sprintf("  median AUC %.3f (IQR %.3f-%.3f); TCGA median %.3f\n", median(b$AUC), quantile(b$AUC, .25),
            quantile(b$AUC, .75), median(b$tcga_AUC)))
cat(sprintf("  Spearman rho TCGA vs GSE96058 AUC: %.2f\n", cor(b$AUC, b$tcga_AUC, method = "spearman")))
print(b %>% group_by(biotype) %>% summarise(n = n(), median_AUC = round(median(AUC), 3), AUC_ge_0.8 = sum(AUC >= 0.8)))
top20 <- c("TTC6", "ERBB4", "MLPH", "SRARP", "ANKRD30A", "DNALI1", "THSD4", "SPDEF", "TPSG1", "CDCA7", "CA12",
           "TFF3", "CCDC170", "FZD9", "SCUBE2", "TBC1D9", "INPP4B", "CAPN8", "ANKRD30B", "ABCC8")
cat("\nTable 2 genes (top 20 by IPS):\n")
print(b %>% filter(gene %in% top20) %>% arrange(match(gene, top20)) %>%
        transmute(gene, tcga_AUC = round(tcga_AUC, 3), AUC = round(AUC, 3), CI = sprintf("%.3f-%.3f", lower, upper)))