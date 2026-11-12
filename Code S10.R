# =============================================================================
# Code S10 - Sensitivity analysis of the integrated prioritization score (IPS)
#
# The IPS of the 181 Basal-like molecular subtype-discriminatory candidates is the
# sum of four z-scored components (cross-validated AUC, absolute mean log2
# expression difference, number of rewired Hallmark pathways, number of significant
# PAM50 correlations), rescaled to 0-100. The ranking is recomputed
#   (a) omitting each component in turn,
#   (b) with percentile ranks instead of z-scores,
#   (c) under 1,000 random weightings (flat Dirichlet, weights summing to four),
# and compared with the original ranking (Spearman rho, top-20 overlap).
#
# Input (working directory): Table S21.csv, Table S25.csv, Table S27.csv,
#                            Table S30.csv, Table S37.csv
# Output: IPS_sensitivity_summary.csv (Table S40A), IPS_all_candidates.csv (Table S40B)
# =============================================================================

library(dplyr)
library(readr)
set.seed(42)

# ---------------------------------------------------------------- 1. components
s21 <- read_csv("Table S21.csv", show_col_types = FALSE)
roc <- s21 %>% filter(!is.na(Gene)) %>% select(Gene, AUC_CV, Mean_Diff)

biotype <- read_csv("Table S27.csv", show_col_types = FALSE)
# Basal-like molecular subtype-discriminatory candidates (non-prognostic rows of Table S27)
candidates <- biotype %>% filter(Subtype == "Basal-like", Biomarker != "Prognostic") %>% pull(Gene) %>% unique()
biotype <- biotype %>% filter(Gene %in% candidates) %>% distinct(Gene, .keep_all = TRUE) %>% select(Gene, Biotype)

s37 <- read_csv("Table S37.csv", show_col_types = FALSE)[, 1:9] %>% filter(!is.na(Gene)) %>%
  transmute(Gene, rewired = Num_Gained + Num_Lost)
s30 <- read_csv("Table S30.csv", show_col_types = FALSE)[, 1:7] %>% filter(!is.na(lncRNA)) %>%
  transmute(Gene = lncRNA, rewired = Num_Gained + Num_Lost)
rewiring <- bind_rows(s37, s30)

s25 <- read_csv("Table S25.csv", show_col_types = FALSE)
pam50 <- s25 %>% count(Gene = .[[1]], name = "pam50")


comp <- tibble(Gene = candidates) %>%
  left_join(roc, by = "Gene") %>%
  left_join(rewiring, by = "Gene") %>%
  left_join(pam50, by = "Gene") %>%
  left_join(biotype, by = "Gene") %>%
  mutate(auc = AUC_CV,
         abs_diff = abs(Mean_Diff),               # |mean log2 expression difference| vs other subtypes
         rewired = coalesce(rewired, 0),          # pseudogenes were not in the rewiring analysis -> 0
         pam50 = coalesce(pam50, 0))
cat("Candidates:", nrow(comp), "\n")

# ---------------------------------------------------------------- 2. the IPS function
M <- as.matrix(comp[, c("auc", "abs_diff", "rewired", "pam50")])
Z <- scale(M)
ips <- function(w = c(1, 1, 1, 1), Zmat = Z) {
  s <- as.numeric(Zmat %*% w)
  100 * (s - min(s)) / (max(s) - min(s))          # rescale to 0-100
}
comp$IPS <- ips()

# top protein-coding candidates (Table 2)
print(comp %>% filter(Biotype == "protein_coding") %>% arrange(desc(IPS)) %>%
        select(Gene, IPS) %>% head(5))

# ---------------------------------------------------------------- 3. helpers
is_pc <- comp$Biotype == "protein_coding"
top20 <- function(score, only_pc = FALSE) {
  keep <- if (only_pc) is_pc else rep(TRUE, length(score))
  g <- comp$Gene[keep]; s <- score[keep]
  g[order(-s)][1:20]
}
base_all <- top20(comp$IPS); base_pc <- top20(comp$IPS, TRUE)     # base_pc = Table 2
ttc6_rank <- function(score) rank(-score, ties.method = "min")[comp$Gene == "TTC6"]
compare <- function(name, score) tibble(
  Scheme = name,
  Spearman_rho = cor(comp$IPS, score, method = "spearman"),
  Top20_overlap_all = length(intersect(top20(score), base_all)),
  Top20_overlap_proteincoding = length(intersect(top20(score, TRUE), base_pc)),
  TTC6_rank = ttc6_rank(score))

# ---------------------------------------------------------------- 4. alternative schemes
res <- compare("Original (equal weights, z-scores)", comp$IPS)
names4 <- c("AUC", "|log2 difference|", "Rewired pathways", "PAM50 correlations")
for (j in 1:4) {                                  # (a) leave one component out
  w <- c(1, 1, 1, 1); w[j] <- 0
  res <- bind_rows(res, compare(paste("Omit", names4[j]), ips(w)))
}
R <- apply(M, 2, function(x) rank(x) / length(x))  # (b) percentile ranks instead of z-scores
res <- bind_rows(res, compare("Percentile ranks instead of z-scores", ips(Zmat = scale(R))))

# (c) 1,000 random weightings (flat Dirichlet = random weights that sum to 4)
B <- 1000
rho <- ov_all <- ov_pc <- ttc <- numeric(B)
count_all <- count_pc <- setNames(numeric(nrow(comp)), comp$Gene)
for (b in 1:B) {
  w <- rgamma(4, 1); w <- 4 * w / sum(w)
  s <- ips(w)
  rho[b] <- cor(comp$IPS, s, method = "spearman")
  ta <- top20(s); tp <- top20(s, TRUE)
  ov_all[b] <- length(intersect(ta, base_all)); ov_pc[b] <- length(intersect(tp, base_pc))
  ttc[b] <- ttc6_rank(s)
  count_all[ta] <- count_all[ta] + 1; count_pc[tp] <- count_pc[tp] + 1
}
res <- bind_rows(res, tibble(
  Scheme = "1,000 random weightings: median (IQR)",
  Spearman_rho = median(rho),
  Top20_overlap_all = median(ov_all), Top20_overlap_proteincoding = median(ov_pc),
  TTC6_rank = median(ttc)))
print(res, width = Inf)
cat(sprintf("rho IQR %.3f-%.3f | TTC6 in top 5: %.1f%% | in top 20: %.1f%%\n",
            quantile(rho, .25), quantile(rho, .75), 100 * mean(ttc <= 5), 100 * mean(ttc <= 20)))

# how often each Table 2 gene stays in the protein-coding top 20
print(round(100 * count_pc[base_pc] / B, 1))

# correlation between components
print(round(cor(M, method = "spearman"), 2))

# ---------------------------------------------------------------- 5. save
write_csv(res, "IPS_sensitivity_summary.csv")
comp %>% mutate(pct_top20_all = 100 * count_all[Gene] / B,
                pct_top20_proteincoding = ifelse(is_pc, 100 * count_pc[Gene] / B, NA)) %>%
  arrange(desc(IPS)) %>% mutate(Rank = row_number(), .before = 1) %>%
  select(Rank, Gene, Biotype, auc, abs_diff, rewired, pam50, IPS, pct_top20_all, pct_top20_proteincoding) %>%
  write_csv("IPS_all_candidates.csv")