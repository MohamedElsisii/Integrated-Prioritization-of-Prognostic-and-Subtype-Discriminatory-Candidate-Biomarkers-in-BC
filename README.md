# Integrated Prioritization of Prognostic and Molecular Subtype-Discriminatory Candidate Biomarkers in Breast Cancer

<img width="3600" height="2120" alt="workflow" src="https://github.com/user-attachments/assets/80d0cc5a-0196-46a9-ae57-58de5fa67764" />

## Overview

This repository contains the complete analysis code for the study
**"Integrated Prioritization of Prognostic and Molecular Subtype-Discriminatory Candidate Biomarkers in Breast Cancer"** (Elsisi and Maged).

Breast cancer (BC) is a heterogeneous disease. Molecular subtyping with the PAM50 classifier has improved its classification, but patient outcomes still vary widely within each subtype. This integrated bioinformatics framework prioritizes subtype-specific candidate biomarkers and characterizes their regulatory and functional context.

The workflow combines differential expression, survival and Cox regression analyses, cross-validated receiver operating characteristic (ROC) analysis, correlation with PAM50 genes, biotype annotation, pseudogene–parent correlation analysis, and pathway-level functional rewiring. It prioritizes **193 candidate biomarkers** across the four intrinsic subtypes (Luminal A, Luminal B, Basal-like, and HER2-enriched):

- **Prognostic candidates (exploratory):** F2RL2 (Luminal A); GLYATL2, ST8SIA2, PID1, and LRIG3 (Luminal B); TRIM59, LINC01615, DLL3, and CIP2A (Basal-like)
- **Molecular subtype-discriminatory candidates:** FOXC1 and STAC (Luminal B); 181 genes (Basal-like); HMGN2P15 (HER2-enriched, low confidence)

The candidates were externally evaluated in the SCAN-B cohort (GEO: GSE96058). The molecular subtype-discriminatory candidates were reproduced, whereas the prognostic candidates were only partly replicated. All findings are hypothesis-generating and require further validation before clinical use.

## Workflow

1. **Data acquisition and subtyping:** TCGA-BRCA RNA-seq counts and clinical data (GDC); PAM50 subtyping with genefu.
2. **Differential expression:** each subtype versus normal breast tissue (edgeR, limma `voomWithQualityWeights`, `treat`, |logFC| ≥ 2, BH < 0.05). Genes significant in exactly one subtype are retained as subtype-specific (unique) DEGs.
3. **Prognostic screening:** Kaplan–Meier analysis at the optimal cut-point, with p-values corrected for cut-point selection (maximally selected rank statistics, permutation test).
4. **Cox regression:** univariate and multivariate models (age, stage); proportional hazards testing with pre-specified remedies; continuous-expression and covariate-adjusted sensitivity analyses (treatment, ER/PR/HER2 status).
5. **Molecular subtype-discriminatory evaluation:** one-versus-rest ROC analysis with repeated stratified 10-fold cross-validation.
6. **Prioritization:** an integrated prioritization score (IPS) for the Basal-like candidates, with a sensitivity analysis of its weighting.
7. **Functional characterization:** correlation with PAM50 genes, cancer hallmark enrichment, pseudogene–parent decoupling, and Hallmark pathway rewiring (GSVA).
8. **External validation:** SCAN-B (GSE96058).

## Scripts

| Script | Purpose | Main outputs |
|---|---|---|
| `Code_S1.R` | PAM50 subtyping and differential expression | Tables S1–S6 |
| `Code_S2.R` | Kaplan–Meier screening of unique DEGs; cut-point-corrected (maxstat) and permutation p-values. Run once per subtype (`SUBTYPE <- "LumA"`, `"LumB"`, `"Basal"`) | Tables S7–S10, Files S1–S3 |
| `Code_S3.R` | Univariate and multivariate Cox regression | Tables S11–S18, Figures S1–S3 |
| `Code_S4.R` | Final Cox models of the prognostic candidates (PH remedies, unadjusted/age/stage models, continuous expression) and ROC curves | Table S38, Figures S4–S9 |
| `Code_S5.R` | One-versus-rest ROC with repeated stratified 10-fold cross-validation, including a fixed-direction sensitivity analysis | Tables S19–S22, Figures 3–4, File S4 |
| `Code_S6.R` | Correlation of candidates with PAM50 genes | Tables S23–S25 |
| `Code_S7.R` | Pseudogene–parent assignment (BLASTN, explicit thresholds) and decoupling analysis (permutation test, bootstrap CIs) | Tables S28–S29, Figure 5 |
| `Code_S8.R` | Pathway-level functional rewiring of lncRNAs (GSVA, MSigDB Hallmark) | Tables S30–S31 |
| `Code_S9.R` | Pathway-level functional rewiring of protein-coding genes (GSVA, MSigDB Hallmark) | Tables S32–S37 |
| `Code_S10.R` | Sensitivity analysis of the integrated prioritization score | Table S40 |
| `Code_S11.R` | Adjustment of the final Cox models for treatment and receptor status | Table S42 |
| `Code_S12.R` | External validation in SCAN-B (GSE96058) | Table S43 |

Cancer hallmark enrichment (Table S41) was performed with the web tool [cancerhallmarks.com](https://www.cancerhallmarks.com) using the integrated cancer hallmark gene set (n = 6,763 genes).

Each script begins with a header describing its input files and outputs. Scripts read their inputs from the working directory.

## Data

- **TCGA-BRCA:** GDC Data Portal (https://portal.gdc.cancer.gov/), STAR – Counts (GENCODE v36), Data Release 45.0, accessed June 2026; retrieved with TCGAbiolinks. Sample identifiers and PAM50 assignments are provided in Table S1.
- **GENCODE v49** transcript sequences (Code S7): https://www.gencodegenes.org/human/
- **SCAN-B:** GEO accession GSE96058 (downloaded automatically by Code S12).

## Requirements

- R ≥ 4.6 with Bioconductor packages, including edgeR, limma, DESeq2, genefu, TCGAbiolinks, biomaRt, org.Hs.eg.db, GSVA, GEOquery, survival, survminer, maxstat, pROC, data.table, tidyverse, and ggpubr. Exact versions are listed in `sessionInfo.txt`.
- BLAST+ 2.17.0 (`blastn`, `makeblastdb`) for Code S7.

## Citation

If you use this code, please cite:

> Elsisi M, Maged M. Integrated Prioritization of Prognostic and Molecular Subtype-Discriminatory Candidate Biomarkers in Breast Cancer. *Network Modeling Analysis in Health Informatics and Bioinformatics* (under review).

## Contact

Mohamad Maged (corresponding author): momaged@nu.edu.eg
