# =============================================================================
# Paths and constants shared by every R script. Each script sources this file
# first; it in turn sources 00_helperFunctions.R.
#
# Three environment variables relocate the tree without editing code:
#   CIRCADIAN_ROOT        project root (input data under Data/)
#   CIRCADIAN_RESULTS     intermediate results (default <root>/Results/)
#   CIRCADIAN_SUBMISSION  figures and supplementary tables (default <root>/Writing/Submission/)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

.envDir <- function(var, default) {
  d <- Sys.getenv(var, default)
  if (!endsWith(d, '/')) d <- paste0(d, '/')
  d
}

ROOT        <- .envDir('CIRCADIAN_ROOT', '/scratch/project_mnt/S0007/solal/circadian/')
RESULTS_DIR <- .envDir('CIRCADIAN_RESULTS', paste0(ROOT, 'Results/'))
SUBMISSION  <- .envDir('CIRCADIAN_SUBMISSION', paste0(ROOT, 'Writing/Submission/'))

# --- Inputs --------------------------------------------------------------------
# GTEx v10 (GTEx_Analysis_v10_eQTL_expression_matrices.tar, _covariates.tar)
GTEX_EXPR_DIR   <- '/QRISdata/Q8106/Public/QTL/eQTL_expression_matrices/'
GTEX_COV_DIR    <- '/QRISdata/Q8106/Public/QTL/eQTL_covariates/'
# GTEx WGS genotypes (plink) and the variant lookup table
GENOTYPE_PLINK  <- paste0('/QRISdata/Q8106/Controlled/Genotype/filtered/',
                          'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_cQTL_subset_inferred')
GENOTYPE_LOOKUP <- paste0('/QRISdata/Q8106/Controlled/Genotype/support_files/',
                          'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.lookup_table.txt.gz')
# Donor inferred phase (CHIRAL, Talamanca et al. 2023): R object `phi`, radians
DIP_FILE        <- paste0(ROOT, 'Data/CHIRAL/DIP.RData')
# GTEx v10 sample size per tissue (Figure 1b)
SAMPLE_SIZE     <- paste0(ROOT, 'Data/PhenotypeFile/GTEx_v10_sample_size.csv')
# Associations named by Chen et al. 2025 (Supplementary Table 14)
CHEN_ASSOC      <- paste0(ROOT, 'Data/published/associationListPublished.csv')
# Chen et al. 2025 Supplementary Data 2: rhythmic eQTL tests, one file per tissue
CHEN_SUPPDATA2  <- paste0(RESULTS_DIR, 'published/SupplementaryData2/')
# Ensembl gene ID -> gene symbol cache, written by 03_geneSymbolMap.py
SYMBOL_MAP      <- paste0(ROOT, 'Data/ref/annot/ensembl_hgnc_map.csv')

# --- Intermediate results ------------------------------------------------------
HARMONIC_DIR <- paste0(RESULTS_DIR, 'inferred/harmonicReg/')
SNP_DIR      <- paste0(RESULTS_DIR, 'inferred/selectedSNPs/')
QTL_DIR      <- paste0(RESULTS_DIR, 'inferred/QTLregression/')
RHYTHM_DIR   <- paste0(RESULTS_DIR, 'inferred/rhythmicityParameters/')
DRYR_DIR     <- paste0(RHYTHM_DIR, 'dryR/')
PHEWAS_DIR   <- paste0(RESULTS_DIR, 'phewas/inferred/')
CHEN_DIR     <- paste0(RESULTS_DIR, 'published/')

# --- Submission outputs --------------------------------------------------------
FIG_DIR    <- paste0(SUBMISSION, 'Figures/')
SUPFIG_DIR <- paste0(SUBMISSION, 'Supplementary_Figures/')
TABLE_DIR  <- paste0(SUBMISSION, 'Supplementary_Tables/')
for (.d in c(RHYTHM_DIR, DRYR_DIR, CHEN_DIR, FIG_DIR, SUPFIG_DIR, TABLE_DIR))
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# --- Analysis constants --------------------------------------------------------
TSS_WINDOW    <- 50000   # bp either side of the TSS (100 kb window)
MAF_THRESHOLD <- 0.1     # minimum minor allele frequency in the cd-eQTL scan
CDEQTL_P      <- 1e-4    # interaction LRT threshold defining a cd-eQTL
MIN_GROUP_N   <- 10      # minimum donors per genotype group (Bingham tests, dryR)
CLASS_FDR     <- 0.1     # BH FDR for the amplitude / acrophase classification
SEX_SPECIFIC  <- c('Ovary', 'Prostate', 'Testis', 'Uterus', 'Vagina')
N_PHENOTYPES  <- 1085    # Pan-UK Biobank phenotypes in the PheWAS
PHEWAS_P      <- 0.05 / N_PHENOTYPES   # 4.6e-5

# Core circadian genes (CCG): BMAL1, CRY1, NFIL3, CIART, PER3, BHLHE41, DBP
CCG <- c('ENSG00000133794.20', 'ENSG00000008405.12', 'ENSG00000165030.4',
         'ENSG00000159208.16', 'ENSG00000049246.15', 'ENSG00000123095.6',
         'ENSG00000105516.11')

# The two example loci of Figures 3-5 and Supplementary Tables 7-12
EXAMPLES <- tibble::tribble(
  ~symbol,  ~gene,                ~variant,     ~tissue,         ~type,
  'BMAL1',  'ENSG00000133794.20', 'rs11022718', 'Artery_Aorta',  'Core Circadian Gene',
  'TMED10', 'ENSG00000170348.9',  'rs76379942', 'Colon_Sigmoid', 'Circadian Regulated Gene')

# --- Shared functions ----------------------------------------------------------
.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
CODE_DIR <- if (length(.f)) dirname(normalizePath(.f)) else '.'
source(file.path(CODE_DIR, '00_helperFunctions.R'))
