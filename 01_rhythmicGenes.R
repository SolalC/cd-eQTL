# =============================================================================
# 01  Rhythmic gene identification (one tissue per SLURM array task)
#
# For every gene of the tissue, compares a cosinor model to a linear model with
# the same covariates (LRT, 2 df), then lists the variants within 50 kb of the
# TSS of the rhythmic genes for the cd-eQTL scan (02).
#
# Rhythmic genes are called downstream (02, 04) at Bonferroni-adjusted p < 0.05
# within tissue. The p.adj column written here is BH and is used only to select
# the variants passed to 02, a superset of the Bonferroni set.
#
# Run:     sbatch 01_rhythmicGenes.sh        (array 1-50)
# Usage:   Rscript 01_rhythmicGenes.R <tissue index> <n cores>
# Outputs: Results/inferred/harmonicReg/<tissue>_LRT_harmonic_linear.csv
#          Results/inferred/selectedSNPs/<tissue>_SNPs.csv
# Used by: Figure 1b, Supplementary Table 1 (via 04)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(broom)
library(lmtest)
library(parallel)

args  <- commandArgs(trailingOnly = TRUE)
i     <- as.numeric(args[1])
nCore <- as.numeric(args[2])

# Tissue names from the expression matrices:
expressionPath <- list.files(GTEX_EXPR_DIR, pattern = '.bed.gz$', full.names = TRUE)
tissueName <- str_remove(basename(expressionPath), '.v10.normalized_expression.bed.gz')[i]
message('Tissue: ', tissueName)

tod <- loadDIP()

# Variant lookup table:
lookup <- fread(GENOTYPE_LOOKUP)
lookup$variant_pos  <- lookup$pos
lookup$variant_pos2 <- lookup$variant_pos

# Harmonic vs linear LRT for every gene of the tissue:
indv  <- selectIndividual(tod, tissueName)
nGene <- ncol(indv[[1]]) - 1
geneTest <- harmonicRegression(indv, nGene, nCore = nCore,
                               sex = !tissueName %in% SEX_SPECIFIC) %>%
  mutate(tissue = tissueName,
         p.adj  = p.adjust(p.value, method = 'BH'))
dir.create(HARMONIC_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(geneTest, paste0(HARMONIC_DIR, tissueName, '_LRT_harmonic_linear.csv'))

# Variants within the TSS window of the rhythmic genes:
geneTest.f <- geneTest %>% filter(p.adj < 0.05)
SNPs <- extractSNPs(indv, geneTest.f, lookup, TSS_WINDOW)
dir.create(SNP_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(SNPs, paste0(SNP_DIR, tissueName, '_SNPs.csv'))
