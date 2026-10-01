# =============================================================================
# 02  cd-eQTL scan (one tissue per SLURM array task)
#
# For every rhythmic gene of the tissue (Bonferroni-adjusted harmonic LRT
# p < 0.05, from 01) and every variant within 50 kb of its TSS with MAF >= 10%,
# fits three nested models and compares them with sequential LRTs:
#   1. gene ~ variant + covariates
#   2. gene ~ variant + covariates + sin + cos
#   3. gene ~ variant + covariates + sin + cos + variant:sin + variant:cos
# The model 3 vs model 2 LRT is the cd-eQTL test (threshold applied in 04).
#
# Run:     sbatch 02_cdeQTLscan.sh           (array 1-50; after 01)
# Usage:   Rscript 02_cdeQTLscan.R <tissue index> <n cores>
# Outputs: Results/inferred/QTLregression/{1_linear,2_harmonic,3_harmonicInteraction,4_LRT}/
# Used by: Figure 2, Supplementary Tables 2-4 (via 04, 05, 06)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)
library(future)
library(future.apply)

args  <- commandArgs(trailingOnly = TRUE)
i     <- as.numeric(args[1])
nCore <- as.numeric(args[2])
plan('multisession', workers = nCore)

expressionPath <- list.files(GTEX_EXPR_DIR, pattern = '.bed.gz$', full.names = TRUE)
tissueName <- str_remove(basename(expressionPath), '.v10.normalized_expression.bed.gz')[i]
message('Tissue: ', tissueName)

tod <- loadDIP()

# Rhythmic genes of the tissue (Bonferroni within tissue):
harmonic <- fread(paste0(HARMONIC_DIR, tissueName, '_LRT_harmonic_linear.csv')) %>%
  mutate(p.adj = p.adjust(p.value, method = 'bonferroni')) %>%
  filter(p.adj < 0.05)
if (nrow(harmonic) == 0) stop(paste0('No significant genes found in ', tissueName))

# Restrict expression to those genes:
indv <- selectIndividual(tod, tissueName)
idx <- which(colnames(indv$tissue) %in% c(harmonic$gene, 'SUBJID'))
indv$tissue <- indv$tissue[, ..idx]
idx <- which(indv$pos$gene_id %in% harmonic$gene)
indv$pos <- indv$pos[idx, ]

plink <- loadGenotypes()
SNPs <- fread(paste0(SNP_DIR, tissueName, '_SNPs.csv')) %>%
  filter(gene_id %in% harmonic$gene)

# Run in chunks of 20 genes:
geneToRun <- unique(SNPs$gene_id)
chunked <- split(seq_along(geneToRun), ceiling(seq_along(geneToRun) / 20))
linear <- data.frame(); harmonic <- data.frame()
harmonicInteraction <- data.frame(); lrt <- data.frame()
for (genes in seq_along(chunked)) {
  SNPs.f <- SNPs %>% filter(gene_id %in% geneToRun[chunked[[genes]]])
  output <- QTLregressionParallel(tissueName, indv, plink, SNPs.f, MAF_THRESHOLD,
                                  sex = !tissueName %in% SEX_SPECIFIC)
  linear              <- rbind(linear, output$linear)
  harmonic            <- rbind(harmonic, output$harmonic)
  harmonicInteraction <- rbind(harmonicInteraction, output$harmonicInteraction)
  lrt                 <- rbind(lrt, output$lrt)
}

out <- list(`1_linear` = linear, `2_harmonic` = harmonic,
            `3_harmonicInteraction` = harmonicInteraction, `4_LRT` = lrt)
suffix <- c(`1_linear` = '_linear', `2_harmonic` = '_harmonic',
            `3_harmonicInteraction` = '_harmonicInteraction', `4_LRT` = '_LRT')
for (d in names(out)) {
  dir.create(paste0(QTL_DIR, d), recursive = TRUE, showWarnings = FALSE)
  fwrite(out[[d]], paste0(QTL_DIR, d, '/', tissueName, suffix[[d]], '.csv'))
}
