# =============================================================================
# 14  Interaction model of each independent cd-eQTL in every tissue
#
# For each of the 54 gene-variant pairs, fits the harmonic interaction model in
# every GTEx tissue expressing the gene, derives the per-genotype amplitude and
# acrophase, and computes the tissue x tissue cosine similarity of the four
# rhythm coefficients (sin, cos, variant:sin, variant:cos) over the 49 tissues
# with at least one rhythmic gene. Loci already fitted are skipped.
#
# Run:     sbatch 14_multiTissueRegression.sh    (after 04)
# Outputs: Results/inferred/rhythmicityParameters/multiTissueReg/<gene>_<variant>.csv
#          Results/inferred/rhythmicityParameters/rhythmicityParam/<gene>_<variant>.csv
#          Results/inferred/rhythmicityParameters/cosineSimilarity/<gene>_<variant>.csv
# Used by: Figures 3c-f and 4, Supplementary Figure 3, Supplementary Tables 11-12
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lsa)

dirs <- paste0(RHYTHM_DIR, c('multiTissueReg/', 'rhythmicityParam/', 'cosineSimilarity/'))
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

tod      <- loadDIP() %>% filter(!is.na(radian))
plink    <- loadGenotypes()$X
cdeqtl   <- readIndependentCdeQTL()
RHYTHMIC <- rhythmicTissues()

done <- str_remove(list.files(dirs[1]), '\\.csv$')
for (i in seq_len(nrow(cdeqtl))) {
  toRun <- cdeqtl[i, ]
  key   <- paste0(toRun$gene, '_', toRun$variant)
  if (key %in% done) next
  message('Fitting ', key)
  rhythmicityParameters <- multiTissueReg(tod, plink, toRun$gene, toRun$variant)
  fwrite(rhythmicityParameters, paste0(dirs[1], key, '.csv'))
  fwrite(calculate_rhythmicity_parameters(rhythmicityParameters), paste0(dirs[2], key, '.csv'))
  similarity <- cosineSimilarity(rhythmicityParameters, RHYTHMIC)
  fwrite(as.data.frame(similarity[[1]]), paste0(dirs[3], key, '.csv'))
}
