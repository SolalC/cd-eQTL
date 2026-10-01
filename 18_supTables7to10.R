# =============================================================================
# 18  Supplementary Tables 7-10: amplitude and acrophase of BMAL1 rs11022718
#     and TMED10 rs76379942 by tissue and genotype group, across the 49
#     rhythmic tissues
#
# From the residualised per-genotype cosinor (the fits compared by the Bingham
# tests, 07); groups with fewer than 10 donors are not fitted. Acrophase is the
# peak time of the fitted curve, atan2(sin, cos), in hours (0-24).
#
# Run:     sbatch 18_supTables7to10.sh     (after 01)
# Outputs: Supplementary_Tables/SupTable_7_BMAL1_amplitude.csv    SupTable_8_BMAL1_phase.csv
#          Supplementary_Tables/SupTable_9_TMED10_amplitude.csv   SupTable_10_TMED10_phase.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)

RHYTHMIC <- rhythmicTissues()
tod   <- loadDIP()
plink <- loadGenotypes()$X

all <- bind_rows(lapply(seq_len(nrow(EXAMPLES)), function(k) {
  ex <- EXAMPLES[k, ]
  bind_rows(lapply(RHYTHMIC, function(tissue) {
    res <- tryCatch(singleTissueInteractionLRT(tod, plink, ex$gene, ex$variant, tissue),
                    error = function(e) { message(ex$symbol, ' ', tissue, ': ', conditionMessage(e)); NULL })
    if (is.null(res)) return(NULL)
    gp <- residualisedGenotypeParameters(res$data)
    if (nrow(gp) == 0) return(NULL)
    gp %>% mutate(gene = ex$gene, hgnc_symbol = ex$symbol, tissue = tissue, variant = ex$variant)
  }))
})) %>%
  mutate(amplitude = sqrt(sin_coef^2 + cos_coef^2),
         phase     = (atan2(sin_coef, cos_coef) %% (2 * pi)) / (2 * pi) * 24) %>%
  select(gene, hgnc_symbol, tissue, variant, genotype, n = n_g, amplitude, phase)

for (k in seq_len(nrow(EXAMPLES))) {
  ex <- EXAMPLES[k, ]
  d  <- filter(all, hgnc_symbol == ex$symbol)
  message(ex$symbol, ': ', n_distinct(d$tissue), ' tissues, ', nrow(d), ' tissue x genotype rows')
  fwrite(select(d, -phase),     paste0(TABLE_DIR, 'SupTable_', 5 + 2 * k, '_', ex$symbol, '_amplitude.csv'))
  fwrite(select(d, -amplitude), paste0(TABLE_DIR, 'SupTable_', 6 + 2 * k, '_', ex$symbol, '_phase.csv'))
}
