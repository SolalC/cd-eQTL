# =============================================================================
# 09  Per-genotype amplitude and acrophase: interaction model vs residualised
#     per-genotype cosinor, for the 54 independent cd-eQTL
#
#   (1) interaction model: b(g) = b + g * b_interaction (additive dosage)
#   (2) residualised cosinor: covariate-residualise once, then an independent
#       cosinor per genotype group (the fits compared by the Bingham tests, 07)
# Acrophase is expressed as the peak time of day (hours) for both.
#
# Run:     sbatch 09_methodComparison.sh   (after 04)
# Output:  Results/inferred/rhythmicityParameters/binghamMethodComparison.csv
# Used by: Supplementary Figure 2 (10)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)

# (2) residualised per-genotype cosinor
groupGenotypeParameters <- function(data) {
  gp <- residualisedGenotypeParameters(data)
  if (nrow(gp) == 0) return(NULL)
  gp %>% transmute(genotype, n_g,
                   amp_group   = sqrt(sin_coef^2 + cos_coef^2),
                   phase_group = atan2(-sin_coef, cos_coef),                        # Bingham convention
                   peakH_group = (atan2(sin_coef, cos_coef) %% (2 * pi)) / (2 * pi) * 24)
}

# (1) interaction model
interactionGenotypeParameters <- function(fit) {
  binghamGenotypeParameters(fit) %>% transmute(
    genotype,
    amp_interaction   = amplitude,
    phase_interaction = phase,
    peakH_interaction = (atan2(sin_coef, cos_coef) %% (2 * pi)) / (2 * pi) * 24)
}

# Smallest signed difference between two times of day on a 24 h circle.
circDiffH <- function(a, b) ((a - b + 12) %% 24) - 12

tod    <- loadDIP()
plink  <- loadGenotypes()$X
cdeqtl <- readIndependentCdeQTL()

cmpList <- vector('list', nrow(cdeqtl))
for (i in seq_len(nrow(cdeqtl))) {
  row <- cdeqtl[i, ]
  message('Running: ', row$tissue, ' ', row$gene, ' ', row$variant)
  res <- tryCatch(singleTissueInteractionLRT(tod, plink, row$gene, row$variant, row$tissue),
                  error = function(e) { message('row ', i, ': ', conditionMessage(e)); NULL })
  if (is.null(res) || is.null(res$fit)) next
  grp <- groupGenotypeParameters(res$data)
  if (is.null(grp)) next
  cmpList[[i]] <- inner_join(grp, interactionGenotypeParameters(res$fit), by = 'genotype') %>%
    mutate(gene = row$gene, variant = row$variant, tissue = row$tissue,
           hgnc_symbol   = row$hgnc_symbol,
           amp_absdiff   = amp_interaction - amp_group,
           amp_reldiff   = (amp_interaction - amp_group) / amp_group,
           peakH_absdiff = circDiffH(peakH_interaction, peakH_group))
}

cmp <- bind_rows(cmpList) %>%
  select(gene, hgnc_symbol, variant, tissue, genotype, n_g,
         amp_interaction, amp_group, amp_absdiff, amp_reldiff,
         peakH_interaction, peakH_group, peakH_absdiff,
         phase_interaction, phase_group)
fwrite(cmp, paste0(RHYTHM_DIR, 'binghamMethodComparison.csv'))
message('Rows (cd-eQTL x genotype): ', nrow(cmp))
