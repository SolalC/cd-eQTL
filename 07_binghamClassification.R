# =============================================================================
# 07  Amplitude / acrophase classification of the 54 independent cd-eQTL
#
# For each locus, expression is residualised on the covariates (once, all
# samples) and an independent cosinor is fitted to each genotype group with
# >= 10 donors. Then:
#   1. Bingham et al. (1982) F-tests for equal amplitude (eq. 49) and equal
#      acrophase (eq. 50) across groups; BH-adjusted over the 54 loci.
#   2. Per-group zero-amplitude test (eq. 34, pooled variance): the acrophase
#      of a group is defined only if its 95% confidence ellipse excludes 0.
#   3. Classification: 'both' if both adjusted p < 0.1, else the parameter that
#      passes; if neither passes, the one with the smaller adjusted p. A phase
#      call then needs >= 2 groups of non-zero amplitude, otherwise
#      'undetermined' (a 'both' locus failing it would become 'amplitude').
#
# Run:     sbatch 07_binghamClassification.sh   (after 04)
# Outputs: Results/inferred/rhythmicityParameters/binghamGroupTest.csv
#          Results/inferred/rhythmicityParameters/binghamZeroAmplitude_perGroup.csv
#          Results/inferred/rhythmicityParameters/binghamGatedClassification.csv
# Used by: Supplementary Table 6 (08)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)

# Zero-amplitude test b' Sigma^-1 b / 2 ~ F(2, nu); the ellipse covers the
# origin exactly when this is not significant at 1 - conf.
zeroAmplitude <- function(b, Sigma, nu, conf = 0.95) {
  q   <- as.numeric(t(b) %*% solve(Sigma) %*% b)
  amp <- sqrt(sum(b^2))
  gA  <- b / amp
  gP  <- c(b[2], -b[1]) / amp^2
  c(p_zeroAmp = pf(q / 2, 2, nu, lower.tail = FALSE),
    covers_origin = q <= 2 * qf(conf, 2, nu),
    amp = amp, amp_se = sqrt(as.numeric(t(gA) %*% Sigma %*% gA)),
    phase_h = (atan2(-b[2], b[1]) %% (2 * pi)) / (2 * pi) * 24,
    phase_se_h = sqrt(as.numeric(t(gP) %*% Sigma %*% gP)) / (2 * pi) * 24)
}

tod      <- loadDIP()
plinkAll <- loadGenotypes()
plink    <- plinkAll$X
bim      <- plinkAll$bim
cdeqtl   <- readIndependentCdeQTL()
message('Classifying ', nrow(cdeqtl), ' independent cd-eQTL')

groupList <- vector('list', nrow(cdeqtl))
perGroup  <- vector('list', nrow(cdeqtl))
for (i in seq_len(nrow(cdeqtl))) {
  row <- cdeqtl[i, ]
  message('Running: ', row$tissue, ' ', row$gene, ' ', row$variant)
  res <- tryCatch(singleTissueInteractionLRT(tod, plink, row$gene, row$variant, row$tissue),
                  error = function(e) { message('row ', i, ': ', conditionMessage(e)); NULL })
  if (is.null(res) || is.null(res$fit)) next

  # 1. Group F-tests
  groupList[[i]] <- binghamGroupTest(res$data) %>%
    mutate(gene = row$gene, variant = row$variant, tissue = row$tissue)

  # 2. Zero-amplitude test per group, with the pooled variance of the F-tests
  gp <- residualisedGenotypeParameters(res$data)
  if (nrow(gp) == 0) next
  nu <- sum(gp$df)
  s2 <- sum(gp$rss) / nu
  perGroup[[i]] <- bind_rows(lapply(seq_len(nrow(gp)), function(j) {
    b     <- c(gp$sin_coef[j], gp$cos_coef[j])
    Sigma <- matrix(c(gp$cov_ss[j], gp$cov_sc[j], gp$cov_sc[j], gp$cov_cc[j]), 2, 2) *
             s2 / (gp$rss[j] / gp$df[j])
    data.frame(gene = row$gene, variant = row$variant, tissue = row$tissue,
               symbol = row$hgnc_symbol, genotype = gp$genotype[j], n_g = gp$n_g[j],
               t(zeroAmplitude(b, Sigma, nu)))
  }))
}

groupDf <- bind_rows(groupList) %>%
  mutate(padj_amplitude_bingham = p.adjust(p_amplitude_bingham, method = 'BH'),
         padj_acrophase_bingham = p.adjust(p_acrophase_bingham, method = 'BH'))
fwrite(groupDf, paste0(RHYTHM_DIR, 'binghamGroupTest.csv'))

perGroup <- bind_rows(perGroup) %>% mutate(covers_origin = as.logical(covers_origin))
fwrite(perGroup, paste0(RHYTHM_DIR, 'binghamZeroAmplitude_perGroup.csv'))

# --- 3. Classification --------------------------------------------------------
# Largest acrophase difference between any two rhythmic groups, wrapped to +-12 h.
maxShift <- function(ph) {
  if (length(ph) < 2) return(NA_real_)
  d <- outer(ph, ph, '-'); d <- (d + 12) %% 24 - 12
  max(abs(d))
}
gate <- perGroup %>%
  group_by(gene, variant, tissue) %>%
  summarise(symbol = first(symbol), n_groups = n(),
            n_rhythmic = sum(!covers_origin),
            rhythmic_groups = paste(genotype[!covers_origin], collapse = '/'),
            group_sizes = paste(n_g, collapse = '/'),
            max_shift_rhythmic_h = maxShift(phase_h[!covers_origin]),
            .groups = 'drop')

chrom <- setNames(as.character(bim$chr), bim$id)
classDf <- groupDf %>%
  left_join(gate, by = c('gene', 'variant', 'tissue')) %>%
  mutate(chr = chrom[variant],
         class_base = case_when(
           padj_amplitude_bingham < CLASS_FDR & padj_acrophase_bingham < CLASS_FDR ~ 'both',
           padj_amplitude_bingham < CLASS_FDR                                       ~ 'amplitude',
           padj_acrophase_bingham < CLASS_FDR                                       ~ 'phase',
           padj_amplitude_bingham <= padj_acrophase_bingham                         ~ 'amplitude',
           TRUE                                                                     ~ 'phase'),
         class_gated = case_when(
           class_base == 'phase' & n_rhythmic < 2 ~ 'undetermined',
           class_base == 'both'  & n_rhythmic < 2 ~ 'amplitude',
           TRUE                                   ~ class_base)) %>%
  arrange(class_gated, padj_acrophase_bingham)
fwrite(classDf, paste0(RHYTHM_DIR, 'binghamGatedClassification.csv'))

print(table(base = classDf$class_base, gated = classDf$class_gated))
