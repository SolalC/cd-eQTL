# =============================================================================
# 10  Supplementary Figure 2: agreement of the per-genotype amplitude (a) and
#     acrophase (b) between the interaction model (x) and the residualised
#     per-genotype cosinor (y), with the identity line and Pearson r
#
# Run:     Rscript 10_suppFigure2_methodComparison.R    (after 09; seconds)
# Outputs: Supplementary_Figures/SuppFigure_2a_EffectSizeComparison_amplitude.pdf
#          Supplementary_Figures/SuppFigure_2b_EffectSizeComparison_acrophase.pdf
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

cmp <- fread(paste0(RHYTHM_DIR, 'binghamMethodComparison.csv')) %>%
  mutate(genotype = factor(genotype, levels = c(0, 1, 2)))

genoCols <- c('0' = '#1b9e77', '1' = '#d95f02', '2' = '#7570b3')

rAmp   <- cor(cmp$amp_interaction,   cmp$amp_group)
rPhase <- cor(cmp$peakH_interaction, cmp$peakH_group)
nPts   <- nrow(cmp)
nQTL   <- dplyr::n_distinct(paste(cmp$gene, cmp$variant, cmp$tissue))

baseTheme <- theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom',
        plot.title = element_text(face = 'bold'),
        aspect.ratio = 1)

# --- Panel A: amplitude ------------------------------------------------------
ampMax <- max(cmp$amp_interaction, cmp$amp_group) * 1.02
pAmp <- ggplot(cmp, aes(amp_interaction, amp_group, colour = genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey50') +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_colour_manual(values = genoCols) +
  coord_equal(xlim = c(0, ampMax), ylim = c(0, ampMax)) +
  labs(title = 'A  Amplitude',
       x = 'Interaction model (per-allele cosinor)',
       y = 'Bingham residualised (per-genotype cosinor)',
       colour = 'genotype') +
  annotate('text', x = 0, y = ampMax, hjust = 0, vjust = 1,
           label = sprintf('r = %.2f\n%d groups, %d cd-eQTL', rAmp, nPts, nQTL)) +
  baseTheme

# --- Panel B: acrophase (peak time of day, hours) ----------------------------
pPhase <- ggplot(cmp, aes(peakH_interaction, peakH_group, colour = genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey50') +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_colour_manual(values = genoCols) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  coord_equal() +
  labs(title = 'B  Acrophase (peak time of day, h)',
       x = 'Interaction model (per-allele cosinor)',
       y = 'Bingham residualised (per-genotype cosinor)',
       colour = 'genotype') +
  annotate('text', x = 0, y = 24, hjust = 0, vjust = 1,
           label = sprintf('r = %.2f', rPhase)) +
  baseTheme

ggsave(paste0(SUPFIG_DIR, 'SuppFigure_2a_EffectSizeComparison_amplitude.pdf'),
       pAmp, width = 6, height = 6.4)
ggsave(paste0(SUPFIG_DIR, 'SuppFigure_2b_EffectSizeComparison_acrophase.pdf'),
       pPhase, width = 6, height = 6.4)
message('Amplitude r = ', round(rAmp, 3), '   Acrophase r = ', round(rPhase, 3))
