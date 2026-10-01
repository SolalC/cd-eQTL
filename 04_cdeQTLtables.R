# =============================================================================
# 04  Rhythmic genes and cd-eQTL tables, Figure 1b
#
#   Rhythmic genes: Bonferroni-adjusted harmonic LRT p < 0.05 within tissue.
#   cd-eQTL: interaction (model 3 vs 2) LRT p < 1e-4.
#   Independent cd-eQTL: the lead (smallest p) variant per gene x tissue.
#
# Run:     sbatch 04_cdeQTLtables.sh          (after 02 and 03)
# Outputs: Figures/Figure_1b_RhythmicGenes.jpeg                  Figure 1b
#          Supplementary_Tables/SupTable_1_HarmonicGenes.csv      Supplementary Table 1
#          Supplementary_Tables/SupTable_2_cdEQTL.csv             Supplementary Table 2
#          Supplementary_Tables/SupTable_3_cdEQTL_54independent.csv  Supplementary Table 3
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(MetBrewer)

# --- Rhythmic genes (Supplementary Table 1) ------------------------------------
harmonic <- bind_rows(lapply(harmonicFiles(), function(file) {
  fread(file) %>%
    mutate(tissue = str_remove(basename(file), '_LRT_harmonic_linear.csv'),
           p.adj  = p.adjust(p.value, 'bonferroni'))
}))
message('Gene x tissue pairs tested: ', nrow(harmonic))
harmonic.f <- harmonic %>% filter(p.adj < 0.05)
message('Rhythmic gene x tissue pairs: ', nrow(harmonic.f),
        ' in ', n_distinct(harmonic.f$tissue), ' tissues')

fwrite(harmonic.f %>% addSymbol() %>%
         select(term, X.Df, LogLik, df, statistic, p.value, gene, hgnc_symbol, tissue, p.adj),
       paste0(TABLE_DIR, 'SupTable_1_HarmonicGenes.csv'))

# --- Figure 1b: sample size and number of rhythmic genes per tissue -------------
sampleSize <- fread(SAMPLE_SIZE)
harmonic.ff <- harmonic.f %>%
  mutate(tissue = str_replace_all(tissue, '_', ' ')) %>%
  group_by(tissue) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  left_join(., sampleSize, by = 'tissue') %>%
  mutate(tissue = fct_reorder(tissue, n))

rhythmicGenesPlot <- ggplot(harmonic.ff, aes(x = n, y = tissue, col = tissue)) +
  geom_point() +
  geom_point(aes(x = -RNAandGenotype)) +
  geom_vline(xintercept = 100, lty = 2) +
  geom_vline(xintercept = 0,   lty = 3) +
  scale_color_manual(values = met.brewer('Renoir', 50)) +
  labs(x = 'Sample size ||| Number of rhythmic genes', y = '') +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(rhythmicGenesPlot, filename = paste0(FIG_DIR, 'Figure_1b_RhythmicGenes.jpeg'),
       width = 6, height = 6, dpi = 300)

# --- cd-eQTL (Supplementary Tables 2 and 3) ----------------------------------------
interactionFull <- bind_rows(lapply(list.files(paste0(QTL_DIR, '4_LRT'), full.names = TRUE),
  function(file) {
    fread(file) %>%
      mutate(tissue = str_remove(basename(file), '_LRT.csv')) %>%
      filter(str_detect(term, ':')) %>%
      mutate(p.adj = p.adjust(p.value, method = 'BH'))
  }))

interaction <- interactionFull %>% filter(p.value < CDEQTL_P) %>% addSymbol()
fwrite(interaction, paste0(TABLE_DIR, 'SupTable_2_cdEQTL.csv'))

interaction.f <- interaction %>%
  group_by(tissue, gene) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup()
fwrite(interaction.f, paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv'))

message('cd-eQTL rows: ', nrow(interaction),
        ' (unique gene/variant/tissue: ', n_distinct(interaction$gene, interaction$variant, interaction$tissue), ')')
message('Independent cd-eQTL: ', nrow(interaction.f), ', ', n_distinct(interaction.f$gene),
        ' genes, ', n_distinct(interaction.f$tissue), ' tissues')
