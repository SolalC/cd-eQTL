# =============================================================================
# 19  Pan-UK Biobank PheWAS of the independent cd-eQTL (European ancestry)
#
# Input: the Pan-UKB associations of the cd-eQTL lead variants across 1,085
# phenotypes (Results/phewas/inferred/phewasClean/*.tsv), the phenotype
# manifest with its PheWAS class (Pan-uk-phenotype.csv) and the Pan-UKB
# variant -> rsID key (Neale_rsID_key.csv). Urine assays are excluded.
# Significance: p < 0.05 / 1,085 = 4.6e-5.
#
# Run:     Rscript 19_figure5_phewas.R    (after 04; seconds)
# Outputs: Figures/Figure_5_PheWAS.{jpeg,pdf}           Figure 5 (rs11022718, BMAL1)
#          Supplementary_Tables/SupTable_13_PheWAS.csv   Supplementary Table 13
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(ggrepel)

FIGURE_VARIANT <- 'rs11022718'
SEED <- 1

# --- Inputs -----------------------------------------------------------------------
pheDf <- rbindlist(lapply(list.files(paste0(PHEWAS_DIR, 'phewasClean'), full.names = TRUE),
                          function(i) fread(i, select = c('variant', 'file', 'neglog10_pval_EUR'))))
variantKey <- fread(paste0(PHEWAS_DIR, 'Neale_rsID_key.csv'))
pheno      <- fread(paste0(PHEWAS_DIR, 'Pan-uk-phenotype.csv'))
cdeqtl     <- readIndependentCdeQTL()

pheDf.f <- pheDf %>%
  left_join(pheno,      by = 'file') %>%
  left_join(variantKey, by = 'variant') %>%
  arrange(phewasClass) %>%
  filter(phewasClass != 'Urine assays')

phewasClassOrder <- c(
  'Blood biochemistry', 'Blood count', 'Physical measures', 'Lifestyle and environment',
  'Congenital malformations, deformations and chromosomal abnormalities',
  'Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism',
  'Diseases of the circulatory system', 'Diseases of the digestive system',
  'Diseases of the ear and mastoid process', 'Diseases of the eye and adnexa',
  'Diseases of the genitourinary system', 'Diseases of the musculoskeletal system and connective tissue',
  'Diseases of the nervous system', 'Diseases of the respiratory system',
  'Diseases of the skin and subcutaneous tissue', 'Endocrine, nutritional and metabolic diseases',
  'Infectious and parasitic diseases',
  'Injury, poisoning and certain other consequences of external causes',
  'Mental and behavioural disorders', 'Neoplasms', 'Pregnancy, childbirth and the puerperium')
pheDf.f$phewasClass <- factor(pheDf.f$phewasClass, levels = phewasClassOrder)

threshold_phewas <- -log10(PHEWAS_P)

# --- Supplementary Table 13: significant associations of all cd-eQTL --------------
geneKey <- cdeqtl %>% distinct(rsID = variant, gene, hgnc_symbol)
supTable13 <- pheDf.f %>%
  filter(rsID %in% cdeqtl$variant, neglog10_pval_EUR >= threshold_phewas) %>%
  arrange(rsID, phewasClass, desc(neglog10_pval_EUR)) %>%
  inner_join(geneKey, by = 'rsID') %>%
  select(gene, hgnc_symbol, everything())
fwrite(supTable13, paste0(TABLE_DIR, 'SupTable_13_PheWAS.csv'))
message('cd-eQTL in the PheWAS: ', n_distinct(pheDf.f$rsID[pheDf.f$rsID %in% cdeqtl$variant]),
        '; significant: ', nrow(supTable13), ' associations, ',
        n_distinct(supTable13$description), ' phenotypes, ',
        n_distinct(supTable13$rsID), ' variants')

# --- Figure 5 -------------------------------------------------------------------------
d   <- pheDf.f %>% filter(rsID == FIGURE_VARIANT)
sig <- d %>% filter(neglog10_pval_EUR > threshold_phewas)
info <- cdeqtl %>% filter(variant == FIGURE_VARIANT)
titleTxt <- paste0('Type: ', ifelse(info$gene[1] %in% CCG, 'Core circadian Gene', 'Circadian regulated Gene'),
                   '\nGene: ', info$hgnc_symbol[1],
                   '\nSignificant: ', paste(str_replace_all(unique(info$tissue), '_', ' '), collapse = ', '),
                   '\nVariant: ', FIGURE_VARIANT)

p <- ggplot(d, aes(x = phewasClass, y = neglog10_pval_EUR, col = phewasClass)) +
  geom_point(position = position_jitter(seed = SEED), size = 1.5) +
  geom_hline(yintercept = threshold_phewas, lty = 2, col = 'darkred',  linewidth = 0.3) +
  geom_hline(yintercept = -log10(5e-4),     lty = 2, col = 'darkblue', linewidth = 0.3) +
  geom_text_repel(data = sig, aes(label = description),
                  size = 3, max.overlaps = Inf, box.padding = 0.3,
                  min.segment.length = 0, segment.size = 0.2,
                  position = position_jitter(seed = SEED), seed = SEED) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 28), drop = FALSE) +
  labs(x = '', y = '-log10(p-value)', title = titleTxt) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, lineheight = 1.2),
        legend.position = 'none',
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, lineheight = 0.85))

for (ext in c('jpeg', 'pdf'))
  ggsave(paste0(FIG_DIR, 'Figure_5_PheWAS.', ext), p, width = 14, height = 11, units = 'in', dpi = 300)
message(FIGURE_VARIANT, ': ', nrow(d), ' phenotypes tested, ', nrow(sig), ' above the threshold')
