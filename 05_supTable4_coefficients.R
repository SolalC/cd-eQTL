# =============================================================================
# 05  Supplementary Table 4: interaction-model coefficients of every cd-eQTL
#
# One row per fitted coefficient of the harmonic interaction model (model 3 of
# 02) for each association of Supplementary Table 2, in fitted order. Rows
# duplicated upstream are kept once.
#
# Run:     Rscript 05_supTable4_coefficients.R     (after 04)
# Output:  Supplementary_Tables/SupTable_4_cdEQTL_coefficients.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

sig <- fread(paste0(TABLE_DIR, 'SupTable_2_cdEQTL.csv')) %>%
  distinct(gene, variant, tissue, .keep_all = TRUE) %>%
  select(gene, hgnc_symbol, variant, tissue) %>%
  mutate(order = row_number())

coef <- bind_rows(lapply(
  list.files(paste0(QTL_DIR, '3_harmonicInteraction'), pattern = '_harmonicInteraction.csv$',
             full.names = TRUE),
  function(f) {
    fread(f, colClasses = 'character') %>%
      filter(term != 'term', gene != '') %>%   # repeated headers from chunked output
      mutate(tissue = str_remove(basename(f), '_harmonicInteraction.csv')) %>%
      semi_join(sig, by = c('gene', 'variant', 'tissue'))
  })) %>%
  distinct(gene, variant, tissue, term, .keep_all = TRUE)   # first block of each term

out <- sig %>%
  inner_join(coef, by = c('gene', 'variant', 'tissue'), relationship = 'one-to-many') %>%
  arrange(order) %>%
  mutate(tissue = str_replace_all(tissue, '_', ' ')) %>%
  select(gene, hgnc_symbol, variant, tissue, term, estimate, std.error, statistic, p.value)
fwrite(out, paste0(TABLE_DIR, 'SupTable_4_cdEQTL_coefficients.csv'))
message('Wrote ', nrow(out), ' rows for ', nrow(sig), ' cd-eQTL')
