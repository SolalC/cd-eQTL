# =============================================================================
# 13  Supplementary Table 5: dryR validation of the 54 independent cd-eQTL
#
# One row per locus: the dryR selected model (from 11) and, for reference, the
# interaction LRT of the regression model (from Supplementary Table 3, prefixed
# regression_LRT_). 'discordant' where dryR selects the model in which every
# genotype group shares one rhythm (no genotype effect on rhythmicity).
#
# Run:     Rscript 13_supTable5_dryR.R     (after 11)
# Output:  Supplementary_Tables/SupTable_5_dryRvalidation.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

sharedModel <- c(`2` = 4L, `3` = 11L)
dry <- fread(paste0(DRYR_DIR, 'dryR_perLocus.csv'))
lrt <- fread(paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv')) %>%
  select(gene, variant, tissue,
         regression_LRT_statistic = statistic,
         regression_LRT_p.value   = p.value,
         regression_LRT_p.adj     = p.adj)
stopifnot(nrow(dry) == 54, nrow(lrt) == 54)

out <- dry %>%
  left_join(lrt, by = c('gene', 'variant', 'tissue')) %>%
  mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == '', gene, hgnc_symbol),
         concordance = ifelse(dryR_chosen_model == sharedModel[as.character(n_groups)],
                              'discordant', 'concordant')) %>%
  select(gene, hgnc_symbol, variant, tissue, n,
         n_genotype_groups = n_groups, genotype_groups,
         regression_LRT_statistic, regression_LRT_p.value, regression_LRT_p.adj,
         dryR_chosen_model, dryR_model_desc, dryR_model_BICW, concordance) %>%
  arrange(tissue, hgnc_symbol)
fwrite(out, paste0(TABLE_DIR, 'SupTable_5_dryRvalidation.csv'))
print(table(out$concordance))
