# =============================================================================
# 08  Supplementary Table 6: amplitude / acrophase classification of the 54
#     independent cd-eQTL (one row per locus)
#
# Run:     Rscript 08_supTable6_classification.R     (after 07)
# Output:  Supplementary_Tables/SupTable_6_cdEQTL_classification.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

symbols <- readIndependentCdeQTL()
zero    <- fread(paste0(RHYTHM_DIR, 'binghamZeroAmplitude_perGroup.csv'))
cls     <- fread(paste0(RHYTHM_DIR, 'binghamGatedClassification.csv'))
stopifnot(nrow(cls) == 54, nrow(semi_join(cls, symbols, by = c('gene', 'variant', 'tissue'))) == 54)

zeroSummary <- zero %>%
  arrange(genotype) %>%
  group_by(gene, variant, tissue) %>%
  summarise(n_z = n(),
            zeroAmp_p_by_group = paste(sprintf('%d:%s', genotype,
                                               formatC(p_zeroAmp, digits = 3, format = 'g')),
                                       collapse = '/'),
            n_groups_nonzero_amplitude = sum(!covers_origin),
            .groups = 'drop')

out <- cls %>%
  left_join(symbols, by = c('gene', 'variant', 'tissue')) %>%
  left_join(zeroSummary, by = c('gene', 'variant', 'tissue')) %>%
  mutate(classOrder = match(class_gated, c('amplitude', 'phase', 'both', 'undetermined')),
         sortP = ifelse(class_gated == 'amplitude', padj_amplitude_bingham, padj_acrophase_bingham)) %>%
  arrange(classOrder, sortP)
stopifnot(all(out$n_z == out$n_groups),
          all(out$n_groups_nonzero_amplitude == out$n_rhythmic),
          !any(out$class_base == 'both' & out$n_groups_nonzero_amplitude < 2))

out <- out %>%
  transmute(gene, hgnc_symbol, variant, tissue,
            gene_type = ifelse(gene %in% CCG, 'CCG', 'CRG'),
            n_genotype_groups = n_groups, group_sizes,
            F_amplitude, p_amplitude = p_amplitude_bingham, padj_amplitude = padj_amplitude_bingham,
            F_acrophase, p_acrophase = p_acrophase_bingham, padj_acrophase = padj_acrophase_bingham,
            df_num = k - 1, df_denom,
            classification = class_gated,
            assigned_by = ifelse(pmin(padj_amplitude, padj_acrophase) < CLASS_FDR,
                                 'FDR', 'smaller adjusted p'),
            zeroAmp_p_by_group, n_groups_nonzero_amplitude)
fwrite(out, paste0(TABLE_DIR, 'SupTable_6_cdEQTL_classification.csv'))
print(table(out$classification, out$gene_type))
