# =============================================================================
# 21  Supplementary Tables workbook
#
# Combines Supplementary Tables 1-14 (CSV) into one .xlsx: a table of contents
# plus one tab per table, each with its legend above the data.
#
# Run:     Rscript 21_supplementaryWorkbook.R    (after every table is built)
# Output:  Writing/Submission/Supplementary_Tables.xlsx
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(openxlsx)

OUT <- paste0(SUBMISSION, 'Supplementary_Tables.xlsx')

TABLES <- tibble::tribble(
  ~tab, ~file, ~legend,
  'S1_HarmonicGenes', 'SupTable_1_HarmonicGenes.csv',
  'Supplementary Table 1. Rhythmic genes identified across GTEx v10 tissues. One row per gene-tissue pair that passed the rhythmicity test (Bonferroni-adjusted likelihood ratio test p < 0.05; n = 4,033 pairs), comparing the harmonic (cosinor) model with the linear model. The Bonferroni correction is applied within tissue, over the genes tested in that tissue; 1,285,139 gene-tissue pairs were tested across 50 tissues. 49 of those 50 tissues contribute at least one rhythmic gene, Cells Cultured fibroblasts contributing none. Columns: term, fitted model; X.Df, model degrees of freedom; LogLik, log-likelihood; df, difference in degrees of freedom between models; statistic, chi-square LRT statistic; p.value, raw LRT p-value; gene, Ensembl gene ID; hgnc_symbol, gene symbol, falling back to the Ensembl gene ID where Ensembl assigns no symbol (resolved from the Ensembl REST lookup endpoint by Code/03_geneSymbolMap.R and cached in Data/ref/annot/ensembl_hgnc_map.csv); tissue, GTEx tissue; p.adj, Bonferroni-adjusted p-value, equal to p.value multiplied by the number of genes tested in that tissue.',
  'S2_cdEQTL', 'SupTable_2_cdEQTL.csv',
  'Supplementary Table 2. All 168 significant cd-eQTL. One row per gene-variant-tissue association passing the interaction-model likelihood ratio test (raw p < 1e-4). Columns: term, fitted interaction model; X.Df, model degrees of freedom; LogLik, log-likelihood; df, difference in degrees of freedom; statistic, chi-square LRT statistic; p.value, raw interaction LRT p-value; gene, Ensembl gene ID; variant, lead SNP; tissue, GTEx tissue; p.adj, adjusted p-value; geneModif, Ensembl gene ID without version suffix; hgnc_symbol, gene symbol, falling back to the Ensembl gene ID where Ensembl assigns no symbol.',
  'S3_cdEQTL_54independent', 'SupTable_3_cdEQTL_54independent.csv',
  'Supplementary Table 3. The 54 independent cd-eQTL obtained after LD pruning (pairwise r2 < 0.1; lead variant retained), associated with 51 genes across 28 tissues. Columns as in Table 2.',
  'S4_cdEQTL_coefficients', 'SupTable_4_cdEQTL_coefficients.csv',
  'Supplementary Table 4. Interaction-model regression coefficients for every significant cd-eQTL. One row per fitted coefficient (n = 4,004 rows across 167 associations), giving the full harmonic interaction model underlying the likelihood ratio tests reported in Table 2. Columns: gene, Ensembl gene ID; hgnc_symbol, gene symbol, falling back to the Ensembl gene ID where Ensembl assigns no symbol; variant, lead SNP; tissue, GTEx tissue; term, model term; estimate, fitted coefficient; std.error; statistic, t statistic; p.value. Terms appear in fitted order: intercept; the variant, named by its rsID and coded as the number of alternate allele copies; the technical covariates sex, platform, pcr, PC1-5 and InferredCov1-10; the cosinor terms sin(radian) and cos(radian); and the genotype-by-time interaction terms rsID:sin(radian) and rsID:cos(radian), whose joint significance defines the cd-eQTL. Genotype-group amplitude is sqrt((b_sin + g*b_variant:sin)^2 + (b_cos + g*b_variant:cos)^2) for g copies of the alternate allele. Four associations in single-sex tissues (Ovary, Testis, Vagina) have 23 rather than 24 coefficients, as sex is constant and dropped from the fit.',
  'S5_dryRvalidation', 'SupTable_5_dryRvalidation.csv',
  'Supplementary Table 5. dryR validation of the 54 independent cd-eQTL. Each locus was re-fitted with dryR on the same covariate-residualised expression used for the cosinor parameter tests, with genotype dosage as the condition and CHIRAL-inferred internal phase as time. dryR fits every candidate rhythm model (5 for two genotype groups, 15 for three) and selects one by Bayesian information criterion (BIC) weight, with no significance threshold. dryR selects a model in which rhythmicity is not shared between genotype groups at 51 of 54 loci (94.4%); at the remaining 3, ZFYVE28 (Brain Cortex), DBP (Nerve Tibial) and ENSG00000236266 (Ovary), it selects one rhythm shared by all genotype groups. The full BIC weight profile of every locus is shown in Supplementary Figure 1. Columns: gene, Ensembl gene ID; hgnc_symbol, gene symbol, falling back to the Ensembl gene ID where Ensembl assigns no symbol; variant, lead SNP; tissue, GTEx tissue; n, donors entering the dryR fit; n_genotype_groups and genotype_groups, the dosage groups retained at n >= 10 (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate); regression_LRT_statistic, regression_LRT_p.value and regression_LRT_p.adj, the chi-square statistic, raw p-value and adjusted p-value of the genotype-by-time interaction likelihood ratio test from the harmonic regression model, carried over from Table 3 and not produced by dryR; dryR_chosen_model, the BIC-selected model index, decoded against the number of groups since the same integer denotes different models for 2 and 3 groups; dryR_model_desc, that model in words, where groups written together share one rhythm and a group with no cosinor term is flat; dryR_model_BICW, the BIC weight of the selected model among the 5 or 15 candidates; concordance, discordant where the selected model is one rhythm shared by all genotype groups (no genotype effect on rhythmicity), concordant otherwise.',
  'S6_cdEQTL_classification', 'SupTable_6_cdEQTL_classification.csv',
  "Supplementary Table 6. Classification of the effect of each of the 54 independent cd-eQTL (Supplementary Table 3) on the amplitude and the acrophase of the rhythm. For every locus a separate cosinor was fitted to each genotype group on covariate-residualised expression, and amplitudes and acrophases were compared across groups with the F-tests of Bingham et al. (1982); genotype groups with fewer than 10 donors were not fitted. Adjusted p-values are Benjamini-Hochberg over the 54 loci, separately for each parameter. These tests describe which rhythmic parameter carries the genetic effect and do not re-assess the association itself, which was established by the interaction likelihood ratio test (Supplementary Table 2). A locus is classified as 'both' if both adjusted p-values are < 0.1, 'amplitude' or 'phase' if only that one is; when neither passes, the locus is assigned to the parameter with the smaller adjusted p-value. Because an acrophase is only defined for a rhythmic group, a phase call further requires at least two genotype groups whose amplitude differs from zero (Bingham zero-amplitude test); phase loci failing this requirement are classified as 'undetermined'. All loci classified as 'both' meet it. Counts: 28 amplitude (5 CCG, 23 CRG), 13 phase (0 CCG, 13 CRG), 10 both (3 CCG, 7 CRG), 3 undetermined (0 CCG, 3 CRG). Columns: gene, Ensembl gene ID; hgnc_symbol, gene symbol, falling back to the Ensembl gene ID where Ensembl assigns no symbol; variant, lead SNP; tissue, GTEx tissue; gene_type, CCG (core circadian gene) or CRG (circadian regulated gene); n_genotype_groups, number of genotype dosage groups fitted; group_sizes, donors per group, in dosage order 0/1/2; F_amplitude, p_amplitude and padj_amplitude, the F statistic, raw and adjusted p-value of the test of equal amplitude; F_acrophase, p_acrophase and padj_acrophase, the same for the test of equal acrophase; df_num and df_denom, the F-test degrees of freedom; classification, amplitude, phase, both or undetermined; assigned_by, FDR when at least one adjusted p-value is < 0.1, smaller adjusted p otherwise; zeroAmp_p_by_group, p-value of the Bingham test that the amplitude of each genotype group is zero, as dosage:p; n_groups_nonzero_amplitude, number of groups whose amplitude differs from zero (95% confidence region excluding the origin).",
  'S7_BMAL1_amplitude', 'SupTable_7_BMAL1_amplitude.csv',
  'Supplementary Table 7. Oscillation amplitude of BMAL1 (variant rs11022718) by tissue and genotype group, across the 49 tissues with at least one rhythmic gene. From a separate cosinor fitted to each genotype group on covariate-residualised expression, the model used for the amplitude and acrophase tests (Supplementary Table 6) and shown in Figure 2b and Figure 3a-b; genotype groups with fewer than 10 donors were not fitted. Columns: gene, Ensembl gene ID; hgnc_symbol, gene symbol; tissue; variant; genotype, number of copies of the rs11022718-T allele; n, donors in the genotype group; amplitude, sqrt(sin^2 + cos^2) of the fitted cosinor coefficients.',
  'S8_BMAL1_phase', 'SupTable_8_BMAL1_phase.csv',
  'Supplementary Table 8. Acrophase of BMAL1 (variant rs11022718) by tissue and genotype group, across the 49 tissues with at least one rhythmic gene. From a separate cosinor fitted to each genotype group on covariate-residualised expression, the model used for the amplitude and acrophase tests (Supplementary Table 6) and shown in Figure 2b and Figure 3a-b; genotype groups with fewer than 10 donors were not fitted. Columns: gene, Ensembl gene ID; hgnc_symbol, gene symbol; tissue; variant; genotype, number of copies of the rs11022718-T allele; n, donors in the genotype group; phase, time of peak expression in hours (0-24) on the CHIRAL-inferred internal time scale, atan2(sin, cos) of the fitted cosinor coefficients.',
  'S9_TMED10_amplitude', 'SupTable_9_TMED10_amplitude.csv',
  'Supplementary Table 9. Oscillation amplitude of TMED10 (variant rs76379942) by tissue and genotype group. Computed as in Table 7; the genotype column counts copies of the rs76379942 deletion allele (AG>A). Columns as in Table 7.',
  'S10_TMED10_phase', 'SupTable_10_TMED10_phase.csv',
  'Supplementary Table 10. Acrophase of TMED10 (variant rs76379942) by tissue and genotype group. Computed as in Table 8; the genotype column counts copies of the rs76379942 deletion allele (AG>A). Columns as in Table 8.',
  'S11_BMAL1_cosineSimilarity', 'SupTable_11_BMAL1_cosineSimilarity.csv',
  'Supplementary Table 11. Pairwise tissue-by-tissue cosine similarity of the rs11022718 effect on BMAL1, computed from the four estimated rhythmic parameters (sin, cos, and their genotype interaction terms). The first two columns give the gene identity (gene, Ensembl gene ID; hgnc_symbol, gene symbol), constant down the table, and the third column is the reference tissue; the remaining 49 columns form the 49 x 49 similarity matrix. A value of 1 indicates an identical effect, -1 an opposite effect, and 0 an orthogonal effect.',
  'S12_TMED10_cosineSimilarity', 'SupTable_12_TMED10_cosineSimilarity.csv',
  'Supplementary Table 12. Pairwise tissue-by-tissue cosine similarity of the rs76379942 effect on TMED10, computed as in Table 11. Columns as in Table 11: gene and hgnc_symbol, then the reference tissue, then the 49 x 49 similarity matrix.',
  'S13_PheWAS', 'SupTable_13_PheWAS.csv',
  'Supplementary Table 13. Significant cd-eQTL-phenotype associations from the Pan-UK Biobank PheWAS (European ancestry). One row per variant-phenotype pair passing p < 4.6e-5 (n = 310 associations spanning 100 distinct phenotypes, 34 variants and 32 genes). The first two columns give the cd-eQTL gene the variant was discovered for (gene, Ensembl gene ID; hgnc_symbol, gene symbol), joined on rsID; every row maps to exactly one gene, spanning 32 genes across the 34 variants. The remaining columns follow the Pan-UK Biobank GWAS manifest; key fields: variant, chromosome:position:reference:alternate; neglog10_pval_EUR, -log10 p-value in Europeans; description / category / phewasClass, phenotype annotation; rsID, dbSNP identifier.',
  'S14_ChenComparison', 'SupTable_14_ChenComparison.csv',
  'Supplementary Table 14. Joint genotype-by-time interaction model applied to 10 cd-eQTL associations reported by Chen et al. (ref. 26) in GTEx v10. Three nested models are fitted per association (linear, harmonic, harmonic interaction), giving three rows each (n = 30). None of the interaction tests reached significance (all p > 0.15). Columns: term, fitted model; X.Df; LogLik; df; statistic; p.value; gene, Ensembl gene ID; variant; tissue; gene_name, gene symbol.')

headerStyle <- createStyle(fgFill = '#1F4E78', fontColour = '#FFFFFF', textDecoration = 'bold', fontSize = 10)
legendStyle <- createStyle(wrapText = TRUE, valign = 'top', fontSize = 10)

wb <- createWorkbook()
addWorksheet(wb, 'TOC')
writeData(wb, 'TOC', 'Supplementary Tables - Table of Contents', startRow = 1)
addStyle(wb, 'TOC', createStyle(textDecoration = 'bold', fontSize = 14), rows = 1, cols = 1)

toc <- data.frame(Tab = TABLES$tab, Table = sub('\\..*', '', TABLES$legend),
                  Description = sub('^[^.]*\\. ', '', TABLES$legend), `Data rows` = NA_integer_,
                  check.names = FALSE)
for (k in seq_len(nrow(TABLES))) {
  d  <- fread(paste0(TABLE_DIR, TABLES$file[k]))
  sh <- TABLES$tab[k]
  toc$`Data rows`[k] <- nrow(d)
  addWorksheet(wb, sh)
  writeData(wb, sh, TABLES$legend[k], startRow = 1)
  mergeCells(wb, sh, cols = 1:min(ncol(d), 12), rows = 1)
  addStyle(wb, sh, legendStyle, rows = 1, cols = 1)
  setRowHeights(wb, sh, rows = 1, heights = 75)
  writeData(wb, sh, d, startRow = 3, headerStyle = headerStyle)
  freezePane(wb, sh, firstActiveRow = 4)
  setColWidths(wb, sh, cols = seq_len(ncol(d)),
               widths = pmin(pmax(nchar(names(d)) + 2, 9,
                                  vapply(head(d, 60), function(x) max(nchar(as.character(x)), 0,
                                                                      na.rm = TRUE), numeric(1)) + 2), 40))
}
writeData(wb, 'TOC', toc, startRow = 3, headerStyle = headerStyle)
for (k in seq_len(nrow(TABLES)))
  writeFormula(wb, 'TOC', startRow = 3 + k, startCol = 1,
               x = makeHyperlinkString(sheet = TABLES$tab[k], row = 1, col = 1, text = TABLES$tab[k]))
addStyle(wb, 'TOC', createStyle(wrapText = TRUE, valign = 'top'),
         rows = 4:(3 + nrow(TABLES)), cols = 3, gridExpand = TRUE)
setColWidths(wb, 'TOC', cols = 1:4, widths = c(28, 34, 90, 11))
freezePane(wb, 'TOC', firstActiveRow = 4)

saveWorkbook(wb, OUT, overwrite = TRUE)
message('Wrote ', OUT)
