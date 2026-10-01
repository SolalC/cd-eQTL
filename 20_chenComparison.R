# =============================================================================
# 20  Comparison with Chen et al. (2025)
#
#   Supplementary Table 14  the joint interaction model (three nested models,
#                           as in 02) applied to the associations named by Chen
#                           et al., in GTEx v10
#   Supplementary Figure 4  AIF1 / rs7740525 / heart left ventricle, expression
#                           against DIP by genotype. The per-panel p is a
#                           rhythmicity LRT within the genotype group
#                           (intercept vs cosinor, no covariates); the title p
#                           is the joint genotype-by-time interaction LRT with
#                           all covariates.
#
# Run:     sbatch 20_chenComparison.sh
# Outputs: Supplementary_Tables/SupTable_14_ChenComparison.csv
#          Supplementary_Figures/SuppFigure_4_AIF1_ChenComparison.{jpeg,pdf}
#          Results/published/publishedReplication_harmonicInteraction.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)

tod   <- loadDIP()
plink <- loadGenotypes()

# --- Supplementary Table 14 --------------------------------------------------------
published <- fread(CHEN_ASSOC)
colnames(published) <- c('SNP', 'Gene', 'Tissue', 'Effect', 'Amplitude', 'pvalue', 'N')

tissueMap <- c(
  'Heart - Left Ventricle'            = 'Heart_Left_Ventricle',
  'Brain - Cerebellar Hemisphere'     = 'Brain_Cerebellar_Hemisphere',
  'Brain - Putamen basal ganglia'     = 'Brain_Putamen_basal_ganglia',
  'Adipose - Visceral'                = 'Adipose_Visceral_Omentum',
  'Muscle - Skeletal (GTEx)'          = 'Muscle_Skeletal',
  'Brain - Anterior cingulate cortex' = 'Brain_Anterior_cingulate_cortex_BA24',
  'Brain - Caudate'                   = 'Brain_Caudate_basal_ganglia',
  'Brain - Frontal Cortex'            = 'Brain_Frontal_Cortex_BA9',
  'Brain - Hippocampus'               = 'Brain_Hippocampus',
  'Brain - Hypothalamus'              = 'Brain_Hypothalamus')
geneMap <- c(
  AIF1   = 'ENSG00000204472', CIART = 'ENSG00000159208', MAP2K3 = 'ENSG00000034152',
  CRY1   = 'ENSG00000008405', TFRC  = 'ENSG00000072274', DISC1  = 'ENSG00000271617',
  NR1D1  = 'ENSG00000126368')

triplets <- published %>%
  select(SNP, Gene, Tissue) %>%
  distinct() %>%
  mutate(tissueName = tissueMap[Tissue], ensembl_gene_id = geneMap[Gene]) %>%
  filter(!is.na(tissueName), !is.na(ensembl_gene_id))

# Variants absent from the genotype subset cannot be tested
missing <- setdiff(unique(triplets$SNP), rownames(plink$X))
if (length(missing)) message('Skipping SNPs not in plink file: ', paste(missing, collapse = ', '))
triplets <- triplets %>% filter(SNP %in% rownames(plink$X))

lrtAll <- data.frame(); coefAll <- data.frame()
for (i in seq_len(nrow(triplets))) {
  t <- triplets[i, ]
  message(sprintf('[%d/%d] %s / %s / %s', i, nrow(triplets), t$Gene, t$SNP, t$tissueName))
  result <- tryCatch(singleTissueInteractionLRT(tod, plink$X, t$ensembl_gene_id, t$SNP, t$tissueName),
                     error = function(e) { message('  Error: ', conditionMessage(e)); NULL })
  if (is.null(result)) next
  lrtAll  <- rbind(lrtAll,  result$lrt  %>% mutate(gene_name = t$Gene))
  coefAll <- rbind(coefAll, result$coef %>% mutate(gene_name = t$Gene, tissue = t$tissueName))
}
fwrite(lrtAll,  paste0(TABLE_DIR, 'SupTable_14_ChenComparison.csv'))
fwrite(coefAll, paste0(CHEN_DIR, 'publishedReplication_harmonicInteraction.csv'))
message('Smallest interaction p: ', signif(min(lrtAll$p.value[str_detect(lrtAll$term, ':')]), 3))

# --- Supplementary Figure 4: AIF1 -----------------------------------------------------
GENE <- 'ENSG00000204472'; VARIANT <- 'rs7740525'; TISSUE <- 'Heart_Left_Ventricle'

indv <- selectIndividual(tod, TISSUE)
variantDf <- data.frame(SUBJID   = paste0('GTEX-', unlist(map(str_split(colnames(plink$X), '-'), 2))),
                        genotype = as.integer(plink$X[VARIANT, ]))
geneCol  <- grep(GENE, colnames(indv$tissue), value = TRUE)
geneExpr <- indv$tissue[, c('SUBJID', geneCol), with = FALSE]
setnames(geneExpr, geneCol, 'expression')

plotData <- indv$cov %>%
  left_join(variantDf, by = 'SUBJID') %>%
  left_join(geneExpr,  by = 'SUBJID') %>%
  filter(!is.na(genotype), !is.na(expression)) %>%
  mutate(hours = radian / (2 * pi) * 24)
nCounts <- plotData %>% count(genotype) %>%
  mutate(facetLabel = paste0('Genotype ', genotype, ' (n=', n, ')'))
plotData <- plotData %>%
  left_join(nCounts %>% select(genotype, facetLabel), by = 'genotype') %>%
  mutate(facetLabel = factor(facetLabel, levels = nCounts$facetLabel))

# Rhythmicity within each genotype group: intercept vs cosinor
gtPvals <- plotData %>%
  group_by(genotype, facetLabel) %>%
  group_map(function(df, key) {
    pval <- tryCatch(
      broom::tidy(lmtest::lrtest(lm(expression ~ 1, data = df),
                                 lm(expression ~ sin(radian) + cos(radian), data = df)))$p.value[2],
      error = function(e) NA_real_)
    data.frame(genotype = key$genotype, facetLabel = key$facetLabel, n = nrow(df),
               p_rhythmicity = pval,
               pLabel = ifelse(is.na(pval), 'p = NA', sprintf('p = %.2e', pval)))
  }) %>%
  bind_rows() %>%
  mutate(facetLabel = factor(facetLabel, levels = levels(plotData$facetLabel)))

# Joint interaction LRT, covariate adjusted
pInter <- singleTissueInteractionLRT(tod, plink$X, GENE, VARIANT, TISSUE)$lrt %>%
  filter(str_detect(term, ':')) %>% pull(p.value)
stopifnot(length(pInter) == 1)
print(gtPvals %>% select(genotype, n, p_rhythmicity))
message('Joint interaction LRT p = ', signif(pInter, 3))

p <- ggplot(plotData, aes(x = hours, y = expression)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = 'lm', formula = y ~ sin(I(x * pi / 12)) + cos(I(x * pi / 12)),
              colour = 'steelblue') +
  geom_text(data = gtPvals, aes(x = 12, y = Inf, label = pLabel),
            vjust = 1.5, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ facetLabel, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 24, by = 6), limits = c(0, 24),
                     labels = paste0(seq(0, 24, by = 6), 'h')) +
  labs(x = 'Time of day (h)', y = 'Normalised expression',
       title = paste0('AIF1  |  Heart Left Ventricle  |  ', VARIANT, '  |  ',
                      sprintf('rhyQTL p = %.2e', pInter))) +
  theme_minimal()
ggsave(paste0(SUPFIG_DIR, 'SuppFigure_4_AIF1_ChenComparison.jpeg'), p, width = 10, height = 4, dpi = 300)
ggsave(paste0(SUPFIG_DIR, 'SuppFigure_4_AIF1_ChenComparison.pdf'),  p, width = 10, height = 4)
