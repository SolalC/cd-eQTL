# =============================================================================
# 15  Cross-tissue similarity of cd-eQTL effects: CCG vs CRG
#
# For each independent cd-eQTL, the tissue x tissue cosine similarity of the
# interaction-model coefficients (from 14) over the 49 rhythmic tissues, under
# two definitions:
#   all four coefficients     sin, cos, variant:sin, variant:cos   (Figure 4)
#   interaction terms only    variant:sin, variant:cos             (Supplementary Figure 3)
# The mean pairwise similarity per locus is compared between the 8 CCG and the
# 46 CRG loci with a Wilcoxon rank-sum test.
#
# Run:     Rscript 15_tissueSimilarity.R    (after 14; minutes)
# Outputs: Results/inferred/rhythmicityParameters/cosineSimilarityInteraction/<gene>_<variant>.csv
#          Results/inferred/rhythmicityParameters/cosineSimilarity_perLocus.csv
#          Results/inferred/rhythmicityParameters/cosineSimilarity_CCGvsCRG_tests.csv
#              mean +- sd per group, W and p (Results text)
#          Supplementary_Tables/SupTable_11_BMAL1_cosineSimilarity.csv    Supplementary Table 11
#          Supplementary_Tables/SupTable_12_TMED10_cosineSimilarity.csv   Supplementary Table 12
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(lsa)

inDir  <- paste0(RHYTHM_DIR, 'multiTissueReg/')
simDir <- paste0(RHYTHM_DIR, 'cosineSimilarity/')
intDir <- paste0(RHYTHM_DIR, 'cosineSimilarityInteraction/')
dir.create(intDir, recursive = TRUE, showWarnings = FALSE)

FULL_TERMS  <- c('sin(radian)', 'cos(radian)', 'variant:sin(radian)', 'variant:cos(radian)')
INTER_TERMS <- c('variant:sin(radian)', 'variant:cos(radian)')
RHYTHMIC    <- rhythmicTissues()

# (length(terms) x nTissue) coefficient matrix, dropping tissues with a missing
# coefficient so both definitions use the same tissues.
coefMatrix <- function(d, terms) {
  wide <- d %>%
    filter(term %in% terms) %>%
    select(tissue, term, estimate) %>%
    pivot_wider(id_cols = term, names_from = tissue, values_from = estimate)
  m <- as.matrix(wide[, -1, drop = FALSE])
  rownames(m) <- wide$term
  m <- m[terms[terms %in% rownames(m)], , drop = FALSE]
  m[, colSums(is.na(m)) == 0, drop = FALSE]
}
meanLower <- function(mat) {
  v <- mat[lower.tri(mat, diag = FALSE)]
  c(mean(v), sd(v))
}

sim <- bind_rows(lapply(list.files(inDir, full.names = TRUE), function(f) {
  key <- str_remove(basename(f), '\\.csv$')
  d   <- fread(f)
  m4  <- coefMatrix(d, FULL_TERMS)
  m2  <- coefMatrix(d, INTER_TERMS)
  common <- intersect(intersect(colnames(m4), colnames(m2)), RHYTHMIC)
  s4 <- lsa::cosine(m4[, common, drop = FALSE])
  s2 <- lsa::cosine(m2[, common, drop = FALSE])
  fwrite(as.data.frame(s2), paste0(intDir, key, '.csv'))
  a4 <- meanLower(s4); a2 <- meanLower(s2)
  data.frame(key = key, gene = unique(d$gene), variant = unique(d$variant),
             circType = ifelse(unique(d$gene) %in% CCG, 'CCG', 'CRG'),
             nTissue = length(common),
             mean_full4 = a4[1], sd_full4 = a4[2],
             mean_inter2 = a2[1], sd_inter2 = a2[2])
}))
fwrite(sim, paste0(RHYTHM_DIR, 'cosineSimilarity_perLocus.csv'))

testOne <- function(col, label) {
  ccg <- sim[[col]][sim$circType == 'CCG']
  crg <- sim[[col]][sim$circType == 'CRG']
  wt  <- wilcox.test(crg, ccg, paired = FALSE)
  data.frame(definition = label,
             n_CCG = length(ccg), n_CRG = length(crg),
             mean_CCG = mean(ccg), sd_CCG = sd(ccg),
             mean_CRG = mean(crg), sd_CRG = sd(crg),
             W = unname(wt$statistic), p = wt$p.value)
}
tests <- bind_rows(testOne('mean_full4',  'all four coefficients'),
                   testOne('mean_inter2', 'interaction terms only'))
fwrite(tests, paste0(RHYTHM_DIR, 'cosineSimilarity_CCGvsCRG_tests.csv'))
print(tests, digits = 3)

# --- Supplementary Tables 11-12: example-locus similarity matrices ---------------
for (k in seq_len(nrow(EXAMPLES))) {
  ex <- EXAMPLES[k, ]
  m  <- as.matrix(fread(paste0(simDir, ex$gene, '_', ex$variant, '.csv')))
  stopifnot(ncol(m) == 49)
  tissues <- str_replace_all(colnames(m), '_', ' ')
  colnames(m) <- tissues
  fwrite(cbind(gene = ex$gene, hgnc_symbol = ex$symbol, tissue = tissues, as.data.frame(m)),
         paste0(TABLE_DIR, 'SupTable_', 10 + k, '_', ex$symbol, '_cosineSimilarity.csv'))
}
