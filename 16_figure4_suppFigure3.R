# =============================================================================
# 16  Figure 4 and Supplementary Figure 3: cross-tissue cosine similarity of the
#     BMAL1 (a) and TMED10 (b) cd-eQTL effects, over the 49 rhythmic tissues
#
#   Figure 4                all four coefficients (matrices from 14)
#   Supplementary Figure 3  interaction terms only (matrices from 15)
# Tissues are ordered by hierarchical clustering.
#
# Run:     Rscript 16_figure4_suppFigure3.R    (after 15; seconds)
# Outputs: Figures/Figure_4_cosineSimilarity.{jpeg,pdf}
#          Supplementary_Figures/SuppFigure_3_cosineSimilarityInteraction.{jpeg,pdf}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(corrplot)
options(bitmapType = 'cairo')

readSim <- function(f) {
  m <- as.matrix(fread(f))
  colnames(m) <- str_replace_all(colnames(m), '_', ' ')
  rownames(m) <- colnames(m)
  m
}

drawPage <- function(mats) {
  par(mfrow = c(1, 2), oma = c(0, 0, 2.5, 0))
  for (k in seq_len(nrow(EXAMPLES))) {
    ex <- EXAMPLES[k, ]
    corrplot(mats[[k]], method = 'square', order = 'hclust', type = 'lower',
             diag = FALSE, col = rev(COL2('RdBu')), tl.cex = 0.7, tl.col = 'black',
             title = sprintf('Type: %s\nGene: %s\nVariant: %s', ex$type, ex$symbol, ex$variant),
             mar = c(0, 0, 5, 0), cex.main = 1.4)
    mtext(letters[k], side = 3, outer = TRUE, at = (k - 1) / 2 + 0.01,
          adj = 0, line = 0.5, font = 2, cex = 2.2)
  }
}

savePage <- function(mats, stem) {
  jpeg(paste0(stem, '.jpeg'), width = 22, height = 12.5, unit = 'in', res = 300, type = 'cairo')
  drawPage(mats); dev.off()
  cairo_pdf(paste0(stem, '.pdf'), width = 22, height = 12.5)
  drawPage(mats); dev.off()
}

for (d in c('cosineSimilarity', 'cosineSimilarityInteraction')) {
  mats <- lapply(seq_len(nrow(EXAMPLES)), function(k)
    readSim(paste0(RHYTHM_DIR, d, '/', EXAMPLES$gene[k], '_', EXAMPLES$variant[k], '.csv')))
  for (k in seq_along(mats))
    message(d, ' ', EXAMPLES$symbol[k], ': ', ncol(mats[[k]]), ' tissues, mean similarity ',
            round(mean(mats[[k]][lower.tri(mats[[k]])]), 3))
  savePage(mats, if (d == 'cosineSimilarity') paste0(FIG_DIR, 'Figure_4_cosineSimilarity')
                 else paste0(SUPFIG_DIR, 'SuppFigure_3_cosineSimilarityInteraction'))
}
