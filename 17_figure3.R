# =============================================================================
# 17  Figure 3 panels: BMAL1 rs11022718 (Artery Aorta, left) and TMED10
#     rs76379942 (Colon Sigmoid, right)
#
#   a-b  expression (covariate-residualised) against DIP by genotype, with the
#        per-genotype cosinor fit and a dashed line at each peak; and the
#        per-genotype rhythm vector with its 95% confidence region (Bingham
#        et al. eq. 33) on a 24 h clock, from the residualised per-genotype
#        cosinors compared by the classification tests (07)
#   c-d  amplitude by tissue and genotype (interaction model, from 14)
#   e-f  acrophase by genotype across tissues (interaction model, from 14)
# The panels were assembled into the final figure by hand.
#
# Run:     sbatch 17_figure3.sh     (after 14)
# Outputs: Figures/Figure_3{a,b}_<gene>_expression.{jpeg,pdf}
#          Figures/Figure_3{a,b}_<gene>_confidenceClock.{jpeg,pdf}
#          Figures/Figure_3{c,d}_<gene>_amplitude.{jpeg,pdf}
#          Figures/Figure_3{e,f}_<gene>_phase.{jpeg,pdf}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)
library(ggExtra)
library(MetBrewer)
options(bitmapType = 'cairo')

genoCols <- setNames(met.brewer(n = 3, name = 'Degas'), c('0', '1', '2'))
panels   <- list(expression = c('3a', '3b'), amplitude = c('3c', '3d'), phase = c('3e', '3f'))

# --- a-b: residualised expression ------------------------------------------------
# Expression residualised on the covariates, as in the Bingham tests.
residualise <- function(data) {
  covTerms <- c(if ('sex' %in% colnames(data)) 'sex',
                'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  covTerms <- covTerms[covTerms %in% colnames(data)]
  data$gene    <- residuals(lm(as.formula(paste('gene ~', paste(covTerms, collapse = '+'))),
                               data = data, na.action = na.exclude))
  data$variant <- factor(round(data$variant), c(0, 1, 2))
  data %>% filter(!is.na(gene), !is.na(variant))
}

plotResidualisedExpression <- function(cov, tissueName, variant, geneTitleName) {
  params <- do.call(rbind, lapply(levels(cov$variant), function(g) {
    d <- cov[cov$variant == g, ]
    if (nrow(d) < 3) return(NULL)
    cc <- coef(lm(gene ~ cos(radian) + sin(radian), data = d))
    sc <- cc[['sin(radian)']]; co <- cc[['cos(radian)']]
    ph <- atan2(sc, co) %% (2 * pi)
    data.frame(variant = factor(g, levels = levels(cov$variant)), n = nrow(d),
               amplitude = sqrt(sc^2 + co^2),
               peak_radian = ph, peak_hour = ph / (2 * pi) * 24)
  }))
  annoCaption <- paste(sprintf('G%s: peak %.1f h, amp %.2f', as.character(params$variant),
                               params$peak_hour, params$amplitude), collapse = '    ')
  p <- ggplot(cov, aes(x = radian, y = gene, col = variant, fill = variant)) +
    theme_minimal() +
    geom_point(aes(fill = variant)) +
    geom_smooth(method = 'lm', formula = y ~ cos(x) + sin(x)) +
    ggtitle(paste0(geneTitleName, '_', tissueName, '_', variant, ' (residualised)')) +
    labs(y = 'residualised expression', caption = annoCaption) +
    theme(title = element_text(hjust = 0.5), legend.position = 'bottom',
          plot.caption = element_text(hjust = 0.5)) +
    scale_x_continuous(name = 'Time of day (h)',
                       breaks = seq(0, 2 * pi, by = pi / 4), labels = seq(0, 24, by = 3)) +
    geom_vline(data = params, aes(xintercept = peak_radian, colour = variant),
               linetype = 'dashed', linewidth = 0.4, show.legend = FALSE) +
    scale_fill_manual(values = genoCols) +
    scale_color_manual(values = genoCols)
  list(plot = ggMarginal(p, type = 'histogram', margins = 'x', groupColour = T, groupFill = T,
                         binwidth = 2 * pi / 12),
       params = params)
}

# --- a-b: confidence clock -----------------------------------------------------------
# Rhythm vector (x = sin coef, y = cos coef) drawn at hour-angle atan2(sin, cos)
# clockwise from the top (the true acrophase), radius = amplitude, with the
# joint 95% confidence region (theta - theta_hat)' Sigma^-1 (theta - theta_hat)
# <= 2 F_0.95(2, nu).
binghamClockPlot <- function(gp, nu, conf = 0.95) {
  rFac <- sqrt(2 * qf(conf, 2, nu))
  ell <- do.call(rbind, lapply(seq_len(nrow(gp)), function(j) {
    mu    <- c(gp$sin_coef[j], gp$cos_coef[j])
    Sigma <- matrix(c(gp$cov_ss[j], gp$cov_sc[j], gp$cov_sc[j], gp$cov_cc[j]), 2, 2)
    e   <- eigen(Sigma, symmetric = TRUE)
    th  <- seq(0, 2 * pi, length.out = 180)
    pts <- mu + rFac * (e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% rbind(cos(th), sin(th)))
    data.frame(genotype = factor(gp$genotype[j]), x = pts[1, ], y = pts[2, ])
  }))
  centres <- data.frame(genotype = factor(gp$genotype), x = gp$sin_coef, y = gp$cos_coef)
  rMax <- max(sqrt(ell$x^2 + ell$y^2))
  clockXY <- function(hour, r) { a <- hour / 24 * 2 * pi; data.frame(x = r * sin(a), y = r * cos(a)) }
  ringR <- pretty(c(0, rMax), n = 4); ringR <- ringR[ringR > 0 & ringR <= rMax]
  rings <- do.call(rbind, lapply(ringR, function(r) {
    d <- clockXY(seq(0, 24, length.out = 200), r); d$r <- r; d }))
  ringLab <- cbind(clockXY(0, ringR), lab = format(ringR, digits = 2))
  hrs    <- seq(0, 21, by = 3)
  spokes <- do.call(rbind, lapply(hrs, function(h) {
    e <- clockXY(h, rMax); data.frame(x = 0, y = 0, xend = e$x, yend = e$y) }))
  hrLab  <- cbind(clockXY(hrs, rMax * 1.13), lab = paste0(hrs, 'h'))

  ggplot() +
    coord_equal(clip = 'off') + theme_void() +
    theme(legend.position = 'none') +
    geom_path(data = rings, aes(x, y, group = r), colour = 'grey85') +
    geom_segment(data = spokes, aes(x = x, y = y, xend = xend, yend = yend), colour = 'grey85') +
    geom_text(data = hrLab, aes(x, y, label = lab), size = 3, colour = 'grey40') +
    geom_text(data = ringLab, aes(x, y, label = lab), size = 2.6, colour = 'grey55',
              hjust = -0.15, vjust = -0.3) +
    geom_polygon(data = ell, aes(x, y, fill = genotype, group = genotype), alpha = 0.25, colour = NA) +
    geom_segment(data = centres, aes(x = 0, y = 0, xend = x, yend = y, colour = genotype),
                 arrow = ggplot2::arrow(length = grid::unit(0.15, 'cm'))) +
    geom_point(data = centres, aes(x, y, colour = genotype), size = 2) +
    scale_colour_manual(values = genoCols) + scale_fill_manual(values = genoCols)
}

tod   <- loadDIP()
plink <- loadGenotypes()$X
for (k in seq_len(nrow(EXAMPLES))) {
  ex  <- EXAMPLES[k, ]
  res <- singleTissueInteractionLRT(tod, plink, ex$gene, ex$variant, ex$tissue)

  expr <- plotResidualisedExpression(residualise(res$data), ex$tissue, ex$variant, ex$symbol)
  message(ex$symbol, ' residualised per-genotype fits:'); print(expr$params, digits = 3)
  gp    <- residualisedGenotypeParameters(res$data)
  clock <- binghamClockPlot(gp, nu = sum(gp$df))

  stem <- paste0(FIG_DIR, 'Figure_', panels$expression[k], '_', ex$symbol)
  for (ext in c('jpeg', 'pdf')) {
    ggsave(paste0(stem, '_expression.', ext), expr$plot, width = 9, height = 5, dpi = 300)
    ggsave(paste0(stem, '_confidenceClock.', ext), clock, width = 4.5, height = 4.5, dpi = 300)
  }
}
rm(plink); gc()

# --- c-f: amplitude and acrophase across tissues (interaction model) ---------------
RHYTHMIC <- rhythmicTissues()
for (k in seq_len(nrow(EXAMPLES))) {
  ex <- EXAMPLES[k, ]
  rhythmicity <- calculate_rhythmicity_parameters(
    fread(paste0(RHYTHM_DIR, 'multiTissueReg/', ex$gene, '_', ex$variant, '.csv')))
  # Tissues ordered by clustering of the per-genotype amplitudes
  amp_wide <- rhythmicity %>%
    dplyr::select(tissue, variant, genotype, amplitude) %>%
    unite('feature', variant, genotype, sep = '_') %>%
    pivot_wider(names_from = feature, values_from = amplitude) %>%
    as.data.frame()
  amp_mat <- as.matrix(amp_wide[, -1]); rownames(amp_mat) <- amp_wide$tissue
  hc <- hclust(dist(amp_mat, method = 'euclidean'), method = 'ward.D2')
  rhythmicity$tissue <- factor(rhythmicity$tissue, levels = hc$label[hc$order])
  r <- rhythmicity %>% filter(tissue %in% RHYTHMIC)
  l <- paste0(ex$symbol, '_', ex$variant)
  if (max(r$amplitude, na.rm = TRUE) > 1.3)
    warning(l, ': amplitude exceeds the 1.3 colour limit; tiles above it render grey')

  ampPlot <- ggplot(r) +
    geom_tile(aes(x = tissue, y = genotype, fill = amplitude)) +
    scale_fill_viridis_c(limits = c(0, 1.3)) +
    scale_y_reverse() +
    ggtitle(paste0(l, '_amplitude')) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  # Acrophase (radians, -pi..pi) as time of day: wrap to [0, 2*pi), rescale to 0-24 h
  phasePlot <- ggplot(r) +
    geom_histogram(aes(x = (phase %% (2 * pi)) / (2 * pi) * 24, fill = tissue)) +
    facet_wrap(~ genotype, nrow = 1) +
    scale_x_continuous(name = 'Acrophase (time of day, h)',
                       breaks = seq(0, 24, by = 6), limits = c(0, 24)) +
    scale_fill_manual(values = met.brewer('Renoir', 50)) +
    ylim(0, 26) +
    ggtitle(paste0(l, '_phase')) +
    theme_minimal() +
    theme(legend.position = 'none')

  for (ext in c('jpeg', 'pdf')) {
    ggsave(paste0(FIG_DIR, 'Figure_', panels$amplitude[k], '_', ex$symbol, '_amplitude.', ext),
           ampPlot, width = 14, height = 7, dpi = 300)
    ggsave(paste0(FIG_DIR, 'Figure_', panels$phase[k], '_', ex$symbol, '_phase.', ext),
           phasePlot, width = 7, height = 7, dpi = 300)
  }
}
