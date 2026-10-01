# =============================================================================
# 12  Supplementary Figure 1: dryR BIC weight of every candidate model at the
#     54 independent cd-eQTL
#
# One row per locus, one column per candidate model, cell colour = BIC weight
# (sums to 1 within a locus), star = selected model. Loci with two usable
# genotype groups (5 models, panel a) and three (15 models, panel b) are drawn
# separately because the same model index means different models in each. The
# model where all groups share one rhythm, and the loci selecting it, are in red.
#
# Run:     sbatch 12_suppFigure1_dryRmodelGrid.sh   (after 11; needs patchwork >= 1.2)
# Outputs: Supplementary_Figures/SuppFigure_1_dryRModelSelection.{jpeg,pdf}
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(patchwork)

perLocus <- fread(paste0(DRYR_DIR, 'dryR_perLocus.csv'))

modelDesc <- list(
  `2` = c('flat in both', 'rhythmic in 0 only', 'rhythmic in 1 only',
          'one rhythm shared by 0 and 1', '0 and 1 rhythmic, separate'),
  `3` = c('flat in all three', 'rhythmic in 0 only', 'rhythmic in 1 only', 'rhythmic in 2 only',
          '0=1 shared, 2 flat', '0 and 1 separate, 2 flat', '0=2 shared, 1 flat',
          '0 and 2 separate, 1 flat', '1=2 shared, 0 flat', '1 and 2 separate, 0 flat',
          'all three share one rhythm', '0 alone, 1=2 shared', '0=2 shared, 1 alone',
          '0=1 shared, 2 alone', 'all three separate'))
nModels     <- c(`2` = 5L, `3` = 15L)
sharedModel <- c(`2` = 4L, `3` = 11L)
SHARED_COL  <- '#CC0000'
XLIM <- c(0.5, 15.5)   # shared by both panels so the columns align

# --- Full BIC weight profile from the saved dryR objects -----------------------
# dryR_raw.rds is in Supplementary Table 3 row order. drylm may have been fed a
# 3-row duplicated matrix, so every row of BICW_rhythm is the same fit.
raw  <- readRDS(paste0(DRYR_DIR, 'dryR_raw.rds'))
loci <- readIndependentCdeQTL() %>% distinct(gene, variant, tissue)
stopifnot(length(raw) == nrow(loci))

bicw <- bind_rows(lapply(seq_along(raw), function(i) {
  o <- raw[[i]]
  if (is.null(o) || is.null(o$BICW_rhythm)) return(NULL)
  w <- as.numeric(o$BICW_rhythm[1, ])
  data.table(gene = loci$gene[i], variant = loci$variant[i], tissue = loci$tissue[i],
             model = seq_along(w), BICW = w)
})) %>%
  left_join(perLocus %>% select(gene, variant, tissue, n_groups, dryR_chosen_model, hgnc_symbol),
            by = c('gene', 'variant', 'tissue')) %>%
  mutate(selected = model == dryR_chosen_model)

# The argmax must be the recorded model and the weights must sum to 1
chk <- bicw %>% group_by(gene, variant, tissue) %>%
  summarise(argmax_ok = model[which.max(BICW)] == first(dryR_chosen_model),
            sums_to_1 = abs(sum(BICW) - 1) < 1e-6,
            nmod_ok   = n() == nModels[as.character(first(n_groups))],
            .groups = 'drop')
stopifnot(all(chk$argmax_ok), all(chk$sums_to_1), all(chk$nmod_ok))

bicw <- bicw %>%
  mutate(label = ifelse(hgnc_symbol == '' | is.na(hgnc_symbol), gene, hgnc_symbol),
         locus = paste0(label, '  ', variant, '  ', str_replace_all(tissue, '_', ' ')))

# --- Panels -------------------------------------------------------------------
ord <- bicw %>% distinct(locus, n_groups, dryR_chosen_model, label) %>%
  arrange(n_groups, dryR_chosen_model, label)

panel <- function(k, tag) {
  d <- bicw %>% filter(n_groups == k) %>%
    mutate(locus = factor(locus, levels = rev(ord$locus[ord$n_groups == k])))
  shared <- sharedModel[as.character(k)]
  xCol <- ifelse(seq_len(nModels[as.character(k)]) == shared, SHARED_COL, 'black')
  yCol <- ifelse(levels(d$locus) %in% d$locus[d$selected & d$model == shared],
                 SHARED_COL, 'black')
  ggplot(d, aes(x = model, y = locus)) +
    geom_tile(aes(fill = BICW), colour = 'white', linewidth = 0.3) +
    geom_text(data = filter(d, selected), label = '*', size = 2.9,
              colour = 'black', vjust = 0.72) +
    scale_fill_viridis_c(name = 'BIC weight   ', limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(position = 'top', limits = XLIM, expand = c(0, 0),
                       breaks = seq_len(nModels[as.character(k)]),
                       labels = paste0(seq_len(nModels[as.character(k)]), '. ',
                                       modelDesc[[as.character(k)]])) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = tag) +
    theme_minimal(base_size = 9) +
    theme(panel.grid      = element_blank(),
          axis.text.y     = element_text(size = 6, family = 'mono', colour = yCol),
          axis.text.x.top = element_text(size = 6.6, angle = 45, hjust = 0, vjust = 0,
                                         colour = xCol),
          axis.ticks      = element_blank(),
          plot.tag        = element_text(size = 12, face = 'bold', hjust = 0),
          plot.tag.position = c(0, 1),
          plot.margin     = margin(14, 74, 10, 2))
}

n2 <- sum(ord$n_groups == 2); n3 <- sum(ord$n_groups == 3)
message('Loci: ', nrow(ord), '  (2-group ', n2, ', 3-group ', n3, ')')

fig <- suppressWarnings(
  (panel(2L, 'a') / panel(3L, 'b')) +
    plot_layout(heights = c(n2, n3), guides = 'collect') &
    theme(legend.position   = 'bottom',
          legend.direction  = 'horizontal',
          legend.title      = element_text(size = 8, vjust = 1),
          legend.text       = element_text(size = 7),
          legend.key.height = unit(8, 'pt'), legend.key.width = unit(58, 'pt')))

ggsave(paste0(SUPFIG_DIR, 'SuppFigure_1_dryRModelSelection.jpeg'), fig,
       width = 11.5, height = 11, dpi = 400)
ggsave(paste0(SUPFIG_DIR, 'SuppFigure_1_dryRModelSelection.pdf'), fig,
       width = 11.5, height = 11)
