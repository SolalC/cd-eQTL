# =============================================================================
# 06  Figure 2: cd-eQTL interaction heatmap (gene x tissue, -log10 interaction p)
#
# Each cell is the smallest interaction LRT p of the gene in that tissue. For
# cd-eQTL genes not rhythmic in a tissue (so not scanned by 02), the lead
# variant is tested there with the same three-model LRT ("backfill"), written to
# QTLregression/missing/ and reused on later runs. Cells whose gene is not
# expressed in the tissue stay white. Stars mark only the 54 cd-eQTL of the main
# scan; backfilled cells reaching p < 1e-4 are not cd-eQTL and are not starred.
# CCG labels are drawn in red.
#
# Run:     sbatch 06_figure2_heatmap.sh      (after 04; loads genotypes only
#          if some backfill combinations have not been computed yet)
# Output:  Figures/Figure_2_cdEQTL_heatmap.jpeg                 Figure 2
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)

readLRT <- function(dir) {
  bind_rows(lapply(list.files(dir, full.names = TRUE), function(file) {
    fread(file) %>%
      mutate(tissue = str_remove(basename(file), '_LRT.csv')) %>%
      filter(str_detect(term, ':'))
  }))
}

# --- Interaction LRTs from the main scan ------------------------------------
interactionFull <- readLRT(paste0(QTL_DIR, '4_LRT'))
interaction.f <- interactionFull %>%
  filter(p.value < CDEQTL_P) %>%
  group_by(tissue, gene) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup()
message('cd-eQTL (gene x tissue, p < ', CDEQTL_P, '): ', nrow(interaction.f))

# Symbols as in Supplementary Table 3
gene_mapping <- readIndependentCdeQTL() %>% distinct(gene, hgnc_symbol)
stopifnot(all(interaction.f$gene %in% gene_mapping$gene))
symbolise <- function(d) {
  d %>% left_join(gene_mapping, by = 'gene') %>%
    mutate(hgnc_symbol = ifelse(hgnc_symbol == '' | is.na(hgnc_symbol), gene, hgnc_symbol))
}
interaction.f <- symbolise(interaction.f)

bestPerCell <- function(d) {
  d %>% filter(gene %in% interaction.f$gene) %>%
    group_by(tissue, gene) %>%
    slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    symbolise()
}
interactionFull.f <- bestPerCell(interactionFull)
allTissues <- unique(interactionFull.f$tissue)

# --- Backfill: lead variant of each cd-eQTL gene in every tissue ---------------
lrtOutDir  <- paste0(QTL_DIR, 'missing/4_LRT/')
coefOutDir <- paste0(QTL_DIR, 'missing/3_harmonicInteraction/')
dir.create(lrtOutDir,  recursive = TRUE, showWarnings = FALSE)
dir.create(coefOutDir, recursive = TRUE, showWarnings = FALSE)

alreadyRun <- bind_rows(
  interactionFull %>% distinct(gene, variant, tissue),
  bind_rows(lapply(list.files(lrtOutDir, full.names = TRUE), function(f)
    fread(f) %>% mutate(tissue = str_remove(basename(f), '_LRT.csv')) %>%
      distinct(gene, variant, tissue))))
missingCombos <- interaction.f %>%
  distinct(gene, variant) %>%
  tidyr::crossing(tissue = allTissues) %>%
  anti_join(alreadyRun, by = c('gene', 'variant', 'tissue'))

# Keep only combinations whose gene is expressed in the tissue; the others
# cannot be tested and stay white.
expressed <- lapply(setNames(unique(missingCombos$tissue), unique(missingCombos$tissue)),
  function(t) sub('\\..*', '', fread(paste0(GTEX_EXPR_DIR, t, '.v10.normalized_expression.bed.gz'),
                                    select = 'gene_id')$gene_id))
missingCombos <- missingCombos %>%
  filter(map2_lgl(gene, tissue, ~ sub('\\..*', '', .x) %in% expressed[[.y]]))
message(nrow(missingCombos), ' backfill combinations to compute')

if (nrow(missingCombos) > 0) {
  plink <- loadGenotypes()
  tod   <- loadDIP() %>% filter(!is.na(radian))
  for (i in seq_len(nrow(missingCombos))) {
    row <- missingCombos[i, ]
    message(sprintf('[%d/%d] %s  %s  %s', i, nrow(missingCombos), row$gene, row$variant, row$tissue))
    singleTissueInteractionLRT(tod, plink$X, gene = row$gene, variant = row$variant,
                               tissueName = row$tissue,
                               lrtOutDir = lrtOutDir, coefOutDir = coefOutDir)
  }
  rm(plink); gc()
}

interactionFullCombined <- bind_rows(interactionFull.f, bestPerCell(readLRT(lrtOutDir))) %>%
  group_by(tissue, gene) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup()

interaction_full <- expand.grid(hgnc_symbol = unique(interactionFullCombined$hgnc_symbol),
                                tissue      = unique(interactionFullCombined$tissue)) %>%
  left_join(interactionFullCombined, by = c('hgnc_symbol', 'tissue')) %>%
  mutate(tissue = str_replace_all(tissue, '_', ' '))

# --- Clustering and plot --------------------------------------------------------
heatmap_matrix <- interaction_full %>%
  select(hgnc_symbol, tissue, p.value) %>%
  spread(key = tissue, value = p.value)
heatmap_matrix[, -1] <- -log10(heatmap_matrix[, -1])
heatmap_matrix[is.na(heatmap_matrix)] <- 0

row_clust <- hclust(dist(heatmap_matrix[, -1]))
col_clust <- hclust(dist(t(heatmap_matrix[, -1])))

interaction_full$hgnc_symbol <- factor(interaction_full$hgnc_symbol,
                                       levels = heatmap_matrix$hgnc_symbol[row_clust$order])
interaction_full$tissue      <- factor(interaction_full$tissue,
                                       levels = colnames(heatmap_matrix)[-1][col_clust$order])
interaction_full <- interaction_full %>% mutate(logP = -log10(p.value))

# Star only the cd-eQTL cells from the main scan.
stars <- interaction_full %>%
  semi_join(interaction.f %>% mutate(tissue = str_replace_all(tissue, '_', ' ')),
            by = c('gene', 'tissue'))
message('Stars drawn: ', nrow(stars))

ccg_symbols  <- interaction.f %>% filter(gene %in% CCG) %>% pull(hgnc_symbol) %>% unique()
gene_order   <- levels(interaction_full$hgnc_symbol)
x_label_cols <- ifelse(gene_order %in% ccg_symbols, 'darkred', 'black')

heatmapPlot <- ggplot(interaction_full, aes(x = hgnc_symbol, y = tissue, fill = logP)) +
  geom_tile(color = 'black', linewidth = 0.3) +
  geom_text(label = '*', data = stars) +
  scale_fill_viridis_c(na.value = 'white') +
  coord_fixed() +
  labs(x = '', y = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = x_label_cols),
        panel.grid  = element_blank())
ggsave(heatmapPlot, filename = paste0(FIG_DIR, 'Figure_2_cdEQTL_heatmap.jpeg'),
       width = 10, height = 10, dpi = 300)
