# =============================================================================
# 23  LD pruning of the cd-eQTL
#
# The cd-eQTL of Supplementary Table 2 (interaction LRT p < 1e-4) are pruned
# within each gene x tissue by greedy clumping on the interaction p-value: the
# variant with the smallest p is kept, every variant with r2 >= threshold with
# it is removed, and the next smallest remaining p is kept, until none is left.
# This is plink --clump with p1 = p2 = 1 and no distance limit. r2 is the
# squared Pearson correlation of alternate-allele dosages over all genotyped
# donors. Pruning is run at r2 < 0.1 and r2 < 0.01.
#
# Run:     sbatch 23_ldPruning.sh          (after 04; loads the full plink set)
# Output:  Results/inferred/ldPruning/cdeQTL_LDpruned.csv
#            one row per cd-eQTL; for each threshold t, whether it is kept
#            (ld<t>_kept), the kept variant it was pruned by (ld<t>_lead,
#            itself if kept) and its r2 with that variant (ld<t>_r2)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)

OUT_DIR <- paste0(RESULTS_DIR, 'inferred/ldPruning/')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
R2_THRESHOLDS <- c(0.1, 0.01)

# Greedy clumping of the variants of one gene x tissue. G: dosage matrix,
# variants in rows. Returns, per variant, the kept variant it belongs to and r2.
clump <- function(variant, p, G, threshold) {
  r2   <- cor(t(G[variant, , drop = FALSE]), use = 'pairwise.complete.obs')^2
  lead <- setNames(rep(NA_character_, length(variant)), variant)
  for (v in variant[order(p)]) {
    if (!is.na(lead[v])) next
    lead[v] <- v
    lead[is.na(lead) & r2[v, names(lead)] >= threshold] <- v
  }
  data.frame(variant = variant, lead = unname(lead),
             r2 = r2[cbind(variant, unname(lead))])
}

cdeqtl <- fread(paste0(TABLE_DIR, 'SupTable_2_cdEQTL.csv')) %>%
  distinct(gene, variant, tissue, .keep_all = TRUE)
message('cd-eQTL: ', nrow(cdeqtl), ' (', n_distinct(cdeqtl$variant), ' variants)')

X <- loadGenotypes()$X
missing <- setdiff(unique(cdeqtl$variant), rownames(X))
if (length(missing)) stop('Variants not in the plink set: ', paste(missing, collapse = ', '))
G <- X[unique(cdeqtl$variant), , drop = FALSE]
rm(X); gc()

out <- cdeqtl
for (t in R2_THRESHOLDS) {
  pruned <- cdeqtl %>%
    group_by(gene, tissue) %>%
    group_modify(~ clump(.x$variant, .x$p.value, G, t)) %>%
    ungroup()
  pfx <- paste0('ld', t, '_')
  out <- out %>%
    left_join(pruned %>% rename_with(~ paste0(pfx, .x), c(lead, r2)),
              by = c('gene', 'tissue', 'variant')) %>%
    mutate(!!paste0(pfx, 'kept') := variant == .data[[paste0(pfx, 'lead')]])
}
fwrite(out, paste0(OUT_DIR, 'cdeQTL_LDpruned.csv'))

for (t in R2_THRESHOLDS) {
  kept <- out %>% filter(.data[[paste0('ld', t, '_kept')]])
  message(sprintf('r2 < %s: %d of %d cd-eQTL kept, %d genes, %d tissues, %d gene x tissue',
                  t, nrow(kept), nrow(out), n_distinct(kept$gene), n_distinct(kept$tissue),
                  n_distinct(kept$gene, kept$tissue)))
}
