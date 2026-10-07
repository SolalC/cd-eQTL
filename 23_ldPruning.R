# =============================================================================
# 23  LD pruning of the cd-eQTL
#
# The cd-eQTL of Supplementary Table 2 (interaction LRT p < 1e-4) are pruned by
# greedy clumping on the interaction p-value: the variant with the smallest p
# is kept, every variant with r2 >= threshold with it is removed, and the next
# smallest remaining p is kept, until none is left. This is plink --clump with
# p1 = p2 = 1 and --clump-kb 1000. r2 is the squared Pearson correlation of
# alternate-allele dosages over all genotyped donors; pairs on different
# chromosomes or more than CLUMP_KB apart are not compared. Pruning is run at
# r2 < 0.1 and r2 < 0.01, in two scopes:
#   gene x tissue  within each gene in each tissue (keeps at least the 54
#                  leads, one per gene x tissue)
#   all            the unique variants of all 167 cd-eQTL together, across
#                  genes and tissues, each variant taking its smallest p. This
#                  is the number of LD-independent signals.
#
# Run:     sbatch 23_ldPruning.sh          (after 04; loads the full plink set)
# Output:  Results/inferred/ldPruning/cdeQTL_LDpruned.csv
#            one row per cd-eQTL; for each threshold t and scope, whether it
#            is kept (ld<t>_kept, ld<t>_all_kept), the kept variant it was
#            pruned by (ld<t>_lead, ld<t>_all_lead; itself if kept) and its r2
#            with that variant (ld<t>_r2, ld<t>_all_r2)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)

OUT_DIR <- paste0(RESULTS_DIR, 'inferred/ldPruning/')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
R2_THRESHOLDS <- c(0.1, 0.01)
CLUMP_KB      <- 1000

# Greedy clumping. G: dosage matrix, variants in rows; chr and pos: named by
# variant. Returns, per variant, the kept variant it belongs to and its r2.
clump <- function(variant, p, G, threshold, chr, pos) {
  r2   <- cor(t(G[variant, , drop = FALSE]), use = 'pairwise.complete.obs')^2
  far  <- outer(chr[variant], chr[variant], '!=') |
          abs(outer(pos[variant], pos[variant], '-')) > CLUMP_KB * 1000
  r2[far] <- 0
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

plink <- loadGenotypes()
missing <- setdiff(unique(cdeqtl$variant), rownames(plink$X))
if (length(missing)) stop('Variants not in the plink set: ', paste(missing, collapse = ', '))
G   <- plink$X[unique(cdeqtl$variant), , drop = FALSE]
chr <- setNames(as.character(plink$bim$chr), plink$bim$id)[rownames(G)]
pos <- setNames(plink$bim$pos, plink$bim$id)[rownames(G)]
rm(plink); gc()

snps <- cdeqtl %>% group_by(variant) %>% summarise(p = min(p.value))

out <- cdeqtl
for (t in R2_THRESHOLDS) {
  pruned <- cdeqtl %>%
    group_by(gene, tissue) %>%
    group_modify(~ clump(.x$variant, .x$p.value, G, t, chr, pos)) %>%
    ungroup()
  prunedAll <- clump(snps$variant, snps$p, G, t, chr, pos)
  pfx <- paste0('ld', t, '_'); pfxAll <- paste0('ld', t, '_all_')
  out <- out %>%
    left_join(pruned %>% rename_with(~ paste0(pfx, .x), c(lead, r2)),
              by = c('gene', 'tissue', 'variant')) %>%
    left_join(prunedAll %>% rename_with(~ paste0(pfxAll, .x), c(lead, r2)), by = 'variant') %>%
    mutate(!!paste0(pfx, 'kept')    := variant == .data[[paste0(pfx, 'lead')]],
           !!paste0(pfxAll, 'kept') := variant == .data[[paste0(pfxAll, 'lead')]])
}
fwrite(out, paste0(OUT_DIR, 'cdeQTL_LDpruned.csv'))

for (t in R2_THRESHOLDS) {
  kept <- out %>% filter(.data[[paste0('ld', t, '_kept')]])
  message(sprintf('r2 < %s, within gene x tissue: %d of %d cd-eQTL kept, %d genes, %d tissues',
                  t, nrow(kept), nrow(out), n_distinct(kept$gene), n_distinct(kept$tissue)))
  keptAll <- unique(out$variant[out[[paste0('ld', t, '_all_kept')]]])
  message(sprintf('r2 < %s, across genes and tissues: %d of %d variants kept (independent signals)',
                  t, length(keptAll), nrow(snps)))
}
