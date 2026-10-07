# =============================================================================
# 23  LD pruning of the cd-eQTL
#
# The cd-eQTL of Supplementary Table 2 (interaction LRT p < 1e-4) are clumped
# with plink 1.9 on the interaction p-value (--clump-p1 1 --clump-p2 1
# --clump-kb 1000), at r2 < 0.1 and r2 < 0.01, on the genotype set of the scan
# (all genotyped donors). Two scopes:
#   gene x tissue  one clump run per gene in each tissue (keeps at least the
#                  54 leads, one per gene x tissue)
#   all            one run on the unique variants of all cd-eQTL, across genes
#                  and tissues, each variant taking its smallest p. This is the
#                  number of LD-independent signals.
#
# Run:     sbatch 23_ldPruning.sh          (after 04; needs plink 1.9 on PATH,
#                                           or its path in the PLINK variable)
# Output:  Results/inferred/ldPruning/cdeQTL_LDpruned.csv
#            one row per cd-eQTL; for each threshold t and scope, whether it
#            is an index variant (ld<t>_kept, ld<t>_all_kept) and the index
#            variant of its clump (ld<t>_lead, ld<t>_all_lead; itself if kept)
#          Results/inferred/ldPruning/plink/   plink inputs and logs
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

OUT_DIR   <- paste0(RESULTS_DIR, 'inferred/ldPruning/')
PLINK_DIR <- paste0(OUT_DIR, 'plink/')
dir.create(PLINK_DIR, recursive = TRUE, showWarnings = FALSE)
PLINK         <- Sys.getenv('PLINK', 'plink')
R2_THRESHOLDS <- c(0.1, 0.01)
CLUMP_KB      <- 1000

runPlink <- function(args, out) {
  status <- system2(PLINK, c(args, '--out', out), stdout = FALSE, stderr = FALSE)
  if (status != 0) stop('plink failed, see ', out, '.log')
}

# plink --clump on one set of variants; returns, per variant, its index variant
clump <- function(variant, p, threshold, tag) {
  out <- paste0(PLINK_DIR, tag)
  fwrite(data.frame(SNP = variant, P = p), paste0(out, '.assoc'), sep = '\t')
  runPlink(c('--bfile', paste0(PLINK_DIR, 'cdeQTL'), '--clump', paste0(out, '.assoc'),
             '--clump-p1', 1, '--clump-p2', 1, '--clump-r2', threshold,
             '--clump-kb', CLUMP_KB), out)
  cl <- fread(paste0(out, '.clumped'))
  members <- strsplit(gsub('\\(1\\)', '', cl$SP2), ',')
  lead <- rbind(data.frame(variant = cl$SNP, lead = cl$SNP),
                data.frame(variant = unlist(members),
                           lead = rep(cl$SNP, lengths(members)))) %>%
    filter(variant != 'NONE')
  if (!setequal(lead$variant, variant)) stop('plink did not return every variant of ', tag)
  lead
}

cdeqtl <- fread(paste0(TABLE_DIR, 'SupTable_2_cdEQTL.csv')) %>%
  distinct(gene, variant, tissue, .keep_all = TRUE)
snps <- cdeqtl %>% group_by(variant) %>% summarise(p = min(p.value))
message('cd-eQTL: ', nrow(cdeqtl), ' (', nrow(snps), ' variants)')

# Genotypes of the cd-eQTL variants only
writeLines(snps$variant, paste0(PLINK_DIR, 'cdeQTL_variants.txt'))
runPlink(c('--bfile', GENOTYPE_PLINK, '--extract', paste0(PLINK_DIR, 'cdeQTL_variants.txt'),
           '--make-bed'), paste0(PLINK_DIR, 'cdeQTL'))
missing <- setdiff(snps$variant, fread(paste0(PLINK_DIR, 'cdeQTL.bim'), header = FALSE)$V2)
if (length(missing)) stop('Variants not in the plink set: ', paste(missing, collapse = ', '))

out <- cdeqtl
for (t in R2_THRESHOLDS) {
  pruned <- cdeqtl %>%
    group_by(gene, tissue) %>%
    group_modify(~ clump(.x$variant, .x$p.value, t,
                         paste0('r2_', t, '_', .y$gene, '_', .y$tissue))) %>%
    ungroup()
  prunedAll <- clump(snps$variant, snps$p, t, paste0('r2_', t, '_all'))
  pfx <- paste0('ld', t, '_'); pfxAll <- paste0('ld', t, '_all_')
  out <- out %>%
    left_join(pruned %>% rename(!!paste0(pfx, 'lead') := lead),
              by = c('gene', 'tissue', 'variant')) %>%
    left_join(prunedAll %>% rename(!!paste0(pfxAll, 'lead') := lead), by = 'variant') %>%
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
