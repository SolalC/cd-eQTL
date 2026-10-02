# =============================================================================
# 21  Are the cd-eQTL in Chen et al. Supplementary Data 2? (check of 21_chenComparison)
#
# Looks up every cd-eQTL of Supplementary Table 2 (gene x variant x tissue) in
# the raw Chen files, independently of the window filter and of the rsID join
# used by 21_chenComparison.R. A cd-eQTL is searched for in four ways:
#   by ID        Chen ID equals the GTEx variant ID (chr_pos_ref_alt_b38) or rsID
#   by position  Chen chromosome and position equal the variant's (any alleles)
#   by gene      Chen gene equals the cd-eQTL gene (unversioned Ensembl ID, or
#                the gene symbol in case the IDs differ)
#   any tissue   the same variant / pair in another Chen tissue file
# Each cd-eQTL gets a verdict: the first of
#   pair in Chen               gene x variant found in the same tissue
#   variant, other gene        variant found in the tissue, but for other genes
#   gene, other variants       gene found in the tissue (Chen did not call this variant)
#   other tissue only          pair or variant found only in other tissues
#   absent from tissue         neither gene nor variant in the tissue file
# and, where 21_chenComparison.R has run, the status it gave that row.
#
# Note: Chen Supplementary Data 2 lists their significant rhyQTLs (G-test p <
# 0.05 and HANOVA BH < 0.05), not every test, so "gene, other variants" and
# "absent" can be genuine non-detections rather than a lookup failure.
#
# Run:     sbatch 21_chenCdeQTLcheck.sh       (after 04; 21_chenComparison optional)
# Outputs: Results/published/chenComparison/cdeQTL_inChen_check.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))

OUT_DIR <- paste0(CHEN_DIR, 'chenComparison/')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
geneBase  <- function(x) sub('\\..*', '', x)
addChr    <- function(x) ifelse(startsWith(x, 'chr'), x, paste0('chr', x))
tissueKey <- function(x) tolower(gsub('[^A-Za-z0-9]', '', x))

# --- cd-eQTL (Supplementary Table 2) and lead flags (Table 3) -------------------------
cd <- fread(paste0(TABLE_DIR, 'SupTable_2_cdEQTL.csv')) %>%
  distinct(gene, variant, tissue, .keep_all = TRUE) %>%
  transmute(tissue, gene, geneBase = geneBase(gene), hgnc_symbol, variant, cdeQTL_p = p.value)
lead <- fread(paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv'))
cd <- cd %>% mutate(lead = paste(tissue, gene, variant) %in% paste(lead$tissue, lead$gene, lead$variant))
message('cd-eQTL: ', nrow(cd), ' (', sum(cd$lead), ' lead)')

# --- Their coordinates from the GTEx lookup ---------------------------------------------
lookupCols <- c('variant_id', 'chr', 'pos', 'ref', 'alt', 'rs_id_dbSNP155_GRCh38p13')
lookup <- fread(GENOTYPE_LOOKUP, select = lookupCols)
lookup <- lookup[rs_id_dbSNP155_GRCh38p13 %in% cd$variant]
setnames(lookup, c('variant_id', 'rs_id_dbSNP155_GRCh38p13'), c('gtex_variant_id', 'variant'))
message('cd-eQTL variants found in the lookup: ', n_distinct(lookup$variant), ' / ', n_distinct(cd$variant),
        if (anyDuplicated(lookup$variant)) ' (some rsIDs have several lookup rows; all are searched)' else '')
cdPos <- cd %>% left_join(lookup, by = 'variant', relationship = 'many-to-many')

# --- Raw Chen files ----------------------------------------------------------------------
bedFiles   <- list.files(GTEX_EXPR_DIR, pattern = '.v10.normalized_expression.bed.gz$')
gtexTissue <- str_remove(bedFiles, '.v10.normalized_expression.bed.gz')
chenFiles  <- list.files(CHEN_SUPPDATA2, pattern = '.txt$', full.names = TRUE)
chenTissue <- str_remove(basename(chenFiles), '.txt$')
tissueMap  <- setNames(gtexTissue[match(tissueKey(chenTissue), tissueKey(gtexTissue))], chenTissue)
stopifnot(!anyNA(tissueMap))
chen <- rbindlist(lapply(seq_along(chenFiles), function(k) {
  fread(chenFiles[k], sep = '\t', select = c('ID', 'Chromosome', 'Position', 'REF', 'ALT',
                                             'rhyGene.ID', 'rhyGene.name', 'pval'),
        colClasses = list(character = c('ID', 'Chromosome', 'REF', 'ALT'))) %>%
    mutate(tissue = tissueMap[[chenTissue[k]]])
})) %>%
  mutate(chr = addChr(Chromosome), geneBase = geneBase(rhyGene.ID))
message('Chen rows: ', nrow(chen), '; IDs look like: ', paste(head(unique(chen$ID), 3), collapse = ', '))

# Restrict to rows touching a cd-eQTL gene or variant (keeps the joins small)
chen <- chen[ID %in% c(cdPos$gtex_variant_id, cdPos$variant) |
             paste(chr, Position) %in% paste(cdPos$chr, cdPos$pos) |
             geneBase %in% cd$geneBase | rhyGene.name %in% cd$hgnc_symbol]

# --- Search ---------------------------------------------------------------------------------
chen <- as.data.frame(chen)
searchOne <- function(r) {
  inTissue <- chen[chen$tissue == r$tissue, ]
  ids      <- c(r$gtex_variant_id, r$variant)
  byID     <- inTissue[inTissue$ID %in% ids, ]
  byPos    <- if (is.na(r$pos)) inTissue[0, ] else inTissue[inTissue$chr == r$chr & inTissue$Position == r$pos, ]
  varRows  <- unique(rbind(byID, byPos))
  isGene   <- function(d) d$geneBase == r$geneBase | d$rhyGene.name == r$hgnc_symbol
  geneRows <- inTissue[isGene(inTissue), ]
  pairRows <- varRows[isGene(varRows), ]
  other    <- chen[chen$tissue != r$tissue &
                     (chen$ID %in% ids | (!is.na(r$pos) & chen$chr == r$chr & chen$Position == r$pos)), ]
  data.frame(
    variant_found_by_ID       = nrow(byID) > 0,
    variant_found_by_position = nrow(byPos) > 0,
    chen_rows_for_gene        = nrow(geneRows),
    pair_in_chen              = nrow(pairRows) > 0,
    chen_ID                   = paste(unique(pairRows$ID), collapse = ';'),
    chen_alleles              = paste(unique(paste0(pairRows$REF, rep('>', nrow(pairRows)), pairRows$ALT)), collapse = ';'),
    chen_pval                 = if (nrow(pairRows)) min(pairRows$pval) else NA_real_,
    chen_genes_for_variant    = paste(unique(varRows$rhyGene.name), collapse = ';'),
    other_tissues             = paste(unique(other$tissue), collapse = ';'),
    verdict = if (nrow(pairRows)) 'pair in Chen' else if (nrow(varRows)) 'variant, other gene' else
      if (nrow(geneRows)) 'gene, other variants' else if (nrow(other)) 'other tissue only' else
      'absent from tissue')
}
check <- bind_cols(cdPos, bind_rows(lapply(seq_len(nrow(cdPos)), function(i) searchOne(cdPos[i, ]))))

# Several lookup rows per rsID: keep the best-matching one per cd-eQTL
rank <- c('pair in Chen', 'variant, other gene', 'gene, other variants', 'other tissue only', 'absent from tissue')
check <- check %>%
  arrange(tissue, gene, variant, match(verdict, rank)) %>%
  distinct(tissue, gene, variant, .keep_all = TRUE)

# --- Status given by 21_chenComparison.R, where it has run ---------------------------------
inWinFile <- paste0(OUT_DIR, 'chen_inWindow.csv')
check$comparison_status <- NA_character_
if (file.exists(inWinFile)) {
  status <- fread(inWinFile, select = c('tissue', 'rhyGene.ID', 'rsID', 'status')) %>%
    transmute(tissue, geneBase = geneBase(rhyGene.ID), variant = rsID, status) %>%
    distinct(tissue, geneBase, variant, .keep_all = TRUE)
  check <- check %>%
    left_join(status, by = c('tissue', 'geneBase', 'variant')) %>%
    mutate(comparison_status = coalesce(status, 'not in chen_inWindow.csv')) %>%
    select(-status)
}

fwrite(check %>% select(-geneBase), paste0(OUT_DIR, 'cdeQTL_inChen_check.csv'))
print(check %>% count(verdict, lead))
print(check %>% filter(pair_in_chen) %>% count(comparison_status))
message('cd-eQTL pairs present in Chen but not "tested in both" in the comparison: ',
        sum(check$pair_in_chen & check$comparison_status %in% c(NA, 'not in chen_inWindow.csv') |
              check$pair_in_chen & !check$comparison_status %in% c(NA, 'tested in both')))
