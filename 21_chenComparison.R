# =============================================================================
# 21  Chen et al. (2025) rhythmic eQTL tests vs the cd-eQTL scan
#
# Part 1  Chen et al. Supplementary Data 2 (one file per tissue) restricted to
#         the variants within the cd-eQTL scan window of their gene: TSS +/-
#         TSS_WINDOW (50 kb), with the TSS from the GTEx v10 expression BED and
#         the same bounds as extractSNPs (start - window <= pos <= end + window).
#         Chen positions are assumed GRCh38 (GTEx v8), like the BED.
# Part 2  Per tissue, the in-window Chen tests merged with the interaction LRT
#         of the cd-eQTL scan (02) on tissue x gene (unversioned Ensembl ID) x
#         variant (rsID). Chen variants are mapped to rsIDs through the GTEx
#         lookup table by chr, position, REF and ALT (REF/ALT swapped as a
#         fallback, flagged in allele_match).
#         Each in-window Chen test gets a status saying whether it was tested in
#         the scan and, if not, the first reason it was not: the variant is not
#         in the lookup, the gene is not in the GTEx v10 tissue, the gene is not
#         rhythmic (Bonferroni, so not scanned), the variant is not in the
#         genotype subset, it is not in the scan's variant list, or it failed
#         the heterozygote-frequency filter (MAF_THRESHOLD).
#         `cdeQTL_call` turns this into a first-failed call for our pipeline:
#         Step 0 lookup, Step 1 tissue/gene (01), Step 2 variant selection (02),
#         Step 3 interaction LRT p >= CDEQTL_P, else cd-eQTL.
#
# Genotype coding: Chen genotype groups 0/1/2 and the scan dosage may count
# different alleles; check allele orientation before comparing effect
# directions. The merged table carries p-values only.
#
# Run:     sbatch 21_chenComparison.sh        (after 01 and 02)
# Inputs:  Results/published/SupplementaryData2/<Chen tissue>.txt
# Outputs: Results/published/chenComparison/
#            chen_inWindow.csv           Chen tests within the window, with status
#            chen_cdeQTL_merged.csv      tests present in both analyses
#            byTissue/<tissue>_merged.csv  the same, one file per tissue
#            chen_cdeQTL_summary.csv     counts per tissue and cdeQTL_call
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
OUT_DIR <- paste0(CHEN_DIR, 'chenComparison/')
dir.create(paste0(OUT_DIR, 'byTissue'), recursive = TRUE, showWarnings = FALSE)

geneBase <- function(x) sub('\\..*', '', x)
addChr   <- function(x) ifelse(startsWith(x, 'chr'), x, paste0('chr', x))
tissueKey <- function(x) tolower(gsub('[^A-Za-z0-9]', '', x))

# --- Tissue names: Chen file names -> GTEx v10 names ------------------------------
bedFiles   <- list.files(GTEX_EXPR_DIR, pattern = '.v10.normalized_expression.bed.gz$', full.names = TRUE)
gtexTissue <- str_remove(basename(bedFiles), '.v10.normalized_expression.bed.gz')
chenFiles  <- list.files(CHEN_SUPPDATA2, pattern = '.txt$', full.names = TRUE)
if (!length(chenFiles)) stop('No Chen files in ', CHEN_SUPPDATA2)
chenTissue <- str_remove(basename(chenFiles), '.txt$')
tissueMap  <- setNames(gtexTissue[match(tissueKey(chenTissue), tissueKey(gtexTissue))], chenTissue)
if (anyNA(tissueMap)) stop('Chen files without a GTEx tissue: ',
                           paste(names(tissueMap)[is.na(tissueMap)], collapse = ', '))
message(length(chenFiles), ' Chen tissue files, all matched to GTEx tissues')

# --- Chen et al. tests ------------------------------------------------------------
chen <- rbindlist(lapply(seq_along(chenFiles), function(k) {
  fread(chenFiles[k], sep = '\t', colClasses = list(character = c('ID', 'Chromosome', 'REF', 'ALT'))) %>%
    mutate(tissue = tissueMap[[chenTissue[k]]])
})) %>%
  mutate(chr = addChr(Chromosome), geneBase = geneBase(rhyGene.ID))
message('Chen tests: ', nrow(chen))

# --- Part 1: TSS window -------------------------------------------------------------
# TSS of every gene in any v10 tissue (identical across tissues). Joining on
# chromosome as well keeps PAR genes (same base ID on chrX and chrY) apart.
tss <- rbindlist(lapply(bedFiles, fread, select = 1:4)) %>%
  distinct() %>%
  rename(tss_chr = `#chr`, tss_start = start, tss_end = end) %>%
  mutate(geneBase = geneBase(gene_id))
stopifnot(!anyDuplicated(tss[, c('geneBase', 'tss_chr')]))

chen <- chen %>%
  left_join(tss, by = c('geneBase', 'chr' = 'tss_chr')) %>%
  mutate(distance_to_TSS = Position - tss_end,
         inWindow = !is.na(tss_end) &
           Position >= tss_start - TSS_WINDOW & Position <= tss_end + TSS_WINDOW)
message('Chen tests whose gene has no TSS on the same chromosome in GTEx v10: ', sum(is.na(chen$tss_end)))
message('Chen tests within TSS +/- ', TSS_WINDOW / 1000, ' kb: ', sum(chen$inWindow))

inWin <- chen %>% filter(inWindow) %>% select(-inWindow)

# --- Variant IDs: chr / pos / REF / ALT -> rsID ---------------------------------------
lookupCols <- c('variant_id', 'chr', 'pos', 'ref', 'alt', 'rs_id_dbSNP155_GRCh38p13')
lookup <- fread(GENOTYPE_LOOKUP, select = lookupCols)
stopifnot(all(lookupCols %in% colnames(lookup)))
lookup <- lookup[paste(chr, pos) %in% paste(inWin$chr, inWin$Position)]
setnames(lookup, c('variant_id', 'rs_id_dbSNP155_GRCh38p13'), c('gtex_variant_id', 'rsID'))

same    <- lookup %>% select(chr, Position = pos, REF = ref, ALT = alt, gtex_variant_id, rsID) %>%
  mutate(allele_match = 'same')
swapped <- lookup %>% select(chr, Position = pos, REF = alt, ALT = ref, gtex_variant_id, rsID) %>%
  mutate(allele_match = 'swapped')
alleles <- bind_rows(same, swapped) %>%
  filter(rsID != '.') %>%
  distinct(chr, Position, REF, ALT, .keep_all = TRUE)   # 'same' kept over 'swapped'

inWin <- inWin %>% left_join(alleles, by = c('chr', 'Position', 'REF', 'ALT'))
message('In-window Chen tests mapped to an rsID: ', sum(!is.na(inWin$rsID)), ' / ', nrow(inWin),
        ' (', sum(inWin$allele_match == 'swapped', na.rm = TRUE), ' with REF/ALT swapped)')
rsChen <- inWin %>% filter(startsWith(ID, 'rs'), !is.na(rsID))
if (nrow(rsChen)) message('Chen rsID differs from the mapped rsID: ',
                          sum(rsChen$ID != rsChen$rsID), ' / ', nrow(rsChen))
rm(lookup, same, swapped, alleles)

# --- cd-eQTL scan -------------------------------------------------------------------
genotyped <- fread(paste0(GENOTYPE_PLINK, '.bim'), header = FALSE, select = 2)[[1]]   # rsIDs

scanTissue <- function(t) {
  harmonicFile <- paste0(HARMONIC_DIR, t, '_LRT_harmonic_linear.csv')
  snpFile      <- paste0(SNP_DIR, t, '_SNPs.csv')
  lrtFile      <- paste0(QTL_DIR, '4_LRT/', t, '_LRT.csv')
  expressed <- rhythmic <- character(0)
  candidates <- lrt <- NULL
  if (file.exists(harmonicFile)) {
    harmonic  <- fread(harmonicFile) %>% mutate(p.adj = p.adjust(p.value, 'bonferroni'))
    expressed <- geneBase(harmonic$gene)
    rhythmic  <- geneBase(harmonic$gene[harmonic$p.adj < 0.05])
  }
  if (file.exists(snpFile))
    candidates <- fread(snpFile) %>% transmute(geneBase = geneBase(gene_id), rsID = rs_id_dbSNP155_GRCh38p13)
  if (file.exists(lrtFile))
    lrt <- fread(lrtFile) %>%
      filter(str_detect(term, ':')) %>%
      distinct(gene, variant, .keep_all = TRUE) %>%
      mutate(p.adj = p.adjust(p.value, method = 'BH')) %>%   # within tissue, as in 04
      transmute(geneBase = geneBase(gene), rsID = variant,
                cdeQTL_gene_id = gene, cdeQTL_LRT_statistic = statistic,
                cdeQTL_p = p.value, cdeQTL_p.adj_BH = p.adj, cdeQTL = p.value < CDEQTL_P)
  list(scanned = file.exists(harmonicFile), expressed = expressed, rhythmic = rhythmic,
       candidates = candidates, lrt = lrt)
}

# --- Part 2: merge per tissue ----------------------------------------------------------
annotated <- lapply(sort(unique(inWin$tissue)), function(t) {
  message('Tissue: ', t)
  s <- scanTissue(t)
  d <- inWin %>% filter(tissue == t)
  if (!is.null(s$lrt)) d <- d %>% left_join(s$lrt, by = c('geneBase', 'rsID'))
  else d <- d %>% mutate(cdeQTL_gene_id = NA_character_, cdeQTL_LRT_statistic = NA_real_,
                         cdeQTL_p = NA_real_, cdeQTL_p.adj_BH = NA_real_, cdeQTL = NA)
  candidatePair <- if (is.null(s$candidates)) rep(FALSE, nrow(d)) else
    paste(d$geneBase, d$rsID) %in% paste(s$candidates$geneBase, s$candidates$rsID)
  d %>% mutate(status = case_when(
    !is.na(cdeQTL_p)               ~ 'tested in both',
    is.na(rsID)                    ~ 'variant not in GTEx lookup',
    !s$scanned                     ~ 'tissue not in 01 output',
    !geneBase %in% s$expressed     ~ 'gene not in GTEx v10 tissue',
    !geneBase %in% s$rhythmic      ~ 'gene not rhythmic (not scanned)',
    !rsID %in% genotyped           ~ 'variant not in genotype subset',
    !candidatePair                 ~ 'not in scan variant list',
    TRUE                           ~ 'heterozygote frequency filter'))
})
inWin <- bind_rows(annotated) %>%
  mutate(cdeQTL_call = case_when(   # our pipeline, labelled by the first step that drops the test
    status == 'variant not in GTEx lookup'      ~ 'Step 0: variant not in GTEx lookup',
    status == 'tissue not in 01 output'         ~ 'Step 1: tissue not in 01 output',
    status == 'gene not in GTEx v10 tissue'     ~ 'Step 1: gene not in GTEx v10 tissue',
    status == 'gene not rhythmic (not scanned)' ~ 'Step 1: gene not rhythmic (Bonferroni)',
    status == 'variant not in genotype subset'  ~ 'Step 2: variant not in genotype subset',
    status == 'not in scan variant list'        ~ 'Step 2: not in scan variant list',
    status == 'heterozygote frequency filter'   ~ 'Step 2: heterozygote frequency filter',
    cdeQTL_p >= CDEQTL_P                        ~ paste0('Step 3: interaction LRT p >= ', CDEQTL_P),
    TRUE                                        ~ paste0('cd-eQTL (interaction LRT p < ', CDEQTL_P, ')')))

fwrite(inWin %>% select(-geneBase), paste0(OUT_DIR, 'chen_inWindow.csv'))

merged <- inWin %>%
  filter(status == 'tested in both') %>%
  select(tissue, gene_id = cdeQTL_gene_id, gene_name = rhyGene.name, rsID, gtex_variant_id,
         chr, Position, REF, ALT, allele_match, distance_to_TSS,
         chen_ID = ID, chen_rhyGene.ID = rhyGene.ID,
         chen_Sample.size.0 = Sample.size.0, chen_Sample.size.1 = Sample.size.1,
         chen_Sample.size.2 = Sample.size.2,
         chen_pval_0 = pval_0, chen_phase_0 = phase_0, chen_amp_0 = amp_0,
         chen_pval_1 = pval_1, chen_phase_1 = phase_1, chen_amp_1 = amp_1,
         chen_Chosen.model = Chosen.model, chen_G.test_pval = G.test_pval, chen_pval = pval,
         cdeQTL_LRT_statistic, cdeQTL_p, cdeQTL_p.adj_BH, cdeQTL, cdeQTL_call)
fwrite(merged, paste0(OUT_DIR, 'chen_cdeQTL_merged.csv'))
for (t in unique(merged$tissue))
  fwrite(merged %>% filter(tissue == t), paste0(OUT_DIR, 'byTissue/', t, '_merged.csv'))

tissueSummary <- chen %>%
  count(tissue, name = 'chen_tests') %>%
  left_join(chen %>% filter(inWindow) %>% count(tissue, name = 'in_window'), by = 'tissue') %>%
  left_join(inWin %>% dplyr::count(tissue, cdeQTL_call) %>%
              pivot_wider(names_from = cdeQTL_call, values_from = n), by = 'tissue') %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0L)))
fwrite(tissueSummary, paste0(OUT_DIR, 'chen_cdeQTL_summary.csv'))

message('Tests in both analyses: ', nrow(merged), ' in ', n_distinct(merged$tissue), ' tissues')
print(inWin %>% dplyr::count(cdeQTL_call, sort = TRUE))
