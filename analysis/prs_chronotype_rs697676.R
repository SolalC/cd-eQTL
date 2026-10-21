# =============================================================================
# Chronotype polygenic score (PGS000336) in the GTEx donors, by rs697676 genotype
#
# Score: PGS Catalog PGS000336 ("GRS313_Chron"; Maukonen et al. 2020, J Biol
#   Rhythms, doi:10.1177/0748730420935328). 313 genome-wide significant variants
#   from the Jones et al. 2019 chronotype GWAS (GCST007565, 23andMe + UK Biobank),
#   weights = ln(OR) of the GWAS binary trait (morning person). A higher score
#   therefore means more morning-associated alleles, if the GWAS coding was kept;
#   check the sign on a variant of known direction before interpreting it.
#   Original build GRCh37. GTEx WGS is GRCh38, so the PGS Catalog harmonised
#   GRCh38 file (hm_chr, hm_pos) is used; the original file is a fallback,
#   matched by rsID only.
#
# Score per donor: PRS = sum_j w_j * d_j, d_j = copies of the effect allele.
#   Variants are matched to the .bim by GRCh38 chr:pos and alleles, then by rsID.
#   Allele orientation: effect allele = bim A1 (counted) -> d = g; = A2 -> d = 2 - g;
#   strand-flipped (complement) alleles are accepted for non-palindromic SNPs.
#   Palindromic SNPs (A/T, C/G) are assumed to be on the forward strand (both the
#   harmonised file and the GTEx calls are GRCh38 forward strand); set
#   DROP_PALINDROMIC=TRUE to exclude them instead. Missing genotypes are replaced
#   by the variant mean (as plink --score). Unmatched variants are listed.
#
# Figure: PRS (y, z-score) by number of copies of the counted allele of rs697676
#   (x, 0/1/2; bim A1, as in eqtl_PER3_DBP_allTissues.R). rs697676 lies in PER3
#   (chr1p36). If rs697676 or a variant in LD with it is in the score, the
#   association is partly circular, so the plot is repeated with a score that
#   excludes every score variant within LOCUS_WINDOW (default 500 kb) of rs697676.
#   The p-value is the slope of PRS ~ dosage (lm, additive).
#
# Caveat: the score was built and evaluated in Europeans; GTEx donors are ~85%
#   European ancestry. Score means differ by ancestry for reasons unrelated to
#   chronotype, so compare within ancestry or adjust for genetic PCs.
#
# Run:     Rscript prs_chronotype_rs697676.R [output_dir]
# Needs:   R >= 4.0, data.table, ggplot2; internet access on first run (the
#          scoring file is downloaded to output_dir and reused afterwards)
# Output (in OUT_DIR, default <ROOT>/Results/PRS/):
#          PGS000336_hmPOS_GRCh38.txt.gz          downloaded scoring file
#          PGS000336_variant_matching.csv         one row per score variant (match, orientation)
#          PGS000336_GTEx_scores.csv              one row per donor (PRS, rs697676 dosage)
#          PGS000336_by_rs697676.pdf              boxplots
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --- Paths (as in eqtl_PER3_DBP_allTissues.R; override with env vars) ---------
ROOT           <- Sys.getenv('CIRCADIAN_ROOT', '/scratch/project_mnt/S0007/solal/circadian/')
GENOTYPE_PLINK <- Sys.getenv('GENOTYPE_PLINK', paste0(
  '/QRISdata/Q8106/Controlled/Genotype/filtered/',
  'GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv_cQTL_subset_inferred'))
DIP_FILE       <- Sys.getenv('DIP_FILE', file.path(ROOT, 'Data/CHIRAL/DIP.RData'))

args    <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args)) args[1] else Sys.getenv('OUT_DIR', file.path(ROOT, 'Results/PRS'))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PGS_ID   <- 'PGS000336'
PGS_FTP  <- sprintf('https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/%s/ScoringFiles/', PGS_ID)
PGS_HM   <- paste0(PGS_FTP, 'Harmonized/', PGS_ID, '_hmPOS_GRCh38.txt.gz')
PGS_ORIG <- paste0(PGS_FTP, PGS_ID, '.txt.gz')

SNP              <- 'rs697676'
LOCUS_WINDOW     <- as.numeric(Sys.getenv('LOCUS_WINDOW', 5e5))
DROP_PALINDROMIC <- as.logical(Sys.getenv('DROP_PALINDROMIC', 'FALSE'))

# --- Download the scoring file (once) -------------------------------------------
# Harmonised GRCh38 first; the original GRCh37 file if that fails.
fetch <- function(url) {
  dest <- file.path(OUT_DIR, basename(url))
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  ok <- tryCatch({ download.file(url, dest, mode = 'wb', quiet = TRUE); TRUE },
                 error = function(e) { message('Download failed: ', url, ' (', conditionMessage(e), ')'); FALSE },
                 warning = function(w) { message('Download failed: ', url, ' (', conditionMessage(w), ')'); FALSE })
  if (ok && file.exists(dest) && file.size(dest) > 0) dest else { unlink(dest); NULL }
}
scoreFile <- fetch(PGS_HM); harmonised <- !is.null(scoreFile)
if (!harmonised) scoreFile <- fetch(PGS_ORIG)
if (is.null(scoreFile)) stop('Could not download ', PGS_ID, '. Download it by hand from ',
                             PGS_FTP, ' into ', OUT_DIR)
message('Scoring file: ', scoreFile, if (harmonised) ' (harmonised GRCh38)' else
        ' (original GRCh37: matching by rsID only)')

# --- Read the scoring file -------------------------------------------------------
# '#' lines are the PGS Catalog header. Harmonised columns: hm_rsID, hm_chr,
# hm_pos, hm_inferOtherAllele (when other_allele was not reported).
readScore <- function(f) {
  l <- readLines(gzfile(f))
  s <- fread(text = l[!startsWith(l, '#')], sep = '\t', colClasses = 'character')
  col <- function(...) { for (x in c(...)) if (x %in% names(s)) return(s[[x]]); rep(NA_character_, nrow(s)) }
  oa <- col('other_allele')
  if (all(is.na(oa) | oa == '')) oa <- col('hm_inferOtherAllele')
  oa[grepl('/', oa) | oa == ''] <- NA   # several possible other alleles: leave unknown
  rs <- col('hm_rsID'); rs[is.na(rs) | rs == ''] <- col('rsID')[is.na(rs) | rs == '']
  data.table(rsID = rs,
             chr = sub('^chr', '', col('hm_chr')),
             pos = suppressWarnings(as.integer(col('hm_pos'))),
             effect_allele = toupper(col('effect_allele')),
             other_allele = toupper(oa),
             weight = as.numeric(col('effect_weight')))
}
score <- readScore(scoreFile)
message(nrow(score), ' variants in ', PGS_ID)
SNP_IN_SCORE <- SNP %in% score$rsID
message(SNP, if (SNP_IN_SCORE) ' IS in the score: the PRS by genotype comparison is partly circular'
        else ' is not in the score (by rsID); LD neighbours are handled by LOCUS_WINDOW')

# --- Match the score to the .bim --------------------------------------------------
bim <- fread(paste0(GENOTYPE_PLINK, '.bim'), header = FALSE,
             col.names = c('chr', 'id', 'posg', 'pos', 'a1', 'a2'),
             colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character'))
bim[, `:=`(row = .I, chr = sub('^chr', '', chr), a1 = toupper(a1), a2 = toupper(a2))]

comp <- function(a) chartr('ACGT', 'TGCA', a)
isPalindromic <- function(x, y) nchar(x) == 1 & nchar(y) == 1 & comp(x) == y

# Orientation of one score variant against candidate bim rows. Returns the first
# compatible row with sign = +1 (effect = A1) or -1 (effect = A2), else NULL.
orient <- function(ea, oa, cand) {
  for (i in seq_len(nrow(cand))) {
    b <- cand[i]
    for (flip in c(FALSE, TRUE)) {
      e <- if (flip) comp(ea) else ea
      o <- if (is.na(oa)) NA else if (flip) comp(oa) else oa
      pal <- isPalindromic(b$a1, b$a2)
      if (flip && (pal || nchar(ea) > 1)) next   # no strand flip for palindromic SNPs or indels
      if (e == b$a1 && (is.na(o) || o == b$a2))
        return(list(row = b$row, sign = 1, flip = flip, palindromic = pal, bim_id = b$id))
      if (e == b$a2 && (is.na(o) || o == b$a1))
        return(list(row = b$row, sign = -1, flip = flip, palindromic = pal, bim_id = b$id))
    }
  }
  NULL
}

setkey(bim, chr, pos)
byId <- bim[id %in% score$rsID]
match <- rbindlist(lapply(seq_len(nrow(score)), function(j) {
  s <- score[j]
  cand <- if (!is.na(s$pos)) bim[.(s$chr, s$pos), nomatch = 0L] else bim[0]
  how  <- 'chr:pos'
  m <- if (nrow(cand)) orient(s$effect_allele, s$other_allele, cand)
  if (is.null(m) && !is.na(s$rsID)) {
    cand <- byId[id == s$rsID]; how <- 'rsID'
    m <- if (nrow(cand)) orient(s$effect_allele, s$other_allele, cand)
  }
  if (is.null(m)) return(cbind(s, matched_by = NA_character_, bim_row = NA_integer_,
                               bim_id = NA_character_, sign = NA_real_, flip = NA,
                               palindromic = NA))
  cbind(s, matched_by = how, bim_row = m$row, bim_id = m$bim_id, sign = m$sign,
        flip = m$flip, palindromic = m$palindromic)
}))
match[, used := !is.na(bim_row) & !(DROP_PALINDROMIC & palindromic %in% TRUE)]
# The same bim row matched twice (duplicated score entries): keep the first
match[used == TRUE & duplicated(bim_row), used := FALSE]
message(sprintf('Matched %d / %d variants (%d by position, %d by rsID, %d strand-flipped, %d palindromic); used %d',
                sum(!is.na(match$bim_row)), nrow(match), sum(match$matched_by == 'chr:pos', na.rm = TRUE),
                sum(match$matched_by == 'rsID', na.rm = TRUE), sum(match$flip %in% TRUE),
                sum(match$palindromic %in% TRUE), sum(match$used)))
if (sum(match$used) < 0.8 * nrow(match))
  warning('Fewer than 80% of the score variants are in ', basename(GENOTYPE_PLINK),
          '. If this plink set is a variant subset, use the full WGS set (GENOTYPE_PLINK).')

# --- Read genotypes of given .bim rows (as readPlinkVariants) ----------------------
# SNP-major .bed: codes 00/10/11/01 = 2/1/0/NA copies of the bim A1 allele.
readBedRows <- function(prefix, rows) {
  fam <- fread(paste0(prefix, '.fam'), header = FALSE)
  n <- nrow(fam); nb <- ceiling(n / 4)
  con <- file(paste0(prefix, '.bed'), 'rb'); on.exit(close(con))
  stopifnot(identical(readBin(con, 'raw', 3), as.raw(c(0x6c, 0x1b, 0x01))))
  lut <- c(2L, NA_integer_, 1L, 0L)
  G <- vapply(rows, function(r) {
    seek(con, 3 + (r - 1) * nb)
    bytes <- as.integer(readBin(con, 'raw', nb))
    codes <- as.vector(rbind(bytes %% 4L, (bytes %/% 4L) %% 4L,
                             (bytes %/% 16L) %% 4L, (bytes %/% 64L) %% 4L))[seq_len(n)]
    lut[codes + 1L]
  }, integer(n))
  G <- matrix(G, nrow = n)
  rownames(G) <- paste0('GTEX-', vapply(strsplit(fam$V2, '-'), `[`, '', 2))
  G[!duplicated(rownames(G)), , drop = FALSE]
}

# --- Score ----------------------------------------------------------------------
use <- match[used == TRUE]
G   <- readBedRows(GENOTYPE_PLINK, use$bim_row)
D   <- sweep(G, 2, ifelse(use$sign > 0, 0, 2), function(g, k) abs(k - g))   # copies of effect allele
nMissing <- rowSums(is.na(D))
D   <- apply(D, 2, function(d) { d[is.na(d)] <- mean(d, na.rm = TRUE); d })
D   <- matrix(D, ncol = nrow(use), dimnames = list(rownames(G), NULL))

# rs697676 (counted allele = bim A1, as in eqtl_PER3_DBP_allTissues.R)
snpRow <- bim[id == SNP]
if (nrow(snpRow) == 0) stop(SNP, ' not in ', basename(GENOTYPE_PLINK), '.bim')
snpRow <- snpRow[1]
snpG   <- readBedRows(GENOTYPE_PLINK, snpRow$row)[, 1]
message(sprintf('%s: chr%s:%d, counted allele %s, other allele %s', SNP, snpRow$chr, snpRow$pos,
                snpRow$a1, snpRow$a2))

# Score variants near rs697676 (same chromosome, within LOCUS_WINDOW): removed
# from the locus-excluded score
bimUse <- bim[order(row)][use$bim_row]
nearLocus <- bimUse$chr == snpRow$chr & abs(bimUse$pos - snpRow$pos) <= LOCUS_WINDOW
match[used == TRUE, near_rs697676 := nearLocus]
message(sum(nearLocus), ' score variant(s) within ', LOCUS_WINDOW / 1e3, ' kb of ', SNP,
        if (any(nearLocus)) paste0(': ', paste(use$rsID[nearLocus], collapse = ', ')))

prs      <- as.vector(D %*% use$weight)
prsLocus <- as.vector(D[, !nearLocus, drop = FALSE] %*% use$weight[!nearLocus])
res <- data.table(SUBJID = rownames(D), prs = prs, prs_excl_locus = prsLocus,
                  n_variants = nrow(use), n_missing = nMissing,
                  rs697676 = snpG[rownames(D)], counted_allele = snpRow$a1)
res[, `:=`(prs_z = as.vector(scale(prs)), prs_excl_locus_z = as.vector(scale(prs_excl_locus)))]
if (file.exists(DIP_FILE)) {
  load(DIP_FILE)   # provides `phi`
  res[, has_DIP := SUBJID %in% paste0('GTEX-', names(phi))[!is.na(phi)]]
}

fwrite(match, file.path(OUT_DIR, paste0(PGS_ID, '_variant_matching.csv')))
fwrite(res,   file.path(OUT_DIR, paste0(PGS_ID, '_GTEx_scores.csv')))
message('Wrote ', OUT_DIR, '/', PGS_ID, '_{variant_matching,GTEx_scores}.csv (', nrow(res), ' donors)')

# --- Test and figure --------------------------------------------------------------
pd <- res[!is.na(rs697676)]
trend <- function(y) summary(lm(y ~ pd$rs697676))$coefficients[2, ]
tests <- data.table(score = c('prs_z', 'prs_excl_locus_z'),
                    label = c(sprintf('Full score (%d variants)', nrow(use)),
                              sprintf('Without %d variant(s) within %g kb of %s',
                                      sum(nearLocus), LOCUS_WINDOW / 1e3, SNP)))
tests[, c('beta', 'se', 't', 'p') := as.list(trend(pd[[score]])), by = score]
print(tests)

# Figure: locus-excluded score only, with the fitted slope. The label gives the
# slope and its standard error (SD of the score per allele), not the SD of the data.
GENO_COLS <- c('0' = '#96410e', '1' = '#e5a335', '2' = '#556219')
tx <- tests[score == 'prs_excl_locus_z']
pd[, z := prs_excl_locus_z]
counts <- pd[, .(n = .N, y = max(z)), by = rs697676]
fitLine <- data.table(x = c(-0.3, 2.3))[, y := mean(pd$z) + tx$beta * (x - mean(pd$rs697676))]
yTop <- max(pd$z) + 0.18 * diff(range(pd$z))

pl <- ggplot(pd, aes(rs697676, z)) +
  geom_boxplot(aes(group = rs697676, fill = factor(rs697676)), outlier.shape = NA,
               width = 0.6, alpha = 0.5) +
  geom_jitter(aes(colour = factor(rs697676)), width = 0.15, height = 0, size = 0.6, alpha = 0.6) +
  geom_line(data = fitLine, aes(x, y), colour = 'grey20', linewidth = 0.7) +
  geom_text(data = counts, aes(y = y, label = paste0('n = ', n)), vjust = -0.6, size = 2.8) +
  annotate('text', x = 1, y = yTop, size = 3.2, parse = TRUE,   # plotmath: no Unicode font needed
           label = sprintf("beta %%+-%% SE == '%.3f' %%+-%% '%.3f' * ',' ~ italic(p) == '%.2g'",
                           tx$beta, tx$se, tx$p)) +
  scale_x_continuous(breaks = 0:2) +
  scale_fill_manual(values = GENO_COLS, guide = 'none') +
  scale_colour_manual(values = GENO_COLS, guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(title = 'PRS: Chronotype\nrs697676 excluded',
       x = sprintf('%s: copies of the %s allele', SNP, snpRow$a1),
       y = 'Chronotype PRS') +
  theme_minimal(base_size = 9) + theme(plot.title = element_text(hjust = 0.5))

pdfOut <- file.path(OUT_DIR, paste0(PGS_ID, '_by_', SNP, '.pdf'))
pdf(pdfOut, width = 4.5, height = 4.2); print(pl); invisible(dev.off())
message('Wrote ', pdfOut)
