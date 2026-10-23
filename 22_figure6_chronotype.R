# =============================================================================
# 22  Figure 6: cd-eQTL of the chronotype loci (PheWAS, 19)
#
#   a  CIART           ~ rs2147324   lead cd-eQTL, in its cd-eQTL tissue
#   b  PER3            ~ rs3135419   lead cd-eQTL, in its cd-eQTL tissue
#   c  ENSG00000236266 ~ rs697676    lead cd-eQTL, in its cd-eQTL tissue (Ovary)
#   d  PER3            ~ rs697676    same variant and tissue as c
# The tissues of a-c are read from Supplementary Table 3; d uses the tissue of c.
#
# Each panel: covariate-adjusted expression against DIP by genotype, with the
# interaction-model (M3) curve of each genotype and its 95% band. Models as in
# the scan (02):
#   M1  expression ~ dosage + covariates
#   M2  M1 + sin(DIP) + cos(DIP)
#   M3  M2 + dosage:sin(DIP) + dosage:cos(DIP)
#   cd-eQTL LRT = M3 vs M2 (2 df)
# Points are the M3 partial residuals (residual + genotype and rhythm terms).
# Dosage = copies of the plink A1 allele (as genio::read_plink); the counted
# allele of each variant is printed and written to the table.
# Only the three variants are read from the .bed (no loadGenotypes()).
#
# Run:     sbatch 22_figure6_chronotype.sh     (after 04; minutes)
# Outputs: Figures/Figure_6_chronotype.{jpeg,pdf}
#          Figures/Figure_6_chronotype_fits.csv   LRT p, n, amplitude and peak hour per genotype
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
options(bitmapType = 'cairo')

PANELS <- tibble::tribble(
  ~tag, ~symbol,           ~gene,             ~variant,
  'a',  'CIART',           'ENSG00000159208', 'rs2147324',
  'b',  'PER3',            'ENSG00000049246', 'rs3135419',
  'c',  'ENSG00000236266', 'ENSG00000236266', 'rs697676',
  'd',  'PER3',            'ENSG00000049246', 'rs697676')

# Genotype colours as in the rs697676 figure (met.brewer('Degas', n = 3))
GENO_COLS  <- c('0' = '#96410e', '1' = '#e5a335', '2' = '#556219')
COVARIATES <- c('sex', 'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
toHour <- function(r) (r %% (2 * pi)) / (2 * pi) * 24

# --- Tissues from Supplementary Table 3 ------------------------------------------
supTable3 <- fread(paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv')) %>%
  mutate(geneId = sub('\\..*', '', gene))
lead <- PANELS %>% filter(tag != 'd') %>%
  left_join(supTable3 %>% select(geneId, variant, tissue, p_table = p.value),
            by = c('gene' = 'geneId', 'variant'))
if (any(is.na(lead$tissue))) stop('Not in Supplementary Table 3: ',
                                  paste(lead$tag[is.na(lead$tissue)], collapse = ', '))
if (anyDuplicated(lead$tag)) {
  message('Lead cd-eQTL in several tissues; plotting the smallest p:')
  print(lead)
  lead <- lead %>% group_by(tag) %>% slice_min(p_table, n = 1, with_ties = FALSE) %>% ungroup()
}
PANELS <- PANELS %>%
  left_join(lead %>% select(tag, tissue, p_table), by = 'tag') %>%
  mutate(tissue = ifelse(tag == 'd', tissue[tag == 'c'], tissue))

# --- Genotypes: read only the requested variants from the plink .bed --------------
# SNP-major .bed: 3 magic bytes, then ceiling(n/4) bytes per variant, 2 bits per
# donor, low bits first. Codes 00/10/11/01 = 2/1/0/NA copies of the A1 allele.
readPlinkVariants <- function(prefix, ids) {
  bim <- fread(paste0(prefix, '.bim'), header = FALSE,
               col.names = c('chr', 'id', 'posg', 'pos', 'a1', 'a2'))
  fam <- fread(paste0(prefix, '.fam'), header = FALSE)
  n   <- nrow(fam)
  nb  <- ceiling(n / 4)
  con <- file(paste0(prefix, '.bed'), 'rb'); on.exit(close(con))
  stopifnot(identical(readBin(con, 'raw', 3), as.raw(c(0x6c, 0x1b, 0x01))))
  donors <- paste0('GTEX-', vapply(strsplit(fam$V2, '-'), `[`, '', 2))
  lut <- c(2L, NA_integer_, 1L, 0L)
  out <- list(); info <- list()
  for (v in ids) {
    idx <- which(bim$id == v)
    if (length(idx) != 1) stop(v, ': ', length(idx), ' rows in the .bim')
    seek(con, 3 + (idx - 1) * nb)
    bytes <- as.integer(readBin(con, 'raw', nb))
    codes <- as.vector(rbind(bytes %% 4L, (bytes %/% 4L) %% 4L,
                             (bytes %/% 16L) %% 4L, (bytes %/% 64L) %% 4L))[seq_len(n)]
    out[[v]]  <- lut[codes + 1L]
    info[[v]] <- bim[idx, .(variant = id, chr, pos, counted_allele = a1, other_allele = a2)]
  }
  geno <- as.data.table(out); geno[, SUBJID := donors]
  list(geno = geno[!duplicated(SUBJID)], info = rbindlist(info))
}

# --- Expression and covariates of one tissue ------------------------------------------
readTissue <- function(tissueName, genes) {
  expr <- fread(paste0(GTEX_EXPR_DIR, tissueName, '.v10.normalized_expression.bed.gz'))
  expr <- expr[sub('\\..*', '', gene_id) %in% genes]
  donorCols <- grep('GTEX', names(expr), value = TRUE)
  ex <- as.data.table(t(as.matrix(expr[, ..donorCols])))
  setnames(ex, sub('\\..*', '', expr$gene_id))
  ex[, SUBJID := donorCols]
  cov <- fread(paste0(GTEX_COV_DIR, tissueName, '.v10.covariates.txt'))
  cv  <- as.data.table(t(as.matrix(cov[, -1])))
  setnames(cv, cov$ID)
  cv[, SUBJID := names(cov)[-1]]
  merge(ex, cv, by = 'SUBJID')
}

# --- M1-M3, LRT, partial residuals and M3 curve per genotype ---------------------
fitPanel <- function(d, gene, variant) {
  d <- d[!is.na(d[[gene]]) & !is.na(d[[variant]]) & !is.na(radian)]
  covs <- COVARIATES[COVARIATES %in% names(d)]
  covs <- covs[vapply(covs, function(x) uniqueN(d[[x]], na.rm = TRUE) > 1, logical(1))]
  d <- copy(d); setnames(d, c(gene, variant), c('gene', 'g'))
  base <- paste('gene ~ g +', paste(covs, collapse = ' + '))
  m2 <- lm(as.formula(paste(base, '+ sin(radian) + cos(radian)')), data = d)
  m3 <- lm(as.formula(paste(base, '+ sin(radian) + cos(radian) + g:sin(radian) + g:cos(radian)')), data = d)
  b  <- coef(m3)
  gm <- mean(d$g)
  curveAt <- function(gg, x) b[['g']] * (gg - gm) +
    (b[['sin(radian)']] + gg * b[['g:sin(radian)']]) * sin(x) +
    (b[['cos(radian)']] + gg * b[['g:cos(radian)']]) * cos(x)
  grid <- seq(0, 2 * pi, length.out = 145)
  vt <- c('g', 'sin(radian)', 'cos(radian)', 'g:sin(radian)', 'g:cos(radian)')
  V  <- vcov(m3)[vt, vt]
  tq <- qt(0.975, df.residual(m3))
  genotypes <- sort(unique(d$g))
  fit <- data.table(n = nrow(d),
                    p_cdeqtl = pchisq(2 * as.numeric(logLik(m3) - logLik(m2)), df = 2, lower.tail = FALSE))
  for (gg in 0:2) {
    bs <- b[['sin(radian)']] + gg * b[['g:sin(radian)']]
    bc <- b[['cos(radian)']] + gg * b[['g:cos(radian)']]
    fit[, paste0('n_g', gg) := sum(d$g == gg)]
    fit[, paste0('amp_g', gg) := sqrt(bs^2 + bc^2)]
    fit[, paste0('peak_h_g', gg) := toHour(atan2(bs, bc))]
  }
  list(fit = fit,
       points = data.table(genotype = d$g, hour = toHour(d$radian),
                           expr_adj = residuals(m3) + curveAt(d$g, d$radian)),
       curves = rbindlist(lapply(genotypes, function(gg) {
         A  <- cbind(gg - gm, sin(grid), cos(grid), gg * sin(grid), gg * cos(grid))
         se <- sqrt(rowSums((A %*% V) * A))
         f  <- curveAt(gg, grid)
         data.table(genotype = gg, hour = grid / (2 * pi) * 24, fit = f,
                    lo = f - tq * se, hi = f + tq * se)
       })))
}

rhythmPlot <- function(fp, title, tag) {
  ggplot(fp$points, aes(hour, expr_adj, colour = factor(genotype))) +
    geom_point(size = 0.6, alpha = 0.5) +
    geom_ribbon(data = fp$curves, aes(hour, ymin = lo, ymax = hi, fill = factor(genotype)),
                inherit.aes = FALSE, alpha = 0.4, colour = NA) +
    geom_line(data = fp$curves, aes(hour, fit), linewidth = 0.9) +
    scale_x_continuous(breaks = seq(0, 24, 3), limits = c(0, 24)) +
    scale_colour_manual(values = GENO_COLS, name = 'Genotype', drop = FALSE) +
    scale_fill_manual(values = GENO_COLS, name = 'Genotype', drop = FALSE) +
    labs(title = title, tag = tag, x = 'Time of day (h)', y = 'Residualised expression') +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = 'bottom',
          plot.tag = element_text(face = 'bold', size = 12))
}

# --- Run ------------------------------------------------------------------------------
load(DIP_FILE)   # provides `phi`, radians
dip <- data.table(SUBJID = paste0('GTEX-', names(phi)), radian = as.numeric(phi))[!is.na(radian)]
gt  <- readPlinkVariants(GENOTYPE_PLINK, unique(PANELS$variant))
print(gt$info)

data <- lapply(setNames(nm = unique(PANELS$tissue)), function(tis)
  merge(merge(readTissue(tis, unique(PANELS$gene)), gt$geno, by = 'SUBJID'), dip, by = 'SUBJID'))

fits <- lapply(seq_len(nrow(PANELS)), function(k)
  with(PANELS[k, ], fitPanel(data[[tissue]], gene, variant)))
res <- bind_cols(PANELS, rbindlist(lapply(fits, `[[`, 'fit'))) %>%
  left_join(gt$info %>% select(variant, counted_allele, other_allele), by = 'variant')
message('cd-eQTL LRT here vs Supplementary Table 3 (should agree for a-c):')
print(res %>% select(tag, symbol, variant, tissue, n, p_cdeqtl, p_table), digits = 3)
fwrite(res, paste0(FIG_DIR, 'Figure_6_chronotype_fits.csv'))

plots <- lapply(seq_len(nrow(res)), function(k) {
  title <- sprintf('%s ~ %s in %s\ncd-eQTL LRT, p = %.2g', res$symbol[k], res$variant[k],
                   str_replace_all(res$tissue[k], '_', ' '), res$p_cdeqtl[k])
  rhythmPlot(fits[[k]], title, res$tag[k])
})
# Shared legend taken from the first panel
g   <- ggplotGrob(plots[[1]])
leg <- g$grobs[[which(grepl('^guide-box', g$layout$name) &
                      vapply(g$grobs, function(x) !inherits(x, 'zeroGrob'), logical(1)))[1]]]
plots <- lapply(plots, function(p) p + theme(legend.position = 'none'))

drawFigure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 2,
                                    heights = grid::unit(c(1, 1, 0.08), 'null'))))
  vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
  for (k in seq_along(plots)) print(plots[[k]], vp = vp((k + 1) %/% 2, 2 - k %% 2))
  grid::pushViewport(vp(3, 1:2)); grid::grid.draw(leg); grid::popViewport()
}
# --- Location of each variant (GENCODE v50, canonical transcripts) --------------------
# GTEx v10 quantified with GENCODE v39; v50 is used for annotation only.
GTF_FILE <- Sys.getenv('GTF_FILE', '/QRISdata/Q9564/GENCODE/v50/GRCh38/gencode.v50.annotation.gtf.gz')
annotateVariant <- function(chr, pos) {
  cmd <- sprintf("zcat '%s' | awk -F'\\t' '$1==\"%s\" && $4<=%d && $5>=%d'",
                 GTF_FILE, chr, as.integer(pos), as.integer(pos))
  g <- fread(cmd = cmd, header = FALSE, sep = '\t', quote = '',
             col.names = c('chr', 'src', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr'))
  if (nrow(g) == 0) return(data.table(feature = 'intergenic'))
  attrOf <- function(key) ifelse(grepl(paste0(key, ' "'), g$attr),
                                 sub(paste0('.*', key, ' "([^"]+)".*'), '\\1', g$attr), NA)
  g[, `:=`(gene_name = attrOf('gene_name'), gene_type = attrOf('gene_type'),
           exon_number = attrOf('exon_number'),
           canonical = grepl('Ensembl_canonical', attr))]
  g[feature == 'gene' | canonical, .(gene_name, gene_type, feature, exon_number, strand, start, end)]
}
nearestTSS <- function(chr, pos, genes) {
  cmd <- sprintf("zcat '%s' | awk -F'\\t' '$1==\"%s\" && $3==\"gene\"'", GTF_FILE, chr)
  g <- fread(cmd = cmd, header = FALSE, sep = '\t', quote = '',
             col.names = c('chr', 'src', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr'))
  g[, gene_id := sub('.*gene_id "([^".]+).*', '\\1', attr)]
  g <- g[gene_id %in% genes]
  g[, .(gene_id, strand, start, end, tss = ifelse(strand == '+', start, end),
        snp_minus_tss = pos - ifelse(strand == '+', start, end),
        snp_in_gene_body = pos >= start & pos <= end)]
}
if (file.exists(GTF_FILE)) {
  for (k in seq_len(nrow(gt$info))) {
    v <- gt$info[k]
    chr <- paste0('chr', sub('^chr', '', v$chr))
    message('\n', v$variant, ' (', chr, ':', v$pos, ') overlaps:')
    print(annotateVariant(chr, v$pos))
    print(nearestTSS(chr, v$pos, unique(PANELS$gene[PANELS$variant == v$variant])))
  }
} else message('GTF not found (', GTF_FILE, '): variant annotation skipped; set GTF_FILE')

stem <- paste0(FIG_DIR, 'Figure_6_chronotype')
pdf(paste0(stem, '.pdf'), width = 8, height = 7.5); drawFigure(); invisible(dev.off())
jpeg(paste0(stem, '.jpeg'), width = 8, height = 7.5, units = 'in', res = 300); drawFigure(); invisible(dev.off())
message('Wrote ', stem, '.{pdf,jpeg}')
