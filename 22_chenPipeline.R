# =============================================================================
# 22  Chen et al. (2025) rhyQTL pipeline applied to the 54 independent cd-eQTL
#
# Re-implements the rhyQTL calling of Chen et al. (github.com/YingChen10/rhyQTL,
# Analysis/0_rhyQTL_mapping, scripts 0_0 to 0_5) for each lead cd-eQTL, and
# reports at which step each locus would be lost. Steps, in their order:
#   0 genotype   (0_0)  at least 2 genotype groups with >= 50 donors
#   1 rhythm     (0_1)  cosinor (lm x ~ cos + sin, LRT vs intercept) in the two
#                       largest genotype groups; kept if either group has
#                       p < 0.01 and either has peak-to-trough amplitude
#                       2*sqrt(a^2 + b^2) > log2(1.5)
#   2 dryR once  (0_2)  largest group downsampled to the size of the second,
#                       dryR on the two groups; kept unless the chosen model is
#                       1 (flat in both) or 4 (shared rhythm)
#   3 dryR x20   (0_3)  the downsampling repeated 20 times; G-test of the chosen
#                       model counts against uniform; kept if p < 0.05 (p = 0
#                       when one model is always chosen)
#   4 HANOVA     (0_4)  ANOVA of cos/sin vs cos/sin + group:cos/sin on the two
#                       largest groups
#   5 combine    (0_5)  cosinor p of the group with the larger amplitude
#                       < 5e-4, G-test p < 0.05, HANOVA BH < 0.05
#
# Differences from Chen et al., forced by the data available here:
#   - Expression is GTEx v10 normalised expression residualised on our
#     covariates (as in 07 and 11), not their log2 CPM with covariates removed.
#     Their amplitude threshold log2(1.5) is on the log2 CPM scale, so it is
#     applied here on a different scale (step 1 reports amplitudes; the call
#     is also given without it, `rhyQTL_noAmpFilter`).
#   - HANOVA: their models have no group intercept, and norm = TRUE (each group
#     divided by its mean) is what removes the genotype main effect. That
#     division is unstable when group means are near 0, as for residuals, so
#     each group is mean-centred instead and HANOVA run with norm = FALSE.
#   - Time is the CHIRAL donor phase in hours (DIP), the same source as the
#     cd-eQTL scan. Chen read `GTEx_donor_time_science.txt`, whose name and
#     `CorrectedTOD` column suggest the Talamanca et al. (Science 2023) donor
#     times; not verified.
#   - BH in step 5 is over all of Chen's tests in the tissue, which cannot be
#     reproduced for 54 loci. The raw HANOVA p is reported; `rhyQTL` uses raw
#     p < 0.05, an upper bound on what Chen would call.
#   - Genotype groups are rounded dosages from our plink set (GTEx v9), not the
#     GTEx v8 VCF. Groups are ordered by size, as in Chen, so allele coding
#     does not matter.
#
# Run:     sbatch 22_chenPipeline.sh         (after 04; needs dryR)
# Output:  Results/published/chenPipeline/chenPipeline_54cdeQTL.csv
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)
library(dryR)

OUT_DIR <- paste0(CHEN_DIR, 'chenPipeline/')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

MIN_SIZE <- 50          # 0_0: donors per genotype group
PERIOD   <- 24
N_DRAWS  <- 20          # 0_3: downsampling repeats

# --- Chen et al. functions ------------------------------------------------------------

# 0_1 harm_reg, unchanged
harm_reg <- function(x, t, period) {
  fit0  <- lm(x ~ 1)
  c     <- cos(2 * pi * t / period)
  s     <- sin(2 * pi * t / period)
  fit1  <- lm(x ~ c + s)
  a     <- coef(fit1)[2]
  b     <- coef(fit1)[3]
  p.val <- lrtest(fit1, fit0)$Pr[2]
  amp   <- 2 * sqrt(a^2 + b^2)
  phase <- period * (atan2(b, a) %% (2 * pi)) / (2 * pi)
  c(pval = p.val, phase = unname(phase), amp = unname(amp))
}

# 0_4 HANOVA for one series per group (their function, reduced to one column)
HANOVA <- function(val1, val2, times1, times2, period, norm = TRUE) {
  if (norm) { val1 <- val1 / mean(val1, na.rm = TRUE); val2 <- val2 / mean(val2, na.rm = TRUE) }
  timevec <- c(times1, times2) / period * 2 * pi
  cosB <- cos(timevec); sinB <- sin(timevec)
  ind  <- c(rep(0, length(times1)), rep(1, length(times2)))
  vals <- c(val1, val2)
  f.r  <- lm(vals ~ cbind(cosB, sinB))
  f.f  <- lm(vals ~ cbind(cosB, sinB, cosB * ind, sinB * ind))
  anova(f.r, f.f)$`Pr(>F)`[2]
}

# dryR on two groups; Chen pass a vector, some dryR versions need a matrix
dryOnce <- function(x, group, time) {
  out <- tryCatch(drylm(x, group, time), error = function(e) NULL)
  if (is.null(out)) {
    m <- matrix(x, nrow = 1, dimnames = list('gene', paste0('s', seq_along(x))))
    out <- drylm(m, group, time)
  }
  p <- as.data.frame(out[['parameters']])[1, ]
  # By name where possible; otherwise Chen's positional layout
  # (mean, a, b, amp, relamp, phase per group, then chosen_model)
  mCol <- grep('^chosen_model$', names(p))
  aCol <- grep('^amp_', names(p))
  list(model = as.integer(p[[if (length(mCol)) mCol else 13]]),
       amp1  = as.numeric(p[[if (length(aCol) >= 2) aCol[1] else 4]]),
       amp2  = as.numeric(p[[if (length(aCol) >= 2) aCol[2] else 10]]))
}

# 0_2 / 0_3 random.sample.test: largest group downsampled to the second's size
randomSampleTest <- function(tmp, num) {
  keep <- sample(which(tmp$genotype == num$Var1[1]), size = num$Freq[2], replace = FALSE)
  r    <- rbind(tmp[keep, ], tmp[tmp$genotype == num$Var1[2], ])
  dryOnce(r$x, r$genotype, r$time)
}

# G-test of observed counts against a uniform distribution (DescTools::GTest default)
gTestUniform <- function(obs) {
  if (length(obs) == 1) return(0)
  e <- sum(obs) / length(obs)
  pchisq(2 * sum(obs * log(obs / e)), df = length(obs) - 1, lower.tail = FALSE)
}

# --- One locus ----------------------------------------------------------------------------
chenLocus <- function(d) {
  covTerms <- c(if ('sex' %in% colnames(d)) 'sex', 'platform', 'pcr',
                paste0('PC', 1:5), paste0('InferredCov', 1:10))
  covTerms <- covTerms[covTerms %in% colnames(d)]
  d$x <- residuals(lm(as.formula(paste('gene ~', paste(covTerms, collapse = '+'))), data = d))
  tmp <- data.frame(genotype = as.character(round(d$variant)),
                    time = d$radian / (2 * pi) * PERIOD, x = d$x)

  num <- as.data.frame(table(factor(tmp$genotype, levels = 0:2)), stringsAsFactors = FALSE)
  out <- list(n_0 = num$Freq[1], n_1 = num$Freq[2], n_2 = num$Freq[3])
  num <- num[order(-num$Freq), ]
  out$group_large <- num$Var1[1]; out$group_second <- num$Var1[2]

  # 0: genotype groups
  out$step0_size <- sum(num$Freq >= MIN_SIZE) >= 2
  if (num$Freq[2] < 3) return(out)   # nothing can be fitted

  # 1: cosinor in the two largest groups
  r1 <- harm_reg(tmp$x[tmp$genotype == num$Var1[1]], tmp$time[tmp$genotype == num$Var1[1]], PERIOD)
  r2 <- harm_reg(tmp$x[tmp$genotype == num$Var1[2]], tmp$time[tmp$genotype == num$Var1[2]], PERIOD)
  out[c('pval_0', 'phase_0', 'amp_0')] <- as.list(r1)
  out[c('pval_1', 'phase_1', 'amp_1')] <- as.list(r2)
  out$step1_pval <- min(r1['pval'], r2['pval']) < 0.01
  out$step1_amp  <- max(r1['amp'], r2['amp']) > log2(1.5)
  out$pval_maxAmp <- if (r1['amp'] >= r2['amp']) r1[['pval']] else r2[['pval']]

  # 2: one downsampled dryR
  tmp2 <- tmp[tmp$genotype %in% num$Var1[1:2], ]
  once <- tryCatch(randomSampleTest(tmp2, num), error = function(e) NULL)
  out$dryR_once_model <- if (is.null(once)) NA_integer_ else once$model
  out$step2_dryR <- !is.na(out$dryR_once_model) & !out$dryR_once_model %in% c(1, 4)

  # 3: 20 downsampled dryR, G-test
  draws <- lapply(seq_len(N_DRAWS), function(j) tryCatch(randomSampleTest(tmp2, num), error = function(e) NULL))
  models <- unlist(lapply(draws, function(z) z$model))
  if (length(models)) {
    obs <- table(models)
    out$gtest_p      <- gTestUniform(as.numeric(obs))
    out$gtest_model  <- as.integer(names(obs)[which.max(obs)])
    out$gtest_times  <- max(obs)
  }
  out$step3_gtest <- isTRUE(out$gtest_p < 0.05)

  # 4: HANOVA (norm = FALSE, see header)
  g1 <- tmp[tmp$genotype == num$Var1[1], ]; g2 <- tmp[tmp$genotype == num$Var1[2], ]
  # Their models have no group intercept; norm = TRUE (x / group mean) is what
  # removes the genotype main effect. On residuals we centre each group instead.
  out$hanova_p <- HANOVA(g1$x - mean(g1$x), g2$x - mean(g2$x), g1$time, g2$time, PERIOD, norm = FALSE)
  out
}

# --- Run ----------------------------------------------------------------------------------
tod    <- loadDIP()
plink  <- loadGenotypes()$X
cdeqtl <- readIndependentCdeQTL()
lead   <- fread(paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv')) %>%
  select(gene, variant, tissue, cdeQTL_p = p.value)
message('Loci: ', nrow(cdeqtl))

res <- vector('list', nrow(cdeqtl))
for (i in seq_len(nrow(cdeqtl))) {
  row <- cdeqtl[i, ]
  message(sprintf('[%2d/%d] %s %s %s', i, nrow(cdeqtl), row$tissue, row$hgnc_symbol, row$variant))
  fit <- tryCatch(singleTissueInteractionLRT(tod, plink, row$gene, row$variant, row$tissue),
                  error = function(e) { message('   refit failed: ', conditionMessage(e)); NULL })
  if (is.null(fit) || is.null(fit$data)) next
  r <- tryCatch(chenLocus(as.data.frame(fit$data)),
                error = function(e) { message('   Chen pipeline failed: ', conditionMessage(e)); NULL })
  if (!is.null(r)) res[[i]] <- as.data.frame(c(as.list(row), r))
}

out <- bind_rows(res) %>%
  left_join(lead, by = c('gene', 'variant', 'tissue')) %>%
  mutate(
    step5_rhythm = pval_maxAmp < 5e-4,
    rhyQTL = step0_size & step1_pval & step1_amp & step2_dryR & step3_gtest &
             step5_rhythm & hanova_p < 0.05,
    rhyQTL_noAmpFilter = step0_size & step1_pval & step2_dryR & step3_gtest &
             step5_rhythm & hanova_p < 0.05,
    first_failed = case_when(
      !step0_size         ~ '0 genotype groups < 50',
      !step1_pval         ~ '1 no group rhythmic at p < 0.01',
      !step1_amp          ~ '1 amplitude < log2(1.5)',
      !step2_dryR         ~ '2 dryR model 1 or 4',
      !step3_gtest        ~ '3 G-test p >= 0.05',
      !step5_rhythm       ~ '5 rhythm p of max-amplitude group >= 5e-4',
      hanova_p >= 0.05    ~ '4 HANOVA p >= 0.05',
      TRUE                ~ 'passes (HANOVA BH not applied)'))

fwrite(out, paste0(OUT_DIR, 'chenPipeline_54cdeQTL.csv'))
print(out %>% count(first_failed, sort = TRUE))
message('Called by the Chen pipeline: ', sum(out$rhyQTL, na.rm = TRUE), ' / ', nrow(out),
        ' (', sum(out$rhyQTL_noAmpFilter, na.rm = TRUE), ' without the amplitude filter)')
