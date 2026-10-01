# =============================================================================
# 11  dryR validation of the 54 independent cd-eQTL
#
# dryR (Weger et al. 2021) cannot take covariates, so each locus is fitted on
# covariate-residualised expression (as for the Bingham tests, 07), with the
# genotype group (>= 10 donors) as condition and DIP in hours as time. dryR fits
# every candidate rhythm model (5 for two groups, 15 for three) and selects one
# by BIC weight. A locus fails validation when the selected model is the one in
# which every genotype group shares one rhythm (model 4 of 5, model 11 of 15).
#
# Run:     sbatch 11_dryRvalidation.sh      (after 04; needs the dryR package)
# Outputs: Results/inferred/rhythmicityParameters/dryR/dryR_raw.rds
#              the dryR object per locus, in Supplementary Table 3 row order
#          Results/inferred/rhythmicityParameters/dryR/dryR_perLocus.csv
#              chosen model, its description and BIC weight per locus
# Used by: Supplementary Figure 1 (12), Supplementary Table 5 (13)
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(genio)
library(broom)
library(lmtest)
library(dryR)

tod    <- loadDIP()
plink  <- loadGenotypes()$X
cdeqtl <- readIndependentCdeQTL()
message('Loci to validate: ', nrow(cdeqtl))

# dryR call, tolerant of the argument naming across dryR versions
runDry <- function(mat, group, time) {
  attempts <- list(
    function() drylm(mat, group, time),
    function() drylm(countData = mat, group = group, time = time),
    function() drylm(data = mat, group = group, time = time)
  )
  for (k in seq_along(attempts)) {
    out <- tryCatch(attempts[[k]](), error = function(e) e)
    if (!inherits(out, 'error')) return(out)
    message('   dryR call style ', k, ' failed: ', conditionMessage(out))
  }
  NULL
}

# Chosen model and its BIC weight from the dryR output
chosenModel <- function(dry) {
  if (is.null(dry)) return(list(model = NA_character_, weight = NA_real_))
  cand <- NULL
  for (nm in c('parameters', 'results', 'ncond_models', 'models')) {
    if (!is.null(dry[[nm]]) && is.data.frame(dry[[nm]])) { cand <- dry[[nm]]; break }
  }
  if (is.null(cand)) return(list(model = NA_character_, weight = NA_real_))
  mcol <- grep('chosen_model$|^chosen_model', colnames(cand), value = TRUE)
  wcol <- grep('BICW|weight',                 colnames(cand), value = TRUE)
  list(model  = if (length(mcol)) as.character(cand[[mcol[1]]][1]) else NA_character_,
       weight = if (length(wcol)) suppressWarnings(as.numeric(cand[[wcol[1]]][1])) else NA_real_)
}

# dryR model semantics, indexed by number of conditions (dryR:::create_matrix_list).
# Groups written together share one rhythm; a group with no cosinor term is flat.
modelDesc <- list(
  `2` = c('flat in both', 'rhythmic in 0 only', 'rhythmic in 1 only',
          'shared rhythm 0=1', 'separate rhythms 0,1'),
  `3` = c('flat in all', 'rhythmic in 0 only', 'rhythmic in 1 only', 'rhythmic in 2 only',
          '0=1 shared, 2 flat', '0,1 separate, 2 flat', '0=2 shared, 1 flat',
          '0,2 separate, 1 flat', '1=2 shared, 0 flat', '1,2 separate, 0 flat',
          'all three shared', '0 alone, 1=2 shared', '0=2 shared, 1 alone',
          '0=1 shared, 2 alone', 'all three separate'))

resList <- vector('list', nrow(cdeqtl))
rawList <- vector('list', nrow(cdeqtl))
for (i in seq_len(nrow(cdeqtl))) {
  row <- cdeqtl[i, ]
  message(sprintf('[%2d/%d] %s %s %s', i, nrow(cdeqtl), row$tissue, row$hgnc_symbol, row$variant))
  fit <- tryCatch(singleTissueInteractionLRT(tod, plink, row$gene, row$variant, row$tissue),
                  error = function(e) { message('   refit failed: ', conditionMessage(e)); NULL })
  if (is.null(fit) || is.null(fit$data)) next
  d <- fit$data

  # Same covariate residualisation as the Bingham classification.
  covTerms <- c(if ('sex' %in% colnames(d)) 'sex',
                'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  covTerms <- covTerms[covTerms %in% colnames(d)]
  d$resid    <- residuals(lm(as.formula(paste('gene ~', paste(covTerms, collapse = '+'))), data = d))
  d$genotype <- round(d$variant)
  d <- d[d$genotype %in% as.numeric(names(which(table(d$genotype) >= MIN_GROUP_N))), ]
  if (length(unique(d$genotype)) < 2) { message('   fewer than two usable genotype groups'); next }

  mat   <- matrix(d$resid, nrow = 1, dimnames = list(row$gene, paste0('s', seq_len(nrow(d)))))
  group <- paste0('g', d$genotype)
  time  <- as.numeric(d$radian) * 24 / (2 * pi)   # radians -> hours

  dry <- runDry(mat, group, time)
  if (is.null(dry)) {
    # drylm drops a one-row matrix to a vector; repeating the gene leaves the
    # per-gene fit unchanged.
    mat2 <- rbind(mat, mat, mat)
    rownames(mat2) <- paste0(row$gene, '_', 1:3)
    dry <- runDry(mat2, group, time)
  }
  rawList[[i]] <- dry
  cm <- chosenModel(dry)
  nGroups <- length(unique(d$genotype))
  model   <- as.integer(cm$model)
  resList[[i]] <- data.frame(
    gene = row$gene, variant = row$variant, tissue = row$tissue, hgnc_symbol = row$hgnc_symbol,
    n = nrow(d), n_groups = nGroups,
    genotype_groups   = paste(sort(unique(d$genotype)), collapse = '/'),
    dryR_chosen_model = model,
    dryR_model_desc   = if (is.na(model)) NA_character_ else modelDesc[[as.character(nGroups)]][model],
    dryR_model_BICW   = cm$weight)
  message('   dryR chosen model: ', model)
}

res <- bind_rows(resList)
fwrite(res, paste0(DRYR_DIR, 'dryR_perLocus.csv'))
saveRDS(rawList, paste0(DRYR_DIR, 'dryR_raw.rds'))
sharedModel <- c(`2` = 4L, `3` = 11L)
message('Loci selecting a shared rhythm (not validated): ',
        sum(res$dryR_chosen_model == sharedModel[as.character(res$n_groups)], na.rm = TRUE),
        ' / ', nrow(res))
