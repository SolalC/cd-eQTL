# =============================================================================
# Shared functions. Sourced by 00_config.R; relies on the paths defined there.
#
#   Data loading     loadDIP, loadGenotypes, selectIndividual, rhythmicTissues,
#                    readIndependentCdeQTL, addSymbol
#   Rhythmic genes   harmonicRegression, extractSNPs
#   cd-eQTL models   QTLregressionParallel, singleTissueInteractionLRT,
#                    multiTissueReg, calculate_rhythmicity_parameters,
#                    cosineSimilarity
#   Bingham et al.   residualisedGenotypeParameters, binghamGroupTest,
#   (1982) tests     binghamGenotypeParameters
# =============================================================================

# --- Data loading --------------------------------------------------------------

# Donor inferred phase (DIP, radians) keyed by GTEx subject ID.
loadDIP <- function() {
  load(DIP_FILE)   # provides `phi`
  data.frame(SUBJID = paste0('GTEX-', names(phi)),
             DIP    = phi,
             radian = phi)
}

# Full plink genotype set (8.25M variants; needs ~100 GB of memory).
loadGenotypes <- function() genio::read_plink(GENOTYPE_PLINK)

# Expression, gene positions and covariates of one tissue, restricted to donors
# with a DIP, expression and covariates.
selectIndividual <- function(tod, tissueName) {
  # Read the tissue of interest:
  tissuePath <- paste0(GTEX_EXPR_DIR, tissueName, '.v10.normalized_expression.bed.gz')
  tissue <- fread(tissuePath)
  tissueIndv <- colnames(tissue)[str_detect(colnames(tissue), 'GTEX')]
  # Gene position:
  pos <- colnames(tissue)[!str_detect(colnames(tissue), 'GTEX')]
  tissuePos <- tissue[,..pos]
  # Covariates:
  covPath <- paste0(GTEX_COV_DIR, tissueName, '.v10.covariates.txt')
  cov <- fread(covPath)
  covColname <- cov$ID
  covIndv <- colnames(cov)[-1]
  cov <- as.data.table(t(cov[,-1]))
  colnames(cov) <- covColname
  rownames(cov) <- covIndv
  # Harmonise the dataset to contain all the same individuals:
  indv <- Reduce(intersect, list(tod$SUBJID,
                                 tissueIndv,
                                 covIndv))
  # Expression:
  tissueExpr <- tissue[,..indv]
  tissueIndvFiltered <- colnames(tissueExpr)
  tissueExpr <- as.data.table(t(tissueExpr))
  colnames(tissueExpr) <- tissuePos$gene_id
  tissueExpr$SUBJID <- tissueIndvFiltered
  # Covariates:
  cov <- cov %>% mutate(SUBJID = rownames(.)) %>%
    filter(SUBJID %in% indv) %>%
    left_join(., tod, by = 'SUBJID')

  output <- list(tissue = tissueExpr, pos = tissuePos, cov = cov)
  return(output)
}

# Per-tissue harmonic LRT files written by 01_rhythmicGenes.R.
harmonicFiles <- function() {
  list.files(HARMONIC_DIR, pattern = '_LRT_harmonic_linear.csv$', full.names = TRUE)
}

# The 49 tissues with at least one rhythmic gene (Bonferroni within tissue < 0.05).
rhythmicTissues <- function() {
  out <- unlist(lapply(harmonicFiles(), function(f) {
    d <- fread(f)
    if (sum(p.adjust(d$p.value, 'bonferroni') < 0.05, na.rm = TRUE) > 0)
      str_remove(basename(f), '_LRT_harmonic_linear.csv')
  }))
  stopifnot(length(out) == 49)
  out
}

# The 54 independent cd-eQTL (Supplementary Table 3, from 04_cdeQTLtables.R).
# Row order matters: dryR objects (11) are stored in this order.
readIndependentCdeQTL <- function() {
  fread(paste0(TABLE_DIR, 'SupTable_3_cdEQTL_54independent.csv')) %>%
    distinct(gene, variant, tissue, hgnc_symbol)
}

# Add the gene symbol from 03_geneSymbolMap.py, falling back to the versioned
# Ensembl ID where Ensembl assigns no symbol.
addSymbol <- function(d) {
  map <- fread(SYMBOL_MAP)
  d %>%
    mutate(geneModif = sub('\\..*', '', gene)) %>%
    left_join(map, by = c('geneModif' = 'ensembl_gene_id')) %>%
    mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == '', gene, hgnc_symbol))
}

# --- Rhythmic genes --------------------------------------------------------------

# Cosinor vs linear model LRT (2 df) for every gene of one tissue, with the same
# 18 covariates in both models (17 in single-sex tissues).
harmonicRegression <- function(tissue, subset = 10, nCore = 1, sex = T) {
  expr <- tissue[[1]]
  pos <- tissue[[2]]
  cov <- tissue[[3]]
  # Get the genes for later iteration:
  genes <- colnames(expr)[str_detect(colnames(expr), 'ENSG')]
  if(subset < length(genes)) {
    genes <- genes[1:subset]
  }
  # Merge the data.frames:
  data <- expr %>% left_join(., cov, by = 'SUBJID')
  # Regression:
  if(sex) {
    regressionResults <- mclapply(genes, function(gene){
      designHarmonic <- paste0(
        gene, "~ sex+platform+pcr+PC1+PC2+PC3+PC4+PC5+InferredCov1+InferredCov2+InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8+InferredCov9+InferredCov10+sin(radian)+cos(radian)")
      designLinear <- paste0(
        gene, "~ sex+platform+pcr+PC1+PC2+PC3+PC4+PC5+InferredCov1+InferredCov2+InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8+InferredCov9+InferredCov10")
      lmHarmonic <- lm(data = data, formula = as.formula(designHarmonic))
      lmLinear <- lm(data = data, formula = as.formula(designLinear))
      lrt <- lmtest::lrtest(lmLinear, lmHarmonic)
      lrt <- broom::tidy(lrt)[2,]
      return(lrt)
    }, mc.cores = nCore)
  }
  if(!sex){
    regressionResults <- mclapply(genes, function(gene){
      designHarmonic <- paste0(
        gene, "~ platform+pcr+PC1+PC2+PC3+PC4+PC5+InferredCov1+InferredCov2+InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8+InferredCov9+InferredCov10+sin(radian)+cos(radian)")
      designLinear <- paste0(
        gene, "~ platform+pcr+PC1+PC2+PC3+PC4+PC5+InferredCov1+InferredCov2+InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8+InferredCov9+InferredCov10")
      lmHarmonic <- lm(data = data, formula = as.formula(designHarmonic))
      lmLinear <- lm(data = data, formula = as.formula(designLinear))
      lrt <- lmtest::lrtest(lmLinear, lmHarmonic)
      lrt <- broom::tidy(lrt)[2,]
      return(lrt)
    }, mc.cores = nCore)
  }

  regressionOutput <- do.call(rbind, regressionResults) %>%
    as.data.table() %>% mutate(gene = genes)

  return(regressionOutput)
}

# Variants within `window` bp of the TSS of each rhythmic gene.
extractSNPs <- function(tissue, regressionResults, lookup, window) {
  pos <- tissue[[2]]
  pos.f <- pos %>% filter(gene_id %in% regressionResults$gene) %>%
    mutate(start = start - window, end = end + window)

  lookupRes <- lookup[
    pos.f,
    on = .(chr = `#chr`, variant_pos2 >= start, variant_pos2 <= end),
    .(gene_id, chr, variant_pos, ref, alt, rs_id_dbSNP155_GRCh38p13,
      num_alt_per_site, variant_id, variant_id_b37),
    nomatch = 0L
  ]
  return(lookupRes)
}

# --- cd-eQTL models --------------------------------------------------------------

# The three nested models of the cd-eQTL scan for every gene-variant pair of one
# tissue, parallelised over genes (future_lapply):
#   1. gene ~ variant + covariates
#   2. gene ~ variant + covariates + sin + cos
#   3. gene ~ variant + covariates + sin + cos + variant:sin + variant:cos
# Variants whose heterozygote frequency is below frqThreshold are skipped.
QTLregressionParallel <- function(tissueName, tissue, plink, SNPs,
                                  frqThreshold = 0.1, sex = T) {
  expr <- tissue[[1]]
  pos <- tissue[[2]]
  cov <- tissue[[3]]
  # Read the variant that will be used for regression:
  genes <- SNPs %>% filter(!duplicated(gene_id)) %>% pull(gene_id)
  SNPs <- SNPs %>% filter(rs_id_dbSNP155_GRCh38p13 != '.',
                          rs_id_dbSNP155_GRCh38p13 %in% rownames(plink$X))
  plinkVariant <- which(rownames(plink$X) %in% SNPs$rs_id_dbSNP155_GRCh38p13)
  col <- str_split(colnames(plink$X), '-')
  colnames(plink$X) <- paste0('GTEX-', unlist(map(col, 2)))
  plinkIndv <- which(colnames(plink$X) %in% cov$SUBJID)
  plink.f <- plink$X[plinkVariant, plinkIndv]
  ID <- colnames(plink.f)
  variant <- rownames(plink.f)
  plink.f <- as.data.table(t(plink.f))
  plink.f$SUBJID <- ID
  # Remove duplicated columns:
  colDup <- duplicated(colnames(plink.f))
  plink.f <- plink.f[,!..colDup]
  # Merge the data.frames:
  data <- expr %>% left_join(., cov, by = 'SUBJID')
  data <- data %>% left_join(., plink.f, by = 'SUBJID')
  # Remove large plink objects no longer needed (prevents them being captured by the closure and
  # serialized to parallel workers by future_lapply):
  rm(plink, plink.f)
  # Consistent formula terms:
  if(sex) {
    all_terms <- c("sex", "platform", "pcr", paste0("PC", 1:5), paste0("InferredCov", 1:10))
  } else {
    all_terms <- c("platform", "pcr", paste0("PC", 1:5), paste0("InferredCov", 1:10))
  }
  designBasis <- paste(all_terms, collapse = "+")
  # Parallel function to process each gene:
  process_gene <- function(gene) {
    innerSNPs <- SNPs %>% filter(gene_id == gene)
    variants <- innerSNPs %>% pull(rs_id_dbSNP155_GRCh38p13)

    linear_results <- data.frame()
    harmonic_results <- data.frame()
    harmonic_interaction_results <- data.frame()
    lrt_results <- data.frame()
    for(variant in variants){
      freq <- data %>% dplyr::select(variant) %>% table()
      freq <- (freq/sum(freq))
      if(is.na(freq[2]) || freq[2] < frqThreshold) {
        next
      }
      designLinear <- paste0(gene, '~', variant,'+', designBasis)
      designHarmonic <- paste0(gene, '~', variant,'+', designBasis, '+sin(radian)+cos(radian)')
      designHarmonicInteraction <- paste0(gene, '~', variant,'+', designBasis, '+sin(radian)+cos(radian)',
                                          '+', variant, ':sin(radian)+', variant, ':cos(radian)')
      # Perform each model:
      lmLinear <- lm(data = data, formula = as.formula(designLinear))
      lmHarmonic <- lm(data = data, formula = as.formula(designHarmonic))
      lmHarmonicInteraction <- lm(data = data, formula = as.formula(designHarmonicInteraction))
      # Check the fit of each model:
      lrtVariant <- lmtest::lrtest(lmLinear, lmHarmonic, lmHarmonicInteraction)
      # Output:
      linear_results <- rbind(linear_results, (lmLinear %>%
                                broom::tidy() %>%
                                mutate(gene = gene, variant = variant)))
      harmonic_results <- rbind(harmonic_results, lmHarmonic %>%
                                  broom::tidy() %>%
                                  mutate(gene = gene, variant = variant))
      harmonic_interaction_results <- rbind(harmonic_interaction_results,
                                            lmHarmonicInteraction %>%
                                              broom::tidy() %>%
                                              mutate(gene = gene, variant = variant))
      lrt_results <- rbind(lrt_results, lrtVariant %>%
                             broom::tidy() %>%
                             mutate(gene = gene, variant = variant))
    }
    return(list(
      linear = linear_results,
      harmonic = harmonic_results,
      harmonicInteraction = harmonic_interaction_results,
      lrt = lrt_results
    ))
  }
  # Run the parallel processing
  results <- future_lapply(genes, process_gene)
  # Filter out NULL results (failed frequency threshold)
  results <- results[!sapply(results, is.null)]
  # Combine results
  regLinear <- bind_rows(lapply(results, function(x) x$linear))
  regHarmonic <- bind_rows(lapply(results, function(x) x$harmonic))
  regHarmonicInteraction <- bind_rows(lapply(results, function(x) x$harmonicInteraction))
  lrt <- bind_rows(lapply(results, function(x) x$lrt))
  return(list(linear = regLinear,
              harmonic = regHarmonic,
              harmonicInteraction = regHarmonicInteraction,
              lrt = lrt))
}

# The same three-model LRT for a single gene / variant / tissue. Optionally
# appends the results to <lrtOutDir>/<tissue>_LRT.csv and
# <coefOutDir>/<tissue>_harmonicInteraction.csv (same format as the scan).
# Returns lrt, coef, the fitted interaction model, nObs and the model data.
singleTissueInteractionLRT <- function(tod, plink, gene, variant, tissueName,
                                       lrtOutDir  = NULL,
                                       coefOutDir = NULL) {
  # --- Load expression -----------------------------------------------------
  tissuePath <- paste0(GTEX_EXPR_DIR, tissueName, '.v10.normalized_expression.bed.gz')
  tissue <- fread(tissuePath)
  tissue <- tissue %>% filter(sub('\\..*', '', gene_id) %in% sub('\\..*', '', gene))
  if (nrow(tissue) == 0) {
    message(tissueName, ': gene not found, skipping')
    return(NULL)
  }
  tissueIndv <- colnames(tissue)[str_detect(colnames(tissue), 'GTEX')]
  pos        <- colnames(tissue)[!str_detect(colnames(tissue), 'GTEX')]
  tissuePos  <- tissue[, ..pos]

  # --- Load covariates -----------------------------------------------------
  covPath <- paste0(GTEX_COV_DIR, tissueName, '.v10.covariates.txt')
  cov        <- fread(covPath)
  covColname <- cov$ID
  covIndv    <- colnames(cov)[-1]
  cov        <- as.data.table(t(cov[, -1]))
  colnames(cov) <- covColname
  rownames(cov) <- covIndv

  # --- Harmonise individuals -----------------------------------------------
  indv <- Reduce(intersect, list(tod$SUBJID, tissueIndv, covIndv))

  tissueExpr              <- tissue[, ..indv]
  tissueIndvFiltered      <- colnames(tissueExpr)
  tissueExpr              <- as.data.table(t(tissueExpr))
  tissueExpr$SUBJID       <- tissueIndvFiltered
  colnames(tissueExpr)    <- c('gene', 'SUBJID')

  cov <- cov %>%
    mutate(SUBJID = rownames(.)) %>%
    filter(SUBJID %in% indv) %>%
    left_join(., tod, by = 'SUBJID')

  # --- Genotype ------------------------------------------------------------
  col           <- str_split(colnames(plink), '-')
  colnames(plink) <- paste0('GTEX-', unlist(map(col, 2)))
  variantDf     <- plink[variant, indv] %>% as.data.frame()
  variantDf$SUBJID <- colnames(plink[, indv])
  colnames(variantDf) <- c('variant', 'SUBJID')

  # --- Merge ---------------------------------------------------------------
  data <- cov %>%
    left_join(., variantDf,  by = 'SUBJID') %>%
    left_join(., tissueExpr, by = 'SUBJID') %>%
    filter(!is.na(variant))

  if (nrow(data) == 0) {
    message(tissueName, ': no overlapping individuals, skipping')
    return(NULL)
  }

  # --- Formula terms (drop sex for sex-specific tissues) -------------------
  all_terms <- if ('sex' %in% colnames(data)) {
    c('sex', 'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  } else {
    c('platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  }
  designBasis <- paste(all_terms, collapse = '+')

  designLinear              <- paste0('gene ~ variant +', designBasis)
  designHarmonic            <- paste0('gene ~ variant +', designBasis,
                                      '+sin(radian)+cos(radian)')
  designHarmonicInteraction <- paste0('gene ~ variant +', designBasis,
                                      '+sin(radian)+cos(radian)',
                                      '+variant:sin(radian)+variant:cos(radian)')

  # --- Three-model LRT -----------------------------------------------------
  lmLinear              <- lm(data = data, formula = as.formula(designLinear))
  lmHarmonic            <- lm(data = data, formula = as.formula(designHarmonic))
  lmHarmonicInteraction <- lm(data = data, formula = as.formula(designHarmonicInteraction))

  lrtResult  <- lmtest::lrtest(lmLinear, lmHarmonic, lmHarmonicInteraction) %>%
    broom::tidy() %>%
    mutate(gene = gene, variant = variant, tissue = tissueName)

  coefResult <- lmHarmonicInteraction %>%
    broom::tidy() %>%
    mutate(gene = gene, variant = variant)

  # --- Write to disk -------------------------------------------------------
  if (!is.null(lrtOutDir)) {
    outFile <- file.path(lrtOutDir, paste0(tissueName, '_LRT.csv'))
    # Append so one file per tissue accumulates all gene/variant rows
    fwrite(lrtResult, outFile, append = file.exists(outFile))
  }
  if (!is.null(coefOutDir)) {
    outFile <- file.path(coefOutDir, paste0(tissueName, '_harmonicInteraction.csv'))
    fwrite(coefResult, outFile, append = file.exists(outFile))
  }

  return(list(lrt = lrtResult, coef = coefResult,
              fit = lmHarmonicInteraction, nObs = nrow(data), data = data))
}

# Interaction model for one gene-variant pair in every tissue expressing the
# gene. Returns the tidy coefficients, one block per tissue.
multiTissueReg <- function(tod, plink, gene, variant, subset=NULL) {
  # Get the list of tissue:
  tissuePaths <- list.files(GTEX_EXPR_DIR, full.names = T)
  tissuePaths <- tissuePaths[!str_detect(tissuePaths, 'tbi')]

  if(!is.null(subset)) {
    tissuePaths <- tissuePaths[str_detect(tissuePaths, str_c(subset, collapse = "|"))]
    if(length(tissuePaths) == 0) stop('Wrong Subset')
  }

  # Get the list of covariates:
  covPaths <- list.files(GTEX_COV_DIR, full.names = T)
  regHarmonicInteraction <- data.frame()
  for(i in tissuePaths){
    # read tissue:
    tissue <- fread(i)
    tissue <- tissue %>% filter(gene_id %in% gene)
    tissueName <- str_remove(basename(i),
    '.v10.normalized_expression.bed.gz')
    if(nrow(tissue) == 0){
      print(paste0(tissueName, ': Skipped'))
      next}

    tissueIndv <- colnames(tissue)[str_detect(colnames(tissue), 'GTEX')]
    # covariate:
    cov <- fread(covPaths[str_detect(covPaths, tissueName)])
    covColname <- cov$ID
    covIndv <- colnames(cov)[-1]
    cov <- as.data.table(t(cov[,-1]))
    colnames(cov) <- covColname
    rownames(cov) <- covIndv
    # Harmonise the tissue to contain the same individuals:
    indv <- Reduce(intersect, list(tod$SUBJID,
                                   tissueIndv,
                                   covIndv))
    # Get the same individual for each of the metrics:
    tissueExpr <- tissue[,..indv]
    tissueIndvFiltered <- colnames(tissueExpr)
    tissueExpr <- as.data.table(t(tissueExpr))
    tissueExpr$SUBJID <- tissueIndvFiltered
    colnames(tissueExpr) <- c('gene', 'SUBJID')
    # Covariates:
    cov <- cov %>% mutate(SUBJID = rownames(.)) %>%
      filter(SUBJID %in% indv) %>%
      left_join(., tod, by = 'SUBJID')
    # Genotype:
    col <- str_split(colnames(plink), '-')
    colnames(plink) <- paste0('GTEX-', unlist(map(col, 2)))
    plinkIndv <- colnames(plink[,indv])
    variantDf <- plink[variant,indv] %>% as.data.frame()
    variantDf$SUBJID <- plinkIndv
    colnames(variantDf) <- c('variant', 'SUBJID')

    cov <- left_join(cov, variantDf, by = 'SUBJID')
    cov <- left_join(cov, tissueExpr, by = 'SUBJID')
    cov <- cov %>% filter(!is.na(variant))
    # Perform the regression:
    all_terms <- c("sex", "platform", "pcr", paste0("PC", 1:5),
                   paste0("InferredCov", 1:10))
    if (!'sex' %in% colnames(cov)) {
      all_terms <- c("platform", "pcr", paste0("PC", 1:5),
                     paste0("InferredCov", 1:10))
      }
    designBasis <- paste(all_terms, collapse = "+")
    designHarmonicInteraction <- paste0('gene ~ variant +',
                                        designBasis, '+sin(radian)+cos(radian)',
                                        '+ variant:sin(radian)+',
                                        'variant:cos(radian)')
    lmHarmonicInteraction <- lm(data = cov,
                                formula = as.formula(designHarmonicInteraction))
    # Output:
    regHarmonicInteraction <- rbind(regHarmonicInteraction,
                                    (lmHarmonicInteraction %>%
                                       broom::tidy() %>%
                                       mutate(gene = gene, variant = variant,
                                              tissue = tissueName)))
    print(paste0(tissueName, ': Done'))
  }
  return(regHarmonicInteraction)
}

# Per-genotype amplitude and acrophase from the interaction model coefficients
# (multiTissueReg output): b(g) = b + g * b_interaction, for g = 0, 1, 2.
calculate_rhythmicity_parameters <- function(df) {
  # Group by gene, variant, and tissue
  result_list <- lapply(split(df, list(df$gene, df$variant, df$tissue)), function(group_df) {
    coef_list <- setNames(group_df$estimate, group_df$term)
    sin_coef_base <- coef_list["sin(radian)"]
    cos_coef_base <- coef_list["cos(radian)"]
    sin_interaction <- coef_list["variant:sin(radian)"]
    cos_interaction <- coef_list["variant:cos(radian)"]

    genotypes <- c(0, 1, 2)
    results <- data.frame(
      gene = unique(group_df$gene),
      variant = unique(group_df$variant),
      tissue = unique(group_df$tissue),
      genotype = genotypes,
      sin_coef = sin_coef_base + genotypes * sin_interaction,
      cos_coef = cos_coef_base + genotypes * cos_interaction
    )
    results$amplitude <- sqrt(results$sin_coef^2 + results$cos_coef^2)
    # True acrophase (peak) in DIP/radian units: the fitted curve
    # sin_coef*sin(x) + cos_coef*cos(x) = A*cos(x - phi) peaks at
    # phi = atan2(sin_coef, cos_coef).
    results$phase <- atan2(results$sin_coef, results$cos_coef)
    return(results)
  })
  rhythmicity_parameters <- do.call(rbind, result_list)
  return(rhythmicity_parameters)
}

# Tissue x tissue cosine similarity of the four rhythm coefficients
# (sin, cos, variant:sin, variant:cos), restricted to the tissues given.
cosineSimilarity <- function(multiRegOutput, tissue) {
  require(lsa)
  # .env: the argument, not the column of the same name
  multiRegOutput <- multiRegOutput %>% filter(tissue %in% .env$tissue)
  # Extract rhythmicity coefficients:
  rhythmCoef <- multiRegOutput %>% filter(str_detect(term, '^sin|^cos')) %>%
    select(term, estimate, std.error, tissue) %>%
    mutate(parameter = ifelse(term == "sin(radian)", "sin", "cos")) %>%
    pivot_wider(
      id_cols = tissue,
      names_from = parameter,
      values_from = c(estimate, std.error),
      names_glue = "{parameter}{ifelse(.value == 'std.error', '.std.error', '')}"
    )
  # Extract interaction coefficients:
  rhythmInter <- multiRegOutput %>% filter(str_detect(term, ':sin|:cos')) %>%
    select(term, estimate, std.error, tissue) %>%
    mutate(parameter = ifelse(term == "variant:sin(radian)", "sinInteraction", "cosInteraction")) %>%
    pivot_wider(
      id_cols = tissue,
      names_from = parameter,
      values_from = c(estimate, std.error),
      names_glue = "{parameter}{ifelse(.value == 'std.error', '.std.error', '')}"
    )
  rhythmInter <- left_join(rhythmInter, rhythmCoef, by = 'tissue')

  interMatrix  <- matrix(multiRegOutput %>%
                           filter(str_detect(term, 'sin|cos')) %>%
                           pull(estimate),
                         nrow=4)
  colnames(interMatrix) <- multiRegOutput %>% pull(tissue) %>% unique()
  cosineMatrix <- (cosine(interMatrix))
  return(list(cosineMatrix, rhythmInter))
}

# --- Bingham et al. (1982) cosinor parameter tests ---------------------------------

# Residualise expression on the covariates once (all samples), then fit an
# independent single cosinor to each genotype group with >= minGroupN donors.
# Returns one row per group: n_g, residual df and RSS, the (sin, cos)
# coefficients and their 2x2 covariance. `data` is singleTissueInteractionLRT()$data.
residualisedGenotypeParameters <- function(data, minGroupN = MIN_GROUP_N) {
  covTerms <- c(if ('sex' %in% colnames(data)) 'sex',
                'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  covTerms <- covTerms[covTerms %in% colnames(data)]
  data$resid    <- residuals(lm(as.formula(paste('gene ~', paste(covTerms, collapse = '+'))),
                                data = data))
  data$genotype <- round(data$variant)
  rows <- lapply(sort(unique(data$genotype)), function(g) {
    d <- data[data$genotype == g, ]
    if (nrow(d) < minGroupN) return(NULL)
    fit <- lm(resid ~ cos(radian) + sin(radian), data = d)
    cf  <- coef(fit); V <- vcov(fit)
    if (any(is.na(cf[c('cos(radian)', 'sin(radian)')]))) return(NULL)
    data.frame(genotype = g, n_g = nrow(d), df = df.residual(fit),
               rss = sum(residuals(fit)^2),
               sin_coef = cf[['sin(radian)']], cos_coef = cf[['cos(radian)']],
               cov_ss = V['sin(radian)', 'sin(radian)'],
               cov_sc = V['sin(radian)', 'cos(radian)'],
               cov_cc = V['cos(radian)', 'cos(radian)'])
  })
  bind_rows(rows[!vapply(rows, is.null, logical(1))])
}

# Exact Bingham et al. (1982, Chronobiologia 9:397) comparison of single
# cosinors across the genotype groups, on covariate-residualised expression.
# Each group contributes a 3-parameter cosinor (M, beta, gamma), matching the
# paper's N - 3k residual degrees of freedom (eq. 48).
#
# Notation follows the paper (eq. 22): per group i,
#   beta_i  = cos coefficient  (=  A_i cos phi_i)
#   gamma_i = sin coefficient  (= -A_i sin phi_i)
#   A_i     = sqrt(beta_i^2 + gamma_i^2);   phi_i = atan2(-gamma_i, beta_i)
#   c22, c23, c33  = (X'X)^-1 elements of the (cos, sin) block
#   c22(phi), c33(phi) = that block rotated into the amplitude / acrophase
#                        directions (eq. 36); Var(A_i) = sigma^2 c22(phi_i)
#
# Returns the amplitude test (F2, eq. 49) and acrophase test (F3, eq. 50), both
# on F(k-1, N-3k). Groups smaller than minGroupN are dropped; with fewer than
# two groups left the tests are NA.
binghamGroupTest <- function(data, minGroupN = MIN_GROUP_N) {
  covTerms <- c(if ('sex' %in% colnames(data)) 'sex',
                'platform', 'pcr', paste0('PC', 1:5), paste0('InferredCov', 1:10))
  covTerms <- covTerms[covTerms %in% colnames(data)]

  naRow <- function(k) data.frame(
    k = k, df_denom = NA_real_, F_amplitude = NA_real_, p_amplitude_bingham = NA_real_,
    F_acrophase = NA_real_, p_acrophase_bingham = NA_real_)

  # 1. Residualise expression on the covariates (once, all samples).
  covFormula <- as.formula(paste('gene ~', paste(covTerms, collapse = '+')))
  data$resid <- residuals(lm(covFormula, data = data))

  # 2. One single cosinor per genotype group on the residuals.
  data$genotype <- round(data$variant)
  per <- lapply(sort(unique(data$genotype)), function(g) {
    d <- data[data$genotype == g, ]
    if (nrow(d) < minGroupN) return(NULL)
    fit <- lm(resid ~ cos(radian) + sin(radian), data = d)
    cf  <- coef(fit)
    if (any(is.na(cf[c('cos(radian)', 'sin(radian)')]))) return(NULL)
    beta_i  <- cf[['cos(radian)']]    #  A cos phi
    gamma_i <- cf[['sin(radian)']]    # -A sin phi
    s2_i    <- sum(residuals(fit)^2) / df.residual(fit)
    Cblock  <- vcov(fit)[c('cos(radian)', 'sin(radian)'),
                         c('cos(radian)', 'sin(radian)')] / s2_i   # (X'X)^-1 block
    list(A   = sqrt(beta_i^2 + gamma_i^2),
         phi = atan2(-gamma_i, beta_i),
         c22 = Cblock[1, 1], c23 = Cblock[1, 2], c33 = Cblock[2, 2],
         rss = sum(residuals(fit)^2), df = df.residual(fit))
  })
  per <- per[!vapply(per, is.null, logical(1))]
  k <- length(per)
  if (k < 2) return(naRow(k))

  A   <- vapply(per, `[[`, numeric(1), 'A')
  phi <- vapply(per, `[[`, numeric(1), 'phi')
  # Rotate each group's (cos, sin) covariance into amplitude / acrophase
  # directions at its own fitted acrophase (eq. 36).
  c22p <- vapply(per, function(p) p$c22 * cos(p$phi)^2 - 2 * p$c23 * cos(p$phi) * sin(p$phi) +
                                  p$c33 * sin(p$phi)^2, numeric(1))
  c33p <- vapply(per, function(p) p$c22 * sin(p$phi)^2 + 2 * p$c23 * cos(p$phi) * sin(p$phi) +
                                  p$c33 * cos(p$phi)^2, numeric(1))
  if (any(!is.finite(c22p)) || any(c22p <= 0) ||
      any(!is.finite(c33p)) || any(c33p <= 0)) return(naRow(k))

  W <- 1 / c22p           # eq. 49 weights
  V <- A^2 / c33p         # eq. 50 weights

  s2       <- sum(vapply(per, `[[`, numeric(1), 'rss')) /
              sum(vapply(per, `[[`, numeric(1), 'df'))   # pooled sigma^2 (eq. 48)
  df_denom <- sum(vapply(per, `[[`, numeric(1), 'df'))   # = N - 3k

  # Amplitude test (H2, eq. 49)
  F2 <- (sum(W * A^2) - sum(W * A)^2 / sum(W)) / ((k - 1) * s2)
  p2 <- pf(F2, k - 1, df_denom, lower.tail = FALSE)

  # Acrophase test (H3, eq. 50): pooled acrophase phiTilde then weighted spread
  phiTilde <- atan2(sum(V * sin(2 * phi)), sum(V * cos(2 * phi))) / 2
  F3 <- sum(V * sin(phi - phiTilde)^2) / ((k - 1) * s2)
  p3 <- pf(F3, k - 1, df_denom, lower.tail = FALSE)

  data.frame(k = k, df_denom = df_denom,
             F_amplitude = F2, p_amplitude_bingham = p2,
             F_acrophase = F3, p_acrophase_bingham = p3)
}

# Per-genotype amplitude and acrophase from the interaction model, with the 2x2
# covariance of b(g) = (beta_s + g*gamma_s, beta_c + g*gamma_c) (delta method).
binghamGenotypeParameters <- function(fit, genotypes = c(0, 1, 2)) {
  rhythmTerms <- c('sin(radian)', 'cos(radian)',
                   'variant:sin(radian)', 'variant:cos(radian)')
  cf <- coef(fit)
  V  <- vcov(fit)[rhythmTerms, rhythmTerms]
  beta_s  <- cf[['sin(radian)']];          beta_c  <- cf[['cos(radian)']]
  gamma_s <- cf[['variant:sin(radian)']];  gamma_c <- cf[['variant:cos(radian)']]

  do.call(rbind, lapply(genotypes, function(g) {
    # Jacobian of b(g) wrt (beta_s, beta_c, gamma_s, gamma_c).
    J <- matrix(c(1, 0, g, 0,
                  0, 1, 0, g), nrow = 2, byrow = TRUE)
    Sigma <- J %*% V %*% t(J)                       # 2x2 covariance of b(g)
    b_s <- beta_s + g * gamma_s
    b_c <- beta_c + g * gamma_c
    amp <- sqrt(b_s^2 + b_c^2)
    gA <- c(b_s, b_c) / amp
    var_amp <- as.numeric(t(gA) %*% Sigma %*% gA)
    denom <- b_s^2 + b_c^2
    gP <- c(b_c / denom, -b_s / denom)
    var_phase <- as.numeric(t(gP) %*% Sigma %*% gP)
    data.frame(
      genotype   = g,
      sin_coef   = b_s,
      cos_coef   = b_c,
      amplitude  = amp,
      amplitude_se = sqrt(var_amp),
      phase      = atan2(-b_c, b_s),
      phase_se   = sqrt(var_phase),
      cov_ss     = Sigma[1, 1],
      cov_sc     = Sigma[1, 2],
      cov_cc     = Sigma[2, 2],
      stringsAsFactors = FALSE
    )
  }))
}
