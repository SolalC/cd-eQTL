# =============================================================================
# 03  Ensembl gene ID -> gene symbol map
#
# Resolves a symbol for every rhythmic gene (Bonferroni-adjusted harmonic LRT
# p < 0.05 within tissue, from 01), which covers every gene in the
# supplementary tables, and caches it. Uses the Ensembl REST lookup endpoint
# (biomaRt is not reachable from the compute nodes). Where Ensembl assigns no
# symbol, addSymbol() falls back to the versioned Ensembl ID.
#
# Run:     Rscript 03_geneSymbolMap.R      (needs internet; after 01)
# Output:  Data/ref/annot/ensembl_hgnc_map.csv   (ensembl_gene_id, hgnc_symbol)
# Used by: gene symbols in Figures 2-5 and Supplementary Tables 1-13
# =============================================================================

.a <- commandArgs(FALSE); .f <- sub('^--file=', '', .a[grep('^--file=', .a)])
source(file.path(if (length(.f)) dirname(normalizePath(.f)) else '.', '00_config.R'))
library(httr)
library(jsonlite)

ids <- unlist(lapply(harmonicFiles(), function(f) {
  d <- fread(f)
  d$gene[p.adjust(d$p.value, 'bonferroni') < 0.05]
}))
ids <- sort(unique(sub('\\..*', '', ids)))
message('IDs to resolve: ', length(ids))

lookup <- function(batch, attempts = 4) {
  for (k in seq_len(attempts)) {
    r <- tryCatch(POST('https://rest.ensembl.org/lookup/id',
                       body = toJSON(list(ids = batch)), encode = 'raw',
                       content_type_json(), accept_json(), timeout(120)),
                  error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200)
      return(fromJSON(content(r, as = 'text', encoding = 'UTF-8'), simplifyVector = FALSE))
    Sys.sleep(5 * k)
  }
  stop('Ensembl REST lookup failed')
}

batches <- split(ids, ceiling(seq_along(ids) / 500))
map <- bind_rows(lapply(batches, function(b) {
  res <- lookup(b)
  Sys.sleep(0.5)
  bind_rows(lapply(b, function(id) {
    name <- res[[id]]$display_name
    if (is.null(name) || name == '') NULL else data.frame(ensembl_gene_id = id, hgnc_symbol = name)
  }))
})) %>% arrange(ensembl_gene_id)

dir.create(dirname(SYMBOL_MAP), recursive = TRUE, showWarnings = FALSE)
fwrite(map, SYMBOL_MAP)
message('Resolved ', nrow(map), '/', length(ids), ' to a symbol; wrote ', SYMBOL_MAP)
