#' Estimates the redundancy of enriched gene sets
#'
#' Wrapper around the `mgsa` method
#'
#' @param scores Dataframe from `r packagedocs::rd_link(testUnidirectionality())` containing the filtered enrichment scores for every gene set
#' @param gsets Named list containing character vectors for each gene set
#' @param Gdiffexp Character vector of differentially expressed genes
#' @return Numeric vector containing for every gene set (in the same order the scores dataframe) its relevance compared with other gene sets, higher is better
#' @export
estimateRedundancy <- function(scores, gsets, Gdiffexp) {
  gsets_filtered = gsets[scores$gsetid]
  G = Reduce(union, gsets_filtered)

  redundancy = mgsa::mgsa(Gdiffexp, gsets_filtered, G)
  redundancy = stats::setNames(redundancy@setsResults$estimate, rownames(redundancy@setsResults))
  redundancy[scores$gsetid]
}
