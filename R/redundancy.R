#' Estimates the redundancy of enriched gene sets
#'
#' Wrapper around the `mgsa` method
#'
#' @param scores Dataframe from `r packagedocs::rd_link(testUnidirectionality())` containing the filtered enrichment scores for every gene set
#' @param gsets Named list containing character vectors for each gene set
#' @param Gdiffexp Character vector of differentially expressed genes
#' @return Numeric vector containing for every gene set (in the same order the scores dataframe) its relevance compared with other gene sets, higher is better
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 4 # the first 100 genes are more upregulated in the first condition
#' barycoords = transformBarycentric(Eoi)
#' Gdiffexp = (1:1000)[barycoords$r > 1]
#'
#' # a and b are redundant, but a is stronger enriched
#' gsets = list(a=1:50, b=c(1:50, 100:110), c=200:500)
#' scores = testUnidirectionality(barycoords, gsets, Gdiffexp=(1:1000)[barycoords$r > 1])
#' scores$redundancy = estimateRedundancy(scores, gsets, Gdiffexp)
#' scores[scores$gsetid == "a", "redundancy"] > scores[scores$gsetid == "b", "redundancy"]
#' @export
estimateRedundancy <- function(scores, gsets, Gdiffexp) {
  gsets_filtered = gsets[scores$gsetid]
  G = Reduce(union, gsets_filtered)

  redundancy = mgsa::mgsa(Gdiffexp, gsets_filtered, G)
  redundancy = stats::setNames(redundancy@setsResults$estimate, rownames(redundancy@setsResults))
  redundancy[scores$gsetid]
}
