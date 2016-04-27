#' @export
estimateRedundancy <- function(scores, gsets, Gdiffexp) {
  scores = scores[order(scores$qval),]
  filteredgsets <- sapply(gsets[scores$gsetid], function(gset) {intersect(gset, Gdiffexp)})
  last <- c()
  prevgsetids = c()
  sapply(scores$gsetid, function(gsetid){
    if (length(prevgsetids) > 0) {
      gsetoi = filteredgsets[[gsetid]]

      jaccards = sapply(filteredgsets[prevgsetids], jaccard, b=gsetoi)
    } else {
      jaccards = c(0)
    }
    prevgsetids <<- c(gsetid, prevgsetids)
    max(jaccards)
  })
}
