# simple: overlap with any other already explained differentially expressed gene
# disadvantage: cell cycle is also informative even if a lot of its subterms already went by
# major disadvantage: if cell cycle has the highest enrichment, all its informative child terms will be gone
Gexplained = c()
scores$overlap = sapply(scores$gsetid, function(gsetid){
  gsetdiffexp =intersect(gsets[[gsetid]], Gdiffexp)
  overlap = length(intersect(gsetdiffexp, Gexplained))/length(gsetdiffexp)

  Gexplained <<- union(Gexplained, gsetdiffexp)

  overlap
})

# looking at the jaccard of every previously explained gset
# two ways: one with the matrix (faster with more gsets?) the other calculates the intersect without any matrix
# I found that in practice the time difference is only minimal
scores = scores[scores$qval < 0.1, ]
gsetmatrix = t(sapply(scores$gsetid, function(gsetid) {Gdiffexp %in% gsets[[gsetid]]}))
colnames(gsetmatrix) = Gdiffexp
sapply(c(1:nrow(gsetmatrix)), function(rowid) {
  if (rowid > 1) {
    row = gsetmatrix[rowid, ]
    max(sapply(c(1:rowid-1), function(rowid2) {sum(gsetmatrix[rowid2, ] & row)}))
  } else {
    return(0)
  }
})

#
filteredgsets = sapply(scores$gsetid, function(gsetid) {intersect(gsets[[gsetid]], Gdiffexp)})
last = c()
gsetidorder = sapply(scores$gsetid, function(gsetid) {
  last <<- c(gsetid, last)
  last
})
scores$overlap = sapply(gsetidorder, function(gsetids){
  if (length(gsetids) > 1) {
    gsetoi = filteredgsets[[gsetids[[1]]]]

    jaccards = sapply(filteredgsets[gsetids[c(2:length(gsetids))]], jaccard, b=gsetoi)

    return(max(jaccards))
  } else {
    return(0)
  }
})
