#' Interactive triwise dotplot
#'
#' Draw an interactive triwise dotplot in which genes can be interactively
#'
#' @import htmlwidgets
#'
#' @param Eoi Expression matrix with the three conditions in the columns
#' @param Gdiffexp Differentially expressed genes
#' @param Goi List with genes of interest
#' @param Glabels Labels for every gene if different from the rownames of `Eoi`
#' @param Gpin Pinned genes
#' @export
interactiveDotplot <- function(Eoi, Gdiffexp=c(), Goi=c(), Glabels=rownames(Eoi), Gpin = c(), plotLocalEnrichment=T, width = NULL, height = NULL) {
  Eoi = Eoi[,c(1,3,2)]

  barycoords = transformBarycentric(Eoi)
  barycoords = addPolar(barycoords)

  if(!is.null(Goi) && plotLocalEnrichment) {
    # calculate local pvalues
    localpvals = testLocality(Goi, Gdiffexp, barycoords)
  } else {
    localpvals = rep(0, 100)
  }

  # automatically choose pinned genes if not given
  if (is.null(Gpin)) {
    Gpin = intersect(Goi, Gdiffexp)
    if (length(Gpin) > 20) {
      Gpin = Gpin[order(barycoords[Gpin,"r"])[1:20]]
    }
  }

  Gmap = setNames(seq(0, nrow(Eoi)), rownames(Eoi))

  params <- list(
    Eoi = list(
      data=Eoi,
      columns=Gmap,
      index=attr(barycoords, "conditions")
    ),
    Gdiffexp = names(Gmap) %in% Gdiffexp,
    Glabels = Glabels,
    Goi = Gmap[Goi],
    Gpin = Gmap[Gpin],
    logpvals=log10(localpvals),
    plotLocalEnrichment=plotLocalEnrichment
  )
  attr(params, 'TOJSON_FUNC') <- function(x) {jsonlite::toJSON(x, matrix="columnmajor")}

  print(height)
  # create widget
  htmlwidgets::createWidget(
    'dotplot',
    params,
    width = width,
    height = height,
    package = 'triwise'
  )
}

## pval plot

#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
interactivePvalplot <- function(scores, gsetlabels, labels, width = NULL, height = NULL) {
  scores$logqval_unidir = log10(scores$qval)
  scores$logqval_unidir[is.infinite(log10(scores$qval))] = "-Inf"
  params <- list(
    scores = scores,
    gsetlabels = gsetlabels,
    labels = labels
  )
  attr(params, 'TOJSON_FUNC') <- function(x) {jsonlite::toJSON(x, matrix="columnmajor")}

  # create widget
  htmlwidgets::createWidget(
    'pvalplot',
    params,
    width = width,
    height = height,
    package = 'triwise'
  )
}
