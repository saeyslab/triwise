#' Interactive triwise dotplot
#'
#' Draw an interactive triwise dotplot using a html widget
#'
#' @inheritParams plotDotplot
#' @param Eoi Expression matrix with the three conditions in the columns
#' @param Glabels Labels for every gene if different from the rownames of `Eoi`
#' @param Gpin Pinned genes, if `NULL` will automatically choose the top 20 most differentially expressed genes
#' @param plotLocalEnrichment Whether to plot local enrichment as a ring around the dotplot
#' @param width Width of the plot
#' @param height Height of the plot
#' @export
interactiveDotplot <- function(Eoi, Gdiffexp=rownames(Eoi), Goi=c(), Glabels=rownames(Eoi), Gpin = c(), Coi=colnames(Eoi), colorvalues=NULL, rmax=5, sizevalues=c(T=2, F=0.5), alphavalues=c(T=0.8, F=0.8), plotLocalEnrichment=F, width = NULL, height = NULL) {
  Eoi = Eoi[,c(1,3,2)] # reorder so that ordering corresponds to the ordering of plotDotplot

  barycoords = transformBarycentric(Eoi)
  barycoords = addPolar(barycoords)

  if(!is.null(Goi) && plotLocalEnrichment) {
    # calculate local pvalues
    localpvals = testLocality(Goi, Gdiffexp, barycoords)
  } else {
    localpvals = rep(1, 48)
  }

  # automatically choose pinned genes if not given
  if (is.null(Gpin)) {
    Gpin = intersect(Goi, Gdiffexp)
    if (length(Gpin) > 20) {
      Gpin = Gpin[order(barycoords[Gpin,"r"], decreasing=T)[1:20]]
    }
  }

  Gmap = stats::setNames(seq(0, nrow(Eoi)), rownames(Eoi))

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
    options = list(plotLocalEnrichment=plotLocalEnrichment)
  )
  attr(params, 'TOJSON_FUNC') <- function(x) {jsonlite::toJSON(x, matrix="columnmajor")}

  # create widget
  htmlwidgets::createWidget(
    'dotplot',
    params,
    width = width,
    height = height,
    package = 'triwise'
  )
}

#' Interactive plot of p-values
#'
#' Draw a dotplot of p-values using a htmlwidget
#'
#' @inheritParams plotPvalplot
#' @param gsetlabels Names of the gene sets, used for hovering
#' @param width Width of the plot
#' @param height Height of the plot
#' @export
interactivePvalplot <- function(scores, gsetlabels, Coi, width = NULL, height = NULL) {
  scores$logqval_unidir = log10(scores$qval)
  scores$logqval_unidir[is.infinite(log10(scores$qval))] = "-Inf"
  params <- list(
    scores = scores,
    gsetlabels = gsetlabels,
    labels = Coi
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
