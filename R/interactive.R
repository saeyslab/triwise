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
#' @examples
#' data(vandelaar)
#' Eoi <- limma::avearrays(vandelaar, Biobase::phenoData(vandelaar)$celltype)
#' Eoi = Eoi[,c("YS_MF", "FL_mono", "BM_mono")]
#' barycoords = transformBarycentric(Eoi)
#' interactiveDotplot(barycoords)
#'
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 1
#' barycoords = transformBarycentric(Eoi)
#' Gdiffexp =(1:1000)[barycoords$r > 1]
#' interactiveDotplot(Eoi)
#' interactiveDotplot(Eoi, Gdiffexp)
#' interactiveDotplot(Eoi, as.character(Gdiffexp), as.character(1:10), as.character(1:1000))
#' interactiveDotplot(Eoi, as.character(Gdiffexp), as.character(1:10), as.character(1:1000), c(50, 200))
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
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 4 # the first 100 genes are more upregulated in the first condition
#' barycoords = transformBarycentric(Eoi)
#'
#' gsets = list(a=1:50, b=80:150, c=200:500)
#' scores = testUnidirectionality(barycoords, gsets, Gdiffexp=(1:1000)[barycoords$r > 1])
#'
#' interactivePvalplot(scores, as.list(setNames(names(gsets), names(gsets))), 1:3)
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
