#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
interactiveDotplot <- function(Eoi, Gdiffexp, Glabels, Goi, Gpin = NULL, width = NULL, height = NULL) {
  barycoords = transformBarycentric(Eoi)

  # calculate local pvalues
  localpvals = testLocality(Goi, Gdiffexp, barycoords)

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
      index=colnames(Eoi)
    ),
    Gdiffexp = names(Gmap) %in% Gdiffexp,
    Glabels = Glabels,
    Goi = Gmap[Goi],
    Gpin = Gmap[Gpin],
    logpvals=log10(localpvals)
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
