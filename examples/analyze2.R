.onLoad <- function(libname, pkgname) {
}

#' Conversion to barycentric coordinates
#'
#' Converts the expression matrix containing 3 biological conditions to barycentric coordinates.
#'
#' @param E Expression matrix
#' @param transfomatrix Transformation matrix
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @export
transformBarycentric = function(E, transfomatrix=NULL) {
  if (is.null(transfomatrix)) {
    anglebase = 0
    transfomatrix = matrix(
      c(
        cos(anglebase), sin(anglebase),
        cos(pi*2/3+anglebase), sin(pi*2/3+anglebase),
        cos(-pi*2/3+anglebase), sin(-pi*2/3+anglebase)
      ),
      nrow=2, ncol=3
    )
  }

  barycoords = as.data.frame(t(transfomatrix %*% t(E)))
  colnames(barycoords) = c("x", "y")
  barycoords$angle = atan2(barycoords$y, barycoords$x)
  barycoords$r = sqrt(barycoords$y**2 + barycoords$x**2)
  barycoords
}

#' @export
transformReverseBarycentric = function(barycoords, transfomatrix=NULL) {
  if (is.null(transfomatrix)) {
    anglebase = 0
    transfomatrix = matrix(
      c(
        cos(anglebase), sin(anglebase),
        cos(pi*2/3+anglebase), sin(pi*2/3+anglebase),
        cos(-pi*2/3+anglebase), sin(-pi*2/3+anglebase)
      ),
      nrow=2, ncol=3
    )
  }

  barycoords = do.call(cbind, list(x=cos(scores$angle), y=sin(scores$angle)))
  enrichmat = barycoords3 %*% t(ginv(transfomatrix))
  rownames(enrichmat) = scores$name

  enrichmat
}

diffexp <- function(expr, conditions=phenoData(expr)$condition, expr_type="counts")
  if (expr_type == "counts") {
    expr = limma::voom(expr, conditions)
  } else if (expr_type == "log2") {

  }

#' Simple test for enrichment in a given set of genes of interest (Goi)
#'
#'
#'
#' @param Goi character vector containing the genes tested for enrichment
#' @param gsets a GeneSetCollection or a
#' @param background character vector
#' @return Dataframe containing for every gene set which passed the filtering: \itemize{
#'  \item p-value of enrichment
#'  \item q-value of enrichment (p-value corrected for multiple testing)
#'  \item estimated odds score of enrichment, how much more likely is a given gene set to be part of
#' }
#' @export
testEnrichment = function(Goi, gsets, background, minknown=2, minfound=2, maxknown=500) {
  scores = do.call(rbind.data.frame, lapply(gsets, function(gset){
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    # much faster than using table
    tp = length(intersect(Goi, gset))
    fn = length(Goi) - tp
    fp = length(gset) - tp
    tn = length(background) - tp - fp - fn

    contingency_table = matrix(c(tp, fp, fn, tn), 2, 2)

    if (contingency_table[1,1] < minfound){
      return(NULL)
    }

    fisher_result = fisher.test(contingency_table, alternative="greater")
    list(pval=fisher_result$p.value, odds=fisher_result$estimate, found=contingency_table[1,1])
  }))
  scores=dplyr::add_rownames(scores, var = "gsetid")
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
  }
  scores
}

#' Test gene sets for unidirectional enrichment
#'
#'
#'
#' @param angles either a dataframe obtained from transformBarycentric or a list of angles with gene ids as names
#' @param gsets
#' @param Gdiffexp differentially expressed genes
#' @return Dataframe containing for every gene set which passed the filtering: \itemize{
#'  \item p-value of unidirectionality
#'  \item q-value of unidirectionality (p-value corrected for multiple testing)
#'  \item average angle
#' }
#' @export
testUnidirectionality = function(angles, gsets, Gdiffexp=NULL, minknown=2, minfound=2, maxknown=500, weight=T, angleweights = NULL) {
  if (is.data.frame(angles)) {
    angles = setNames(angles$angle, rownames(angles))
  }

  if (!is.null(Gdiffexp)) {
    angles = angles[Gdiffexp]
  }
  background = names(angles)

  scores = do.call(rbind.data.frame, lapply(gsets, function(gset){
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    gset_filtered = intersect(gset, background)

    if (length(gset_filtered) < minfound){
      return(NULL)
    }


    angles_gset = circular::circular(angles[gset_filtered], type="angles", units="radians")

    angle = as.numeric(circularMean(angles_gset)) %% (2*pi)

    if (weight) {
      angleweight = angleweights[[floor(angle / (2*pi) * length(angleweights))+1]]
    } else {
      angleweight = 1
    }
    #rayleigh_result = circular::rayleigh.test(angles_gset)

    list(pval=testRayleigh(angles_gset, angleweight), angle=angle, n=length(gset_filtered))
  }))

  scores=dplyr::add_rownames(scores, var = "gsetid")
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
  }
  scores
}

rescale = function(weights, total=length(weights)) {weights*total/sum(weights)}
#' @export
testRayleigh = function(angles, weight=1) {
  n = length(angles)

  z = (n * (sqrt(sum(sin(angles))**2 + sum(cos(angles) )**2)/n)**2) * weight
  pval = exp(-z)

  tmp = 1 + (2 * z - z**2)/(4 * n) - (24 * z - 132 * z**2 + 76 * z**3 - 9 * z**4)/(288 * n**2)

  max(0, min(pval * tmp, 1))
}

testLocality = function(Goi, Gdiffexp, angles, deltangle=pi/24, bandwidth=pi/3) {
  if (is.data.frame(angles)) {
    angles = setNames(angles$angle, rownames(angles))
  }

  Gdiffexp = names(angles) %in% Gdiffexp
  Goi = names(angles) %in% Goi

  localpvals = lapply(seq(0, pi*2-0.001, deltangle), function(angle) {
    angle1 = angle - bandwidth/2
    angle2 = angle + bandwidth/2

    p = betweenCircular(angles, angle1, angle2) & Gdiffexp
    n = !p

    contingency = table(p, Goi)

    fisher_result = fisher.test(contingency, alternative="greater")
    fisher_result$p.value
  })

  p.adjust(localpvals, method="fdr")
}
