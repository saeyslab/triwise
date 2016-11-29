#' Conversion to barycentric coordinates
#'
#' Converts the expression matrix containing three biological conditions to barycentric coordinates, reducing its dimensionality by one while retaining information of differential expression.
#'
#' @param E Expression matrix
#' @param transfomatrix \code{NULL} (default) or a numeric matrix containing the transformation matrix
#' @return Dataframe containing for every original point its new coordinates in `d-1` dimensions
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 1
#' barycoords = transformBarycentric(Eoi)
#' hist(barycoords$angle)
#' @export
transformBarycentric = function(E, transfomatrix=NULL) {
  if (is.null(transfomatrix)) {
    transfomatrix = getTransformationMatrix()
  }

  barycoords = as.data.frame(t(transfomatrix %*% t(E)))
  colnames(barycoords) = c("x", "y")
  barycoords$angle = atan2(barycoords$y, barycoords$x)
  barycoords$r = sqrt(barycoords$y**2 + barycoords$x**2)

  attr(barycoords, "conditions") = colnames(E)

  barycoords
}

#' Convert from barycentric coordinates back to original expression values
#'
#' Converts a dataframe contain barycentric coordinates in the x and y columns back to original coordinates (apart from a constant for each gene).
#'
#' @param barycoords Dataframe of barycentric coordinates with x and y in separate columns
#' @param transfomatrix \code{NULL} (default) or a numeric matrix containing the transformation matrix
#' @return Matrix containing for every gene original expression values centered around 0
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 1
#' barycoords = transformBarycentric(Eoi)
#' reversebarycoords = transformReverseBarycentric(barycoords)
#' all.equal(scale(t(Eoi)), scale(t(reversebarycoords)), check.attributes=FALSE)
#' @export
transformReverseBarycentric = function(barycoords, transfomatrix=NULL) {
  if (is.null(transfomatrix)) {
    transfomatrix = getTransformationMatrix()
  }

  barycoords = do.call(cbind, list(x=cos(barycoords$angle), y=sin(barycoords$angle)))
  barycoords %*% t(MASS::ginv(transfomatrix))
}

#' Rayleigh z-test implementation
#' @param angles Numeric vector containing angles in radians
#' @return P-value of unidirectionality under the uniformity null hypothesis
#' @examples
#' testRayleigh(as.numeric(circular::rvonmises(10, circular::circular(1), 2)))
#' testRayleigh(seq(0, 2*pi, 0.1))
#' @export
testRayleigh = function(angles) {
  n = length(angles)

  z = (n * (sqrt(sum(sin(angles))**2 + sum(cos(angles) )**2)/n)**2)
  pval = exp(-z)

  tmp = 1 + (2 * z - z**2)/(4 * n) - (24 * z - 132 * z**2 + 76 * z**3 - 9 * z**4)/(288 * n**2)

  max(0, min(pval * tmp, 1))
}

#' Simple fisher's exact test for enrichment in a given set of genes of interest (Goi)
#'
#' @param Goi character vector containing the genes tested for enrichment
#' @param background character vector containing background gene identifiers
#' @inheritParams testUnidirectionality
#' @return Dataframe containing for every gene set which passed the filtering: \itemize{
#'  \item p-value of enrichment
#'  \item q-value of enrichment (p-value corrected for multiple testing)
#'  \item estimated odds score of enrichment, how much more likely is a given gene set to be part of
#' }
#' @examples
#' Goi = 1:50
#' gsets = list(a= 20:70, b = 45:80)
#' background = 1:100
#' testEnrichment(Goi, gsets, background)
#' @export
testEnrichment = function(Goi, gsets, background, minknown=2, mindiffexp=2, maxknown=500) {
  scores = do.call(rbind.data.frame, lapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    # much faster than using table
    tp = length(intersect(Goi, gset))
    fn = length(Goi) - tp
    fp = length(gset) - tp
    tn = length(background) - tp - fp - fn

    contingency_table = matrix(c(tp, fp, fn, tn), 2, 2)

    if (contingency_table[1,1] < mindiffexp){
      return(NULL)
    }

    fisher_result = stats::fisher.test(contingency_table, alternative="greater")
    list(pval=fisher_result$p.value, odds=fisher_result$estimate, found=contingency_table[1,1], gsetid=gsetid)
  }))
  if (nrow(scores) > 0) {
    scores$qval = stats::p.adjust(scores$pval, method="fdr")
  }
  scores
}

#' Generate background models
#' @description Generates a background model by randomly resampling genes at different `n` (number of genes) and angles and calculating z distributions
#' @inheritParams testUnidirectionality
#' @param barycoords Dataframe containing for every gene its barycentric coordinates, as returned by \code{\link[triwise]{transformBarycentric}}. Will use the \code{z} column as test statistic, or if this column is not given the `r` column
#' @param noi Integer vector denoting the number of genes at which to sample, the larger the more accurate the p-values
#' @param anglesoi Numeric vector denoting the angles (in radians) at which to pre-calculate null distribution, the larger the more accurate the p-values
#' @param nsamples Number of samples, higher for more accurate and stable p-values
#' @param bw Bandwidth of the von-mises distribution for weighing the samples. A higher bandwidth leads to a more accurate p-value estimate as long as `nsamples` is high enough
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 1
#' barycoords = transformBarycentric(Eoi)
#'
#' hist(barycoords$angle)
#'
#' bm = generateBackgroundModel(barycoords)
#'
#' # the distribution of mean angle of the samples is not uniform due to the
#' # non-uniform distribution of the angles of individual genes
#' hist(bm$backmodels[[1]]$angles)
#'
#' # the whole distribution (and therefore also the p-value) also depends on the mean angle
#' plotdata = data.frame(angle = cut(bm$backmodels[[1]]$angles, 10), z = bm$backmodels[[1]]$z)
#' ggplot2::ggplot(plotdata) + ggplot2::geom_violin(ggplot2::aes(angle, z))
#' @return A list containing: \itemize{
#' \item noi: number of genes, in the same order as the elements in \code{backmodels}
#' \item anglesoi: angles at which weights were calculated using the von-mises distribution
#' \item nsamples: number of samples for each \code{n}
#' \item bw: bandwidth
#' \item backmodels: for each \code{n} a second list containing: \itemize{
#'   \item angles: mean angle of a sample
#'   \item z: strength of unidirectional upregulation
#'   \item weights: weight from von-mises distribution for every sample and angle in \code{anglesoi}, these weights will be used to calculate the p-value
#'   }
#' }
#' @export
generateBackgroundModel <- function(barycoords, noi = seq(5, 100, 5), anglesoi = seqClosed(0, 2*pi, 24), nsamples=100000, bw=20, mc.cores=getOption("mc.cores", default = 1)) {
  if (length(noi) == 1) {
    noi = seq(5, 100, length=noi)
  }

  if (length(anglesoi) == 1) {
    anglesoi = seqClosed(0, 2*pi, anglesoi)
  }

  if (!("z" %in% colnames(barycoords))) {
    barycoords$z = barycoords$r
  }

  backmodels = parallel::mclapply(noi, function(n) {
    backmodel = backgroundModel2(barycoords$angle, barycoords$z, nsamples, n, anglesoi, bw)
  }, mc.cores = mc.cores)

  named.list(noi, anglesoi, nsamples, bw, backmodels)
}


#' Test gene sets for unidirectional enrichment
#'
#' @param barycoords Dataframe containing for every gene its barycentric coordinates, as returned by \code{\link[triwise]{transformBarycentric}}
#' @param gsets List of character vectors, each containing a set of genes (gene identifiers)
#' @param Gdiffexp Differentially expressed genes
#' @param statistic A string denoting the measure used for the strength of upregulation of a particular gene. \itemize{
#'   \item diffexp: Whether a gene is differentially expressed (1 versus 0)
#'   \item rank: The rank of the maximal log fold-change
#'   \item r: The maximal log fold-change
#'   \item z: Custom using the z-column within barycoords
#'   \item angle: Uses a rayleigh z-test ignoring non-differentially expressed genes within the gene set
#' }
#' @param bm Previously calculated background model using the \code{\link[triwise]{generateBackgroundModel}}
#' @param minknown Minimal number of genes within a gene set for it to be considered for enrichment
#' @param mindiffexp Minimal number of genes differentially expressed within a gene set for it to be considered for enrichment
#' @param maxknown Maximal number of genes within a gene set for it to be considered for enrichment
#' @param mc.cores Number of processor cores to use. Due to limitations of the parallel package, this does not work on Windows
#' @param ... Parameters for \code{\link[triwise]{generateBackgroundModel}}
#' @return Dataframe containing for every gene set which passed the filtering: \itemize{
#'  \item p-value of unidirectionality
#'  \item q-value of unidirectionality (p-value corrected for multiple testing)
#'  \item average angle
#' }
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 4 # the first 100 genes are more upregulated in the first condition
#' barycoords = transformBarycentric(Eoi)
#'
#' gsets = list(a=1:50, b=80:150, c=200:500)
#' testUnidirectionality(barycoords, gsets, Gdiffexp=(1:1000)[barycoords$r > 1])
#' @export
testUnidirectionality = function(barycoords, gsets, Gdiffexp=NULL, statistic="diffexp", bm=NULL, minknown=5, mindiffexp=0, maxknown=1500, mc.cores=getOption("mc.cores", default=1), ...) {
  if(!is.data.frame(barycoords)) stop("barycoords should be a data.frame")
  if(!all(c("x", "y", "angle", "r") %in% colnames(barycoords))) stop("barycoords should contain x, y, angle and r columns")

  if(!is.null(Gdiffexp)) {
    barycoords$diffexp = rownames(barycoords) %in% Gdiffexp
  }

  if(statistic == "diffexp" || statistic == "angle") {
    if(is.null(Gdiffexp)) stop("Gdiffexp should be given if statistic == \"diffexp\" or statistic == \"angle\"")
    barycoords$z = rownames(barycoords) %in% Gdiffexp
  } else if(statistic == "rank") {
    barycoords$z = rank(barycoords$r)
  } else if(statistic == "r"){
    barycoords$z = barycoords$r
  } else if(statistic == "z"){
    # otherwise use the z column from the user
  } else {
    stop("no valid statistic provided (diffexp, rank, r, z or angle)")
  }

  if(is.null(bm) && statistic != "angle") {
    bm = generateBackgroundModel(barycoords, mc.cores=mc.cores, ...)
  }

  background = rownames(barycoords)
  scores = dplyr::bind_rows(parallel::mclapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    gset_filtered = intersect(gset, background)

    subbarycoords = barycoords[gset_filtered, ]

    if (length(gset_filtered) == 0) return(NULL)

    angles_gset = subbarycoords$angle %% (2*pi)
    rs_gset = as.numeric(subbarycoords$z)
    angle = as.numeric(circularMean(angles_gset, rs_gset)) %% (2*pi)

    if (sum(subbarycoords$diffexp) < mindiffexp) {
      return(NULL)
    }

    if(statistic != "angle") {
      pval = empiricalPvalue(angles_gset,rs_gset, bm)
    } else {
      angles_gset_filtered = angles_gset[subbarycoords$z > 0]
      pval = testRayleigh(angles_gset_filtered)
    }

    data.frame(pval=pval, angle=angle, n=length(rs_gset), gsetid=gsetid,z=circularZ(angles_gset, rs_gset), stringsAsFactors=FALSE)
  }, mc.cores = mc.cores))
  if (nrow(scores) > 0) {
    scores$qval = stats::p.adjust(scores$pval, method="fdr")
  }
  scores
}


empiricalPvalue = function(angles, rs, bm) {
  nid = which.min(abs(bm$noi-length(angles)))
  backmodel = bm$backmodels[[nid]]

  if (max(rs) == 0) {
    return(1)
  }

  z = circularZ(angles, rs)

  angle = circularMean(angles, rs) %% (2*pi)
  angleid = which.min(abs(sapply(bm$anglesoi, diffCircular, angle)))

  higher = backmodel$z > z
  sum(backmodel$weights[angleid,higher]) * backmodel$basepval# * backmodel$anglesp[[angleid]]
}

#' Test local upregulation
#' @description Tests for local upregulation locally in certain directions
#' @param Goi Gene set for which to test enrichment
#' @param Gdiffexp Differentially expressed genes
#' @param angles Numeric vector with angles of all genes or a dataframe as returned by `transformBarycentric`
#' @param deltangle Stepsize of angles
#' @param bandwidth Bandwidth of angles
#' @return Corrected p-values for enrichment around every angle
#' @examples
#' Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
#' Eoi[1:100,1] = Eoi[1:100,1] + 4 # the first 100 genes are more upregulated in the first condition
#' barycoords = transformBarycentric(Eoi)
#'
#' Goi = 1:50
#' qvals = testLocality(Goi, Gdiffexp=(1:1000)[barycoords$r > 1], barycoords)
#' plot(attr(qvals, "anglesoi"), qvals)
#' @export
testLocality = function(Goi, Gdiffexp, angles, deltangle=pi/24, bandwidth=pi/3) {
  if (is.data.frame(angles)) {
    angles = stats::setNames(angles$angle, rownames(angles))
  }

  Gdiffexp = names(angles) %in% Gdiffexp
  Goi = names(angles) %in% Goi

  anglesoi = seq(0, pi*2-0.001, deltangle)

  localpvals = lapply(anglesoi, function(angle) {
    angle1 = angle - bandwidth/2
    angle2 = angle + bandwidth/2

    p = betweenCircular(angles, angle1, angle2) & Gdiffexp
    n = !p

    contingency = table(p, Goi)

    fisher_result = stats::fisher.test(contingency, alternative="greater")
    fisher_result$p.value
  })

  qvals = stats::p.adjust(localpvals, method="fdr")
  attr(qvals, "anglesoi") = anglesoi

  qvals
}
