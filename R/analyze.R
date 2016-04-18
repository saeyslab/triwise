#' Conversion to barycentric coordinates
#'
#' Converts the expression matrix containing two or more biological conditions to barycentric coordinates, reducing its dimensionality by one while retaining information of differential expression.
#'
#' @param E Expression matrix
#' @return Dataframe containing for every original point its new coordinates in `d-1` dimensions
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
#' Converts a dataframe contain barycentric coordinates in the x and y columns back to original coordinates (apart from a constant shift for each gene).
#'
#' @param barycoords Dataframe of barycentric coordinates with x and y in separate columns
#' @param transfomatrix \code{NULL} (default) or a numeric matrix containing the transformation
#' @return Matrix containing for every gene
#' @export
transformReverseBarycentric = function(scores, transfomatrix=NULL) {
  if (is.null(transfomatrix)) {
    transfomatrix = getTransformationMatrix()
  }

  barycoords = do.call(cbind, list(x=cos(barycoords$angle), y=sin(barycoords$angle)))
  barycoords %*% t(MASS::ginv(transfomatrix))
}



#' Test gene sets for unidirectional enrichment
#'
#' @param angles either a dataframe obtained from transformBarycentric or a list of angles with gene identifiers as names
#' @param gsets list of character vectors, each containing a set of genes (gene identifiers)
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

  scores = do.call(rbind.data.frame, lapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    gset_filtered = intersect(gset, background)

    if (length(gset_filtered) < minfound){
      return(NULL)
    }

    angles_gset = circular::circular(angles[gset_filtered], type="angles", units="radians")

    angle = as.numeric(circularMean(angles_gset)) %% (2*pi)

    rayleigh_result = circular::rayleigh.test(angles_gset)

    list(pval=testRayleigh(angles_gset), angle=angle, n=length(gset_filtered), gsetid=gsetid)
  }))

  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
  }
  scores
}

#' Rayleigh z-test
#'
#' @param angles Numeric vector containing angles in radians
#' @return P-value of unidirectionality under the uniformity null hypothesis
#' @export
testRayleigh = function(angles) {
  n = length(angles)

  z = (n * (sqrt(sum(sin(angles))**2 + sum(cos(angles) )**2)/n)**2)
  pval = exp(-z)

  tmp = 1 + (2 * z - z**2)/(4 * n) - (24 * z - 132 * z**2 + 76 * z**3 - 9 * z**4)/(288 * n**2)

  max(0, min(pval * tmp, 1))
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

    if (contingency_table[1,1] < minfound){
      return(NULL)
    }

    fisher_result = fisher.test(contingency_table, alternative="greater")
    list(pval=fisher_result$p.value, odds=fisher_result$estimate, found=contingency_table[1,1], gsetid=gsetid)
  }))
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
  }
  scores
}







#####

#' Generate background models
#' @description Generates a background model by randomly resampling genes at different n and calculating z distributions at different n and angles
#' @export
generateBackgroundModel <- function(barycoords, noi = seq(5, 100, 5), anglesoi = seqClosed(0, 2*pi, 24), nsamples=100000, bw=20, mc.cores=options("mc.cores")) {
  barycoords$z = barycoords$r#rank(barycoords$r)

  if (length(noi) == 1) {
    noi = seq(5, 100, length=noi)
  }

  if (length(anglesoi) == 1) {
    anglesoi = seqClosed(0, 2*pi, anglesoi)
  }

  backmodels = mclapply(noi, function(n) {
    backmodel = triwise::backgroundModel2(barycoords$angle, barycoords$z, nsamples, n, anglesoi, bw)
  }, mc.cores = mc.cores)

  named.list(noi, anglesoi, nsamples, bw, backmodels)
}



#' @export
testUnidirectionality = function(barycoords, gsets, bm=NULL, minknown=5, minfound=5, maxknown=500, mc.cores=options("mc.cores")) {
  barycoords$z = barycoords$r#rank(barycoords$r)

  if(is.null(bm)) {
    bm = generateBackgroundModel(barycoords)
  }

  background = rownames(barycoords)
  scores = dplyr::bind_rows(mclapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]
    if (length(gset) < minknown | length(gset) > maxknown) {
      return(NULL)
    }

    gset_filtered = intersect(gset, background)

    subbarycoords = barycoords[gset_filtered, ]

    angles_gset = subbarycoords$angle %% (2*pi)
    rs_gset = as.numeric(subbarycoords$z)
    angle = as.numeric(circularMean(angles_gset, rs_gset)) %% (2*pi)

    if (length(rs_gset) < minfound) {
      return(NULL)
    }

    data.frame(pval=empiricalPvalue(angles_gset,rs_gset, bm), angle=angle, n=sum(rs_gset), gsetid=gsetid,z=circularZ(angles_gset, rs_gset), stringsAsFactors=F)
  }, mc.cores = mc.cores))
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="fdr")
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
  angleid = which.min(abs(sapply(bm$anglesoi, triwise::diffCircular, angle)))

  higher = backmodel$z > z
  sum(backmodel$weights[angleid,higher]) * backmodel$basepval# * backmodel$anglesp[[angleid]]
}

