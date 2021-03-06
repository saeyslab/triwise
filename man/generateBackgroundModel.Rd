% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{generateBackgroundModel}
\alias{generateBackgroundModel}
\title{Generate background models}
\usage{
generateBackgroundModel(barycoords, noi = seq(5, 100, 5),
  anglesoi = seqClosed(0, 2 * pi, 24), nsamples = 1e+05, bw = 20,
  mc.cores = getOption("mc.cores", default = 1))
}
\arguments{
\item{barycoords}{Dataframe containing for every gene its barycentric coordinates, as returned by \code{\link[triwise]{transformBarycentric}}. Will use the \code{z} column as test statistic, or if this column is not given the `r` column}

\item{noi}{Integer vector denoting the number of genes at which to sample, the larger the more accurate the p-values}

\item{anglesoi}{Numeric vector denoting the angles (in radians) at which to pre-calculate null distribution, the larger the more accurate the p-values}

\item{nsamples}{Number of samples, higher for more accurate and stable p-values}

\item{bw}{Bandwidth of the von-mises distribution for weighing the samples. A higher bandwidth leads to a more accurate p-value estimate as long as `nsamples` is high enough}

\item{mc.cores}{Number of processor cores to use. Due to limitations of the parallel package, this does not work on Windows}
}
\value{
A list containing: \itemize{
\item noi: number of genes, in the same order as the elements in \code{backmodels}
\item anglesoi: angles at which weights were calculated using the von-mises distribution
\item nsamples: number of samples for each \code{n}
\item bw: bandwidth
\item backmodels: for each \code{n} a second list containing: \itemize{
  \item angles: mean angle of a sample
  \item z: strength of unidirectional upregulation
  \item weights: weight from von-mises distribution for every sample and angle in \code{anglesoi}, these weights will be used to calculate the p-value
  }
}
}
\description{
Generates a background model by randomly resampling genes at different `n` (number of genes) and angles and calculating z distributions
}
\examples{
Eoi = matrix(rnorm(1000*3, sd=0.5), 1000, 3, dimnames=list(1:1000, c(1,2,3)))
Eoi[1:100,1] = Eoi[1:100,1] + 1
barycoords = transformBarycentric(Eoi)

hist(barycoords$angle)

bm = generateBackgroundModel(barycoords)

# the distribution of mean angle of the samples is not uniform due to the
# non-uniform distribution of the angles of individual genes
hist(bm$backmodels[[1]]$angles)

# the whole distribution (and therefore also the p-value) also depends on the mean angle
plotdata = data.frame(angle = cut(bm$backmodels[[1]]$angles, 10), z = bm$backmodels[[1]]$z)
ggplot2::ggplot(plotdata) + ggplot2::geom_violin(ggplot2::aes(angle, z))
}
