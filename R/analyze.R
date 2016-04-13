#' Transform barycentric coordinates in `d-1` dimensions back to regular coordinates in `d` dimensions (centered around 0)
#' @param barypoints Barycentric coordinates
#' @return E Expression matrix centered around zero
#' @export
transformReverseBarycentric = function(barypoints) {
  d = length(barypoints)+1
  transfomatrix = t(regularsimplex(d))

  barypoints %*% t(MASS::ginv(transfomatrix))
}

#' Coordinates of a regular simplex
#'
#' Get the coordinates of a regular simplex with `d+1` corners in `d` dimensions
#'
#' @param E Expression matrix
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @export
# see http://mathoverflow.net/questions/38724/coordinates-of-vertices-of-regular-simplex?answertab=votes#tab-top
# answer by Mariano Suárez-Alvarez and others
regularsimplex = function(d) {
  points = t(sapply(c(1:d), function(i) {
    if (i < d) {
      point = rep(0, d-1)
      point[i]=1
    } else {
      point = rep(-1/(1+sqrt(d)), d-1)
    }
    point
  }))
  points = points - apply(points, 2, sum)/(d) # centering (equal distance to center)
  points = points/(sqrt(sum(points[1,]^2))) # scaling (distance to center = 1)

  points
}

#' Conversion to barycentric coordinates
#'
#' Converts the expression matrix containing two or more biological conditions to barycentric coordinates, reducing its dimensionality by one while retaining information of differential expression.
#'
#' @param E Expression matrix
#' @return Dataframe containing for every original point its new coordinates in `d-1` dimensions
#' @export
transformBarycentric = function(E) {
  transfomatrix = t(regularsimplex(ncol(E)))
  barycoords = t(transfomatrix %*% t(E))
  barycoords
}
