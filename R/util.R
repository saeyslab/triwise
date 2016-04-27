#' @importFrom Rcpp evalCpp
#' @useDynLib triwise


#' @title Barycentric transformation matrix
#'
#' @description Get the matrix for the barycentric transformation
#'
#' @param anglebase Number of radians the first barycentric direction should be rotated anticlockwise
#' @export
getTransformationMatrix <- function(anglebase=0) {
  matrix(
    c(
      cos(anglebase), sin(anglebase),
      cos(pi*2/3+anglebase), sin(pi*2/3+anglebase),
      cos(-pi*2/3+anglebase), sin(-pi*2/3+anglebase)
    ),
    nrow=2, ncol=3
  )
}

addPolar <- function(barypoints) {
  barypoints$angle = atan2(barypoints$y, barypoints$x)
  barypoints$r = sqrt(barypoints$x^2 + barypoints$y^2)

  barypoints
}

hexagonPolar <- function(angle, radius=1) {
  delta <- 2*pi/6
  cos(delta/2)/cos((angle %% delta)-delta/2) * radius
}

clipHexagon <- function(barypoints, rmax) {
  barypoints["rclip"] = mapply(function(angle, r) {min(hexagonPolar(angle, rmax), r)}, barypoints$angle, barypoints$r)
  barypoints["xclip"] = cos(barypoints["angle"]) * barypoints["rclip"]
  barypoints["yclip"] = sin(barypoints["angle"]) * barypoints["rclip"]

  barypoints
}

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

betweenCircular <- function(angle, angle1, angle2) {
  (angle-angle1)%%(pi*2) < (angle2-angle1)
}

#' @export
diffCircular <- function(a, b) {
  diff = (a-b)%%(2*pi)
  if (diff > pi) {
    return(-(2*pi - diff))
  } else {
    return(diff)
  }
}

#' @export
circularMean <- function(angles, rs=1){
  atan2(mean(sin(angles) * rs), mean(cos(angles) * rs))
}

#' @export
circularZ <- function(angles) {
  sqrt(sum(cos(angles))**2 + sum(sin(angles))**2)/length(angles)
}

#' @export
circularZ <- function(angles, rs) {
  sqrt(sum(cos(angles) * rs)**2 + sum(sin(angles) * rs)**2)/length(angles)
}



#' @export
jaccard= function(a, b) {length(intersect(a, b))/length(union(a, b))}


#' make a named list using the variable names (avoids repeating the same name twice)
#' @export
named.list <- function(...) {
  l <- list(...)
  names(l) <- sapply(substitute(list(...)), deparse)[-1]
  l
}

#'
#' @export
seqClosed <- function(a=0, b, length) {
  head(seq(a, b, length=length+1), -1)
}
