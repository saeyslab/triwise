#' @importFrom Rcpp evalCpp
#' @useDynLib triwise, .registration=TRUE
x <- 1



#' Barycentric transformation matrix
#'
#' Get the matrix for the barycentric transformation
#'
#' @param anglebase Number of radians the first barycentric direction should be rotated anticlockwise
#' @examples
#' plot(t(getTransformationMatrix(0)), asp=1)
#' plot(t(getTransformationMatrix(pi/2)), asp=1)
#' plot(t(getTransformationMatrix(pi)), asp=1)
#' @return 2 by 3 matrix
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

hexagonPolar <- function(angle, radius=1, baseangle=0) {
  delta <- 2*pi/6
  cos(delta/2)/cos((angle %% delta)-delta/2) * radius
}

# will also rotate points given the baseangle
clipHexagon <- function(barypoints, rmax, baseangle=0) {
  barypoints["rclip"] = mapply(function(angle, r) {min(hexagonPolar(angle, rmax, baseangle), r)}, barypoints$angle, barypoints$r)
  barypoints["xclip"] = cos(barypoints["angle"]+baseangle) * barypoints["rclip"]
  barypoints["yclip"] = sin(barypoints["angle"]+baseangle) * barypoints["rclip"]

  barypoints
}

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)),
             error = function(e) FALSE)
  })
}

betweenCircular <- function(angle, angle1, angle2) {
  (angle-angle1)%%(pi*2) < (angle2-angle1)
}

diffCircular <- function(a, b) {
  diff = (a-b)%%(2*pi)
  if (diff > pi) {
    return(-(2*pi - diff))
  } else {
    return(diff)
  }
}

circularMean <- function(angles, rs=1){
  atan2(mean(sin(angles) * rs), mean(cos(angles) * rs))
}

circularZ <- function(angles) {
  sqrt(sum(cos(angles))**2 + sum(sin(angles))**2)/length(angles)
}

circularZ <- function(angles, rs=1) {
  sqrt(sum(cos(angles) * rs)**2 + sum(sin(angles) * rs)**2)/length(angles)
}

jaccard= function(a, b) {length(intersect(a, b))/length(union(a, b))}

named.list <- function(...) {
  l <- list(...)
  names(l) <- sapply(substitute(list(...)), deparse)[-1]
  l
}

seqClosed <- function(a=0, b, length) {
  utils::head(seq(a, b, length=length+1), -1)
}

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

itself = function(x) x
