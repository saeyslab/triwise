hexagonPolar <- function(angle, radius=1) {
  delta <- 2*pi/6
  cos(delta/2)/cos((angle %% delta)-delta/2) * radius
}

clipHexagon <- function(barycoords, rmax) {
  barycoords["rclip"] = apply(barycoords, 1, function(row) min(hexagonPolar(row["angle"], rmax), row["r"]))
  barycoords["xclip"] = cos(barycoords["angle"]) * barycoords["rclip"]
  barycoords["yclip"] = sin(barycoords["angle"]) * barycoords["rclip"]

  barycoords
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

circularMean <- function(angles, weights){
  atan2(sum(sin(angles) * weights)/sum(weights), sum(cos(angles) * weights)/sum(weights))
}


circularMean <- function(angles, weights){
  atan2(mean(sin(angles)), mean(cos(angles)))
}

circularZ <- function(angles) {
  sqrt(sum(cos(angles))**2 + sum(sin(angles))**2)
}

jaccard = function(a, b) {length(intersect(a, b))/length(union(a, b))}
