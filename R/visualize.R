#' @import ggplot2
drawHexagonGrid <- function(rmax=5) {
  hexagonpoints <- plyr::ldply(seq(0, pi*2-0.001, pi/3), function(angle) data.frame(x=cos(angle), y=sin(angle)))
  gridData <- dplyr::bind_rows(lapply(seq_len(rmax), function(r) data.frame(r=r, x=hexagonpoints$x * r, y =hexagonpoints$y * r)))


  ggplot() +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.ticks=element_blank(),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background=element_blank()
    ) +
    geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour="black")
}

#' @import ggplot2
drawDirections <- function(rmax, labels, labelmargin=0.05) {
  directionsData = plyr::ldply(seq(0, pi*2-0.001, pi/3*2), function(angle) {
    x = cos(angle) * rmax
    y = sin(angle) * rmax

    if (y<0){
      va=1
    } else if (y == 0) {
      va=0.5
    } else{
      va=0
    }
    if (x<0){
      ha = 1
    } else if (x < 0.1 & x > -0.1) {
      ha=0.5
    } else{
      ha=0
    }

    data.frame(x=x, y=y, va=va, ha=ha, xlabel=x*(1+labelmargin),ylabel=y*(1+labelmargin))
  })
  directionsData$label = labels

  c(
    geom_text(aes(label=label, x=xlabel, y=ylabel, vjust=va, hjust=ha), data=directionsData),
    geom_segment(aes(xend=0, x=x, yend=0, y=y, group=label), data=directionsData)
  )
}

#' @import ggplot2
drawDotplot <- function(barycoords, rmax=5, color="red", alpha=1, order=NULL) {
  barycoords = clipHexagon(barycoords, rmax)
  barycoords$color = color
  barycoords$alpha = alpha

  if(!is.null(order)) {
    barycoords = barycoords[order,]
  }

  dotplot = geom_point(aes(x=xclip, y=yclip, color=color), data=barycoords)

  if (all(areColors(barycoords$color))) {
    dotplot = c(dotplot, scale_color_identity())
  }

  dotplot
}

#' Creating a dotplot
#'
#' Converts the expression matrix containing 3 biological conditions to barycentric coordinates.
#'
#' @param E Expression matrix
#' @param Gdiffexp Differentially expressed genes
#' @param Goi Genes of interest
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @import ggplot2
#' @export
#'
plotDotplot <- function(E, Gdiffexp=NULL, Goi=NULL)  {
  barycoords = transformBarycentric(E)

  barycoords$goi = rownames(barycoords) %in% Goi
  barycoords$diffexp = rownames(barycoords) %in% Gdiffexp

  order = with(barycoords, order(goi, diffexp))

  alpha = rep(0.5, nrow(E))
  alpha[barycoords$goi] = 1
  alpha[barycoords$diffexp] = 1
  alpha[barycoords$goi & barycoords$diffexp] = 1

  color = rep("#333333", nrow(E))
  color[barycoords$goi] = "green"
  color[barycoords$diffexp] = "red"
  color[barycoords$goi & barycoords$diffexp] = "blue"

  plot = drawHexagonGrid() +
    drawDirections(5, colnames(E)) +
    drawDotplot(barycoords, color=color, order=order, alpha=alpha)

  plot
}
