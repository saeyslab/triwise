#' @import ggplot2
#' @export
drawHexagonGrid <- function(rmax=5, color="#999999") {
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
    geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color) +
    scale_x_continuous(expand = c(0.1, 0.1))
}

#' @import ggplot2
#' @export
drawCircleGrid <- function(rmax=5, color="#999999") {
  circlepoints <- plyr::ldply(seq(0, pi*2-0.001, pi/100), function(angle) data.frame(x=cos(angle), y=sin(angle)))
  gridData <- dplyr::bind_rows(lapply(seq_len(rmax), function(r) data.frame(r=r, x=circlepoints$x * r, y =circlepoints$y * r)))

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
    geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color)
}

#' @import ggplot2
#' @export
drawDirections <- function(rmax, labels, labelmargin=0.05, type="hexagonal") {
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

    rot = 0
    if(cos(angle) > 0.8) {
      rot = -90
      ha2 = ha
      ha=va
      va = ha2
    } else {
      if(type == "hexagonal") {
         if((cos(angle) < 0) & (sin(angle) > 0)) {
          rot = 60
        } else if ((cos(angle) < 0) & (sin(angle) < 0)) {
          rot = -60
        }
      } else if (type == "circular") {
        if((cos(angle) < 0) & (sin(angle) > 0)) {
          rot = 30
          va = 0
          ha = 0.5
        } else if ((cos(angle) < 0) & (sin(angle) < 0)) {
          rot = -30
          va = 1
          ha = 0.5
        }
      }
    }

    data.frame(x=x, y=y, va=va, ha=ha, xlabel=x*(1+labelmargin),ylabel=y*(1+labelmargin), rot=rot)
  })
  directionsData$label = labels

  c(
    geom_text(aes(label=label, x=xlabel, y=ylabel, vjust=va, hjust=ha, angle=rot), data=directionsData),
    geom_segment(aes(xend=0, x=x, yend=0, y=y, group=label), data=directionsData, alpha=0.5)
  )
}

#' @import ggplot2
#' @export
drawDotplot <- function(barycoords, rmax=5, color="black", alpha=1, order=NULL) {
  barycoords = clipHexagon(barycoords, rmax)
  barycoords$color = color
  barycoords$alpha = alpha

  if(!is.null(order)) {
    barycoords = barycoords[order,]
  }

  dotplot = geom_point(aes(x=xclip, y=yclip, color=color), data=barycoords, size=0.5)

  if (all(areColors(barycoords$color))) {
    dotplot = c(dotplot, scale_color_identity())
  }

  dotplot
}

#' Creating a dotplot
#'
#' Plot a dotplot
#'
#' @param E Expression matrix
#' @param Gdiffexp Differentially expressed genes
#' @param Goi Genes of interest
#' @param colorby Color by differential expression ("diffexp") or by log fold-change ("z")
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @import ggplot2
#' @export
plotDotplot <- function(E, barycoords=NULL, Gdiffexp=NULL, Goi=NULL, Coi=colnames(E), colorby="diffexp") {
  if (is.null(barycoords)) {
    barycoords = transformBarycentric(E)
  }

  barycoords$goi = rownames(barycoords) %in% Goi
  barycoords$diffexp = rownames(barycoords) %in% Gdiffexp

  order = with(barycoords, order(goi, diffexp))

  alpha = rep(0.5, nrow(E))
  alpha[barycoords$goi] = 1
  alpha[barycoords$diffexp] = 1
  alpha[barycoords$goi & barycoords$diffexp] = 1

  rmax=5

  colorpal = colorRamp(brewer.pal(9, "YlOrRd"))

  if (colorby == "diffexp") {
    color = rep("#888888", nrow(E))
    color[barycoords$goi] = "green"
    color[barycoords$diffexp] = "#111111"
    color[barycoords$goi & barycoords$diffexp] = "blue"
  } else if(colorby == "z") {
    color = rgb(colorpal(barycoords$z/max(barycoords$z)), max=255)
  }

  plot = drawHexagonGrid() +
    drawDirections(rmax, Coi) +
    drawDotplot(barycoords, color=color, order=order, alpha=alpha)

  plot
}

#' @export
plotPvalplot <- function(scores, Coi=c("", "", "")) {
  plot = drawCircleGrid() + drawDirections(5, Coi, type="circular")

  scores$r = pmin(-log10(scores$qval), 5)
  scores$x = cos(scores$angle) * scores$r
  scores$y = sin(scores$angle) * scores$r

  plot + geom_point(aes(x=x, y=y), data=scores, size=0.5)
}
