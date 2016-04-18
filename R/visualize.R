#' @import ggplot2
#' @export
drawHexagonGrid <- function(rmax=5, color="#999999") {
  hexagonpoints <- plyr::ldply(seq(0, pi*2-0.001, pi/3), function(angle) data.frame(x=cos(angle), y=sin(angle)))
  gridData <- dplyr::bind_rows(lapply(seq_len(rmax), function(r) data.frame(r=r, x=hexagonpoints$x * r, y =hexagonpoints$y * r)))

    drawGridBasis() +
    geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color) +
    scale_x_continuous(expand = c(0.1, 0.1))
}

#' @import ggplot2
#' @export
drawCircleGrid <- function(rmax=1, rstep=0.2, rbase=0, color="#999999", squared=F) {
  circlepoints <- plyr::ldply(seq(0, pi*2-0.001, pi/100), function(angle) data.frame(x=cos(angle), y=sin(angle)))

  radii = seq(rbase, rmax+rbase, rstep)
  if(squared) radii = sqrt(radii)
  gridData <- dplyr::bind_rows(lapply(radii, function(r) data.frame(r=r, x=circlepoints$x * r, y =circlepoints$y * r)))

  drawGridBasis() + geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color)
}

#' @import ggplot2
#' @export
drawGridBasis <- function() {
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
    )
}

#' @import ggplot2
#' @title plotRoseplot
#' @name plotRoseplot
#' @description A rose plot shows the distribution of the differentially expressed genes in different directions.
#' @export
plotRoseplot = function(barycoords, Gdiffexp, Goi=rownames(barycoords), conditions=attr(barycoords, "conditions"), nbins=12, bincolors=rainbow(nbins, start=0, v=0.8, s=0.6)) {
  barycoords = triwise::transformBarycentric(Eoi)

  deltaalpha = pi*2/nbins

  Goidiffexp = intersect(Goi, Gdiffexp)
  percnochange = (1- length(Goidiffexp)/length(Goi))
  if(is.na(percnochange)) {percnochange=0}

  bins = seq(0, pi*2-0.0001, by=deltaalpha)
  if (length(Goidiffexp) > 0) {
    binned = sapply(barycoords[Goidiffexp, "angle"], function(x) {
      binid = which.min(abs(sapply(bins, diffCircular, x)))
      if (is.na(binid)) {
        binid = 1
      }
      binid
    })

    bincounts = tabulate(binned, nbins=nbins)

  } else {
    bincounts = rep(0,nbins)
  }
  percschange = bincounts/length(Goidiffexp)

  rmax = 1#ceiling(max(percschange) * 10)/10

  plot = drawCircleGrid(rmax, 0.1, squared=T)

  circlesect = function(angle1=0, angle2=pi*2, r=1) {
    points = plyr::ldply(seq(angle1, angle2, (angle2-angle1)/100), function(angle) data.frame(x=cos(angle)*r, y=sin(angle)*r))
    points = rbind(points, data.frame(x=0,y=0))
  }

  for (binid in 1:length(bins)) {
    angle1 = bins[binid] - deltaalpha/2
    angle2 = (angle1 + deltaalpha)

    perchange = percschange[binid]

    radius = sqrt(perchange)

    sector = circlesect(angle1, angle2, r = radius)

    plot = plot + geom_polygon(data=sector, aes(x,y), fill=bincolors[binid])
  }
  plot + drawDirections(rmax, conditions, type="circular")
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
drawDotplot <- function(barypoints, rmax=5, color=scale_color_grey(), alpha=scale_alpha(), size=scale_size(), order=NULL) {
  barypoints = clipHexagon(barypoints, rmax)

  if(!is.null(order)) {
    barypoints = barypoints[order,]
  }

  dotplot = geom_point(aes(x=xclip, y=yclip, color=colorby, alpha=alphaby, size=sizeby), data=barypoints)
  dotplot = c(dotplot, color, alpha, size)

  dotplot
}

#' Creating a dotplot
#'
#' Plot a dotplot
#'
#' @param barycoords Dataframe containing for every gene its barycentric coordinates, as return by \code{transformBarycentric}
#' @param Gdiffexp Differentially expressed genes
#' @param gsets Genes of interest, a character or numeric vector to plot one set of genes, a named list to plot multiple gene lists
#' @param conditions Character vector specifying the names of the three biological conditions
#' @param colorby Color by differential expression ("diffexp") or by log fold-change ("z")
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @import ggplot2
#' @export
plotDotplot <- function(barycoords, diffexp=NULL, gsets=NULL, conditions=attr(barycoords, "conditions"), colorby="diffexp", colors=NULL, rmax=5) {
  if (!is.list(gsets)) {
    gsets = list(gset=gsets)
  }

  barypoints = as.data.frame(barycoords)
  colnames(barypoints) = c("x", "y")
  barypoints = addPolar(barypoints)

  barypoints$diffexp = rownames(barypoints) %in% diffexp
  barypoints$ingset = F
  barypoints$gsetname = "all"
  for (gsetname in names(gsets)) {
    barypoints[gsets[[gsetname]], "ingset"] = T
    barypoints[gsets[[gsetname]], "gsetname"] = gsetname
  }
  barypoints$type = paste0(ifelse(barypoints$diffexp, "diff", "nodiff"), barypoints$gsetname)

  if (is.null(colors)) {
    if (colorby == "diffexp") {
      # make colorvalues in two steps so that extra colors (eg. if 1 gene set) are filled with nas
      
      # different palletes for diffexp and nodiffexp
      colorvalues = setNames(c("#333333", RColorBrewer::brewer.pal(max(length(gsets), 3), "Set1")), c("diffall", paste0("diff",names(gsets))))
      colorvalues = c(colorvalues, setNames(c("#AAAAAA", RColorBrewer::brewer.pal(max(length(gsets), 3), "Pastel1")), c("nodiffall", paste0("nodiff", names(gsets)))))
      
      # sample palette
      #colorvalues = setNames(c("#222222", RColorBrewer::brewer.pal(max(length(gsets), 3), "Set1")), c("all", names(gsets)))
      color = scale_colour_manual(
        values=colorvalues
      ,name="type")
      barypoints$colorby = barypoints$type
      barypoints$alphaby = barypoints$diffexp
      barypoints$sizeby = barypoints$ingset
      alpha = scale_alpha_manual(values=setNames(c(0.8,0.8), c(F,T)), name="diffexp")
      size = scale_size_manual(values=setNames(c(0.5, 2), c(F,T)), name="ingset")
    } else if (colorby == "z") {
      color = scale_colour_continuous()
      barypoints$colorby = barypoints$z
      barypoints$alphaby = T
      alpha = scale_alpha_discrete(breaks=c(F,T), limits=c(0.4,1))
    } else {
      color = scale_colour_continuous()
      barypoints$colorby = barypoints$r
      barypoints$alphaby = T
      alpha = scale_alpha_discrete(breaks=c(F,T), limits=c(0.4,1))
    }
  }

  order = with(barypoints, order(ingset, diffexp))

  plot = drawHexagonGrid() +
    drawDirections(rmax, conditions) +
    drawDotplot(barypoints, rmax, color=color, order=order, alpha=alpha, size=size)

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
