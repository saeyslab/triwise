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
#' @param barycoords Dataframe contain barycentric coordinates as returned by `transformBarycentric`
#' @param Gdiffexp List of differentially expressed genes
#' @param Goi List of genes of interest
#' @param size Should the `radius` or the `surface` of a circle sector denote the number of genes differentially expressed in a particular direction
#' @param Coi Names of the three biological conditions
#' @param nbins Number of bins, should be a multiple of 3 to make sense
#' @param bincolors Colors of every bin
#' @export
plotRoseplot = function(barycoords, Gdiffexp=rownames(barycoords), Goi=rownames(barycoords), size="surface", Coi=attr(barycoords, "conditions"), nbins=12, bincolors=rainbow(nbins, start=0, v=0.8, s=0.6)) {
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


  plot = drawCircleGrid(rmax, 0.2, squared=(size=="surface"))

  circlesect = function(angle1=0, angle2=pi*2, r=1) {
    points = plyr::ldply(seq(angle1, angle2, (angle2-angle1)/100), function(angle) data.frame(x=cos(angle)*r, y=sin(angle)*r))
    points = rbind(points, data.frame(x=0,y=0))
  }

  for (binid in 1:length(bins)) {
    angle1 = bins[binid] - deltaalpha/2
    angle2 = (angle1 + deltaalpha)

    perchange = percschange[binid]

    radius = ifelse(size=="radius", perchange, sqrt(perchange))

    sector = circlesect(angle1, angle2, r = radius)

    plot = plot + geom_polygon(data=sector, aes(x,y), fill=bincolors[binid])
  }
  plot + drawDirections(rmax, Coi, type="circular")
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

#' @import ggplot2
#' @export
drawConnectionplot <- function(barypoints, barypoints2, rmax=5, order=NULL) {
  barypoints = clipHexagon(barypoints, rmax)
  barypoints2 = clipHexagon(barypoints2, rmax)

  colnames(barypoints2) = paste0(colnames(barypoints2), "2")

  allbarypoints =data.frame(barypoints, barypoints2)

  if(!is.null(order)) {
    allbarypoints = allbarypoints[order,]
  }

  dotplot = geom_segment(aes(x=xclip, y=yclip, xend=xclip2, yend=yclip2, color=colorby), data=allbarypoints, arrow=arrow(type="closed", length=unit(0.05, "inches")))

  dotplot
}

#' Creating a dotplot
#'
#' Plot a dotplot
#'
#' @param barycoords Dataframe containing for every gene its barycentric coordinates, as return by \code{transformBarycentric}
#' @param Gdiffexp Differentially expressed genes
#' @param Goi Genes of interest, a character or numeric vector to plot one set of genes, a named list to plot multiple gene lists
#' @param Coi Character vector specifying the names of the three biological conditions
#' @param colorby Color by differential expression ("diffexp") or by log fold-change ("z")
#' @param colorvalues Colors used according to colorby
#' @param sizevalues List with the size of each dot if differentially expressed (`T`) or not (`F`)
#' @param alphavalues List with the alpha value of each dot if differentially expressed or not
#' @param barycoords2 Dataframe containing for every gene a second set of barycentric coordinates, as returned by \code{transformBarycentric}. An arrow will be drawn from the coordinates in `barycoords` to those in `barycoords2`
#' @return Dataframe containing for every original point its new x and y coordinates in barycentric space
#' @import ggplot2
#' @export
plotDotplot <- function(barycoords, Gdiffexp=rownames(barycoords), Goi=NULL, Coi=attr(barycoords, "conditions"), colorby="diffexp", colorvalues=NULL, rmax=5, sizevalues=c(T=2, F=0.5), alphavalues=c(T=0.8, F=0.8), barycoords2=NULL) {
  if (!is.list(Goi)) {
    Goi = list(gset=Goi)
  }

  barypoints = as.data.frame(barycoords)
  #colnames(barypoints) = c("x", "y")
  #barypoints = addPolar(barypoints)

  barypoints$diffexp = rownames(barypoints) %in% Gdiffexp
  barypoints$ingset = F
  barypoints$gsetname = "all"
  for (gsetname in names(Goi)) {
    barypoints[Goi[[gsetname]], "ingset"] = T
    barypoints[Goi[[gsetname]], "gsetname"] = gsetname
  }
  barypoints$type = paste0(ifelse(barypoints$diffexp, "diff", "nodiff"), barypoints$gsetname)

  if (colorby == "diffexp") {
    if (is.null(colorvalues)) {
      # make colorvalues in two steps so that extra colors (eg. if 1 gene set) are filled with nas

      # different palletes for diffexp and nodiffexp
      colorvalues = setNames(c("#333333", RColorBrewer::brewer.pal(max(length(Goi), 3), "Set1")), c("diffall", paste0("diff",names(Goi))))
      colorvalues = c(colorvalues, setNames(c("#AAAAAA", RColorBrewer::brewer.pal(max(length(Goi), 3), "Pastel1")), c("nodiffall", paste0("nodiff", names(Goi)))))
    }

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

  order = with(barypoints, order(ingset, diffexp))

  plot = drawHexagonGrid() +
    drawDirections(rmax, Coi) +
    drawDotplot(barypoints, rmax, color=color, order=order, alpha=alpha, size=size)

  if(!is.null(barycoords2) && sum(barypoints$ingset) > 0) {
    barypoints2 = as.data.frame(barycoords2)
    #colnames(barypoints2) = c("x", "y")
    #barypoints2 = addPolar(barypoints2)

    plot = plot + drawConnectionplot(barypoints[barypoints$ingset, ], barypoints2[barypoints$ingset, ], rmax)
  }

  plot
}

#' Plot results from unidirectional enrichment
#' @param scores Dataframe as returned by `testUnidirectionality` containing for every gene set a q-value and an associated angle
#' @param Coi Names of the three biological conditions
#' @param colorby Column in scores used for coloring
#' @export
plotPvalplot <- function(scores, Coi=c("", "", ""), colorby=NULL) {
  plot = drawCircleGrid(5, 1) + drawDirections(5, Coi, type="circular")

  scores$r = pmin(-log10(scores$qval), 5)
  scores$x = cos(scores$angle) * scores$r
  scores$y = sin(scores$angle) * scores$r

  if (!is.null(colorby)) {
    scores$colorby = scores[,colorby]
  } else {
    scores$colorby = 1
  }

  plot = plot + geom_point(aes(x=x, y=y, color=colorby), data=scores, size=5)

  plot
}
