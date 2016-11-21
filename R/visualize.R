#' @import ggplot2
drawHexagonGrid <- function(rmax=5, color="#999999", showlabels=T, baseangle=0) {
  radii = seq_len(rmax)
  hexagonpoints <- plyr::ldply(seq(0, pi*2-0.001, pi/3), function(angle) data.frame(x=cos(angle+baseangle), y=sin(angle+baseangle)))
  gridData <- dplyr::bind_rows(lapply(radii, function(r) data.frame(r=r, x=hexagonpoints$x * r, y =hexagonpoints$y * r)))

  plot =  drawGridBasis() +
    geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color) +
    scale_x_continuous(expand = c(0.1, 0.1))

  if (showlabels) {
    labelData <- data.frame(r=radii, label=2^radii, y=sapply(radii, function(x) hexagonPolar(pi/2, x)), x=0)
    plot = plot + geom_label(aes(x=x, y=y, label=label), hjust=0.5, data=labelData)
  }
  plot
}

#' @import ggplot2
drawCircleGrid <- function(rmax=1, rstep=0.2, rbase=rstep, color="#999999", squared=F, showlabels=T, labeller=function(x) x) {
  circlepoints <- plyr::ldply(seq(0, pi*2-0.001, pi/100), function(angle) data.frame(x=cos(angle), y=sin(angle)))

  radii = seq(rbase, rmax, rstep)
  labels = labeller(radii)
  if(squared) radii = sqrt(radii)
  gridData <- dplyr::bind_rows(lapply(radii, function(r) data.frame(r=r, x=circlepoints$x * r, y =circlepoints$y * r)))

  plot = drawGridBasis() + geom_polygon(aes(x=x, y=y, group=r), gridData, fill=NA, colour=color)

  if (showlabels) {
    labelData <- data.frame(r=radii, label=unlist(labels), y=radii, x=0)
    plot = plot + geom_label(aes(x=x, y=y, label=label), hjust=0.5, data=labelData)
  }
  plot
}

#' @import ggplot2
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
#' @description A rose plot shows the distribution of a given set of genes in different directions.
#' @param barycoords Dataframe containing barycentric coordinates as returned by `transformBarycentric`
#' @param Gdiffexp List of differentially expressed genes
#' @param Goi List of genes of interest
#' @param size Should the `radius` or the `surface` of a circle sector denote the number of genes differentially expressed in a particular direction
#' @param relative Whether to show the relative number of genes or the absolute number of genes
#' @param showlabels Whether to label the grid
#' @param Coi Names of the three biological conditions, used for labelling
#' @param nbins Number of bins, should be a multiple of 3 to make sense
#' @param bincolors Colors of every bin, defaults to a rainbow palette
#' @param rmax Number or "auto" (default), denotes the maximal radius of the grid.
#' @param baseangle The angle by which to rotate the whole plot (default to 0)
#' @return A ggplot2 plot, which can be used to further customize the plot
#' @export
plotRoseplot = function(barycoords, Gdiffexp=rownames(barycoords), Goi=rownames(barycoords), size="surface", relative=T, showlabels=T, Coi=attr(barycoords, "conditions"), nbins=12, bincolors=rainbow(nbins, start=0, v=0.8, s=0.6), rmax="auto", baseangle=0) {
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
  if(relative) {
    percschange = bincounts/length(Goidiffexp)
  } else {
    percschange = bincounts
  }

  if (rmax == "auto") {
    #rmax = 10^(ceiling(log10(max(percschange))))
    rmax = roundUpNice(max(percschange))
    ticks = grDevices::axisTicks(c(0, rmax), F, nint=4)
    rmax = ticks[length(ticks)]
  }
  ticks = grDevices::axisTicks(c(0, rmax), F, nint=4)

  labeller = if(relative) {scales::percent} else {itself}
  plot = drawCircleGrid(rmax, ticks[[2]] - ticks[[1]], squared=(size=="surface"), showlabels=showlabels, labeller=labeller)

  circlesect = function(angle1=0, angle2=pi*2, r=1) {
    points = plyr::ldply(seq(angle1, angle2, (angle2-angle1)/100), function(angle) data.frame(x=cos(angle)*r, y=sin(angle)*r))
    points = rbind(points, data.frame(x=0,y=0))
  }

  for (binid in 1:length(bins)) {
    angle1 = bins[binid] - deltaalpha/2
    angle2 = (angle1 + deltaalpha)

    perchange = percschange[binid]

    radius = ifelse(size=="radius", perchange, sqrt(perchange))

    sector = circlesect(angle1+baseangle, angle2+baseangle, r = radius)

    plot = plot + geom_polygon(data=sector, aes(x,y), fill=bincolors[binid])
  }
  plot + drawDirections(ifelse(size=="surface", sqrt(rmax), rmax), Coi, type="circular", baseangle=baseangle)
}

#' @import ggplot2
drawDirections <- function(rmax, labels, labelmargin=0.05, type="hexagonal", baseangle=0) {
  directionsData = plyr::ldply(seq(0, pi*2-0.001, pi/3*2), function(angle) {
    angle = angle + baseangle
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
    if(sin(angle) >= 0) {
      rot = angle - pi/2
      va = 0
      ha = 0.5
    } else {
      rot = angle + pi/2
      va = 1
      ha = 0.5
    }

    rot = rot / pi * 180

    data.frame(x=x, y=y, va=va, ha=ha, xlabel=x*(1+labelmargin),ylabel=y*(1+labelmargin), rot=rot)
  })
  directionsData$label = labels

  c(
    geom_text(aes(label=label, x=xlabel, y=ylabel, vjust=va, hjust=ha, angle=rot), data=directionsData),
    geom_segment(aes(xend=0, x=x, yend=0, y=y, group=label), data=directionsData, alpha=0.5)
  )
}

#' @import ggplot2
drawDotplot <- function(barypoints, rmax=5, color=scale_color_grey(), alpha=scale_alpha(), size=scale_size(), order=NULL, baseangle=0) {
  barypoints = clipHexagon(barypoints, rmax, baseangle)

  if(!is.null(order)) {
    barypoints = barypoints[order,]
  }

  dotplot = geom_point(aes(x=xclip, y=yclip, color=colorby, alpha=alphaby, size=sizeby), data=barypoints)
  dotplot = c(dotplot, color, alpha, size)

  dotplot
}

#' @import ggplot2
drawConnectionplot <- function(barypoints, barypoints2, rmax=5, order=NULL, baseangle=0) {
  barypoints = clipHexagon(barypoints, rmax, baseangle)
  barypoints2 = clipHexagon(barypoints2, rmax, baseangle)

  colnames(barypoints2) = paste0(colnames(barypoints2), "2")

  allbarypoints = data.frame(barypoints, barypoints2)

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
#' @param Goi Genes of interest, a character or numeric vector to plot one set of genes, a named list containing different such vectors to plot multiple gene sets
#' @param Coi Character vector specifying the names of the three biological conditions, used for labelling
#' @param colorby Color by differential expression ("diffexp") or by log fold-change ("z")
#' @param colorvalues Color values, different syntax depending on the colorby parameter: \itemize{
#'   \item diffexp: Named list with colors. First part of the name denotes whether a gene is differentially expressed (`diffexp` or `nodiffexp`). The second part denotes the name of the gene set. Genes not within a gene set are denoted by `all`. \cr
#'   For example: list(diffexpall="#000000", nodiffexpall="#AAAAAA", nodiffexpgset="#FFAAAA", diffexpgset="#FF0000")
#'   \item z: A character vector of color values, used to generate the color gradient
#' }
#' @param rmax Number denoting the maximal radius of the grid. All points outside of the grid will be clipped on the boundaries.
#' @param showlabels Whether to show labels on the grid
#' @param sizevalues Named list with the size of each dot if differentially expressed (TRUE) or not (FALSE)
#' @param alphavalues Named list with the alpha value of each dot if differentially expressed (TRUE) or not (FALSE)
#' @param barycoords2 Dataframe containing for every gene a second set of barycentric coordinates, as returned by \code{transformBarycentric}. An arrow will be drawn from the coordinates in `barycoords` to those in `barycoords2`
#' @param baseangle The angle by which to rotate the whole plot (default to 0)
#' @examples
#' data(vandelaar)
#' Eoi_replicates <- vandelaar[, phenoData(vandelaar)$celltype %in% c("BM_mono", "FL_mono", "YS_MF")]
#' Eoi <- limma::avearrays(Eoi_replicates, phenoData(Eoi_replicates)$celltype)
#' Eoi = Eoi[,c("YS_MF", "FL_mono", "BM_mono")]
#' barycoords = transformBarycentric(Eoi)
#' plotDotplot(barycoords)
#'
#' @return A ggplot2 plot, which can be used for further customization
#' @import ggplot2
#' @export
plotDotplot <- function(barycoords, Gdiffexp=rownames(barycoords), Goi=NULL, Coi=attr(barycoords, "conditions"), colorby="diffexp", colorvalues=NULL, rmax=5, showlabels=T, sizevalues=setNames(c(0.5,2), c(F,T)), alphavalues=setNames(c(0.8, 0.8), c(F, T)), barycoords2=NULL, baseangle=0) {
  if (!is.list(Goi)) {
    Goi = list(gset=Goi)
  }

  barypoints = as.data.frame(barycoords)
  #colnames(barypoints) = c("x", "y")
  #barypoints = addPolar(barypoints)

  barypoints$diffexp = rownames(barypoints) %in% Gdiffexp
  barypoints$ingset = F
  barypoints$gsetname = "all"
  if(is.null(names(Goi))) names(Goi) = 1:length(Goi)
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
    color = scale_colour_manual(values=colorvalues,name="type")
    barypoints$colorby = barypoints$type
    barypoints$alphaby = barypoints$diffexp
    barypoints$sizeby = barypoints$ingset
  } else if (colorby == "z") {
    if(!("z" %in% colnames(barypoints))) stop("z column not defined")

    color = scale_colour_continuous()
    barypoints$colorby = barypoints$z
    barypoints$alphaby = T
    barypoints$sizeby = T
  } else {
    color = scale_colour_continuous()
    barypoints$colorby = barypoints$r
    barypoints$alphaby = T
    barypoints$sizeby = T
  }
  alpha = scale_alpha_manual(values=setNames(c(0.8,0.8), c(F,T)), name="diffexp")
  size = scale_size_manual(values=setNames(c(0.5,2), c(F,T)), name="ingset")

  order = with(barypoints, order(ingset, diffexp))

  plot = drawHexagonGrid(rmax, showlabels=showlabels, baseangle=baseangle) +
    drawDirections(rmax, Coi, baseangle=baseangle) +
    drawDotplot(barypoints, rmax, color=color, order=order, alpha=alpha, size=size, baseangle=baseangle)

  if(!is.null(barycoords2) && sum(barypoints$ingset) > 0) {
    barypoints2 = as.data.frame(barycoords2)[rownames(barypoints),]
    #colnames(barypoints2) = c("x", "y")
    #barypoints2 = addPolar(barypoints2)

    plot = plot + drawConnectionplot(barypoints[barypoints$ingset, ], barypoints2[barypoints$ingset, ], rmax, baseangle=baseangle)
  }

  plot
}

#' Plot results from unidirectional enrichment
#' @param scores Dataframe as returned by `testUnidirectionality` containing for every gene set a q-value and an associated angle
#' @param Coi Names of the three biological conditions
#' @param colorby Column in scores used for coloring
#' @param showlabels Whether to show labels on the grid
#' @param baseangle The angle by which to rotate the whole plot (default to 0)
#' @export
plotPvalplot <- function(scores, Coi=c("", "", ""), colorby=NULL, showlabels=T, baseangle=0) {
  labeller = function(x) paste0("10^-", x, "")
  plot = drawCircleGrid(5, 1, showlabels=showlabels, labeller=labeller) + drawDirections(5, Coi, type="circular", baseangle=baseangle)

  scores$r = pmin(-log10(scores$qval), 5)
  scores$x = cos(scores$angle+baseangle) * scores$r
  scores$y = sin(scores$angle+baseangle) * scores$r

  if (!is.null(colorby)) {
    scores$colorby = scores[,colorby]
    color = scale_color_continuous()
  } else {
    scores$colorby = factor(1)
    color = scale_color_manual(values=c(`1`="#333333"))
  }

  plot = plot + geom_point(aes(x=x, y=y, color=colorby), data=scores, size=1) + color

  plot
}
