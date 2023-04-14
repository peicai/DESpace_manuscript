
library(assertthat)
library(scales)
library(concaveman)
library(ggforce)
library(grid)
library(gridExtra)
library(ggnewscale)
#' Spatial plotting functions
#'
#' @param color Optional hex code to set color of borders around spots. Set to
#'   \code{NA} to remove borders.
#' @param ... Additional arguments for \code{geom_polygon()}. \code{size}, to
#'   specify the linewidth of these borders, is likely the most useful.
#' @param platform Spatial sequencing platform. If "Visium", the hex spot layout
#'   will be used, otherwise square spots will be plotted.\cr
#'   NOTE: specifying this argument is only necessary if \code{sce} was not
#'   created by \code{spatialCluster()} or \code{spatialEnhance()}.
#' @param is.enhanced True if \code{sce} contains subspot-level data instead of
#'   spots. Spatial sequencing platform. If true, the respective subspot lattice
#'   for each platform will be plotted.\cr
#'   NOTE: specifying this argument is only necessary if \code{sce} was not
#'   created by \code{spatialCluster()} or \code{spatialEnhance()}.
#'
#' @keywords internal
#' @name spatialPlot
NULL

## Use Seaborn colorblind palette as default
palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
                      "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9")
                      
#' Plot spatial cluster assignments. 
#' 
#' @param sce SingleCellExperiment. If \code{fill} is specified and is a string,
#'   it must exist as a column in \code{colData(sce)}.
#' @param label Labels used to color each spot. May be the name of a column in
#'   \code{colData(sce)}, or a vector of discrete values.
#' @param palette Optional vector of hex codes to use for discrete spot values.
#' @inheritParams spatialPlot
#' 
#' @return Returns a ggplot object.
#' 
#' @examples
#' sce <- exampleSCE()
#' clusterPlot(sce)
#'
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_equal labs theme_void
#' @export
clusterPlot <- function(sce, label="spatial.cluster",
                        palette=NULL, color=NULL,
                        platform=NULL, is.enhanced=NULL,
                        ...) {
  
  if (is.null(platform))
    platform <- .bsData(sce, "platform", "Visium")
  if (is.null(is.enhanced))
    is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
  
  vertices <- .make_vertices(sce, label, platform, is.enhanced)
  
  ## No borders around subspots by default
  if (is.null(color)) {
    color <- ifelse(is.enhanced, NA, "#d8dcd6")
  }
  
  splot <- ggplot(data=vertices, 
                  aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +
    geom_polygon(color=color, ...) +
    labs(fill="Cluster") +
    coord_equal() +
    theme_void()
  
  if (!is.null(palette))
    splot <- splot + scale_fill_manual(values=palette)
  
  splot
}

#' Plot spatial gene expression.
#' 
#' @param sce SingleCellExperiment. If \code{feature} is specified and is a 
#'   string, it must exist as a row in the specified assay of \code{sce}.
#' @param feature Feature vector used to color each spot. May be the name of a
#'   gene/row in an assay of \code{sce}, or a vector of continuous values.
#' @param assay.type String indicating which assay in \code{sce} the expression
#'   vector should be taken from.
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous spot values.
#' @param diverging If true, use a diverging color gradient in
#'   \code{featurePlot()} (e.g. when plotting a fold change) instead of a
#'   sequential gradient (e.g. when plotting expression).
#' @param layer_col Column name of clusters in \code{colData(sce)}
#' @param cluster Cluster names used for drawing a boundary around a group of points (belong to the specify cluster) to drive attention.
#' Can be NULL, "all"/"ALL", and a vector of cluster names.
#' @param legend TRUE of FALSE. Show the legend for the shaped layers 
#' @param label TRUE of FALSE. Adding a label and an arrow pointing to a group.
#' @param ncol The dimensions of the grid to create. By default, if the length of feature equals to 1, set the dimension as 1. Otherwise, set the dimension as 3.
#' @param title TRUE or FALSE. If true, the title name of each (subplot) is the gene name.
#' @inheritParams spatialPlot
#' 
#' @return Returns a ggplot object.
#' 
#' @examples
#' sce <- exampleSCE()
#' featurePlot(sce, "gene_2")
#' 
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void
#' @importFrom scales muted
#' @importFrom assertthat assert_that
#' @export
#' 
#' 
FeaturePlot <- function(sce, feature,
                        assay.type="logcounts", 
                        diverging=FALSE,
                        low=NULL, high=NULL, mid=NULL,
                        color=NULL,
                        platform= "Visium", is.enhanced=FALSE,
                        layer_col = NULL, cluster = NULL,legend_layer = FALSE,
                        label = FALSE, ncol = 3, title = F,
                        legend_exprs = FALSE,title_size = 5,
                        ...) {
  #if (is.null(platform))
  #platform <- .bsData(sce, "platform", "Visium")
  #if (is.null(is.enhanced))
  #is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
  
  ## extract expression from logcounts if a gene name is passed.
  ## otherwise, assume a vector of counts was passed and let
  ## .make_vertices helpers check validity
  if (is.character(feature)) {
    assert_that(all(feature %in% rownames(sce)),
                msg=sprintf("Feature %s not in SCE.", feature))
    fill <- assay(sce, assay.type)[feature, ]
    fill.name <-  feature
  } else {
    fill <- feature
    x <- fill
    ## this could be an argument, but it's easily overwritten with labs()
    ## and we should encourage composing ggplot functions instead
    fill.name <- "Expression"
  }
  
  ## No borders around subspots by default
  if (is.null(color)) {
    color <- ifelse(is.enhanced, NA, "#d8dcd6")
  }
  
  MoreArgs = list(sce, diverging, low, high, mid, color, legend_exprs,
                  platform, is.enhanced, layer_col, cluster, label, title,title_size)
  
  ## if feature is a vector of gene names
  if(length(feature) > 1){
    x <- asplit(as.data.frame(as.matrix(fill)), 1)
    plot.list <- mapply(.geneExprsPlot, x, fill.name, MoreArgs = MoreArgs, SIMPLIFY = FALSE)
    expression.plots <- lapply(plot.list, `[[`, 1) 
    vertices <- lapply(plot.list, `[[`, 2)[[1]]
    widths <- NULL
  }else{
    x <- fill
    plot.list <- .geneExprsPlot(x, fill.name, sce, diverging, low, high, mid, color, legend_exprs,
                                platform, is.enhanced, layer_col, cluster, label, title)
    expression.plots <- list(plot.list$plot)
    vertices <- plot.list$vertices
    ncol=1
    widths <- c(5, 1)
  }
  
  
  ## if 'layer_col' and 'cluster' are specified, and legend_layer == T,  draw the shape of the outline of the group and add the legend of layers
  if(!is.null(layer_col) && !is.null(cluster) && legend_layer){
    ## create the legend of layer
    if(cluster %notin% c("all", "ALL")){
      filter <- vertices$Layer %in% (cluster)
    }else{
      filter <- NULL
    }
    
    my_hist <- vertices %>% 
      ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
      ggforce::geom_mark_hull( aes(x=x.vertex, y=y.vertex, 
                                   color = Layer, fill=Layer, filter = filter))  +
      theme(legend.position = "right") 
    
    legend <- cowplot::get_legend(my_hist)
    plots <- append(expression.plots, list(ggpubr::as_ggplot(legend)) )
    plots <- patchwork::wrap_plots(plots, ncol=min(length(expression.plots),ncol)+1,widths = widths) & 
      theme(legend.position = "bottom") 
    plots
  }else{
    plots <- patchwork::wrap_plots(expression.plots, ncol=min(length(expression.plots),ncol))
    plots
  }
}



#' Make gene expression plots 
#' 
#' @param x fill; A vector of values to use as fill for each spot
#' @param y fill.name; expression legend name 
#' @param sce SingleCellExperiment with row/col in colData
#' @param platform "Visium" or "ST", used to determine spot layout
#' @param is.enhanced If true, \code{sce} contains enhanced subspot data instead
#'   of spot-level expression. Used to determine spot layout.
#' @param layer_col Column name of clusters in \code{colData(sce)}
#' @param cluster Cluster names used for drawing a boundary around a group of points (belong to the specify cluster) to drive attention.
#' Can be NULL, "all"/"ALL", and a vector of cluster names.
#' @param label TRUE of FALSE. Adding a label and an arrow pointing to a group.
#' @param title TRUE or FALSE. If true, the title name of each (subplot) is the gene name.

#' @return Returns a list containing a (list of) ggplot object and a vertices (table of (x.pos, y.pos, spot, fill, Layer)).
#' 
#' @keywords internal
.geneExprsPlot <- function( x, y, sce,
                            diverging,
                            low, high, mid,legend_exprs,
                            color,
                            platform, is.enhanced,
                            layer_col, cluster,
                            label, title,title_size=5,
                            ...){
  fill = as.vector(x)
  fill.name = y
  ## fill.name -> expression legend name
  ## title.name -> subtitle name
  ## if fill.name == 'Expression', don't provide title (title <- FALSE)
  ## if fill.name is gene.name: 
  #### if title == T, title name -> gene name, and the legend name (fill.name) is NULL
  #### if title != T, title name -> NULL, and fill.name is gene name
  if(fill.name == 'Expression'){
    title <-  FALSE
  }
  
  if(title){
    title.name <-  fill.name
    fill.name <- NULL
  }else{
    title.name <- NULL
  }
  
  
  vertices <- .make_vertices(sce, fill, platform, is.enhanced)
  
  if (diverging) {
    low <- ifelse(is.null(low), "#F0F0F0", low)
    mid <- NULL
    high <- ifelse(is.null(high), muted("red"), high)
  } else {
    low <- ifelse(is.null(low), muted("blue"), low)
    mid <- ifelse(is.null(mid), "#F0F0F0", mid)
    high <- ifelse(is.null(high), muted("red"), high)
  }
  
  ## if 'layer_col' and 'cluster' are specified, draw the shape of the outline of the group
  if(!is.null(layer_col) && !is.null(cluster)){
    cdata <- data.frame(colData(sce))
    Layer <- as.character(cdata[[layer_col]])
    if(cluster %notin% c("all", "ALL")){
      Layer[Layer %notin% (cluster)] <- 'Others'}
    Layer = as.factor(Layer)
    vertices <- cbind(vertices, Layer)
    vertices <- vertices %>% filter(!is.na(Layer))
    
    ## annotate clusters
    if(label){
      label <- Layer
    }else{
      label <- NULL
    }
    
    if(cluster %notin% c("all", "ALL")){
      ## Add filter
      splot <- vertices %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
        geom_polygon(aes(group=spot, fill=fill), color=color) +
        labs(fill=fill.name) + coord_equal() +
        theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high)+ 
        new_scale_fill() + new_scale_color()  + 
        ggforce::geom_mark_hull( aes(x=x.vertex, y=y.vertex, 
                                     color = Layer, fill=Layer,
                                     label = label, filter = Layer %in% (cluster)),
                                 alpha=0, expand = unit(0.1, "mm"), show.legend = FALSE) 
    }else{
      splot <- vertices %>% 
        ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
        geom_polygon(aes(group=spot, fill=fill), color=color) +
        labs(fill=fill.name) + coord_equal() +
        theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high)+ 
        new_scale_fill() + new_scale_color()  + 
        ggforce::geom_mark_hull( aes(x=x.vertex, y=y.vertex, 
                                     color = Layer, fill=Layer,
                                     label = label),
                                 alpha=0, expand = unit(0.1, "mm"), show.legend = FALSE) 
    }
    
    
    splot <- splot + scale_alpha(guide="none") + theme_void() + 
      scale_fill_discrete() +
      scale_color_discrete() 
    
  }else{
    # if (color == F) {
    #   color <- ifelse(is.enhanced, NA, "#d8dcd6")
    # }
    splot <-  vertices %>% 
      ggplot(mapping = aes(x=x.vertex, y=y.vertex)) +
      geom_polygon(aes(group=spot, fill=fill), color=color) +
      labs(fill=fill.name) +
      coord_equal() +
      theme_void() + scale_fill_gradient2(low=low, mid=mid, high=high)
    
  }
  if(is.na(legend_exprs)){
    legend_exprs = FALSE
  }
  if(legend_exprs == F){
    splot <- splot + theme(legend.position = "none")
  }
  if(title){
    splot <- splot + labs(title=title.name)+ theme(plot.title = element_text(size=title_size))
  }
  return(list(plot = splot, vertices = vertices))
}




#' Make vertices outlining spots/subspots for geom_polygon()
#' 
#' @param sce SingleCellExperiment with row/col in colData
#' @param fill Name of a column in \code{colData(sce)} or a vector of values to
#'   use as fill for each spot
#' @param platform "Visium" or "ST", used to determine spot layout
#' @param is.enhanced If true, \code{sce} contains enhanced subspot data instead
#'   of spot-level expression. Used to determine spot layout.
#'   
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_vertices <- function(sce, fill, platform, is.enhanced) {
  cdata <- data.frame(colData(sce))
  
  if (platform == "Visium") {
    if (is.enhanced) {
      vertices <- .make_triangle_subspots(cdata, fill)
    } else {
      vertices <- .make_hex_spots(cdata, fill)
    }
  } else if (platform == "ST") {
    if (is.enhanced) {
      vertices <- .make_square_spots(cdata, fill, scale.factor=(1/3))
    } else {
      vertices <- .make_square_spots(cdata, fill)
    }
  } else {
    stop("Unsupported platform: \"", platform, "\". Cannot create spot layout.")
  }
  
  vertices
}

#' Helper to extract x, y, fill ID from colData
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
#' @importFrom assertthat assert_that
.select_spot_positions <- function(cdata, x="col", y="row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  ## I think this is the best way to check if something is a string
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")    
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y)]
    colnames(spot_positions) <- c("x.pos", "y.pos")    
    spot_positions$fill <- fill
  }
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

#' Compute vertex coordinates for each spot in frame of plot
#'
#' @param spot_positions Center for hex, top left for square
#' @param vertex_offsets Data frame of (x, y) offsets wrt spot position for each
#'   vertex of spot
#' 
#' @return Cartesian product of positions and offsets, with coordinates
#'   computed as (pos + offset)
#'
#' @keywords internal
.make_spot_vertices <- function(spot_positions, vertex_offsets) {
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  
  as.data.frame(spot_vertices)
}

#' Make vertices for each hex spot
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_hex_spots <- function(cdata, fill) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  spot_positions <- .adjust_hex_centers(spot_positions)
  
  ## vertices of each hex (with respect to center coordinates)
  ## start at top center, loop clockwise
  vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                               y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  
  spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  
  ## Flip to match image orientation
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  
  spot_vertices
}

#' Adjust hex spot positions so hexagons are adjacent to each other in plot
#'
#' Spots are regular hexagons with one unit of horizontal distance
#' between centers
#' 
#' @return Shifted spot centers
#' 
#' @keywords internal
.adjust_hex_centers <- function(spot_positions) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  ## Start at (1-indexed origin)
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
  
  ## Shift centers up so rows are adjacent
  spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
  
  ## Spot columns are offset by row
  ## (i.e. odd rows have odd numbered columns, even rows have even)
  ## Shift centers to the left so columns are adjacent (but hexes stay offset)
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
  
  spot_positions
}

#' Make vertices for each square spot
#'
#' Squares are simple, just mae a unit square at each array coordinate
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_square_spots <- function(cdata, fill="spatial.cluster", scale.factor=1) {
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  
  vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                               y.offset=c(0, 0, 1, 1))
  vertex_offsets <- vertex_offsets * scale.factor
  
  .make_spot_vertices(spot_positions, vertex_offsets)
}

#' Helper to pull out subspot position columns
#' Probably redundant with select_spot_positions above, but we need subspot.idx
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, "subspot.idx", fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx", "fill")
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y, "subspot.idx")]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx")    
    spot_positions$fill <- fill
  }
  
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

#' Make vertices for each triangle subspot of a hex
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#'
#' @keywords internal
.make_triangle_subspots <- function(cdata, fill="spatial.cluster") {
  spot_positions <- .select_subspot_positions(cdata, x="spot.col", y="spot.row", fill=fill)
  spot_positions <- .adjust_hex_centers(spot_positions)
  
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  ## Make lists of triangle vertices (with respect to hex center)
  ## subspot.idx is same ordering as `shift` in spatialEnhance
  ## that is, beginning in top right and proceeding clockwise, (1, 5, 3, 4, 6, 2)
  ## NOTE: however, we need to reflect over x-axis to match global negation of y-coordinate
  vertex_offsets <- do.call(rbind, list(
    data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=3),
    data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
    data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=1),
    data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=2),
    data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
    data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=4)
  ))
  
  ## note that instead of cartesian product, `merge()` does an outer join
  ## on subspot.idx here
  spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  spot_vertices$y.vertex <- -spot_vertices$y.vertex
  
  spot_vertices
}

.bsData <- function(sce, name, default=NULL, warn=FALSE) {
  if (!exists("BayesSpace.data", metadata(sce)))
    stop("BayesSpace metadata not present in this object.")
  
  bsData <- metadata(sce)[["BayesSpace.data"]]
  if (exists(name, bsData)) {
    bsData[[name]]
  } else {
    if (warn) {
      default.name <- ifelse(is.null(default), "NULL", default)
      warning(name, " not found in BayesSpace metadata. Using default: ", default.name)
    }
    default
  }
}

