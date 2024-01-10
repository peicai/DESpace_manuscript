#https://derekogle.com/NCGraphing/resources/colors
library(RColorBrewer); library(ggplot2)
all_colours = c(
  "#08306B", # "Manual_DESpace"
  "#A6CEE3", # "BayesSpace_DESpace"
  "#1F78B4", # "StLearn_DESpace"
  "#FFD92F", # "SPARK"
  "#A6761D", # "SPARK-X"
  "#B2DF8A", # "SpatialDE"
  "#33A02C", # "SpatialDE2"
  "#F4A582", # "MERINGUE"
  "#BABABA",#"#666666", # "SpaGCN"
  "#C51B7D",  # "nnSVG"
  "plum2",# BayesSpace_scran_FindMarkers
  "mediumpurple",# StLearn_scran_FindMarkers
  "turquoise",# BayesSpace_seurat_FindAllMarkers
  "turquoise4"# BayesSpace_seurat_FindAllMarkers
)
colours_CompuCost = c(
  "BayesSpace" = "#A63603",#"#A6CEE3",
  "StLearn" = "#A63603"#"#1F78B4"
)

methods_all = c(
  "Manual_DESpace",
  "BayesSpace_DESpace",
  "StLearn_DESpace",
  "SPARK",
  "SPARK-X",
  "SpatialDE",
  "SpatialDE2",
  "MERINGUE",
  "SpaGCN",
  "nnSVG",
  "BayesSpace_findMarkers",
  "StLearn_findMarkers",
  "BayesSpace_FindAllMarkers",
  "StLearn_FindAllMarkers"
)

colors_method = all_colours
names(colors_method) = methods_all

# methods names:
methods_order = c(
  "SV_edgeR_counts",
  "StLearn_edgeR",
  "spark",
  "spark_x",
  "spatialDE",
  "spatialDE2",
  "meringue",
  "SpaGCN",
  "nnSVG",
  "BayesSpace_findMarkers",
  "StLearn_findMarkers",
  "BayesSpace_FindAllMarkers",
  "StLearn_FindAllMarkers"
)
colors_method2 = all_colours[-1]
names(colors_method2) = methods_order
# points, borders:
shape_border = c(0, 0, 0,1, 2, 5, 6, 3, 4, 8, 10, 10, 12, 12)
# points, fill:
shape_fill = c(15, 15, 15, 16, 17, 23, 25, 3, 4, 8, 10, 10, 12, 12)

.prettify <- function(theme = NULL, ...) {
  if (is.null(theme)) theme <- "classic"
  base <- paste0("theme_", theme)
  base <- getFromNamespace(base, "ggplot2")
  base(base_size = 8) + theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(2, "mm"),
    strip.background = element_rect(fill = NA),
    plot.margin = unit(rep(1, 4), "mm"),
    panel.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,1,0,"mm"),
    ...)}

my_theme <- theme(aspect.ratio = 0.67,
                  #legend.box.just = "left",
                  panel.grid = element_blank(),
                  panel.spacing = unit(1, "mm"),
                  panel.border = element_rect(color = "grey"),
                  strip.text = element_text(size = 4),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.5)),
                  axis.text.y=element_text(size=rel(1.5)),
                  axis.title.y = element_text(size=rel(1.5)),
                  axis.title.x = element_text(size=rel(1.5)),
                  axis.title.x.top = element_text(size=rel(1.5)),
                  legend.title=element_text(size=rel(1.5)),
                  legend.text=element_text(size=rel(1.5)),
                  legend.key.width=unit(0.4, "cm"),
                  legend.position="bottom",
                  legend.direction = "horizontal",
                  legend.box="horizontal",
                  strip.text.x = element_text(size = rel(2)),
                  strip.text.y = element_text(size = rel(2)),
                  legend.margin=margin())