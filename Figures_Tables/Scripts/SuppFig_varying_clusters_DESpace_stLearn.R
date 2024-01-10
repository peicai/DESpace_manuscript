rm(list = ls())
path <- "./Simulation/Output/151507/split_mixture_patch_"
library(clue)
library(mclust)
library(RcppHungarian)
library(RColorBrewer)
library(ggplot2)
library(SpatialExperiment)

## use the hungarian algorithm to match BayesSpace clusters to the simulated clusters

getPredLabels <- function(ref_labels, pred_clusters) {
  cost_matrix <- .computeCostMatrix(ref_labels, pred_clusters)
  cluster_map <- .getClusterMapping(cost_matrix)
  
  pred_labels <- unlist(cluster_map[pred_clusters])
  names(pred_labels) <- NULL
  
  return(pred_labels)
}

.computeCostMatrix <- function(ref_labels, pred_clusters) {
  # Create a matrix to store the cost
  unique_ref_labels <- unique(ref_labels)
  unique_pred_clusters <- unique(pred_clusters)
  count_matrix <- matrix(0, nrow = length(unique_ref_labels), ncol = length(unique_pred_clusters), dimnames = list(unique_ref_labels, unique_pred_clusters))
  
  # Iterate over the indices and update the matrix
  for (i in seq_along(ref_labels)) {
    count_matrix[ref_labels[i], pred_clusters[i]] <- count_matrix[ref_labels[i], pred_clusters[i]] + 1
  }
  
  if (ncol(count_matrix) > 1) {
    cost_matrix <- apply(count_matrix, 1, function(row) max(row) - row)
  } else {
    cost_matrix <- t(max(count_matrix) - count_matrix)
  }
  
  return(cost_matrix)
}

.getClusterMapping <- function(cost_matrix) {
  solved <- HungarianSolver(cost_matrix)
  
  cluster_map <- list()
  for (i in 1:nrow(solved$pairs)) {
    from <- rownames(cost_matrix)[solved$pairs[i, 1]]
    to <- if(solved$pairs[i, 2] == 0) from else colnames(cost_matrix)[solved$pairs[i, 2]] 
    cluster_map[[from]] <- to
  }
  
  return(cluster_map)
}
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))

p_all <- lapply(c(2,4,6,8,10,12),function(x){
  res_path <- paste0(path, x, "/probs_0.5_0.9FALSE/stLearn_results.csv")
  stLearn_results <- read.csv(res_path)
  CD1 <- as.data.frame(stLearn_results)
  res_path <- paste0(path, x, "/probs_0.5_0.9FALSE/sce_edgeR.rda")
  load(res_path)
  CD2 <- as.data.frame(colData(sce))
  merged_CD <- merge(CD1, CD2, by.x = "X", by.y = "barcode")
  
  ## match clusters' colour across different number of clusters
  if(x == 2){cols = mycolors}
  if(x == 4){cols = mycolors[c(3,2,4,1)]}
  if(x == 6){cols = mycolors[c(1,3,2,6,4,5)]}
  if(x == 8){cols = mycolors[c(4,8,7,5,1,6,2)]}
  if(x == 10){cols = mycolors[c(8,5,1,7,4,10,2,6,9)]}
  if(x == 12){cols = mycolors[c(8,6,2,10,5,12,7,9,11)]}
  
  p <- ggplot(merged_CD, aes(x = imagecol.x, imagerow.x)) +
    scale_x_reverse() +
    geom_point(aes(colour = factor(louvain)),size = 0.5) +
    theme_classic() +
    scale_color_manual(values = cols) +
    labs(title = paste0("Number of Clusters: ", x)) +
    theme(axis.title.x=element_blank(),
          line = element_blank(),
          text=element_text(size=10),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none")
  return(p)
}
)
pp <- gridExtra::grid.arrange(grobs = p_all,
                              ncol = 3)
path_save = "./Figures/Figures/Supplementary/"
ggsave(filename = paste0(path_save,"SuppFig_clusters_DESpace_stLearn.pdf"),
       plot = pp,
       width = 6,
       height = 5)
