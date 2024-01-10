rm(list = ls())
path <- "./Simulation/Output/151507/split_mixture_patch_"
library(clue)
library(mclust)
library(RcppHungarian)

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

p_all <- lapply(c(2,4,6,8,10,12),function(x){
  print(x)
  res_path <- paste0(path, x, "/probs_0.5_0.9FALSE/sce_edgeR.rda")
  load(res_path)
  CD <- as.data.frame(colData(sce))
  mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
  if(x != 10) {
    CD$mapped_clusters <- getPredLabels(CD$pattern, CD$spatial.cluster)
  }else {
    CD$mapped_clusters <- CD$spatial.cluster
    mycolors <- c("#D95F02", "#A6CEE3", "#66A61E", "#1B9E77", "#E6AB02", "#666666", 
                  "#1F78B4", "#A6761D", "#E7298A" )
  }
  p <- ggplot(CD, aes(x = imagecol, imagerow)) +
    geom_point(aes(colour = factor(mapped_clusters)),size = 0.5) +
    scale_x_reverse() +
    theme_classic() +
    scale_color_manual(values = mycolors) +
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
ggsave(filename = paste0(path_save,"SuppFig_clusters_DESpace_BayesSpace.pdf"),
       plot = pp,
       width = 6,
       height = 5)
