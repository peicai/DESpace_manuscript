rm(list = ls())
source("./Analyses/03_simulations/02_simulation_function_data_specific.R")
## download FULL data:
LIBD_files <- list.files(path = "./DESpace_data/Data/LIBD",
                         pattern = "^15 *.*rda$", full.names = TRUE)

sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", LIBD_files)),1,6)
List_sce_one = list()

# load filtered LIBD data
set.seed(169612)
setwd("./DESpace_data/Data/LIBD")
for(ii in 1:length(sample_names)){
  load(LIBD_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$layer_guess_reordered
  List_sce_one[[ii]] = sce_one
}

library(doParallel)
detectCores()
registerDoParallel(3)

pattern_list <- c('bottom_patch','circle_patch','Manual_clusters_patch','MixClusters_patch','MixClusters_reverse_patch')

########################################### Main simulation (BayesSpace) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','Manual_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,rep(0.9,5))
    SV_percent = c(0.5,rep(0.1,5))
    cluster_index = list(3,4,5,6,7)
  }
  
  default = F
  for(i in c(1,5,9)){ # 3 biological samples
    LIBD_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts','spark_x','meringue','spark','nnSVG'
                    ),
                    spat_methods_params = list(NA,NA,NA,NA,NA
                    ),
                    spat_methods_names = c('SV_edgeR_counts','spark_x',
                                           'meringue','spark','nnSVG'
                    ),
                    layer_names = "spatial.cluster", # Use re-computed BayesSpace clusters
                    original_layer_names = "layer_guess_reordered", # Use manual annotations to construct pattern
                    default = default,
                    reps =1
    )
  }
}

################################# Weak vs. strong (BayesSpace) ##############################################################
## Simulation with 50% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','Manual_clusters_patch') ){
    spatial_probs = c(0.6,0.9)
    SV_percent = 0.5
    cluster_index = NULL
  }else{
    spatial_probs = c(rep(0.6,5),rep(0.9,5))
    SV_percent = c(rep(0.1,10))
    cluster_index = list(3,4,5,6,7,3,4,5,6,7)
  }
  
  default = F
  for(i in c(1,5,9)){ # 3 biological samples
    LIBD_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts','spark_x','meringue','spark','nnSVG'
                    ),
                    spat_methods_params = list(NA,NA,NA,NA,NA
                    ),
                    spat_methods_names = c('SV_edgeR_counts','spark_x',
                                           'meringue','spark','nnSVG'
                    ),
                    layer_names = "spatial.cluster",
                    original_layer_names = "layer_guess_reordered", # Use manual annotations
                    default = default,
                    reps =1
    )
  }
}

########################################### Main simulation (stLearn) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','Manual_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,rep(0.9,5))
    SV_percent = c(0.5,rep(0.1,5))
    cluster_index = list(3,4,5,6,7)
  }
  
  default = F
  for(i in c(1,5,9)){ # 3 biological samples
    LIBD_sub_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts'),
                    spat_methods_params = list(NA),
                    spat_methods_names = c('SV_edgeR_counts'),
                    layer_names =  "stLearn_pca_kmeans", # Use re-computed stLearn clusters
                    original_layer_names = "layer_guess_reordered", # Use manual annotations to construct pattern
                    default = default,
                    reps =1
    )
  }
}

################################# Weak vs. strong (stLearn) ##############################################################
## Simulation with 50% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','Manual_clusters_patch') ){
    spatial_probs = c(0.6,0.9)
    SV_percent = 0.5
    cluster_index = NULL
  }else{
    spatial_probs = c(rep(c(0.6,0.9), each = 5))
    SV_percent = c(rep(0.1,10))
    cluster_index = list(3,4,5,6,7,3,4,5,6,7)
  }
  
  default = F
  for(i in c(1,5,9)){ # 3 biological samples
    LIBD_sub_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/LIBD/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts'),
                    spat_methods_params = list(NA),
                    spat_methods_names = c('SV_edgeR_counts'),
                    layer_names =  "stLearn_pca_kmeans",
                    original_layer_names = "layer_guess_reordered", # Use manual annotations
                    default = default,
                    reps =1
    )
  }
}
