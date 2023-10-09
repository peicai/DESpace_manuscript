rm(list = ls())
source("./Analyses/03_simulations/02_simulation_function_data_specific.R")
## download FULL data:
melanoma_files <- list.files(path = "./DESpace_data/Data/melanoma",
                         pattern = "^mel *.*rda$", full.names = TRUE)
sample_names <- substr(gsub(".rda.*$", "", gsub("^.*/", "", melanoma_files)),1,9)
List_sce_one = list()

# load filtered melanoma data
set.seed(169612)
setwd("./DESpace_data/Data/melanoma")
for(ii in 1:length(sample_names)){
  load(melanoma_files[ii])
  sample <- sample_names[ii]
  layer = colData(sce_one)$spatial.cluster
  List_sce_one[[ii]] = sce_one
}

library(doParallel)
detectCores()
registerDoParallel(3)

pattern_list <- c('right_patch','circle_patch','BayesSpace_clusters_patch','MixClusters_patch','MixClusters_reverse_patch')

########################################### Main simulation (BayesSpace) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 4 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('right_patch','circle_patch','BayesSpace_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,0.9,0.9,0.9)
    SV_percent = c(0.5,1/6,1/6,1/6)
    cluster_index = list(1,c(2,3),4)
  }
  
  default = F
  for(i in c(1,3,5,7)){ # 4 biological samples
    melanoma_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts','spark_x','meringue','spark','nnSVG',
                                     'scranMarker','seuratMarker'
                    ),
                    spat_methods_params = list(NA,NA,NA,NA,NA,NA,NA
                    ),
                    spat_methods_names = c('SV_edgeR_counts','spark_x',
                                           'meringue','spark','nnSVG',
                                           'scranMarker','seuratMarker'
                    ),
                    layer_names = "spatial.cluster", # Use re-computed BayesSpace clusters
                    original_layer_names = "spatial.cluster", # Use manual annotations to construct pattern
                    default = default,
                    reps =1
    )
  }
}

################################# Weak vs. strong (BayesSpace) ##############################################################
## Simulation with 50% of genes are SVGs
## 5 patterns and 4 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('right_patch','circle_patch','BayesSpace_clusters_patch') ){
    spatial_probs = c(0.6,0.9)
    SV_percent = 0.5
    cluster_index = NULL
  }else{
    spatial_probs = c(rep(0.6,3),rep(0.9,3))
    SV_percent = c(rep(1/6,6))
    cluster_index = list(1,c(2,3),4,1,c(2,3),4) # combine BayesSpace cluster 2 and 3 as they're mixed
  }
  
  default = F
  for(i in c(1,3,5,7)){ # 4 biological samples
    melanoma_simulation(i,List_sce_one,
                    working_directory = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                    pattern_type = pattern,
                    my_dir = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                    spatial_probs = spatial_probs,
                    SV_percent = SV_percent,
                    spat_methods = c('SV_edgeR_counts','spark_x','meringue','spark','nnSVG',
                                     'scranMarker','seuratMarker'
                    ),
                    spat_methods_params = list(NA,NA,NA,NA,NA,NA,NA
                    ),
                    spat_methods_names = c('SV_edgeR_counts','spark_x',
                                           'meringue','spark','nnSVG',
                                           'scranMarker','seuratMarker'
                    ),
                    layer_names = "spatial.cluster",
                    original_layer_names = "spatial.cluster", # Use manual annotations
                    default = default,
                    reps =1
    )
  }
}

########################################### Main simulation (stLearn) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 4 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('right_patch','circle_patch','BayesSpace_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,rep(0.9,3))
    SV_percent = c(0.5,rep(1/6,3))
    cluster_index = list(1,c(2,3),4)
  }
  
  default = F
  for(i in c(1,3,5,7)){ # 4 biological samples
    melanoma_sub_simulation(i,List_sce_one,
                        working_directory = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                        pattern_type = pattern,
                        my_dir = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                        spatial_probs = spatial_probs,
                        SV_percent = SV_percent,
                        spat_methods = c('SV_edgeR_counts', 'scranMarker','seuratMarker'),
                        spat_methods_params = list(NA,NA,NA),
                        spat_methods_names = c('SV_edgeR_counts','scranMarker','seuratMarker'),
                        layer_names =  "stLearn_pca_kmeans", # Use re-computed stLearn clusters
                        original_layer_names = "spatial.cluster", # Use manual annotations to construct pattern
                        default = default,
                        reps =1
    )
  }
}

################################# Weak vs. strong (stLearn) ##############################################################
## Simulation with 50% of genes are SVGs
## 5 patterns and 4 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('right_patch','circle_patch','BayesSpace_clusters_patch') ){
    spatial_probs = c(0.6,0.9)
    SV_percent = 0.5
    cluster_index = NULL
  }else{
    spatial_probs = c(rep(c(0.6,0.9), each = 3))
    SV_percent = c(rep(1/6,3))
    cluster_index = list(1,c(2,3),4,1,c(2,3),4)
  }
  
  default = F
  for(i in c(1,3,5,7)){ # 4 biological samples
    melanoma_sub_simulation(i,List_sce_one,
                        working_directory = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                        pattern_type = pattern,
                        my_dir = paste0("./DESpace_data/Simulation/Output/melanoma/",sample_names[i]),
                        spatial_probs = spatial_probs,
                        SV_percent = SV_percent,
                        spat_methods = c('SV_edgeR_counts','scranMarker','seuratMarker'),
                        spat_methods_params = list(NA,NA,NA),
                        spat_methods_names = c('SV_edgeR_counts','scranMarker','seuratMarker'),
                        layer_names =  "stLearn_pca_kmeans",
                        original_layer_names = "spatial.cluster", # Use manual annotations
                        default = default,
                        reps =1
    )
  }
}
