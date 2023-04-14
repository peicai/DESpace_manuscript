rm(list = ls())
source("./Analyses/03_simulations/02_simulation_function_data_specific.R")
## download FULL data:
load("./DESpace_data/Data/mouse_cerebellum/cerebellum_filtered.rda")

library(doParallel)
detectCores()
registerDoParallel(3)

pattern_list <- c('bottom_patch','circle_patch','StLearn_clusters_patch','MixClusters_patch','MixClusters_reverse_patch')

########################################### Main simulation (BayesSpace) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','StLearn_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,rep(0.9,4))
    SV_percent = c(0.5,rep(0.5/4,4))
    cluster_index = list(0,1,2,3,4)
  }
  
  default = F
    Slideseq_simulation(i,List_sce_one,
                        working_directory = paste0("./DESpace_data/Simulation/Output/mouse_cerebellum"),
                        pattern_type = pattern,
                        my_dir = paste0("./DESpace_data/Simulation/Output/mouse_cerebellum"),
                        spatial_probs = spatial_probs,
                        SV_percent = SV_percent,
                        spat_methods = c('SV_edgeR_counts','spark_x','meringue','spark','nnSVG'
                        ),
                        spat_methods_params = list(NA,NA,NA,NA,NA
                        ),
                        spat_methods_names = c('SV_edgeR_counts','spark_x',
                                               'meringue','spark','nnSVG'
                        ),
                        layer_names = "spatial.cluster",# Use re-computed BayesSpace clusters
                        original_layer_names = "stLearn_results.stLearn_pca_kmeans", # Use stLearn clusters to construct pattern
                        default = default,
                        reps =1
    )
}


########################################### Main simulation (stLearn) ###############################################################
## Main simulation with 33% of genes are SVGs
## 5 patterns and 3 samples

foreach(j = c(1:5)) %dopar% {
  pattern <- pattern_list[j]
  if(pattern %in% c('bottom_patch','circle_patch','StLearn_clusters_patch') ){
    spatial_probs = c(0.5,0.9)
    SV_percent = 0.33
    cluster_index = NULL
  }else{
    ## For the other 2 more complicated patterns, i.e., Mixture and inverted mixture
    spatial_probs = c(0.5,rep(0.9,4))
    SV_percent = c(0.5,rep(0.5/4,4))
    cluster_index = list(0,1,2,3)
  }
  
  default = F
    Slideseq_sub_simulation(i,List_sce_one,
                            working_directory = paste0("./DESpace_data/Simulation/Output/mouse_cerebellum"),
                            pattern_type = pattern,
                            my_dir = paste0("./DESpace_data/Simulation/Output/mouse_cerebellum"),
                            spatial_probs = spatial_probs,
                            SV_percent = SV_percent,
                            spat_methods = c('SV_edgeR_counts'),
                            spat_methods_params = list(NA),
                            spat_methods_names = c('SV_edgeR_counts'),
                            layer_names =  "stLearn_pca_kmeans", # Use re-computed stLearn clusters
                            original_layer_names = "stLearn_results.stLearn_pca_kmeans", # Use stLearn clusters to construct pattern
                            default = default,
                            reps =1
    )
}


